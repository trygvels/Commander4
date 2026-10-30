"""Tests of compsep/index_field.py: amplitude priors as pseudo-bands, region maps, per-region
mixing matrices and the prior-aware per-pixel amplitude draw (M2.4 stage 1)."""
from copy import deepcopy

import healpy as hp
import numpy as np
import pytest
from mpi4py import MPI
from pixell import curvedsky
from pixell.bunch import Bunch

from commander4.backend import compsep as cpp_compsep
from commander4.compsep import index_field as ifd
from commander4.compsep.perpix_solver import solve_compsep_perpix, solve_perpix_with_priors
from commander4.data_models.detector_map import DetectorMap
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust

NSIDE, NSIDE_BETA = 4, 1
NPIX = 12 * NSIDE**2
FREQS = np.array([30.0, 44.0, 70.0, 100.0, 143.0, 353.0])


def _component(cls, shortname, **sed):
    params = Bunch(polarization="I", shortname=shortname, spatially_varying_MM=False,
                   Cl_prior_amplitude=None, lmax=2, units="uK_RJ", **sed)
    object.__setattr__(params, "_name", shortname)
    return cls(params, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True,
               eval_pol="I", comp_name=shortname)


def _components():
    return [_component(CMB, "cmb"),
            _component(Synchrotron, "sync", beta=-3.0, nu_ref=30.0),
            _component(ThermalDust, "dust", beta=1.55, T=20.0, nu_ref=353.0)]


def _stream(seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(FREQS.size, NPIX)), rng.uniform(0.25, 4.0, (FREQS.size, NPIX))


# --- AmpPrior ---------------------------------------------------------------------------------

def test_amp_prior_validation() -> None:
    with pytest.raises(ValueError):
        ifd.AmpPrior("flat")
    with pytest.raises(ValueError):
        ifd.AmpPrior("gaussian", mean=0.0)
    with pytest.raises(ValueError):
        ifd.AmpPrior("gaussian", mean=0.0, rms=0.0)
    with pytest.raises(ValueError):
        ifd.AmpPrior("hier", template=np.array([1.0, -1.0]))
    with pytest.raises(ValueError):
        ifd.AmpPrior("hier")


def test_amp_prior_pseudo_bands() -> None:
    m, w = ifd.AmpPrior("none").pseudo_band(5)
    np.testing.assert_array_equal(m, 0.0)
    np.testing.assert_array_equal(w, 0.0)
    m, w = ifd.AmpPrior("gaussian", mean=3.0, rms=0.5).pseudo_band(5)
    np.testing.assert_array_equal(m, 3.0)
    np.testing.assert_allclose(w, 4.0)
    g = np.array([1.0, 2.0, 4.0])
    prior = ifd.AmpPrior("hier", template=g, mu=1.5, log_s=np.log(0.25))
    m, w = prior.pseudo_band(3)
    np.testing.assert_allclose(m, 1.5 * g)
    np.testing.assert_allclose(w, 1.0 / (0.25 * g)**2)
    with pytest.raises(ValueError):          # not a HEALPix map: cannot be regridded
        prior.pseudo_band(4)
    # A template at another nside is ud_graded to the stream (children take the parent value).
    g12 = 1.0 + np.arange(12.0)
    prior = ifd.AmpPrior("hier", template=g12, mu=1.0, log_s=0.0)
    m, w = prior.pseudo_band(48)
    np.testing.assert_allclose(m, hp.ud_grade(g12, 2))
    np.testing.assert_allclose(w, 1.0 / hp.ud_grade(g12, 2)**2)
    used = np.ones(48, dtype=bool)
    grad = prior.hyper_grad(np.ones(48), np.zeros(48), used)
    assert np.isfinite(grad).all()
    m48, w48 = ifd.AmpPrior("gaussian", mean=g12, rms=g12).pseudo_band(48)
    np.testing.assert_allclose(m48, hp.ud_grade(g12, 2))


def test_amp_prior_trailing_terms_and_hyperprior() -> None:
    g = np.array([1.0, 2.0, 4.0, 8.0])
    prior = ifd.AmpPrior("hier", template=g, mu=0.8, log_s=-1.0)
    used = np.array([True, False, True, True])
    m, inv_T = prior.pseudo_band(4)
    expected = -0.5 * np.sum((m**2 * inv_T)[used]) + 0.5 * np.sum(np.log(inv_T[used]))
    assert prior.trailing_terms(used) == pytest.approx(expected, rel=1e-14)
    assert ifd.AmpPrior("none").trailing_terms(used) == 0.0
    assert prior.log_hyperprior() == pytest.approx(
        -0.5 * ((0.8 - 1.0) / 0.5)**2 - 0.5 * ((-1.0 - np.log(0.3)) / 1.0)**2)
    prior.log_s = ifd.HIER_LOGS_MIN - 1e-9
    assert prior.log_hyperprior() == -np.inf
    assert ifd.AmpPrior("none").log_hyperprior() == 0.0
    h = 1e-6
    prior = ifd.AmpPrior("hier", template=g, mu=0.8, log_s=-1.0)
    for name, idx in (("mu", 0), ("log_s", 1)):
        value = getattr(prior, name)
        setattr(prior, name, value + h)
        up = prior.log_hyperprior()
        setattr(prior, name, value - h)
        down = prior.log_hyperprior()
        setattr(prior, name, value)
        assert prior.log_hyperprior_grad()[idx] == pytest.approx((up - down) / (2 * h), rel=1e-6)


def test_hier_hyper_grad_matches_finite_differences_of_the_augmented_target() -> None:
    """dF/dmu and dF/dlog s from the kernel extras equal central differences of
    F = sum_p [1/2 b'^T A'^-1 b' - 1/2 log det A'] + trailing terms, with every pixel used."""
    maps_sky, inv_var = _stream(seed=3)
    comps = _components()
    rng = np.random.default_rng(4)
    g = np.exp(rng.normal(size=NPIX))
    comps[1].amp_prior = ifd.AmpPrior("hier", template=g, mu=1.2, log_s=-0.9)
    comps[2].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=50.0)
    comps[0].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=70.0)
    region = ifd.region_map(0, NSIDE)
    M_aug, dM_aug = ifd.augment_mixing(*ifd.mixing_regions(FREQS, comps, 1, True))

    def target(extras=False):
        sky_aug, w_aug = ifd.augment_stream(maps_sky, inv_var, comps)
        out = cpp_compsep.marginal_loglike_regions(sky_aug, w_aug, region, M_aug, dM_aug, 1,
                                                   want_extras=extras)
        used = np.ones(NPIX, dtype=bool) if not extras else out[5]
        value = float(out[0].sum()) + comps[1].amp_prior.trailing_terms(used)
        return (value, out) if extras else value

    value, out = target(extras=True)
    assert out[2] == 0 and out[5].all()
    grad = comps[1].amp_prior.hyper_grad(out[3][1], out[4][1], out[5])
    h = 1e-5
    for name, idx in (("mu", 0), ("log_s", 1)):
        base = getattr(comps[1].amp_prior, name)
        setattr(comps[1].amp_prior, name, base + h)
        up = target()
        setattr(comps[1].amp_prior, name, base - h)
        down = target()
        setattr(comps[1].amp_prior, name, base)
        assert grad[idx] == pytest.approx((up - down) / (2 * h), rel=1e-6, abs=1e-6)


# --- regions and mixing ------------------------------------------------------------------------

def test_region_map_gives_parent_ring_index_and_masks() -> None:
    region = ifd.region_map(NSIDE_BETA, NSIDE)
    assert region.dtype == np.int64 and region.shape == (NPIX,)
    nest_parent = hp.ring2nest(NSIDE, np.arange(NPIX)) // (NSIDE // NSIDE_BETA)**2
    np.testing.assert_array_equal(region, hp.nest2ring(NSIDE_BETA, nest_parent))
    keep = np.ones(NPIX, dtype=bool)
    keep[:10] = False
    masked = ifd.region_map(NSIDE_BETA, NSIDE, keep)
    assert np.all(masked[:10] == -1) and np.all(masked[10:] == region[10:])
    np.testing.assert_array_equal(ifd.region_map(0, NSIDE), 0)
    with pytest.raises(ValueError):
        ifd.region_map(8, NSIDE)


def test_field_nside_and_needs_region_solver() -> None:
    comps = _components()
    assert ifd.field_nside(comps) == 0 and not ifd.needs_region_solver(comps)
    comps[1].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=1.0)
    assert ifd.needs_region_solver(comps)
    comps[1].amp_prior = None
    comps[1].set_beta_map(np.full(12, -3.0), nside=1)
    assert ifd.field_nside(comps) == 1 and ifd.needs_region_solver(comps)
    comps[2].set_beta_map(np.full(48, 1.5), nside=2)
    with pytest.raises(ValueError):
        ifd.field_nside(comps)


def test_mixing_regions_uses_the_field_where_present() -> None:
    comps = _components()
    beta_map = -3.0 + 0.1 * np.arange(12)
    comps[1].set_beta_map(beta_map, nside=1)
    M_reg, dM_reg = ifd.mixing_regions(FREQS, comps, 12, with_derivative=True)
    assert M_reg.shape == dM_reg.shape == (12, FREQS.size, 3)
    for r in range(12):
        np.testing.assert_array_equal(M_reg[r, :, 0], comps[0].get_sed(FREQS))
        np.testing.assert_allclose(M_reg[r, :, 1], (FREQS / 30.0)**beta_map[r], rtol=1e-14)
        np.testing.assert_array_equal(M_reg[r, :, 2], comps[2].get_sed(FREQS))
        np.testing.assert_allclose(dM_reg[r, :, 1], M_reg[r, :, 1] * np.log(FREQS / 30.0))
        np.testing.assert_array_equal(dM_reg[r, :, [0, 2]], 0.0)
    with pytest.raises(ValueError):
        ifd.mixing_regions(FREQS, comps, 48)


def test_augment_stream_and_mixing_shapes() -> None:
    maps_sky, inv_var = _stream()
    comps = _components()
    comps[2].amp_prior = ifd.AmpPrior("gaussian", mean=2.0, rms=4.0)
    sky_aug, w_aug = ifd.augment_stream(maps_sky, inv_var, comps)
    assert sky_aug.shape == w_aug.shape == (FREQS.size + 3, NPIX)
    np.testing.assert_array_equal(sky_aug[:FREQS.size], maps_sky)
    np.testing.assert_array_equal(w_aug[FREQS.size:FREQS.size + 2], 0.0)
    np.testing.assert_array_equal(sky_aug[-1], 2.0)
    np.testing.assert_allclose(w_aug[-1], 1.0 / 16.0)
    M_aug = ifd.augment_mixing(np.ones((2, FREQS.size, 3)))
    assert M_aug.shape == (2, FREQS.size + 3, 3) and M_aug.flags.c_contiguous
    np.testing.assert_array_equal(M_aug[:, FREQS.size:], np.broadcast_to(np.eye(3), (2, 3, 3)))


# --- the amplitude draw ------------------------------------------------------------------------

def test_solve_perpix_with_priors_matches_numpy() -> None:
    """The draw is A'^-1 (b' + fluctuation) with the region's mixing matrix, the Gaussian and
    hier priors as pseudo-bands and a flat cmb; reproduced in NumPy with the same random
    numbers."""
    maps_sky, inv_var = _stream(seed=5)
    comps = _components()
    beta_map = -3.0 + 0.1 * np.arange(12)
    comps[1].set_beta_map(beta_map, nside=NSIDE_BETA)
    g = np.exp(np.random.default_rng(6).normal(size=NPIX))
    comps[1].amp_prior = ifd.AmpPrior("hier", template=g, mu=1.1, log_s=-1.2)
    comps[2].amp_prior = ifd.AmpPrior("gaussian", mean=3.0, rms=0.5)
    draw = solve_perpix_with_priors(FREQS, maps_sky, inv_var, comps, NSIDE,
                                    np.random.default_rng(7), 1)
    rand = np.random.default_rng(7).standard_normal((NPIX, FREQS.size + 3))
    region = ifd.region_map(NSIDE_BETA, NSIDE)
    M = np.stack([np.broadcast_to(comps[0].get_sed(FREQS), (NPIX, FREQS.size)),
                  comps[1].get_sed(FREQS, beta=beta_map[region]),
                  np.broadcast_to(comps[2].get_sed(FREQS), (NPIX, FREQS.size))], axis=2)
    A = np.einsum("pbi,bp,pbj->pij", M, inv_var, M)
    b = np.einsum("pbi,bp->pi", M, inv_var * maps_sky + np.sqrt(inv_var) * rand[:, :FREQS.size].T)
    inv_T = np.stack([np.zeros(NPIX), 1.0 / (np.exp(-1.2) * g)**2, np.full(NPIX, 4.0)], axis=1)
    m = np.stack([np.zeros(NPIX), 1.1 * g, np.full(NPIX, 3.0)], axis=1)
    A += np.einsum("pi,ij->pij", inv_T, np.eye(3))
    b += inv_T * m + np.sqrt(inv_T) * rand[:, FREQS.size:]
    expected = np.linalg.solve(A, b[..., None])[..., 0].T
    np.testing.assert_allclose(draw, expected, rtol=1e-10, atol=1e-12)


def test_solve_compsep_perpix_without_fields_or_priors_is_todays_draw() -> None:
    """No field and flat priors: the plain `solve_perpix` path, bit for bit, with the random
    stream consumed as before."""
    maps_sky, inv_var = _stream(seed=8)
    comps = _components()
    data = DetectorMap(maps_sky[0], inv_var[0], nu=FREQS[0], fwhm=60.0, nside=NSIDE,
                       double_precision=True)
    # One band on COMM_SELF: the solver sees a single-band stream.
    solved = solve_compsep_perpix(MPI.COMM_SELF, data, deepcopy(comps), double_precision=True,
                                  nthreads=1, rng=np.random.default_rng(9))
    M = np.stack([c.get_sed(np.array([FREQS[0]])) for c in comps], axis=1)
    rand = np.random.default_rng(9).standard_normal((NPIX, 1))
    comp_maps = cpp_compsep.solve_perpix(maps_sky[:1], inv_var[:1], M, rand, 1)
    for icomp, comp in enumerate(solved):
        expected = curvedsky.map2alm_healpix(comp_maps[icomp][None].astype(np.float64), niter=3,
                                             spin=0, lmax=comp.lmax)
        np.testing.assert_array_equal(comp.alms, expected.astype(np.complex128))


def test_solve_compsep_perpix_routes_fields_and_priors_through_the_region_draw() -> None:
    maps_sky, inv_var = _stream(seed=10)
    comps = _components()
    comps[1].set_beta_map(-3.0 + 0.1 * np.arange(12), nside=NSIDE_BETA)
    comps[2].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=2.0)
    data = DetectorMap(maps_sky[0], inv_var[0], nu=FREQS[0], fwhm=60.0, nside=NSIDE,
                       double_precision=True)
    solved = solve_compsep_perpix(MPI.COMM_SELF, data, deepcopy(comps), double_precision=True,
                                  nthreads=1, rng=np.random.default_rng(11))
    comp_maps = solve_perpix_with_priors(np.array([FREQS[0]]), maps_sky[:1], inv_var[:1], comps,
                                         NSIDE, np.random.default_rng(11), 1)
    for icomp, comp in enumerate(solved):
        expected = curvedsky.map2alm_healpix(comp_maps[icomp][None].astype(np.float64), niter=3,
                                             spin=0, lmax=comp.lmax)
        np.testing.assert_array_equal(comp.alms, expected.astype(np.complex128))


# --- RegionIndexLikelihood: the joint hier target ---------------------------------------------

import importlib
import sys
from pathlib import Path

from commander4.math_utils import gmrf_field

_PROTO_DIR = Path(__file__).resolve().parents[1] / "nextgen" / "prototypes"


@pytest.fixture(scope="module")
def proto():
    """The prototype package imported as `indexfield` from nextgen/prototypes."""
    if not (_PROTO_DIR / "indexfield" / "posterior.py").exists():
        pytest.skip("prototype nextgen/prototypes/indexfield not present")
    if str(_PROTO_DIR) not in sys.path:
        sys.path.insert(0, str(_PROTO_DIR))
    return Bunch(field=importlib.import_module("indexfield.field"),
                 posterior=importlib.import_module("indexfield.posterior"),
                 model=importlib.import_module("indexfield.model"),
                 scalar=importlib.import_module("indexfield.scalar_coverage"))


NU_PROTO = np.array([23.0, 30.0, 44.0, 70.0, 100.0, 143.0, 353.0])
BETA_D, T_D = 1.6, 20.0
MU0, S_MU = -3.0, 0.2


def _proto_like_components(g):
    """Commander4 components with the prototype's mixing columns (cmb, sync, dust at 30 and 545
    GHz pivots) and amplitude priors: cmb ~ N(0, 70^2), sync and dust hier with template g."""
    comps = [_component(CMB, "cmb", nu_ref=1e-3),   # cmb_to_rj(1e-3 GHz) = 1 to 3e-11
             _component(Synchrotron, "sync", beta=-3.0, nu_ref=30.0),
             _component(ThermalDust, "dust", beta=BETA_D, T=T_D, nu_ref=545.0)]
    comps[0].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=70.0)
    comps[1].amp_prior = ifd.AmpPrior("hier", template=g[:, 1])
    comps[2].amp_prior = ifd.AmpPrior("hier", template=g[:, 2])
    return comps


def _proto_dataset(proto, seed, nside_beta):
    """(PixelData, g, truth) of a template-like sky with white noise, as in the E016 toy."""
    rng = np.random.default_rng(seed)
    region = ifd.region_map(nside_beta, NSIDE) if nside_beta else np.zeros(NPIX, dtype=np.int64)
    n_beta = 12 * nside_beta**2 if nside_beta else 1
    beta = MU0 + 0.15 * rng.standard_normal(n_beta)
    g = np.stack([np.ones(NPIX), 10**(1.0 + 0.5 * rng.random(NPIX)),
                  10**(1.5 + 0.5 * rng.random(NPIX))], axis=1)
    a = np.stack([70.0 * rng.standard_normal(NPIX), g[:, 1] * 10**(0.1 * rng.standard_normal(NPIX)),
                  g[:, 2] * 10**(0.1 * rng.standard_normal(NPIX))], axis=1)
    M, _ = proto.model.mixing(NU_PROTO, beta[region], beta_d=BETA_D, T_d=T_D)
    rms = 5.0 * (0.5 + rng.random((NPIX, NU_PROTO.size)))
    d = np.einsum("pbc,pc->pb", M, a) + rms * rng.standard_normal((NPIX, NU_PROTO.size))
    data = proto.model.PixelData(nu=NU_PROTO, d=d, inv_var=1.0 / rms**2, region=region,
                                 n_beta=n_beta, beta_d=BETA_D, T_d=T_D)
    return data, g


def _c4_field_target(proto, data, g, keep=None, centred=False):
    """RegionIndexLikelihood for the sync field (nside_beta 1) with hier sync and dust."""
    comps = _proto_like_components(g)
    field = gmrf_field.IndexFieldPrior(
        gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(NSIDE_BETA)),
        gmrf_field.HyperPrior(mu0=MU0, s_mu=S_MU), centred=centred)
    layout = ifd.StateLayout([ifd.IndexSpec("sync", field, MU0, S_MU)], ["sync", "dust"])
    like = ifd.RegionIndexLikelihood.from_arrays(NU_PROTO, data.d.T.copy(), data.inv_var.T.copy(),
                                                 comps, layout, NSIDE, keep_mask=keep)
    return like, layout, field


def _to_proto_order(x, K):
    """Commander4 packing (block, mu, ls, u, mu_s, ls_s, mu_d, ls_d) -> prototype
    (block, mu, ls, u, mu_s, mu_d, ls_s, ls_d)."""
    return np.concatenate([x[:K + 3], x[[K + 3, K + 5, K + 4, K + 6]]])


def test_state_layout_packs_and_unpacks() -> None:
    field = gmrf_field.IndexFieldPrior(
        gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(1)),
        gmrf_field.HyperPrior(mu0=-3.0, s_mu=0.3))
    layout = ifd.StateLayout([ifd.IndexSpec("sync", field, -3.0, 0.3),
                              ifd.IndexSpec("dust", None, 1.5, 0.1)], ["sync"])
    assert layout.dim == field.dim + 1 + 2 and layout.nreg == 12 and layout.nside_beta == 1
    rng = np.random.default_rng(0)
    blocks = {"sync": rng.normal(size=field.dim), "dust": np.array([1.52])}
    x = layout.pack(blocks, {"sync": (0.9, -1.1)})
    betas = layout.betas(x)
    np.testing.assert_allclose(betas["sync"], field.beta(blocks["sync"]))
    np.testing.assert_array_equal(betas["dust"], 1.52)
    assert layout.hypers(x) == {"sync": (0.9, -1.1)}
    with pytest.raises(ValueError):
        ifd.StateLayout([ifd.IndexSpec("sync", field, -3.0, 0.3),
                         ifd.IndexSpec("sync", None, -3.0, 0.3)], [])
    field2 = gmrf_field.IndexFieldPrior(
        gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(2)),
        gmrf_field.HyperPrior(mu0=1.5, s_mu=0.3))
    with pytest.raises(ValueError):
        ifd.StateLayout([ifd.IndexSpec("sync", field, -3.0, 0.3),
                         ifd.IndexSpec("dust", field2, 1.5, 0.3)], [])


@pytest.mark.parametrize("centred", [False, True])
def test_region_index_likelihood_matches_prototype_hier_posterior(proto, centred) -> None:
    """On identical arrays the Commander4 target equals the prototype HierNonCentredPosterior /
    HierPosterior up to a state-independent constant (the fixed cmb prior's trailing terms, which
    the prototype keeps), and the gradients agree to 1e-8."""
    data, g = _proto_dataset(proto, seed=21, nside_beta=NSIDE_BETA)
    like, layout, field = _c4_field_target(proto, data, g, centred=centred)
    pfield = proto.field.MaternField(proto.field.GraphLaplacianBasis(NSIDE_BETA))
    phyper = proto.posterior.HyperPrior(mu0=MU0, s_mu=S_MU)
    pcls = proto.posterior.HierPosterior if centred else proto.posterior.HierNonCentredPosterior
    ppost = pcls(pfield, phyper, data, g)
    rng = np.random.default_rng(22)
    K = field.K
    lps, plps = [], []
    for _ in range(3):
        xp = ppost.random_init(rng)
        x = np.concatenate([xp[:K + 3], xp[[K + 3, K + 5, K + 4, K + 6]]])
        lp, grad = like.logp_and_grad(x)
        plp, pgrad = ppost.logp_and_grad(xp)
        assert np.isfinite(lp) and np.isfinite(plp)
        lps.append(lp)
        plps.append(plp)
        np.testing.assert_allclose(_to_proto_order(grad, K), pgrad, rtol=1e-8, atol=1e-8)
    diffs = np.diff(lps) - np.diff(plps)
    np.testing.assert_allclose(diffs, 0.0, atol=1e-7 * max(1.0, abs(lps[0])))
    # The constant is the cmb prior's trailing term the prototype keeps: 1/2 N log(1/70^2).
    assert lps[0] - plps[0] == pytest.approx(-0.5 * NPIX * np.log(1.0 / 70.0**2), rel=1e-6)


def test_region_index_likelihood_gradient_matches_finite_differences(proto) -> None:
    """4th-order central differences in every block of x (field, mu, log sigma, u, mu_c, log
    s_c)."""
    data, g = _proto_dataset(proto, seed=23, nside_beta=NSIDE_BETA)
    like, layout, field = _c4_field_target(proto, data, g)
    rng = np.random.default_rng(24)
    x = np.concatenate([rng.normal(size=field.K), [MU0 + 0.05, np.log(0.2), 0.3],
                        [1.05, -1.3, 0.95, -1.1]])
    lp, grad = like.logp_and_grad(x)
    assert np.isfinite(lp)
    h = 1e-4
    fd = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h
        f = [like.logp_and_grad(x + k * e)[0] for k in (-2, -1, 1, 2)]
        fd[i] = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * h)
    np.testing.assert_allclose(grad, fd, rtol=1e-6, atol=1e-4)


def test_region_index_likelihood_masked_pixels_contribute_nothing(proto) -> None:
    """With a keep mask the target equals the unmasked target on the data set restricted to the
    kept pixels (the trailing terms run over the kept pixels only)."""
    data, g = _proto_dataset(proto, seed=25, nside_beta=NSIDE_BETA)
    keep = np.random.default_rng(26).random(NPIX) > 0.3
    like, layout, field = _c4_field_target(proto, data, g, keep=keep)
    sub = proto.model.PixelData(nu=NU_PROTO, d=data.d[keep], inv_var=data.inv_var[keep],
                                region=data.region[keep], n_beta=data.n_beta, beta_d=BETA_D,
                                T_d=T_D)
    comps = _proto_like_components(g[keep])
    like_sub = ifd.RegionIndexLikelihood(MPI.COMM_SELF, None, None, comps, layout, None, 1)
    like_sub.is_solver_root = True
    like_sub.region = data.region[keep]
    like_sub.streams = [(NU_PROTO, sub.d.T.copy(), sub.inv_var.T.copy())]
    rng = np.random.default_rng(27)
    x = np.concatenate([rng.normal(size=field.K), [MU0, np.log(0.2), 0.0], [1.0, -1.2, 1.0, -1.2]])
    lp, grad = like.logp_and_grad(x)
    lp_sub, grad_sub = like_sub.logp_and_grad(x)
    assert lp == pytest.approx(lp_sub, rel=1e-12)
    np.testing.assert_allclose(grad, grad_sub, rtol=1e-10, atol=1e-10)


def test_region_index_likelihood_rejects_guard_and_non_finite_states(proto) -> None:
    data, g = _proto_dataset(proto, seed=28, nside_beta=NSIDE_BETA)
    like, layout, field = _c4_field_target(proto, data, g)
    x = np.concatenate([np.zeros(field.K), [MU0, np.log(0.2), 0.0], [1.0, -1.2, 1.0, -1.2]])
    assert np.isfinite(like.logp_and_grad(x)[0])
    bad = x.copy()
    bad[layout.hier_slices["dust"]] = (1.0, ifd.HIER_LOGS_MIN - 1e-6)
    lp, grad = like.logp_and_grad(bad)
    assert lp == -np.inf and np.all(grad == 0.0)
    bad = x.copy()
    bad[0] = np.nan
    assert like.logp_and_grad(bad)[0] == -np.inf
    assert like.n_rejected == 2


def test_scalar_index_target_matches_prototype_scalar_hier_posterior(proto) -> None:
    """nside_beta 0: one index with prior N(mu0, sd0^2) plus the hier hyperparameters equals the
    E017 ScalarHierPosterior up to the cmb constant."""
    data, g = _proto_dataset(proto, seed=29, nside_beta=0)
    comps = _proto_like_components(g)
    layout = ifd.StateLayout([ifd.IndexSpec("sync", None, MU0, S_MU)], ["sync", "dust"])
    like = ifd.RegionIndexLikelihood.from_arrays(NU_PROTO, data.d.T.copy(), data.inv_var.T.copy(),
                                                 comps, layout, NSIDE)
    ppost = proto.scalar.ScalarHierPosterior(data, g)
    rng = np.random.default_rng(30)
    lps, plps = [], []
    for _ in range(3):
        xp = ppost.random_init(rng)                     # (beta, mu_s, mu_d, ls_s, ls_d)
        x = xp[[0, 1, 3, 2, 4]]                         # (beta, mu_s, ls_s, mu_d, ls_d)
        lp, grad = like.logp_and_grad(x)
        plp, pgrad = ppost.logp_and_grad(xp)
        lps.append(lp)
        plps.append(plp)
        np.testing.assert_allclose(grad[[0, 1, 3, 2, 4]], pgrad, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(np.diff(lps), np.diff(plps), rtol=1e-9, atol=1e-7)


# --- independent covariance-form check (from the 2026-10-27 stats review) ----------------------

def _covariance_form_target(comps, d, w, region, layout, g_s, g_d):
    """log p(x | d) written from scratch: per pixel the Gaussian marginal
    log N(d_p; M m, N + M T M^T) over the observed bands, plus the hyperpriors and the field
    prior, with the pixel set {region >= 0, at least one observed band} (proper priors)."""
    spec = layout.specs[0]
    fld = spec.field

    def sync_sed(nu, b):
        return (nu / 30.0)**b

    def dust_sed(nu, b=BETA_D, T=T_D):
        h_over_k = 0.0479924466
        x, x0 = h_over_k * nu / T, h_over_k * 545.0 / T
        return (nu / 545.0)**(b + 1) * np.expm1(x0) / np.expm1(x)

    def logp(x):
        if fld is None:
            beta = np.full(1, x[0])
            lp_field = -0.5 * ((x[0] - spec.prior_mean) / spec.prior_sd)**2
        else:
            K = fld.K
            z, mu, ls, u = x[:K], x[K], x[K + 1], x[K + 2]
            sigma = np.exp(ls)
            L = np.log(np.pi / fld.basis.h)
            s = 1.0 / (1.0 + np.exp(-u))
            rho = np.exp(np.log(fld.basis.h) + L * s)
            spectrum = (8.0 / rho**2 + fld.basis.lam)**(-1.0)
            wgt = spectrum / np.sqrt(np.sum(spectrum**2) / fld.n_beta)
            beta = mu + sigma * fld.basis.V @ (wgt * z)
            lp_field = -0.5 * z @ z - 0.5 * ((mu - MU0) / S_MU)**2 - 0.5 * (sigma / 0.3)**2 + ls \
                + np.log(L) + np.log(s) + np.log(1 - s)
        mu_s, ls_s, mu_d, ls_d = x[-4:]
        if ls_s < -8 or ls_d < -8:
            return -np.inf
        lp = 0.0
        for p in range(NPIX):
            r = region[p]
            obs = w[:, p] > 0
            if r < 0 or obs.sum() == 0:
                continue
            M = np.stack([comps[0].get_sed(NU_PROTO[obs]), sync_sed(NU_PROTO[obs], beta[r]),
                          dust_sed(NU_PROTO[obs])], axis=1)
            m = np.array([0.0, mu_s * g_s[p], mu_d * g_d[p]])
            T = np.diag([70.0**2, (np.exp(ls_s) * g_s[p])**2, (np.exp(ls_d) * g_d[p])**2])
            C = np.diag(1.0 / w[obs, p]) + M @ T @ M.T
            rsd = d[obs, p] - M @ m
            lp += -0.5 * rsd @ np.linalg.solve(C, rsd) - 0.5 * np.linalg.slogdet(C)[1]
        lp += -0.5 * ((mu_s - 1) / 0.5)**2 - 0.5 * ((mu_d - 1) / 0.5)**2
        lp += -0.5 * (ls_s - np.log(0.3))**2 - 0.5 * (ls_d - np.log(0.3))**2
        return lp + lp_field

    return logp


@pytest.mark.parametrize("scalar,mask", [(False, False), (False, True), (True, False)])
def test_region_index_likelihood_matches_covariance_form_marginal(scalar, mask) -> None:
    """Code minus an independently written covariance-form target is constant over random states
    (to 1e-9 relative on values of order 1e5) and the code gradient matches 4th-order finite
    differences of that target, with pixels observed in 0, 1 and 2 bands and under a mask."""
    rng = np.random.default_rng(123)
    g_s = 10**(1 + 0.5 * rng.random(NPIX))
    g_d = 10**(1.5 + 0.5 * rng.random(NPIX))
    comps = _proto_like_components(np.stack([np.ones(NPIX), g_s, g_d], axis=1))
    d = rng.normal(size=(NU_PROTO.size, NPIX)) * 30 + 50
    w = 1.0 / (5.0 * (0.5 + rng.random((NU_PROTO.size, NPIX))))**2
    w[:, 0] = 0.0
    w[1:, 1] = 0.0
    w[2:, 2] = 0.0
    keep = rng.random(NPIX) > 0.3 if mask else None
    if scalar:
        field = None
    else:
        field = gmrf_field.IndexFieldPrior(
            gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(NSIDE_BETA)),
            gmrf_field.HyperPrior(mu0=MU0, s_mu=S_MU))
    layout = ifd.StateLayout([ifd.IndexSpec("sync", field, MU0, S_MU)], ["sync", "dust"])
    like = ifd.RegionIndexLikelihood.from_arrays(NU_PROTO, d, w, comps, layout, NSIDE,
                                                 keep_mask=keep)
    region = ifd.region_map(0 if scalar else NSIDE_BETA, NSIDE, keep)
    mine = _covariance_form_target(comps, d, w, region, layout, g_s, g_d)
    rng2 = np.random.default_rng(7)
    xs = []
    for _ in range(3):
        core = [MU0 + 0.1 * rng2.normal()] if scalar else np.concatenate(
            [rng2.normal(size=field.K), [MU0 + 0.1 * rng2.normal(), np.log(0.2) + 0.2 * rng2.normal(),
                                         rng2.normal()]])
        xs.append(np.concatenate([core, [1 + 0.2 * rng2.normal(), np.log(0.3) + 0.3 * rng2.normal(),
                                         1 + 0.2 * rng2.normal(), np.log(0.3) + 0.3 * rng2.normal()]]))
    diff = np.array([like.logp_and_grad(x)[0] - mine(x) for x in xs])
    assert np.all(np.isfinite(diff))
    np.testing.assert_allclose(diff - diff[0], 0.0, atol=1e-9 * abs(diff[0]) + 1e-6)
    x = xs[0]
    grad = like.logp_and_grad(x)[1]
    h = 1e-4
    fd = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h
        f = [mine(x + k * e) for k in (-2, -1, 1, 2)]
        fd[i] = (f[0] - 8 * f[1] + 8 * f[2] - f[3]) / (12 * h)
    np.testing.assert_allclose(grad, fd, rtol=1e-6, atol=1e-4)


def test_region_index_likelihood_rejects_singular_and_overflowing_states() -> None:
    """A state at which A' is singular in a used pixel (two flat components with collinear columns
    at some beta) or at which the SED overflows has zero density; the pixel is never dropped."""
    rng = np.random.default_rng(31)
    d = rng.normal(size=(2, NPIX))
    w = np.ones((2, NPIX))
    comps = [_component(Synchrotron, "sync", beta=-3.0, nu_ref=30.0),
             _component(Synchrotron, "sync2", beta=-2.0, nu_ref=30.0)]   # both flat
    layout = ifd.StateLayout([ifd.IndexSpec("sync", None, -3.0, 1.0)], [])
    like = ifd.RegionIndexLikelihood.from_arrays(np.array([30.0, 44.0]), d, w, comps, layout, NSIDE)
    assert np.isfinite(like.logp_and_grad(np.array([-3.0]))[0])
    lp, grad = like.logp_and_grad(np.array([-2.0]))   # collinear with sync2: A singular
    assert lp == -np.inf and np.all(grad == 0.0) and like.nskipped_last == NPIX
    lp, _ = like.logp_and_grad(np.array([5000.0]))    # (44/30)^5000 overflows
    assert lp == -np.inf
    assert like.n_rejected == 2
    # The restore of hyperparameters after an evaluation: the components keep their own values.
    g = np.ones(NPIX)
    comps3 = _proto_like_components(np.stack([g, g, g], axis=1))
    layout3 = ifd.StateLayout([ifd.IndexSpec("sync", None, -3.0, 0.2)], ["sync", "dust"])
    like3 = ifd.RegionIndexLikelihood.from_arrays(NU_PROTO, rng.normal(size=(NU_PROTO.size, NPIX)),
                                                  np.ones((NU_PROTO.size, NPIX)), comps3, layout3,
                                                  NSIDE)
    comps3[1].amp_prior.mu, comps3[1].amp_prior.log_s = 1.3, -0.7
    like3.logp_and_grad(np.array([-3.0, 0.8, -1.5, 1.1, -1.1]))
    assert (comps3[1].amp_prior.mu, comps3[1].amp_prior.log_s) == (1.3, -0.7)


# --- hier prior with a mu map (template multiplier per region, E025 / D023) -------------------

def test_hier_mu_map_pseudo_band_and_hyperprior() -> None:
    rng = np.random.default_rng(5)
    g = np.exp(rng.normal(size=NPIX))
    prior = ifd.AmpPrior("hier", template=g, mu=1.2, log_s=-1.0, mu_nside=1)
    assert prior.n_mu == 12 and prior.hyper_dim == 13
    np.testing.assert_array_equal(prior.mu, np.full(12, 1.2))
    mu = 1.0 + 0.3 * rng.standard_normal(12)
    prior.mu = mu
    m, inv_T = prior.pseudo_band(NPIX)
    region = ifd.region_map(1, NSIDE)
    np.testing.assert_allclose(m, mu[region] * g)
    np.testing.assert_allclose(inv_T, 1.0 / (np.exp(-1.0) * g)**2)
    expected = -0.5 * np.sum(((mu - 1.0) / 0.5)**2) - 0.5 * ((-1.0 - np.log(0.3)) / 1.0)**2
    assert prior.log_hyperprior() == pytest.approx(expected)
    grad = prior.log_hyperprior_grad()
    assert grad.shape == (13,)
    h = 1e-6
    for k in range(13):
        if k < 12:
            prior.mu[k] += h
            up = prior.log_hyperprior()
            prior.mu[k] -= 2 * h
            down = prior.log_hyperprior()
            prior.mu[k] += h
        else:
            prior.log_s += h
            up = prior.log_hyperprior()
            prior.log_s -= 2 * h
            down = prior.log_hyperprior()
            prior.log_s += h
        assert grad[k] == pytest.approx((up - down) / (2 * h), rel=1e-5, abs=1e-8)
    with pytest.raises(ValueError):
        ifd.AmpPrior("hier", template=g, mu=np.ones(5), mu_nside=1)
    scalar = ifd.AmpPrior("hier", template=g, mu=0.9)
    assert scalar.n_mu == 1 and scalar.hyper_dim == 2 and isinstance(scalar.mu, float)


def test_hier_mu_map_hyper_grad_matches_finite_differences() -> None:
    """With a mu map the per-region dF/dmu_r and dF/dlog s from the kernel extras equal central
    differences of the augmented target."""
    maps_sky, inv_var = _stream(seed=6)
    comps = _components()
    rng = np.random.default_rng(7)
    g = np.exp(rng.normal(size=NPIX))
    comps[1].amp_prior = ifd.AmpPrior("hier", template=g, mu=1.0 + 0.2 * rng.standard_normal(12),
                                      log_s=-0.8, mu_nside=1)
    comps[2].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=50.0)
    comps[0].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=70.0)
    region = ifd.region_map(0, NSIDE)
    M_aug, dM_aug = ifd.augment_mixing(*ifd.mixing_regions(FREQS, comps, 1, True))

    def target(extras=False):
        sky_aug, w_aug = ifd.augment_stream(maps_sky, inv_var, comps)
        out = cpp_compsep.marginal_loglike_regions(sky_aug, w_aug, region, M_aug, dM_aug, 1,
                                                   want_extras=extras)
        used = np.ones(NPIX, dtype=bool) if not extras else out[5]
        value = float(out[0].sum()) + comps[1].amp_prior.trailing_terms(used)
        return (value, out) if extras else value

    value, out = target(extras=True)
    assert out[5].all()
    grad = comps[1].amp_prior.hyper_grad(out[3][1], out[4][1], out[5])
    assert grad.shape == (13,)
    h = 1e-5
    prior = comps[1].amp_prior
    for k in range(13):
        if k < 12:
            prior.mu[k] += h
            up = target()
            prior.mu[k] -= 2 * h
            down = target()
            prior.mu[k] += h
        else:
            prior.log_s += h
            up = target()
            prior.log_s -= 2 * h
            down = target()
            prior.log_s += h
        assert grad[k] == pytest.approx((up - down) / (2 * h), rel=1e-6, abs=1e-6)


def test_state_layout_hier_dims_round_trip() -> None:
    specs = [ifd.IndexSpec("sync", None, -3.0, 0.2)]
    layout = ifd.StateLayout(specs, ["sync", "dust"], {"sync": 13})
    assert layout.dim == 1 + 13 + 2 and layout.hier_dims == {"sync": 13, "dust": 2}
    mu = np.linspace(0.8, 1.2, 12)
    x = layout.pack({"sync": np.array([-2.9])}, {"sync": (mu, -1.1), "dust": (1.05, -1.3)})
    hyp = layout.hypers(x)
    np.testing.assert_array_equal(hyp["sync"][0], mu)
    assert hyp["sync"][1] == -1.1 and hyp["dust"] == (1.05, -1.3)
    assert isinstance(hyp["dust"][0], float)
