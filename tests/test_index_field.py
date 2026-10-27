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
    with pytest.raises(ValueError):
        prior.pseudo_band(4)


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
