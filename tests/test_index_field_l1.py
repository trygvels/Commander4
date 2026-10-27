"""L1 (E018): the Commander4 index-field target sampled with NUTS against references that share
no code with it, from the E012 prototype's reference.py (plain NumPy loops, own basis
construction): (a) a scalar index against a 1-D grid, (b) an nside-1 field with fixed amplitude
priors against the MH-within-Gibbs centred reference chain.

These are chain-length tests (a few minutes); they run only with C4_SLOW_TESTS=1 in the
environment, as in E018:

    C4_SLOW_TESTS=1 python -m pytest -q -p no:cacheprovider tests/test_index_field_l1.py
"""
import importlib
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from pixell.bunch import Bunch

from commander4.compsep import index_field as ifd
from commander4.math_utils import gmrf_field
from commander4.math_utils.nuts import NUTSSampler
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust

_PROTO_DIR = Path(__file__).resolve().parents[1] / "nextgen" / "prototypes"
pytestmark = pytest.mark.skipif(not os.environ.get("C4_SLOW_TESTS"),
                                reason="chain-length L1 tests; set C4_SLOW_TESTS=1")
NSIDE = 4
NPIX = 12 * NSIDE**2
NU = np.array([23.0, 30.0, 44.0, 70.0, 100.0, 143.0, 353.0])
RMS = np.array([70.0, 60.0, 60.0, 70.0, 24.0, 16.0, 80.0])
MU0, S_MU = -3.0, 0.2


@pytest.fixture(scope="module")
def ref():
    if not (_PROTO_DIR / "indexfield" / "reference.py").exists():
        pytest.skip("prototype nextgen/prototypes/indexfield not present")
    if str(_PROTO_DIR) not in sys.path:
        sys.path.insert(0, str(_PROTO_DIR))
    return importlib.import_module("indexfield.reference")


def _component(cls, name, **sed):
    p = Bunch(polarization="I", shortname=name, spatially_varying_MM=False, Cl_prior_amplitude=None,
              lmax=2, units="uK_RJ", **sed)
    object.__setattr__(p, "_name", name)
    return cls(p, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True,
               eval_pol="I", comp_name=name)


def _dataset(seed, nside_beta, noise_scale=1.0, beta_sd=0.1, sync_boost=1.0):
    """Template sky with a sync index per nside_beta pixel (scatter beta_sd around MU0), white
    noise (RMS times noise_scale), and the fixed Gaussian amplitude prior (m = 1.03 g,
    T = (0.23 g)^2 for sync and dust, cmb N(0, 70^2))."""
    rng = np.random.default_rng(seed)
    lat = np.abs(90.0 - np.degrees(np.arctan2(1.0, 1.0)) * 0 - np.degrees(
        __import__("healpy").pix2ang(NSIDE, np.arange(NPIX))[0]))
    t = np.exp(-lat / 15.0)
    g = np.stack([np.zeros(NPIX), sync_boost * 10**(0.3 + 1.0 * t), 10**(0.5 + 0.8 * t)], axis=1)
    m = 1.03 * g
    sd = 0.23 * g
    m[:, 0], sd[:, 0] = 0.0, 70.0
    region = ifd.region_map(nside_beta, NSIDE) if nside_beta else np.zeros(NPIX, dtype=np.int64)
    n_beta = 12 * nside_beta**2 if nside_beta else 1
    beta = MU0 + beta_sd * rng.standard_normal(n_beta)
    a = m + sd * rng.standard_normal(m.shape)
    comps = [_component(CMB, "cmb", nu_ref=1e-3), _component(Synchrotron, "sync", beta=MU0, nu_ref=30.0),
             _component(ThermalDust, "dust", beta=1.6, T=20.0, nu_ref=545.0)]
    M = np.stack([np.broadcast_to(comps[0].get_sed(NU), (NPIX, NU.size)),
                  comps[1].get_sed(NU, beta=beta[region]),
                  np.broadcast_to(comps[2].get_sed(NU), (NPIX, NU.size))], axis=2)
    rms = noise_scale * RMS
    d = np.einsum("pbc,pc->pb", M, a) + rms[None, :] * rng.standard_normal((NPIX, NU.size))
    inv_var = np.broadcast_to(1.0 / rms**2, d.shape).copy()
    comps[0].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=70.0)
    comps[1].amp_prior = ifd.AmpPrior("gaussian", mean=m[:, 1], rms=sd[:, 1])
    comps[2].amp_prior = ifd.AmpPrior("gaussian", mean=m[:, 2], rms=sd[:, 2])
    return comps, d, inv_var, region, n_beta, m, 1.0 / sd**2, beta


def _run_nuts(like, x0, rng, n_warmup, n_samples):
    sampler = NUTSSampler(x0.size, rng, n_warmup=n_warmup, max_depth=8)
    x = x0.copy()
    out = np.empty((n_samples, x0.size))
    n_div = 0
    for it in range(n_warmup + n_samples):
        x, info = sampler.transition(x, like.logp_and_grad)
        if it >= n_warmup:
            out[it - n_warmup] = x
            n_div += int(info.divergent)
    return out, n_div


def _ess(x):
    """Geyer initial positive sequence ESS of a 1-D chain."""
    x = np.asarray(x, dtype=float)
    n = x.size
    x = x - x.mean()
    acf = np.correlate(x, x, "full")[n - 1:] / (np.arange(n, 0, -1) * x.var() + 1e-300)
    s = 0.0
    for k in range(0, n - 1, 2):
        pair = acf[k] + (acf[k + 1] if k + 1 < n else 0.0)
        if pair <= 0 and k > 0:
            break
        s += pair
    return float(n / max(2 * s - 1.0, 1.0))


def test_scalar_index_nuts_matches_reference_grid(ref) -> None:
    comps, d, inv_var, region, n_beta, m, inv_T, _ = _dataset(seed=41, nside_beta=0)
    layout = ifd.StateLayout([ifd.IndexSpec("sync", None, MU0, S_MU)], [])
    like = ifd.RegionIndexLikelihood.from_arrays(NU, d.T.copy(), inv_var.T.copy(), comps, layout,
                                                 NSIDE)
    data = ref.from_arrays(NU, d, inv_var, region, 1, beta_d=1.6, T_d=20.0, m=m, inv_T=inv_T)
    grid = np.linspace(-3.6, -2.4, 2401)
    lp = ref.grid_1d(data, "gauss", grid, theta_fixed=(MU0, S_MU))
    p = np.exp(lp - lp.max())
    p /= p.sum()
    mean = float(np.sum(p * grid))
    sd = float(np.sqrt(np.sum(p * (grid - mean)**2)))
    assert sd > 5 * (grid[1] - grid[0])
    chain, n_div = _run_nuts(like, np.array([MU0]), np.random.default_rng(42), 300, 2000)
    ess = _ess(chain[:, 0])
    assert n_div == 0 and ess > 200
    assert abs(chain[:, 0].mean() - mean) < 4 * sd / np.sqrt(ess)
    assert 0.9 < chain[:, 0].std() / sd < 1.1


def test_field_nuts_matches_centred_gibbs_reference(ref) -> None:
    """nside_beta 1 (K = 11) sync field with fixed Gaussian amplitude priors: Commander4 NUTS
    (non-centred) against the prototype's MH-within-Gibbs reference in centred coordinates.

    The synchrotron is made 30 times brighter than the template so that the field scatter sigma
    is identified (per-pixel index error 0.2 against a scatter of 0.15); at the template
    brightness log sigma has a long tail towards small sigma (a funnel) that the centred
    reference chain does not explore (E018, 2026-10-27: reference mean -1.52 drifting towards the
    non-centred NUTS value -1.83 as it runs longer), so there it is not a usable reference."""
    comps, d, inv_var, region, n_beta, m, inv_T, beta_true = _dataset(seed=43, nside_beta=1,
                                                                      beta_sd=0.15, sync_boost=30.0)
    field = gmrf_field.IndexFieldPrior(
        gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(1)),
        gmrf_field.HyperPrior(mu0=MU0, s_mu=S_MU, sigma_scale=0.3))
    layout = ifd.StateLayout([ifd.IndexSpec("sync", field, MU0, S_MU)], [])
    like = ifd.RegionIndexLikelihood.from_arrays(NU, d.T.copy(), inv_var.T.copy(), comps, layout,
                                                 NSIDE)
    x0 = field.pack(np.zeros(field.K), MU0, np.log(0.1), 0.0)
    chain, n_div = _run_nuts(like, x0, np.random.default_rng(44), 500, 4000)
    assert n_div <= 40
    betas = np.array([field.beta(x) for x in chain])                       # (n, 12)
    K = field.K
    mu, ls, u = chain[:, K], chain[:, K + 1], chain[:, K + 2]
    log_rho = np.log(field.transform.rho(u))

    data = ref.from_arrays(NU, d, inv_var, region, n_beta, beta_d=1.6, T_d=20.0, m=m, inv_T=inv_T)
    fspec = ref.field_spec(1, "graph_laplacian")
    hyper = ref.default_hyper(fspec, mu0=MU0, s_mu=S_MU, sigma_scale=0.3, rho_prior="loguniform")
    out = ref.centred_gibbs_reference(data, fspec, hyper, "gauss", np.random.default_rng(45),
                                      n_iter=20000, n_burn=5000, n_mh_sub=5)
    for name, ours, theirs in (("beta", betas, out["beta"]), ("mu", mu, out["mu"]),
                               ("log_sigma", ls, out["log_sigma"]),
                               ("log_rho", log_rho, out["log_rho"])):
        ours = np.atleast_2d(ours.T).T
        theirs = np.atleast_2d(np.asarray(theirs).T).T
        for j in range(ours.shape[1]):
            e1, e2 = _ess(ours[:, j]), _ess(theirs[:, j])
            se = np.sqrt(ours[:, j].var() / e1 + theirs[:, j].var() / e2)
            diff = ours[:, j].mean() - theirs[:, j].mean()
            assert abs(diff) < 4 * se, f"{name}[{j}]: diff {diff:.4f}, se {se:.4f}"
            ratio = ours[:, j].std() / theirs[:, j].std()
            assert 0.8 < ratio < 1.25, f"{name}[{j}]: sd ratio {ratio:.3f}"
    # The field is recovered where the data constrain it (sanity, single realization).
    z = (betas.mean(0) - beta_true) / betas.std(0)
    assert np.all(np.abs(z) < 4.5)
