"""Correctness tests for the NUTS sampler on targets with known moments, and an integration test
with the per-region marginal-likelihood kernel (spatially varying spectral indices)."""
import numpy as np
import pytest

from commander4.backend import compsep as cpp_compsep
from commander4.math_utils.nuts import NUTSSampler, _warmup_windows


def _run(logp_and_grad, theta0, n_warmup, n_samples, seed):
    sampler = NUTSSampler(len(theta0), np.random.default_rng(seed), n_warmup=n_warmup)
    theta = np.asarray(theta0, dtype=float)
    draws, infos = [], []
    for _ in range(n_warmup + n_samples):
        theta, info = sampler.transition(theta, logp_and_grad)
        if not info.adapting:
            draws.append(theta)
            infos.append(info)
    return np.array(draws), infos, sampler


def _ess(x):
    """Crude effective sample size from the initial positive autocorrelation sequence."""
    x = x - x.mean()
    n = len(x)
    acf = np.correlate(x, x, mode="full")[n - 1:] / (x.var() * n)
    tau = 1.0
    for k in range(1, n // 2, 2):
        pair = acf[k] + acf[k + 1]
        if pair < 0:
            break
        tau += 2 * pair
    return n / tau


def test_warmup_windows_follow_stan_schedule():
    assert _warmup_windows(1000) == [100, 150, 250, 450, 950]
    assert _warmup_windows(10) == []
    ends = _warmup_windows(100)
    assert ends and ends[-1] <= 90


def test_correlated_gaussian_moments():
    cov = np.array([[1.0, 0.95 * 10.0], [0.95 * 10.0, 100.0]])
    prec = np.linalg.inv(cov)
    mean = np.array([3.0, -20.0])

    def f(theta):
        d = theta - mean
        return -0.5 * d @ prec @ d, -prec @ d

    draws, infos, _ = _run(f, [0.0, 0.0], 1000, 4000, seed=1)
    assert not any(i.divergent for i in infos)
    for j in range(2):
        se = np.sqrt(cov[j, j] / _ess(draws[:, j]))
        assert abs(draws[:, j].mean() - mean[j]) < 4 * se
    emp = np.cov(draws.T)
    np.testing.assert_allclose(np.diag(emp), np.diag(cov), rtol=0.15)
    assert abs(emp[0, 1] / np.sqrt(emp[0, 0] * emp[1, 1]) - 0.95) < 0.02


def test_mass_matrix_adapts_to_widely_varying_scales():
    scales = np.logspace(-3, 3, 20)

    def f(theta):
        return -0.5 * np.sum((theta / scales)**2), -theta / scales**2

    draws, infos, sampler = _run(f, np.zeros(20), 1000, 2000, seed=2)
    ratio = draws.std(axis=0) / scales
    assert np.all((ratio > 0.8) & (ratio < 1.25))
    # The adapted metric follows the target variances, so trajectories stay short. (For variances
    # far below 1e-3, Stan's regularization towards a small isotropic metric dominates by design.)
    big = scales >= 1e-2
    np.testing.assert_allclose(np.sqrt(sampler.inv_mass[big]) / scales[big], 1.0, rtol=0.5)
    assert np.mean([i.tree_depth for i in infos]) < 5


def test_non_gaussian_log_gamma():
    """theta = ln x with x ~ Gamma(k, 1): p(theta) ~ exp(k theta - e^theta), E[x] = Var[x] = k."""
    k = 3.0

    def f(theta):
        return float(k * theta[0] - np.exp(theta[0])), np.array([k - np.exp(theta[0])])

    draws, _, _ = _run(f, [0.0], 500, 6000, seed=3)
    x = np.exp(draws[:, 0])
    se = np.sqrt(k / _ess(x))
    assert abs(x.mean() - k) < 4 * se
    assert abs(x.var() / k - 1.0) < 0.15


def test_non_finite_density_is_rejected():
    """A hard boundary (zero density outside) is respected: no draw ever lands outside it."""
    def f(theta):
        if theta[0] <= 0.0:
            return -np.inf, np.zeros(1)
        return float(-theta[0]), np.array([-1.0])  # Exponential(1) on theta > 0

    draws, _, _ = _run(f, [1.0], 300, 3000, seed=4)
    assert np.all(draws > 0.0)
    assert abs(draws.mean() - 1.0) < 0.15


def test_regional_spectral_indices_match_exact_per_region_posteriors():
    """NUTS over the indices of many regions jointly, using the region kernel's analytic
    gradient, reproduces each region's exact posterior (here the regions are independent, so
    each one's exact posterior follows from a 1D grid)."""
    rng = np.random.default_rng(5)
    nu = np.array([30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0])
    nreg, npix = 12, 600
    region = np.arange(npix, dtype=np.int64) % nreg
    beta_true = -3.0 + 0.1 * rng.normal(size=nreg)
    amps = np.stack([rng.normal(0, 20, npix), 10**rng.uniform(0.5, 1.5, npix)])
    M_true = np.stack([np.ones_like(nu)[None].repeat(npix, 0),
                       (nu[None] / 30.0)**beta_true[region][:, None]], axis=2)  # (npix, nb, 2)
    d = np.einsum("pbc,cp->bp", M_true, amps) + rng.normal(size=(nu.size, npix))
    w = np.ones_like(d)
    prior_mean, prior_sd = -3.0, 0.3

    def mixing(beta):
        M = np.empty((nreg, nu.size, 2))
        M[:, :, 0] = 1.0
        M[:, :, 1] = (nu[None] / 30.0)**beta[:, None]
        dM = np.zeros_like(M)
        dM[:, :, 1] = M[:, :, 1] * np.log(nu / 30.0)[None]
        return M, dM

    def f(beta):
        ll, grad, _ = cpp_compsep.marginal_loglike_regions(d, w, region, *mixing(beta), 1)
        lp = ll.sum() - 0.5 * np.sum(((beta - prior_mean) / prior_sd)**2)
        return float(lp), grad[:, 1] - (beta - prior_mean) / prior_sd**2

    draws, infos, _ = _run(f, np.full(nreg, prior_mean), 500, 2000, seed=6)
    assert np.mean([i.divergent for i in infos]) < 0.01

    for r in range(nreg):
        grid = np.linspace(draws[:, r].mean() - 0.5, draws[:, r].mean() + 0.5, 2001)
        logp = []
        for b in grid:
            beta = np.full(nreg, -3.0)
            beta[r] = b
            ll = cpp_compsep.marginal_loglike_regions(d, w, region, *mixing(beta), 1,
                                                      True, False)[0][r]
            logp.append(ll - 0.5 * ((b - prior_mean) / prior_sd)**2)
        p = np.exp(np.array(logp) - max(logp))
        p /= p.sum()
        mean = np.sum(p * grid)
        sd = np.sqrt(np.sum(p * (grid - mean)**2))
        se = draws[:, r].std() / np.sqrt(_ess(draws[:, r]))
        assert abs(draws[:, r].mean() - mean) < 4 * se + 1e-3 * sd
        assert abs(draws[:, r].std() / sd - 1.0) < 0.2


def test_window_floor_is_relative_to_the_metric_in_force():
    """With a scaled starting metric (as `init_inv_mass: fd_hessian` gives), a coordinate whose
    posterior variance is far below Stan's absolute floor (1e-3 * 5/(n+5)) keeps its small inverse
    mass after a window update; from the unit metric the update is Stan's (floor 1e-3 * 5/(n+5)
    with n the window length, the last window here holding 50 draws)."""
    rng = np.random.default_rng(5)
    scales = np.array([1e-4, 1e-4, 1.0, 1.0])

    def logp_and_grad(x):
        return -0.5 * float(np.sum((x / scales) ** 2)), -x / scales**2

    sampler = NUTSSampler(4, rng, n_warmup=200)
    sampler.inv_mass = scales**2
    x = np.zeros(4)
    for _ in range(200):
        x, _ = sampler.transition(x, logp_and_grad)
    assert np.all(sampler.inv_mass[:2] < 1e-6)
    np.testing.assert_allclose(np.sqrt(sampler.inv_mass) / scales, 1.0, rtol=0.5)
    unit = NUTSSampler(4, np.random.default_rng(6), n_warmup=200)
    x = np.zeros(4)
    for _ in range(200):
        x, _ = unit.transition(x, logp_and_grad)
    assert np.all(unit.inv_mass[:2] >= 1e-3 * 5.0 / (50 + 5.0) * 0.99)   # Stan's floor stands


def test_fixed_metric_matches_exact_target_and_is_never_updated():
    """``adapt_mass=False``: the inverse mass given at construction survives the whole warm-up
    (no window updates) and the sampler still reproduces an exact target."""
    scales = np.array([1e-3, 1.0, 30.0])
    prec = np.diag(1.0 / scales**2)
    mean = np.array([0.5, -2.0, 10.0])

    def f(theta):
        d = theta - mean
        return -0.5 * d @ prec @ d, -prec @ d

    rng = np.random.default_rng(7)
    sampler = NUTSSampler(3, rng, n_warmup=300, inv_mass=scales**2, adapt_mass=False)
    theta = mean.copy()
    draws, grads = [], []
    for it in range(300 + 3000):
        theta, info = sampler.transition(theta, f)
        if it >= 300:
            draws.append(theta.copy())
            grads.append(info.n_leapfrog)
    draws = np.array(draws)
    np.testing.assert_array_equal(sampler.inv_mass, scales**2)
    assert sampler._window_ends == []
    for j in range(3):
        se = scales[j] / np.sqrt(_ess(draws[:, j]))
        assert abs(draws[:, j].mean() - mean[j]) < 4 * se
    np.testing.assert_allclose(draws.std(axis=0), scales, rtol=0.15)
    assert np.mean(grads) < 16   # a matched metric needs short trajectories
