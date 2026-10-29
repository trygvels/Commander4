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
    """``adapt_mass=False``: an inverse mass assigned after construction (the fd_hessian path)
    survives the whole warm-up untouched, including a caller that later edits its own array, and
    the sampler reproduces an exact non-Gaussian target even with a metric that is wrong by a
    factor 3 to 10 in every coordinate."""
    # Target: x0 = log of a Gamma(3) variable (mean digamma(3) = 0.92278, variance trigamma(3)
    # = 0.39493) times a rho 0.8 Gaussian in (x1, x2) with variances 1 and 9.
    cov = np.array([[1.0, 0.8 * 3.0], [0.8 * 3.0, 9.0]])
    prec = np.linalg.inv(cov)
    true_var = np.array([0.39493407, 1.0, 9.0])

    def f(theta):
        x, y = theta[0], theta[1:]
        logp = 3.0 * x - np.exp(x) - 0.5 * y @ prec @ y
        return logp, np.concatenate([[3.0 - np.exp(x)], -prec @ y])

    rng = np.random.default_rng(7)
    sampler = NUTSSampler(3, rng, n_warmup=400, adapt_mass=False)
    metric = np.array([0.13, 3.0, 0.9])     # 3x too small, 3x too large, 10x too small
    sampler.inv_mass = metric
    metric_at_start = metric.copy()
    theta = np.array([np.log(3.0), 0.0, 0.0])
    draws, steps = [], []
    for it in range(400 + 6000):
        theta, info = sampler.transition(theta, f)
        if it == 10:
            metric[:] = 1e6                  # a caller editing its array must not reach the kernel
        if it >= 400:
            draws.append(theta.copy())
            steps.append(info.step_size)
    draws = np.array(draws)
    np.testing.assert_array_equal(sampler.inv_mass, metric_at_start)
    assert sampler._window_ends == []
    assert len(set(steps)) == 1              # one fixed post-warm-up kernel
    means = np.array([0.92278434, 0.0, 0.0])
    for j in range(3):
        ess = _ess(draws[:, j])
        assert ess > 300
        assert abs(draws[:, j].mean() - means[j]) < 4 * np.sqrt(true_var[j] / ess)
        # Var of the sample variance is about 2 var^2 / ess for a Gaussian; 4 se.
        assert abs(draws[:, j].var() - true_var[j]) < 4 * true_var[j] * np.sqrt(2.0 / ess)
    emp = np.corrcoef(draws[:, 1], draws[:, 2])[0, 1]
    assert abs(emp - 0.8) < 0.05


def test_reset_metric_during_warmup_replaces_metric_and_keeps_target():
    """``reset_metric`` (the fd_hessian refresh): a metric that is wrong by 3 to 10x is replaced
    at transition 30 by the true variances; the step size is re-found and the adaptation
    restarted; the post-warm-up kernel then reproduces the exact target, and the call is refused
    after warm-up."""
    cov = np.array([[1.0, 0.8 * 3.0], [0.8 * 3.0, 9.0]])
    prec = np.linalg.inv(cov)
    true_var = np.array([0.39493407, 1.0, 9.0])

    def f(theta):
        x, y = theta[0], theta[1:]
        logp = 3.0 * x - np.exp(x) - 0.5 * y @ prec @ y
        return logp, np.concatenate([[3.0 - np.exp(x)], -prec @ y])

    rng = np.random.default_rng(11)
    sampler = NUTSSampler(3, rng, n_warmup=400, adapt_mass=False)
    sampler.inv_mass = np.array([0.13, 3.0, 0.9])
    theta = np.array([np.log(3.0), 0.0, 0.0])
    draws, steps = [], []
    for it in range(400 + 6000):
        if it == 30:
            old_step = sampler.step_size
            sampler.reset_metric(true_var, theta, f)
            np.testing.assert_array_equal(sampler.inv_mass, true_var)
            assert sampler._window_samples == []
            assert sampler.step_size != old_step
        theta, info = sampler.transition(theta, f)
        if it >= 400:
            draws.append(theta.copy())
            steps.append(info.step_size)
    draws = np.array(draws)
    np.testing.assert_array_equal(sampler.inv_mass, true_var)
    assert len(set(steps)) == 1
    means = np.array([0.92278434, 0.0, 0.0])
    for j in range(3):
        ess = _ess(draws[:, j])
        assert ess > 300
        assert abs(draws[:, j].mean() - means[j]) < 4 * np.sqrt(true_var[j] / ess)
        assert abs(draws[:, j].var() - true_var[j]) < 4 * true_var[j] * np.sqrt(2.0 / ess)
    assert abs(np.corrcoef(draws[:, 1], draws[:, 2])[0, 1] - 0.8) < 0.05
    with pytest.raises(ValueError, match="warm-up"):
        sampler.reset_metric(true_var, theta, f)


def test_two_refreshes_with_windows_keep_target():
    """Two ``reset_metric`` calls under ``adapt_mass=True`` (one inside a window, one on a window
    end) on a correlated Gaussian: the post-warm-up kernel is fixed and reproduces the target."""
    cov = np.array([[1.0, 0.8 * 3.0], [0.8 * 3.0, 9.0]])
    prec = np.linalg.inv(cov)

    def f(theta):
        return -0.5 * theta @ prec @ theta, -prec @ theta

    rng = np.random.default_rng(5)
    sampler = NUTSSampler(2, rng, n_warmup=300, adapt_mass=True)
    ends = sampler._window_ends
    theta = np.zeros(2)
    draws, steps = [], []
    for it in range(300 + 5000):
        if it in (ends[0] - 10, ends[1]):
            sampler.reset_metric(np.array([0.5, 20.0]), theta, f)
        theta, info = sampler.transition(theta, f)
        if it >= 300:
            draws.append(theta.copy())
            steps.append(info.step_size)
    draws = np.array(draws)
    assert len(set(steps)) == 1
    for j, var in enumerate((1.0, 9.0)):
        ess = _ess(draws[:, j])
        assert ess > 300
        assert abs(draws[:, j].mean()) < 4 * np.sqrt(var / ess)
        assert abs(draws[:, j].var() - var) < 4 * var * np.sqrt(2.0 / ess)
    assert abs(np.corrcoef(draws[:, 0], draws[:, 1])[0, 1] - 0.8) < 0.05


# ---------------------------------------------------------------- dense metric
def _gauss(cov, mean):
    prec = np.linalg.inv(cov)

    def f(theta):
        d = theta - mean
        return -0.5 * d @ prec @ d, -prec @ d
    return f


def test_dense_metric_equals_diagonal_metric_when_diagonal():
    """A dense inverse mass that is diagonal gives the same chain as the diagonal one from the
    same random stream (momentum draw, kinetic energy, leapfrog and U-turn test agree)."""
    cov = np.diag([1.0, 25.0, 0.04])
    f = _gauss(cov, np.array([1.0, -2.0, 0.5]))
    diag = np.array([0.7, 20.0, 0.05])
    chains = []
    for metric in (diag, np.diag(diag)):
        sampler = NUTSSampler(3, np.random.default_rng(11), n_warmup=50, adapt_mass=False,
                              inv_mass=metric)
        theta = np.zeros(3)
        out = []
        for _ in range(120):
            theta, _ = sampler.transition(theta, f)
            out.append(theta.copy())
        chains.append(np.array(out))
    assert chains[1].shape == chains[0].shape
    np.testing.assert_allclose(chains[1], chains[0], rtol=1e-9, atol=1e-12)


def _correlated_gaussian_runs():
    """An 8-dimensional Gaussian with correlations of 0.9 and scales spanning a factor 100,
    sampled with the exact covariance as a dense inverse mass and with the best diagonal metric
    (the marginal variances); 3000 post-warm-up draws each."""
    rng = np.random.default_rng(3)
    dim = 8
    scales = np.logspace(-1, 1, dim)
    corr = 0.9 ** np.abs(np.subtract.outer(np.arange(dim), np.arange(dim)))
    cov = corr * np.outer(scales, scales)
    mean = rng.normal(size=dim)
    f = _gauss(cov, mean)
    results = {}
    for name, metric in (("dense", cov), ("diag", np.diag(cov))):
        sampler = NUTSSampler(dim, np.random.default_rng(5), n_warmup=300, adapt_mass=False,
                              inv_mass=metric)
        theta = np.zeros(dim)
        draws, leaps = [], []
        for it in range(300 + 3000):
            theta, info = sampler.transition(theta, f)
            if it >= 300:
                assert not info.divergent
                draws.append(theta.copy())
                leaps.append(info.n_leapfrog)
        results[name] = (np.array(draws), np.array(leaps))
    return cov, corr, mean, results


def test_dense_metric_matches_the_moments_of_a_correlated_gaussian():
    cov, corr, mean, results = _correlated_gaussian_runs()
    draws, _ = results["dense"]
    dim = cov.shape[0]
    ess = np.array([_ess(draws[:, j]) for j in range(dim)])
    assert ess.min() > 500
    for j in range(dim):
        assert abs(draws[:, j].mean() - mean[j]) < 4 * np.sqrt(cov[j, j] / ess[j])
        assert abs(draws[:, j].var() - cov[j, j]) < 4 * cov[j, j] * np.sqrt(2.0 / ess[j])
    emp = np.corrcoef(draws.T)
    for i in range(dim):
        for j in range(i + 1, dim):
            # SE of a sample correlation, (1 - rho^2) / sqrt(n_eff), at the smaller ESS; 4 SE.
            se = (1.0 - corr[i, j] ** 2) / np.sqrt(min(ess[i], ess[j]))
            assert abs(emp[i, j] - corr[i, j]) < 4 * se + 0.01


def test_dense_metric_needs_fewer_leapfrogs_than_the_diagonal_one():
    _, _, _, results = _correlated_gaussian_runs()
    assert np.median(results["dense"][1]) < 0.5 * np.median(results["diag"][1])


def test_dense_metric_on_a_non_gaussian_target_with_a_reset_during_warmup():
    """A log-Gamma coordinate times a rho 0.95 Gaussian pair with scales 1 and 20, sampled with
    a dense metric that is not aligned with the coordinates, reset to a second dense metric at
    warm-up transition 150: the moments of the exact target are reproduced."""
    cov = np.array([[1.0, 0.95 * 20.0], [0.95 * 20.0, 400.0]])
    prec = np.linalg.inv(cov)
    true_mean = np.array([0.92278434, 0.0, 0.0])          # digamma(3)
    true_var = np.array([0.39493407, 1.0, 400.0])          # trigamma(3)

    def f(theta):
        x, y = theta[0], theta[1:]
        return 3.0 * x - np.exp(x) - 0.5 * y @ prec @ y, np.concatenate([[3.0 - np.exp(x)],
                                                                            -prec @ y])
    rng = np.random.default_rng(21)
    A = rng.normal(size=(3, 3))
    m0 = A @ A.T + 0.5 * np.eye(3)                          # a random SPD metric
    full = np.zeros((3, 3)); full[0, 0] = true_var[0]; full[1:, 1:] = cov
    sampler = NUTSSampler(3, np.random.default_rng(9), n_warmup=400, adapt_mass=False,
                          inv_mass=m0)
    theta = np.array([np.log(3.0), 0.0, 0.0])
    draws = []
    for it in range(400 + 6000):
        if it == 150:
            sampler.reset_metric(0.7 * full, theta, f)    # dense reset, 30 percent too tight
            assert sampler.dense_metric
        theta, info = sampler.transition(theta, f)
        if it >= 400:
            draws.append(theta.copy())
    draws = np.array(draws)
    for j in range(3):
        ess = _ess(draws[:, j])
        assert ess > 300
        assert abs(draws[:, j].mean() - true_mean[j]) < 4 * np.sqrt(true_var[j] / ess)
        assert abs(draws[:, j].var() - true_var[j]) < 4 * true_var[j] * np.sqrt(2.0 / ess)
    assert abs(np.corrcoef(draws[:, 1], draws[:, 2])[0, 1] - 0.95) < 0.03


def test_dense_metric_rejects_adaptation_and_bad_matrices():
    sampler = NUTSSampler(2, np.random.default_rng(0), n_warmup=10, adapt_mass=False)
    with pytest.raises(ValueError):
        sampler.inv_mass = np.array([[1.0, 2.0], [2.0, 1.0]])       # not positive definite
    with pytest.raises(ValueError):
        sampler.inv_mass = np.array([[1.0, 0.5], [0.0, 1.0]])       # not symmetric
    with pytest.raises(ValueError):
        NUTSSampler(2, np.random.default_rng(0), n_warmup=10, adapt_mass=True,
                    inv_mass=np.eye(2))
    sampler.inv_mass = np.array([[2.0, 0.5], [0.5, 1.0]])
    assert sampler.dense_metric
    sampler.inv_mass = np.array([1.0, 2.0])
    assert not sampler.dense_metric


def test_fd_inverse_mass_dense_recovers_the_covariance():
    """On a quadratic log density the dense finite-difference Hessian gives the covariance; an
    eigendirection of negative curvature keeps its scale (1/|lam|), curvatures below 1/cap are
    capped, and a coordinate whose stencil point is rejected (non-finite density, zero
    gradient) is decoupled at unit inverse mass instead of corrupting the matrix."""
    from commander4.compsep.index_field_sampler import fd_inverse_mass_dense
    rng = np.random.default_rng(2)
    A = rng.normal(size=(5, 5))
    Q, _ = np.linalg.qr(A)
    lam = np.array([50.0, 4.0, 0.5, 0.05, -2.0])      # 0.05 < 1/cap = 0.1; -2 is a saddle
    H = (Q * lam) @ Q.T
    f = lambda x: (-0.5 * x @ H @ x, -H @ x)
    inv, summ = fd_inverse_mass_dense(f, np.zeros(5))
    expected = (Q * np.array([1 / 50.0, 1 / 4.0, 2.0, 10.0, 0.5])) @ Q.T
    np.testing.assert_allclose(inv, expected, rtol=1e-6, atol=1e-8)
    assert summ["n_negative"] == 1 and summ["n_capped"] == 1 and summ["n_rejected"] == 0
    assert summ["asymmetry"] < 1e-6 and summ["condition"] == pytest.approx(500.0)
    sampler = NUTSSampler(5, np.random.default_rng(1), n_warmup=100, adapt_mass=False,
                          inv_mass=inv)
    assert sampler.dense_metric
    # Rejected stencil: the density is -inf for x[0] above 5e-4 (the +h point of coordinate 0).
    Hc = np.array([[4.0, 1.0], [1.0, 2.0]])

    def g(x):
        if x[0] > 5e-4:
            return -np.inf, np.zeros(2)
        return -0.5 * x @ Hc @ x, -Hc @ x
    inv2, summ2 = fd_inverse_mass_dense(g, np.zeros(2))
    assert summ2["n_rejected"] == 1
    np.testing.assert_allclose(inv2, np.diag([1.0, 0.5]), rtol=1e-6, atol=1e-9)
    # Same contract for the diagonal estimate.
    from commander4.compsep.index_field_sampler import fd_inverse_mass
    np.testing.assert_allclose(fd_inverse_mass(g, np.zeros(2)), [1.0, 0.5], rtol=1e-6)
