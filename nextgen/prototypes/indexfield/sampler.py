"""NUTS driver, Geyer ESS and E-BFMI for the index-field posteriors.

run_nuts wraps commander4.math_utils.nuts.NUTSSampler (Hoffman and Gelman 2014, Algorithm 6, with
dual averaging and Stan-style windowed diagonal mass adaptation). The sampler adapts during its
first n_warmup transitions; run_nuts performs n_warmup + n_samples transitions and keeps the last
n_samples states.

Effective sample size per coordinate uses Geyer's initial positive sequence estimator (Geyer 1992,
Statistical Science 7, 473): with autocovariances gamma_k of the chain and the pair sums
Gamma_m = gamma_2m + gamma_2m+1, the integrated autocorrelation time is

    tau = -1 + 2 sum_{m=0}^{M} Gamma_m / gamma_0,

where M is the largest index such that Gamma_0, ..., Gamma_M are all positive, and ESS = n / tau.
For antithetic chains tau < 1 and the ESS exceeds n; it is not capped.

E-BFMI (Betancourt 2016, arXiv:1604.00695) is sum_n (E_n+1 - E_n)^2 / sum_n (E_n - mean E)^2 over
the sampling phase, with E_n = -log pi(x_n) + 1/2 r_n^T M^-1 r_n the Hamiltonian at the start of
transition n and r_n the momentum the sampler draws for it. NUTSSampler does not return r_n, so
run_nuts hands it a thin wrapper around the random generator that records the last standard-normal
draw; the sampler forms r_n = draw / sqrt(inv_mass), so the kinetic energy is 1/2 |draw|^2.
Values below about 0.3 indicate that the momentum resampling cannot explore the energy
distribution, the funnel symptom of hierarchical models.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from commander4.math_utils.nuts import NUTSSampler


@dataclass
class Result:
    """Output of run_nuts (sampling phase only unless the name says warmup)."""

    samples: np.ndarray          # (n_samples, dim)
    logp: np.ndarray             # (n_samples,)
    n_divergent: int
    step_size: float
    ebfmi: float
    ess: np.ndarray              # (dim,)
    accept_stat: float           # mean acceptance statistic
    n_leapfrog: int              # total gradient evaluations
    tree_depth_mean: float
    tree_depth_max: int
    n_divergent_warmup: int
    n_leapfrog_warmup: int
    inv_mass: np.ndarray         # (dim,) adapted diagonal inverse mass matrix
    energy: np.ndarray           # (n_samples,) Hamiltonian at the start of each transition
    n_loglike_errors: int        # posterior evaluations that returned zero density, if tracked
    n_loglike_errors_sampling: int = 0   # the same, counted in the sampling phase only

    def summary(self) -> str:
        return (f"n={len(self.logp)} step={self.step_size:.3g} accept={self.accept_stat:.3f} "
                f"div={self.n_divergent} (warmup {self.n_divergent_warmup}) "
                f"ebfmi={self.ebfmi:.3f} ess[min/med/max]={self.ess.min():.0f}/"
                f"{np.median(self.ess):.0f}/{self.ess.max():.0f} "
                f"leapfrog={self.n_leapfrog} depth[mean/max]={self.tree_depth_mean:.1f}/"
                f"{self.tree_depth_max} loglike_errors={self.n_loglike_errors} "
                f"(sampling {self.n_loglike_errors_sampling})")


class _RecordingRNG:
    """Forwards every call to a numpy Generator and remembers the latest normal() draw."""

    def __init__(self, rng: np.random.Generator):
        self._rng = rng
        self.last_normal = None

    def normal(self, *args, **kwargs):
        x = self._rng.normal(*args, **kwargs)
        self.last_normal = x
        return x

    def __getattr__(self, name):
        return getattr(self._rng, name)


def ess_geyer(x: np.ndarray) -> np.ndarray | float:
    """Geyer initial-positive-sequence ESS of a chain (n,) or of each column of (n, dim)."""
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        return float(_ess_1d(x))
    return np.array([_ess_1d(x[:, j]) for j in range(x.shape[1])])


def _ess_1d(x: np.ndarray) -> float:
    n = x.size
    if n < 4:
        return float("nan")
    xc = x - x.mean()
    gamma0 = float(xc @ xc) / n
    if gamma0 <= 0.0 or not np.isfinite(gamma0):
        return float("nan")
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(xc, nfft)
    acov = np.fft.irfft(f * np.conj(f), nfft)[:n] / n       # biased estimator, lags 0..n-1
    m_max = n // 2
    big_gamma = acov[0:2 * m_max:2] + acov[1:2 * m_max:2]   # Gamma_m, m = 0..m_max-1
    neg = np.nonzero(big_gamma <= 0.0)[0]
    M = neg[0] if neg.size else big_gamma.size             # sum over m < M
    tau = -1.0 + 2.0 * float(big_gamma[:M].sum()) / gamma0
    if tau <= 0.0:
        tau = 1e-12
    return n / tau


def ebfmi(energy: np.ndarray) -> float:
    """Energy Bayesian fraction of missing information of a Hamiltonian series."""
    e = np.asarray(energy, dtype=float)
    var = float(np.sum((e - e.mean()) ** 2))
    if var <= 0.0:
        return float("nan")
    return float(np.sum(np.diff(e) ** 2) / var)


def run_nuts(post, rng: np.random.Generator, n_warmup: int, n_samples: int, x0=None,
             target_accept: float = 0.8, max_depth: int = 10, step_size: float | None = None,
             inv_mass: np.ndarray | None = None, progress: bool = False) -> Result:
    """Run NUTS on ``post`` (an object with dim, logp_and_grad(x), random_init(rng)).

    The sampler adapts the step size (dual averaging towards target_accept) and the diagonal mass
    matrix during the first n_warmup transitions and is then frozen. The returned samples are the
    n_samples states after warm-up. log pi of every state is taken from a per-transition cache of
    the trajectory evaluations, so no extra gradient evaluations are spent.
    """
    dim = int(post.dim)
    x = np.asarray(post.random_init(rng) if x0 is None else x0, dtype=float)
    if x.shape != (dim,):
        raise ValueError(f"x0 must have shape ({dim},)")
    rec = _RecordingRNG(rng)
    sampler = NUTSSampler(dim, rec, n_warmup=n_warmup, step_size=step_size, inv_mass=inv_mass,
                          target_accept=target_accept, max_depth=max_depth)
    cache: dict[bytes, float] = {}

    def logp_and_grad(theta):
        lp, g = post.logp_and_grad(theta)
        cache[np.asarray(theta).tobytes()] = float(lp)
        return lp, g

    logp_x, _ = post.logp_and_grad(x)
    if not np.isfinite(logp_x):
        raise ValueError("initial state has zero density")
    samples = np.empty((n_samples, dim))
    logp = np.empty(n_samples)
    energy = np.empty(n_samples)
    accept, n_leap, depth, n_div = 0.0, 0, [], 0
    n_div_w, n_leap_w = 0, 0
    n_total = n_warmup + n_samples
    errors_after_warmup = int(getattr(post, "n_loglike_errors", 0))
    for it in range(n_total):
        if it == n_warmup:
            errors_after_warmup = int(getattr(post, "n_loglike_errors", 0))
        cache.clear()
        x_new, info = sampler.transition(x, logp_and_grad)
        kinetic = 0.5 * float(np.sum(rec.last_normal ** 2))
        lp_new = cache.get(x_new.tobytes())
        if lp_new is None:
            lp_new = float(post.logp_and_grad(x_new)[0])
        if it >= n_warmup:
            k = it - n_warmup
            energy[k] = -logp_x + kinetic          # Hamiltonian at the start of this transition
            samples[k] = x_new
            logp[k] = lp_new
            accept += info.accept_stat
            n_leap += info.n_leapfrog
            depth.append(info.tree_depth)
            n_div += int(info.divergent)
        else:
            n_div_w += int(info.divergent)
            n_leap_w += info.n_leapfrog
        x, logp_x = x_new, lp_new
        if progress and (it + 1) % max(1, n_total // 10) == 0:
            print(f"    nuts {it + 1}/{n_total} step={sampler.step_size:.3g} "
                  f"depth={info.tree_depth} div={n_div_w + n_div}")
    return Result(samples=samples, logp=logp, n_divergent=n_div, step_size=float(sampler.step_size),
                  ebfmi=ebfmi(energy), ess=ess_geyer(samples),
                  accept_stat=accept / n_samples, n_leapfrog=n_leap,
                  tree_depth_mean=float(np.mean(depth)), tree_depth_max=int(np.max(depth)),
                  n_divergent_warmup=n_div_w, n_leapfrog_warmup=n_leap_w,
                  inv_mass=np.array(sampler.inv_mass, dtype=float), energy=energy,
                  n_loglike_errors=int(getattr(post, "n_loglike_errors", 0)),
                  n_loglike_errors_sampling=int(getattr(post, "n_loglike_errors", 0))
                  - errors_after_warmup)


# ------------------------------------------------------------------------------------ self-tests
class _Gaussian:
    """Correlated Gaussian test target N(m, Sigma)."""

    def __init__(self, mean, cov):
        self.mean = np.asarray(mean, dtype=float)
        self.cov = np.asarray(cov, dtype=float)
        self.prec = np.linalg.inv(self.cov)
        self.dim = self.mean.size

    def logp_and_grad(self, x):
        r = x - self.mean
        pr = self.prec @ r
        return -0.5 * float(r @ pr), -pr

    def random_init(self, rng):
        return self.mean + rng.normal(size=self.dim) * 3.0


def _selftest():
    rng = np.random.default_rng(7)
    ok = True

    def report(name, value, thresh, passed):
        nonlocal ok
        ok &= passed
        print(f"  {'PASS' if passed else 'FAIL'}  {name}: {value:.4g} (threshold {thresh:g})")

    # 1. ESS estimator on AR(1) chains with known tau = (1 + phi) / (1 - phi), and on iid noise.
    print("Geyer ESS on AR(1)")
    n = 400_000
    for phi in (0.0, 0.5, 0.9, -0.5):
        eps = rng.normal(size=n)
        x = np.empty(n)
        x[0] = eps[0]
        for t in range(1, n):
            x[t] = phi * x[t - 1] + eps[t]
        tau = (1 + phi) / (1 - phi)
        ratio = ess_geyer(x) / (n / tau)
        report(f"phi={phi}: ESS / (n/tau)", ratio, 0.1, abs(ratio - 1) < 0.1)

    # 2. NUTS on a 20-D correlated Gaussian: means and covariances within 4 Monte Carlo SE.
    print("NUTS on a 20-D correlated Gaussian")
    dim = 20
    scales = np.logspace(-1, 1, dim)
    corr = 0.9 ** np.abs(np.subtract.outer(np.arange(dim), np.arange(dim)))
    cov = corr * np.outer(scales, scales)
    mean = rng.normal(size=dim) * 2.0
    target = _Gaussian(mean, cov)
    res = run_nuts(target, rng, n_warmup=1000, n_samples=6000)
    print("        " + res.summary())
    s = res.samples
    ess = res.ess
    report("min ESS", ess.min(), 100, ess.min() > 100)
    # Means: SE_i = sd_i / sqrt(ESS_i).
    zmean = (s.mean(axis=0) - mean) / (np.sqrt(np.diag(cov)) / np.sqrt(ess))
    print(f"        mean z-scores: max |z| = {np.abs(zmean).max():.2f}")
    report("max |z| of means", np.abs(zmean).max(), 4.0, np.abs(zmean).max() < 4.0)
    # Covariances: the estimator of Sigma_ij is the mean of p_n = (x_i,n - m_i)(x_j,n - m_j)
    # (true means), whose SE is sd(p) / sqrt(ESS(p)) with ESS of the product series itself.
    r = s - mean
    zcov = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(i, dim):
            p = r[:, i] * r[:, j]
            se = p.std(ddof=1) / np.sqrt(ess_geyer(p))
            zcov[i, j] = zcov[j, i] = (p.mean() - cov[i, j]) / se
    rms = np.sqrt(np.mean(zcov[np.triu_indices(dim)] ** 2))
    print(f"        covariance z-scores: max |z| = {np.abs(zcov).max():.2f} over "
          f"{dim * (dim + 1) // 2} elements, rms {rms:.2f}")
    report("max |z| of covariances", np.abs(zcov).max(), 4.0, np.abs(zcov).max() < 4.0)
    sd_ratio = s.std(axis=0, ddof=1) / np.sqrt(np.diag(cov))
    print(f"        sd ratios: min {sd_ratio.min():.3f} max {sd_ratio.max():.3f}")
    # E-BFMI of a Gaussian target with a well adapted metric is close to 1.
    report("E-BFMI (Gaussian, expect ~1)", res.ebfmi, 0.5, 0.5 < res.ebfmi < 1.5)
    report("divergences", res.n_divergent, 0, res.n_divergent == 0)
    print("        adapted sqrt(inv_mass) / true sd: "
          f"min {np.sqrt(res.inv_mass / np.diag(cov)).min():.2f} "
          f"max {np.sqrt(res.inv_mass / np.diag(cov)).max():.2f}")
    # The cached logp equals a direct evaluation.
    direct = np.array([target.logp_and_grad(v)[0] for v in s[:200]])
    d = float(np.max(np.abs(res.logp[:200] - direct)))
    report("cached logp = direct logp", d, 1e-12, d < 1e-12)

    print("sampler.py self-test:", "ALL PASS" if ok else "FAILURES")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest() else 1)
