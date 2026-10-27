"""Metropolis-Adjusted Microcanonical Sampler (MAMS) in NumPy, for the M2.3 / Q4 comparison.

Algorithm: Robnik, Cohn-Gordon and Seljak, arXiv 2503.01707 (MAMS), on the isokinetic dynamics of
arXiv 2212.08549; extraction with equation numbers in nextgen/research/2026-10-27_mams_algorithm.md.

One transition from x (with cached log density lp and gradient of -log p, g):
  1. draw u uniformly on the unit sphere S^(d-1) (full refresh, Sec. 4);
  2. n = max(1, ceil(h rescale(L/eps))) leapfrog steps B_{eps/2} A_eps B_{eps/2} (Eq. 5), h ~ U(0,1),
     with the exact sphere kick (Eq. 28) in preconditioned coordinates y = x / sigma;
  3. accept with probability min(1, exp(-W)), W = L(x_n) - L(x_0) + sum_kicks (d-1) log[cosh(delta)
     + (e.u) sinh(delta)] (Lemma 4.2, Eqs. 6 and 8); on rejection x stays, u is discarded.
Tuning (Sec. 6): dual averaging of eps to the target acceptance (0.9 as in the paper's benchmarks),
diagonal preconditioner sigma from the sample variances of a warm-up window, trajectory length
L = sqrt(d) in sigma units then the ALBA rule L = 0.3 eps mean_i(N / ESS_i) from a warm-up window
(the code form of Eq. 12 to 14). Cost: n gradient evaluations per transition.

Interface as sampler.run_nuts: run_mams(post, rng, n_warmup, n_samples) -> Result with samples,
n_leapfrog (gradients in the sampling phase), n_leapfrog_warmup, n_divergent (W above the
divergence threshold or non-finite), step_size, ess, accept_stat, L.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from indexfield.sampler import ess_geyer

TARGET_ACCEPT = 0.9
ALBA_FACTOR = 0.3
DIVERGENCE_W = 1000.0


@dataclass
class MAMSResult:
    samples: np.ndarray
    logp: np.ndarray
    n_divergent: int
    step_size: float
    L: float
    sigma: np.ndarray
    ess: np.ndarray
    accept_stat: float
    n_leapfrog: int
    n_leapfrog_warmup: int
    n_divergent_warmup: int
    ebfmi: float = float("nan")
    tree_depth_mean: float = float("nan")
    tree_depth_max: int = 0


def rescale(mu: float) -> float:
    """s such that n = ceil(h s), h ~ U(0, 1), has E[n] = mu (BlackJAX `rescale`, footnote 5 of
    2503.01707 as derived in the research note)."""
    k = np.floor(2.0 * mu - 1.0)
    if k < 1:
        return max(mu, 1e-12) * 2.0
    return float(k + k * (mu - 0.5 * (k + 1.0)) / (k + 1.0 - mu))


def kick(u: np.ndarray, g: np.ndarray, h: float) -> tuple[np.ndarray, float]:
    """Exact isokinetic velocity update B_h for the force -g (g = grad of -log p in the
    preconditioned coordinates), Eq. 28 of 2503.01707, and its weight (d-1) log[cosh(delta) +
    (e.u) sinh(delta)] (Eq. 8), written with zeta = exp(-delta) for stability."""
    d = u.size
    gn = float(np.linalg.norm(g))
    if gn == 0.0:
        return u, 0.0
    e = -g / gn
    delta = h * gn / (d - 1)
    zeta = np.exp(-delta)
    c = float(e @ u)
    u_new = e * (1.0 - zeta) * (1.0 + zeta + c * (1.0 - zeta)) + 2.0 * zeta * u
    u_new /= np.linalg.norm(u_new)
    dk = (d - 1) * (delta - np.log(2.0) + np.log((1.0 + zeta * zeta) + c * (1.0 - zeta * zeta)))
    return u_new, float(dk)


def _logp_grad(post, x):
    lp, grad_lp = post.logp_and_grad(x)
    return float(lp), -np.asarray(grad_lp, dtype=float)      # g = grad L, L = -log p


def transition(post, x, lp, g, eps: float, L: float, sigma: np.ndarray, rng: np.random.Generator):
    """One MAMS transition. Returns (x, lp, g, accepted, W, n_steps, divergent)."""
    d = x.size
    u = rng.standard_normal(d)
    u /= np.linalg.norm(u)
    n = max(1, int(np.ceil(rng.uniform() * rescale(L / eps))))
    x0, lp0, g0 = x, lp, g
    dK = 0.0
    xc, lpc, gc = x, lp, g
    finite = True
    for _ in range(n):
        u, dk = kick(u, sigma * gc, 0.5 * eps)
        dK += dk
        xc = xc + eps * sigma * u
        lpc, gc = _logp_grad(post, xc)
        if not (np.isfinite(lpc) and np.all(np.isfinite(gc))):
            finite = False
            break
        u, dk = kick(u, sigma * gc, 0.5 * eps)
        dK += dk
    if not finite:
        return x0, lp0, g0, False, np.inf, n, True
    W = (-lpc) - (-lp0) + dK
    divergent = not np.isfinite(W) or W > DIVERGENCE_W
    if not divergent and np.log(rng.uniform()) < -W:
        return xc, lpc, gc, True, W, n, False
    return x0, lp0, g0, False, W, n, divergent


class _DualAveraging:
    """Step-size adaptation of Hoffman and Gelman (2014), section 3.2, on the MH acceptance."""

    def __init__(self, step_size: float, target: float, gamma=0.05, t0=10.0, kappa=0.75):
        self.target, self.gamma, self.t0, self.kappa = target, gamma, t0, kappa
        self.mu = np.log(10.0 * step_size)
        self.h_bar = 0.0
        self.log_eps_bar = np.log(step_size)
        self.m = 0

    def update(self, accept: float) -> float:
        self.m += 1
        w = 1.0 / (self.m + self.t0)
        self.h_bar = (1 - w) * self.h_bar + w * (self.target - accept)
        log_eps = self.mu - self.h_bar * np.sqrt(self.m) / self.gamma
        eta = self.m ** (-self.kappa)
        self.log_eps_bar = eta * log_eps + (1 - eta) * self.log_eps_bar
        return float(np.exp(log_eps))

    def final(self) -> float:
        return float(np.exp(self.log_eps_bar))


def run_mams(post, rng: np.random.Generator, n_warmup: int, n_samples: int, x0=None,
             target_accept: float = TARGET_ACCEPT, step_size: float | None = None,
             L: float | None = None, sigma: np.ndarray | None = None,
             alba_factor: float = ALBA_FACTOR, preconditioning: bool = True,
             progress: bool = False, **_ignored) -> MAMSResult:
    """Run MAMS on `post` (dim, logp_and_grad(x), random_init(rng)).

    Warm-up (n_warmup transitions, adapted state kept, samples discarded), in the order of Sec. 6:
    the first 40 percent adapt eps by dual averaging at L = sqrt(d) (sigma units) with sigma = 1;
    the next 30 percent keep adapting eps while collecting samples from which sigma (per-dimension
    sd) and the ALBA L are set at the end of the window; the last 30 percent re-adapt eps for the
    new (sigma, L). Sampling then runs at the dual-averaged step size.
    """
    d = int(post.dim)
    x = np.asarray(post.random_init(rng) if x0 is None else x0, dtype=float)
    lp, g = _logp_grad(post, x)
    if not np.isfinite(lp):
        raise ValueError("initial state has zero density")
    sig = np.ones(d) if sigma is None else np.asarray(sigma, dtype=float)
    Lc = float(np.sqrt(d)) if L is None else float(L)
    eps = 0.1 if step_size is None else float(step_size)
    n1 = int(0.4 * n_warmup)
    n2 = int(0.7 * n_warmup)
    da = _DualAveraging(eps, target_accept)
    window = []
    n_leap_w = n_div_w = 0
    for it in range(n_warmup):
        x, lp, g, acc, W, n, div = transition(post, x, lp, g, eps, Lc, sig, rng)
        n_leap_w += n
        n_div_w += int(div)
        eps = da.update(0.0 if not np.isfinite(W) else min(1.0, np.exp(-W)))
        if n1 <= it < n2:
            window.append(x.copy())
        if it == n2 - 1 and len(window) > 10:
            w = np.array(window)
            if preconditioning and sigma is None:
                sd = w.std(axis=0)
                sig = np.where(sd > 0, sd, 1.0)
            if L is None:
                ess = ess_geyer(w / sig)
                ess = np.where(np.isfinite(ess) & (ess > 0), ess, 1.0)
                Lc = float(alba_factor * eps * np.mean(w.shape[0] / ess))
                Lc = min(max(Lc, eps), 500 * eps)
            da = _DualAveraging(eps, target_accept)
        if progress and (it + 1) % max(1, n_warmup // 5) == 0:
            print(f"    mams warmup {it + 1}/{n_warmup} eps={eps:.3g} L={Lc:.3g}")
    eps = da.final() if n_warmup > 0 else eps
    samples = np.empty((n_samples, d))
    logp = np.empty(n_samples)
    acc_sum = 0.0
    n_leap = n_div = 0
    for k in range(n_samples):
        x, lp, g, acc, W, n, div = transition(post, x, lp, g, eps, Lc, sig, rng)
        samples[k] = x
        logp[k] = lp
        acc_sum += 0.0 if not np.isfinite(W) else min(1.0, np.exp(-W))
        n_leap += n
        n_div += int(div)
    return MAMSResult(samples=samples, logp=logp, n_divergent=n_div, step_size=eps, L=Lc, sigma=sig,
                      ess=ess_geyer(samples), accept_stat=acc_sum / max(1, n_samples),
                      n_leapfrog=n_leap, n_leapfrog_warmup=n_leap_w, n_divergent_warmup=n_div_w)


# ------------------------------------------------------------------------------------ self-tests
def _selftest():
    """Exact-target checks: (1) the kick weight equals the log-Jacobian form (d-1) log[cosh(delta)
    + (e.u) sinh(delta)] and preserves |u| = 1; (2) rescale gives E[n] = mu; (3) a 20-d correlated
    Gaussian: means within 4 MC SE, variances within 10 percent, with ESS per gradient against NUTS
    on the same target."""
    from indexfield.sampler import _Gaussian, run_nuts
    rng = np.random.default_rng(3)
    ok = True
    # (1) kick weight
    d = 7
    u = rng.standard_normal(d)
    u /= np.linalg.norm(u)
    g = rng.standard_normal(d) * 3.0
    u2, dk = kick(u, g, 0.3)
    gn = np.linalg.norm(g)
    e = -g / gn
    delta = 0.3 * gn / (d - 1)
    ref = (d - 1) * np.log(np.cosh(delta) + (e @ u) * np.sinh(delta))
    print(f"kick: |u'| - 1 = {abs(np.linalg.norm(u2) - 1):.1e}, weight {dk:.10f} vs {ref:.10f}")
    ok &= abs(dk - ref) < 1e-10 and abs(np.linalg.norm(u2) - 1) < 1e-12
    # (2) rescale
    for mu in (1.3, 2.0, 7.5, 40.0):
        s = rescale(mu)
        n = np.ceil(rng.uniform(size=400000) * s)
        print(f"rescale({mu}) = {s:.4f}: E[n] = {n.mean():.4f}")
        ok &= abs(n.mean() - mu) < 0.02 * mu
    # (3) Gaussian
    d = 20
    A = rng.standard_normal((d, d))
    cov = A @ A.T / d + np.diag(np.linspace(0.1, 10.0, d))
    mean = rng.standard_normal(d)
    target = _Gaussian(mean, cov)
    res = run_mams(target, rng, n_warmup=2000, n_samples=20000)
    z = (res.samples.mean(0) - mean) / np.sqrt(np.diag(cov) / res.ess)
    vr = res.samples.var(0) / np.diag(cov)
    print(f"mams: eps {res.step_size:.3f} L {res.L:.2f} accept {res.accept_stat:.2f} grads "
          f"{res.n_leapfrog} ESS min/med {res.ess.min():.0f}/{np.median(res.ess):.0f} "
          f"max|z| {np.abs(z).max():.2f} var ratio {vr.min():.3f} to {vr.max():.3f} "
          f"div {res.n_divergent}; ESS/grad min {res.ess.min() / res.n_leapfrog:.4f}")
    ok &= np.abs(z).max() < 4.0 and np.all((vr > 0.9) & (vr < 1.1)) and res.n_divergent == 0
    nres = run_nuts(target, np.random.default_rng(4), n_warmup=2000, n_samples=20000)
    print(f"nuts: grads {nres.n_leapfrog} ESS min/med {nres.ess.min():.0f}/"
          f"{np.median(nres.ess):.0f}; ESS/grad min {nres.ess.min() / nres.n_leapfrog:.4f}; "
          f"MAMS/NUTS ESS per gradient (min) = "
          f"{(res.ess.min() / res.n_leapfrog) / (nres.ess.min() / nres.n_leapfrog):.2f}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(_selftest())
