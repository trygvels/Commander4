"""Independent toy for the legacy resolve_chisq move (stats review 2026-11-09).

Toy: chi2_min(beta) = beta^2, so the profile posterior is N(0, 1) (variance 1).
The excess chi2 at freshly drawn amplitudes is X ~ chi2_k.
Legacy move from (beta, X): beta' = beta + sig*N(0,1), X' ~ chi2_k,
accept with min(1, exp(-(beta'^2 - beta^2 + X' - X)/2)), keep (beta', X') on accept.
Composite: redraw X ~ chi2_k, then n legacy moves.
"""
from math import lgamma, log, exp, pi, sqrt
import numpy as np
from scipy.signal import lfilter
from scipy.special import gammaln, ndtr, betainc


# ---------- analytic predictions (reviewer's own derivation) ----------
def c_k(k):
    """c_k = E[exp(-Z/2); Z>0], Z = X'-X with X, X' iid chi2_k. Closed form 2^-(k+1) C(k, k/2)."""
    return exp(-(k + 1) * log(2.0) + lgamma(k + 1.0) - 2.0 * lgamma(k / 2.0 + 1.0))


def T_small(k):
    """Small-step temperature of the composite chain, n=1: stationary law is profile^(1/T)."""
    return 0.5 + 1.0 / (4.0 * c_k(k))


def T_lead(k):
    return sqrt(pi * k / 8.0)


def record_sums(n, ngrid=200001, lim=12.0):
    """N_acc = sum_{j=1..n} 1/(j+1); N_phi = sum_{j=1..n} E[phi(min of j iid N(0,1))]."""
    t = np.linspace(-lim, lim, ngrid)
    ph = np.exp(-0.5 * t * t) / sqrt(2 * pi)
    S = ndtr(-t)
    # sum_j j S^(j-1) for j = 1..n, closed form of the finite arithmetico-geometric sum
    j = np.arange(1, n + 1)
    tot = np.zeros_like(t)
    Sp = np.ones_like(t)
    for jj in j:
        tot += jj * Sp
        Sp = Sp * S
    N_phi = np.trapezoid(ph * ph * tot, t)
    N_acc = float(np.sum(1.0 / (j + 1.0)))
    return N_acc, N_phi


def var_scaled_n(n):
    """Leading-order Var(beta)/sqrt(k) for n moves per redraw, small steps."""
    N_acc, N_phi = record_sums(n)
    return N_acc / (2.0 * sqrt(2.0) * N_phi)


# ---------- exact effective acceptance for n = 1 (FFT convolution, no Gaussian approx) ----------
def alpha_table(k, h=None):
    """Return (z, alpha(D=-z)) with alpha(D) = E[min(1, exp(-(D + X' - X)/2))], X, X' iid chi2_k."""
    sd = sqrt(2.0 * k)
    if h is None:
        h = min(0.05, sd / 200.0)
    lo = max(0.0, k - 14.0 * sd - 40.0)
    hi = k + 14.0 * sd + 80.0
    N = int((hi - lo) / h) + 1
    x = lo + h * np.arange(N)
    with np.errstate(divide="ignore"):
        logf = (k / 2.0 - 1.0) * np.log(x) - x / 2.0 - (k / 2.0) * log(2.0) - gammaln(k / 2.0)
    f = np.exp(logf)
    f[~np.isfinite(f)] = 0.0
    f /= f.sum() * h
    # f_Z(z) = int f(x) f(x+z) dx  (autocorrelation), z_j = (j-(N-1)) h
    M = 1 << int(np.ceil(np.log2(2 * N)))
    F = np.fft.rfft(f, M)
    ac = np.fft.irfft(F * np.conj(F), M) * h
    fz = np.concatenate([ac[-(N - 1):], ac[:N]])
    fz[fz < 0] = 0.0
    z = h * (np.arange(2 * N - 1) - (N - 1))
    fz /= fz.sum() * h
    # G1(z_j) = P(Z < z_j)  (midpoint convention: half of own cell)
    G1 = np.cumsum(fz) * h - 0.5 * fz * h
    # G2(z_j) = int_{z_j}^inf fZ(z) exp(-(z - z_j)/2) dz, backward recursion
    a = exp(-h / 2.0)
    w_same = (1.0 - a) / (h / 2.0) / 2.0 * h  # weight of fz over a half-open cell, approx
    # trapezoid with exponential weight: contribution of cell [z_j, z_{j+1}]
    # int_0^h (f_j + (f_{j+1}-f_j) s/h) exp(-s/2) ds
    I0 = 2.0 * (1.0 - a)
    I1 = (4.0 * (1.0 - a) - 2.0 * h * a) / h
    cell = fz[:-1] * (I0 - I1) + fz[1:] * I1
    b = np.concatenate([cell, [0.0]])[::-1]
    G2 = lfilter([1.0], [1.0, -a], b)[::-1]
    # lfilter gives y_m = b_m + a y_{m-1}; reversed index => G2_j = cell_j + a G2_{j+1}
    return z, G1 + G2, G2


def stationary_grid(k, sig, nb=1201, Lfac=7.0, tab=None):
    """Exact (discretized) stationary law of beta for the n=1 composite chain. Returns var, kurt, acc."""
    if tab is None:
        tab = alpha_table(k)
    z, alp, _ = tab
    T = T_small(k)
    L = Lfac * sqrt(T)
    b = np.linspace(-L, L, nb)
    hb = b[1] - b[0]
    D = b[None, :] ** 2 - b[:, None] ** 2          # D_ij for move i -> j
    A = np.interp(-D, z, alp, left=0.0, right=1.0)
    q = np.exp(-0.5 * ((b[None, :] - b[:, None]) / sig) ** 2) / (sig * sqrt(2 * pi)) * hb
    K = q * A
    accrow = K.sum(axis=1)                 # acceptance incl. the same-cell proposal mass
    np.fill_diagonal(K, 0.0)
    rowacc = K.sum(axis=1)
    Kfull = K.copy()
    Kfull[np.diag_indices(nb)] = 1.0 - rowacc
    Mx = Kfull.T - np.eye(nb)
    Mx[-1, :] = 1.0
    rhs = np.zeros(nb)
    rhs[-1] = 1.0
    p = np.linalg.solve(Mx, rhs)
    m2 = np.sum(p * b * b)
    m4 = np.sum(p * b ** 4)
    return m2, m4 / m2 ** 2, float(np.sum(p * accrow)), (b, p / hb)


# ---------- simulation ----------
def simulate(k, n, sig, nch, nit, burn, rng, redraw=1.0, init_sd=1.0, x0_scale=1.0):
    """redraw: 1.0 = X ~ chi2_k (true conditional); 0.5 = X ~ chi2_k/2 (half covariance);
    None = no redraw (legacy move alone, X carried across iterations)."""
    beta = init_sd * rng.standard_normal(nch)
    half = nch // 2
    beta[:half] = 0.0                      # half the chains start under-dispersed, half over
    beta[half:] *= 2.0
    X = x0_scale * rng.chisquare(k, nch) if k > 0 else np.zeros(nch)
    s2 = np.zeros(nch)
    s4 = np.zeros(nch)
    sX = np.zeros(nch)
    nacc = 0
    for it in range(nit):
        if redraw is not None and k > 0:
            X = redraw * rng.chisquare(k, nch)
        for _ in range(n):
            bp = beta + sig * rng.standard_normal(nch)
            Xp = rng.chisquare(k, nch) if k > 0 else X
            la = -0.5 * (bp * bp - beta * beta + Xp - X)
            acc = np.log(rng.random(nch)) < la
            beta = np.where(acc, bp, beta)
            X = np.where(acc, Xp, X)
            if it >= burn:
                nacc += int(acc.sum())
        if it >= burn:
            b2 = beta * beta
            s2 += b2
            s4 += b2 * b2
            sX += X
    m = nit - burn
    v = s2 / m
    var = v.mean()
    se = v.std(ddof=1) / sqrt(nch)
    v_lo, v_hi = v[:half].mean(), v[half:].mean()
    kurt = (s4 / m).mean() / var ** 2
    return dict(var=var, se=se, var_under=v_lo, var_over=v_hi, kurt=kurt,
                acc=nacc / (m * n * nch), meanX=sX.mean() / m)
