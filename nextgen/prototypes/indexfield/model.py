"""Per-pixel data model and collapsed (amplitude-marginalized) likelihood of the index field.

Mathematics: E012 v2, sections "Data model", "Amplitude treatment" and "Gradient". Interface and
shapes: SPEC.md, section model.py. Units: GHz, uK_RJ, radians.

For each data pixel p with mixing matrix M_p (B x C) evaluated at beta_{r(p)}, noise precision
N_p^-1 = diag(inv_var_p) and data d_p:

    A_p = M_p^T N_p^-1 M_p  (C x C),    b_p = M_p^T N_p^-1 d_p  (C,).

The three collapsed log-likelihoods implemented by `collapsed_loglike` are, per pixel,

    "w1"    :  1/2 b_p^T A_p^-1 b_p - 1/2 log det A_p            (flat amplitude prior, marginal)
    "w0"    :  1/2 b_p^T A_p^-1 b_p                              (profile / ridge likelihood)
    "gauss" :  1/2 b'_p^T A'_p^-1 b'_p - 1/2 log det A'_p,       (proper Gaussian prior a ~ N(m, T))
               A'_p = A_p + T_p^-1,  b'_p = b_p + T_p^-1 m_p

summed over p. Dropped constants: see the docstring of `collapsed_loglike`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Physical constants (CODATA 2018), in units that make h nu / (k T) dimensionless with nu in GHz.
H_OVER_K_GHZ = 6.62607015e-34 / 1.380649e-23 * 1e9      # h/k in K/GHz = 0.0479924...
T_CMB = 2.7255                                           # K

COMPONENTS = ("cmb", "sync", "dust")
SYNC_COL = 1                                             # column of M carrying beta_s


def cmb_to_rj(nu_ghz):
    """uK_CMB to uK_RJ conversion factor x^2 e^x / (e^x - 1)^2 with x = h nu / (k T_CMB).

    Tends to 1 for nu -> 0. Accepts a float or an ndarray of frequencies in GHz.
    """
    x = np.asarray(nu_ghz, dtype=float) * (H_OVER_K_GHZ / T_CMB)
    ex = np.exp(x)
    out = x**2 * ex / (ex - 1.0)**2
    return float(out) if out.ndim == 0 else out


def mixing(nu_ghz, beta_s, beta_d=1.6, T_d=20.0, nu0_s=30.0, nu0_d=545.0):
    """Mixing matrix M (B, C) in uK_RJ and the derivative of its synchrotron column.

    Columns, in the order (cmb, sync, dust):
        cmb  = cmb_to_rj(nu)
        sync = (nu / nu0_s)^beta_s
        dust = (nu / nu0_d)^(beta_d + 1) (exp(h nu0_d / k T_d) - 1) / (exp(h nu / k T_d) - 1)
    (a modified black body in RJ units normalized to 1 at nu0_d). Returns (M, dM_dbeta_s) with
    dM_dbeta_s = log(nu / nu0_s) sync, the derivative of the sync column, shape (B,).

    beta_s may also be an array of shape (N,): then M has shape (N, B, C) and dM_dbeta_s has
    shape (N, B). This is the vectorized form used by `collapsed_loglike` (SPEC: "stack M as
    (N, B, C) from mixing evaluated at beta[region]").
    """
    nu = np.asarray(nu_ghz, dtype=float)
    beta_s = np.asarray(beta_s, dtype=float)
    cmb = cmb_to_rj(nu)
    dust = ((nu / nu0_d)**(beta_d + 1.0) * np.expm1(H_OVER_K_GHZ * nu0_d / T_d)
            / np.expm1(H_OVER_K_GHZ * nu / T_d))
    log_ratio = np.log(nu / nu0_s)
    if beta_s.ndim == 0:
        sync = np.exp(beta_s * log_ratio)
        M = np.stack([cmb, sync, dust], axis=1)
        return M, log_ratio * sync
    sync = np.exp(beta_s[:, None] * log_ratio[None, :])                  # (N, B)
    M = np.empty(sync.shape + (3,))
    M[:, :, 0] = cmb
    M[:, :, 1] = sync
    M[:, :, 2] = dust
    return M, log_ratio[None, :] * sync


@dataclass
class PixelData:
    """Per-pixel data set (SPEC.md). All maps in RING order at nside_d.

    nu: (B,) GHz. d: (N, B) uK_RJ. inv_var: (N, B) uK_RJ^-2. region: (N,) int, index pixel of
    each data pixel in [0, n_beta). beta_d, T_d: fixed dust parameters. m, inv_T: (N, C) mean and
    diagonal precision of the Gaussian amplitude prior (mode "gauss"), or None.
    """
    nu: np.ndarray
    d: np.ndarray
    inv_var: np.ndarray
    region: np.ndarray
    n_beta: int
    beta_d: float
    T_d: float
    m: np.ndarray | None = None
    inv_T: np.ndarray | None = None


def _collapsed_core(M, dM, d, inv_var, region, n_beta, mode, m=None, inv_T=None,
                    sync_col=SYNC_COL):
    """Collapsed log-likelihood and per-index-pixel gradient for a given mixing stack.

    M: (N, B, C), dM: (N, B) derivative of column `sync_col` with respect to beta. The component
    count C is whatever M carries (the C = 2 quadrature self-test uses this entry point). See
    `collapsed_loglike` for the formulas and the dropped constants.

    Gradient (E012 "Gradient", envelope theorem): with x = A^-1 b (or A'^-1 b'),
        g_p = x_c g^T N^-1 (d - M x) - w (A^-1 M^T N^-1 g)_c,    g = dM_p[:, c] / dbeta,
    w = 1 for "w1" and "gauss", w = 0 for "w0". The first term is the derivative of the maximized
    quadratic form at fixed maximizer; the second is -1/2 tr(A^-1 dA/dbeta) with dA = e_c g^T N^-1
    M + M^T N^-1 g e_c^T.
    """
    if mode not in ("gauss", "w0", "w1"):
        raise ValueError(f"mode must be 'gauss', 'w0' or 'w1', got {mode!r}")
    N, B, C = M.shape
    NinvM = M * inv_var[:, :, None]                                    # N^-1 M        (N, B, C)
    A = np.einsum("nbc,nbk->nck", M, NinvM)                             # M^T N^-1 M    (N, C, C)
    b = np.einsum("nbc,nb->nc", NinvM, d)                               # M^T N^-1 d    (N, C)
    if mode == "gauss":
        if m is None or inv_T is None:
            raise ValueError("mode 'gauss' needs the amplitude prior (m, inv_T)")
        idx = np.arange(C)
        A[:, idx, idx] += inv_T
        b = b + inv_T * m
    try:
        L = np.linalg.cholesky(A)
    except np.linalg.LinAlgError as err:
        bad = np.flatnonzero(np.linalg.eigvalsh(A)[:, 0] <= 0.0)
        raise ValueError(f"A (M^T N^-1 M{' + T^-1' if mode == 'gauss' else ''}) is not positive "
                         f"definite in {bad.size} pixel(s), first: {bad[:5].tolist()}") from err
    # Solves through the Cholesky factor: A x = b  <=>  L y = b, L^T x = y. For C <= 3 the
    # batched general solve is as fast as a triangular one in NumPy, so use it on L and L^T.
    x = np.linalg.solve(np.swapaxes(L, 1, 2), np.linalg.solve(L, b[:, :, None]))[:, :, 0]
    quad = 0.5 * np.einsum("nc,nc->n", b, x)                            # 1/2 b^T A^-1 b
    logdet = 2.0 * np.sum(np.log(np.diagonal(L, axis1=1, axis2=2)), axis=1)
    w = 0.0 if mode == "w0" else 1.0
    logL_p = quad - 0.5 * w * logdet

    resid = d - np.einsum("nbc,nc->nb", M, x)                           # d - M x       (N, B)
    gNr = np.einsum("nb,nb->n", dM * inv_var, resid)                    # g^T N^-1 (d - M x)
    grad_p = x[:, sync_col] * gNr
    if w != 0.0:
        MNg = np.einsum("nbc,nb->nc", NinvM, dM)                        # M^T N^-1 g    (N, C)
        y = np.linalg.solve(np.swapaxes(L, 1, 2), np.linalg.solve(L, MNg[:, :, None]))[:, :, 0]
        grad_p = grad_p - w * y[:, sync_col]
    grad = np.bincount(region, weights=grad_p, minlength=n_beta)
    return float(np.sum(logL_p)), grad


def collapsed_loglike(beta, data: PixelData, mode: str):
    """Collapsed log-likelihood log L(beta) and its gradient with respect to the index field.

    beta: (Nb,) synchrotron index per index pixel; evaluated at beta[data.region] per data pixel.
    mode: "gauss", "w0" or "w1". Returns (logL: float, grad: ndarray (Nb,)) with
    grad[r] = sum_{p in r} g_p (E012 "Gradient"; pixels grouped by region with np.bincount).

    Per pixel, with A = M^T N^-1 M, b = M^T N^-1 d, A' = A + T^-1, b' = b + T^-1 m:
        "w1"   : 1/2 b^T A^-1 b   - 1/2 log det A
        "w0"   : 1/2 b^T A^-1 b
        "gauss": 1/2 b'^T A'^-1 b' - 1/2 log det A'

    Dropped constants. Relative to the per-mode formulas in SPEC.md the only term dropped is
    -1/2 d^T N^-1 d, in all modes (so "w1" and "w0" return -1/2 chi2min + 1/2 d^T N^-1 d
    - (w/2) log det A). Nothing else is dropped: every grid, reference or quadrature check must
    use this same convention. For orientation, the exact log marginal
    log int N(d; M a, N) p(a) da equals the returned value plus the beta-independent
        -1/2 d^T N^-1 d - (B/2) log 2 pi - 1/2 log det N + (C/2) log 2 pi          ("w1"),
        -1/2 d^T N^-1 d - (B/2) log 2 pi - 1/2 log det N - 1/2 m^T T^-1 m
            - 1/2 log det T                                                         ("gauss"),
    which the quadrature self-test below adds back. "w0" is a profile, not a marginal.

    Raises ValueError if any A (or A') is not positive definite (np.linalg.cholesky on the
    (N, C, C) stack), or if mode is "gauss" without (m, inv_T).
    """
    beta = np.asarray(beta, dtype=float)
    if beta.shape != (data.n_beta,):
        raise ValueError(f"beta must have shape ({data.n_beta},), got {beta.shape}")
    M, dM = mixing(data.nu, beta[data.region], beta_d=data.beta_d, T_d=data.T_d)
    return _collapsed_core(M, dM, data.d, data.inv_var, data.region, data.n_beta, mode,
                           m=data.m, inv_T=data.inv_T, sync_col=SYNC_COL)


# ----------------------------------------------------------------------------------------------
# Self-tests
# ----------------------------------------------------------------------------------------------

BANDS_GHZ_TEST = np.array([23., 30., 33., 41., 44., 61., 70., 94., 100., 143., 217., 353.])


def _toy_data(rng, N, B, Nb, beta_true=-3.0, with_prior=True):
    """Synthetic PixelData: template-like amplitudes, white noise, a Gaussian amplitude prior."""
    nu = BANDS_GHZ_TEST[:B]
    region = np.sort(rng.integers(0, Nb, size=N))
    beta = beta_true + 0.2 * rng.standard_normal(Nb)
    a = np.stack([70.0 * rng.standard_normal(N),
                  10**(0.3 + 2.0 * rng.random(N)),
                  10**(0.5 + 1.5 * rng.random(N))], axis=1)
    M, _ = mixing(nu, beta[region])
    rms = np.full((N, B), 30.0) * (0.5 + rng.random((N, B)))
    d = np.einsum("nbc,nc->nb", M, a) + rms * rng.standard_normal((N, B))
    m = inv_T = None
    if with_prior:
        m = np.stack([np.zeros(N), 10**(0.3 + rng.random(N)), 10**(0.5 + rng.random(N))], axis=1)
        inv_T = 1.0 / (10.0 * np.maximum(np.abs(m), 10.0))**2
    data = PixelData(nu=nu, d=d, inv_var=1.0 / rms**2, region=region, n_beta=Nb, beta_d=1.6,
                     T_d=20.0, m=m, inv_T=inv_T)
    return data, beta


def _check_gradient(rng, N=200, B=12, Nb=4, h=3e-5):
    """Gradient against central finite differences per index pixel, all three modes."""
    data, beta_true = _toy_data(rng, N, B, Nb)
    beta0 = beta_true + 0.1 * rng.standard_normal(Nb)       # off the optimum, large gradients
    ok_all = True
    for mode in ("gauss", "w0", "w1"):
        logL, grad = collapsed_loglike(beta0, data, mode)
        g_fd = np.empty(Nb)
        for r in range(Nb):
            e = np.zeros(Nb)
            e[r] = h
            g_fd[r] = (collapsed_loglike(beta0 + e, data, mode)[0]
                       - collapsed_loglike(beta0 - e, data, mode)[0]) / (2 * h)
        err = np.abs(g_fd - grad)
        tol = 1e-6 * np.abs(grad) + 1e-7 * abs(logL)
        ok = bool(np.all(err < tol))
        ok_all &= ok
        print(f"  gradient {mode:5s}: logL={logL:.6e} |g|max={np.abs(grad).max():.3e} "
              f"max|g_fd-g|={err.max():.3e} max rel={np.max(err / np.abs(grad)):.2e} "
              f"tol min={tol.min():.3e}  {'PASS' if ok else 'FAIL'}")
    return ok_all


def _check_mixing_derivative(h=1e-6):
    """dM/dbeta_s against finite differences and cmb_to_rj(nu -> 0) -> 1."""
    nu = BANDS_GHZ_TEST
    M, dM = mixing(nu, -3.0)
    fd = (mixing(nu, -3.0 + h)[0][:, SYNC_COL] - mixing(nu, -3.0 - h)[0][:, SYNC_COL]) / (2 * h)
    err = np.max(np.abs(fd - dM)) / np.max(np.abs(dM))     # dM = 0 exactly at nu = nu0_s
    Mv, dMv = mixing(nu, np.array([-3.0, -2.5]))
    stack_ok = np.allclose(Mv[0], M) and np.allclose(dMv[0], dM)
    ok = err < 1e-8 and abs(cmb_to_rj(1e-3) - 1.0) < 1e-6 and stack_ok and \
        abs(M[np.argmin(np.abs(nu - 30.0)), SYNC_COL] - 1.0) < 1e-12
    print(f"  mixing: max rel err dM/dbeta = {err:.2e}, cmb_to_rj(1e-3 GHz) = "
          f"{cmb_to_rj(1e-3):.8f}, cmb_to_rj(100) = {cmb_to_rj(100.0):.6f}, "
          f"stacked == scalar: {stack_ok}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def _check_not_pd(rng):
    """A pixel with zero weight gives A = 0, which must raise."""
    data, beta = _toy_data(rng, 20, 12, 2, with_prior=False)
    data.inv_var[3] = 0.0
    try:
        collapsed_loglike(beta, data, "w1")
    except ValueError as err:
        print(f"  non-PD raises ValueError: PASS ({str(err)[:60]}...)")
        return True
    print("  non-PD did not raise: FAIL")
    return False


def _check_timing(rng, N=3072, B=12, Nb=192, n_rep=20):
    """Wall time of one evaluation at the E012 size (target under 20 ms)."""
    import time
    data, beta = _toy_data(rng, N, B, Nb)
    ok = True
    for mode in ("gauss", "w0", "w1"):
        collapsed_loglike(beta, data, mode)
        t0 = time.perf_counter()
        for _ in range(n_rep):
            collapsed_loglike(beta, data, mode)
        dt = (time.perf_counter() - t0) / n_rep * 1e3
        ok &= dt < 20.0
        print(f"  timing {mode:5s} N={N} B={B} Nb={Nb}: {dt:.2f} ms per evaluation "
              f"{'PASS' if dt < 20.0 else 'FAIL (target 20 ms)'}")
    return ok


def _check_quadrature(rng, B=6, n_grid=801, half_width_sd=9.0):
    """Mode "gauss" (and "w1") against brute-force 2-D trapezoid integration over a_p, C = 2.

    A C = 2 variant (cmb, sync) built inside this test goes through `_collapsed_core`, the same
    code path as `collapsed_loglike`. The brute-force integral of
        exp(-1/2 (d - M a)^T N^-1 (d - M a)) p(a),   p(a) = N(a; m, T)  or  1,
    is evaluated on a grid centred at the posterior mean and spanning +-half_width_sd posterior
    sd per component; the trapezoid rule converges exponentially for a Gaussian integrand. The
    dropped constants listed in `collapsed_loglike` are added back before comparing.
    """
    nu = BANDS_GHZ_TEST[:B]
    beta_s = -2.9
    M3, dM3 = mixing(nu, beta_s)
    M = M3[None, :, :2].copy()                               # (1, B, 2): cmb and sync columns
    dM = dM3[None, :]
    a_true = np.array([30.0, 60.0])
    rms = 20.0 * (0.5 + rng.random(B))
    d = (M[0] @ a_true + rms * rng.standard_normal(B))[None, :]
    inv_var = (1.0 / rms**2)[None, :]
    m = np.array([[0.0, 40.0]])
    inv_T = 1.0 / np.array([[100.0, 400.0]])**2
    region = np.zeros(1, dtype=int)

    ok_all = True
    for mode in ("gauss", "w1"):
        logL, _ = _collapsed_core(M, dM, d, inv_var, region, 1, mode, m=m, inv_T=inv_T)
        # Posterior mean and covariance of a for the grid placement (only used to place the grid).
        A = M[0].T @ (M[0] * inv_var[0, :, None])
        b = M[0].T @ (inv_var[0] * d[0])
        if mode == "gauss":
            A = A + np.diag(inv_T[0])
            b = b + inv_T[0] * m[0]
        x = np.linalg.solve(A, b)
        sd = np.sqrt(np.diag(np.linalg.inv(A)))
        a0 = x[0] + half_width_sd * sd[0] * np.linspace(-1, 1, n_grid)
        a1 = x[1] + half_width_sd * sd[1] * np.linspace(-1, 1, n_grid)
        A0, A1 = np.meshgrid(a0, a1, indexing="ij")
        model = M[0][:, 0] * A0[..., None] + M[0][:, 1] * A1[..., None]           # (n, n, B)
        log_f = -0.5 * np.sum(inv_var[0] * (d[0] - model)**2, axis=-1)
        if mode == "gauss":
            log_f = log_f - 0.5 * (inv_T[0, 0] * (A0 - m[0, 0])**2
                                   + inv_T[0, 1] * (A1 - m[0, 1])**2)
        peak = log_f.max()
        wgt = np.ones(n_grid)
        wgt[0] = wgt[-1] = 0.5
        integral = np.sum(np.exp(log_f - peak) * wgt[:, None] * wgt[None, :])
        log_int = peak + np.log(integral) + np.log(a0[1] - a0[0]) + np.log(a1[1] - a1[0])
        # The integrand above omits the prior normalization (2 pi)^(-C/2) det(T)^(-1/2) and the
        # noise normalization, which are beta independent on both sides. Completing the square,
        #   log_int = -1/2 dNd [- 1/2 mTm] + 1/2 b'^T A'^-1 b' - 1/2 log det A' + (C/2) log 2 pi,
        # so the code value (which drops only -1/2 dNd relative to SPEC) is expected to equal
        #   log_int + 1/2 dNd [+ 1/2 mTm] - (C/2) log 2 pi,    C = 2.
        dNd = np.sum(inv_var[0] * d[0]**2)
        expected = log_int + 0.5 * dNd - np.log(2 * np.pi)
        if mode == "gauss":
            expected += 0.5 * np.sum(inv_T[0] * m[0]**2)
        err = abs(logL - expected)
        ok = err < 1e-8 * max(1.0, abs(logL))
        ok_all &= ok
        print(f"  quadrature {mode:5s} C=2: code={logL:.10e} brute={expected:.10e} "
              f"|diff|={err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok_all


if __name__ == "__main__":
    rng = np.random.default_rng(12345)
    print("model.py self-tests")
    results = [
        _check_mixing_derivative(),
        _check_gradient(rng, N=200, B=12, Nb=4),
        _check_not_pd(rng),
        _check_quadrature(rng),
        _check_timing(rng, N=3072, B=12, Nb=192),
    ]
    print("ALL PASS" if all(results) else "SOME FAILED")
    raise SystemExit(0 if all(results) else 1)
