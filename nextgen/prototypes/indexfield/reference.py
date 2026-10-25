"""Independent reference path for the index-field prototype (SPEC.md, section reference.py).

Everything in this module is written from E012 v2 and SPEC.md alone and shares no code with
model.py, field.py, posterior.py or sampler.py, so that an error in those modules cannot cancel
against an error here. All functions take plain arrays (see from_arrays) and never need the
PixelData class; any object exposing the attribute names listed in from_arrays works by duck typing.

Contents
  loglike_pixel_loops          plain-loop collapsed likelihood, modes gauss, w0, w1
  region_loglike_table, grid_1d, grid_2d
                               exact grids for N_beta = 1 and N_beta = 2
  region_map, field_spec, matern_spectrum, logdet_K
                               own construction of the basis V, eigenvalues and Matern spectrum
  rho_from_u, u_from_rho, log_dlogrho_du, RefHyper, default_hyper, log_hyperprior
  hyper_conditional_grid, draw_from_grid
                               exact 2-D grid draw of (log sigma, u) | y
  centred_gibbs_reference      MH-within-Gibbs reference chain in centred coordinates
  linear_gaussian_closed_form  closed form for the fixed-amplitude linear model
  ess_geyer                    effective sample size (Geyer initial positive sequence)

Dropped constants. In every likelihood mode the term -1/2 d^T N^-1 d is dropped (it does not depend
on beta) and nothing else is dropped. With A = M^T N^-1 M and b = M^T N^-1 d per data pixel,

  w1:    1/2 b^T A^-1 b - 1/2 log det A
  w0:    1/2 b^T A^-1 b
  gauss: 1/2 b'^T A'^-1 b' - 1/2 log det A',   A' = A + T^-1,   b' = b + T^-1 m,

summed over data pixels. Self-test (i) shows that the gauss mode equals the exact log marginal
log N(d; M m, M T M^T + N) up to the constant -1/2 d^T N^-1 d - 1/2 m^T T^-1 m - 1/2 log det N
- 1/2 log det T - B/2 log 2 pi, by brute-force quadrature over the amplitudes.

Units: frequencies in GHz, data and amplitudes in uK_RJ, angles (rho, h) in radians.
"""
from __future__ import annotations

import math
import sys
import time
from typing import Callable, NamedTuple

import healpy as hp
import numpy as np

try:
    import numba

    _HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    numba = None
    _HAVE_NUMBA = False

H_PLANCK = 6.62607015e-34
K_BOLTZ = 1.380649e-23
T_CMB = 2.7255
NU0_S = 30.0
NU0_D = 545.0
MODE_INT = {"gauss": 0, "w0": 1, "w1": 2}


# ----------------------------------------------------------------------------------------------
# Data container
# ----------------------------------------------------------------------------------------------
class RefData(NamedTuple):
    """Plain-array data container with the attribute names of SPEC PixelData (plus F, offset)."""

    nu: np.ndarray
    d: np.ndarray
    inv_var: np.ndarray
    region: np.ndarray
    n_beta: int
    beta_d: float
    T_d: float
    m: np.ndarray | None
    inv_T: np.ndarray | None
    F: np.ndarray | None
    offset: np.ndarray | None


def from_arrays(nu, d, inv_var, region, n_beta, beta_d=1.6, T_d=20.0, m=None, inv_T=None,
                F=None, offset=None) -> RefData:
    """Build the reference data container from plain arrays.

    nu (B,) GHz; d (N, B) uK_RJ; inv_var (N, B) uK_RJ^-2 (1/sigma^2, white noise); region (N,) int
    in [0, n_beta) giving the index pixel of each data pixel; beta_d, T_d the fixed dust SED
    parameters; m (N, C) and inv_T (N, C) the mean and diagonal precision of the Gaussian amplitude
    prior (mode "gauss" only; C = 3 in the order cmb, sync, dust). For the fixed-amplitude linear
    model of linear_gaussian_closed_form and mode "linear" of centred_gibbs_reference, F (N, B) is
    the per-pixel response d = offset + F beta_region + n and offset (N, B) the known signal.

    A SPEC PixelData object can be passed directly to every function here instead of this
    container (duck typing on the attribute names); this function exists so that no PixelData
    class is needed. Nothing is copied or validated beyond dtype and contiguity.
    """
    nu = np.ascontiguousarray(nu, dtype=np.float64)
    d = np.ascontiguousarray(d, dtype=np.float64)
    inv_var = np.ascontiguousarray(inv_var, dtype=np.float64)
    region = np.ascontiguousarray(region, dtype=np.int64)
    if d.shape != inv_var.shape or d.shape[1] != nu.shape[0] or region.shape[0] != d.shape[0]:
        raise ValueError("inconsistent shapes of nu, d, inv_var, region")
    if region.min() < 0 or region.max() >= n_beta:
        raise ValueError("region entries must lie in [0, n_beta)")
    m = None if m is None else np.ascontiguousarray(m, dtype=np.float64)
    inv_T = None if inv_T is None else np.ascontiguousarray(inv_T, dtype=np.float64)
    F = None if F is None else np.ascontiguousarray(F, dtype=np.float64)
    offset = None if offset is None else np.ascontiguousarray(offset, dtype=np.float64)
    return RefData(nu, d, inv_var, region, int(n_beta), float(beta_d), float(T_d), m, inv_T, F,
                   offset)


# ----------------------------------------------------------------------------------------------
# Plain-loop collapsed likelihood
# ----------------------------------------------------------------------------------------------
def _maybe_jit(fn):
    return numba.njit(cache=False)(fn) if _HAVE_NUMBA else fn


def _fixed_columns(nu, beta_d, T_d):
    """CMB and dust columns of the mixing matrix (shared by all pixels).

    cmb: uK_CMB -> uK_RJ factor x^2 e^x / (e^x - 1)^2 with x = h nu / (k T_cmb).
    dust: modified black body in RJ units normalized to one at nu0_d,
          (nu/nu0_d)^(beta_d + 1) (exp(h nu0_d / k T_d) - 1) / (exp(h nu / k T_d) - 1).
    """
    nu = np.asarray(nu, dtype=np.float64)
    x = H_PLANCK * nu * 1e9 / (K_BOLTZ * T_CMB)
    cmb = x * x * np.exp(x) / (np.exp(x) - 1.0) ** 2
    xd = H_PLANCK * nu * 1e9 / (K_BOLTZ * T_d)
    xd0 = H_PLANCK * NU0_D * 1e9 / (K_BOLTZ * T_d)
    dust = (nu / NU0_D) ** (beta_d + 1.0) * (np.exp(xd0) - 1.0) / (np.exp(xd) - 1.0)
    return np.ascontiguousarray(cmb), np.ascontiguousarray(dust)


def _pixel_term_py(nu, cmb, dust, dp, wp, beta, mode, mp, tp, M, A, bvec, L, z):
    """One pixel's contribution; explicit loops and a hand-written C x C Cholesky.

    Builds M (B x 3), A = M^T N^-1 M, b = M^T N^-1 d, adds T^-1 and T^-1 m in mode 0 (gauss),
    factorizes A = L L^T and returns (1/2 ||L^-1 b||^2 - w/2 log det A, True) with w = 0 for mode 1
    and w = 1 otherwise, or (0, False) if A is not positive definite.
    """
    B = nu.shape[0]
    C = 3
    for b in range(B):
        M[b, 0] = cmb[b]
        M[b, 1] = (nu[b] / NU0_S) ** beta
        M[b, 2] = dust[b]
    for i in range(C):
        bvec[i] = 0.0
        for j in range(C):
            A[i, j] = 0.0
    for b in range(B):
        w = wp[b]
        for i in range(C):
            bvec[i] += M[b, i] * w * dp[b]
            for j in range(C):
                A[i, j] += M[b, i] * w * M[b, j]
    if mode == 0:
        for i in range(C):
            A[i, i] += tp[i]
            bvec[i] += tp[i] * mp[i]
    ok = True
    logdet = 0.0
    for i in range(C):
        for j in range(i + 1):
            s = A[i, j]
            for k in range(j):
                s -= L[i, k] * L[j, k]
            if i == j:
                if s <= 0.0:
                    ok = False
                    s = 1.0
                L[i, i] = math.sqrt(s)
                logdet += 2.0 * math.log(L[i, i])
            else:
                L[i, j] = s / L[j, j]
    if not ok:
        return 0.0, False
    for i in range(C):
        s = bvec[i]
        for k in range(i):
            s -= L[i, k] * z[k]
        z[i] = s / L[i, i]
    quad = 0.0
    for i in range(C):
        quad += z[i] * z[i]
    term = 0.5 * quad
    if mode != 1:
        term -= 0.5 * logdet
    return term, True


_pixel_term = _maybe_jit(_pixel_term_py)


def _make_cores(pixel_term):
    def loglike_core(nu, cmb, dust, d, inv_var, beta_data, mode, m, inv_T):
        N = d.shape[0]
        B = d.shape[1]
        M = np.empty((B, 3))
        A = np.empty((3, 3))
        bvec = np.empty(3)
        L = np.zeros((3, 3))
        z = np.empty(3)
        total = 0.0
        nfail = 0
        for p in range(N):
            term, ok = pixel_term(nu, cmb, dust, d[p], inv_var[p], beta_data[p], mode, m[p],
                                  inv_T[p], M, A, bvec, L, z)
            if ok:
                total += term
            else:
                nfail += 1
        return total, nfail

    def table_core(nu, cmb, dust, d, inv_var, region, n_beta, beta_grid, mode, m, inv_T):
        N = d.shape[0]
        B = d.shape[1]
        G = beta_grid.shape[0]
        M = np.empty((B, 3))
        A = np.empty((3, 3))
        bvec = np.empty(3)
        L = np.zeros((3, 3))
        z = np.empty(3)
        table = np.zeros((n_beta, G))
        nfail = 0
        for j in range(G):
            for p in range(N):
                term, ok = pixel_term(nu, cmb, dust, d[p], inv_var[p], beta_grid[j], mode, m[p],
                                      inv_T[p], M, A, bvec, L, z)
                if ok:
                    table[region[p], j] += term
                else:
                    nfail += 1
        return table, nfail

    return loglike_core, table_core


_loglike_core_py, _table_core_py = _make_cores(_pixel_term_py)
_loglike_core, _table_core = (_maybe_jit(f) for f in _make_cores(_pixel_term))


def _prep(data, mode):
    if mode not in MODE_INT:
        raise ValueError(f"mode must be one of {sorted(MODE_INT)}, got {mode!r}")
    nu = np.ascontiguousarray(data.nu, dtype=np.float64)
    d = np.ascontiguousarray(data.d, dtype=np.float64)
    inv_var = np.ascontiguousarray(data.inv_var, dtype=np.float64)
    N = d.shape[0]
    if mode == "gauss":
        if data.m is None or data.inv_T is None:
            raise ValueError("mode 'gauss' needs m and inv_T")
        m = np.ascontiguousarray(data.m, dtype=np.float64)
        inv_T = np.ascontiguousarray(data.inv_T, dtype=np.float64)
        if m.shape != (N, 3) or inv_T.shape != (N, 3):
            raise ValueError("m and inv_T must have shape (N, 3)")
    else:
        m = np.zeros((N, 3))
        inv_T = np.zeros((N, 3))
    cmb, dust = _fixed_columns(nu, data.beta_d, data.T_d)
    return nu, cmb, dust, d, inv_var, MODE_INT[mode], m, inv_T


def loglike_pixel_loops(beta_field_data, data, mode, _core=None) -> float:
    """Collapsed log likelihood summed over data pixels, plain loops, no shared code with model.py.

    beta_field_data: (N,) synchrotron index per DATA pixel (beta[region]). Modes and the dropped
    constant are stated in the module docstring. Raises LinAlgError if any A (or A') is not
    positive definite, as E012 v2 requires (a skipped pixel is an error, not a change of target).
    """
    core = _loglike_core if _core is None else _core
    nu, cmb, dust, d, inv_var, mode_int, m, inv_T = _prep(data, mode)
    beta_data = np.ascontiguousarray(beta_field_data, dtype=np.float64)
    if beta_data.shape != (d.shape[0],):
        raise ValueError("beta_field_data must have shape (N,)")
    total, nfail = core(nu, cmb, dust, d, inv_var, beta_data, mode_int, m, inv_T)
    if nfail:
        raise np.linalg.LinAlgError(f"{nfail} pixels with a non positive definite A")
    return float(total)


def region_loglike_table(data, mode, beta_grid, _core=None) -> np.ndarray:
    """Table T[r, j] = sum over data pixels p in index pixel r of the pixel log likelihood at
    beta = beta_grid[j]. Since the collapsed likelihood is a sum over pixels and each pixel sees
    only beta_{r(p)}, log L(beta) = sum_r T[r, beta_r] exactly, which grid_1d and grid_2d use."""
    core = _table_core if _core is None else _core
    nu, cmb, dust, d, inv_var, mode_int, m, inv_T = _prep(data, mode)
    region = np.ascontiguousarray(data.region, dtype=np.int64)
    grid = np.ascontiguousarray(beta_grid, dtype=np.float64)
    table, nfail = core(nu, cmb, dust, d, inv_var, region, int(data.n_beta), grid, mode_int, m,
                        inv_T)
    if nfail:
        raise np.linalg.LinAlgError(f"{nfail} (pixel, grid point) pairs with non PD A")
    return table


def _log_gauss_prior(beta_points, theta_fixed):
    """Optional Gaussian prior on the index pixels for the grids.

    theta_fixed None: flat prior (returns 0). (mu, sd) scalars: iid N(mu, sd^2) on every index
    pixel, normalized. (mu_vec, cov): full multivariate normal N(mu_vec, cov), normalized.
    beta_points has shape (..., Nb).
    """
    beta_points = np.asarray(beta_points, dtype=np.float64)
    if theta_fixed is None:
        return np.zeros(beta_points.shape[:-1])
    mu, scale = theta_fixed
    nb = beta_points.shape[-1]
    if np.ndim(scale) == 0:
        r = (beta_points - float(mu)) / float(scale)
        return -0.5 * np.sum(r * r, axis=-1) - nb * math.log(float(scale)) \
            - 0.5 * nb * math.log(2.0 * math.pi)
    cov = np.asarray(scale, dtype=np.float64)
    res = beta_points - np.asarray(mu, dtype=np.float64)
    sol = np.linalg.solve(cov, res[..., None])[..., 0]
    _, logdet = np.linalg.slogdet(cov)
    return -0.5 * np.sum(res * sol, axis=-1) - 0.5 * logdet - 0.5 * nb * math.log(2.0 * math.pi)


def grid_1d(data, mode, beta_grid, theta_fixed=None) -> np.ndarray:
    """Unnormalized log posterior of the scalar index on beta_grid for N_beta = 1.

    log post(beta_j) = log L(beta_j) (dropped constant as in the module docstring) + log prior,
    where theta_fixed is None (flat) or (mu, sd) for N(mu, sd^2); see _log_gauss_prior.
    """
    if int(data.n_beta) != 1:
        raise ValueError("grid_1d needs n_beta == 1")
    grid = np.asarray(beta_grid, dtype=np.float64)
    table = region_loglike_table(data, mode, grid)
    return table[0] + _log_gauss_prior(grid[:, None], theta_fixed)


def grid_2d(data, mode, b1_grid, b2_grid, theta_fixed=None, brute=False) -> np.ndarray:
    """Unnormalized log posterior on the (b1, b2) grid for N_beta = 2, shape (G1, G2).

    log L(b1, b2) = T[0, b1] + T[1, b2] with the per-region tables of region_loglike_table (the
    pixels of the two regions do not interact in the likelihood); the prior couples them only
    when theta_fixed gives a full covariance. brute=True evaluates loglike_pixel_loops at every
    grid pair instead (self-test of the separable path).
    """
    if int(data.n_beta) != 2:
        raise ValueError("grid_2d needs n_beta == 2")
    g1 = np.asarray(b1_grid, dtype=np.float64)
    g2 = np.asarray(b2_grid, dtype=np.float64)
    if brute:
        region = np.asarray(data.region)
        ll = np.empty((g1.size, g2.size))
        for i, b1 in enumerate(g1):
            for j, b2 in enumerate(g2):
                ll[i, j] = loglike_pixel_loops(np.array([b1, b2])[region], data, mode)
    else:
        ll = region_loglike_table(data, mode, g1)[0][:, None] \
            + region_loglike_table(data, mode, g2)[1][None, :]
    pts = np.stack(np.meshgrid(g1, g2, indexing="ij"), axis=-1)
    return ll + _log_gauss_prior(pts, theta_fixed)


# ----------------------------------------------------------------------------------------------
# Basis, spectrum, rho transform, hyperpriors
# ----------------------------------------------------------------------------------------------
def region_map(nside_d, nside_beta) -> np.ndarray:
    """Index pixel r(p) (RING, nside_beta) containing data pixel p (RING, nside_d): nested
    parent pixel, as SPEC 'Shapes and conventions' prescribes."""
    if nside_beta > nside_d or nside_d % nside_beta:
        raise ValueError("nside_beta must divide nside_d")
    shift = 2 * int(round(math.log2(nside_d // nside_beta)))
    p = np.arange(12 * nside_d * nside_d)
    return np.asarray(hp.nest2ring(nside_beta, hp.ring2nest(nside_d, p) >> shift), dtype=np.int64)


class FieldSpec(NamedTuple):
    """Orthonormal basis V (Nb, K), eigenvalues lam (K,) >= 0 of the (continuum-scaled)
    Laplacian, mean pixel spacing h = sqrt(4 pi / Nb), K, Nb, nside_beta, family and, for the
    harmonic family, the multipole ell of each column (None for the graph Laplacian)."""

    V: np.ndarray
    lam: np.ndarray
    h: float
    K: int
    n_beta: int
    nside_beta: int
    family: str
    ell: np.ndarray | None


def _real_ylm(ell, m, theta, phi):
    """Real orthonormal spherical harmonic from the associated Legendre function lpmv:
    m = 0: sqrt((2l+1)/4pi) P_l(cos theta); m > 0: sqrt(2) N_lm P_l^m cos(m phi);
    m < 0: sqrt(2) N_l|m| P_l^|m| sin(|m| phi), N_lm^2 = (2l+1)/(4pi) (l-m)!/(l+m)!."""
    from scipy.special import gammaln, lpmv

    am = abs(m)
    norm = math.sqrt((2 * ell + 1) / (4 * math.pi)
                     * math.exp(gammaln(ell - am + 1) - gammaln(ell + am + 1)))
    P = lpmv(am, ell, np.cos(theta))
    if m == 0:
        return norm * P
    if m > 0:
        return math.sqrt(2.0) * norm * P * np.cos(m * phi)
    return math.sqrt(2.0) * norm * P * np.sin(am * phi)


def field_spec(nside_beta, family, lmax=None, remove_constant=True) -> FieldSpec:
    """Own construction of the basis and eigenvalues of E012 v2.

    family "graph_laplacian" (alias "gmrf"): W = symmetrized 8-neighbour adjacency of the
    HEALPix pixels at nside_beta (hp.get_all_neighbours, -1 entries dropped), G = D - W, eigh(G);
    the eigenvector with eigenvalue 0 (constant mode) is dropped, so K = Nb - 1; lam = eigenvalue
    / (3 h^2) with h = sqrt(4 pi / Nb) (D009 continuum scale).
    family "harmonic": real Y_lm for 1 <= l <= lmax (default 2 nside_beta) at the pixel centres,
    columns ordered by l then m = -l..l. With remove_constant=True (SPEC decision "Constant mode
    in the harmonic basis") every pixel-evaluated column y is first replaced by
    y - (1^T y / Nb) 1, since the l >= 1 harmonics sampled on the HEALPix grid are not exactly
    orthogonal to 1 (|1^T Y| / sqrt(Nb) up to 0.17 at nside 1), so dropping only l = 0 would
    leave part of the constant mode in V; the columns are then orthonormalized by QR (the QR
    keeps the span of the leading columns, so column (l, m) mixes only with lower or equal l),
    which gives 1^T V = 0 exactly and range(V) orthogonal to 1. remove_constant=False keeps the
    previous behaviour (QR of the raw columns). lam = l(l+1) of the column's l;
    K = (lmax+1)^2 - 1, which must be < Nb.
    """
    n_beta = 12 * nside_beta * nside_beta
    h = math.sqrt(4.0 * math.pi / n_beta)
    if family in ("graph_laplacian", "gmrf"):
        nb = hp.get_all_neighbours(nside_beta, np.arange(n_beta))
        W = np.zeros((n_beta, n_beta))
        for i in range(n_beta):
            for q in nb[:, i]:
                if q >= 0 and q != i:
                    W[i, q] = 1.0
        W = np.maximum(W, W.T)
        G = np.diag(W.sum(axis=1)) - W
        evals, evecs = np.linalg.eigh(G)
        if not (abs(evals[0]) < 1e-9 and evals[1] > 1e-9):
            raise RuntimeError("graph Laplacian: zero eigenvalue not unique or not found")
        V = np.ascontiguousarray(evecs[:, 1:])
        lam = evals[1:] / (3.0 * h * h)
        return FieldSpec(V, lam, h, n_beta - 1, n_beta, int(nside_beta), "graph_laplacian", None)
    if family == "harmonic":
        if lmax is None:
            lmax = 2 * nside_beta
        K = (lmax + 1) ** 2 - 1
        if K >= n_beta:
            raise ValueError("harmonic basis needs (lmax+1)^2 - 1 < Nb")
        theta, phi = hp.pix2ang(nside_beta, np.arange(n_beta))
        cols, ells = [], []
        for ell in range(1, lmax + 1):
            for m in range(-ell, ell + 1):
                cols.append(_real_ylm(ell, m, theta, phi))
                ells.append(ell)
        Y = np.stack(cols, axis=1)
        if remove_constant:
            Y = Y - Y.mean(axis=0, keepdims=True)
        Q, R = np.linalg.qr(Y)
        rdiag = np.abs(np.diag(R))
        if rdiag.min() < 1e-8 * rdiag.max():
            raise RuntimeError("harmonic basis: columns linearly dependent after projection")
        Q = Q * np.sign(np.diag(R))[None, :]
        ell_arr = np.array(ells, dtype=np.int64)
        return FieldSpec(np.ascontiguousarray(Q), (ell_arr * (ell_arr + 1)).astype(np.float64), h,
                         K, n_beta, int(nside_beta), "harmonic", ell_arr)
    raise ValueError(f"unknown family {family!r}")


def matern_spectrum(fspec, rho, alpha=2.0):
    """Matern spectrum s_i = (kappa^2 + lam_i)^(-alpha/2), kappa = sqrt(8 (alpha - 1)) / rho, and
    the normalization c with c^2 = (1/Nb) sum_i s_i^2, so that beta - mu 1 = sigma V diag(s/c) z
    has pixel-averaged variance exactly sigma^2. rho may be an array (leading axes)."""
    rho = np.asarray(rho, dtype=np.float64)
    kappa2 = 8.0 * (alpha - 1.0) / (rho * rho)
    s = (kappa2[..., None] + fspec.lam) ** (-0.5 * alpha)
    c = np.sqrt(np.sum(s * s, axis=-1) / fspec.n_beta)
    return s, c


def logdet_K(fspec, rho, alpha=2.0):
    """log det of K(rho) = V diag((s/c)^2) V^T restricted to range(V): sum_i 2 log(s_i / c)."""
    s, c = matern_spectrum(fspec, rho, alpha)
    return np.sum(2.0 * np.log(s), axis=-1) - 2.0 * fspec.K * np.log(c)


def _expit(u):
    return 0.5 * (1.0 + np.tanh(0.5 * np.asarray(u, dtype=np.float64)))


def rho_from_u(u, rho_min, rho_max):
    """rho = exp(log rho_min + L s(u)), L = log(rho_max / rho_min), s the logistic function."""
    L = math.log(rho_max / rho_min)
    return np.exp(math.log(rho_min) + L * _expit(u))


def u_from_rho(rho, rho_min, rho_max):
    """u = logit((log rho - log rho_min) / L)."""
    L = math.log(rho_max / rho_min)
    t = (np.log(rho) - math.log(rho_min)) / L
    return np.log(t) - np.log1p(-t)


def log_dlogrho_du(u, rho_min, rho_max):
    """log |d log rho / d u| = log L + log s(u) + log(1 - s(u)), computed stably."""
    u = np.asarray(u, dtype=np.float64)
    return math.log(math.log(rho_max / rho_min)) - np.logaddexp(0.0, -u) - np.logaddexp(0.0, u)


class RefHyper(NamedTuple):
    """Hyperprior settings of E012 v2 (same attribute names as SPEC HyperPrior).

    mu ~ N(mu0, s_mu^2); sigma ~ half-Normal(0, sigma_scale) (rho_prior "loguniform") or
    Exponential with P(sigma > pc_sigma0) = 0.05 (rho_prior "pc", the joint PC prior of Fuglstad et
    al. 2019); rho log-uniform on [rho_min, rho_max] or PC with P(rho < pc_rho0) = 0.05, truncated
    to the same bounds. Both are expressed on log sigma and on u with their Jacobians.
    """

    mu0: float = -3.0
    s_mu: float = 0.5
    sigma_scale: float = 0.3
    rho_prior: str = "loguniform"
    rho_min: float = 1.0
    rho_max: float = math.pi
    pc_rho0: float | None = None
    pc_sigma0: float | None = None


def default_hyper(fspec, **overrides) -> RefHyper:
    """E012 defaults for a basis: rho_min = h, rho_max = pi, pc_rho0 = 2 h, pc_sigma0 = 0.5."""
    kw = dict(rho_min=fspec.h, rho_max=math.pi, pc_rho0=2.0 * fspec.h, pc_sigma0=0.5)
    kw.update(overrides)
    return RefHyper(**kw)


def _lp_mu(mu, hyper):
    """log N(mu; mu0, s_mu^2) up to a constant."""
    r = (np.asarray(mu, dtype=np.float64) - hyper.mu0) / hyper.s_mu
    return -0.5 * r * r


def _lp_logsigma(log_sigma, hyper):
    """log p(log sigma) = log p(sigma) + log sigma (Jacobian), constants dropped.
    half-normal: -sigma^2 / (2 sigma_scale^2) + log sigma; PC (exponential): -lambda_2 sigma +
    log sigma, lambda_2 = -log(0.05) / pc_sigma0."""
    log_sigma = np.asarray(log_sigma, dtype=np.float64)
    sigma = np.exp(log_sigma)
    if getattr(hyper, "rho_prior", "loguniform") == "pc":
        lam2 = -math.log(0.05) / float(hyper.pc_sigma0)
        return -lam2 * sigma + log_sigma
    return -0.5 * (sigma / hyper.sigma_scale) ** 2 + log_sigma


def _lp_u(u, hyper):
    """log p(u) = log p(log rho) + log |d log rho / d u|, constants dropped.
    log-uniform: p(log rho) = 1/L, so log p(u) = log s(u) + log(1 - s(u)) (+ const).
    PC (d = 2): p(rho) = lambda_1 rho^-2 exp(-lambda_1 / rho), lambda_1 = -pc_rho0 log(0.05),
    p(log rho) = p(rho) rho, truncated to [rho_min, rho_max] (normalization dropped)."""
    u = np.asarray(u, dtype=np.float64)
    base = -np.logaddexp(0.0, -u) - np.logaddexp(0.0, u)
    if getattr(hyper, "rho_prior", "loguniform") == "pc":
        rho = rho_from_u(u, hyper.rho_min, hyper.rho_max)
        lam1 = -float(hyper.pc_rho0) * math.log(0.05)
        base = base + math.log(math.log(hyper.rho_max / hyper.rho_min)) \
            + math.log(lam1) - np.log(rho) - lam1 / rho
    return base


def log_hyperprior(mu, log_sigma, u, hyper):
    """log p(mu) + log p(log sigma) + log p(u) in the sampled coordinates (broadcasting)."""
    return _lp_mu(mu, hyper) + _lp_logsigma(log_sigma, hyper) + _lp_u(u, hyper)


# ----------------------------------------------------------------------------------------------
# Exact 2-D grid draw of (log sigma, u) | y
# ----------------------------------------------------------------------------------------------
class HyperGrid(NamedTuple):
    ls: np.ndarray      # (n_ls,) cell centres in log sigma
    u: np.ndarray       # (n_u,) cell centres in u
    dls: float
    du: float
    prob: np.ndarray    # (n_ls, n_u) normalized cell probabilities
    logdens: np.ndarray


def hyper_conditional_logdensity(ls, uu, y, fspec, hyper, alpha=2.0):
    """log p(log sigma, u | y) on the product grid ls (n_ls,) x uu (n_u,), constants dropped:
    -K log sigma - 1/2 logdet K(rho) - 1/2 y^T diag((c/s)^2) y / sigma^2
    + log p(log sigma) + log p(u), where y = V^T (beta - mu 1) are the basis coefficients and
    K = number of basis modes. The likelihood does not enter (beta is fixed)."""
    y = np.asarray(y, dtype=np.float64)
    ls = np.asarray(ls, dtype=np.float64)
    uu = np.asarray(uu, dtype=np.float64)
    rho = rho_from_u(uu, hyper.rho_min, hyper.rho_max)
    s, c = matern_spectrum(fspec, rho, alpha)                       # (n_u, K), (n_u,)
    logdetK = np.sum(2.0 * np.log(s), axis=-1) - 2.0 * fspec.K * np.log(c)
    Q = np.sum((y * y)[None, :] * (c * c)[:, None] / (s * s), axis=-1)
    inv_s2 = np.exp(-2.0 * ls)
    return (-fspec.K * ls[:, None] - 0.5 * logdetK[None, :]
            - 0.5 * Q[None, :] * inv_s2[:, None]
            + _lp_logsigma(ls, hyper)[:, None] + _lp_u(uu, hyper)[None, :])


def hyper_conditional_grid(y, fspec, hyper, alpha=2.0, n_coarse=(128, 128), n_fine=(256, 256),
                           ls_range=(-14.0, 3.0), u_range=(-16.0, 16.0), drop=30.0) -> HyperGrid:
    """Two-pass grid of p(log sigma, u | y): a coarse pass over the full ranges finds the box
    where the log density is within `drop` of its maximum (extended by one coarse cell), then a
    fine midpoint grid in that box gives cell probabilities proportional to the density at the
    cell centre. The draw from this grid (draw_from_grid) is exact for the piecewise-constant
    approximation of the conditional, with O(h^2) discretization error in the density."""
    ls_c = np.linspace(ls_range[0], ls_range[1], n_coarse[0])
    u_c = np.linspace(u_range[0], u_range[1], n_coarse[1])
    ld = hyper_conditional_logdensity(ls_c, u_c, y, fspec, hyper, alpha)
    mask = ld > ld.max() - drop
    ii, jj = np.where(mask)
    dls_c = ls_c[1] - ls_c[0]
    du_c = u_c[1] - u_c[0]
    lo_ls = max(ls_range[0], ls_c[ii.min()] - dls_c)
    hi_ls = min(ls_range[1], ls_c[ii.max()] + dls_c)
    lo_u = max(u_range[0], u_c[jj.min()] - du_c)
    hi_u = min(u_range[1], u_c[jj.max()] + du_c)
    dls = (hi_ls - lo_ls) / n_fine[0]
    du = (hi_u - lo_u) / n_fine[1]
    ls_f = lo_ls + dls * (np.arange(n_fine[0]) + 0.5)
    u_f = lo_u + du * (np.arange(n_fine[1]) + 0.5)
    ldf = hyper_conditional_logdensity(ls_f, u_f, y, fspec, hyper, alpha)
    prob = np.exp(ldf - ldf.max())
    prob /= prob.sum()
    return HyperGrid(ls_f, u_f, dls, du, prob, ldf)


def draw_from_grid(grid: HyperGrid, rng, n=1):
    """Inverse-CDF draw of n cells from grid.prob, with a uniform jitter inside each cell.
    Returns (log_sigma (n,), u (n,), i (n,), j (n,)) with the cell indices."""
    cdf = np.cumsum(grid.prob.ravel())
    cdf /= cdf[-1]
    flat = np.searchsorted(cdf, rng.random(n), side="right")
    flat = np.minimum(flat, cdf.size - 1)
    i, j = np.unravel_index(flat, grid.prob.shape)
    ls = grid.ls[i] + (rng.random(n) - 0.5) * grid.dls
    uu = grid.u[j] + (rng.random(n) - 0.5) * grid.du
    return ls, uu, i, j


# ----------------------------------------------------------------------------------------------
# Effective sample size
# ----------------------------------------------------------------------------------------------
def ess_geyer(x) -> float:
    """ESS = n / tau with tau = -1 + 2 sum_k Gamma_k / gamma_0, Gamma_k = gamma_{2k} +
    gamma_{2k+1} the paired autocovariances, truncated at the first non-positive Gamma_k
    (Geyer's initial positive sequence) and made non-increasing (initial monotone sequence)."""
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    x = x - x.mean()
    if n < 4 or not np.any(x != 0.0):
        return float("nan")
    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    f = np.fft.rfft(x, nfft)
    acov = np.fft.irfft(f * np.conj(f), nfft)[:n] / n
    gamma0 = acov[0]
    npairs = n // 2
    Gamma = acov[0:2 * npairs:2] + acov[1:2 * npairs:2]
    total = 0.0
    prev = np.inf
    for g in Gamma:
        if g <= 0.0:
            break
        g = min(g, prev)
        prev = g
        total += g
    tau = -1.0 + 2.0 * total / gamma0
    return float(n / max(tau, 1e-12))


# ----------------------------------------------------------------------------------------------
# Centred Gibbs reference
# ----------------------------------------------------------------------------------------------
def centred_gibbs_reference(data, fspec, hyper, mode, rng, n_iter, n_burn=None, n_mh_sub=5,
                            alpha=2.0, init=None, init_sd=0.05, grid_kwargs=None, verbose=False):
    """MH-within-Gibbs reference chain in centred coordinates (E012 v2, Reference 1c).

    Target. With the constant mode dropped from the basis (E012 v2) the field prior
    N(beta; mu 1, sigma^2 K(rho)), K = V diag((s/c)^2) V^T, is supported on the affine subspace
    mu 1 + range(V); it is enforced by the parameterization beta = mu 1 + V y, y in R^K, and the
    target is

      p(y, mu, log sigma, u | d) propto L(mu 1 + V y) N(y; 0, sigma^2 diag((s/c)^2))
                                        p(mu) p(log sigma) p(u),

    with log N(y; 0, sigma^2 D) = -K log sigma - 1/2 logdet K(rho) - 1/2 y^T D^-1 y / sigma^2
    + const (K = number of basis modes), p(log sigma) the half-normal on sigma times the log sigma
    Jacobian and p(u) the rho prior times the logit Jacobian (log_hyperprior). This is the same
    target as SPEC CentredPosterior and the same joint as the non-centred posterior.

    Steps per iteration.
    1. (y, mu) | sigma, rho, d: random-walk Metropolis on w = (y, mu) in R^(K+1), n_mh_sub
       proposals per iteration. The proposal covariance (scaled empirical covariance, Haario) and
       the scale (Robbins-Monro towards acceptance 0.234) adapt during the burn-in only; the
       burn-in is discarded. mu is moved together with y because beta pins mu: given beta the
       decomposition beta = mu 1 + V y is unique (1 is not in range(V)), so the "exact Gaussian
       mu | beta" step of a full-rank prior does not exist here (it would be a delta). This
       deviation from the SPEC wording is deliberate and keeps the sampler valid.
    2. (log sigma, u) | y, mu: exact inverse-CDF draw from the fine 2-D grid of
       hyper_conditional_grid (the likelihood does not enter).

    mode: "gauss", "w0", "w1" (SED likelihood through loglike_pixel_loops, raising on a non PD
    A as E012 requires) or "linear" (d = offset + F beta_region + n, needs data.F and data.offset).

    Returns a dict with samples after burn-in: beta (n_iter, Nb), y (n_iter, K), mu, log_sigma,
    log_rho, u (n_iter,), accept_rate, ess (dict with mu, log_sigma, log_rho, u, beta (Nb,)),
    and n_burn, K, seconds.
    """
    Nb, K, V = fspec.n_beta, fspec.K, fspec.V
    if int(data.n_beta) != Nb:
        raise ValueError("data.n_beta does not match fspec.n_beta")
    region = np.asarray(data.region, dtype=np.int64)
    if n_burn is None:
        n_burn = max(1000, n_iter // 2)
    grid_kwargs = {} if grid_kwargs is None else dict(grid_kwargs)

    if mode == "linear":
        if data.F is None or data.offset is None:
            raise ValueError("mode 'linear' needs data.F and data.offset")
        d, F, off, iv = data.d, data.F, data.offset, data.inv_var

        def loglike(beta):
            r = d - off - F * beta[region][:, None]
            return -0.5 * float(np.sum(r * r * iv))
    else:
        nu, cmb, dust, d, inv_var, mode_int, m, inv_T = _prep(data, mode)

        def loglike(beta):
            total, nfail = _loglike_core(nu, cmb, dust, d, inv_var, beta[region], mode_int, m,
                                         inv_T)
            if nfail:
                raise np.linalg.LinAlgError(f"{nfail} pixels with a non positive definite A")
            return float(total)

    def prior_w(w, dinv_s2):
        y, mu = w[:K], w[K]
        return -0.5 * float(np.sum(y * y * dinv_s2)) + float(_lp_mu(mu, hyper))

    def dinv_over_sigma2(ls, u):
        s, c = matern_spectrum(fspec, rho_from_u(u, hyper.rho_min, hyper.rho_max), alpha)
        return (c * c) / (s * s) * math.exp(-2.0 * ls)

    # initial state
    init = {} if init is None else dict(init)
    mu = float(init.get("mu", hyper.mu0))
    y = np.zeros(K)
    if "beta" in init:
        sol, *_ = np.linalg.lstsq(np.column_stack([np.ones(Nb), V]), np.asarray(init["beta"]),
                                  rcond=None)
        mu, y = float(sol[0]), sol[1:].copy()
    if "y" in init:
        y = np.asarray(init["y"], dtype=np.float64).copy()
    ls = float(init.get("log_sigma", math.log(0.1)))
    u = float(init.get("u", 0.0))
    w = np.concatenate([y, [mu]])
    dim = K + 1
    dinv_s2 = dinv_over_sigma2(ls, u)
    cur_ll = loglike(mu + V @ y)
    cur_lp = cur_ll + prior_w(w, dinv_s2)

    # adaptation state
    log_scale = math.log(2.38 / math.sqrt(dim))
    chol = np.eye(dim) * init_sd
    mean_w = w.copy()
    cov_acc = np.zeros((dim, dim))
    count = 0
    t_adapt = 0
    reset_at = n_burn // 3

    out_beta = np.empty((n_iter, Nb))
    out_y = np.empty((n_iter, K))
    out_mu = np.empty(n_iter)
    out_ls = np.empty(n_iter)
    out_u = np.empty(n_iter)
    n_acc = 0
    n_prop = 0
    t0 = time.time()
    for it in range(n_burn + n_iter):
        burn = it < n_burn
        if it == reset_at:
            mean_w = w.copy()
            cov_acc[:] = 0.0
            count = 0
        for _ in range(n_mh_sub):
            prop = w + math.exp(log_scale) * (chol @ rng.standard_normal(dim))
            ll_prop = loglike(prop[K] + V @ prop[:K])
            lp_prop = ll_prop + prior_w(prop, dinv_s2)
            log_r = lp_prop - cur_lp
            acc_prob = min(1.0, math.exp(min(0.0, log_r)))
            if math.log(rng.random()) < log_r:
                w, cur_ll, cur_lp = prop, ll_prop, lp_prop
                if not burn:
                    n_acc += 1
            if not burn:
                n_prop += 1
            if burn:
                t_adapt += 1
                log_scale += (acc_prob - 0.234) / (t_adapt + 20.0) ** 0.6
                count += 1
                delta = w - mean_w
                mean_w = mean_w + delta / count
                cov_acc += np.outer(delta, w - mean_w)
        if burn and count > 2 * dim + 10 and it % 20 == 0:
            cov = cov_acc / (count - 1) + 1e-12 * np.eye(dim)
            chol = np.linalg.cholesky(cov)
        # exact conditional draw of (log sigma, u) given y
        grid = hyper_conditional_grid(w[:K], fspec, hyper, alpha, **grid_kwargs)
        ls_draw, u_draw, _, _ = draw_from_grid(grid, rng, 1)
        ls, u = float(ls_draw[0]), float(u_draw[0])
        dinv_s2 = dinv_over_sigma2(ls, u)
        cur_lp = cur_ll + prior_w(w, dinv_s2)
        if not burn:
            k = it - n_burn
            out_y[k] = w[:K]
            out_mu[k] = w[K]
            out_beta[k] = w[K] + V @ w[:K]
            out_ls[k] = ls
            out_u[k] = u
        if verbose and (it + 1) % max(1, (n_burn + n_iter) // 10) == 0:
            print(f"  gibbs it {it + 1}/{n_burn + n_iter}  acc {n_acc / max(n_prop, 1):.3f}  "
                  f"scale {math.exp(log_scale):.3g}", file=sys.stderr)
    log_rho = np.log(rho_from_u(out_u, hyper.rho_min, hyper.rho_max))
    ess = dict(mu=ess_geyer(out_mu), log_sigma=ess_geyer(out_ls), log_rho=ess_geyer(log_rho),
               u=ess_geyer(out_u), beta=np.array([ess_geyer(out_beta[:, r]) for r in range(Nb)]))
    return dict(beta=out_beta, y=out_y, mu=out_mu, log_sigma=out_ls, log_rho=log_rho, u=out_u,
                accept_rate=n_acc / max(n_prop, 1), ess=ess, n_burn=n_burn, K=K,
                seconds=time.time() - t0)


# ----------------------------------------------------------------------------------------------
# Linear-Gaussian closed form
# ----------------------------------------------------------------------------------------------
class ClosedForm(NamedTuple):
    mu_grid: np.ndarray
    logsig_grid: np.ndarray
    u_grid: np.ndarray
    logpost: np.ndarray          # (n_mu, n_ls, n_u) unnormalized log p(theta | d)
    weights: np.ndarray          # normalized grid weights
    theta_mean: np.ndarray       # (mu, log sigma, u, log rho) posterior means
    theta_sd: np.ndarray
    beta_mean: np.ndarray        # (Nb,) marginal posterior mean of beta
    beta_sd: np.ndarray
    beta_given_theta: Callable   # (mu, log_sigma, u) -> (mean (Nb,), cov (Nb, Nb))
    log_norm: float              # log sum of exp(logpost) times the cell volume


def linear_gaussian_closed_form(d, F, offset, inv_var, region, fspec, hyper, mu_grid, logsig_grid,
                                u_grid, alpha=2.0) -> ClosedForm:
    """Closed form for the fixed-amplitude linear model d_pb = offset_pb + F_pb beta_r(p) + n_pb,
    n ~ N(0, N), N = diag(1/inv_var) (E012 v2, Reference 1d).

    Write the data as a vector of length N B, X (N B x Nb) the design with X[(p,b), r(p)] = F_pb
    and U = X V. With beta = mu 1 + V y, y ~ N(0, sigma^2 D), D = diag((s/c)^2):

      y | theta, d ~ N(P^-1 U^T N^-1 r, P^-1),   P = U^T N^-1 U + D^-1 / sigma^2,
                                                 r = d - offset - mu X 1,
      beta | theta, d: mean mu 1 + V E[y], covariance V P^-1 V^T  (equal to E012's
                       Sigma (F^T N^-1 (d - m) + Q mu 1) with Q = (sigma^2 K)^-1 on range(V)),
      theta | d propto p(theta) N(d - offset; mu X 1, sigma^2 X K X^T + N),

    the marginal evaluated in K dimensions by the Woodbury identity and the matrix determinant
    lemma: log N = -1/2 [r^T N^-1 r - (U^T N^-1 r)^T P^-1 (U^T N^-1 r)] - 1/2 [log det N
    + log det(sigma^2 D) + log det P] - (N B / 2) log 2 pi. p(theta) is log_hyperprior on
    (mu, log sigma, u). The posterior of beta is the grid mixture of the conditionals; the grid
    quadrature is the only approximation (no Monte Carlo error). mu_grid, logsig_grid, u_grid are
    1-D arrays of grid points (treated as midpoints of equal cells).
    """
    d = np.asarray(d, dtype=np.float64)
    F = np.asarray(F, dtype=np.float64)
    offset = np.asarray(offset, dtype=np.float64)
    inv_var = np.asarray(inv_var, dtype=np.float64)
    region = np.asarray(region, dtype=np.int64)
    N, B = d.shape
    Nb, K, V = fspec.n_beta, fspec.K, fspec.V
    mu_grid = np.asarray(mu_grid, dtype=np.float64)
    ls_grid = np.asarray(logsig_grid, dtype=np.float64)
    u_grid = np.asarray(u_grid, dtype=np.float64)

    ivv = inv_var.ravel()
    r0 = (d - offset).ravel()
    Fvec = F.ravel()
    region_rep = np.repeat(region, B)
    U = Fvec[:, None] * V[region_rep]                          # (N B, K)
    UtNU = U.T @ (ivv[:, None] * U)
    a0 = U.T @ (ivv * r0)
    a1 = U.T @ (ivv * Fvec)
    r0Nr0 = float(np.sum(ivv * r0 * r0))
    x1Nr0 = float(np.sum(ivv * Fvec * r0))
    x1Nx1 = float(np.sum(ivv * Fvec * Fvec))
    logdetN = -float(np.sum(np.log(ivv)))
    const = -0.5 * logdetN - 0.5 * N * B * math.log(2.0 * math.pi)

    def _P_and_logdet(ls_vec, u_scalar):
        rho = rho_from_u(u_scalar, hyper.rho_min, hyper.rho_max)
        s, c = matern_spectrum(fspec, rho, alpha)
        logdetK = float(np.sum(2.0 * np.log(s)) - 2.0 * K * np.log(c))
        dinv = (c * c) / (s * s)
        P = UtNU[None] + np.einsum("j,k->jk", np.exp(-2.0 * ls_vec), dinv)[:, :, None] * np.eye(K)
        _, logdetP = np.linalg.slogdet(P)
        return P, logdetP, logdetK

    n_mu, n_ls, n_u = mu_grid.size, ls_grid.size, u_grid.size
    logpost = np.empty((n_mu, n_ls, n_u))
    lp_mu = _lp_mu(mu_grid, hyper)
    lp_ls = _lp_logsigma(ls_grid, hyper)
    lp_u = _lp_u(u_grid, hyper)
    rhs = np.stack([a0, a1], axis=1)                          # (K, 2)
    store = []
    for k in range(n_u):
        P, logdetP, logdetK = _P_and_logdet(ls_grid, u_grid[k])
        sol = np.linalg.solve(P, np.broadcast_to(rhs, (n_ls, K, 2)))   # (n_ls, K, 2)
        q00 = np.einsum("j,jk->k", a0, sol[:, :, 0].T)        # a0^T P^-1 a0 per ls
        q01 = np.einsum("j,jk->k", a0, sol[:, :, 1].T)
        q11 = np.einsum("j,jk->k", a1, sol[:, :, 1].T)
        mu = mu_grid[:, None]
        rNr = r0Nr0 - 2.0 * mu * x1Nr0 + mu * mu * x1Nx1
        quad = q00[None, :] - 2.0 * mu * q01[None, :] + mu * mu * q11[None, :]
        logdet_total = logdetN + 2.0 * K * ls_grid + logdetK + logdetP
        loglik = -0.5 * (rNr - quad) - 0.5 * logdet_total[None, :] + const + 0.5 * logdetN
        logpost[:, :, k] = loglik + lp_mu[:, None] + lp_ls[None, :] + lp_u[k]
        store.append((P, sol))

    lmax = logpost.max()
    weights = np.exp(logpost - lmax)
    cell = (mu_grid[1] - mu_grid[0]) * (ls_grid[1] - ls_grid[0]) * (u_grid[1] - u_grid[0]) \
        if min(n_mu, n_ls, n_u) > 1 else 1.0
    log_norm = float(lmax + math.log(weights.sum()) + math.log(cell))
    weights /= weights.sum()

    # theta moments
    logrho_grid = np.log(rho_from_u(u_grid, hyper.rho_min, hyper.rho_max))
    w_mu = weights.sum(axis=(1, 2))
    w_ls = weights.sum(axis=(0, 2))
    w_u = weights.sum(axis=(0, 1))

    def _msd(g, wg):
        m1 = float(np.sum(wg * g))
        return m1, math.sqrt(max(float(np.sum(wg * (g - m1) ** 2)), 0.0))

    mm, ms = zip(*[_msd(mu_grid, w_mu), _msd(ls_grid, w_ls), _msd(u_grid, w_u),
                   _msd(logrho_grid, w_u)])
    theta_mean, theta_sd = np.array(mm), np.array(ms)

    # beta marginal moments: mixture over the grid of N(mu 1 + V (p0 - mu p1), V P^-1 V^T)
    e1 = np.zeros(Nb)
    e2 = np.zeros(Nb)
    for k in range(n_u):
        P, sol = store[k]
        Pinv = np.linalg.inv(P)                                      # (n_ls, K, K)
        var = np.einsum("ri,jik,rk->jr", V, Pinv, V)                 # (n_ls, Nb)
        Avec = sol[:, :, 0] @ V.T                                    # (n_ls, Nb)  V p0
        Bvec = 1.0 - sol[:, :, 1] @ V.T                              # (n_ls, Nb)  1 - V p1
        wk = weights[:, :, k]                                        # (n_mu, n_ls)
        w0 = wk.sum(axis=0)
        w1 = (wk * mu_grid[:, None]).sum(axis=0)
        w2 = (wk * (mu_grid * mu_grid)[:, None]).sum(axis=0)
        e1 += (w0[:, None] * Avec + w1[:, None] * Bvec).sum(axis=0)
        e2 += (w0[:, None] * (Avec * Avec + var) + 2.0 * w1[:, None] * Avec * Bvec
               + w2[:, None] * Bvec * Bvec).sum(axis=0)
    beta_sd = np.sqrt(np.maximum(e2 - e1 * e1, 0.0))

    def beta_given_theta(mu, log_sigma, u):
        P, _, _ = _P_and_logdet(np.array([float(log_sigma)]), float(u))
        Pinv = np.linalg.inv(P[0])
        yv = Pinv @ (a0 - float(mu) * a1)
        return float(mu) + V @ yv, V @ Pinv @ V.T

    return ClosedForm(mu_grid, ls_grid, u_grid, logpost, weights, theta_mean, theta_sd, e1,
                      beta_sd, beta_given_theta, log_norm)


def simulate_linear(fspec, hyper, nside_d, B, rng, mu, sigma, rho, f_scale=50.0, noise_sd=3.0,
                    alpha=2.0):
    """Draw beta = mu 1 + sigma V diag(s/c) z and data d = offset + F beta_region + n for the
    linear model, F ~ U(0.5, 1.5) f_scale per (pixel, band), offset ~ N(0, 10^2), white noise
    noise_sd. Returns (RefData, dict(beta, z))."""
    Nb, K = fspec.n_beta, fspec.K
    s, c = matern_spectrum(fspec, rho, alpha)
    z = rng.standard_normal(K)
    beta = mu + sigma * (fspec.V @ (s / c * z))
    region = region_map(nside_d, fspec.nside_beta)
    N = region.size
    F = f_scale * rng.uniform(0.5, 1.5, size=(N, B))
    offset = 10.0 * rng.standard_normal((N, B))
    d = offset + F * beta[region][:, None] + noise_sd * rng.standard_normal((N, B))
    inv_var = np.full((N, B), 1.0 / noise_sd**2)
    data = from_arrays(np.linspace(20.0, 40.0, B), d, inv_var, region, Nb, F=F, offset=offset)
    return data, dict(beta=beta, z=z)


# ----------------------------------------------------------------------------------------------
# Self-tests
# ----------------------------------------------------------------------------------------------
def _selftest():
    from scipy import stats

    rng = np.random.default_rng(20261025)
    fails = []

    def check(name, ok, detail):
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
        if not ok:
            fails.append(name)

    # --- (a) pixel-loop likelihood against a direct NumPy evaluation, all modes, jit vs python
    nu = np.array([23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 94.0, 100.0, 143.0, 217.0, 353.0])
    B = nu.size
    N = 40
    cmb, dust = _fixed_columns(nu, 1.6, 20.0)
    a_true = np.column_stack([70.0 * rng.standard_normal(N),
                              10 ** rng.uniform(0.5, 2.5, N), 10 ** rng.uniform(0.5, 2.0, N)])
    beta_true = -3.0 + 0.2 * rng.standard_normal(N)
    sig = 2.0 + 3.0 * rng.random((N, B))
    M = np.empty((N, B, 3))
    M[:, :, 0] = cmb
    M[:, :, 1] = (nu[None, :] / NU0_S) ** beta_true[:, None]
    M[:, :, 2] = dust
    d = np.einsum("pbc,pc->pb", M, a_true) + sig * rng.standard_normal((N, B))
    inv_var = 1.0 / sig**2
    m_prior = np.column_stack([np.zeros(N), 10 ** rng.uniform(0.5, 2.5, N),
                               10 ** rng.uniform(0.5, 2.0, N)])
    inv_T = 1.0 / (10.0 * np.maximum(np.abs(m_prior), 10.0)) ** 2
    region = np.zeros(N, dtype=np.int64)
    data = from_arrays(nu, d, inv_var, region, 1, 1.6, 20.0, m_prior, inv_T)
    beta_eval = -2.9 + 0.1 * rng.standard_normal(N)
    Me = M.copy()
    Me[:, :, 1] = (nu[None, :] / NU0_S) ** beta_eval[:, None]
    A = np.einsum("pbi,pb,pbj->pij", Me, inv_var, Me)
    bvec = np.einsum("pbi,pb,pb->pi", Me, inv_var, d)
    for mode in ("gauss", "w0", "w1"):
        if mode == "gauss":
            Ap = A + np.einsum("pi,ij->pij", inv_T, np.eye(3))
            bp = bvec + inv_T * m_prior
        else:
            Ap, bp = A, bvec
        quad = np.einsum("pi,pi->p", bp, np.linalg.solve(Ap, bp[:, :, None])[:, :, 0])
        logdet = np.linalg.slogdet(Ap)[1]
        ref = 0.5 * quad.sum() - (0.0 if mode == "w0" else 0.5 * logdet.sum())
        val = loglike_pixel_loops(beta_eval, data, mode)
        val_py = loglike_pixel_loops(beta_eval, data, mode, _core=_loglike_core_py)
        check(f"loglike_pixel_loops {mode} vs numpy", abs(val - ref) < 1e-8 * abs(ref) + 1e-9,
              f"loops {val:.10g} numpy {ref:.10g} diff {val - ref:.2e}; jit-python diff "
              f"{val - val_py:.1e}")
    # non PD must raise
    bad = from_arrays(nu[:2], d[:, :2], inv_var[:, :2], region, 1)
    try:
        loglike_pixel_loops(beta_eval, bad, "w1")
        raised = False
    except np.linalg.LinAlgError:
        raised = True
    check("non PD A raises (B=2 < C=3)", raised, "LinAlgError raised" if raised else "no error")

    # --- (i) grid_1d mode gauss vs brute-force 3-D quadrature over a_p, 1 pixel, C = 3
    nu1 = np.array([23.0, 30.0, 44.0, 70.0, 100.0, 143.0])
    cmb1, dust1 = _fixed_columns(nu1, 1.6, 20.0)
    m1 = np.array([0.0, 60.0, 30.0])
    T1 = np.array([200.0, 600.0, 300.0]) ** 2
    a1 = m1 + np.sqrt(T1) * rng.standard_normal(3)
    beta0 = -3.1
    M1 = np.column_stack([cmb1, (nu1 / NU0_S) ** beta0, dust1])
    sig1 = np.array([3.0, 2.5, 2.5, 2.0, 2.0, 2.5])
    d1 = M1 @ a1 + sig1 * rng.standard_normal(6)
    data1 = from_arrays(nu1, d1[None, :], (1.0 / sig1**2)[None, :], np.zeros(1, dtype=np.int64),
                        1, 1.6, 20.0, m1[None, :], (1.0 / T1)[None, :])
    bgrid = np.linspace(-3.6, -2.6, 9)
    g1 = grid_1d(data1, "gauss", bgrid)
    Ninv = np.diag(1.0 / sig1**2)
    const = (-0.5 * d1 @ Ninv @ d1 - 0.5 * np.sum(m1 * m1 / T1) - 0.5 * np.sum(np.log(sig1**2))
             - 0.5 * np.sum(np.log(T1)) - 0.5 * 6 * math.log(2.0 * math.pi))
    brute = np.empty_like(bgrid)
    n_q = 121
    for j, bj in enumerate(bgrid):
        Mj = np.column_stack([cmb1, (nu1 / NU0_S) ** bj, dust1])
        Ap = Mj.T @ Ninv @ Mj + np.diag(1.0 / T1)
        bp = Mj.T @ Ninv @ d1 + m1 / T1
        centre = np.linalg.solve(Ap, bp)
        sd = np.sqrt(np.diag(np.linalg.inv(Ap)))
        axes = [np.linspace(centre[i] - 9 * sd[i], centre[i] + 9 * sd[i], n_q) for i in range(3)]
        vol = np.prod([ax[1] - ax[0] for ax in axes])
        G = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
        res = d1[None, :] - G @ Mj.T
        expo = (-0.5 * np.sum(res * res / sig1**2, axis=1)
                - 0.5 * np.sum((G - m1) ** 2 / T1, axis=1))
        expo += (-0.5 * np.sum(np.log(sig1**2)) - 0.5 * np.sum(np.log(T1))
                 - 0.5 * 9 * math.log(2.0 * math.pi))
        mx = expo.max()
        brute[j] = mx + math.log(np.sum(np.exp(expo - mx))) + math.log(vol)
        edge = np.exp(expo.reshape(n_q, n_q, n_q)[[0, -1], :, :].max() - mx)
        assert edge < 1e-12, edge
    diff = g1 + const - brute
    check("(i) grid_1d gauss = brute-force quadrature + stated constant",
          np.max(np.abs(diff)) < 1e-6,
          f"max |grid + const - quadrature| = {np.max(np.abs(diff)):.2e} over {bgrid.size} "
          f"beta values; log L range {g1.max() - g1.min():.3g}")
    # grid_2d separable path vs brute force
    region2 = (np.arange(N) % 2).astype(np.int64)
    data2 = from_arrays(nu, d, inv_var, region2, 2, 1.6, 20.0, m_prior, inv_T)
    b1g = np.linspace(-3.2, -2.8, 5)
    b2g = np.linspace(-3.3, -2.7, 4)
    for mode in ("gauss", "w1"):
        fast = grid_2d(data2, mode, b1g, b2g, theta_fixed=(-3.0, 0.3))
        slow = grid_2d(data2, mode, b1g, b2g, theta_fixed=(-3.0, 0.3), brute=True)
        check(f"grid_2d separable = brute {mode}", np.max(np.abs(fast - slow)) < 1e-8,
              f"max diff {np.max(np.abs(fast - slow)):.2e}")

    # --- (b) field_spec: orthonormality, zero mode dropped, eigen residual, harmonic span
    for nside in (1, 2):
        fs = field_spec(nside, "graph_laplacian")
        orth = np.max(np.abs(fs.V.T @ fs.V - np.eye(fs.K)))
        one = np.max(np.abs(fs.V.T @ np.ones(fs.n_beta)))
        nb = hp.get_all_neighbours(nside, np.arange(fs.n_beta))
        W = np.zeros((fs.n_beta, fs.n_beta))
        for i in range(fs.n_beta):
            for q in nb[:, i]:
                if q >= 0 and q != i:
                    W[i, q] = W[q, i] = 1.0
        G = np.diag(W.sum(1)) - W
        resid = np.max(np.abs(G @ fs.V - fs.V * (3 * fs.h**2 * fs.lam)[None, :]))
        check(f"graph Laplacian nside {nside}: K={fs.K}=Nb-1, orthonormal, V^T 1 = 0, G V = lam V",
              fs.K == fs.n_beta - 1 and orth < 1e-12 and one < 1e-12 and resid < 1e-10
              and fs.lam.min() > 0, f"orth {orth:.1e} one {one:.1e} resid {resid:.1e} "
              f"lam/(l(l+1)) at l=1 ~ {fs.lam[:3].mean() / 2:.3f}")
    fs = field_spec(1, "harmonic", lmax=1)
    xyz = np.column_stack(hp.pix2vec(1, np.arange(12)))
    proj = fs.V @ (fs.V.T @ xyz)
    check("harmonic nside 1 lmax 1: K=3, orthonormal, spans pixel-centre (x,y,z), lam = 2",
          fs.K == 3 and np.max(np.abs(fs.V.T @ fs.V - np.eye(3))) < 1e-12
          and np.max(np.abs(proj - xyz)) < 1e-12 and np.all(fs.lam == 2.0),
          f"span residual {np.max(np.abs(proj - xyz)):.1e}")
    fs4 = field_spec(4, "harmonic")
    check("harmonic nside 4 default lmax 8: K = 80 < 192, orthonormal, ell ordered",
          fs4.K == 80 and np.max(np.abs(fs4.V.T @ fs4.V - np.eye(80))) < 1e-12
          and np.all(np.diff(fs4.ell) >= 0), f"K={fs4.K}, ell[-1]={fs4.ell[-1]}")
    # constant mode projected out (SPEC decision): |1^T V| / sqrt(Nb) < 1e-12 at nside 1, 2, 4;
    # the raw basis (remove_constant=False) must show the known leakage so the test has teeth
    for nside in (1, 2, 4):
        fsh = field_spec(nside, "harmonic")
        one = np.linalg.norm(fsh.V.T @ np.ones(fsh.n_beta)) / math.sqrt(fsh.n_beta)
        fsr = field_spec(nside, "harmonic", remove_constant=False)
        one_raw = np.linalg.norm(fsr.V.T @ np.ones(fsr.n_beta)) / math.sqrt(fsr.n_beta)
        orth = np.max(np.abs(fsh.V.T @ fsh.V - np.eye(fsh.K)))
        same_span_ell1 = np.max(np.abs(fsh.V[:, :3] @ fsh.V[:, :3].T
                                       - fsr.V[:, :3] @ fsr.V[:, :3].T))
        check(f"harmonic nside {nside}: |1^T V|/sqrt(Nb) < 1e-12 with constant removed",
              one < 1e-12 and orth < 1e-12 and one_raw > 1e-3 and fsh.K == fsr.K,
              f"removed {one:.1e}, raw {one_raw:.3f}, orth {orth:.1e}, K={fsh.K}, "
              f"l=1 projector diff vs raw {same_span_ell1:.1e}")
    # rho transform round trip and Jacobian by finite differences
    hyp = default_hyper(fs)
    uu = np.linspace(-6, 6, 13)
    rt = u_from_rho(rho_from_u(uu, hyp.rho_min, hyp.rho_max), hyp.rho_min, hyp.rho_max)
    eps = 1e-6
    fd = (np.log(rho_from_u(uu + eps, hyp.rho_min, hyp.rho_max))
          - np.log(rho_from_u(uu - eps, hyp.rho_min, hyp.rho_max))) / (2 * eps)
    check("rho transform round trip and d log rho / du",
          np.max(np.abs(rt - uu)) < 1e-9
          and np.max(np.abs(np.log(fd) - log_dlogrho_du(uu, hyp.rho_min, hyp.rho_max))) < 1e-7,
          f"round trip {np.max(np.abs(rt - uu)):.1e}")
    # spectrum normalization: pixel-averaged prior variance equals sigma^2
    s, c = matern_spectrum(fs4, 0.3)
    Kmat = fs4.V @ np.diag((s / c) ** 2) @ fs4.V.T
    check("c(rho) normalization: trace(K)/Nb = 1", abs(np.trace(Kmat) / fs4.n_beta - 1) < 1e-12,
          f"{np.trace(Kmat) / fs4.n_beta:.12f}")

    # --- (c) ESS sanity
    n = 20000
    iid = rng.standard_normal(n)
    ar = np.empty(n)
    ar[0] = 0.0
    for i in range(1, n):
        ar[i] = 0.9 * ar[i - 1] + rng.standard_normal()
    e_iid, e_ar = ess_geyer(iid), ess_geyer(ar)
    check("ess_geyer: iid ~ n, AR(1) phi=0.9 ~ n/19",
          0.7 * n < e_iid < 1.4 * n and 0.6 * n / 19 < e_ar < 1.6 * n / 19,
          f"iid {e_iid:.0f} (n={n}), AR(1) {e_ar:.0f} (expected {n / 19:.0f})")

    # --- (iii) the 2-D grid conditional draw is exact: histogram vs grid probabilities
    fs1 = field_spec(1, "graph_laplacian")
    hyp1 = default_hyper(fs1)
    s, c = matern_spectrum(fs1, 0.8)
    y_test = 0.2 * (s / c) * rng.standard_normal(fs1.K)
    grid = hyper_conditional_grid(y_test, fs1, hyp1)
    n_draw = 20000
    ls_d, u_d, i_d, j_d = draw_from_grid(grid, rng, n_draw)
    # within-cell jitter must map back to the same cell
    i_back = np.floor((ls_d - (grid.ls[0] - 0.5 * grid.dls)) / grid.dls).astype(int)
    j_back = np.floor((u_d - (grid.u[0] - 0.5 * grid.du)) / grid.du).astype(int)
    cells_ok = np.all(i_back == i_d) and np.all(j_back == j_d)
    # aggregate fine cells into blocks of 16 x 16, merge blocks with expectation < 5
    nb_i, nb_j = grid.prob.shape[0] // 16, grid.prob.shape[1] // 16
    expected = grid.prob.reshape(nb_i, 16, nb_j, 16).sum(axis=(1, 3)) * n_draw
    counts = np.zeros_like(expected)
    np.add.at(counts, (i_d // 16, j_d // 16), 1)
    big = expected >= 5.0
    e_bins = np.concatenate([expected[big], [expected[~big].sum()]])
    c_bins = np.concatenate([counts[big], [counts[~big].sum()]])
    keep = e_bins > 0
    chi2 = float(np.sum((c_bins[keep] - e_bins[keep]) ** 2 / e_bins[keep]))
    dof = int(keep.sum()) - 1
    pval = float(stats.chi2.sf(chi2, dof))
    check("(iii) grid draw exact: chi-square of 20000 draws vs grid probabilities (p > 0.01)",
          cells_ok and pval > 0.01,
          f"chi2 {chi2:.1f} dof {dof} p {pval:.3f}; bins {dof + 1}; mean log sigma draws "
          f"{ls_d.mean():.4f} vs grid {np.sum(grid.prob.sum(1) * grid.ls):.4f}; cells "
          f"consistent {cells_ok}")

    # --- (d) closed form: Woodbury path vs direct N B-dimensional Gaussian at random theta
    fs3 = field_spec(1, "harmonic", lmax=1)
    hyp3 = default_hyper(fs3)
    data3, truth3 = simulate_linear(fs3, hyp3, 4, 2, rng, mu=-3.0, sigma=0.15, rho=1.2)
    Nd, Bd = data3.d.shape
    for trial in range(3):
        mu_t = -3.0 + 0.01 * rng.standard_normal()
        ls_t = math.log(0.2) + 0.5 * rng.standard_normal()
        u_t = 2.0 * rng.standard_normal()
        cf1 = linear_gaussian_closed_form(data3.d, data3.F, data3.offset, data3.inv_var,
                                          data3.region, fs3, hyp3, np.array([mu_t]),
                                          np.array([ls_t]), np.array([u_t]))
        # direct
        s, c = matern_spectrum(fs3, rho_from_u(u_t, hyp3.rho_min, hyp3.rho_max))
        Kmat = fs3.V @ np.diag((s / c) ** 2) @ fs3.V.T
        X = np.zeros((Nd * Bd, fs3.n_beta))
        X[np.arange(Nd * Bd), np.repeat(data3.region, Bd)] = data3.F.ravel()
        Sigma = math.exp(2 * ls_t) * X @ Kmat @ X.T + np.diag(1.0 / data3.inv_var.ravel())
        r = (data3.d - data3.offset).ravel() - mu_t * X @ np.ones(fs3.n_beta)
        logdet_sigma = 2.0 * float(np.sum(np.log(np.diag(np.linalg.cholesky(Sigma)))))
        direct = (-0.5 * r @ np.linalg.solve(Sigma, r) - 0.5 * logdet_sigma
                  - 0.5 * Nd * Bd * math.log(2 * math.pi)
                  + float(log_hyperprior(mu_t, ls_t, u_t, hyp3)))
        # beta | theta direct in y space
        U = X @ fs3.V
        P = U.T @ (data3.inv_var.ravel()[:, None] * U) + np.diag((c / s) ** 2) * math.exp(-2 * ls_t)
        ymean = np.linalg.solve(P, U.T @ (data3.inv_var.ravel() * r))
        bmean_direct = mu_t + fs3.V @ ymean
        bmean, bcov = cf1.beta_given_theta(mu_t, ls_t, u_t)
        check(f"closed form Woodbury = direct (trial {trial})",
              abs(cf1.logpost[0, 0, 0] - direct) < 1e-7 * abs(direct) + 1e-8
              and np.max(np.abs(bmean - bmean_direct)) < 1e-10
              and np.max(np.abs(bcov - fs3.V @ np.linalg.inv(P) @ fs3.V.T)) < 1e-12,
              f"logpost {cf1.logpost[0, 0, 0]:.8f} direct {direct:.8f}")

    # --- (ii) centred Gibbs reference vs closed form, tiny linear-Gaussian problem, K = 3
    info = np.bincount(np.repeat(data3.region, Bd), weights=(data3.F**2 * data3.inv_var).ravel(),
                       minlength=fs3.n_beta)
    bhat = np.bincount(np.repeat(data3.region, Bd), weights=(data3.F * (data3.d - data3.offset)
                                                             * data3.inv_var).ravel(),
                       minlength=fs3.n_beta) / info
    mu_hat = bhat.mean()
    sd_hat = math.sqrt(1.0 / info.sum())
    mu_grid = np.linspace(mu_hat - 10 * sd_hat, mu_hat + 10 * sd_hat, 61)
    ls_grid = np.linspace(-5.5, 1.5, 85)
    u_grid = np.linspace(-11.0, 11.0, 89)
    t0 = time.time()
    cf = linear_gaussian_closed_form(data3.d, data3.F, data3.offset, data3.inv_var, data3.region,
                                     fs3, hyp3, mu_grid, ls_grid, u_grid)
    t_cf = time.time() - t0
    edge = max(cf.weights[[0, -1], :, :].sum(), cf.weights[:, [0, -1], :].sum(),
               cf.weights[:, :, [0, -1]].sum())
    print(f"  closed form: {t_cf:.1f} s, edge mass {edge:.1e}, theta mean "
          f"{np.round(cf.theta_mean, 4)}, sd {np.round(cf.theta_sd, 4)}, truth beta mean "
          f"{truth3['beta'].mean():.4f}")
    n_it = 20000
    res = centred_gibbs_reference(data3, fs3, hyp3, "linear", rng, n_it, n_burn=6000,
                                  n_mh_sub=5, init=dict(mu=mu_hat, log_sigma=math.log(0.2)))
    names = ["mu", "log_sigma", "u", "log_rho"]
    chains = [res["mu"], res["log_sigma"], res["u"], res["log_rho"]]
    worst = 0.0
    lines = []
    for nm, ch, ref_m, ref_s in zip(names, chains, cf.theta_mean, cf.theta_sd):
        ess = ess_geyer(ch)
        se = ch.std() / math.sqrt(ess)
        zscore = (ch.mean() - ref_m) / se
        worst = max(worst, abs(zscore))
        lines.append(f"{nm}: gibbs {ch.mean():.4f} cf {ref_m:.4f} (sd {ref_s:.4f}) ESS "
                     f"{ess:.0f} z {zscore:+.2f}")
    zb = []
    for r in range(fs3.n_beta):
        ch = res["beta"][:, r]
        se = ch.std() / math.sqrt(res["ess"]["beta"][r])
        zb.append((ch.mean() - cf.beta_mean[r]) / se)
    zb = np.array(zb)
    worst = max(worst, np.max(np.abs(zb)))
    sd_ratio = res["beta"].std(axis=0) / cf.beta_sd
    check("(ii) centred Gibbs = closed form (all means within 4 SE)", worst < 4.0,
          "; ".join(lines) + f"; beta pixels max |z| {np.max(np.abs(zb)):.2f}, beta ESS min "
          f"{res['ess']['beta'].min():.0f}, sd ratio range {sd_ratio.min():.3f} to "
          f"{sd_ratio.max():.3f}; accept {res['accept_rate']:.3f}; {res['seconds']:.0f} s; "
          f"{n_it} iterations")

    print(f"\n{len(fails)} failure(s)" + (": " + ", ".join(fails) if fails else ""))
    return len(fails) == 0


if __name__ == "__main__":
    sys.exit(0 if _selftest() else 1)
