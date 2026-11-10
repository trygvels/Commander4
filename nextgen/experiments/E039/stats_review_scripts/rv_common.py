"""Independent reference implementation for the stats review of nextgen/lowres (reviewer code).

Everything here is built from explicit stacked matrices (M, N, priors) and does not call the code under
review, except `banddata()` which wraps the same numbers into the BandData objects that code expects.
Spectra for the reference are computed here from scratch (delta bands), not through spectra.py.
"""
import sys

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, "/Users/trygvels/code/Commander4/nextgen")
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402,F401

HK = 0.04799243073  # h/k in K/GHz (CODATA 2018), typed independently of spectra.py
TCMB = 2.7255


def g_ref(nu):
    x = HK * np.asarray(nu, float) / TCMB
    return x ** 2 * np.exp(x) / (np.exp(x) - 1.0) ** 2


def sync_ref(nu, beta, nu_ref):
    """delta band: thermodynamic response to (nu/nu_ref)^beta, and d/dbeta."""
    beta = np.atleast_1d(np.asarray(beta, float))
    s = (nu / nu_ref) ** beta / g_ref(nu)
    return s, s * np.log(nu / nu_ref)


def dust_ref(nu, bd=1.55, Td=19.6, nu0=353.0):
    x, x0 = HK * nu / Td, HK * nu0 / Td
    return (nu / nu0) ** (bd + 1.0) * (np.exp(x0) - 1.0) / (np.exp(x) - 1.0) / g_ref(nu)


class Problem:
    """Small synthetic problem. Noise covariance per band is dense for bands in `dense`."""

    def __init__(self, n=48, nreg=3, nus=(22.8, 28.4, 33.0, 44.1, 70.4, 353.0), dense=(0, 2, 5), seed=0,
                 noise=1.0, noise353=0.2, nu_ref=23.0, beta_true=None, amp=8.0, amp_d=0.5, amp_c=0.3,
                 random_regions=True, corr=0.5, noise_slope=0.0):
        rng = np.random.default_rng(seed)
        self.n, self.nreg, self.nus, self.nu_ref = n, nreg, np.array(nus, float), nu_ref
        if random_regions:
            # unequal, shuffled partition (every region non-empty)
            reg = np.concatenate([np.arange(nreg), rng.integers(0, nreg, n - nreg)])
            self.region = rng.permutation(reg)
        else:
            self.region = np.arange(n) % nreg
        self.beta_true = np.full(nreg, -3.1) if beta_true is None else np.asarray(beta_true, float)
        self.a_s = amp * rng.standard_normal(n) * (1.0 + 3.0 * rng.random(n))
        self.a_d = amp_d * rng.standard_normal(n)
        self.a_c = amp_c * rng.standard_normal(n)
        self.cov, self.chol, self.ninv_dense, self.is_dense = [], [], [], []
        for i, nu in enumerate(self.nus):
            # noise_slope > 0 makes the noise fall with frequency so that several bands carry weight
            sig = noise353 * noise if nu > 300 else noise * rng.uniform(1.0, 2.5) * (nu / 23.0) ** (-noise_slope)
            if i in dense:
                X = rng.standard_normal((n, n))
                C = sig ** 2 * (np.eye(n) + corr * X @ X.T / n)
            else:
                C = np.diag(sig ** 2 * rng.uniform(0.6, 1.6, n))
            self.cov.append(C)
            self.chol.append(np.linalg.cholesky(C))
            self.ninv_dense.append(np.linalg.inv(C))
            self.is_dense.append(i in dense)
        self.dust = np.array([dust_ref(nu) for nu in self.nus])
        self.nb = len(self.nus)

    # ---- reference matrices -------------------------------------------------------------------
    def M(self, beta, use_dust=True, use_cmb=True, nu_ref=None):
        nu_ref = self.nu_ref if nu_ref is None else nu_ref
        n, blocks = self.n, []
        for b, nu in enumerate(self.nus):
            s, _ = sync_ref(nu, beta, nu_ref)
            row = [np.diag(s[self.region])]
            if use_dust:
                row.append(self.dust[b] * np.eye(n))
            if use_cmb:
                row.append(np.eye(n))
            blocks.append(np.hstack(row))
        return np.vstack(blocks)

    def dM(self, beta, r, ncomp, nu_ref=None):
        """dM/dbeta_r, stacked (B n, ncomp n)."""
        nu_ref = self.nu_ref if nu_ref is None else nu_ref
        n = self.n
        out = np.zeros((self.nb * n, ncomp * n))
        sel = (self.region == r).astype(float)
        for b, nu in enumerate(self.nus):
            _, ds = sync_ref(nu, beta, nu_ref)
            out[b * n:(b + 1) * n, :n] = np.diag(ds[self.region] * sel)
        return out

    def N(self):
        return sl.block_diag(*self.cov)

    def W(self):
        return sl.block_diag(*self.ninv_dense)

    def prior_precision(self, ncomp, use_dust=True, use_cmb=True, cmb_sd=None, log_tau=None):
        n = self.n
        lam = np.zeros(ncomp * n)
        if log_tau is not None:
            lam[:n] = np.exp(-2.0 * np.asarray(log_tau))[self.region]
        if use_cmb and cmb_sd:
            k = 1 + int(use_dust)
            lam[k * n:(k + 1) * n] = 1.0 / cmb_sd ** 2
        return lam

    def mean(self, beta=None, a_s=None, a_d=None, a_c=None, nu_ref=None):
        """Noise-free data per band (B, n), generated with all three components."""
        beta = self.beta_true if beta is None else beta
        a_s = self.a_s if a_s is None else a_s
        a_d = self.a_d if a_d is None else a_d
        a_c = self.a_c if a_c is None else a_c
        nu_ref = self.nu_ref if nu_ref is None else nu_ref
        out = []
        for b, nu in enumerate(self.nus):
            s, _ = sync_ref(nu, beta, nu_ref)
            out.append(s[self.region] * a_s + self.dust[b] * a_d + a_c)
        return np.array(out)

    def noise(self, rng):
        return np.array([self.chol[b] @ rng.standard_normal(self.n) for b in range(self.nb)])

    # ---- wrap for the code under review -------------------------------------------------------
    def banddata(self, d, nu_ref=None, sync="ref"):
        """List of BandData for data d (B, n). sync='ref': reviewer's delta-band spectrum;
        sync='code': Band.delta(...).power_law from spectra.py."""
        nu_ref = self.nu_ref if nu_ref is None else nu_ref
        out = []
        for b, nu in enumerate(self.nus):
            ninv = self.ninv_dense[b] if self.is_dense[b] else np.diag(self.ninv_dense[b]).copy()
            if sync == "ref":
                f = (lambda be, nu=nu: sync_ref(nu, be, nu_ref))
            else:
                from lowres.spectra import Band
                B = Band.delta(f"b{b}", nu)
                f = (lambda be, B=B: B.power_law(be, nu_ref))
            out.append(BandData(f"b{b}", d[b], ninv, f, float(self.dust[b])))
        return out


def ref_solve(prob, d, beta, use_dust=True, use_cmb=True, cmb_sd=None, log_tau=None, nu_ref=None):
    """Reference solution by QR of the whitened, prior-augmented least-squares problem (no normal equations,
    uses N not N^-1). Returns dict with chi2 (min of 2Q), logdetA, x."""
    n = prob.n
    ncomp = 1 + int(use_dust) + int(use_cmb)
    M = prob.M(beta, use_dust, use_cmb, nu_ref)
    lam = prob.prior_precision(ncomp, use_dust, use_cmb, cmb_sd, log_tau)
    Mw = np.vstack([sl.solve_triangular(prob.chol[b], M[b * n:(b + 1) * n], lower=True) for b in range(prob.nb)])
    dw = np.concatenate([sl.solve_triangular(prob.chol[b], d[b], lower=True) for b in range(prob.nb)])
    keep = lam > 0
    aug = np.diag(np.sqrt(lam))[keep]
    Maug = np.vstack([Mw, aug])
    daug = np.concatenate([dw, np.zeros(aug.shape[0])])
    Q, R = np.linalg.qr(Maug)
    x = sl.solve_triangular(R, Q.T @ daug)
    res = daug - Maug @ x
    return {"chi2": float(res @ res), "logdetA": 2.0 * float(np.sum(np.log(np.abs(np.diag(R))))), "x": x, "M": M,
            "lam": lam}


def ref_info(prob, beta, a_s, use_dust=True, use_cmb=True, cmb_sd=None, log_tau=None, nu_ref=None):
    """Reference J, T, and the matrices A^-1, P, B from explicit stacked algebra."""
    ncomp = 1 + int(use_dust) + int(use_cmb)
    M = prob.M(beta, use_dust, use_cmb, nu_ref)
    W = prob.W()
    lam = prob.prior_precision(ncomp, use_dust, use_cmb, cmb_sd, log_tau)
    A = M.T @ W @ M + np.diag(lam)
    Ai = np.linalg.inv(A)
    B = Ai @ M.T @ W
    P = W - W @ M @ B
    x = np.zeros(ncomp * prob.n)
    x[:prob.n] = a_s
    dMs = [prob.dM(beta, r, ncomp, nu_ref) for r in range(prob.nreg)]
    J = np.array([[(dMs[r] @ x) @ P @ (dMs[q] @ x) for q in range(prob.nreg)] for r in range(prob.nreg)])
    T = np.array([[np.trace(Ai @ dMs[r].T @ P @ dMs[q]) for q in range(prob.nreg)] for r in range(prob.nreg)])
    return {"J": J, "T": T, "Ai": Ai, "P": P, "B": B, "M": M, "W": W, "dM": dMs, "lam": lam, "A": A}


def fd_grad(f, x, h=1e-5):
    x = np.asarray(x, float)
    g = np.zeros_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g


def fd_hess(f, x, h=1e-3):
    x = np.asarray(x, float)
    k = x.size
    H = np.zeros((k, k))
    f0 = f(x)
    for i in range(k):
        ei = np.zeros(k); ei[i] = h
        H[i, i] = (f(x + ei) - 2 * f0 + f(x - ei)) / h ** 2
        for j in range(i + 1, k):
            ej = np.zeros(k); ej[j] = h
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h ** 2)
    return H
