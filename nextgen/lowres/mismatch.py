"""What a wrong covariance for the CMB amplitudes does to the profile estimate of the region indices, to first order.

Model as in `likelihood.py`: d = F a + n, F = [M_s(beta) M_d U] (synchrotron, dust, CMB columns stacked over bands),
A = F^T N^-1 F + Lam on the CMB block, a_hat = G d with G = A^-1 F^T N^-1, residual r = R d with R = 1 - F G. The
score of the profile objective is g_r = -(F_r a_hat)^T N^-1 r with F_r = dF / dbeta_r.

Let the CMB amplitudes have covariance S_true while the fit assumes S_w = Lam^-1, and write dS = S_true - S_w. With
eps = U c + n and Z = A^-1 F^T N^-1 U = E_c - A^-1[:, cmb] Lam (so that G U = Z):

    E[g_r]        = -tr( dS  (F_r Z)^T  N^-1 R U )                         (zero when the prior is the truth)
    Cov(g_r, g_s) = J_rs + T_rs + q_r^T dS q_s ,   q_r = (N^-1 R U)^T F_r a   (first order in the amplitude errors)

so the profile optimum is biased by -J^-1 E[g] and has covariance J^-1 (J + T + q^T dS q) J^-1. The first is exact
to the order of the sandwich; the second neglects the change of the noise-noise term T under the mismatch.
N^-1 R U is formed band by band as Q_b = N_b^-1 (cmb_b 1 - diag(s_b) Z_s - dust_b Z_d - cmb_b Z_c).

Suggested by the second referee of Paper A1 (2026-11-11) in place of a simulation at another resolution.
"""
import numpy as np

from .fit import sandwich_plugin


def cmb_prior_mismatch(like, beta, lam, d_cov):
    """First-order bias and covariance of the profile optimum when the CMB covariance differs from the prior.

    like: RegionLikelihood built with the prior of the fit. beta: where to evaluate (the optimum).
    lam: the prior precision on the CMB block, (n, n) or a scalar (white prior 1 / sd^2).
    d_cov: S_true - S_prior, (n, n).
    Returns dict: bias (nreg,), cov_matched (nreg, nreg) (the debiased sandwich), cov_true (nreg, nreg),
    eg (nreg,) the expected score, dK (nreg, nreg)."""
    assert like.use_cmb
    n, nreg = like.n, like.nreg
    cov_c, cov_s, Jp, T, sol, ok = sandwich_plugin(like, beta)
    Ai = sol["Ainv"]
    a_s = sol["x"][:n]
    k = like.ic
    Ac = Ai[:, k * n:(k + 1) * n]                                   # (ncomp n, n)
    Zfull = -(Ac / np.asarray(lam) if np.ndim(lam) == 0 else Ac @ lam)
    if np.ndim(lam) == 0:
        Zfull = -Ac * float(lam)
    Zfull[k * n:(k + 1) * n][np.diag_indices(n)] += 1.0
    Zs = Zfull[:n]
    Zd = Zfull[like.id * n:(like.id + 1) * n] if like.use_dust else None
    Zc = Zfull[k * n:(k + 1) * n]
    row = np.zeros(n)                                              # sum_b ds_b * rowsum((Q_b dS) * Z_s)
    V = np.zeros((n, nreg))
    Qq = np.zeros((n, nreg))
    for i, b in enumerate(like.bands):
        inner = b.cmb * np.eye(n) - sol["s"][i][:, None] * Zs - b.cmb * Zc
        if Zd is not None:
            inner -= b.dust * Zd
        Q = b.ninv @ inner if b.dense else b.ninv[:, None] * inner
        row += sol["ds"][i] * np.sum((Q @ d_cov) * Zs, axis=1)
        V[:] = 0.0
        V[np.arange(n), like.region] = sol["ds"][i] * a_s
        Qq += Q.T @ V
        del inner, Q
    eg = -np.bincount(like.region, weights=row, minlength=nreg)
    dK = Qq.T @ d_cov @ Qq
    Jc = Jp - T
    Ji = np.linalg.inv(Jc)
    return {"bias": -Ji @ eg, "cov_matched": cov_s, "cov_true": Ji @ (Jp + dK) @ Ji, "eg": eg, "dK": 0.5 * (dK + dK.T), "ok": ok}
