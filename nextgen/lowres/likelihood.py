"""Collapsed likelihood for one synchrotron index per sky region with dense noise covariance.

Data model, per band b (thermodynamic units), elements p = 1..n (Stokes Q then U, any pixel order):

    d_b = s_b(beta_r(p)) a_s + f_b a_d + a_c + n_b,      n_b ~ N(0, N_b),

with a_s, a_d, a_c the synchrotron, dust and CMB amplitude vectors (n each), s_b the band-averaged
power law of index beta at the reference frequency nu_ref, f_b the band-averaged dust spectrum, and
the CMB column equal to one. N_b^-1 is dense (n x n) or diagonal (n). One index per region.

Amplitude priors: flat on dust; flat or N(0, cmb_prior_sd^2) per element on the CMB; flat or
N(0, tau_r^2) per element on synchrotron (tau per region). With x = (a_s, a_d, a_c),

    A(beta) = sum_b M_b^T N_b^-1 M_b + prior precisions,     rhs(beta) = sum_b M_b^T N_b^-1 d_b,
    chi2_min(beta) = sum_b d_b^T N_b^-1 d_b - rhs^T A^-1 rhs.

Three objectives (to be minimized), all exact:

    profile        chi2_min / 2
    marginal       chi2_min / 2 + (1/2) ln det A [+ sum_p ln tau_r(p) with a synchrotron prior]

With flat synchrotron priors `marginal` is the flat-prior marginal, which depends on nu_ref. With the
zero-mean prior it is the proper marginal. With flat priors on every amplitude the profile score is
unbiased for any fixed true amplitudes and Gaussian noise of known covariance; the covariance of the
optimum is the sandwich J^-1 K J^-1 with

    J_rr' = (M_r a)^T P (M_r' a),     K = J + T,     T_rr' = tr(A^-1 M_r^T P M_r'),

P = N^-1 - N^-1 M A^-1 M^T N^-1 and M_r = dM / dbeta_r (derivation in the E039 record).

With a Gaussian CMB prior the same statements hold over CMB skies drawn from that prior, with the CMB
sample variance counted in K. They do not hold conditional on one sky. If the true CMB variance is
kappa^2 times the prior variance, the expected score is (kappa^2 - 1) t with

    t_r = tr(Lam A^-1 M_r^T B^T),     B = A^-1 M^T N^-1,   Lam the prior precision,

and the optimum is biased by -(kappa^2 - 1) J^-1 t to first order (2026-11-10 stats review, checked by
Monte Carlo). `information` returns t in `sol["t_cmb"]` so that the sensitivity can be quoted.
"""
import os

import numpy as np
import scipy.linalg as sl

# Which triangle the Cholesky factor of the (C-ordered) normal matrix uses. On macOS 27.0 (build 26A5388g)
# Accelerate's dpotrf with uplo = "U" writes up to 768 kB past an internal buffer for orders of about 1500 to
# 6000 (E043). scipy maps a C-ordered array with lower=True to uplo = "U", so the factor is taken with
# lower=False (uplo = "L"), which does not overflow, costs the same and gives the same numbers. The
# environment variable restores the old path, for reproducing the overflow under the canary only.
CHO_LOWER = os.environ.get("LOWRES_CHOLESKY_UPLO_U", "0") == "1"


class BandData:
    """One band: data vector d (n,), inverse covariance (n, n) dense or (n,) diagonal, and its spectra:
    `sync(beta) -> (s, ds/dbeta)` for an array of beta, scalars `dust` and `cmb`."""

    def __init__(self, name, d, ninv, sync, dust, cmb=1.0):
        self.name, self.d, self.ninv, self.sync, self.dust, self.cmb = name, np.asarray(d, float), ninv, sync, dust, cmb
        self.dense = (np.ndim(ninv) == 2)

    def apply(self, v):
        return self.ninv @ v if self.dense else self.ninv * v

    def scaled(self, left, right):
        """diag(left) N^-1 diag(right) as a dense (n, n) array."""
        if self.dense:
            return (left[:, None] * self.ninv) * right[None, :]
        return np.diag(left * self.ninv * right)


class RegionLikelihood:
    def __init__(self, bands, region, use_dust=True, use_cmb=True, cmb_prior_sd=None, sync_prior=False,
                 cmb_prior_prec=None):
        """cmb_prior_sd: white Gaussian prior per element. cmb_prior_prec: dense prior precision (n, n) for the
        CMB amplitudes, for example the inverse of an LCDM pixel covariance; give one or the other."""
        assert not (cmb_prior_sd and cmb_prior_prec is not None)
        self.bands = bands
        self.region = np.asarray(region, dtype=int)
        self.n = self.region.size
        self.nreg = int(self.region.max()) + 1
        self.use_dust, self.use_cmb, self.sync_prior = use_dust, use_cmb, sync_prior
        self.cmb_prior_sd = cmb_prior_sd
        n = self.n
        self.w = [b.apply(b.d) for b in bands]                              # N_b^-1 d_b
        self.dNd = float(sum(b.d @ w for b, w in zip(bands, self.w)))
        self.ncomp = 1 + int(use_dust) + int(use_cmb)
        self.id, self.ic = (1 if use_dust else None), (1 + int(use_dust) if use_cmb else None)
        # constant blocks of A and rhs
        one = np.ones(n)
        self.A_const = np.zeros((self.ncomp * n, self.ncomp * n))
        self.rhs_const = np.zeros(self.ncomp * n)
        cols = self._nuisance_columns()
        for (i, ci) in cols:
            self.rhs_const[i * n:(i + 1) * n] = sum(getattr(b, ci) * w for b, w in zip(bands, self.w))
            for (j, cj) in cols:
                blk = sum(b.scaled(one * getattr(b, ci), one * getattr(b, cj)) for b in bands)
                self.A_const[i * n:(i + 1) * n, j * n:(j + 1) * n] = blk
        if use_cmb and cmb_prior_sd:
            k = self.ic
            self.A_const[k * n:(k + 1) * n, k * n:(k + 1) * n][np.diag_indices(n)] += 1.0 / cmb_prior_sd ** 2
        if use_cmb and cmb_prior_prec is not None:
            k = self.ic
            self.A_const[k * n:(k + 1) * n, k * n:(k + 1) * n] += cmb_prior_prec

    def _nuisance_columns(self):
        cols = []
        if self.use_dust:
            cols.append((self.id, "dust"))
        if self.use_cmb:
            cols.append((self.ic, "cmb"))
        return cols

    # ---- assembly -------------------------------------------------------------------------------
    def sync_columns(self, beta):
        """s (B, n) and ds/dbeta (B, n) for region indices beta (nreg,)."""
        s, ds = np.empty((len(self.bands), self.n)), np.empty((len(self.bands), self.n))
        for i, b in enumerate(self.bands):
            sr, dsr = b.sync(beta)
            s[i], ds[i] = sr[self.region], dsr[self.region]
        return s, ds

    def assemble(self, beta, log_tau=None):
        n = self.n
        s, ds = self.sync_columns(beta)
        A = self.A_const.copy()
        rhs = self.rhs_const.copy()
        A[:n, :n] = sum(b.scaled(s[i], s[i]) for i, b in enumerate(self.bands))
        rhs[:n] = sum(s[i] * w for i, w in enumerate(self.w))
        one = np.ones(n)
        for (j, cj) in self._nuisance_columns():
            blk = sum(b.scaled(s[i], one * getattr(b, cj)) for i, b in enumerate(self.bands))
            A[:n, j * n:(j + 1) * n] = blk
            A[j * n:(j + 1) * n, :n] = blk.T
        if self.sync_prior:
            A[:n, :n][np.diag_indices(n)] += np.exp(-2.0 * np.asarray(log_tau))[self.region]
        return A, rhs, s, ds

    def solve(self, beta, log_tau=None, need_inverse=False):
        A, rhs, s, ds = self.assemble(beta, log_tau)
        c = sl.cho_factor(A, lower=CHO_LOWER, overwrite_a=True, check_finite=False)
        x = sl.cho_solve(c, rhs, check_finite=False)
        logdet = 2.0 * float(np.sum(np.log(np.diag(c[0]))))
        chi2 = self.dNd - float(rhs @ x)
        out = {"x": x, "chi2": chi2, "logdet": logdet, "s": s, "ds": ds, "chol": c}
        if need_inverse:
            # the identity is made in Fortran order and overwritten, so that the inverse costs one matrix of memory, not two
            out["Ainv"] = sl.cho_solve(c, np.eye(A.shape[0], order="F"), overwrite_b=True, check_finite=False)
        return out

    # ---- objectives and gradients ---------------------------------------------------------------
    def residuals(self, sol):
        n = self.n
        x = sol["x"]
        a_s = x[:n]
        res = []
        for i, b in enumerate(self.bands):
            m = sol["s"][i] * a_s
            if self.use_dust:
                m = m + b.dust * x[self.id * n:(self.id + 1) * n]
            if self.use_cmb:
                m = m + b.cmb * x[self.ic * n:(self.ic + 1) * n]
            res.append(b.d - m)
        return res

    def _grad_chi2_half(self, sol):
        """d(chi2_min / 2)/dbeta_r = -sum_b (N_b^-1 r_b)^T (ds_b a_s) over the elements of region r."""
        a_s = sol["x"][:self.n]
        g = np.zeros(self.n)
        for i, (b, r) in enumerate(zip(self.bands, self.residuals(sol))):
            g -= b.apply(r) * sol["ds"][i] * a_s
        return np.bincount(self.region, weights=g, minlength=self.nreg)

    def _grad_logdet_half(self, sol):
        """d((1/2) ln det A)/dbeta_r = sum_b sum_{p in r} ds_b(p) [A^-1_{s.} M_b^T N_b^-1]_{pp}."""
        n, Ai = self.n, sol["Ainv"]
        g = np.zeros(n)
        for i, b in enumerate(self.bands):
            Y = Ai[:n, :n] * sol["s"][i][None, :]
            if self.use_dust:
                Y = Y + b.dust * Ai[:n, self.id * n:(self.id + 1) * n]
            if self.use_cmb:
                Y = Y + b.cmb * Ai[:n, self.ic * n:(self.ic + 1) * n]
            t = (Y * b.ninv.T).sum(1) if b.dense else np.diag(Y) * b.ninv
            g += sol["ds"][i] * t
        return np.bincount(self.region, weights=g, minlength=self.nreg)

    def profile(self, beta, grad=False):
        sol = self.solve(beta)
        f = 0.5 * sol["chi2"]
        return (f, self._grad_chi2_half(sol)) if grad else f

    def marginal(self, beta, log_tau=None, grad=False):
        """Flat-prior marginal (sync_prior False) or zero-mean-prior marginal (sync_prior True).
        With grad: gradient in beta, and in log_tau when the prior is on."""
        sol = self.solve(beta, log_tau, need_inverse=grad)
        f = 0.5 * sol["chi2"] + 0.5 * sol["logdet"]
        if self.sync_prior:
            f += float(np.sum(np.asarray(log_tau)[self.region]))
        if not grad:
            return f
        g = self._grad_chi2_half(sol) + self._grad_logdet_half(sol)
        if not self.sync_prior:
            return f, g
        # d/dlog_tau_r: the prior term a_s^2 / tau^2 inside chi2 and tau^-2 on the diagonal of A
        n = self.n
        a_s = sol["x"][:n]
        e = np.exp(-2.0 * np.asarray(log_tau))[self.region]
        gt = -a_s ** 2 * e - np.diag(sol["Ainv"])[:n] * e + 1.0
        return f, g, np.bincount(self.region, weights=gt, minlength=self.nreg)

    # ---- information and sandwich ---------------------------------------------------------------
    def information(self, beta, log_tau=None):
        """J (profiled Fisher information at the fitted amplitudes), T (noise-noise term), and the
        fitted amplitudes. K = J + T; the estimator covariance is J^-1 K J^-1."""
        n = self.n
        sol = self.solve(beta, log_tau, need_inverse=True)
        Ai, a_s = sol["Ainv"], sol["x"][:n]
        E = sum(b.scaled(sol["ds"][i], sol["ds"][i]) for i, b in enumerate(self.bands))
        C = np.zeros((n, self.ncomp * n))
        C[:, :n] = sum(b.scaled(sol["ds"][i], sol["s"][i]) for i, b in enumerate(self.bands))
        one = np.ones(n)
        for (j, cj) in self._nuisance_columns():
            C[:, j * n:(j + 1) * n] = sum(b.scaled(sol["ds"][i], one * getattr(b, cj)) for i, b in enumerate(self.bands))
        X = E - C @ Ai @ C.T
        R = np.zeros((self.nreg, n))
        R[self.region, np.arange(n)] = 1.0
        J = R @ (a_s[:, None] * X * a_s[None, :]) @ R.T
        T = R @ (Ai[:n, :n] * X) @ R.T
        if self.use_cmb and self.cmb_prior_sd:
            # t_r = sigma_c^-2 sum_{p in r} sum_k (C A^-1)[p, cmb k] (A^-1)[cmb k, sync p]
            k = self.ic
            CA = C @ Ai[:, k * n:(k + 1) * n]                                 # (n, n): rows sync p, columns cmb k
            sol["t_cmb"] = (R @ np.sum(CA * Ai[:n, k * n:(k + 1) * n], axis=1)) / self.cmb_prior_sd ** 2
        return J, T, sol
