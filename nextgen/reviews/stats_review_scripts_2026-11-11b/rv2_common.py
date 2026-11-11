"""Reviewer's reference implementation, independent of nextgen/lowres/likelihood.py.

Stacked model: d = M_sd(beta) (a_s, a_d) + U a_c + n over the bands; with a Gaussian prior a_c ~ N(0, C) the CMB is
integrated into the noise, Sigma = N + U C U^T (dense across bands), and everything is ordinary generalized least
squares in (a_s, a_d):

    At = M_sd^T Sigma^-1 M_sd,   Wt = Sigma^-1 - Sigma^-1 M_sd At^-1 M_sd^T Sigma^-1,
    profile = d^T Wt d / 2,      marginal = profile + ln det At / 2 (+ a constant in beta),
    J_rr' = (D_r a_s)^T Wt (D_r' a_s),   T_rr' = tr([At^-1]_ss D_r^T Wt D_r'),   D_r = dM_sd/dbeta_r (synchrotron column).

Bands may also enter through a linear map L_b (rows = coefficients): data L_b d_b, covariance L_b N_b L_b^T.
Nothing here calls the code under review.
"""
import numpy as np

FREQS = np.array([23.0, 33.0, 41.0, 61.0, 353.0])
NU_REF = 23.0


def sync_fun(nu):
    lr = np.log(nu / NU_REF)

    def f(beta):
        beta = np.atleast_1d(np.asarray(beta, float))
        s = np.exp(beta * lr)
        return s, s * lr
    return f


def dust_coef(nu, beta_d=1.5):
    return float((nu / 353.0) ** beta_d)


class Reference:
    """Dense reference. covs: list of noise covariances (n, n) per band; maps: list of L_b (k_b, n) or None;
    cmb_cov: prior covariance (n, n) of the CMB amplitudes, or None for no CMB component."""

    def __init__(self, freqs, covs, region, cmb_cov=None, maps=None, use_dust=True):
        self.freqs, self.region = np.asarray(freqs, float), np.asarray(region, int)
        self.n, self.nreg, self.nb = self.region.size, int(self.region.max()) + 1, len(freqs)
        n = self.n
        self.maps = [np.eye(n) if (maps is None or maps[i] is None) else maps[i] for i in range(self.nb)]
        self.use_dust = use_dust
        self.Lstack = np.zeros((sum(m.shape[0] for m in self.maps), self.nb * n))
        r0 = 0
        self.rows = []
        for i, m in enumerate(self.maps):
            self.Lstack[r0:r0 + m.shape[0], i * n:(i + 1) * n] = m
            self.rows.append(slice(r0, r0 + m.shape[0]))
            r0 += m.shape[0]
        N = np.zeros((self.nb * n, self.nb * n))
        for i, c in enumerate(covs):
            N[i * n:(i + 1) * n, i * n:(i + 1) * n] = c
        U = np.tile(np.eye(n), (self.nb, 1))
        Sig = N if cmb_cov is None else N + U @ cmb_cov @ U.T
        self.Sigma = self.Lstack @ Sig @ self.Lstack.T                         # covariance of the stacked (mapped) data
        self.Sigma = 0.5 * (self.Sigma + self.Sigma.T)
        self.Si = np.linalg.inv(self.Sigma)
        self.Si = 0.5 * (self.Si + self.Si.T)
        self.sync = [sync_fun(f) for f in self.freqs]
        self.fd = [dust_coef(f) for f in self.freqs]

    def design(self, beta):
        """M_sd (rows, 2n or n) and the per-element derivative columns Dm (rows, n): d(model)/d(beta at element p) per unit a_s(p)."""
        n = self.n
        cols = 2 if self.use_dust else 1
        M = np.zeros((self.nb * n, cols * n))
        Dm = np.zeros((self.nb * n, n))
        for i in range(self.nb):
            s, ds = self.sync[i](beta)
            M[i * n:(i + 1) * n, :n] = np.diag(s[self.region])
            Dm[i * n:(i + 1) * n, :] = np.diag(ds[self.region])
            if self.use_dust:
                M[i * n:(i + 1) * n, n:] = self.fd[i] * np.eye(n)
        return self.Lstack @ M, self.Lstack @ Dm

    def parts(self, beta):
        M, Dm = self.design(beta)
        SM = self.Si @ M
        At = M.T @ SM
        Ati = np.linalg.inv(At)
        Wt = self.Si - SM @ Ati @ SM.T
        return M, Dm, At, Ati, 0.5 * (Wt + Wt.T)

    def stack(self, dvecs):
        return np.concatenate([m @ d for m, d in zip(self.maps, dvecs)])

    def profile(self, beta, dvecs):
        M, Dm, At, Ati, Wt = self.parts(beta)
        d = self.stack(dvecs)
        return 0.5 * float(d @ Wt @ d)

    def marginal(self, beta, dvecs):
        M, Dm, At, Ati, Wt = self.parts(beta)
        d = self.stack(dvecs)
        return 0.5 * float(d @ Wt @ d) + 0.5 * np.linalg.slogdet(At)[1]

    def info(self, beta, a_s):
        """J(a_s), T, and the per-element matrix X = Dm^T Wt Dm."""
        M, Dm, At, Ati, Wt = self.parts(beta)
        X = Dm.T @ Wt @ Dm
        R = np.zeros((self.nreg, self.n))
        R[self.region, np.arange(self.n)] = 1.0
        J = R @ (a_s[:, None] * X * a_s[None, :]) @ R.T
        T = R @ (Ati[:self.n, :self.n] * X) @ R.T
        return J, T

    def score_moments(self, beta, a_s, a_d=None):
        """Exact mean and covariance of the profile score at beta when the data are N(M_sd(beta) a, Sigma):
        g_r = -d^T B_r d with B_r = Wt D_r [At^-1]_{s.} M^T Sigma^-1 (D_r the region-r derivative columns)."""
        M, Dm, At, Ati, Wt = self.parts(beta)
        a = a_s if a_d is None else np.concatenate([a_s, a_d])
        m = M @ a
        G = Ati[:self.n] @ M.T @ self.Si                                      # a_hat_s = G d
        mean, lin, Bs = np.zeros(self.nreg), [], []
        for r in range(self.nreg):
            Dr = Dm * (self.region == r)[None, :]
            B = Wt @ Dr @ G
            Bsym = 0.5 * (B + B.T)
            Bs.append(Bsym)
            mean[r] = -(np.trace(B @ self.Sigma) + m @ B @ m)
        cov = np.zeros((self.nreg, self.nreg))
        for r in range(self.nreg):
            for q in range(self.nreg):
                cov[r, q] = 2 * np.trace(Bs[r] @ self.Sigma @ Bs[q] @ self.Sigma) + 4 * m @ Bs[r] @ self.Sigma @ Bs[q] @ m
        return mean, cov

    def expected_objective(self, beta, beta0, a_s, a_d=None):
        """E[profile(beta)] when the data are N(M_sd(beta0) a, Sigma)."""
        M0, _ = self.design(beta0)
        a = a_s if a_d is None else np.concatenate([a_s, a_d])
        m = M0 @ a
        _, _, _, _, Wt = self.parts(beta)
        return 0.5 * (m @ Wt @ m + np.trace(Wt @ self.Sigma))


def fd_grad(f, x, h=1e-5):
    x = np.asarray(x, float)
    g = np.zeros(x.size)
    for i in range(x.size):
        e = np.zeros(x.size)
        e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g


def fd_hess(f, x, h=1e-3):
    x = np.asarray(x, float)
    n = x.size
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei, ej = np.zeros(n), np.zeros(n)
            ei[i], ej[j] = h, h
            H[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
    return H


def random_spd(n, rng, cond=30.0, scale=1.0):
    q, _ = np.linalg.qr(rng.standard_normal((n, n)))
    w = scale * np.exp(rng.uniform(0, np.log(cond), n))
    return (q * w[None, :]) @ q.T
