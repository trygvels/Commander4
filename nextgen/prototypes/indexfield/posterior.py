"""Log posterior densities and gradients of the hierarchical index-field model (E012 v2).

Data enter only through a collapsed log-likelihood log L(beta) and its gradient g = dlog L/dbeta,
by default model.collapsed_loglike(beta, data, mode); any callable beta -> (logL, grad) can be
substituted (used by the self-test and the negative controls).

Hyperparameters theta = (mu, log sigma, u), rho = rho(u) through field.RhoTransform, with the
proper hyperpriors of E012 v2 (all Jacobians included):

    log N(mu; mu_0, s_mu^2)
    log halfN(sigma; sigma_scale) + log sigma                 (sampling log sigma)
    log-uniform rho on [rho_min, rho_max]: log |dlog rho / du| = log[L s(u)(1 - s(u))]
    or the PC prior of Fuglstad et al. 2019 (rho_prior = "pc"):
    p(rho) = lam_1 rho^-2 exp(-lam_1/rho),
    lam_1 = -log(0.05) rho_0 (so P(rho < rho_0) = 0.05), jointly with sigma ~ Exp(lam_2),
    lam_2 = -log(0.05) / sigma_0 (P(sigma > sigma_0) = 0.05), both mapped to the sampled
    coordinates: + log rho + log |dlog rho / du| for u, + log sigma for log sigma. The PC prior is
    truncated to [rho_min, rho_max] by the u map; its normalization constant is dropped.

Non-centred (NonCentredPosterior): x = (z in R^K, mu, log sigma, u), beta = mu 1 + sigma S(rho) z,

    log pi(x) = log L(beta) - 1/2 z^T z + hyperprior and Jacobian terms.

Centred (CentredPosterior): x = (y in R^K, mu, log sigma, u), beta = mu 1 + V y, so that the
component of beta - mu 1 outside range(V) (zero prior variance) is excluded by construction,

    log pi(x) = log L(beta) - 1/2 y^T diag((c/s)^2) y / sigma^2 - K log sigma - 1/2 log det K(rho)
                + hyperprior and Jacobian terms,

with log det K = sum_i 2 log(s_i/c) (field.MaternField.logdet_K). For y = sigma (s/c) z the two
densities differ exactly by the change-of-variables Jacobian -K log sigma - 1/2 log det K.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from .field import MaternField, RhoTransform

try:
    from .model import collapsed_loglike
except ImportError:                      # model.py written by another agent; optional here
    collapsed_loglike = None

PC_ALPHA = 0.05      # tail probabilities of the PC prior, P(rho < rho_0) = P(sigma > sigma_0)


@dataclass
class HyperPrior:
    """Hyperprior settings of E012 v2. rho_min, pc_rho0 and pc_sigma0 default to h, 2 h and 0.5
    once the pixel spacing h is known (resolved(h))."""

    mu0: float = -3.0
    s_mu: float = 0.5
    sigma_scale: float = 0.3
    rho_prior: str = "loguniform"
    rho_min: float | None = None
    rho_max: float = float(np.pi)
    pc_rho0: float | None = None
    pc_sigma0: float | None = None

    def resolved(self, h: float) -> "HyperPrior":
        if self.rho_prior not in ("loguniform", "pc"):
            raise ValueError(f"unknown rho_prior {self.rho_prior!r}")
        return replace(self, rho_min=h if self.rho_min is None else self.rho_min,
                       pc_rho0=2.0 * h if self.pc_rho0 is None else self.pc_rho0,
                       pc_sigma0=0.5 if self.pc_sigma0 is None else self.pc_sigma0)


class _FieldPosterior:
    """Shared machinery: likelihood call, hyperprior terms, packing conventions.

    x = (field block of length K, mu, log sigma, u); dim = K + 3.
    """

    centred: bool

    def __init__(self, field: MaternField, hyper: HyperPrior, data=None, mode: str = "gauss",
                 loglike=None):
        self.field = field
        self.basis = field.basis
        self.n_beta = field.n_beta
        self.K = field.K
        self.hyper = hyper.resolved(self.basis.h)
        self.transform = RhoTransform(self.hyper.rho_min, self.hyper.rho_max)
        self.mode = mode
        if loglike is None:
            if collapsed_loglike is None:
                raise ImportError("indexfield.model is not importable; pass loglike=")
            if data is None:
                raise ValueError("data (model.PixelData) is required when loglike is None")
            loglike = lambda beta: collapsed_loglike(beta, data, mode)   # noqa: E731
        self._loglike = loglike
        self.n_loglike_errors = 0
        if self.hyper.rho_prior == "pc":
            self.pc_lam1 = -np.log(PC_ALPHA) * self.hyper.pc_rho0
            self.pc_lam2 = -np.log(PC_ALPHA) / self.hyper.pc_sigma0

    @property
    def dim(self) -> int:
        return self.K + 3

    def loglike(self, beta: np.ndarray):
        """(log L, dlog L/dbeta). A Cholesky failure in the collapsed likelihood (amplitude normal
        matrix A, or A + T^-1, not positive definite, which happens only for absurd beta reached
        during warm-up step-size search) is reported as zero density, log L = -inf, and counted in
        n_loglike_errors so that it is never silent (SPEC "Zero-density evaluations"). Both
        exception types are caught: model._collapsed_core converts the LinAlgError of
        np.linalg.cholesky into a ValueError, and a user-supplied loglike may raise either."""
        try:
            logL, g = self._loglike(beta)
        except (ValueError, np.linalg.LinAlgError):
            self.n_loglike_errors += 1
            return -np.inf, np.zeros(self.n_beta)
        if not np.isfinite(logL):
            self.n_loglike_errors += 1
            return -np.inf, np.zeros(self.n_beta)
        return float(logL), np.asarray(g, dtype=float)

    def _split(self, x: np.ndarray):
        x = np.asarray(x, dtype=float)
        if x.shape != (self.dim,):
            raise ValueError(f"x must have shape ({self.dim},)")
        return x[:self.K], float(x[self.K]), float(x[self.K + 1]), float(x[self.K + 2])

    def _hyper_terms(self, mu: float, log_sigma: float, u: float):
        """Hyperprior log density in the sampled coordinates and its gradient
        (d/dmu, d/dlog sigma, d/du), Jacobians included."""
        hp_ = self.hyper
        sigma = np.exp(log_sigma)
        tr = self.transform
        # mu ~ N(mu0, s_mu^2)
        lp = -0.5 * ((mu - hp_.mu0) / hp_.s_mu) ** 2
        dmu = -(mu - hp_.mu0) / hp_.s_mu ** 2
        # u: Jacobian of the logit map, common to both rho priors
        lp += tr.log_jacobian(u)
        du = tr.dlog_jacobian_du(u)
        if hp_.rho_prior == "loguniform":
            # sigma ~ halfN(0, sigma_scale), + log sigma Jacobian
            lp += -0.5 * (sigma / hp_.sigma_scale) ** 2 + log_sigma
            dls = 1.0 - (sigma / hp_.sigma_scale) ** 2
            # log-uniform rho: constant density in log rho (its -log L constant is dropped)
        else:
            rho = tr.rho(u)
            # sigma ~ Exp(lam2), + log sigma
            lp += np.log(self.pc_lam2) - self.pc_lam2 * sigma + log_sigma
            dls = 1.0 - self.pc_lam2 * sigma
            # rho: log p(rho) = log lam1 - 2 log rho - lam1/rho, times drho/du = rho dlog rho/du
            lp += np.log(self.pc_lam1) - 2.0 * np.log(rho) - self.pc_lam1 / rho + np.log(rho)
            du += (self.pc_lam1 / rho - 1.0) * tr.dlogrho_du(u)
        return float(lp), float(dmu), float(dls), float(du)

    def theta(self, x: np.ndarray):
        """(mu, sigma, rho) of a state."""
        _, mu, ls, u = self._split(x)
        return mu, float(np.exp(ls)), float(self.transform.rho(u))

    def _random_theta(self, rng: np.random.Generator):
        hp_ = self.hyper
        mu = rng.normal(hp_.mu0, hp_.s_mu)
        if hp_.rho_prior == "loguniform":
            sigma = max(abs(rng.normal(0.0, hp_.sigma_scale)), 0.05 * hp_.sigma_scale)
            u = float(self.transform.u(np.exp(rng.uniform(0.1, 0.9) * self.transform.L
                                              + self.transform.log_rho_min)))
        else:
            sigma = max(rng.exponential(1.0 / self.pc_lam2), 0.01 * hp_.pc_sigma0)
            rho = self.pc_lam1 / rng.exponential(1.0)
            rho = np.clip(rho, 1.01 * hp_.rho_min, 0.99 * hp_.rho_max)
            u = float(np.clip(self.transform.u(rho), -4.0, 4.0))
        return float(mu), float(np.log(sigma)), u

    def pack(self, block: np.ndarray, mu: float, log_sigma: float, u: float) -> np.ndarray:
        return np.concatenate([np.asarray(block, dtype=float), [mu, log_sigma, u]])


class NonCentredPosterior(_FieldPosterior):
    """x = (z[K], mu, log_sigma, u), beta = mu 1 + sigma S(rho) z (see module docstring).

    Gradient (chain rule of E012 "Gradient"), with g = dlog L/dbeta and w = s/c:
        d/dz = sigma w * (V^T g) - z
        d/dmu = 1^T g + hyperprior
        d/dlog sigma = sigma (S z)^T g + hyperprior and Jacobian
        d/du = sigma z^T (dS/dlog rho)^T g dlog rho/du + hyperprior and Jacobian.
    """

    centred = False

    def unpack(self, x: np.ndarray) -> dict:
        z, mu, ls, u = self._split(x)
        sigma, rho = float(np.exp(ls)), float(self.transform.rho(u))
        beta = mu + sigma * self.field.S_times(rho, z)
        return dict(z=z, mu=mu, log_sigma=ls, sigma=sigma, u=u, rho=rho, beta=beta)

    def beta(self, x: np.ndarray) -> np.ndarray:
        return self.unpack(x)["beta"]

    def logp_and_grad(self, x: np.ndarray):
        z, mu, ls, u = self._split(x)
        sigma = np.exp(ls)
        rho = float(self.transform.rho(u))
        w, dw = self.field.weights(rho)
        V = self.basis.V
        Sz = V @ (w * z)
        beta = mu + sigma * Sz
        logL, g = self.loglike(beta)
        lp_h, dmu_h, dls_h, du_h = self._hyper_terms(mu, ls, u)
        lp = logL - 0.5 * float(z @ z) + lp_h
        if not np.isfinite(lp):
            return -np.inf, np.zeros(self.dim)
        Vtg = V.T @ g
        dz = sigma * w * Vtg - z
        dmu = float(g.sum()) + dmu_h
        dls = sigma * float(Sz @ g) + dls_h
        du = sigma * float((dw * z) @ Vtg) * float(self.transform.dlogrho_du(u)) + du_h
        return float(lp), self.pack(dz, dmu, dls, du)

    def random_init(self, rng: np.random.Generator) -> np.ndarray:
        mu, ls, u = self._random_theta(rng)
        return self.pack(rng.normal(size=self.K), mu, ls, u)


class CentredPosterior(_FieldPosterior):
    """x = (y[K], mu, log_sigma, u), beta = mu 1 + V y (see module docstring).

    The prior N(beta; mu 1, sigma^2 K(rho)) is supported on mu 1 + range(V) only; parameterizing
    by y in R^K enforces this (the orthogonal complement has zero prior variance and is never
    sampled). Gradient, with w = s/c, dw = dw/dlog rho, Q = sum_i y_i^2 / w_i^2:
        d/dy = V^T g - y / (w^2 sigma^2)
        d/dmu = 1^T g + hyperprior
        d/dlog sigma = Q / sigma^2 - K + hyperprior and Jacobian
        d/du = [sum_i y_i^2 dw_i / (w_i^3 sigma^2) - sum_i dw_i / w_i] dlog rho/du + hyperprior.
    """

    centred = True

    def unpack(self, x: np.ndarray) -> dict:
        y, mu, ls, u = self._split(x)
        sigma, rho = float(np.exp(ls)), float(self.transform.rho(u))
        beta = mu + self.basis.V @ y
        return dict(y=y, mu=mu, log_sigma=ls, sigma=sigma, u=u, rho=rho, beta=beta)

    def beta(self, x: np.ndarray) -> np.ndarray:
        return self.unpack(x)["beta"]

    def logp_and_grad(self, x: np.ndarray):
        y, mu, ls, u = self._split(x)
        sigma = np.exp(ls)
        rho = float(self.transform.rho(u))
        w, dw = self.field.weights(rho)
        V = self.basis.V
        beta = mu + V @ y
        logL, g = self.loglike(beta)
        lp_h, dmu_h, dls_h, du_h = self._hyper_terms(mu, ls, u)
        yw2 = (y / w) ** 2
        Q = float(yw2.sum())
        logdet_K = float(2.0 * np.sum(np.log(w)))
        lp = logL - 0.5 * Q / sigma ** 2 - self.K * ls - 0.5 * logdet_K + lp_h
        if not np.isfinite(lp):
            return -np.inf, np.zeros(self.dim)
        Vtg = V.T @ g
        dy = Vtg - y / (w * w * sigma ** 2)
        dmu = float(g.sum()) + dmu_h
        dls = Q / sigma ** 2 - self.K + dls_h
        dlogrho = float(np.sum(yw2 * dw / w) / sigma ** 2 - np.sum(dw / w))
        du = dlogrho * float(self.transform.dlogrho_du(u)) + du_h
        return float(lp), self.pack(dy, dmu, dls, du)

    def random_init(self, rng: np.random.Generator) -> np.ndarray:
        mu, ls, u = self._random_theta(rng)
        w, _ = self.field.weights(float(self.transform.rho(u)))
        return self.pack(np.exp(ls) * w * rng.normal(size=self.K), mu, ls, u)

    def from_noncentred(self, x_nc: np.ndarray) -> np.ndarray:
        """Map a non-centred state (z, theta) to the centred state (y = sigma w z, theta)."""
        z, mu, ls, u = self._split(x_nc)
        w, _ = self.field.weights(float(self.transform.rho(u)))
        return self.pack(np.exp(ls) * w * z, mu, ls, u)


# ------------------------------------------------------------------------------------ self-tests
def _fd4_grad(f, x, h):
    """Fourth-order central-difference gradient of a scalar function of a vector."""
    g = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h
        g[i] = (8.0 * (f(x + e) - f(x - e)) - (f(x + 2 * e) - f(x - 2 * e))) / (12.0 * h)
    return g


def _selftest():
    from .field import GraphLaplacianBasis, HarmonicBasis
    rng = np.random.default_rng(3)
    ok = True

    def report(name, value, thresh, passed):
        nonlocal ok
        ok &= passed
        print(f"  {'PASS' if passed else 'FAIL'}  {name}: {value:.3e} (threshold {thresh:g})")

    # Gaussian stand-in for the collapsed likelihood (NOT part of the package): strong per-pixel
    # constraints of the size a real nside 16 / nside_beta 2 run gives, so log L is O(1e3).
    def make_gauss_loglike(nb):
        prec = 10.0 ** rng.uniform(2.0, 4.0, size=nb)
        b0 = -3.0 + 0.2 * rng.normal(size=nb)

        def loglike(beta):
            r = beta - b0
            return -0.5 * float(np.sum(prec * r * r)), -prec * r
        return loglike

    h_fd = 1e-3
    for basis in (GraphLaplacianBasis(2), HarmonicBasis(2)):
        loglike = make_gauss_loglike(basis.n_beta)
        field = MaternField(basis)
        for rho_prior in ("loguniform", "pc"):
            hyper = HyperPrior(rho_prior=rho_prior)
            posts = dict(noncentred=NonCentredPosterior(field, hyper, loglike=loglike),
                         centred=CentredPosterior(field, hyper, loglike=loglike))
            for name, post in posts.items():
                worst = 0.0
                worst_rel = 0.0
                for trial in range(3):
                    x = post.random_init(rng)
                    lp, grad = post.logp_and_grad(x)
                    g_fd = _fd4_grad(lambda v: post.logp_and_grad(v)[0], x, h_fd)
                    err = np.abs(g_fd - grad)
                    worst = max(worst, float(np.max(err / (1e-6 * np.abs(grad) + 1e-8))))
                    worst_rel = max(worst_rel, float(np.max(err / np.maximum(np.abs(grad),
                                                                             1e-8))))
                tag = f"{basis.kind} {rho_prior} {name} (K={post.K}, logp~{abs(lp):.0e})"
                report(f"gradient FD |g_fd - g| / (1e-6|g| + 1e-8), {tag}", worst, 1.0,
                       worst < 1.0)
                print(f"        max relative error {worst_rel:.2e}")
            # Consistency of the two parameterizations: same (beta, theta), density ratio is
            # exactly the change-of-variables Jacobian -K log sigma - 1/2 log det K(rho).
            x_nc = posts["noncentred"].random_init(rng)
            x_c = posts["centred"].from_noncentred(x_nc)
            lp_nc = posts["noncentred"].logp_and_grad(x_nc)[0]
            lp_c = posts["centred"].logp_and_grad(x_c)[0]
            mu, sigma, rho = posts["noncentred"].theta(x_nc)
            expected = -post.K * np.log(sigma) - 0.5 * field.logdet_K(rho)
            d = abs((lp_c - lp_nc) - expected) / max(1.0, abs(expected))
            report(f"centred - noncentred = -K log sigma - 1/2 logdet K, {basis.kind} {rho_prior}",
                   d, 1e-10, d < 1e-10)
            db = np.max(np.abs(posts["noncentred"].beta(x_nc) - posts["centred"].beta(x_c)))
            report(f"same beta in both parameterizations, {basis.kind} {rho_prior}", db, 1e-12,
                   db < 1e-12)
            # Centred quadratic form against field.Kinv_quadform on beta - mu 1.
            y = x_c[:post.K]
            q_direct = float(np.sum((y / field.weights(rho)[0]) ** 2))
            q_field = field.Kinv_quadform(rho, posts["centred"].beta(x_c) - mu)
            report(f"Q(y) = Kinv_quadform(beta - mu 1), {basis.kind} {rho_prior}",
                   abs(q_direct - q_field) / q_field, 1e-10,
                   abs(q_direct - q_field) / q_field < 1e-10)

    # Hyperprior terms alone against finite differences (no likelihood), both rho priors.
    basis = GraphLaplacianBasis(1)
    field = MaternField(basis)
    for rho_prior in ("loguniform", "pc"):
        post = NonCentredPosterior(field, HyperPrior(rho_prior=rho_prior),
                                   loglike=lambda beta: (0.0, np.zeros_like(beta)))
        for _ in range(3):
            mu, ls, u = post._random_theta(rng)
            t = np.array([mu, ls, u])
            lp, *grad = post._hyper_terms(*t)
            g_fd = _fd4_grad(lambda v: post._hyper_terms(*v)[0], t, 1e-4)
            err = np.max(np.abs(g_fd - np.array(grad)) / (1e-6 * np.abs(grad) + 1e-8))
            report(f"hyperprior gradient FD, {rho_prior}", err, 1.0, err < 1.0)
    # PC prior tail probabilities by quadrature on rho (untruncated): P(rho < rho0) = 0.05.
    post = NonCentredPosterior(field, HyperPrior(rho_prior="pc"),
                               loglike=lambda beta: (0.0, np.zeros_like(beta)))
    rho0 = post.hyper.pc_rho0
    from scipy.integrate import quad
    pdf = lambda r: post.pc_lam1 / r ** 2 * np.exp(-post.pc_lam1 / r)      # noqa: E731
    p_tail = quad(pdf, 1e-12, rho0)[0]
    report("PC prior P(rho < rho0)", abs(p_tail - PC_ALPHA), 1e-8, abs(p_tail - PC_ALPHA) < 1e-8)
    p_sig = np.exp(-post.pc_lam2 * post.hyper.pc_sigma0)
    report("PC prior P(sigma > sigma0)", abs(p_sig - PC_ALPHA), 1e-12,
           abs(p_sig - PC_ALPHA) < 1e-12)
    # Zero-density path: a likelihood that fails as model._collapsed_core does on a singular A
    # (ValueError) or as np.linalg.cholesky does (LinAlgError) gives log pi = -inf, a zero
    # gradient and one counted error per call, in both parameterizations.
    for exc in (ValueError, np.linalg.LinAlgError):
        def failing(beta, _exc=exc):
            raise _exc("A is not positive definite")
        for Post in (NonCentredPosterior, CentredPosterior):
            post = Post(field, HyperPrior(), loglike=failing)
            lp, grad = post.logp_and_grad(post.random_init(rng))
            good = lp == -np.inf and np.all(grad == 0.0) and post.n_loglike_errors == 1
            report(f"{exc.__name__} in loglike -> -inf, counted, {Post.__name__}",
                   float(post.n_loglike_errors), 1, good)

    print("posterior.py self-test:", "ALL PASS" if ok else "FAILURES")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest() else 1)
