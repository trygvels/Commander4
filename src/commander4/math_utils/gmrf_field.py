"""Hierarchical Gaussian prior for a spectral-index field on HEALPix pixels (port of the E012 v2
prototype nextgen/prototypes/indexfield/field.py and the hyperprior terms of posterior.py).

The index field beta in R^Nb lives on the HEALPix index pixels (nside_beta, Nb = 12 nside_beta^2)
and is a Gaussian field in a fixed orthonormal basis V (Nb x K) with the constant mode removed:

    beta = mu 1 + sigma S(rho) z,     S(rho) = V diag(s_i(rho)) / c(rho),     z ~ N(0, I_K),

    s_i(rho) = (kappa^2 + lam_i)^(-alpha/2),   kappa = sqrt(8 (alpha - 1)) / rho,

    c(rho)^2 = (1/Nb) sum_i s_i(rho)^2,

so that trace(S S^T) = Nb and the pixel-averaged marginal variance of beta - mu 1 is sigma^2
(D009). V is the eigenbasis of the HEALPix 8-neighbour graph Laplacian with its eigenvalues
rescaled to the continuum Laplace-Beltrami spectrum, lam = eigenvalue / (3 h^2),
h = sqrt(4 pi / Nb). rho, kappa and h are in radians. K = Nb - 1.

Hyperparameters theta = (mu, log sigma, u) with rho = rho(u) through RhoTransform,

    log rho = log rho_min + L s(u),   s(u) = 1 / (1 + e^-u),   L = log(rho_max / rho_min),

and the proper hyperpriors of E012 v2 with every Jacobian included (HyperPrior, IndexFieldPrior):

    log N(mu; mu0, s_mu^2),
    log halfN(sigma; sigma_scale) + log sigma            (sampling log sigma),
    log-uniform rho on [rho_min, rho_max]: log |dlog rho / du| = log[L s(u) (1 - s(u))].

Two parameterizations of the field block share one target. Non-centred: x = (z[K], mu, log sigma,
u), prior log density -1/2 z^T z + hyperprior terms. Centred: x = (y[K], mu, log sigma, u),
beta = mu 1 + V y, prior log density

    -1/2 y^T diag((c/s)^2) y / sigma^2 - K log sigma - 1/2 log det C(rho) + hyperprior terms,

C(rho) = S S^T on range(V), log det C = 2 sum_i log(s_i / c). For y = sigma (s/c) z the two
densities differ by exactly the change-of-variables Jacobian -K log sigma - 1/2 log det C(rho)
(IndexFieldPrior.centred_log_jacobian). The data enter only through a collapsed log likelihood
log L(beta) and its gradient dlog L/dbeta, which IndexFieldPrior.logp_and_grad maps to the sampled
coordinates by the chain rule; the likelihood itself lives in compsep.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import healpy as hp
import numpy as np
from numpy.typing import NDArray

NSIDE_BETA_MAX = 16   # Dense eigh of Nb x Nb = 3072 x 3072; above this a sparse solver is needed.


class Basis:
    """Orthonormal basis V (Nb, K) with non-negative spectrum lam (K,) and pixel spacing h.

    Attributes:
        nside_beta: HEALPix resolution of the index pixels.
        n_beta: Nb = 12 nside_beta^2.
        V: (Nb, K) matrix with orthonormal columns, V^T V = I_K; the constant mode is excluded.
        lam: (K,) eigenvalue of the (continuum-scaled) Laplacian associated with each column.
        h: mean pixel spacing sqrt(4 pi / Nb) in radians.
        K: number of modes.
    """

    kind = "abstract"

    def __init__(self, nside_beta: int, V: NDArray, lam: NDArray):
        self.nside_beta = int(nside_beta)
        self.n_beta = hp.nside2npix(self.nside_beta)
        self.V = np.ascontiguousarray(V, dtype=float)
        self.lam = np.ascontiguousarray(lam, dtype=float)
        self.h = float(np.sqrt(4.0 * np.pi / self.n_beta))
        if self.V.shape[0] != self.n_beta or self.V.shape[1] != self.lam.shape[0]:
            raise ValueError("V must be (Nb, K) and lam (K,)")
        if np.any(self.lam < 0):
            raise ValueError("spectrum lam must be non-negative")

    @property
    def K(self) -> int:
        return self.V.shape[1]

    def orthonormality_error(self) -> float:
        """max |V^T V - I|, zero for an exactly orthonormal basis."""
        return float(np.max(np.abs(self.V.T @ self.V - np.eye(self.K))))


class GraphLaplacianBasis(Basis):
    """Eigenbasis of the HEALPix 8-neighbour graph Laplacian G = D - W at nside_beta.

    W is the symmetrized adjacency of hp.get_all_neighbours (an edge exists if either pixel lists
    the other; -1 entries dropped), D = diag(W 1). eigh(G) gives orthonormal eigenvectors; the
    eigenvector with eigenvalue 0 (the constant mode, unique because the graph is connected) is
    dropped. The eigenvalues are rescaled to the continuum Laplace-Beltrami spectrum,
    lam = eigenvalue / (3 h^2), h = sqrt(4 pi / Nb) (D009; the unweighted 8-neighbour Laplacian
    of a plane grid of spacing h approximates -3 h^2 times the continuum Laplacian).

    The neighbour graph is built in the pixel ordering of hp.get_all_neighbours with nest=False,
    so V indexes RING-ordered index pixels, as in the prototype.

    nside_beta is capped at NSIDE_BETA_MAX (dense eigh); larger values raise ValueError.

    Attributes in addition to Basis: G (Nb, Nb) the graph Laplacian and graph_eigenvalues (K,)
    the unscaled eigenvalues of the kept modes.
    """

    kind = "graph"

    def __init__(self, nside_beta: int):
        nside_beta = int(nside_beta)
        if nside_beta < 1 or not hp.isnsideok(nside_beta, nest=True):
            raise ValueError(f"nside_beta must be a power of two >= 1, got {nside_beta}")
        if nside_beta > NSIDE_BETA_MAX:
            raise ValueError(f"nside_beta = {nside_beta} exceeds NSIDE_BETA_MAX = {NSIDE_BETA_MAX}"
                             f" (dense eigh of {hp.nside2npix(nside_beta)}^2 is not supported)")
        nb = hp.nside2npix(nside_beta)
        nbrs = hp.get_all_neighbours(nside_beta, np.arange(nb))      # (8, Nb), -1 if missing
        src = np.broadcast_to(np.arange(nb)[None, :], nbrs.shape)
        ok = nbrs >= 0
        W = np.zeros((nb, nb))
        W[src[ok], nbrs[ok]] = 1.0
        W = np.maximum(W, W.T)
        np.fill_diagonal(W, 0.0)
        G = np.diag(W.sum(axis=1)) - W
        evals, evecs = np.linalg.eigh(G)
        if not (abs(evals[0]) < 1e-9 * max(1.0, evals[-1]) and evals[1] > 1e-9):
            raise RuntimeError("graph Laplacian: expected exactly one zero eigenvalue "
                               f"(got {evals[:3]})")
        h = np.sqrt(4.0 * np.pi / nb)
        super().__init__(nside_beta, evecs[:, 1:], evals[1:] / (3.0 * h * h))
        self.G = G
        self.graph_eigenvalues = evals[1:].copy()


class MaternField:
    """Matern spectrum in a Basis with range parameter rho (radians).

    s_i(rho) = (kappa^2 + lam_i)^(-alpha/2), kappa = sqrt(8 (alpha - 1)) / rho (the Lindgren et
    al. 2011 practical range with smoothness nu = alpha - 1 in two dimensions; alpha = 2 is the
    E012 default), c(rho)^2 = (1/Nb) sum_i s_i^2, S(rho) = V diag(s/c). Derivatives with respect
    to log rho, using d kappa^2 / d log rho = -2 kappa^2:

        ds_i/dlog rho = alpha kappa^2 s_i / (kappa^2 + lam_i)   (= 2 kappa^2 s_i^2 for alpha = 2),
        dc/dlog rho   = (1 / (c Nb)) sum_i s_i ds_i/dlog rho,
        d(s_i/c)/dlog rho = (ds_i - s_i dc/c) / c.

    In the centred parameterization, beta - mu 1 = V y has prior N(0, sigma^2 C(rho)) restricted
    to range(V) with C = S S^T, log det C = sum_i 2 log(s_i / c), and
    v^T C^-1 v = sum_i (c / s_i)^2 (V^T v)_i^2 for v in range(V).
    """

    def __init__(self, basis: Basis, alpha: float = 2.0):
        if alpha <= 1.0:
            raise ValueError("alpha must exceed 1 (kappa = sqrt(8 (alpha - 1)) / rho)")
        self.basis = basis
        self.alpha = float(alpha)
        self.n_beta = basis.n_beta
        self.K = basis.K

    def kappa(self, rho: float) -> float:
        return float(np.sqrt(8.0 * (self.alpha - 1.0)) / rho)

    def spectrum(self, rho: float):
        """Returns (s (K,), c, ds_dlogrho (K,), dc_dlogrho)."""
        if not rho > 0:
            raise ValueError("rho must be positive")
        k2 = self.kappa(rho) ** 2
        q = k2 + self.basis.lam
        s = q ** (-0.5 * self.alpha)
        ds = self.alpha * k2 * s / q
        c = float(np.sqrt(np.sum(s * s) / self.n_beta))
        dc = float(np.sum(s * ds) / (c * self.n_beta))
        return s, c, ds, dc

    def weights(self, rho: float):
        """Returns (w, dw) with w = s / c the normalized spectrum and dw = dw / dlog rho."""
        s, c, ds, dc = self.spectrum(rho)
        w = s / c
        dw = (ds - s * dc / c) / c
        return w, dw

    def S_times(self, rho: float, z: NDArray) -> NDArray:
        """V diag(s/c) z, shape (Nb,)."""
        w, _ = self.weights(rho)
        return self.basis.V @ (w * z)

    def ST_times(self, rho: float, g: NDArray) -> NDArray:
        """diag(s/c) V^T g, shape (K,)."""
        w, _ = self.weights(rho)
        return w * (self.basis.V.T @ g)

    def dS_dlogrho_times(self, rho: float, z: NDArray) -> NDArray:
        """V diag(d(s/c)/dlog rho) z, shape (Nb,)."""
        _, dw = self.weights(rho)
        return self.basis.V @ (dw * z)

    def logdet_K(self, rho: float) -> float:
        """log det of C(rho) = S S^T on range(V): sum_i 2 log(s_i / c). The prototype's name is
        kept; the plan (M2.4 section 1) calls this matrix C(rho)."""
        w, _ = self.weights(rho)
        return float(2.0 * np.sum(np.log(w)))

    def Kinv_quadform(self, rho: float, v: NDArray) -> float:
        """v^T V diag((c/s)^2) V^T v for v in R^Nb (the component outside range(V) is ignored)."""
        w, _ = self.weights(rho)
        coeff = self.basis.V.T @ v
        return float(np.sum((coeff / w) ** 2))

    def pixel_variances(self, rho: float) -> NDArray:
        """diag(S S^T), the per-pixel marginal variance of (beta - mu 1) / sigma; its mean is 1."""
        w, _ = self.weights(rho)
        return (self.basis.V ** 2) @ (w * w)


class RhoTransform:
    """Logit map between the bounded correlation length rho in (rho_min, rho_max) and u in R:

        log rho = log rho_min + L s(u),   s(u) = 1 / (1 + e^-u),   L = log(rho_max / rho_min),
        u = logit((log rho - log rho_min) / L),
        dlog rho / du = L s(u) (1 - s(u)),   log |dlog rho / du| = log L + log s + log(1 - s).
    """

    def __init__(self, rho_min: float, rho_max: float):
        if not (0 < rho_min < rho_max):
            raise ValueError("need 0 < rho_min < rho_max")
        self.rho_min = float(rho_min)
        self.rho_max = float(rho_max)
        self.log_rho_min = float(np.log(rho_min))
        self.L = float(np.log(rho_max / rho_min))

    @staticmethod
    def _sigmoid(u):
        return 1.0 / (1.0 + np.exp(-np.asarray(u, dtype=float)))

    def rho(self, u):
        return np.exp(self.log_rho_min + self.L * self._sigmoid(u))

    def u(self, rho):
        f = (np.log(np.asarray(rho, dtype=float)) - self.log_rho_min) / self.L
        return np.log(f) - np.log1p(-f)

    def dlogrho_du(self, u):
        s = self._sigmoid(u)
        return self.L * s * (1.0 - s)

    def log_jacobian(self, u):
        """log |dlog rho / du| = log L - softplus(-u) - softplus(u)."""
        u = np.asarray(u, dtype=float)
        return np.log(self.L) - np.logaddexp(0.0, -u) - np.logaddexp(0.0, u)

    def dlog_jacobian_du(self, u):
        """d/du log |dlog rho / du| = 1 - 2 s(u)."""
        return 1.0 - 2.0 * self._sigmoid(u)


@dataclass
class HyperPrior:
    """Hyperprior settings of E012 v2 (M2.4 section 1): mu ~ N(mu0, s_mu^2) with (mu0, s_mu)
    from the component's spectral_index_prior, sigma ~ halfN(0, sigma_scale), rho log-uniform on
    [rho_min, rho_max]. rho_min defaults to the pixel spacing h once it is known (resolved(h)).
    Only rho_prior = "loguniform" is ported; the prototype's PC prior is not part of the target.
    """

    mu0: float
    s_mu: float
    sigma_scale: float = 0.3
    rho_prior: str = "loguniform"
    rho_min: float | None = None
    rho_max: float = float(np.pi)

    def resolved(self, h: float) -> "HyperPrior":
        if self.rho_prior != "loguniform":
            raise ValueError(f"rho_prior {self.rho_prior!r} is not supported (only 'loguniform')")
        if not self.s_mu > 0 or not self.sigma_scale > 0:
            raise ValueError("s_mu and sigma_scale must be positive")
        return replace(self, rho_min=h if self.rho_min is None else self.rho_min)


class IndexFieldPrior:
    """Prior of the index field in the sampled coordinates, both parameterizations.

    x = (field block of length K, mu, log sigma, u); dim = K + 3. The field block is z (non-centred,
    beta = mu 1 + sigma S(rho) z) or y (centred, beta = mu 1 + V y), see the module docstring.

    Non-centred prior and chain rule, with g = dlog L/dbeta and w = s/c (posterior.py
    NonCentredPosterior):
        log pi = log L(beta) - 1/2 z^T z + hyperprior terms
        d/dz = sigma w * (V^T g) - z
        d/dmu = 1^T g + hyperprior
        d/dlog sigma = sigma (S z)^T g + hyperprior and Jacobian
        d/du = sigma z^T (dS/dlog rho)^T g dlog rho/du + hyperprior and Jacobian.

    Centred, with dw = dw/dlog rho and Q = sum_i y_i^2 / w_i^2 (posterior.py CentredPosterior):
        log pi = log L(beta) - 1/2 Q / sigma^2 - K log sigma - 1/2 log det C(rho) + hyperprior
        d/dy = V^T g - y / (w^2 sigma^2)
        d/dmu = 1^T g + hyperprior
        d/dlog sigma = Q / sigma^2 - K + hyperprior and Jacobian
        d/du = [sum_i y_i^2 dw_i / (w_i^3 sigma^2) - sum_i dw_i / w_i] dlog rho/du + hyperprior.
    """

    def __init__(self, field: MaternField, hyper: HyperPrior, centred: bool = False):
        self.field = field
        self.basis = field.basis
        self.n_beta = field.n_beta
        self.K = field.K
        self.hyper = hyper.resolved(self.basis.h)
        self.transform = RhoTransform(self.hyper.rho_min, self.hyper.rho_max)
        self.centred = bool(centred)

    @property
    def dim(self) -> int:
        return self.K + 3

    def split(self, x: NDArray):
        """(block (K,), mu, log sigma, u) of a packed state."""
        x = np.asarray(x, dtype=float)
        if x.shape != (self.dim,):
            raise ValueError(f"x must have shape ({self.dim},)")
        return x[:self.K], float(x[self.K]), float(x[self.K + 1]), float(x[self.K + 2])

    def pack(self, block: NDArray, mu: float, log_sigma: float, u: float) -> NDArray:
        return np.concatenate([np.asarray(block, dtype=float), [mu, log_sigma, u]])

    def theta(self, x: NDArray):
        """(mu, sigma, rho) of a state."""
        _, mu, ls, u = self.split(x)
        return mu, float(np.exp(ls)), float(self.transform.rho(u))

    def hyper_terms(self, mu: float, log_sigma: float, u: float):
        """Hyperprior log density in the sampled coordinates and its gradient
        (d/dmu, d/dlog sigma, d/du), Jacobians included (posterior.py _hyper_terms):

            -1/2 ((mu - mu0) / s_mu)^2 - 1/2 (sigma / sigma_scale)^2 + log sigma
            + log L + log s(u) + log(1 - s(u)),

        the log-uniform density of rho being constant in log rho (its normalization dropped).
        """
        hp_ = self.hyper
        sigma = np.exp(log_sigma)
        tr = self.transform
        # mu ~ N(mu0, s_mu^2)
        lp = -0.5 * ((mu - hp_.mu0) / hp_.s_mu) ** 2
        dmu = -(mu - hp_.mu0) / hp_.s_mu ** 2
        # u: Jacobian of the logit map
        lp += tr.log_jacobian(u)
        du = tr.dlog_jacobian_du(u)
        # sigma ~ halfN(0, sigma_scale), + log sigma Jacobian
        lp += -0.5 * (sigma / hp_.sigma_scale) ** 2 + log_sigma
        dls = 1.0 - (sigma / hp_.sigma_scale) ** 2
        return float(lp), float(dmu), float(dls), float(du)

    def beta(self, x: NDArray) -> NDArray:
        """The index field beta (Nb,) of a packed state in the active parameterization."""
        block, mu, ls, u = self.split(x)
        if self.centred:
            return mu + self.basis.V @ block
        return mu + np.exp(ls) * self.field.S_times(float(self.transform.rho(u)), block)

    def centred_log_jacobian(self, x: NDArray) -> float:
        """-K log sigma - 1/2 log det C(rho): the log density ratio centred minus non-centred at
        the same (beta, theta), i.e. the Jacobian of y = sigma (s/c) z."""
        _, sigma, rho = self.theta(x)
        return float(-self.K * np.log(sigma) - 0.5 * self.field.logdet_K(rho))

    def from_noncentred(self, x_nc: NDArray) -> NDArray:
        """Map a non-centred state (z, theta) to the centred state (y = sigma w z, theta)."""
        z, mu, ls, u = self.split(x_nc)
        w, _ = self.field.weights(float(self.transform.rho(u)))
        return self.pack(np.exp(ls) * w * z, mu, ls, u)

    def to_noncentred(self, x_c: NDArray) -> NDArray:
        """Map a centred state (y, theta) to the non-centred state (z = y / (sigma w), theta)."""
        y, mu, ls, u = self.split(x_c)
        w, _ = self.field.weights(float(self.transform.rho(u)))
        return self.pack(y / (np.exp(ls) * w), mu, ls, u)

    def prior_logp_and_grad(self, x: NDArray):
        """Prior log density (field block and hyperpriors, Jacobians included) and its gradient
        in the sampled coordinates; the likelihood is not included (logp_and_grad adds it)."""
        return self.logp_and_grad(x, loglike=None)

    def logp_and_grad(self, x: NDArray, loglike=None):
        """Log posterior density and gradient in the sampled coordinates.

        loglike(beta) -> (log L, dlog L/dbeta) is the collapsed likelihood; None means prior only.
        A non-finite density gives (-inf, zeros) (the nuts.py contract); a loglike raising
        ValueError or LinAlgError (a non-positive-definite amplitude normal matrix) is treated the
        same way.
        """
        block, mu, ls, u = self.split(x)
        sigma = np.exp(ls)
        rho = float(self.transform.rho(u))
        w, dw = self.field.weights(rho)
        V = self.basis.V
        if self.centred:
            beta = mu + V @ block
        else:
            Sz = V @ (w * block)
            beta = mu + sigma * Sz
        if loglike is None:
            logL, g = 0.0, np.zeros(self.n_beta)
        else:
            try:
                logL, g = loglike(beta)
            except (ValueError, np.linalg.LinAlgError):
                return -np.inf, np.zeros(self.dim)
            if not np.isfinite(logL):
                return -np.inf, np.zeros(self.dim)
            logL, g = float(logL), np.asarray(g, dtype=float)
        lp_h, dmu_h, dls_h, du_h = self.hyper_terms(mu, ls, u)
        Vtg = V.T @ g
        if self.centred:
            y = block
            yw2 = (y / w) ** 2
            Q = float(yw2.sum())
            logdet_K = float(2.0 * np.sum(np.log(w)))
            lp = logL - 0.5 * Q / sigma ** 2 - self.K * ls - 0.5 * logdet_K + lp_h
            if not np.isfinite(lp):
                return -np.inf, np.zeros(self.dim)
            dblock = Vtg - y / (w * w * sigma ** 2)
            dmu = float(g.sum()) + dmu_h
            dls = Q / sigma ** 2 - self.K + dls_h
            dlogrho = float(np.sum(yw2 * dw / w) / sigma ** 2 - np.sum(dw / w))
            du = dlogrho * float(self.transform.dlogrho_du(u)) + du_h
        else:
            z = block
            lp = logL - 0.5 * float(z @ z) + lp_h
            if not np.isfinite(lp):
                return -np.inf, np.zeros(self.dim)
            dblock = sigma * w * Vtg - z
            dmu = float(g.sum()) + dmu_h
            dls = sigma * float(Sz @ g) + dls_h
            du = sigma * float((dw * z) @ Vtg) * float(self.transform.dlogrho_du(u)) + du_h
        return float(lp), self.pack(dblock, dmu, dls, du)


class _RegionGrid:
    """The HEALPix pixelization of an index field without a basis: `nside_beta`, `n_beta`, `h`."""

    def __init__(self, nside_beta: int):
        self.nside_beta = int(nside_beta)
        self.n_beta = 12 * self.nside_beta ** 2
        self.h = float(np.sqrt(4.0 * np.pi / self.n_beta))


class _IndependentHyper:
    """(mu0, s_mu) holder so that `IndexSpec` can check its prior against the field's."""

    def __init__(self, mu0: float, s_mu: float):
        self.mu0 = float(mu0)
        self.s_mu = float(s_mu)


class IndependentIndexPrior:
    """Independent region indices at nside_beta: the fixed-nside ladder of parametric fits.

    Target in the sampled coordinates x = beta (n_beta,):
        log pi = log L(beta) + sum_r [-1/2 ((beta_r - mu0) / sd0)^2]
        d/dbeta_r = dlog L/dbeta_r - (beta_r - mu0) / sd0^2,
    no hyperparameters, no spatial coupling. Same interface as `IndexFieldPrior` where the
    sampler needs it (`dim`, `n_beta`, `K`, `basis.nside_beta`, `hyper.mu0/s_mu`, `beta`,
    `logp_and_grad`, `prior_logp_and_grad`); `centred` is False (beta itself is sampled).
    """

    centred = False

    def __init__(self, nside_beta: int, mu0: float, sd0: float):
        if nside_beta < 1:
            raise ValueError("IndependentIndexPrior needs nside_beta >= 1")
        if not sd0 > 0.0:
            raise ValueError("IndependentIndexPrior needs sd0 > 0")
        self.basis = _RegionGrid(nside_beta)
        self.n_beta = self.basis.n_beta
        self.K = self.n_beta
        self.hyper = _IndependentHyper(mu0, sd0)
        self.mu0 = float(mu0)
        self.sd0 = float(sd0)

    @property
    def dim(self) -> int:
        return self.n_beta

    def beta(self, x: NDArray) -> NDArray:
        x = np.asarray(x, dtype=float)
        if x.shape != (self.dim,):
            raise ValueError(f"x must have shape ({self.dim},)")
        return x.copy()

    def prior_logp_and_grad(self, x: NDArray):
        return self.logp_and_grad(x, loglike=None)

    def logp_and_grad(self, x: NDArray, loglike=None):
        """Log posterior density and gradient; loglike(beta) -> (log L, dlog L/dbeta) or None for
        the prior alone. Non-finite values and raising likelihoods give (-inf, zeros)."""
        beta = self.beta(x)
        if loglike is None:
            logL, g = 0.0, np.zeros(self.n_beta)
        else:
            try:
                logL, g = loglike(beta)
            except (ValueError, np.linalg.LinAlgError):
                return -np.inf, np.zeros(self.dim)
            if not np.isfinite(logL):
                return -np.inf, np.zeros(self.dim)
            logL, g = float(logL), np.asarray(g, dtype=float)
        r = (beta - self.mu0) / self.sd0
        lp = logL - 0.5 * float(r @ r)
        if not np.isfinite(lp):
            return -np.inf, np.zeros(self.dim)
        return lp, g - r / self.sd0
