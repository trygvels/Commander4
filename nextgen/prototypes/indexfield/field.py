"""Basis expansions, Matern spectra and the correlation-length transform for the index field.

The spectral-index field beta in R^Nb lives on the HEALPix index pixels (nside_beta, RING order).
E012 v2 models it as a hierarchical Gaussian field in a fixed orthonormal basis V (Nb x K) with
the constant mode removed:

    beta = mu 1 + sigma S(rho) z,     S(rho) = V diag(s_i(rho)) / c(rho),     z ~ N(0, I_K),

    s_i(rho) = (kappa^2 + lam_i)^(-alpha/2),   kappa = sqrt(8 (alpha - 1)) / rho,

    c(rho)^2 = (1/Nb) sum_i s_i(rho)^2,

so that trace(S S^T) = Nb and the pixel-averaged marginal variance of beta - mu 1 is sigma^2.
Two bases are provided: the eigenvectors of the HEALPix 8-neighbour graph Laplacian with the
eigenvalues rescaled to the continuum Laplace-Beltrami operator (D009), and real spherical
harmonics 1 <= l <= lmax orthonormalized by QR with lam = l (l + 1). rho and kappa are in radians.

The correlation length is sampled through u = logit((log rho - log rho_min) / L),
L = log(rho_max / rho_min), with log Jacobian log |d log rho / d u| = log[L s(u) (1 - s(u))].
"""
from __future__ import annotations

import numpy as np
import healpy as hp
from scipy.special import sph_harm_y


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

    def __init__(self, nside_beta: int, V: np.ndarray, lam: np.ndarray):
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

    Attributes in addition to Basis: G (Nb, Nb) the graph Laplacian and graph_eigenvalues (K,)
    the unscaled eigenvalues of the kept modes.
    """

    kind = "graph"

    def __init__(self, nside_beta: int):
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


def real_spherical_harmonics(nside: int, lmax: int):
    """Real spherical harmonics 1 <= l <= lmax at the pixel centres of a HEALPix map.

    Columns ordered by l, then m = -l..l. With Y_l^m the complex harmonics (Condon-Shortley
    phase) the real harmonics are Y_l0 for m = 0, sqrt(2) (-1)^m Re Y_l^m for m > 0 and
    sqrt(2) (-1)^m Im Y_l^|m| for m < 0. The columns are scaled by h = sqrt(4 pi / Npix) so that
    they are close to orthonormal under the pixel sum (exactly orthonormal in the continuum).

    Returns (Y (Npix, K), ells (K,), ems (K,)) with K = (lmax + 1)^2 - 1.
    """
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    cols, ells, ems = [], [], []
    for ell in range(1, lmax + 1):
        for m in range(-ell, ell + 1):
            y = sph_harm_y(ell, abs(m), theta, phi)
            if m == 0:
                col = y.real
            elif m > 0:
                col = np.sqrt(2.0) * (-1) ** m * y.real
            else:
                col = np.sqrt(2.0) * (-1) ** m * y.imag
            cols.append(col)
            ells.append(ell)
            ems.append(m)
    Y = np.column_stack(cols) * np.sqrt(4.0 * np.pi / npix)
    return Y, np.array(ells), np.array(ems)


class HarmonicBasis(Basis):
    """Real spherical harmonics 1 <= l <= lmax (default 2 nside_beta) at the index-pixel centres,
    orthonormalized by a thin QR decomposition (columns ordered by l then m, each Q column signed
    so that the corresponding diagonal of R is positive, which keeps Q_k close to the k-th
    harmonic). lam = l (l + 1) of the column's l. K = (lmax + 1)^2 - 1 < Nb.

    Attributes in addition to Basis: ells, ems (K,), and qr_departure = max |R - I| as a measure of
    how far the pixel-evaluated harmonics were from orthonormal before the QR.
    """

    kind = "harmonic"

    def __init__(self, nside_beta: int, lmax: int | None = None):
        if lmax is None:
            lmax = 2 * nside_beta
        nb = hp.nside2npix(nside_beta)
        if (lmax + 1) ** 2 - 1 >= nb:
            raise ValueError(f"lmax = {lmax} gives K >= Nb = {nb}; the QR would be rank deficient")
        Y, ells, ems = real_spherical_harmonics(nside_beta, lmax)
        Q, R = np.linalg.qr(Y, mode="reduced")
        sign = np.sign(np.diag(R))
        sign[sign == 0] = 1.0
        Q = Q * sign[None, :]
        R = R * sign[:, None]
        if np.min(np.diag(R)) < 1e-8:
            raise RuntimeError("pixelized harmonics are rank deficient at this nside/lmax")
        super().__init__(nside_beta, Q, ells * (ells + 1.0))
        self.lmax = int(lmax)
        self.ells = ells
        self.ems = ems
        self.qr_departure = float(np.max(np.abs(R - np.eye(R.shape[0]))))


class MaternField:
    """Matern spectrum in a Basis with range parameter rho (radians).

    s_i(rho) = (kappa^2 + lam_i)^(-alpha/2), kappa = sqrt(8 (alpha - 1)) / rho (the Lindgren et
    al. 2011 practical range with smoothness nu = alpha - 1 in two dimensions; alpha = 2 is the
    E012 default), c(rho)^2 = (1/Nb) sum_i s_i^2, S(rho) = V diag(s/c). Derivatives with respect
    to log rho, using d kappa^2 / d log rho = -2 kappa^2:

        ds_i/dlog rho = alpha kappa^2 s_i / (kappa^2 + lam_i)   (= 2 kappa^2 s_i^2 for alpha = 2),
        dc/dlog rho   = (1 / (c Nb)) sum_i s_i ds_i/dlog rho,
        d(s_i/c)/dlog rho = (ds_i - s_i dc/c) / c.

    In the centred parameterization, beta - mu 1 = V y has prior N(0, sigma^2 K(rho)) restricted
    to range(V) with K = S S^T, log det K = sum_i 2 log(s_i / c), and
    v^T K^-1 v = sum_i (c / s_i)^2 (V^T v)_i^2 for v in range(V).
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

    def S_times(self, rho: float, z: np.ndarray) -> np.ndarray:
        """V diag(s/c) z, shape (Nb,)."""
        w, _ = self.weights(rho)
        return self.basis.V @ (w * z)

    def ST_times(self, rho: float, g: np.ndarray) -> np.ndarray:
        """diag(s/c) V^T g, shape (K,)."""
        w, _ = self.weights(rho)
        return w * (self.basis.V.T @ g)

    def dS_dlogrho_times(self, rho: float, z: np.ndarray) -> np.ndarray:
        """V diag(d(s/c)/dlog rho) z, shape (Nb,)."""
        _, dw = self.weights(rho)
        return self.basis.V @ (dw * z)

    def logdet_K(self, rho: float) -> float:
        """log det of K = S S^T on range(V): sum_i 2 log(s_i / c)."""
        w, _ = self.weights(rho)
        return float(2.0 * np.sum(np.log(w)))

    def Kinv_quadform(self, rho: float, v: np.ndarray) -> float:
        """v^T V diag((c/s)^2) V^T v for v in R^Nb (the component outside range(V) is ignored)."""
        w, _ = self.weights(rho)
        coeff = self.basis.V.T @ v
        return float(np.sum((coeff / w) ** 2))

    def pixel_variances(self, rho: float) -> np.ndarray:
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


# ------------------------------------------------------------------------------------ self-tests
def _fd4(f, x, h):
    """Fourth-order central difference of a scalar or vector function of a scalar."""
    return (8.0 * (f(x + h) - f(x - h)) - (f(x + 2 * h) - f(x - 2 * h))) / (12.0 * h)


def _selftest():
    rng = np.random.default_rng(12)
    ok = True

    def report(name, value, thresh, passed):
        nonlocal ok
        ok &= passed
        print(f"  {'PASS' if passed else 'FAIL'}  {name}: {value:.3e} (threshold {thresh:g})")

    # 1. Orthonormality of both bases at several resolutions.
    print("orthonormality of V")
    for nside in (1, 2, 4, 8):
        for cls in (GraphLaplacianBasis, HarmonicBasis):
            b = cls(nside)
            report(f"{b.kind} nside {nside} K={b.K}", b.orthonormality_error(), 1e-10,
                   b.orthonormality_error() < 1e-10)
            if b.kind == "harmonic":
                print(f"        pre-QR departure from orthonormality max|R - I| = "
                      f"{b.qr_departure:.3e}")
            const = np.abs(b.V.sum(axis=0)).max() / np.sqrt(b.n_beta)
            report(f"{b.kind} nside {nside} constant-mode leakage |1^T V|/sqrt(Nb)", const,
                   1e-10, const < 1e-10)

    # 2. Graph Laplacian eigenvalue scale against l(l+1), nside 8, by projecting onto harmonics.
    print("graph Laplacian spectrum against l(l+1) at nside 8")
    g, hb = GraphLaplacianBasis(8), HarmonicBasis(8)
    proj = hb.V.T @ g.V                        # (K_h, K_g): harmonic content of each graph mode
    power = (proj ** 2).sum(axis=0)
    l_eff = (proj ** 2 * hb.lam[:, None]).sum(axis=0) / power
    ratio = g.lam / l_eff
    n_low = 15                                  # l = 1, 2, 3: 3 + 5 + 7 modes
    print(f"        lowest {n_low} graph modes: lam/(l(l+1)) = "
          f"{np.array2string(ratio[:n_low], precision=3)}")
    print(f"        harmonic power captured (lmax 16): min {power[:n_low].min():.4f}")
    report("min ratio, lowest modes", ratio[:n_low].min(), 0.8, ratio[:n_low].min() >= 0.8)
    report("max ratio, lowest modes", ratio[:n_low].max(), 1.1, ratio[:n_low].max() <= 1.1)
    sel = (hb.ells == 8)
    ray = np.einsum("ik,ij,jk->k", hb.V[:, sel], g.G, hb.V[:, sel]) / (3 * g.h ** 2) / 72.0
    print(f"        Rayleigh quotient of l = 8 harmonics / 72: mean {ray.mean():.3f} "
          f"(E012 review quoted 2.65/3 = 0.88)")

    # 3. Spectrum derivatives and the S / K operators.
    print("Matern spectrum and operators")
    for basis in (GraphLaplacianBasis(2), HarmonicBasis(2)):
        for alpha in (2.0, 3.0):
            fld = MaternField(basis, alpha=alpha)
            for rho in (basis.h, 0.3, np.pi):
                lr = np.log(rho)
                s, c, ds, dc = fld.spectrum(rho)
                ds_fd = _fd4(lambda t: fld.spectrum(np.exp(t))[0], lr, 1e-3)
                dc_fd = _fd4(lambda t: fld.spectrum(np.exp(t))[1], lr, 1e-3)
                e1 = np.max(np.abs(ds_fd - ds) / (1e-6 * np.abs(ds) + 1e-12))
                e2 = abs(dc_fd - dc) / (1e-6 * abs(dc) + 1e-12)
                z = rng.normal(size=basis.K)
                dSz_fd = _fd4(lambda t: fld.S_times(np.exp(t), z), lr, 1e-3)
                e3 = np.max(np.abs(dSz_fd - fld.dS_dlogrho_times(rho, z))
                            / (1e-6 * np.abs(dSz_fd) + 1e-12))
                tag = f"{basis.kind} alpha {alpha} rho {np.degrees(rho):.1f} deg"
                report(f"ds/dlogrho FD, {tag}", e1, 1.0, e1 < 1.0)
                report(f"dc/dlogrho FD, {tag}", e2, 1.0, e2 < 1.0)
                report(f"dS/dlogrho z FD, {tag}", e3, 1.0, e3 < 1.0)
                # trace(S S^T) = Nb, Kinv_quadform(S z) = |z|^2, logdet matches explicit product.
                tr = fld.pixel_variances(rho).sum() / basis.n_beta
                report(f"mean pixel variance = 1, {tag}", abs(tr - 1), 1e-12, abs(tr - 1) < 1e-12)
                q = fld.Kinv_quadform(rho, fld.S_times(rho, z))
                report(f"Kinv_quadform(S z) = |z|^2, {tag}", abs(q - z @ z) / (z @ z), 1e-10,
                       abs(q - z @ z) / (z @ z) < 1e-10)
                Kmat = fld.basis.V @ np.diag((s / c) ** 2) @ fld.basis.V.T
                ld = np.sum(np.log(np.linalg.eigvalsh(Kmat)[-basis.K:]))
                report(f"logdet_K vs eigvalsh, {tag}", abs(ld - fld.logdet_K(rho)), 1e-8,
                       abs(ld - fld.logdet_K(rho)) < 1e-8)
                # ST_times is the transpose of S_times.
                gvec = rng.normal(size=basis.n_beta)
                adj = abs(gvec @ fld.S_times(rho, z) - fld.ST_times(rho, gvec) @ z)
                report(f"<g, S z> = <S^T g, z>, {tag}", adj, 1e-9, adj < 1e-9)

    # 4. RhoTransform round trip, Jacobian and its derivative.
    print("RhoTransform")
    tr = RhoTransform(0.05, np.pi)
    for u in (-6.0, -1.0, 0.0, 2.5, 7.0):
        rt = abs(tr.u(tr.rho(u)) - u)
        report(f"round trip u={u}", rt, 1e-9, rt < 1e-9)
        d_fd = _fd4(lambda t: np.log(tr.rho(t)), u, 1e-3)
        e = abs(d_fd - tr.dlogrho_du(u)) / (1e-6 * abs(d_fd) + 1e-12)
        report(f"dlogrho/du FD u={u}", e, 1.0, e < 1.0)
        lj = abs(np.log(abs(tr.dlogrho_du(u))) - tr.log_jacobian(u))
        report(f"log_jacobian = log|dlogrho/du| u={u}", lj, 1e-10, lj < 1e-10)
        dj_fd = _fd4(tr.log_jacobian, u, 1e-3)
        e = abs(dj_fd - tr.dlog_jacobian_du(u)) / (1e-6 * abs(dj_fd) + 1e-10)
        report(f"dlog_jacobian/du FD u={u}", e, 1.0, e < 1.0)

    print("field.py self-test:", "ALL PASS" if ok else "FAILURES")
    return ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if _selftest() else 1)
