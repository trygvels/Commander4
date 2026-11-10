"""A scale-dependent noise model for the Planck PR4 LFI Q/U maps at low resolution.

The map files carry a white-noise covariance W (2 by 2 per pixel). The ring-half difference shows
that an nside-16 pixel average is noisier than W says at every scale (E038 b). Model:

    N = W^(1/2) R W^(1/2),      R = r I + Y diag(A l^-gamma) Y^T,

with W^(1/2) the symmetric square root per pixel, Y the real spin-2 harmonics of E and B for
2 <= l <= lmax on the grid, scaled so that Y^T Y is close to the identity, and (r, A, gamma) fitted by
maximum likelihood to the whitened half difference z = [(W1 + W2) / 4]^(-1/2) (d1 - d2) / 2:

    -2 ln L = z^T R^-1 z + ln det R.

Everything is NESTED, (Q elements, then U elements).
"""
import healpy as hp
import numpy as np
from scipy.optimize import minimize

from . import data as D


def spin2_basis(nside, lmax):
    """Y (2 npix, nmodes) and the multipole of each column. Columns are E then B; within each, (l, m)
    with the real part for m >= 0 and the imaginary part for m > 0."""
    npix = hp.nside2npix(nside)
    cols, ells = [], []
    size = hp.Alm.getsize(lmax)
    zero = np.zeros(size, dtype=complex)
    for which in (1, 2):
        for l in range(2, lmax + 1):
            for m in range(0, l + 1):
                for part in ((1.0,) if m == 0 else (1.0, 1.0j)):
                    alm = np.zeros(size, dtype=complex)
                    alm[hp.Alm.getidx(lmax, l, m)] = part
                    alms = [zero, alm, zero] if which == 1 else [zero, zero, alm]
                    _, q, u = hp.alm2map(alms, nside, lmax=lmax, pol=True)
                    v = np.concatenate([hp.reorder(q, r2n=True), hp.reorder(u, r2n=True)])
                    cols.append(v * np.sqrt(4 * np.pi / npix) / (1.0 if m == 0 else np.sqrt(2.0)))
                    ells.append(l)
    return np.array(cols).T, np.array(ells)


def sqrt_blocks(a, b, c):
    """Symmetric square root and inverse square root of the 2 by 2 matrices [[a, c], [c, b]], elementwise.
    Returns (s_a, s_b, s_c), (i_a, i_b, i_c)."""
    tr, det = a + b, a * b - c * c
    s = np.sqrt(det)
    t = np.sqrt(tr + 2 * s)
    sa, sb, sc = (a + s) / t, (b + s) / t, c / t
    return (sa, sb, sc), (sb / s, sa / s, -sc / s)


def cov_blocks(ninv):
    n = ninv.shape[0] // 2
    i = np.arange(n)
    a, b, c = ninv[i, i], ninv[n + i, n + i], ninv[i, n + i]
    det = a * b - c * c
    return b / det, a / det, -c / det


def apply_blocks(blk, v):
    a, b, c = blk
    n = a.size
    return np.concatenate([a * v[:n] + c * v[n:], c * v[:n] + b * v[n:]])


def whitened_half_difference(freq, nside=16):
    """z, and the white covariance blocks of the full-mission map."""
    d1, n1 = D.load_planck_qu(D.PR4_LFI.format(freq=freq, split="full-ringhalf-1"), nside)
    d2, n2 = D.load_planck_qu(D.PR4_LFI.format(freq=freq, split="full-ringhalf-2"), nside)
    q1, u1, x1 = cov_blocks(n1)
    q2, u2, x2 = cov_blocks(n2)
    _, isq = sqrt_blocks((q1 + q2) / 4, (u1 + u2) / 4, (x1 + x2) / 4)
    return apply_blocks(isq, (d1 - d2) / 2)


class NoiseModel:
    def __init__(self, nside=16, lmax=None):
        self.nside = nside
        self.lmax = 3 * nside - 1 if lmax is None else lmax
        self.Y, self.ell = spin2_basis(nside, self.lmax)
        self.G = self.Y.T @ self.Y

    def dl(self, A, gamma):
        return A * (self.ell / 10.0) ** (-gamma)

    def neg2loglike(self, theta, z):
        """-2 ln L for theta = (ln r, ln A, gamma), by the Woodbury identity."""
        r, A, gamma = np.exp(theta[0]), np.exp(theta[1]), theta[2]
        dl = self.dl(A, gamma)
        # R = r I + Y D Y^T; R^-1 = (I - Y (r D^-1 + G)^-1 Y^T) / r; det R = r^n det(I + D G / r)
        M = np.diag(r / dl) + self.G
        c = np.linalg.cholesky(M)
        yz = self.Y.T @ z
        w = np.linalg.solve(c, yz)
        quad = (z @ z - w @ w) / r
        logdet = z.size * np.log(r) + 2 * np.sum(np.log(np.diag(c))) + np.sum(np.log(dl / r))
        return quad + logdet

    def fit(self, z, theta0=(0.2, 0.0, 1.0)):
        res = minimize(self.neg2loglike, np.array(theta0), args=(z,), method="Nelder-Mead",
                       options=dict(xatol=1e-3, fatol=1e-3, maxiter=400))
        return res

    def R(self, theta):
        r, A, gamma = np.exp(theta[0]), np.exp(theta[1]), theta[2]
        return r * np.eye(self.Y.shape[0]) + (self.Y * self.dl(A, gamma)[None, :]) @ self.Y.T


def lfi_ninv(freq, theta, model, nside=16):
    """Dense N^-1 of the full-mission PR4 LFI Q/U map at `nside` under the fitted model, and the map."""
    d, nw = D.load_planck_qu(D.PR4_LFI.format(freq=freq, split="full"), nside)
    _, isq = sqrt_blocks(*cov_blocks(nw))
    Rinv = np.linalg.inv(model.R(theta))
    n = d.size // 2
    a, b, c = isq
    # W^(-1/2) as a sparse 2 by 2 block operator applied on both sides
    left = np.concatenate([a[:, None] * Rinv[:n] + c[:, None] * Rinv[n:], c[:, None] * Rinv[:n] + b[:, None] * Rinv[n:]])
    out = np.concatenate([left[:, :n] * a[None, :] + left[:, n:] * c[None, :],
                          left[:, :n] * c[None, :] + left[:, n:] * b[None, :]], axis=1)
    return d, 0.5 * (out + out.T)
