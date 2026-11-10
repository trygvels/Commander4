"""Beam differences between bands at low resolution.

The released low-resolution maps are pixel averages of maps with different beams. To multipole
lmax = 3 nside - 1 the difference between two beams is a multipole-dependent factor, which acts on a
(Q, U) vector at nside 16 as

    S = 1 + Y diag(f_l - 1) Y^+,      f_l = b_l(target) / b_l(band),

with Y the real spin-2 harmonics of E and B on the grid and Y^+ = (Y^T Y)^-1 Y^T. S is invertible
(S^-1 has 1 / f_l), so d' = S d with N'^-1 = S^-T N^-1 S^-1 loses nothing: it is the same likelihood for
a model in which the band sees the target-resolution sky through S^-1. Structure above lmax is not
touched, which leaves a residual that the E038 (c) simulation measures.
"""
import healpy as hp
import numpy as np

from .lfi_noise import spin2_basis

FWHM_ARCMIN = {"K": 52.8, "Ka": 39.6, "Q": 30.6, "V": 21.0, "W": 13.2, "30": 32.3, "44": 27.9, "70": 13.2, "353": 4.9}


class WindowMatcher:
    def __init__(self, nside=16, lmax=None, target_arcmin=60.0):
        self.nside = nside
        self.lmax = 3 * nside - 1 if lmax is None else lmax
        self.target = target_arcmin
        self.Y, self.ell = spin2_basis(nside, self.lmax)
        self.Yp = np.linalg.solve(self.Y.T @ self.Y, self.Y.T)          # pseudo-inverse, (nmodes, 2 npix)

    def factor(self, fwhm_arcmin):
        bt = hp.gauss_beam(np.radians(self.target / 60.0), lmax=self.lmax, pol=True)[:, 1]
        bb = hp.gauss_beam(np.radians(fwhm_arcmin / 60.0), lmax=self.lmax, pol=True)[:, 1]
        return (bt / bb)[self.ell]

    def operators(self, fwhm_arcmin):
        """S and S^-1 as dense (2 npix, 2 npix) arrays."""
        f = self.factor(fwhm_arcmin)
        n = self.Y.shape[0]
        S = np.eye(n) + (self.Y * (f - 1.0)[None, :]) @ self.Yp
        Si = np.eye(n) + (self.Y * (1.0 / f - 1.0)[None, :]) @ self.Yp
        return S, Si

    def apply(self, d, ninv, fwhm_arcmin):
        """Data and inverse covariance of the band brought to the target window."""
        S, Si = self.operators(fwhm_arcmin)
        out = Si.T @ ninv @ Si
        return S @ d, 0.5 * (out + out.T)
