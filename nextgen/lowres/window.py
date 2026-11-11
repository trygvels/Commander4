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



def lcdm_cmb_precision(matcher, floor_uk=0.01, tau=0.0544):
    """Prior precision (2 npix, 2 npix) for CMB Q/U amplitudes on the grid from the LCDM E and B spectra
    (CAMB, Planck 2018 base parameters), seen through the matcher's target Gaussian window and the pixel
    window, plus a white floor that keeps the matrix invertible outside the harmonic basis. With the basis
    columns y = Y_real sqrt(Omega_pix), the covariance of a field with spectrum C_l is Y diag(C_l / Omega_pix) Y^T.
    Returns the precision and the covariance."""
    import camb
    lmax = matcher.lmax
    pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=tau, As=2.1e-9, ns=0.9649, lmax=max(300, lmax))
    cls = camb.get_results(pars).get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]
    bl = hp.gauss_beam(np.radians(matcher.target / 60.0), lmax=lmax, pol=True)[:, 1]
    pw = hp.pixwin(matcher.nside, pol=True, lmax=lmax)[1]
    w2 = (bl * pw) ** 2
    npix = hp.nside2npix(matcher.nside)
    omega = 4 * np.pi / npix
    nE = matcher.Y.shape[1] // 2
    cl = np.concatenate([(cls[: lmax + 1, 1] * w2)[matcher.ell[:nE]], (cls[: lmax + 1, 2] * w2)[matcher.ell[nE:]]])
    cov = (matcher.Y * (cl / omega)[None, :]) @ matcher.Y.T + floor_uk ** 2 * np.eye(matcher.Y.shape[0])
    cov = 0.5 * (cov + cov.T)
    return np.linalg.inv(cov), cov
