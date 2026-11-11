"""Referee check: how does a white CMB prior of 0.2 uK per Q/U element at nside 16 compare, mode by mode,
with the LCDM EE spectrum at the multipoles that carry the index in the low-pass fits (2 to 20)?

White prior: every Q and U element N(0, sigma^2), so C_l^EE = C_l^BB = sigma^2 * Omega_pix for all l on the
grid. LCDM: CAMB, Planck 2018 base parameters, tau = 0.054, with a 60 arcmin Gaussian window and the
nside-16 pixel window (the working window of the paper's fits). A few seconds, one thread.
"""
import numpy as np
import camb
import healpy as hp

sigma = 0.2                                   # uK per element
omega = 4 * np.pi / (12 * 16 ** 2)
cl_white = sigma ** 2 * omega                 # uK^2 sr, E and B each

pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=300)
res = camb.get_results(pars)
cls = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]   # TT, EE, BB, TE
ell = np.arange(cls.shape[0])
bl = hp.gauss_beam(np.radians(1.0), lmax=ell.max(), pol=True)[:, 1]
pw = hp.pixwin(16, pol=True, lmax=min(ell.max(), 3 * 16 - 1))[1]
w2 = bl ** 2
w2[: pw.size] *= pw ** 2
w2[pw.size:] = 0.0
ee = cls[:, 1] * w2

print(f"white prior: C_l = {cl_white:.3e} uK^2 sr for E and for B at every multipole")
print("  l   D_l^EE(LCDM, windowed)   D_l(white prior)   ratio kappa^2 = true / prior")
for l in (2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 30, 40, 47):
    d_true = l * (l + 1) * ee[l] / (2 * np.pi)
    d_white = l * (l + 1) * cl_white / (2 * np.pi)
    print(f"{l:3d}   {d_true:10.4e}              {d_white:10.4e}         {ee[l] / cl_white:8.2f}")
var_true = np.sum((2 * ell[2:48] + 1) / (4 * np.pi) * ee[2:48]) / 2      # per Q or U element, E only
print(f"rms per Q/U element from EE, l = 2 to 47, 60 arcmin and pixel window: {np.sqrt(var_true):.3f} uK (prior: {sigma})")
for L in (6, 12, 20, 47):
    sel = slice(2, L + 1)
    k2 = np.sum((2 * ell[sel] + 1) * ee[sel]) / np.sum((2 * ell[sel] + 1) * cl_white)
    print(f"mode-count-weighted kappa^2 of the E modes at l = 2 to {L}: {k2:6.2f}   (B modes: 0)")
