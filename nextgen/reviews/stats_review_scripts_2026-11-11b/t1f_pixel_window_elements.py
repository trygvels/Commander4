"""Addition 1: where the isotropic pixel-window model differs from the exact pixel average (nside 4, modes of t1b).
Spectrum: white in l up to the band limit (the hardest case), and the LCDM spectrum stretched by 4 in l (so that the
nside-4 band limit sits where the nside-16 one does on the real spectrum)."""
import sys

import camb
import healpy as hp
import numpy as np

sys.path.insert(0, "nextgen")
from lowres.window import WindowMatcher  # noqa: E402

z = np.load("nextgen/reviews/stats_review_scripts_2026-11-11b/pixavg_modes_n4.npz")
PY, ells, kinds = z["PY"], z["ells"], z["kinds"]
WM = WindowMatcher(4, target_arcmin=240.0)
Y, ell = WM.Y, WM.ell
n = Y.shape[0]
npix = n // 2
omega = 4 * np.pi / npix
pw = hp.pixwin(4, pol=True, lmax=16)[1]
th = np.degrees(hp.pix2ang(4, np.arange(npix), nest=True)[0])
colat = np.tile(np.minimum(th, 180 - th), 2)
pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=600)
cls = camb.get_results(pars).get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]
bl = hp.gauss_beam(np.radians(240.0 / 60.0), lmax=64, pol=True)[:, 1]
nE = Y.shape[1] // 2
for label, cE, cB in (("white to l = 11", np.ones(65), np.ones(65)),
                      ("LCDM stretched: C'(l) = 16 C(4 l), window 240 arcmin", 16 * cls[4 * np.arange(65), 1] * bl ** 2, 16 * cls[4 * np.arange(65), 2] * bl ** 2)):
    cmode = np.where(kinds == 1, cE[ells], cB[ells])
    inb = ells <= WM.lmax
    C_in = (PY[:, inb] * cmode[inb][None, :]) @ PY[:, inb].T                    # exact pixel average, l <= 11
    C_all = (PY * cmode[None, :]) @ PY.T                                        # exact pixel average, l <= 64
    cl = np.concatenate([(cE[:12] * pw[:12] ** 2)[ell[:nE]], (cB[:12] * pw[:12] ** 2)[ell[nE:]]])
    C_mod = (Y * (cl / omega)[None, :]) @ Y.T
    r_in = np.diag(C_mod) / np.diag(C_in)
    r_all = np.diag(C_mod) / np.diag(C_all)
    print(f"## {label}")
    print(f"   variance per element, model / exact pixel average of the same multipoles (l <= 11): mean {r_in.mean():.3f}, "
          f"10 to 90 percent {np.percentile(r_in, 10):.3f} to {np.percentile(r_in, 90):.3f}, extremes {r_in.min():.3f} {r_in.max():.3f}; "
          f"elements off by more than 10 percent: {int(np.sum(np.abs(r_in - 1) > 0.1))} of {n}, all within {colat[np.abs(r_in - 1) > 0.1].max() if np.any(np.abs(r_in - 1) > 0.1) else 0:.0f} deg of a pole")
    print(f"   model / exact pixel average of all multipoles to 64: mean {r_all.mean():.3f} (the sky has {1 / r_all.mean():.2f} times the prior variance)")
    # how wrong is the prior as a whole: eigenvalues of C_model^-1/2 C_exact C_model^-1/2 with the code's floor
    for floor in (0.01, None):
        fl = floor ** 2 if floor else (np.trace(C_all) - np.trace(C_mod)) / n
        Lm = np.linalg.cholesky(C_mod + fl * np.eye(n))
        ev = np.linalg.eigvalsh(np.linalg.solve(Lm, np.linalg.solve(Lm, C_all).T))
        print(f"   true / prior variance along the prior's principal directions, floor variance {fl:.2e}: median {np.median(ev):.2f}, "
              f"90th percentile {np.percentile(ev, 90):.1f}, maximum {ev.max():.1f}; directions with more than twice the prior variance: {int(np.sum(ev > 2))} of {n}")
