"""E036 addendum: the sandwich factor of the profile likelihood, region by region, from the truth alone.

For one amplitude per pixel and Stokes parameter with known diagonal noise, the profile score for the
index has curvature u V and variance (u + 1) V per element (u = a^2 A the true amplitude S/N squared,
A = sum_b w f^2, V the variance of ln(nu_b / nu_ref) under the weights w f^2). So for a region the sd of
z = (estimate - truth) / (profile posterior sd) is predicted to be sqrt(1 + sum V / sum (u V)).
Prints the quadrature mean of that prediction over the regions of each Fisher-sd class of E036.

    python -I nextgen/scripts/e036_sandwich_factor.py [--noise-scale 1 2 4]
"""
import argparse
import os
import sys

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e031_noise_regime_fisher import real_noise  # noqa: E402
from e036_qu_amplitude_priors import CENTRE, DR1  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--noise-scale", type=float, nargs="+", default=[1.0, 2.0, 4.0])
ap.add_argument("--beta", type=float, default=-3.1); ap.add_argument("--nu-ref", type=float, default=23.0)
a = ap.parse_args()
nside, npix = 16, hp.nside2npix(16)
bands = [23.0, 30.0, 33.0, 41.0, 44.0]
dr1 = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
dr1 = hp.smoothing(hp.ud_grade(dr1, 64), fwhm=np.radians(5.0), pol=True)
a30 = hp.ud_grade(dr1[1:], nside).reshape(2 * npix)
nus = np.array([CENTRE[b] for b in bands])
rn = real_noise(nside)
sig0 = np.stack([np.stack([rn[b][1][1], rn[b][1][2]]) for b in bands], axis=2).reshape(2 * npix, len(nus))
aref = a30 * (a.nu_ref / 30.0) ** a.beta
ft = (nus / a.nu_ref) ** a.beta; L = np.log(nus / a.nu_ref)
classes = [("< 0.03", 0.0, 0.03), ("0.03 to 0.1", 0.03, 0.1), ("0.1 to 0.25", 0.1, 0.25), ("> 0.25", 0.25, np.inf)]
print("# predicted sd of z of the profile, sqrt(mean over regions of 1 + sum V / sum u V), per Fisher-sd class")
for ns in a.noise_scale:
    w = 1.0 / (sig0 * ns) ** 2
    A = (w * ft ** 2).sum(1)
    Lbar = (w * ft ** 2 * L).sum(1) / A
    V = (w * ft ** 2 * (L - Lbar[:, None]) ** 2).sum(1) / A
    u = aref ** 2 * A
    print(f"\n## noise scale {ns:g}: median u per element {np.median(u):.2f} (S/N {np.sqrt(np.median(u)):.2f})")
    for ns_r in (2, 4):
        reg = np.tile(hp.ud_grade(np.arange(hp.nside2npix(ns_r), dtype=float), nside).astype(int), 2)
        nreg = hp.nside2npix(ns_r)
        J = np.bincount(reg, weights=u * V, minlength=nreg)            # curvature = Fisher information
        sV = np.bincount(reg, weights=V, minlength=nreg)
        fsd = 1.0 / np.sqrt(J)
        fac2 = 1.0 + sV / J
        out = []
        for cname, lo, hi in classes:
            sel = (fsd >= lo) & (fsd < hi)
            if sel.any():
                out.append(f"{cname}: {np.sqrt(fac2[sel].mean()):.3f} ({int(sel.sum())} regions, range "
                           f"{np.sqrt(fac2[sel].min()):.2f} to {np.sqrt(fac2[sel].max()):.2f})")
        print(f"nside {ns_r}: " + " | ".join(out))
