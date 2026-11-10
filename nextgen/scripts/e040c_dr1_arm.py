"""E040, stage 2 arm: the same region fits on the Cosmoglobe DR1 reprocessed maps (WMAP K to V and LFI
30 to 70 GHz from the DR1 release, PR4 353 GHz for dust).

The DR1 frequency maps carry a per-pixel white-noise covariance only (Q_RMS, U_RMS and a QU term at
nside 512), so the errors of this arm are optimistic at large scales and the arm is for the point estimates:
does the index between K band and 30 GHz change when both maps come from the joint reprocessing?

    python -I nextgen/scripts/e040c_dr1_arm.py --bands K 30 353 [--match-window 60] [--noise-scale 1.0]
"""
import argparse
import os
import sys
import time

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, lfi_band, wmap_band  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1"
FILES = {"K": ("CG_023-WMAP_K_IQU_n0512_v1.fits", 1e3), "Ka": ("CG_030-WMAP_Ka_IQU_n0512_v1.fits", 1e3),
         "30": ("CG_030_IQU_n0512_v1.fits", 1.0), "44": ("CG_044_IQU_n0512_v1.fits", 1.0)}
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def load_dr1(band, nside_out=16):
    """Q/U of a DR1 frequency map at nside_out with its white covariance (QQ, UU and QU per pixel), by the
    inverse-covariance weighted average used for the Planck maps. Returns d (uK) and the dense N^-1."""
    fn, scale = FILES[band]
    q, u, qr, ur, quc = [np.asarray(m, dtype=np.float64) for m in
                         hp.read_map(f"{DR1}/{fn}", field=(1, 2, 4, 5, 6), nest=True, dtype=np.float64)]
    q, u = q * scale, u * scale
    qq, uu, qu = (qr * scale) ** 2, (ur * scale) ** 2, quc * scale ** 2
    bad = ~np.isfinite(q) | (qq <= 0) | (uu <= 0)
    det = qq * uu - qu ** 2
    bad |= det <= 0
    wqq, wuu, wqu = uu / det, qq / det, -qu / det
    for w in (wqq, wuu, wqu):
        w[bad] = 0.0
    q[bad], u[bad] = 0.0, 0.0
    k = (hp.npix2nside(q.size) // nside_out) ** 2
    n = hp.nside2npix(nside_out)
    red = lambda x: x.reshape(n, k).sum(1)  # noqa: E731
    Wqq, Wuu, Wqu = red(wqq), red(wuu), red(wqu)
    bq, bu = red(wqq * q + wqu * u), red(wqu * q + wuu * u)
    Dd = Wqq * Wuu - Wqu ** 2
    dq, du = (Wuu * bq - Wqu * bu) / Dd, (Wqq * bu - Wqu * bq) / Dd
    ninv = np.zeros((2 * n, 2 * n))
    i = np.arange(n)
    ninv[i, i], ninv[n + i, n + i], ninv[i, n + i], ninv[n + i, i] = Wqq, Wuu, Wqu, Wqu
    return np.concatenate([dq, du]), ninv, float(bad.mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["K", "30", "353"])
    ap.add_argument("--match-window", type=float, default=0.0)
    ap.add_argument("--noise-scale", type=float, default=1.0, help="factor on the DR1 white-noise sd")
    ap.add_argument("--cmb-prior-sd", type=float, default=0.2)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    WM = None
    if a.match_window > 0:
        from lowres.window import FWHM_ARCMIN, WindowMatcher
        WM = WindowMatcher(16, target_arcmin=a.match_window)
    bands = []
    for b in a.bands:
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=True)
        else:
            d, ninv, fbad = load_dr1(b)
            ninv = ninv / a.noise_scale ** 2
            B = lfi_band(int(b)) if b in ("30", "44") else wmap_band(b, D.WMAP_BANDPASS)
            print(f"# DR1 {b}: fraction of bad pixels {fbad:.4f}; white sd per element, median "
                  f"{np.median(1 / np.sqrt(np.diag(ninv))):.3f} uK")
        if WM is not None:
            d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        bands.append(BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF)))
    region = D.bpxv_regions()
    n, nreg = region.size, 4
    L = RegionLikelihood(bands, region, cmb_prior_sd=a.cmb_prior_sd, use_dust=("353" in a.bands))
    bh, res = fit_profile(L, np.full(nreg, -3.1))
    Ji, Sw, J, T, sol, ok = sandwich_plugin(L, bh)
    sdn, sds = np.sqrt(np.diag(Ji)), np.sqrt(np.diag(Sw))
    rb = L.residuals(sol)
    print(f"# E040 DR1 arm: bands {a.bands}, window matched to {a.match_window}, DR1 noise scale {a.noise_scale}; "
          f"{res.nit} iterations, converged {res.success}; chi2_min {2 * res.fun:.1f}; {time.time() - t0:.0f} s")
    print("# chi2 per band r^T N^-1 r / n: " + "  ".join(f"{b.name} {float(r @ b.apply(r)) / n:.3f}" for b, r in zip(L.bands, rb)))
    print(f"{'region':16s} {'profile':>9s} {'sd naive':>9s} {'sd sand':>9s}")
    for r in range(nreg):
        print(f"{NAMES[r]:16s} {bh[r]:+9.4f} {sdn[r]:9.4f} {sds[r]:9.4f}")
    if a.out:
        np.savez(a.out, beta=bh, sd_naive=sdn, sd_sand=sds, J=J, T=T, ok=ok)


if __name__ == "__main__":
    main()
