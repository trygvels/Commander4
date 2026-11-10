"""E042 (f): the Planck 2018 (PR3) 30 GHz low-resolution map with its released noise covariance matrix.

The released nside-16 map is the noise-weighted downgrade of the full-resolution map, smoothed in polarization
with a cosine window (1 to multipole 16, zero at 48), with a regularization noise of 0.02 uK; the released
`offset_covmat` is the covariance of that map including the regularization. Here:

  1. the transfer function of the released map against my own noise-weighted downgrade of the PR3 map is
     measured from their cross spectrum (E and B together);
  2. the map is expanded in the spin-2 harmonics of the grid, the coefficients to multipole L are kept and
     divided by that transfer function and by the ratio of the 32 arcmin beam to the common 60 arcmin window,
     and their covariance follows from the released matrix;
  3. the band enters the region likelihood through d' = Y c' and N'^-1 = Y+^T Cov(c')^-1 Y+, an inverse
     covariance of rank K that gives no weight to anything else;
  4. the covariance is tested on the ring-half difference of the PR3 maps put through the same steps;
  5. fits: K, 30, 353; 30, Ka, 353; K, Ka, 30, 353. For comparison the PR4 map with its noise model is cut to
     the same multipoles by projecting the others out.

    python -I nextgen/scripts/e042d_pr3_ncvm.py --lmax 24 --out nextgen/experiments/E042/pr3_ncvm.npz
"""
import argparse
import gc
import os
import sys
import time

import healpy as hp
import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, lfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher  # noqa: E402

X = "/Users/trygvels/data/planck_pr3_lfi_ncvm/extracted/"
PR3 = "/Users/trygvels/data/planck_pr3_lfi/LFI_SkyMap_030-BPassCorrected_1024_R3.00_{split}.fits"
T0 = time.time()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=24)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    L = a.lmax
    WM = WindowMatcher(16, target_arcmin=60.0)
    Y, Yp, ell = WM.Y, WM.Yp, WM.ell
    keep = ell <= L
    K = int(keep.sum())
    t = fits.open(X + "LFI_SkyMap_030_0016_coswin_DX12_Corrected_full_regnoise.fits")[1].data
    d_rel = np.concatenate([np.asarray(t["Q_STOKES"], float).ravel(), np.asarray(t["U_STOKES"], float).ravel()]) * 1e6
    C = np.fromfile(X + "offset_covmat_toast_nside64T16_coswin_025sec_30GHz_DX12_full_bin.dat", dtype=np.float64)
    N9 = int(np.sqrt(C.size))
    C = C.reshape(N9, N9)[3072:, 3072:] * 1e12
    C = 0.5 * (C + C.T)
    # 1. transfer function of the released map against my downgrade, per multipole, from harmonic coefficients on the grid
    mine, _ = D.load_planck_qu(PR3.format(split="full"), 16)
    c_rel, c_mine = Yp @ d_rel, Yp @ mine
    tl = np.ones(48)
    for l in range(2, 48):
        s = ell == l
        tl[l] = np.sum(c_rel[s] * c_mine[s]) / np.sum(c_mine[s] ** 2)
    cosw = np.where(np.arange(48) <= 16, 1.0, 0.5 * (1 + np.cos(np.pi * (np.arange(48) - 16) / 32.0)))
    print("# transfer function of the released map against a plain noise-weighted downgrade (cosine window in brackets): "
          + "  ".join(f"l={l} {tl[l]:.3f} ({cosw[l]:.3f})" for l in (2, 8, 16, 20, 24, 28, 32)), flush=True)
    # smooth the measured transfer function: 1 below l = 12, the cosine window times a constant fitted over 17 to L above
    kfac = np.median(tl[17:max(L, 20) + 1] / cosw[17:max(L, 20) + 1])
    tsm = np.where(np.arange(48) <= 16, 1.0, cosw * kfac)
    print(f"# adopted transfer function: 1 to l = 16, cosine window times {kfac:.3f} above", flush=True)
    b60 = hp.gauss_beam(np.radians(1.0), lmax=47, pol=True)[:, 1]
    b32 = hp.gauss_beam(np.radians(FWHM_ARCMIN["30"] / 60.0), lmax=47, pol=True)[:, 1]
    f = (b60 / b32 / tsm)[ell[keep]]
    Ypk, Yk = Yp[keep], Y[:, keep]
    cc = f * (Ypk @ d_rel)
    Cc = (f[:, None] * (Ypk @ C @ Ypk.T)) * f[None, :]
    Cci = np.linalg.inv(0.5 * (Cc + Cc.T))
    d30 = Yk @ cc
    n30 = Ypk.T @ Cci @ Ypk
    n30 = 0.5 * (n30 + n30.T)
    # 4. the released covariance against the ring-half difference put through the same steps
    h1, _ = D.load_planck_qu(PR3.format(split="full-ringhalf-1"), 16)
    h2, _ = D.load_planck_qu(PR3.format(split="full-ringhalf-2"), 16)
    ch = f * tsm[ell[keep]] * (Ypk @ ((h1 - h2) / 2))          # same window as the released map, then the same correction
    chi = float(ch @ Cci @ ch) / K
    print(f"# ring-half difference of PR3 through the same steps against the released covariance: chi2 per mode {chi:.3f} for {K} modes "
          f"(sd {np.sqrt(2.0 / K):.3f}); by multipole range: " + "  ".join(
              f"l {lo}-{hi} {float(np.mean(ch[(ell[keep] >= lo) & (ell[keep] <= hi)] ** 2 / np.diag(Cc)[(ell[keep] >= lo) & (ell[keep] <= hi)])):.2f}"
              for lo, hi in ((2, 3), (4, 6), (7, 10), (11, 16), (17, L))), flush=True)
    region = D.bpxv_regions()
    n = region.size

    def band(B, dd, ni):
        return BandData(B.name, dd, ni, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF))

    raw = {"30pr3": band(lfi_band(30), d30, n30)}
    for b in ("K", "Ka"):
        dd, ni = D.load_wmap(b)
        dd, ni = WM.apply(dd, ni, FWHM_ARCMIN[b])
        raw[b] = band(wmap_band(b, D.WMAP_BANDPASS), dd, ni)
    dd, ni = D.load_planck_qu(D.PR4_353, 16)
    dd, ni = WM.apply(dd, ni, FWHM_ARCMIN["353"])
    raw["353"] = band(hfi_band(353, psb=True), dd, ni)
    # PR4 30 GHz with its noise model, cut to the same multipoles by projecting out everything above L in the basis
    d4, n4 = D.load_lfi(30, 16)
    d4, n4 = WM.apply(d4, n4, FWHM_ARCMIN["30"])
    V = Y[:, ~keep]
    NV = n4 @ V
    n4c = n4 - NV @ np.linalg.solve(V.T @ NV, NV.T)
    raw["30pr4cut"] = band(lfi_band(30), d4, 0.5 * (n4c + n4c.T))
    store = {"tl": tl, "kfac": kfac, "chi_halfring": chi}
    print(f"# profile index and debiased sandwich error (High Latitude, Spur, Galactic Center, Plane); 30 GHz restricted to l <= {L}", flush=True)
    for label, bs in (("K 30(PR3 NCVM) 353", ("K", "30pr3", "353")), ("30(PR3 NCVM) Ka 353", ("30pr3", "Ka", "353")),
                      ("K Ka 30(PR3 NCVM) 353", ("K", "Ka", "30pr3", "353")), ("K 30(PR4, same l) 353", ("K", "30pr4cut", "353")),
                      ("30(PR4, same l) Ka 353", ("30pr4cut", "Ka", "353"))):
        Lk = RegionLikelihood([raw[b] for b in bs], region, cmb_prior_sd=0.2, use_dust=True)
        bh, r = fit_profile(Lk, np.full(4, -3.1))
        _, Sw, J, T, sol, ok = sandwich_plugin(Lk, bh)
        sd = np.sqrt(np.diag(Sw))
        rb = Lk.residuals(sol)
        print(f"{label:24s} " + "  ".join(f"{x:+.4f} +- {s:.4f}" for x, s in zip(bh, sd)) + f" | chi2 {2 * r.fun:.1f} | "
              + "  ".join(f"{b.name} {float(v @ b.apply(v)) / n:.3f}" for b, v in zip(Lk.bands, rb)) + f" | {time.time() - T0:.0f} s", flush=True)
        store[label.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")] = np.concatenate([bh, sd])
        del Lk, sol
        gc.collect()
    if a.out:
        np.savez(a.out, **store)


if __name__ == "__main__":
    main()
