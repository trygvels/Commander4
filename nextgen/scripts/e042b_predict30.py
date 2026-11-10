"""E042 (g): predict one band from a fit to the others and study observed minus predicted.

    python -I nextgen/scripts/e042b_predict30.py --target 30 --fit K Ka Q V 353 --out FILE.npz
    python -I nextgen/scripts/e042b_predict30.py --target Ka --fit K Q V 353 --out FILE.npz      (control)

The fit is the E040 profile fit (four BP XV regions, nside 16, window-matched, CMB prior 0.2 uK). The target
band is not in the fit. D = d_target - (s_target(beta_hat) a_s + f_target a_d + a_c). Printed: per region,
the weighted regression slope of D on the predicted synchrotron with its error under the target's noise
model (prediction noise not included, so the error is a lower bound); chi-square of D per element; E and B
power of the whitened D by multipole bin against the white expectation of the noise model, before and after
removing the per-region slopes.
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D_  # noqa: E402
from lowres.fit import fit_profile  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, lfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher  # noqa: E402

NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
BINS = ((2, 3), (4, 6), (7, 10), (11, 20), (21, 29), (30, 47))


def load(b, WM):
    if b == "353":
        d, ninv = D_.load_planck_qu(D_.PR4_353, 16)
        B = hfi_band(353, psb=True)
    elif b in ("30", "44", "70"):
        d, ninv = D_.load_lfi(int(b), 16)
        B = lfi_band(int(b))
    else:
        d, ninv = D_.load_wmap(b)
        B = wmap_band(b, D_.WMAP_BANDPASS)
    d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
    return B, d, ninv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="30")
    ap.add_argument("--fit", nargs="+", default=["K", "Ka", "Q", "V", "353"])
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    WM = WindowMatcher(16, target_arcmin=60.0)
    region = D_.bpxv_regions()
    n = region.size
    bands = []
    for b in a.fit:
        B, d, ninv = load(b, WM)
        bands.append(BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D_.DUST_BETA, D_.DUST_T, D_.DUST_REF)))
    L = RegionLikelihood(bands, region, cmb_prior_sd=0.2, use_dust=True)
    bh, res = fit_profile(L, np.full(4, -3.1))
    sol = L.solve(bh)
    x = sol["x"]
    a_s, a_d, a_c = x[:n], x[n:2 * n], x[2 * n:3 * n]
    del L, sol
    Bt, dt, nt = load(a.target, WM)
    st = Bt.power_law(bh, 23.0)[0][region]
    ft = Bt.mbb(D_.DUST_BETA, D_.DUST_T, D_.DUST_REF)
    psync = st * a_s
    pred = psync + ft * a_d + a_c
    Dm = dt - pred
    print(f"# E042 (g): target {a.target} predicted from {a.fit}; indices {np.round(bh, 4)}; {time.time() - t0:.0f} s")
    print(f"# rms of the target {np.sqrt(np.mean(dt ** 2)):.2f}, of the prediction {np.sqrt(np.mean(pred ** 2)):.2f}, of D {np.sqrt(np.mean(Dm ** 2)):.2f} uK; "
          f"chi2 of D per element under the target's noise model {Dm @ nt @ Dm / n:.3f}")
    # per-region slope of D on the predicted synchrotron, generalized least squares with the target's N^-1, all regions jointly
    X = np.zeros((n, 4))
    X[np.arange(n), region] = psync
    F = X.T @ nt @ X
    slope = np.linalg.solve(F, X.T @ nt @ Dm)
    sd = np.sqrt(np.diag(np.linalg.inv(F)))
    print("# slope of D on the predicted synchrotron (a fractional excess of the target), per region, with its error from the target's noise alone:")
    for k in range(4):
        print(f"#   {NAMES[k]:16s} {slope[k]:+.4f} +- {sd[k]:.4f}   (as an index shift against the lowest fitted band: see record)")
    # one global slope
    den = float(psync @ nt @ psync)
    sg = float(psync @ nt @ Dm) / den
    print(f"#   one slope for the sky {sg:+.4f} +- {1 / np.sqrt(den):.4f}")
    D2 = Dm - X @ slope
    print(f"# chi2 per element after removing the four slopes: {D2 @ nt @ D2 / n:.3f}")
    # power by multipole of D whitened with the symmetric square root of N^-1
    w, U = np.linalg.eigh(nt)
    Wh = (U * np.sqrt(np.clip(w, 0, None))[None, :]) @ U.T
    Y, ell = WM.Y, WM.ell
    nE = Y.shape[1] // 2
    isE = np.arange(Y.shape[1]) < nE
    G = np.einsum("ij,ij->j", Y, Y)
    print("# power of the whitened D by multipole bin over the white expectation: E, B before removing slopes | E, B after")
    for lo, hi in BINS:
        sel = (ell >= lo) & (ell <= hi)
        row = []
        for v in (Dm, D2):
            c = Y.T @ (Wh @ v)
            row += [np.mean(c[sel & isE] ** 2 / G[sel & isE]), np.mean(c[sel & ~isE] ** 2 / G[sel & ~isE])]
        nm = int((sel & isE).sum())
        print(f"#   l {lo:2d} to {hi:2d} ({nm:4d} modes each): {row[0]:6.2f} {row[1]:6.2f} | {row[2]:6.2f} {row[3]:6.2f}   (sd of a white ratio {np.sqrt(2.0 / nm):.2f})")
    if a.out:
        np.savez(a.out, beta=bh, D=Dm, D2=D2, pred=pred, psync=psync, slope=slope, slope_sd=sd, target=dt, region=region)


if __name__ == "__main__":
    main()
