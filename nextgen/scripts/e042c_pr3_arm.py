"""E042 (f): the K band against 30 GHz fits with the Planck 2018 (PR3) 30 GHz map in place of PR4.

Noise: the same construction as for PR4 (E038 b2): white covariance of the map file times a scale-dependent
factor fitted to the ring-half difference of the PR3 maps. Fits: K, 30, 353 and 30, Ka, 353, four BP XV
regions, nside 16, window-matched, profile likelihood with debiased sandwich errors.

    python -I nextgen/scripts/e042c_pr3_arm.py --out nextgen/experiments/E042/pr3_arm.npz
"""
import argparse
import gc
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres import lfi_noise as LN  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, lfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher  # noqa: E402

PR3 = "/Users/trygvels/data/planck_pr3_lfi/LFI_SkyMap_030-BPassCorrected_1024_R3.00_{split}.fits"
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    # noise model from the PR3 ring halves
    d1, n1 = D.load_planck_qu(PR3.format(split="full-ringhalf-1"), 16)
    d2, n2 = D.load_planck_qu(PR3.format(split="full-ringhalf-2"), 16)
    q1, u1, x1 = LN.cov_blocks(n1)
    q2, u2, x2 = LN.cov_blocks(n2)
    _, isq = LN.sqrt_blocks((q1 + q2) / 4, (u1 + u2) / 4, (x1 + x2) / 4)
    z = LN.apply_blocks(isq, (d1 - d2) / 2)
    M = LN.NoiseModel(16)
    res = M.fit(z)
    th = res.x
    print(f"# PR3 30 GHz: z^T z / n = {z @ z / z.size:.4f}; noise model r = {np.exp(th[0]):.3f}, A(l=10) = {np.exp(th[1]):.3f}, gamma = {th[2]:.3f}; "
          f"converged {res.success}; {time.time() - t0:.0f} s", flush=True)
    d, nw = D.load_planck_qu(PR3.format(split="full"), 16)
    _, isqf = LN.sqrt_blocks(*LN.cov_blocks(nw))
    Rinv = np.linalg.inv(M.R(th))
    n = d.size
    h = n // 2
    aa, bb, cc = isqf
    left = np.concatenate([aa[:, None] * Rinv[:h] + cc[:, None] * Rinv[h:], cc[:, None] * Rinv[:h] + bb[:, None] * Rinv[h:]])
    n30 = np.concatenate([left[:, :h] * aa[None, :] + left[:, h:] * cc[None, :], left[:, :h] * cc[None, :] + left[:, h:] * bb[None, :]], axis=1)
    n30 = 0.5 * (n30 + n30.T)
    del Rinv, left
    # sign and units against PR4
    d4, _ = D.load_planck_qu(D.PR4_LFI.format(freq=30, split="full"), 16)
    big = np.abs(d4) > 10
    print(f"# PR3 against PR4 30 GHz at nside 16: correlation {np.corrcoef(d, d4)[0, 1]:.4f}; slope of PR3 on PR4 over elements above 10 uK "
          f"{np.sum(d[big] * d4[big]) / np.sum(d4[big] ** 2):.4f}; rms of the difference {np.sqrt(np.mean((d - d4) ** 2)):.3f} uK", flush=True)
    WM = WindowMatcher(16, target_arcmin=60.0)
    region = D.bpxv_regions()

    def band(B, dd, ni):
        return BandData(B.name, dd, ni, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF))

    raw = {}
    d30, n30 = WM.apply(d, n30, FWHM_ARCMIN["30"])
    raw["30"] = band(lfi_band(30), d30, n30)
    for b in ("K", "Ka"):
        dd, ni = D.load_wmap(b)
        dd, ni = WM.apply(dd, ni, FWHM_ARCMIN[b])
        raw[b] = band(wmap_band(b, D.WMAP_BANDPASS), dd, ni)
    dd, ni = D.load_planck_qu(D.PR4_353, 16)
    dd, ni = WM.apply(dd, ni, FWHM_ARCMIN["353"])
    raw["353"] = band(hfi_band(353, psb=True), dd, ni)
    store = {"theta": th}
    for label, bs in (("K 30(PR3) 353", ("K", "30", "353")), ("30(PR3) Ka 353", ("30", "Ka", "353"))):
        L = RegionLikelihood([raw[b] for b in bs], region, cmb_prior_sd=0.2, use_dust=True)
        bh, r = fit_profile(L, np.full(4, -3.1))
        _, Sw, J, T, sol, ok = sandwich_plugin(L, bh)
        sd = np.sqrt(np.diag(Sw))
        rb = L.residuals(sol)
        print(f"{label:18s} " + "  ".join(f"{x:+.4f} +- {s:.4f}" for x, s in zip(bh, sd)) + f" | chi2 {2 * r.fun:.1f} | "
              + "  ".join(f"{b.name} {float(v @ b.apply(v)) / n:.3f}" for b, v in zip(L.bands, rb)) + f" | {time.time() - t0:.0f} s", flush=True)
        store[label.replace(" ", "_").replace("(", "").replace(")", "")] = np.concatenate([bh, sd])
        del L, sol
        gc.collect()
    if a.out:
        np.savez(a.out, **store)


if __name__ == "__main__":
    main()
