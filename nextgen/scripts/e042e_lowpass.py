"""E042, low-pass arm: region fits that use only the multipoles 2 to L of every band.

Motivation (E042 negative control, 2026-11-11): with large scales projected out, the K, Ka index of the
bright regions flattens by up to 0.23, and every amplitude ratio to K band rises by 6 to 10 percent at
multipoles 35 to 47. The pixel scale of the nside-16 maps is where the window matching cannot work
(structure above multipole 47 aliases differently for different beams). This arm asks the opposite
question: how much do the indices change when the pixel scale is not used at all?

Each band (after window matching) is expanded in the spin-2 harmonics of the grid; the coefficients at
multipoles 2 to L are kept with their covariance Y+ N Y+^T, and the band enters through d' = Y c and
N'^-1 = Y+^T Cov(c)^-1 Y+, which gives no weight to anything else. The amplitude system is then rank
deficient, so two fictitious bands with weight 1e-8 uK^-2 hold the unconstrained synchrotron and dust
directions at zero; their effect on the constrained ones is of relative order 1e-8.

    python -I nextgen/scripts/e042e_lowpass.py --lmax 30 --out nextgen/experiments/E042/lowpass30.npz
"""
import argparse
import gc
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, lfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher  # noqa: E402

T0 = time.time()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=30)
    ap.add_argument("--sets", nargs="+", default=["K,Ka,Q,V,353", "K,Ka,353", "K,30,353", "30,Ka,353"])
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    WM = WindowMatcher(16, target_arcmin=60.0)
    keep = WM.ell <= a.lmax
    Yk, Ypk = WM.Y[:, keep], WM.Yp[keep]
    region = D.bpxv_regions()
    n = region.size
    need = sorted({b for s in a.sets for b in s.split(",")})
    raw = {}
    for b in need:
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=True)
        elif b in ("30", "44", "70"):
            d, ninv = D.load_lfi(int(b), 16)
            B = lfi_band(int(b))
        else:
            d, ninv = D.load_wmap(b)
            B = wmap_band(b, D.WMAP_BANDPASS)
        d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        cov = np.linalg.inv(ninv)
        Cc = Ypk @ cov @ Ypk.T
        Cci = np.linalg.inv(0.5 * (Cc + Cc.T))
        nlp = Ypk.T @ Cci @ Ypk
        raw[b] = BandData(B.name, Yk @ (Ypk @ d), 0.5 * (nlp + nlp.T), (lambda be, B=B: B.power_law(be, 23.0)),
                          B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF))
        del cov, Cc, Cci, nlp
        gc.collect()
        print(f"# {B.name}: low-passed to l <= {a.lmax} ({int(keep.sum())} modes); {time.time() - T0:.0f} s", flush=True)
    eps = 1e-8
    zero = np.zeros(n)
    fake_s = BandData("reg_sync", zero, np.full(n, eps), (lambda be: (np.ones(np.size(be)), np.zeros(np.size(be)))), 0.0, cmb=0.0)
    fake_d = BandData("reg_dust", zero, np.full(n, eps), (lambda be: (np.zeros(np.size(be)), np.zeros(np.size(be)))), 1.0, cmb=0.0)
    store = {}
    print(f"# E042 low-pass arm: every band restricted to multipoles 2 to {a.lmax}; profile index and debiased sandwich error "
          "(High Latitude, Spur, Galactic Center, Plane)", flush=True)
    for s in a.sets:
        bs = s.split(",")
        L = RegionLikelihood([raw[b] for b in bs] + [fake_s, fake_d], region, cmb_prior_sd=0.2, use_dust=True)
        bh, res = fit_profile(L, np.full(4, -3.1))
        _, Sw, J, T, sol, ok = sandwich_plugin(L, bh)
        sd = np.sqrt(np.diag(Sw))
        print(f"{s:16s} " + "  ".join(f"{x:+.4f} +- {e:.4f}" for x, e in zip(bh, sd)) + f" | chi2 {2 * res.fun:.1f} | T/J_c "
              + " ".join(f"{T[k, k] / max(J[k, k] - T[k, k], 1e-30):.3f}" for k in range(4)) + f" | converged {res.success} | {time.time() - T0:.0f} s", flush=True)
        store[s.replace(",", "_")] = np.concatenate([bh, sd])
        del L, sol
        gc.collect()
    if a.out:
        np.savez(a.out, **store)


if __name__ == "__main__":
    main()
