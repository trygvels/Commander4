"""E042, low-pass arm: region fits in which the index is read from multipoles 2 to L only.

Motivation (E042 negative control, 2026-11-11): with large scales projected out, the K, Ka index of the
bright regions flattens by up to 0.37, and every amplitude ratio to K band rises by 6 to 10 percent at
multipoles 35 to 47. The pixel scale of the nside-16 maps is where the window matching cannot work
(structure above multipole 47 aliases differently for different beams). This arm asks the opposite
question: how much do the indices change when the pixel scale does not inform them?

One band of each set, the anchor (the lowest frequency), keeps all its modes, and so does 353 GHz: they pin
the synchrotron and dust amplitudes in every pixel. Every other band is expanded in the spin-2 harmonics of
the grid and enters only through its coefficients at multipoles 2 to L with their covariance Y+ N Y+^T:
d' = Y c, N'^-1 = Y+^T Cov(c)^-1 Y+, which gives no weight to anything else. The index is then informed by
those multipoles alone.

A first version low-passed every band and held the unconstrained amplitude directions with a weak prior. It
was ill-posed: with per-pixel amplitudes and band-limited data any index fits (T / J_c above 1.6 in every
region, two regions on the bounds; kept as `lowpass20_first_ill_posed.txt`).

    python -I nextgen/scripts/e042e_lowpass.py --lmax 20 --out nextgen/experiments/E042/lowpass20.npz
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
    ap.add_argument("--cmb-prior-lcdm", action="store_true", help="LCDM pixel covariance as the CMB prior instead of white 0.2 uK")
    a = ap.parse_args()
    WM = WindowMatcher(16, target_arcmin=60.0)
    kw = dict(cmb_prior_sd=0.2)
    if a.cmb_prior_lcdm:
        from lowres.window import lcdm_cmb_precision
        kw = dict(cmb_prior_sd=None, cmb_prior_prec=lcdm_cmb_precision(WM)[0])
    keep = WM.ell <= a.lmax
    Yk, Ypk = WM.Y[:, keep], WM.Yp[keep]
    region = D.bpxv_regions()
    n = region.size
    need = sorted({b for s in a.sets for b in s.split(",")})
    full, low = {}, {}
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
        mk = lambda dd, ni, B=B: BandData(B.name, dd, ni, (lambda be, B=B: B.power_law(be, 23.0)),  # noqa: E731
                                          B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF))
        full[b] = mk(d, ninv)
        if b != "353":
            cov = np.linalg.inv(ninv)
            Cc = Ypk @ cov @ Ypk.T
            Cci = np.linalg.inv(0.5 * (Cc + Cc.T))
            nlp = Ypk.T @ Cci @ Ypk
            low[b] = mk(Yk @ (Ypk @ d), 0.5 * (nlp + nlp.T))
            del cov, Cc, Cci, nlp
            gc.collect()
        print(f"# {B.name}: loaded; low-passed copy to l <= {a.lmax} ({int(keep.sum())} modes); {time.time() - T0:.0f} s", flush=True)
    store = {}
    print(f"# CMB prior: {'LCDM covariance' if a.cmb_prior_lcdm else 'white 0.2 uK'}", flush=True)
    print(f"# E042 low-pass arm: the first band of each set and 353 GHz keep all modes; the others enter through multipoles 2 to {a.lmax}; "
          "profile index and debiased sandwich error (High Latitude, Spur, Galactic Center, Plane)", flush=True)
    for s in a.sets:
        bs = s.split(",")
        bands = [full[bs[0]]] + [low[b] if b != "353" else full[b] for b in bs[1:]]
        L = RegionLikelihood(bands, region, use_dust=True, **kw)
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
