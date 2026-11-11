"""E040 stage 5 on the maps, second version: block errors from `lowres.block` with the diagnostics asked for by the
independent review of 2026-11-11 (blocks with negative information, cross-region share of the block Jacobian,
no leverage-corrected error above a leverage of 0.9), the block scores and informations saved for a figure,
and an option for the white floor of the LCDM prior (0.11 uK is the LCDM power above the band limit of the
nside-16 basis as seen through the pixel average, by the review's calculation).

    python -I nextgen/scripts/e040i_block_maps.py [--prior lcdm|white] [--floor-uk 0.01] --out FILE.npz
"""
import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.block import KINDS, block_errors  # noqa: E402
from lowres.fit import fit_profile  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher, lcdm_cmb_precision  # noqa: E402

NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", default="lcdm", choices=["lcdm", "white"])
    ap.add_argument("--floor-uk", type=float, default=0.01)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    WM = WindowMatcher(16, target_arcmin=60.0)
    bands = []
    for b in ("K", "Ka", "Q", "V", "353"):
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=True)
        else:
            d, ninv = D.load_wmap(b)
            B = wmap_band(b, D.WMAP_BANDPASS)
        d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        bands.append(BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF)))
    region = D.bpxv_regions()
    n = region.size
    if a.prior == "lcdm":
        prec, cov = lcdm_cmb_precision(WM, floor_uk=a.floor_uk)
        rms = float(np.sqrt(np.mean(np.diag(cov))))
        del cov
        kw = dict(cmb_prior_sd=None, cmb_prior_prec=prec)
    else:
        kw, rms = dict(cmb_prior_sd=0.2), 0.2
    L = RegionLikelihood(bands, region, use_dust=True, **kw)
    bh, res = fit_profile(L, np.full(4, -3.1))
    block = np.tile(np.arange(n // 2) // 64, 2)
    nblk = np.array([np.unique(block[region == k]).size for k in range(4)])
    out = block_errors(L, bh, block, 48)
    print(f"# E040 stage 5, maps: K, Ka, Q, V, 353, window-matched, nside 16, CMB prior {a.prior}"
          + (f" with a white floor of {a.floor_uk} uK (rms per element {rms:.3f} uK)" if a.prior == "lcdm" else "")
          + f"; converged {res.success}; chi2_min {2 * res.fun:.1f}; {time.time() - t0:.0f} s")
    print(f"# blocks with negative information per region: {out['n_neg'].tolist()}; cross-region share of the block Jacobian: {out['cross']:.4f}")
    print(f"{'region':16s} {'blocks':>6s} {'n_eff':>6s} {'max h':>6s} {'sum D/J_c':>9s} {'index':>9s} {'sandwich':>9s} " + " ".join(f"{k:>8s}" for k in KINDS)
          + "  ratio to sandwich: " + " ".join(f"{k:>6s}" for k in KINDS))
    for k_ in range(4):
        print(f"{NAMES[k_]:16s} {nblk[k_]:6d} {out['neff'][k_]:6.1f} {out['lev_max'][k_]:6.2f} {out['d_over_jc'][k_]:9.3f} {bh[k_]:+9.4f} {out['sd_model'][k_]:9.4f} "
              + " ".join(f"{out[k][k_]:8.4f}" for k in KINDS) + "                     " + " ".join(f"{out[k][k_] / out['sd_model'][k_]:6.2f}" for k in KINDS))
    if a.out:
        np.savez(a.out, beta=bh, sd_model=out["sd_model"], neff=out["neff"], lev_max=out["lev_max"], n_neg=out["n_neg"], cross=out["cross"],
                 d_over_jc=out["d_over_jc"], G=out["G"], D=out["D"], block=block, region=region, **{k: out[k] for k in KINDS})


if __name__ == "__main__":
    main()
