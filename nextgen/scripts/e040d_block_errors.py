"""E040: errors of the profile index that do not assume the model fits, from the scatter of the score over
sky blocks.

At the profile optimum the score of region r is a sum over elements, sum_p g_p = 0, with
g_p = -sum_b (N_b^-1 r_b)_p (ds_b a_s)_p. Under the model its covariance is J + T. If the model misfits (the
index varies inside a region, residual systematics), the scatter of the score between blocks of sky says
so. With blocks B (the nside-2 pixels, 48 over the sky, both Stokes parameters of a pixel in one block),
G_{B,r} = sum_{p in B and r} g_p, and K_block = sum_B G_B G_B^T, the block-robust covariance is
J_c^-1 K_block J_c^-1. Blocks are treated as independent, which the dense noise covariance violates
mildly, and a region needs many blocks for the estimate to mean anything (the count is printed).

    python -I nextgen/scripts/e040d_block_errors.py [--bands K Ka Q V 353] [--cmb-prior-sd 0.2] [--out FILE.npz]
"""
import argparse
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

NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["K", "Ka", "Q", "V", "353"])
    ap.add_argument("--cmb-prior-sd", type=float, default=0.2)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    WM = WindowMatcher(16, target_arcmin=60.0)
    bands = []
    for b in a.bands:
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
        bands.append(BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF)))
    region = D.bpxv_regions()
    n = region.size
    L = RegionLikelihood(bands, region, cmb_prior_sd=(a.cmb_prior_sd if a.cmb_prior_sd > 0 else None), use_dust=True)
    bh, res = fit_profile(L, np.full(4, -3.1))
    cov_c, cov_s, Jp, T, sol, ok = sandwich_plugin(L, bh)
    a_s = sol["x"][:n]
    g = np.zeros(n)
    for i, (b, r) in enumerate(zip(L.bands, L.residuals(sol))):
        g -= b.apply(r) * sol["ds"][i] * a_s
    block = np.tile(np.arange(n // 2) // 64, 2)                 # nside-2 pixel of each nside-16 pixel (NESTED), Q and U together
    G = np.zeros((48, 4))
    np.add.at(G, (block, region), g)
    nblk = np.array([np.unique(block[region == k]).size for k in range(4)])
    K = G.T @ G
    Jc = Jp - T
    Ji = np.linalg.inv(Jc)
    cov_b = Ji @ K @ Ji
    sd_m, sd_b = np.sqrt(np.diag(cov_s)), np.sqrt(np.diag(cov_b))
    print(f"# E040 block-robust errors: bands {a.bands}, CMB prior sd {a.cmb_prior_sd}; {time.time() - t0:.0f} s")
    print("# sum of the score over all elements per region (should be zero): " + " ".join(f"{x:+.3e}" for x in G.sum(0)))
    print("# score variance: model J(a_hat) diagonal " + " ".join(f"{x:.1f}" for x in np.diag(Jp)) + "; from blocks " + " ".join(f"{x:.1f}" for x in np.diag(K)))
    print(f"{'region':16s} {'blocks':>6s} {'index':>9s} {'sd model':>9s} {'sd blocks':>9s} {'ratio':>6s}")
    for k in range(4):
        print(f"{NAMES[k]:16s} {nblk[k]:6d} {bh[k]:+9.4f} {sd_m[k]:9.4f} {sd_b[k]:9.4f} {sd_b[k] / sd_m[k]:6.2f}")
    if a.out:
        np.savez(a.out, beta=bh, sd_model=sd_m, sd_block=sd_b, G=G, K=K, J=Jp, T=T, nblk=nblk)


if __name__ == "__main__":
    main()
