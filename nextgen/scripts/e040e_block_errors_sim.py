"""E040 stage 5, check D1: the block-robust error on simulations where the model is true.

The E039 nside-8 simulations (four BP XV regions, K, Ka, Q, V, 353; one true index per region). For each
realization: profile optimum, debiased sandwich error, and the block-robust error from the scatter of the
score over the 48 nside-2 pixels, as in `e040d_block_errors.py`. Under the model the ratio is 1.

    python -I nextgen/scripts/e040e_block_errors_sim.py --nreal 60
"""
import argparse
import os
import sys
import time

import healpy as hp
import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres.data import (DUST_BETA, DUST_REF, DUST_T, PR4_353, WMAP_BANDPASS, bpxv_regions, degrade,  # noqa: E402
                         load_planck_qu, load_wmap)
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band, g_rj, wmap_band  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
CMB_SD = 0.2
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nreal", type=int, default=60)
    ap.add_argument("--seed", type=int, default=39)
    ap.add_argument("--spread", type=float, default=0.0,
                    help="sd of a true index offset per nside-2 block inside each region (0: the model is true)")
    a = ap.parse_args()
    ns = 8
    n = 2 * hp.nside2npix(ns)
    raw = []
    for b in ("K", "Ka", "Q", "V"):
        d, ninv = load_wmap(b)
        d, ninv = degrade(d, ninv, 16, ns)
        raw.append((f"WMAP_{b}", wmap_band(b, WMAP_BANDPASS), d, ninv))
    d3, n3 = load_planck_qu(PR4_353, ns)
    raw.append(("Planck_353", Band.delta("Planck_353", 353.0), d3, n3))
    reg16 = bpxv_regions()
    r = reg16[:reg16.size // 2].reshape(-1, 4)
    region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in r]), 2)
    block = np.tile(np.arange(n // 2) // 16, 2)                 # nside-2 pixel of each nside-8 pixel (NESTED)
    nblk = np.array([np.unique(block[region == k]).size for k in range(4)])
    beta_reg = np.array([-3.15, -3.25, -2.85, -2.95])
    a30 = truth_sync_30(ns)
    a_d = d3 * g_rj(353.0)
    fd = [B.mbb(DUST_BETA, DUST_T, DUST_REF) for _, B, _, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]

    def simulate(r_):
        """With a spread, every realization draws its own index offset per nside-2 block (a random effect), so
        that the scatter of the estimate about the region mean includes it, which is what the block error
        estimates. The noise and CMB draws are those of E039."""
        off = a.spread * np.random.default_rng([a.seed, 3, r_]).standard_normal(48)
        beta_el = beta_reg[region] + off[block]                 # true index per element
        a_s = a30 * (23.0 / 30.0) ** beta_el
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = CMB_SD * rng.standard_normal(n)
        return [B.power_law(beta_el, 23.0)[0] * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False)
                for (_, B, _, _), f, L in zip(raw, fd, chol)]

    print(f"# E040 stage 5 check: nside 8, blocks per region {nblk.tolist()}, true index spread between blocks {a.spread}", flush=True)
    t0 = time.time()
    B_, Sm, Sb = [], [], []
    for r_ in range(a.nreal):
        dv = simulate(r_)
        bands = [BandData(nm, d_, ninv, (lambda be, B=B: B.power_law(be, 23.0)), f) for (nm, B, _, ninv), d_, f in zip(raw, dv, fd)]
        L = RegionLikelihood(bands, region, cmb_prior_sd=CMB_SD)
        bh, _ = fit_profile(L, beta_reg)
        _, cs, Jp, T, sol, ok = sandwich_plugin(L, bh)
        g = np.zeros(n)
        for i, (b, res) in enumerate(zip(L.bands, L.residuals(sol))):
            g -= b.apply(res) * sol["ds"][i] * sol["x"][:n]
        G = np.zeros((48, 4))
        np.add.at(G, (block, region), g)
        Ji = np.linalg.inv(Jp - T)
        cb = Ji @ (G.T @ G) @ Ji
        B_.append(bh); Sm.append(np.sqrt(np.diag(cs))); Sb.append(np.sqrt(np.diag(cb)))
        if r_ % 10 == 0:
            print("real", r_, "index", " ".join(f"{x:+.4f}" for x in bh), "| block / sandwich", " ".join(f"{x:.2f}" for x in Sb[-1] / Sm[-1]),
                  f"[{time.time() - t0:.0f} s]", flush=True)
    B_, Sm, Sb = np.array(B_), np.array(Sm), np.array(Sb)
    # the index the region fit should return when the index varies: the information-weighted mean is not computed here;
    # the scatter of the estimate over noise realizations is what both errors are meant to match
    print(f"## over {len(B_)} realizations")
    print(f"{'region':16s} {'blocks':>6s} {'scatter':>8s} {'sandwich':>9s} {'block':>8s} {'block/sand median':>18s} {'16% - 84%':>14s} {'expected sd of ratio':>21s}")
    for k in range(4):
        q = np.percentile(Sb[:, k] / Sm[:, k], [16, 50, 84])
        print(f"{NAMES[k]:16s} {nblk[k]:6d} {B_[:, k].std(ddof=1):8.4f} {np.median(Sm[:, k]):9.4f} {np.median(Sb[:, k]):8.4f} {q[1]:18.3f} "
              f"{q[0]:6.2f} - {q[2]:.2f} {1 / np.sqrt(2 * max(nblk[k] - 1, 1)):21.2f}")
    print("mean estimate minus region truth: " + " ".join(f"{x:+.4f}" for x in B_.mean(0) - beta_reg))


if __name__ == "__main__":
    main()
