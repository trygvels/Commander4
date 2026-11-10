"""E038 (b2): fit and test the scale-dependent LFI noise model on the ring-half difference, and store
the fitted parameters.

    python -I nextgen/scripts/e038_lfi_noise_model.py --freqs 30 44 70 --out nextgen/experiments/E038/lfi_noise_model.npz
"""
import argparse
import os
import sys
import time

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.lfi_noise import NoiseModel, whitened_half_difference  # noqa: E402

BINS = ((2, 3), (4, 6), (7, 10), (11, 20), (21, 29), (30, 47))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freqs", nargs="+", type=int, default=[30, 44, 70])
    ap.add_argument("--nsim", type=int, default=200)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    M = NoiseModel(16)
    npix = hp.nside2npix(16)
    print(f"# basis: {M.Y.shape[1]} modes on {M.Y.shape[0]} elements; Y^T Y diagonal {np.diag(M.G).min():.3f} to "
          f"{np.diag(M.G).max():.3f}, largest off-diagonal {np.abs(M.G - np.diag(np.diag(M.G))).max():.3f}; {time.time() - t0:.0f} s")
    lat = np.abs(hp.pix2ang(16, np.arange(npix), nest=True, lonlat=True)[1])
    zone = np.tile(np.digitize(lat, [10, 20, 30, 45, 60]), 2)
    reg = D.bpxv_regions()
    nE = M.Y.shape[1] // 2
    isE = np.arange(M.Y.shape[1]) < nE
    rng = np.random.default_rng(3802)
    store = {}
    for f in a.freqs:
        z = whitened_half_difference(f)
        res = M.fit(z)
        th = res.x
        r, A, g = np.exp(th[0]), np.exp(th[1]), th[2]
        print(f"## LFI {f} GHz: z^T z / n = {z @ z / z.size:.4f}; fit r = {r:.3f}, A(l=10) = {A:.3f}, gamma = {g:.3f}; "
              f"-2 ln L = {res.fun:.1f} against {M.neg2loglike(np.array([0.0, -20.0, 1.0]), z):.1f} for white; "
              f"{res.nit} iterations, converged {res.success}; {time.time() - t0:.0f} s")
        R = M.R(th)
        L = np.linalg.cholesky(R)
        Rinv = np.linalg.inv(R)
        print(f"   z^T R^-1 z / n = {z @ Rinv @ z / z.size:.4f}; model variance per element over white: median "
              f"{np.median(np.diag(R)):.3f}; multipole 2, 5, 10, 20, 40 variance over white: "
              + ", ".join(f"{r + A * (l / 10.0) ** (-g):.2f}" for l in (2, 5, 10, 20, 40)))

        def stats(v):
            u = v * (Rinv @ v)
            c = M.Y.T @ v
            out = [np.sum(u[zone == k]) for k in range(6)] + [np.sum(u[reg == k]) for k in range(4)]
            for lo, hi in BINS:
                sel = (M.ell >= lo) & (M.ell <= hi)
                out += [np.mean(c[sel & isE] ** 2), np.mean(c[sel & ~isE] ** 2)]
            return np.array(out)

        obs = stats(z)
        sims = np.array([stats(L @ rng.standard_normal(z.size)) for _ in range(a.nsim)])
        dev = (obs - sims.mean(0)) / sims.std(0)
        names = [f"zone{k}" for k in range(6)] + ["HighLat", "Spur", "GC", "Plane"] + \
                [f"{x}{lo}-{hi}" for lo, hi in BINS for x in ("E", "B")]
        print("   deviation of the real half difference from model simulations, in sd:")
        print("   " + "  ".join(f"{nm} {dv:+.2f}" for nm, dv in zip(names, dev)))
        print(f"   cells beyond 3 sd: {int(np.sum(np.abs(dev) > 3))} of {dev.size}; beyond 2 sd: {int(np.sum(np.abs(dev) > 2))}; "
              f"rms deviation {np.sqrt(np.mean(dev ** 2)):.2f}; {time.time() - t0:.0f} s")
        store[f"theta_{f}"] = th
    if a.out:
        np.savez(a.out, **store)


if __name__ == "__main__":
    main()
