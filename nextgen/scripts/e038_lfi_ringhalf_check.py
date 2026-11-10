"""E038 (b): ring-half difference of the PR4 LFI Q/U maps at nside 16 against the white-noise covariance
of the map files, and a sign and units check against WMAP K band.

    python -I nextgen/scripts/e038_lfi_ringhalf_check.py [--freqs 30 44 70] [--nside 16]
"""
import argparse
import os
import sys

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.spectra import g_rj, lfi_band, wmap_band  # noqa: E402


def blocks(ninv, n):
    i = np.arange(n)
    return ninv[i, i], ninv[n + i, n + i], ninv[i, n + i]


def cov_blocks(ninv, n):
    a, b, c = blocks(ninv, n)
    det = a * b - c * c
    return b / det, a / det, -c / det          # QQ, UU, QU of the covariance


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freqs", nargs="+", type=int, default=[30, 44, 70])
    ap.add_argument("--nside", type=int, default=16)
    a = ap.parse_args()
    ns = a.nside
    n = hp.nside2npix(ns)
    lmax = 3 * ns - 1
    rng = np.random.default_rng(38)
    for f in a.freqs:
        d1, n1 = D.load_planck_qu(D.PR4_LFI.format(freq=f, split="full-ringhalf-1"), ns)
        d2, n2 = D.load_planck_qu(D.PR4_LFI.format(freq=f, split="full-ringhalf-2"), ns)
        df, nf = D.load_planck_qu(D.PR4_LFI.format(freq=f, split="full"), ns)
        q1, u1, x1 = cov_blocks(n1, n)
        q2, u2, x2 = cov_blocks(n2, n)
        cq, cu, cx = (q1 + q2) / 4, (u1 + u2) / 4, (x1 + x2) / 4
        dq, du = (d1[:n] - d2[:n]) / 2, (d1[n:] - d2[n:]) / 2
        det = cq * cu - cx * cx
        chi = (cu * dq ** 2 + cq * du ** 2 - 2 * cx * dq * du) / det
        qf, uf, xf = cov_blocks(nf, n)
        print(f"## LFI {f} GHz, nside {ns}")
        print(f"white-noise sd per element, full map: median {np.sqrt(np.median(np.r_[qf, uf])):.3f} uK_CMB "
              f"({np.sqrt(np.median(np.r_[qf, uf])) * g_rj({30: 28.4, 44: 44.1, 70: 70.4}[f]):.3f} uK_RJ); "
              f"half-difference covariance over full-map covariance, median {np.median(np.r_[cq / qf, cu / uf]):.3f}")
        print(f"half difference: chi2 per element {chi.sum() / (2 * n):.4f} (expected 1, sd {np.sqrt(1.0 / n):.4f}); "
              f"rms {np.sqrt(np.mean(np.r_[dq, du] ** 2)):.3f} uK against white expectation "
              f"{np.sqrt(np.mean(np.r_[cq, cu])):.3f} uK")
        # power by multipole against white-noise simulations with the same covariance
        def spec(q, u):
            m = [np.zeros(n), hp.reorder(q, n2r=True), hp.reorder(u, n2r=True)]
            cl = hp.anafast(m, lmax=lmax, pol=True)
            return cl[1], cl[2]
        ee, bb = spec(dq, du)
        sims = []
        l11, l22 = np.sqrt(cq), np.sqrt(np.clip(cu - cx ** 2 / cq, 0, None))
        l21 = cx / np.sqrt(cq)
        for _ in range(200):
            z1, z2 = rng.standard_normal(n), rng.standard_normal(n)
            sims.append(spec(l11 * z1, l21 * z1 + l22 * z2))
        sims = np.array(sims)                       # (nsim, 2, lmax+1)
        me, sd = sims.mean(0), sims.std(0)
        for lo, hi in ((2, 3), (4, 6), (7, 10), (11, 20), (21, 29), (30, lmax)):
            sl = slice(lo, hi + 1)
            w = 2 * np.arange(lo, hi + 1) + 1
            re = np.sum(w * ee[sl]) / np.sum(w * me[0, sl])
            rb = np.sum(w * bb[sl]) / np.sum(w * me[1, sl])
            se = np.std([np.sum(w * s[0, sl]) / np.sum(w * me[0, sl]) for s in sims])
            sb = np.std([np.sum(w * s[1, sl]) / np.sum(w * me[1, sl]) for s in sims])
            print(f"  l {lo:2d} to {hi:2d}: EE ratio {re:6.3f} (white sd {se:.3f})   BB ratio {rb:6.3f} (white sd {sb:.3f})")
        if f == 30:
            dk, _ = D.load_wmap("K")
            if ns != 16:
                raise SystemExit("K-band comparison is at nside 16")
            BK, BL = wmap_band("K", D.WMAP_BANDPASS), lfi_band(30)
            pred = BL.power_law(-3.1, 23.0)[0][0] / BK.power_law(-3.1, 23.0)[0][0]
            for name, sl in (("Q", slice(0, n)), ("U", slice(n, 2 * n))):
                x, y = dk[sl], df[sl]
                big = np.abs(x) > 10.0
                slope = np.sum(x[big] * y[big]) / np.sum(x[big] ** 2)
                print(f"  K band against LFI 30, {name}: {big.sum()} elements with abs(K) > 10 uK, correlation "
                      f"{np.corrcoef(x[big], y[big])[0, 1]:.3f}, slope {slope:.3f} (predicted for index -3.1: {pred:.3f}; "
                      f"implied index {-3.1 + np.log(slope / pred) / np.log(BL.effective_frequency(-3.1) / BK.effective_frequency(-3.1)):.2f})")


if __name__ == "__main__":
    main()
