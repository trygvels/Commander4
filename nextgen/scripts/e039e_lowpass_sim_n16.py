"""E039 addendum B, spot check at the working resolution: the large-scale-only estimator at nside 16.

The E039 simulation at nside 16 (four BP XV regions, K, Ka, Q, V, 353; same synchrotron and dust truth), with
the CMB drawn from the LCDM pixel covariance and that covariance as the CMB prior (the central configuration
of the paper after E040 stage 4). K band and 353 GHz keep all their modes; Ka, Q and V enter through their
harmonic coefficients at multipoles 2 to L with their covariance (the construction of `e042e_lowpass.py`).
For each realization: the all-scales profile fit and the low-pass fit, with debiased sandwich errors.

    python -I nextgen/scripts/e039e_lowpass_sim_n16.py --lmax 20 --nreal 8
"""
import argparse
import gc
import os
import sys
import time

import healpy as hp
import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres.data import (DUST_BETA, DUST_REF, DUST_T, PR4_353, WMAP_BANDPASS, bpxv_regions, load_planck_qu, load_wmap)  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band, g_rj, wmap_band  # noqa: E402
from lowres.window import WindowMatcher, lcdm_cmb_precision  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=20)
    ap.add_argument("--nreal", type=int, default=8)
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--seed", type=int, default=39)
    a = ap.parse_args()
    ns = 16
    n = 2 * hp.nside2npix(ns)
    raw = []
    for b in ("K", "Ka", "Q", "V"):
        d, ninv = load_wmap(b)
        raw.append((f"WMAP_{b}", wmap_band(b, WMAP_BANDPASS), d, ninv))
    d3, n3 = load_planck_qu(PR4_353, ns)
    raw.append(("Planck_353", Band.delta("Planck_353", 353.0), d3, n3))
    region = bpxv_regions()
    beta_true = np.array([-3.15, -3.25, -2.85, -2.95])
    a30 = truth_sync_30(ns)
    a_s = a30 * (23.0 / 30.0) ** beta_true[region]
    a_d = d3 * g_rj(353.0)
    fd = [B.mbb(DUST_BETA, DUST_T, DUST_REF) for _, B, _, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]
    WM = WindowMatcher(ns, target_arcmin=60.0)
    prec, cov_c = lcdm_cmb_precision(WM)
    Lc = np.linalg.cholesky(cov_c)
    del cov_c
    keep = WM.ell <= a.lmax
    Yk, Ypk = WM.Y[:, keep], WM.Yp[keep]
    low_ninv = {}
    for (nm, B, _, ninv) in raw[1:4]:
        Cc = Ypk @ np.linalg.inv(ninv) @ Ypk.T
        nl = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
        low_ninv[nm] = 0.5 * (nl + nl.T)
        del Cc, nl
        gc.collect()
    s_true = [B.power_law(beta_true, 23.0)[0][region] for _, B, _, _ in raw]

    def simulate(r_):
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = Lc @ np.random.default_rng([a.seed, 11, r_]).standard_normal(n)
        return [s0 * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False)
                for s0, f, L in zip(s_true, fd, chol)]

    def bands_for(dvecs, lowpass):
        out = []
        for (nm, B, _, ninv), dv, f in zip(raw, dvecs, fd):
            if lowpass and nm in low_ninv:
                out.append(BandData(nm, Yk @ (Ypk @ dv), low_ninv[nm], (lambda be, B=B: B.power_law(be, 23.0)), f))
            else:
                out.append(BandData(nm, dv, ninv, (lambda be, B=B: B.power_law(be, 23.0)), f))
        return out

    print(f"# E039 addendum B at nside 16: K and 353 at all scales, Ka, Q, V through multipoles 2 to {a.lmax} ({int(keep.sum())} modes of "
          f"{WM.Y.shape[1]}); LCDM CMB drawn and LCDM prior; truth {beta_true}", flush=True)
    t0 = time.time()
    rows = []
    for r_ in range(a.first, a.first + a.nreal):
        dv = simulate(r_)
        res = []
        for lp in (False, True):
            L = RegionLikelihood(bands_for(dv, lp), region, cmb_prior_sd=None, cmb_prior_prec=prec)
            bh, _ = fit_profile(L, beta_true)
            _, cs, Jp, T, _, ok = sandwich_plugin(L, bh)
            res.append((bh, np.sqrt(np.diag(cs))))
            del L
            gc.collect()
        rows.append(res)
        print("real", r_, "all scales", " ".join(f"{x:+.4f}" for x in res[0][0]), "+-", " ".join(f"{x:.4f}" for x in res[0][1]),
              "| low-pass", " ".join(f"{x:+.4f}" for x in res[1][0]), "+-", " ".join(f"{x:.4f}" for x in res[1][1]), f"[{time.time() - t0:.0f} s]", flush=True)
    m = len(rows)
    for j, name in enumerate(("all scales", f"index from l <= {a.lmax}")):
        B_ = np.array([x[j][0] for x in rows])
        S = np.array([x[j][1] for x in rows])
        z = (B_ - beta_true) / S
        print(f"## {name}: {m} realizations")
        print("   mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (B_ - beta_true).mean(0)) + "  +- "
              + " ".join(f"{x:.4f}" for x in B_.std(0, ddof=1) / np.sqrt(m)))
        print(f"   pooled over regions: mean z {z.mean():+.3f} +- {1 / np.sqrt(z.size):.3f}, sd z {z.std(ddof=1):.3f}, coverage of +-1 sd {np.mean(np.abs(z) < 1):.3f}")
        print("   median quoted error   " + " ".join(f"{x:.4f}" for x in np.median(S, axis=0)) + "; scatter of the estimate "
              + " ".join(f"{x:.4f}" for x in B_.std(0, ddof=1)))
    dlt = np.array([x[1][0] - x[0][0] for x in rows])
    print("## low-pass minus all scales, paired: mean " + " ".join(f"{x:+.4f}" for x in dlt.mean(0)) + "  +- "
          + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1) / np.sqrt(m)) + "; sd " + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1)))


if __name__ == "__main__":
    main()
