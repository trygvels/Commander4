"""E039 addendum: the estimator that reads the index from large scales only, on simulations with known truth.

Same truth, noise and seeds as `e039_lowres_validation.py` at nside 8 (four BP XV regions, K, Ka, Q, V, 353).
K band and 353 GHz keep all their modes; Ka, Q and V enter through their harmonic coefficients at
multipoles 2 to L with their covariance (the construction of `e042e_lowpass.py`). For each realization: the
profile optimum and its debiased sandwich error.

    python -I nextgen/scripts/e039c_lowpass_sim.py --lmax 12 --nreal 100
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
from lowres.lfi_noise import spin2_basis  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band, g_rj, wmap_band  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
CMB_SD = 0.2


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmax", type=int, default=12)
    ap.add_argument("--nreal", type=int, default=100)
    ap.add_argument("--seed", type=int, default=39)
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
    beta_true = np.array([-3.15, -3.25, -2.85, -2.95])
    a30 = truth_sync_30(ns)
    a_s = a30 * (23.0 / 30.0) ** beta_true[region]
    a_d = d3 * g_rj(353.0)
    fd = [B.mbb(DUST_BETA, DUST_T, DUST_REF) for _, B, _, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]
    Y, ell = spin2_basis(ns, 3 * ns - 1)
    Yp = np.linalg.solve(Y.T @ Y, Y.T)
    keep = ell <= a.lmax
    Yk, Ypk = Y[:, keep], Yp[keep]
    low_ninv = {}
    for (nm, B, _, ninv) in raw[1:4]:
        Cc = Ypk @ np.linalg.inv(ninv) @ Ypk.T
        nl = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
        low_ninv[nm] = 0.5 * (nl + nl.T)

    def simulate(r_):
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = CMB_SD * rng.standard_normal(n)
        out = []
        for (nm, B, _, _), f, L in zip(raw, fd, chol):
            s0, _ = B.power_law(beta_true, 23.0)
            out.append(s0[region] * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False))
        return out

    def bands_for(dvecs, lowpass):
        out = []
        for (nm, B, _, ninv), dv, f in zip(raw, dvecs, fd):
            if lowpass and nm in low_ninv:
                out.append(BandData(nm, Yk @ (Ypk @ dv), low_ninv[nm], (lambda be, B=B: B.power_law(be, 23.0)), f))
            else:
                out.append(BandData(nm, dv, ninv, (lambda be, B=B: B.power_law(be, 23.0)), f))
        return out

    print(f"# E039 low-pass addendum: nside 8, four BP XV regions, K and 353 at all scales, Ka, Q, V through multipoles 2 to {a.lmax} "
          f"({int(keep.sum())} modes of {Y.shape[1]}); truth {beta_true}", flush=True)
    t0 = time.time()
    rows = []
    for r_ in range(a.nreal):
        dv = simulate(r_)
        res = []
        for lp in (False, True):
            L = RegionLikelihood(bands_for(dv, lp), region, cmb_prior_sd=CMB_SD)
            bh, _ = fit_profile(L, beta_true)
            _, cs, _, _, _, ok = sandwich_plugin(L, bh)
            res.append((bh, np.sqrt(np.diag(cs))))
        rows.append(res)
        if r_ % 10 == 0:
            print("real", r_, "all scales", " ".join(f"{x:+.4f}" for x in res[0][0]), "| low-pass", " ".join(f"{x:+.4f}" for x in res[1][0]),
                  f"[{time.time() - t0:.0f} s]", flush=True)
    for j, name in enumerate(("all scales", f"index from l <= {a.lmax}")):
        B = np.array([x[j][0] for x in rows])
        S = np.array([x[j][1] for x in rows])
        z = (B - beta_true) / S
        print(f"## {name}: over {len(rows)} realizations")
        print("   mean shift          " + " ".join(f"{x:+.4f}" for x in (B - beta_true).mean(0)))
        print("   mean z +- se        " + " ".join(f"{x:+.3f}" for x in z.mean(0)) + f"  (se {1 / np.sqrt(len(rows)):.3f})")
        print("   sd of z             " + " ".join(f"{x:.3f}" for x in z.std(0, ddof=1)))
        print("   coverage of +-1 sd  " + " ".join(f"{x:.3f}" for x in np.mean(np.abs(z) < 1, axis=0)))
        print("   median quoted error " + " ".join(f"{x:.4f}" for x in np.median(S, axis=0)) + "; scatter of the estimate "
              + " ".join(f"{x:.4f}" for x in B.std(0, ddof=1)))
    D_ = np.array([x[1][0] - x[0][0] for x in rows])
    print("## low-pass minus all scales: mean " + " ".join(f"{x:+.4f}" for x in D_.mean(0)) + "; sd " + " ".join(f"{x:.4f}" for x in D_.std(0, ddof=1)))


if __name__ == "__main__":
    main()
