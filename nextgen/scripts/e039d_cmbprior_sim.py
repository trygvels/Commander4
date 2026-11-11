"""E040 stage 4, check C4: the three CMB priors of the central table on simulations whose CMB is LCDM.

The E039 simulations (four BP XV regions, K, Ka, Q, V, 353; same synchrotron and dust truth, same noise
seeds), but the CMB is drawn from the LCDM pixel covariance of `lowres.window.lcdm_cmb_precision` (E and B
spectra through the 60 arcmin window and the pixel window) and not from a white 0.2 uK field. Each
realization is fitted three times with the profile likelihood: with the LCDM covariance as the CMB prior
(matched), with the white 0.2 uK prior (right total variance, wrong shape), and with a flat CMB prior.

    python -I nextgen/scripts/e039d_cmbprior_sim.py --nside 8 --nreal 60
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
from lowres.window import WindowMatcher, lcdm_cmb_precision  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nside", type=int, default=8, choices=[8, 16])
    ap.add_argument("--nreal", type=int, default=60)
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--seed", type=int, default=39)
    ap.add_argument("--priors", nargs="+", default=["lcdm", "white", "flat"], choices=["lcdm", "white", "flat"])
    a = ap.parse_args()
    ns = a.nside
    n = 2 * hp.nside2npix(ns)
    raw = []
    for b in ("K", "Ka", "Q", "V"):
        d, ninv = load_wmap(b)
        if ns != 16:
            d, ninv = degrade(d, ninv, 16, ns)
        raw.append((f"WMAP_{b}", wmap_band(b, WMAP_BANDPASS), d, ninv))
    d3, n3 = load_planck_qu(PR4_353, ns)
    raw.append(("Planck_353", Band.delta("Planck_353", 353.0), d3, n3))
    reg16 = bpxv_regions()
    if ns == 16:
        region = reg16
    else:
        r = reg16[:reg16.size // 2].reshape(-1, (16 // ns) ** 2)
        region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in r]), 2)
    beta_true = np.array([-3.15, -3.25, -2.85, -2.95])
    a30 = truth_sync_30(ns)
    a_s = a30 * (23.0 / 30.0) ** beta_true[region]
    a_d = d3 * g_rj(353.0)
    fd = [B.mbb(DUST_BETA, DUST_T, DUST_REF) for _, B, _, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]
    WM = WindowMatcher(ns, target_arcmin=60.0)
    prec, cov_c = lcdm_cmb_precision(WM)
    Lc = np.linalg.cholesky(cov_c)
    s_true = [B.power_law(beta_true, 23.0)[0][region] for _, B, _, _ in raw]
    kws = {"lcdm": dict(cmb_prior_sd=None, cmb_prior_prec=prec), "white": dict(cmb_prior_sd=0.2), "flat": dict(cmb_prior_sd=None)}

    def simulate(r_):
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = Lc @ np.random.default_rng([a.seed, 11, r_]).standard_normal(n)        # LCDM CMB, uK
        return [s0 * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False)
                for s0, f, L in zip(s_true, fd, chol)]

    print(f"# E040 C4: nside {ns}, CMB drawn from the LCDM covariance (rms per element {np.sqrt(np.mean(np.diag(cov_c))):.3f} uK), "
          f"priors {a.priors}, truth {beta_true}", flush=True)
    t0 = time.time()
    est = {k: [] for k in a.priors}
    err = {k: [] for k in a.priors}
    for r_ in range(a.first, a.first + a.nreal):
        dv = simulate(r_)
        line = [f"real {r_}"]
        for k in a.priors:
            bands = [BandData(nm, d_, ninv, (lambda be, B=B: B.power_law(be, 23.0)), f) for (nm, B, _, ninv), d_, f in zip(raw, dv, fd)]
            L = RegionLikelihood(bands, region, **kws[k])
            bh, res = fit_profile(L, beta_true)
            _, cs, _, _, _, ok = sandwich_plugin(L, bh)
            est[k].append(bh)
            err[k].append(np.sqrt(np.diag(cs)))
            line.append(f"{k} " + " ".join(f"{x:+.4f}" for x in bh))
        print(" | ".join(line) + f" [{time.time() - t0:.0f} s]", flush=True)
    est = {k: np.array(v) for k, v in est.items()}
    err = {k: np.array(v) for k, v in err.items()}
    m = len(next(iter(est.values())))
    for k in a.priors:
        z = (est[k] - beta_true) / err[k]
        print(f"## CMB prior {k}: {m} realizations")
        print("   mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (est[k] - beta_true).mean(0))
              + "  +- " + " ".join(f"{x:.4f}" for x in est[k].std(0, ddof=1) / np.sqrt(m)))
        print("   mean z                " + " ".join(f"{x:+.3f}" for x in z.mean(0)) + f"  (se {1 / np.sqrt(m):.3f})")
        print("   sd of z               " + " ".join(f"{x:.3f}" for x in z.std(0, ddof=1)))
        print("   coverage of +-1 sd    " + " ".join(f"{x:.3f}" for x in np.mean(np.abs(z) < 1, axis=0)))
        print("   median quoted error   " + " ".join(f"{x:.4f}" for x in np.median(err[k], axis=0))
              + "; scatter of the estimate " + " ".join(f"{x:.4f}" for x in est[k].std(0, ddof=1)))
    if "lcdm" in est:
        for k in a.priors:
            if k == "lcdm":
                continue
            dlt = est[k] - est["lcdm"]
            print(f"## {k} minus lcdm, paired: mean " + " ".join(f"{x:+.4f}" for x in dlt.mean(0)) + "  +- "
                  + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1) / np.sqrt(m)) + "; sd " + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1)))


if __name__ == "__main__":
    main()
