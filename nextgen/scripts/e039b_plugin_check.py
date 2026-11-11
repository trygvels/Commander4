"""E039 addendum: the errors that one data set gives (plug-in information), on the E039 simulations.

Same truth, noise and seeds as `e039_lowres_validation.py`. For each realization: the profile optimum and
its errors from `fit.sandwich_plugin` (curvature J_c^-1 and sandwich J_c^-1 J(a_hat) J_c^-1 with
J_c = J(a_hat) - T), and the plain plug-in sandwich for comparison.

Original docstring: E039: simulations with the real WMAP noise covariances for the dense region-index estimator.

Truth: DR1 polarized synchrotron amplitudes (5 deg), PR4 353 GHz dust amplitudes, CMB N(0, 0.2^2), one
true index per region. Noise: drawn from the released WMAP covariances and the PR4 353 GHz Q/U covariance.
Fits: profile for every realization (with naive and sandwich errors at the truth), and for the first
`--full` realizations also the flat-prior marginal at references of 23 and 30 GHz and the zero-mean prior.

    python -I nextgen/scripts/e039b_plugin_check.py --nside 8 --regions hp2 --nreal 40
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
                         healpix_regions, load_planck_qu, load_wmap)
from lowres.fit import fit_profile, sandwich, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band, g_rj, wmap_band  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
CMB_SD = 0.2


def truth_sync_30(nside):
    """DR1 polarized synchrotron at 30 GHz (uK_RJ), 5 deg, NESTED (Q, U) at nside."""
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nside", type=int, default=8)
    ap.add_argument("--regions", choices=["bpxv", "hp2"], default="bpxv")
    ap.add_argument("--nreal", type=int, default=200)
    ap.add_argument("--full", type=int, default=0, help="realizations that also get the marginal fits")
    ap.add_argument("--seed", type=int, default=39)
    ap.add_argument("--first", type=int, default=0, help="index of the first realization (for parallel runs)")
    ap.add_argument("--no-drift", action="store_true",
                    help="WMAP bands without the nine-year frequency drift (the band model of the runs before 2026-11-10 22:04)")
    a = ap.parse_args()
    ns, n = a.nside, 2 * hp.nside2npix(a.nside)
    names = ["K", "Ka", "Q", "V"]
    raw = []
    for b in names:
        d, ninv = load_wmap(b)
        if ns != 16:
            d, ninv = degrade(d, ninv, 16, ns)
        raw.append((f"WMAP_{b}", wmap_band(b, WMAP_BANDPASS, drift=not a.no_drift), d, ninv))
    d3, n3 = load_planck_qu(PR4_353, ns)
    raw.append(("Planck_353", Band.delta("Planck_353", 353.0), d3, n3))
    if a.regions == "bpxv":
        reg16 = bpxv_regions()
        if ns == 16:
            region = reg16
        else:      # majority vote of the four nside-16 children, per Stokes block
            k = (16 // ns) ** 2
            r = reg16[:reg16.size // 2].reshape(-1, k)
            region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in r]), 2)
        beta_true = np.array([-3.15, -3.25, -2.85, -2.95])
    else:
        region = healpix_regions(2, ns)
        beta_true = -3.1 + 0.12 * np.random.default_rng(2).standard_normal(48)
    nreg = beta_true.size
    # truth
    a30 = truth_sync_30(ns)
    a_s = a30 * (23.0 / 30.0) ** beta_true[region]
    a_d = d3 * g_rj(353.0)                                   # 353 GHz data as the dust amplitude, uK_RJ
    rng_c = np.random.default_rng([a.seed, 1])
    fd = [B.mbb(DUST_BETA, DUST_T, DUST_REF) for _, B, _, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]   # N^-1 = L L^T; noise = L^-T xi

    def bands_for(dvecs, nu_ref):
        return [BandData(nm, dv, ninv, (lambda be, B=B: B.power_law(be, nu_ref)), f)
                for (nm, B, _, ninv), dv, f in zip(raw, dvecs, fd)]

    def simulate(r):
        rng = np.random.default_rng([a.seed, 7, r])
        a_c = CMB_SD * rng.standard_normal(n)
        out = []
        for (nm, B, _, _), f, L in zip(raw, fd, chol):
            s0, _ = B.power_law(beta_true, 23.0)
            out.append(s0[region] * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False))
        return out

    # expected information at the truth, from a noiseless data vector
    clean = [B.power_law(beta_true, 23.0)[0][region] * a_s + f * a_d for (nm, B, _, _), f in zip(raw, fd)]
    Lc = RegionLikelihood(bands_for(clean, 23.0), region, cmb_prior_sd=CMB_SD)
    Ji, Sw, J, T, sol = sandwich(Lc, beta_true)
    sd_naive, sd_sand = np.sqrt(np.diag(Ji)), np.sqrt(np.diag(Sw))
    nper = np.bincount(region, minlength=nreg)
    print(f"# E039 nside {ns}, regions {a.regions} ({nreg}), bands {[x[0] for x in raw]}, CMB prior sd {CMB_SD} uK")
    print(f"# truth amplitudes at 23 GHz: median |a_s| {np.median(np.abs(a_s)):.2f} uK_RJ; white-noise sd per element (K) "
          f"{np.median(1 / np.sqrt(np.diag(raw[0][3]))):.2f} uK")
    print("# region  n_el  beta_true  sd_naive  sd_sandwich  ratio")
    for r in range(nreg):
        print(f"# {r:3d} {nper[r]:6d} {beta_true[r]:+.3f} {sd_naive[r]:9.4f} {sd_sand[r]:9.4f} {sd_sand[r] / sd_naive[r]:6.3f}")
    sys.stdout.flush()
    t0 = time.time()
    Z = {k: [] for k in ("curv", "plain", "deb")}
    est, oks = [], []
    for r in range(a.first, a.first + a.nreal):
        dv = simulate(r)
        L23 = RegionLikelihood(bands_for(dv, 23.0), region, cmb_prior_sd=CMB_SD)
        bh, res = fit_profile(L23, beta_true)
        cc, cs, Jp, Tp, _, ok = sandwich_plugin(L23, bh)
        Jpi = np.linalg.inv(Jp)
        sd_plain = np.sqrt(np.diag(Jpi @ (Jp + Tp) @ Jpi))
        est.append(bh)
        oks.append(ok)
        Z["curv"].append((bh - beta_true) / np.sqrt(np.diag(cc)))
        Z["plain"].append((bh - beta_true) / sd_plain)
        Z["deb"].append((bh - beta_true) / np.sqrt(np.diag(cs)))
        print("real", r, "nok", int(ok.sum()), "z_deb", " ".join(f"{x:+.2f}" for x in Z["deb"][-1]), f"[{time.time() - t0:.0f} s]", flush=True)
    est, oks = np.array(est), np.array(oks)
    Z = {k: np.array(v) for k, v in Z.items()}
    meas = sd_sand < 0.25
    t = np.diag(T) / np.diag(J)
    print(f"\n## plug-in errors over {len(est)} realizations; T/J at the truth: min {t.min():.4f}, median {np.median(t):.4f}, max {t.max():.4f}; "
          f"regions with sandwich sd < 0.25 at the truth: {int(meas.sum())} of {nreg}; realizations with every region measured: {int(oks.all(1).sum())}")
    print("class of T/J at truth | regions | sd z: curvature, plain plug-in sandwich, debiased sandwich | coverage of +-1 sd: the same three | mean z debiased")
    for lo, hi in ((0, 0.02), (0.02, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 1e9)):
        sel = meas & (t >= lo) & (t < hi)
        if not sel.any():
            continue
        m = oks[:, sel]
        row = []
        for k in ("curv", "plain", "deb"):
            z = Z[k][:, sel][m]
            row.append((np.std(z, ddof=1), np.mean(np.abs(z) < 1)))
        zd = Z["deb"][:, sel][m]
        print(f"{lo:g} to {hi:g} | {int(sel.sum())} | {row[0][0]:.3f} {row[1][0]:.3f} {row[2][0]:.3f} | {row[0][1]:.3f} {row[1][1]:.3f} {row[2][1]:.3f} | "
              f"{zd.mean():+.3f} +- {zd.std(ddof=1) / np.sqrt(zd.size):.3f} (n = {zd.size})")
    z = Z["deb"][:, meas][oks[:, meas]]
    zc = Z["curv"][:, meas][oks[:, meas]]
    print(f"all measured regions: sd z curvature {zc.std(ddof=1):.3f}, debiased {z.std(ddof=1):.3f}; coverage curvature "
          f"{np.mean(np.abs(zc) < 1):.3f}, debiased {np.mean(np.abs(z) < 1):.3f}; mean z debiased {z.mean():+.3f}; n = {z.size}")


if __name__ == "__main__":
    main()
