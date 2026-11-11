"""E039: simulations with the real WMAP noise covariances for the dense region-index estimator.

Truth: DR1 polarized synchrotron amplitudes (5 deg), PR4 353 GHz dust amplitudes, CMB N(0, 0.2^2), one
true index per region. Noise: drawn from the released WMAP covariances and the PR4 353 GHz Q/U covariance.
Fits: profile for every realization (with naive and sandwich errors at the truth), and for the first
`--full` realizations also the flat-prior marginal at references of 23 and 30 GHz and the zero-mean prior.

    python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions bpxv --nreal 200 --full 40
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
from lowres.fit import fit_marginal, fit_profile, sandwich  # noqa: E402
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
    rows, t0 = [], time.time()
    for r in range(a.first, a.first + a.nreal):
        dv = simulate(r)
        L23 = RegionLikelihood(bands_for(dv, 23.0), region, cmb_prior_sd=CMB_SD)
        bh, res = fit_profile(L23, beta_true)
        dof = len(raw) * n - 2 * n      # CMB carries a prior of the size of its signal: about one constraint per element
        row = {"r": r, "prof": bh, "chi2": 2.0 * res.fun, "nit": res.nit, "dof": dof}
        if r < a.full:
            f23, _, _ = fit_marginal(L23, bh)
            L30 = RegionLikelihood(bands_for(dv, 30.0), region, cmb_prior_sd=CMB_SD)
            f30, _, _ = fit_marginal(L30, bh)
            Lz = RegionLikelihood(bands_for(dv, 23.0), region, cmb_prior_sd=CMB_SD, sync_prior=True)
            lt0 = np.log(np.array([np.sqrt(np.mean(a_s[region == k] ** 2)) for k in range(nreg)]))
            bz, tz, _ = fit_marginal(Lz, bh, lt0)
            row.update(flat23=f23, flat30=f30, zm=bz, zm_tau=np.exp(tz))
        rows.append(row)
        print("real", r, "prof", " ".join(f"{x:+.4f}" for x in bh), f"chi2/dof {row['chi2'] / dof:.4f} it {res.nit}"
              + ("" if r >= a.full else " | flat23 " + " ".join(f"{x:+.4f}" for x in row["flat23"]) + " | flat30 "
                 + " ".join(f"{x:+.4f}" for x in row["flat30"]) + " | zm " + " ".join(f"{x:+.4f}" for x in row["zm"])),
              f"[{time.time() - t0:.0f} s]", flush=True)
    P = np.array([x["prof"] for x in rows])
    zs, zn = (P - beta_true) / sd_sand, (P - beta_true) / sd_naive
    print("\n## profile over", len(rows), "realizations: per region mean shift, mean z (sandwich) +- se, sd z sandwich, sd z naive, "
          "coverage of +-1 sandwich sd")
    for r in range(nreg):
        print(f"{r:3d} {np.mean(P[:, r] - beta_true[r]):+.4f} {zs[:, r].mean():+.3f} {zs[:, r].std(ddof=1) / np.sqrt(len(rows)):.3f} "
              f"{zs[:, r].std(ddof=1):.3f} {zn[:, r].std(ddof=1):.3f} {np.mean(np.abs(zs[:, r]) < 1):.3f}")
    ok = sd_sand < 0.25
    print(f"pooled over the {int(ok.sum())} regions with sandwich sd < 0.25: mean z {zs[:, ok].mean():+.3f}, sd z sandwich "
          f"{zs[:, ok].std(ddof=1):.3f}, sd z naive {zn[:, ok].std(ddof=1):.3f}, coverage {np.mean(np.abs(zs[:, ok]) < 1):.3f}")
    c2 = np.array([x["chi2"] / x["dof"] for x in rows])
    print(f"chi2_min / (n_band - 2) n: mean {c2.mean():.4f}, sd {c2.std(ddof=1):.4f}")
    full = [x for x in rows if "flat23" in x]
    if full:
        F23, F30, Z = (np.array([x[k] for x in full]) for k in ("flat23", "flat30", "zm"))
        Pf = np.array([x["prof"] for x in full])
        print(f"\n## marginal targets over {len(full)} realizations: mean of (target - profile) per region, in index units and "
              "in sandwich sd; predictions for flat30 - flat23: scalar ln(30/23) n sd_naive^2, and coupled J^-1 (n ln(30/23))")
        predv = Ji @ (nper * np.log(30.0 / 23.0))
        for r in range(nreg):
            pred = np.log(30.0 / 23.0) * nper[r] * sd_naive[r] ** 2
            print(f"{r:3d} flat23 {np.mean(F23[:, r] - Pf[:, r]):+.4f} ({np.mean(F23[:, r] - Pf[:, r]) / sd_sand[r]:+.2f} sd)  "
                  f"flat30 {np.mean(F30[:, r] - Pf[:, r]):+.4f} ({np.mean(F30[:, r] - Pf[:, r]) / sd_sand[r]:+.2f} sd)  "
                  f"flat30 - flat23 {np.mean(F30[:, r] - F23[:, r]):+.4f} predicted scalar {pred:+.4f} coupled {predv[r]:+.4f}  "
                  f"zm {np.mean(Z[:, r] - Pf[:, r]):+.4f} ({np.mean(Z[:, r] - Pf[:, r]) / sd_sand[r]:+.2f} sd, rms "
                  f"{np.sqrt(np.mean((Z[:, r] - Pf[:, r]) ** 2)) / sd_sand[r]:.2f})")


if __name__ == "__main__":
    main()
