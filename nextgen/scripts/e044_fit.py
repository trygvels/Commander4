"""E044: the estimator comparison of E040 on regions that are what they are named, with a processing mask.

WMAP K, Ka, Q, V and Planck 353 GHz polarization at nside 16 with the released dense covariances, all bands brought
to a common window, masked elements removed exactly (`lowres.mask`), regions and mask from `lowres.regions`.
One fit per run, one likelihood object alive at a time.

    python -I nextgen/scripts/e044_fit.py [--prior lcdm|white|flat] [--estimators profile flat zm] [--block]
        [--mismatch] [--lowpass L] [--bands K Ka Q V 353] [--dust-beta 1.53] [--no-drift] [--delta-bands]
        [--bp353 psb|full] [--geometry named|bpxv] [--no-mask] [--out FILE.npz]

Prints: the profile optimum with its debiased sandwich error; chi-square in total, per band and per region; with
`flat`, the optima of the flat-prior marginal with the synchrotron amplitude defined at 23 and at 30 GHz and two
first-order predictions of their difference (with the curvature of the fitted marginal, and with J_c, which needs
no flat-prior fit); with `zm`, the zero-mean prior with a free scale; with `--block`, block errors from the joint
score scatter and from each region's own scores (`lowres.block`); with `--mismatch` (white prior), the exact
first-order bias and error of this fit if the CMB has the LCDM covariance (`lowres.mismatch`).
`--lowpass L`: the first band and 353 GHz keep all modes, the others enter through multipoles 2 to L.
`--nside 8` exists for an end-to-end test of this script only.
"""
import argparse
import gc
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.block import KINDS, block_errors  # noqa: E402
from lowres.fit import fit_marginal, fit_profile, marginal_hessian, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.mask import marginalize, restrict_prior  # noqa: E402
from lowres.regions import NAMES, named_regions  # noqa: E402
from lowres.spectra import WMAP_NOMINAL, Band, hfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher, lcdm_cmb_precision  # noqa: E402

BPXV_NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["K", "Ka", "Q", "V", "353"])
    ap.add_argument("--prior", default="lcdm", choices=["lcdm", "white", "flat"])
    ap.add_argument("--white-sd", type=float, default=0.2)
    ap.add_argument("--floor-uk", type=float, default=0.01)
    ap.add_argument("--estimators", nargs="+", default=["profile"], choices=["profile", "flat", "zm"])
    ap.add_argument("--block", action="store_true")
    ap.add_argument("--mismatch", action="store_true", help="white-prior fit judged under the LCDM covariance")
    ap.add_argument("--lowpass", type=int, default=0)
    ap.add_argument("--dust-beta", type=float, default=D.DUST_BETA)
    ap.add_argument("--no-drift", action="store_true")
    ap.add_argument("--delta-bands", action="store_true")
    ap.add_argument("--bp353", default="psb", choices=["psb", "full"])
    ap.add_argument("--geometry", default="named", choices=["named", "bpxv"])
    ap.add_argument("--no-mask", action="store_true")
    ap.add_argument("--tau-a-only", action="store_true", help="mask the disc around Tau A but not the brightest 5 percent")
    ap.add_argument("--bright-fraction", type=float, default=0.05, help="fraction of pixels masked by polarized amplitude (referee read 3: 0.03 and 0.10 as arms)")
    ap.add_argument("--mask-extra", default="", help="npy file, one boolean per nside-16 pixel in NESTED order, True for pixels to mask in addition (H16)")
    ap.add_argument("--beta0", type=float, default=-3.1)
    ap.add_argument("--nside", type=int, default=16, choices=[8, 16])
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    ns = a.nside
    if a.geometry == "named":
        region, keep, block = named_regions(ns, mask=not a.no_mask, bright_fraction=(0.0 if a.tau_a_only else a.bright_fraction))
        names = NAMES
    else:
        assert ns == 16
        region = D.bpxv_regions()
        keep = np.ones(region.size, dtype=bool)
        if not a.no_mask:
            _, keep, _ = named_regions(16, bright_fraction=(0.0 if a.tau_a_only else a.bright_fraction))
        block = np.tile(np.arange(region.size // 2) // 64, 2)
        names = BPXV_NAMES
    if a.mask_extra:
        assert ns == 16
        extra = np.load(a.mask_extra)
        assert extra.dtype == bool and extra.size == keep.size // 2
        keep = keep & ~np.tile(extra, 2)
    WM16 = WindowMatcher(16, target_arcmin=60.0)
    WMn = WM16 if ns == 16 else WindowMatcher(ns, target_arcmin=60.0)
    lp = a.lowpass > 0
    if lp:
        sel = WMn.ell <= a.lowpass
        Yk, Ypk = WMn.Y[:, sel], WMn.Yp[sel]
        assert np.abs(Ypk @ Yk - np.eye(Yk.shape[1])).max() < 1e-8
    raw = []
    for ib, b in enumerate(a.bands):
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=(a.bp353 == "psb"))
        else:
            d, ninv = D.load_wmap(b)
            B = Band.delta(f"WMAP_{b}", WMAP_NOMINAL[b]) if a.delta_bands else wmap_band(b, D.WMAP_BANDPASS, drift=not a.no_drift)
        d, ninv = WM16.apply(d, ninv, FWHM_ARCMIN[b])
        if ns != 16:
            d, ninv = D.degrade(d, ninv, 16, ns)
        low = lp and ib > 0 and b != "353"
        if low:
            Cc = Ypk @ np.linalg.inv(ninv) @ Ypk.T
            nl = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
            d, ninv = Yk @ (Ypk @ d), 0.5 * (nl + nl.T)
            del Cc, nl
        d, ninv = marginalize(d, ninv, keep, pinv=low)
        raw.append((B, d, ninv, B.mbb(a.dust_beta, D.DUST_T, D.DUST_REF)))
        gc.collect()
    region, block = region[keep], block[keep]
    n, nreg = region.size, int(region.max()) + 1
    nper = np.bincount(region, minlength=nreg)
    cov_l = None
    if a.prior == "lcdm" or a.mismatch:
        _, cov_full = lcdm_cmb_precision(WMn, floor_uk=a.floor_uk)
        cov_l = cov_full[np.ix_(keep, keep)]
        del cov_full
    if a.prior == "lcdm":
        kw = dict(cmb_prior_sd=None, cmb_prior_prec=restrict_prior(cov_l, np.ones(n, dtype=bool)))
    elif a.prior == "white":
        kw = dict(cmb_prior_sd=a.white_sd)
    else:
        kw = dict(cmb_prior_sd=None)
    if not a.mismatch:
        cov_l = None
    del WM16
    if not lp:
        del WMn
    gc.collect()

    def bands_for(nu_ref):
        return [BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, nu_ref)), f) for B, d, ninv, f in raw]

    L = RegionLikelihood(bands_for(23.0), region, **kw)
    bh, res = fit_profile(L, np.full(nreg, a.beta0))
    Ji, Sw, J, T, sol, okreg = sandwich_plugin(L, bh)
    sdn, sds = np.sqrt(np.diag(Ji)), np.sqrt(np.diag(Sw))
    chi2 = 2.0 * res.fun
    dof = len(raw) * n - L.ncomp * n - nreg + (n if a.prior != "flat" else 0)
    print(f"# E044 nside {ns}, geometry {a.geometry}, mask {'none' if a.no_mask else ('Tau A disc' if a.tau_a_only else f'brightest {100 * a.bright_fraction:g} percent and Tau A disc')} "
          f"({int(keep.sum())} of {keep.size} elements kept), bands {a.bands}, CMB prior {a.prior}"
          + (f" (floor {a.floor_uk} uK)" if a.prior == "lcdm" else (f" ({a.white_sd} uK)" if a.prior == "white" else ""))
          + (f", extra mask {os.path.basename(a.mask_extra)}" if a.mask_extra else "")
          + f", dust beta {a.dust_beta}, 353 bandpass {a.bp353}, {'delta bands' if a.delta_bands else 'band-averaged spectra'}, "
          f"drift correction {not a.no_drift}, low-pass {a.lowpass if lp else 'none'}")
    print(f"# profile fit: {res.nit} iterations, {res.nfev} evaluations, converged {res.success}; chi2_min {chi2:.1f}; nominal dof {dof} "
          f"(a CMB prior counted as one constraint per element); ratio {chi2 / dof:.4f}; {time.time() - t0:.0f} s")
    rb = L.residuals(sol)
    if not lp:
        print("# chi2 per band r^T N^-1 r / n: " + "  ".join(f"{b.name} {float(r @ b.apply(r)) / n:.3f}" for b, r in zip(L.bands, rb)))
        creg = np.zeros(nreg)
        for b, r in zip(L.bands, rb):
            creg += np.bincount(region, weights=r * b.apply(r), minlength=nreg)
        creg = creg / (nper * len(raw))
        print("# chi2 by region, sum_b r_R^T (N^-1 r)_R / (n_R n_bands): " + "  ".join(f"{names[k]} {creg[k]:.3f}" for k in range(nreg)))
    else:
        creg = np.full(nreg, np.nan)
    a_s = sol["x"][:n]
    out = {"beta": bh, "sd_naive": sdn, "sd_sand": sds, "J": J, "T": T, "region": region, "keep": keep, "block": block, "a_s": a_s,
           "chi2": chi2, "dof": dof, "chi2_region": creg, "resid": np.array(rb), "ok": okreg, "cov_sand": Sw, "nper": nper,
           "converged": np.array(bool(res.success))}
    sol = None
    gc.collect()
    Jc = J - T
    pred_j = np.linalg.solve(Jc, nper * np.log(30.0 / 23.0)) if okreg.all() else np.full(nreg, np.nan)
    out["pred_jc"] = pred_j
    blk = None
    if a.block:
        blk = block_errors(L, bh, block, 48)
        for k_ in list(KINDS) + [k + "_own" for k in KINDS] + ["neff", "lev_max", "n_neg", "cross", "d_over_jc", "G", "D"]:
            out["blk_" + k_] = blk[k_]
        print(f"# block errors: blocks with negative information per region {blk['n_neg'].tolist()}; cross-region share of the block Jacobian "
              f"{blk['cross']:.3f}; sum D / J_c " + " ".join(f"{x:.3f}" for x in blk["d_over_jc"]) + f"; {time.time() - t0:.0f} s")
    if a.mismatch:
        assert a.prior == "white"
        from lowres.mismatch import cmb_prior_mismatch
        mm = cmb_prior_mismatch(L, bh, 1.0 / a.white_sd ** 2, cov_l - a.white_sd ** 2 * np.eye(n))
        out.update(mm_bias=mm["bias"], mm_cov_true=mm["cov_true"], mm_dK=mm["dK"], mm_eg=mm["eg"])
        sdt = np.sqrt(np.diag(mm["cov_true"]))
        print("# white prior judged under the LCDM covariance (first order): bias of the index " + " ".join(f"{x:+.4f}" for x in mm["bias"])
              + "; error " + " ".join(f"{x:.4f}" for x in sdt) + " where the fit quotes " + " ".join(f"{x:.4f}" for x in sds)
              + " (ratio " + " ".join(f"{x / y:.3f}" for x, y in zip(sdt, sds)) + f"); {time.time() - t0:.0f} s")
        del mm, cov_l
        gc.collect()
    full = {}
    if "flat" in a.estimators:
        f23, _, r23 = fit_marginal(L, bh)
        Hm = marginal_hessian(L, f23)
        del L
        gc.collect()
        L30 = RegionLikelihood(bands_for(30.0), region, **kw)
        f30, _, r30 = fit_marginal(L30, bh)
        del L30
        gc.collect()
        full["flat23"], full["flat30"] = f23, f30
        out.update(flat23=f23, flat30=f30, Hmarg=Hm, pred_marg=np.linalg.solve(Hm, nper * np.log(30.0 / 23.0)),
                   flat_converged=np.array([bool(r23.success), bool(r30.success)]))
        print(f"# flat-prior marginal fits: ref 23 {r23.nit} it (converged {r23.success}), ref 30 {r30.nit} it (converged {r30.success}); {time.time() - t0:.0f} s")
        L = None
    if "zm" in a.estimators:
        Lz = RegionLikelihood(bands_for(23.0), region, sync_prior=True, **kw)
        lt0 = np.log(np.array([np.sqrt(np.mean(a_s[region == k] ** 2)) for k in range(nreg)]))
        bz, tz, rz = fit_marginal(Lz, bh, lt0)
        del Lz
        gc.collect()
        full["zm"] = bz
        out.update(zm=bz, zm_tau=np.exp(tz), zm_converged=np.array(bool(rz.success)))
        print(f"# zero-mean prior fit: {rz.nit} it (converged {rz.success}); {time.time() - t0:.0f} s")
    if np.any(np.abs(bh + 4.6) < 1e-6) or np.any(np.abs(bh + 1.6) < 1e-6):
        print("# WARNING: a region index sits on a bound of the allowed range")
    print("# T / J_c per region: " + "  ".join(f"{names[k]} {T[k, k] / max(Jc[k, k], 1e-30):.4f}" for k in range(nreg)))
    hdr = f"{'region':16s} {'n_el':>6s} {'rms a_s':>8s} {'profile':>9s} {'sd sand':>9s}"
    if blk is not None:
        hdr += f" {'CR3':>8s} {'CR3 own':>8s} {'n_eff':>6s} {'max h':>6s}"
    if "flat23" in full:
        hdr += f" {'flat23-p':>9s} {'flat30-p':>9s} {'30-23':>8s} {'pred J_c':>8s} {'pred H':>8s}"
    if "zm" in full:
        hdr += f" {'zm-p':>8s} {'tau':>7s}"
    print(hdr)
    for r in range(nreg):
        line = f"{names[r]:16s} {nper[r]:6d} {np.sqrt(np.mean(a_s[region == r] ** 2)):8.2f} {bh[r]:+9.4f} {sds[r]:9.4f}"
        if blk is not None:
            line += f" {blk['CR3'][r]:8.4f} {blk['CR3_own'][r]:8.4f} {blk['neff'][r]:6.1f} {blk['lev_max'][r]:6.2f}"
        if "flat23" in full:
            line += (f" {full['flat23'][r] - bh[r]:+9.4f} {full['flat30'][r] - bh[r]:+9.4f} {full['flat30'][r] - full['flat23'][r]:+8.4f} "
                     f"{pred_j[r]:+8.4f} {out['pred_marg'][r]:+8.4f}")
        if "zm" in full:
            line += f" {full['zm'][r] - bh[r]:+8.4f} {out['zm_tau'][r]:7.2f}"
        print(line)
    if a.out:
        np.savez(a.out, **out)


if __name__ == "__main__":
    main()
