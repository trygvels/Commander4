"""E040: polarized synchrotron index per sky region from WMAP K, Ka, Q, V and Planck 353 GHz at low
resolution with the released dense WMAP noise covariances; several estimators on the same maps.

    python -I nextgen/scripts/e040_wmap_region_indices.py --regions bpxv [--nside 16] [--bands K Ka Q V 353]
        [--cmb-prior-sd 0.2] [--delta-bands] [--dust-beta 1.53] [--estimators profile flat zm] [--out FILE.npz]

Prints, per region: profile optimum with curvature and sandwich errors (both from the debiased
information J(a_hat) - T); the flat-prior marginal optimum at
references of 23 and 30 GHz with the coupled prediction J^-1 (n ln(30/23)) for their difference; the
zero-mean-prior optimum. Also chi-square at the profile optimum, in total and per band.
"""
import argparse
import os
import sys
import time

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.fit import fit_marginal, fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import WMAP_NOMINAL, Band, hfi_band, lfi_band, wmap_band  # noqa: E402

LFI_NOMINAL = {"30": 28.4, "44": 44.1, "70": 70.4}
REGION_NAMES = {"bpxv": ["High Latitude", "Spur", "Galactic Center", "Plane"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nside", type=int, default=16)
    ap.add_argument("--regions", default="bpxv", choices=["bpxv", "hp2", "hp1", "lat", "dr1"])
    ap.add_argument("--bands", nargs="+", default=["K", "Ka", "Q", "V", "353"])
    ap.add_argument("--cmb-prior-sd", type=float, default=0.2, help="uK thermodynamic; 0 means flat")
    ap.add_argument("--delta-bands", action="store_true", help="delta bandpasses at the nominal WMAP frequencies")
    ap.add_argument("--dust-beta", type=float, default=D.DUST_BETA)
    ap.add_argument("--bp353", default="psb", choices=["psb", "full", "delta"],
                    help="353 GHz bandpass: polarization-sensitive bolometers, all bolometers, or a delta at 353.0")
    ap.add_argument("--lfi-split", default="full", help="PR4 LFI map: full, full-ringhalf-1 or full-ringhalf-2")
    ap.add_argument("--lfi-noise-scale", type=float, default=1.0, help="factor on the LFI noise sd")
    ap.add_argument("--lfi-noise", default="model", choices=["model", "white"],
                    help="scale-dependent model fitted to the ring-half difference (E038 b2), or the map-file columns")
    ap.add_argument("--lfi-a-scale", type=float, default=1.0, help="factor on the power-law term of the LFI noise model")
    ap.add_argument("--estimators", nargs="+", default=["profile"], choices=["profile", "flat", "zm"])
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    ns, npix = a.nside, hp.nside2npix(a.nside)
    t0 = time.time()
    raw = []
    for b in a.bands:
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, ns)
            B = (Band.delta("Planck_353", 353.0) if (a.delta_bands or a.bp353 == "delta")
                 else hfi_band(353, psb=(a.bp353 == "psb")))
        elif b in LFI_NOMINAL:
            d, ninv = D.load_lfi(int(b), ns, noise=a.lfi_noise, split=a.lfi_split, a_scale=a.lfi_a_scale)
            if a.lfi_noise_scale != 1.0:
                ninv = ninv / a.lfi_noise_scale ** 2
            B = Band.delta(f"Planck_{int(b):03d}", LFI_NOMINAL[b]) if a.delta_bands else lfi_band(int(b))
        else:
            d, ninv = D.load_wmap(b)
            if ns != 16:
                d, ninv = D.degrade(d, ninv, 16, ns)
            B = Band.delta(f"WMAP_{b}", WMAP_NOMINAL[b]) if a.delta_bands else wmap_band(b, D.WMAP_BANDPASS)
        raw.append((B, d, ninv, B.mbb(a.dust_beta, D.DUST_T, D.DUST_REF)))

    def bands_for(nu_ref):
        return [BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, nu_ref)), f) for B, d, ninv, f in raw]

    if a.regions == "bpxv":
        reg16 = D.bpxv_regions()
        if ns == 16:
            region = reg16
        else:
            k = (16 // ns) ** 2
            r = reg16[:reg16.size // 2].reshape(-1, k)
            region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in r]), 2)
    elif a.regions == "hp2":
        region = D.healpix_regions(2, ns)
    elif a.regions == "hp1":
        region = D.healpix_regions(1, ns)
    elif a.regions == "lat":
        region = D.latitude_regions((0, 20, 45, 90), ns)
    else:
        region, _ = D.label_regions("/Users/trygvels/data/cosmoglobe_DR1/synch_beta/region_labels.fits", 1, ns)
    nreg = int(region.max()) + 1
    names = REGION_NAMES.get(a.regions, [str(i) for i in range(nreg)])
    n = region.size
    prior = a.cmb_prior_sd if a.cmb_prior_sd > 0 else None
    has_dust = "353" in a.bands
    kw = dict(cmb_prior_sd=prior, use_dust=has_dust)
    L = RegionLikelihood(bands_for(23.0), region, **kw)
    bh, res = fit_profile(L, np.full(nreg, -3.1))
    # errors from this one data set: J(a_hat) - T estimates J, and J(a_hat) the score covariance (fit.sandwich_plugin)
    Ji, Sw, J, T, sol, okreg = sandwich_plugin(L, bh)
    sdn, sds = np.sqrt(np.diag(Ji)), np.sqrt(np.diag(Sw))
    nper = np.bincount(region, minlength=nreg)
    chi2 = 2.0 * res.fun
    ncomp = L.ncomp
    dof = len(raw) * n - ncomp * n - nreg + (n if prior else 0)
    print(f"# E040 nside {ns}, regions {a.regions} ({nreg}), bands {a.bands}, CMB prior sd {a.cmb_prior_sd}, "
          f"{'delta bands at nominal frequencies' if a.delta_bands else 'band-averaged spectra'}, dust beta {a.dust_beta}, "
          f"353 bandpass {a.bp353}, LFI split {a.lfi_split}, LFI noise {a.lfi_noise} x {a.lfi_noise_scale} (A x {a.lfi_a_scale})")
    print("# effective frequencies for index -3.1: " + "  ".join(
        f"{B.name} {B.effective_frequency(-3.1):.2f}" for B, _, _, _ in raw if B.nu.size > 1 and B.nu.max() < 200))
    print(f"# profile fit: {res.nit} iterations, {res.nfev} evaluations, converged {res.success}; chi2_min {chi2:.1f}; "
          f"nominal dof {dof} (a CMB prior counted as one constraint per element); ratio {chi2 / dof:.4f}; {time.time() - t0:.0f} s")
    rb = L.residuals(sol)
    print("# chi2 per band r^T N^-1 r / n: " + "  ".join(f"{b.name} {float(r @ b.apply(r)) / n:.3f}" for b, r in zip(L.bands, rb)))
    # chi-square by region: sum over bands of r_R^T (N^-1 r)_R, which adds up to the total without the prior term
    creg = np.zeros(nreg)
    for b, r in zip(L.bands, rb):
        creg += np.bincount(region, weights=r * b.apply(r), minlength=nreg)
    print("# chi2 by region, sum_b r_R^T (N^-1 r)_R / (n_R n_bands): " + "  ".join(
        f"{names[k]} {creg[k] / (nper[k] * len(raw)):.3f}" for k in range(nreg)))
    a_s = sol["x"][:n]
    out = {"beta": bh, "sd_naive": sdn, "sd_sand": sds, "J": J, "T": T, "region": region, "a_s": a_s, "chi2": chi2, "dof": dof,
           "resid": np.array(rb), "ok": okreg, "cov_sand": Sw, "nper": nper}
    full = {}
    if "flat" in a.estimators:
        f23, _, r23 = fit_marginal(L, bh)
        L30 = RegionLikelihood(bands_for(30.0), region, **kw)
        f30, _, r30 = fit_marginal(L30, bh)
        full["flat23"], full["flat30"] = f23, f30
        out.update(flat23=f23, flat30=f30)
        print(f"# flat-prior marginal fits: ref 23 {r23.nit} it, ref 30 {r30.nit} it; {time.time() - t0:.0f} s")
    if "zm" in a.estimators:
        Lz = RegionLikelihood(bands_for(23.0), region, sync_prior=True, **kw)
        lt0 = np.log(np.array([np.sqrt(np.mean(a_s[region == k] ** 2)) for k in range(nreg)]))
        bz, tz, rz = fit_marginal(Lz, bh, lt0)
        full["zm"] = bz
        out.update(zm=bz, zm_tau=np.exp(tz))
        print(f"# zero-mean prior fit: {rz.nit} it; {time.time() - t0:.0f} s")
    Jfin = np.where(np.isfinite(Ji), Ji, 0.0)
    predv = Jfin @ (nper * okreg * np.log(30.0 / 23.0))
    print("# T / J_c per region (noise share of the score variance): " + "  ".join(
        f"{names[k]} {T[k, k] / max(J[k, k] - T[k, k], 1e-30):.4f}" for k in range(nreg)))
    hdr = f"{'region':16s} {'n_el':>6s} {'rms a_s':>8s} {'profile':>9s} {'sd naive':>9s} {'sd sand':>9s}"
    if "flat23" in full:
        hdr += f" {'flat23-p':>9s} {'flat30-p':>9s} {'30-23':>8s} {'pred':>8s}"
    if "zm" in full:
        hdr += f" {'zm-p':>8s} {'tau':>7s}"
    print(hdr)
    for r in range(nreg):
        line = (f"{names[r]:16s} {nper[r]:6d} {np.sqrt(np.mean(a_s[region == r] ** 2)):8.2f} {bh[r]:+9.4f} {sdn[r]:9.4f} "
                f"{sds[r]:9.4f}")
        if "flat23" in full:
            line += (f" {full['flat23'][r] - bh[r]:+9.4f} {full['flat30'][r] - bh[r]:+9.4f} "
                     f"{full['flat30'][r] - full['flat23'][r]:+8.4f} {predv[r]:+8.4f}")
        if "zm" in full:
            line += f" {full['zm'][r] - bh[r]:+8.4f} {out['zm_tau'][r]:7.2f}"
        print(line)
    if a.out:
        np.savez(a.out, **out)


if __name__ == "__main__":
    main()
