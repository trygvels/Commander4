"""E036: amplitude priors for signed Q/U amplitudes at real WMAP/LFI depth; coverage of a region index.

Truth: Cosmoglobe DR1 polarized synchrotron posterior mean (Q, U at 30 GHz, uK_RJ) smoothed to 5 deg
at nside 16, one index over the sky. Noise: real Q/U white-noise rms per nside-16 pixel (DR1 rms
maps, through the E031 helper). Estimators of one index per region, exact on a grid:

  flat    flat-prior marginal            -(S - b^2/A)/2 - (1/2) log A      (reference dependent)
  prof    profile likelihood             -(S - b^2/A)/2
  zm      zero-mean prior, a ~ N(0, tau_r^2), tau_r flat in log, marginalized on a grid
  zm_tpl  the same with tau_p = s_r g_p, g_p a polarized-intensity template with a 30 percent error
  zm_orc  the same with the oracle scale map g_p = |a_p|

with S = sum_b w d^2, b = sum_b w f d, A = sum_b w f^2, f_b = (nu_b/nu_ref)^beta, w = 1/sigma^2,
all per pixel and Stokes parameter. See the E036 record for the predictions.

    python -I nextgen/scripts/e036_qu_amplitude_priors.py [--nreal 100] [--beta -3.1]
"""
import argparse
import os
import sys

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e031_noise_regime_fisher import cmb_sed, dust_sed, real_noise  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
DR1_DUST = "/Users/trygvels/data/cosmoglobe_DR1/CG_dust_IQU_n1024_v1.fits"   # Q, U referenced at 70 GHz (F13)
CENTRE = {23.0: 22.8, 30.0: 28.4, 33.0: 33.0, 41.0: 40.6, 44.0: 44.1, 61.0: 60.8, 70.0: 70.1,
          353.0: 353.0}
DUST_BETA, DUST_T, DUST_REF = 1.55, 19.6, 70.0


def summarize(logL, betas, truth):
    """Posterior mean, sd and CDF at the truth for each column of logL (nbeta, nreg), flat prior."""
    p = np.exp(logL - logL.max(axis=0, keepdims=True))
    p /= p.sum(axis=0, keepdims=True)
    mean = (betas[:, None] * p).sum(axis=0)
    sd = np.sqrt(((betas[:, None] - mean[None, :]) ** 2 * p).sum(axis=0))
    i = int(round((truth - betas[0]) / (betas[1] - betas[0])))
    cdf = p[:i].sum(axis=0) + 0.5 * p[i]
    edge = p[:3].sum(axis=0) + p[-3:].sum(axis=0)
    return mean, sd, cdf, edge


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nreal", type=int, default=100)
    ap.add_argument("--beta", type=float, default=-3.1)
    ap.add_argument("--bands", type=float, nargs="+", default=[23.0, 30.0, 33.0, 41.0, 44.0])
    ap.add_argument("--nu-ref", type=float, default=23.0)
    ap.add_argument("--seed", type=int, default=36)
    ap.add_argument("--ntau", type=int, default=32)
    ap.add_argument("--noise-scale", type=float, default=1.0)
    ap.add_argument("--nuisance", choices=["none", "cmb_dust"], default="none",
                    help="cmb_dust: a free CMB and a free dust amplitude per pixel and Stokes parameter, SEDs "
                         "known, flat priors (the CMB optionally Gaussian), projected out exactly")
    ap.add_argument("--cmb-prior-sd", type=float, default=0.0, help="uK_CMB; 0 means a flat CMB prior")
    ap.add_argument("--global-arm", action="store_true",
                    help="add zero-mean priors with ONE scale shared by all regions times a scale map (the true "
                         "local rms, and the same 30 percent high at |b| > 30 deg), stated at 23 GHz, the pivot and 30 GHz")
    a = ap.parse_args()
    if a.nuisance == "cmb_dust" and a.bands == [23.0, 30.0, 33.0, 41.0, 44.0]:
        a.bands = [23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0]
    nside, npix = 16, hp.nside2npix(16)
    rng = np.random.default_rng(a.seed)

    # truth amplitudes at 30 GHz: DR1 polarized synchrotron, 5 deg, nside 16
    dr1 = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    dr1 = hp.ud_grade(dr1, 64)
    dr1 = hp.smoothing(dr1, fwhm=np.radians(5.0), pol=True)
    qu30 = hp.ud_grade(dr1[1:], nside)                                   # (2, npix)
    nus = np.array([CENTRE[b] for b in a.bands])
    rn = real_noise(nside)
    sig = np.stack([np.stack([rn[b][1][1], rn[b][1][2]]) for b in a.bands], axis=2)  # (2, npix, nb)
    sig = sig.reshape(2 * npix, len(nus)) * a.noise_scale
    a30 = qu30.reshape(2 * npix)                                         # pixel-Stokes vector
    nps = a30.size
    w = 1.0 / sig ** 2
    truth_sig = a30[:, None] * (nus[None, :] / 30.0) ** a.beta
    # per-element weight matrix Wp (nps, nb, nb): diag(w), or diag(w) with the CMB and dust columns projected out
    Wp = np.zeros((nps, nus.size, nus.size))
    Wp[:, np.arange(nus.size), np.arange(nus.size)] = w
    if a.nuisance == "cmb_dust":
        fd = dust_sed(nus, DUST_BETA, DUST_T, DUST_REF)[0]
        M = np.stack([cmb_sed(nus), fd], axis=1)                         # (nb, 2), uK_RJ per unit amplitude
        dd = np.array(hp.read_map(DR1_DUST, field=(0, 1, 2)), dtype=float)
        dd = hp.smoothing(hp.ud_grade(dd, 64), fwhm=np.radians(5.0), pol=True)
        d70 = hp.ud_grade(dd[1:], nside).reshape(nps)
        truth_sig = truth_sig + d70[:, None] * fd[None, :]                # the CMB truth is zero
        G = np.einsum("pb,bi,bj->pij", w, M, M)
        if a.cmb_prior_sd > 0:
            G[:, 0, 0] += 1.0 / a.cmb_prior_sd ** 2
        WM = w[:, :, None] * M[None, :, :]
        Wp -= np.einsum("pbi,pij,pcj->pbc", WM, np.linalg.inv(G), WM)
        i353 = int(np.argmax(nus))
        print(f"# nuisance: CMB (prior sd {a.cmb_prior_sd} uK, 0 = flat) and dust (MBB {DUST_BETA}, {DUST_T} K) per "
              f"element; truth dust |a| at {nus[i353]:g} GHz median {np.median(np.abs(d70) * fd[i353]):.2f} uK_RJ, "
              f"at 70 GHz {np.median(np.abs(d70)):.3f}")
    aref = a30 * (a.nu_ref / 30.0) ** a.beta                              # amplitude at nu_ref

    betas = np.round(np.arange(-4.6, -1.6 + 1e-9, 0.01), 2)
    lognu = np.log(nus / a.nu_ref)
    F = np.exp(betas[:, None] * lognu[None, :])                           # (nbeta, nb)

    # scale maps for the zero-mean variants (uK at nu_ref)
    P = np.hypot(qu30[0], qu30[1]) * (a.nu_ref / 30.0) ** a.beta
    np.random.seed(11)
    cl = np.zeros(48); cl[1:] = 1.0
    f = hp.smoothing(hp.synfast(cl, nside, new=True), fwhm=np.radians(10.0))
    f *= 0.3 / f.std()
    g_tpl = np.tile(P * np.exp(f), 2)
    g_orc = np.maximum(np.abs(aref), 1e-3)
    variants = {"zm": (np.ones(nps), np.geomspace(0.05, 2000.0, a.ntau)),
                "zm_tpl": (g_tpl, np.geomspace(0.03, 30.0, a.ntau)),
                "zm_orc": (g_orc, np.geomspace(0.03, 30.0, a.ntau))}

    # global-scale arms: scale maps at nu_ref and the references at which the prior is stated
    gvariants = {}
    if a.global_arm:
        q23 = qu30 * (a.nu_ref / 30.0) ** a.beta
        rms2 = hp.smoothing(0.5 * (q23[0] ** 2 + q23[1] ** 2), fwhm=np.radians(10.0))
        g_true = np.tile(np.sqrt(np.maximum(rms2, 1e-6)), 2)
        blat = np.tile(np.abs(hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]), 2)
        g_hi = g_true * np.where(blat > 30.0, 1.3, 1.0)
        wf = w * ((nus / a.nu_ref) ** a.beta) ** 2
        nu_pivot = a.nu_ref * np.exp((wf * np.log(nus / a.nu_ref)).sum() / wf.sum())
        for mname, g in (("true", g_true), ("hi", g_hi)):
            for rname, ref in (("23", 23.0), ("pv", float(nu_pivot)), ("30", 30.0)):
                gvariants[f"zmG_{mname}_{rname}"] = (g, ref)
        print(f"# global-scale arms: information pivot {nu_pivot:.2f} GHz; scale maps: true local rms (10 deg), and the same "
              f"times 1.3 at |b| > 30 deg", flush=True)
    sgrid_g = np.geomspace(0.1, 10.0, a.ntau)

    # regions and their Fisher sd at the truth (amplitudes profiled)
    ft = (nus / a.nu_ref) ** a.beta
    ut = ft * lognu
    A0 = np.einsum("b,pbc,c->p", ft, Wp, ft)
    I_p = aref ** 2 * (np.einsum("b,pbc,c->p", ut, Wp, ut) - np.einsum("b,pbc,c->p", ut, Wp, ft) ** 2 / A0)
    A = np.einsum("kpc,kc->kp", np.einsum("kb,pbc->kpc", F, Wp), F)      # (nbeta, nps), noise only
    REG = {}
    for ns_r in (2, 4):
        reg = np.tile(hp.ud_grade(np.arange(hp.nside2npix(ns_r), dtype=float), nside).astype(int), 2)
        order = np.argsort(reg, kind="stable")
        starts = np.searchsorted(reg[order], np.arange(hp.nside2npix(ns_r)))
        fsd = 1.0 / np.sqrt(np.add.reduceat(I_p[order], starts))
        nper = np.diff(np.append(starts, nps))
        REG[ns_r] = (order, starts, fsd, nper)

    names = ["flat23", "flat30", "flat60", "prof", "zm", "zm_tpl", "zm_orc"] + list(gvariants)
    res = {ns_r: {n: [] for n in names} for ns_r in REG}
    for r in range(a.nreal):
        d = truth_sig + sig * rng.standard_normal(sig.shape)
        Pd = np.einsum("pbc,pc->pb", Wp, d)
        S = (d * Pd).sum(1)
        b = F @ Pd.T                                                     # (nbeta, nps)
        prof = -0.5 * (S[None, :] - b ** 2 / A)
        flat = prof - 0.5 * np.log(A)
        per = {"prof": prof, "flat23": flat}
        for ns_r, (order, starts, fsd, nper) in REG.items():
            Lp = np.add.reduceat(prof[:, order], starts, axis=1)
            Lf = np.add.reduceat(flat[:, order], starts, axis=1)
            res[ns_r]["prof"].append(summarize(Lp, betas, a.beta))
            for ref, nm in ((23.0, "flat23"), (30.0, "flat30"), (60.0, "flat60")):
                tilt = betas[:, None] * np.log(ref / a.nu_ref) * nper[None, :]
                res[ns_r][nm].append(summarize(Lf + tilt, betas, a.beta))
        for nm, (g, sgrid) in variants.items():
            acc = {ns_r: np.empty((sgrid.size, betas.size, REG[ns_r][1].size)) for ns_r in REG}
            for it, s in enumerate(sgrid):
                t2 = (s * g) ** 2
                q = t2[None, :] * A
                L = -0.5 * (S[None, :] - b ** 2 * t2[None, :] / (1.0 + q)) - 0.5 * np.log1p(q)
                for ns_r, (order, starts, fsd, nper) in REG.items():
                    acc[ns_r][it] = np.add.reduceat(L[:, order], starts, axis=1)
            for ns_r in REG:
                m = acc[ns_r].max(axis=0)
                Lz = m + np.log(np.exp(acc[ns_r] - m[None]).sum(axis=0))
                res[ns_r][nm].append(summarize(Lz, betas, a.beta))
        for nm, (g, ref) in gvariants.items():
            # prior N(0, (s g')^2) on the amplitude at `ref`, g' the map moved with the true index: on the nu_ref
            # amplitude that is a variance (s g)^2 (ref/nu_ref)^(2 (beta_true - beta))
            tilt2 = (ref / a.nu_ref) ** (2.0 * (a.beta - betas))            # (nbeta,)
            acc = {ns_r: np.empty((sgrid_g.size, betas.size, REG[ns_r][1].size)) for ns_r in REG}
            for it, s in enumerate(sgrid_g):
                t2 = (s * g) ** 2 * tilt2[:, None]                          # (nbeta, nps)
                q = t2 * A
                L = -0.5 * (S[None, :] - b ** 2 * t2 / (1.0 + q)) - 0.5 * np.log1p(q)
                for ns_r, (order, starts, fsd, nper) in REG.items():
                    acc[ns_r][it] = np.add.reduceat(L[:, order], starts, axis=1)
            for ns_r in REG:
                m = acc[ns_r].max(axis=1)                                   # (ns, nreg)
                logZ = m + np.log(np.exp(acc[ns_r] - m[:, None, :]).sum(axis=1))
                tot = logZ.sum(axis=1)                                      # (ns,), one scale for the whole sky
                J = acc[ns_r] + (tot[:, None] - logZ)[:, None, :]           # region r given s, the others marginalized
                mj = J.max(axis=0)
                res[ns_r][nm].append(summarize(mj + np.log(np.exp(J - mj[None]).sum(axis=0)), betas, a.beta))
        if r == 0 or (r + 1) % 20 == 0:
            print(f"# realization {r + 1}/{a.nreal} done", file=sys.stderr, flush=True)

    print(f"# E036: truth beta {a.beta}, bands {[float(x) for x in nus]} GHz, nu_ref {a.nu_ref}, {a.nreal} noise "
          f"realizations, noise scale {a.noise_scale}, nuisance {a.nuisance}")
    print(f"# median Q/U noise rms per nside-16 pixel (uK_RJ): {np.round(np.median(sig, axis=0), 3).tolist()}")
    print(f"# truth |a| at nu_ref per pixel-Stokes: median {np.median(np.abs(aref)):.2f}, 16/84 pct "
          f"{np.percentile(np.abs(aref), 16):.2f}/{np.percentile(np.abs(aref), 84):.2f} uK_RJ; "
          f"median S/N per element at the lowest band {np.median(np.abs(aref) * ft[0] / sig[:, 0]):.2f}")
    classes = [("fisher sd < 0.03", 0.0, 0.03), ("0.03 to 0.1", 0.03, 0.1), ("0.1 to 0.25", 0.1, 0.25),
               ("> 0.25", 0.25, np.inf)]
    for ns_r, (order, starts, fsd, nper) in REG.items():
        print(f"\n## regions at nside {ns_r}: {fsd.size} regions, {int(nper[0])} amplitudes each; Fisher sd median "
              f"{np.median(fsd):.3f}, 16/84 pct {np.percentile(fsd, 16):.3f}/{np.percentile(fsd, 84):.3f}")
        for cname, lo, hi in classes:
            sel = (fsd >= lo) & (fsd < hi)
            if not sel.any():
                continue
            print(f"\n# class {cname}: {int(sel.sum())} regions x {a.nreal} realizations; median Fisher sd "
                  f"{np.median(fsd[sel]):.3f}")
            print(f"{'estimator':9s} {'shift':>8s} {'+-':>7s} {'med sd':>7s} {'mean z':>7s} {'+-':>6s} {'sd z':>6s} "
                  f"{'cov68':>6s} {'cov95':>6s} {'edge':>6s}")
            for nm in names:
                arr = np.array(res[ns_r][nm])                             # (nreal, 4, nreg)
                mean, sd, cdf, edge = (arr[:, k][:, sel] for k in range(4))
                z = (mean - a.beta) / sd
                n_eff = z.size
                print(f"{nm:9s} {np.mean(mean - a.beta):+8.4f} {np.std(mean - a.beta) / np.sqrt(n_eff):7.4f} "
                      f"{np.median(sd):7.4f} {z.mean():+7.3f} {z.std() / np.sqrt(n_eff):6.3f} {z.std():6.3f} "
                      f"{np.mean((cdf > 0.16) & (cdf < 0.84)):6.3f} {np.mean((cdf > 0.025) & (cdf < 0.975)):6.3f} "
                      f"{np.mean(edge > 0.01):6.3f}")
        # latitude zones, regions with Fisher sd below 0.25 only
        lat_r = np.abs(hp.pix2ang(ns_r, np.arange(fsd.size), lonlat=True)[1])
        for zname, zsel in (("|b| > 30", (lat_r > 30) & (fsd < 0.25)), ("|b| <= 30", (lat_r <= 30) & (fsd < 0.25))):
            if not zsel.any():
                continue
            print(f"\n# zone {zname}, Fisher sd < 0.25: {int(zsel.sum())} regions; shift (mean posterior mean - truth), its "
                  f"error over realizations, mean z")
            for nm in names:
                arr = np.array(res[ns_r][nm])
                mean, sd = arr[:, 0][:, zsel], arr[:, 1][:, zsel]
                per_real = (mean - a.beta).mean(axis=1)
                print(f"{nm:11s} {per_real.mean():+8.4f} {per_real.std() / np.sqrt(len(per_real)):7.4f} "
                      f"{((mean - a.beta) / sd).mean():+7.3f}")


if __name__ == "__main__":
    main()
