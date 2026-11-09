"""E031: the noise regime of the simulation ladder against real maps.

Compares the per-pixel white noise of the PySM intensity ladder (`pysm_intensity_sim.py`, noise
scale 1 and 0.03) with the white-noise rms of the real WMAP, LFI and HFI maps averaged to nside
64, and evaluates the per-region Fisher sd of the spectral indices at the truth for

    d_pb = sum_c a_pc f_c(nu_b; theta_p) + n_pb,   n_pb ~ N(0, sigma_pb^2),

with the amplitudes either known or profiled under flat priors: per pixel the Fisher matrix of
(amplitudes, indices) is reduced to the indices by the Schur complement
S_p = F_bb - F_ab^T F_aa^{-1} F_ab, summed over the pixels of a region and inverted. For uniform
noise this is `fisher_sd` / `fisher_sd_profiled` of `pysm_intensity_sim.py`, whose stored values
the script reproduces first. In polarization Q and U share one index per region and their Schur
complements add.

    python -I nextgen/scripts/e031_noise_regime_fisher.py --out nextgen/experiments/E031/
"""
import argparse
import os

import healpy as hp
import numpy as np
from astropy.io import fits

H_K = 0.0479924466   # K / GHz
TCMB = 2.7255
RMS16 = np.array([70.0, 60.0, 60.0, 55.0, 60.0, 60.0, 70.0, 40.0, 24.0, 16.0, 24.0, 80.0])
FREQS = np.array([23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 94.0, 100.0, 143.0, 217.0, 353.0])
Q1IN = "/Users/trygvels/data/q1_inputs"
# label, ladder band (GHz), band centre (GHz), rms files (differencing assemblies)
REAL = [
    ("WMAP K", 23.0, 22.8, [f"{Q1IN}/CG_023-WMAP_K_IQU_rms_n0512_nodip_v1.fits"]),
    ("LFI 30", 30.0, 28.4, [f"{Q1IN}/CG_030_IQU_rms_n0512_nodip_v1.fits"]),
    ("WMAP Ka", 33.0, 33.0, [f"{Q1IN}/CG_030-WMAP_Ka_IQU_rms_n0512_nodip_v1.fits"]),
    ("WMAP Q", 41.0, 40.6, [f"{Q1IN}/CG_040-WMAP_Q{i}_IQU_rms_n0512_nodip_v1.fits" for i in (1, 2)]),
    ("LFI 44", 44.0, 44.1, [f"{Q1IN}/CG_044_IQU_rms_n0512_nodip_v1.fits"]),
    ("WMAP V", 61.0, 60.8, [f"{Q1IN}/CG_060-WMAP_V{i}_IQU_rms_n0512_nodip_v1.fits" for i in (1, 2)]),
    ("LFI 70", 70.0, 70.1, [f"{Q1IN}/CG_070_IQU_rms_n1024_nodip_v1.fits"]),
    ("WMAP W", 94.0, 93.5, [f"{Q1IN}/CG_090-WMAP_W{i}_IQU_rms_n0512_nodip_v1.fits"
                            for i in (1, 2, 3, 4)]),
    ("HFI 100", 100.0, 100.0, ["/Users/trygvels/data/planck_pr3/n512/pr3_100_IQU_rms_n0512.fits"]),
    ("HFI 353 (NPIPE A)", 353.0, 353.0,
     ["/Users/trygvels/data/beyondplanck/npipe/n512/npipe_353_IQU_rms_n0512.fits"]),
]
# Planck 2018 I Table 4 temperature sensitivities, uK_CMB deg (quoted from memory; 100 and 353 are
# checked against the local maps before 143 and 217 are used).
PLANCK18_T = {100.0: 1.29, 143.0: 0.55, 217.0: 0.78, 353.0: 2.56}
LFI = (30.0, 44.0, 70.0)


def cmb_sed(nu):
    """uK_CMB to uK_RJ per band (the CMB mixing column)."""
    x = H_K * np.asarray(nu, dtype=float) / TCMB
    return x ** 2 * np.exp(x) / np.expm1(x) ** 2


def dust_sed(nu, beta, T, nu0):
    """(npix, B) modified blackbody in RJ units at pivot nu0: (nu/nu0)^(beta+1) expm1(x0)/expm1(x)."""
    nu = np.asarray(nu, dtype=float)
    beta, T = np.atleast_1d(beta).astype(float), np.atleast_1d(T).astype(float)
    x, x0 = H_K * nu[None, :] / T[:, None], H_K * nu0 / T[:, None]
    return (nu[None, :] / nu0) ** (beta[:, None] + 1.0) * np.expm1(x0) / np.expm1(x)


def white_rms(path, nside_out=64):
    """(3, npix_out) white-noise rms in uK_CMB of the nside_out average of an (I, Q, U) rms map:
    sigma^2 = sum_i sigma_i^2 / n^2 over the n observed sub-pixels (rms > 0); inf if none."""
    unit = fits.open(path)[1].header.get("TUNIT1", "")
    scale = 1000.0 if "mK" in unit else 1.0
    rms = np.array(hp.read_map(path, field=(0, 1, 2)), dtype=float) * scale
    n = (hp.get_nside(rms[0]) // nside_out) ** 2
    out = np.empty((3, hp.nside2npix(nside_out)))
    for i in range(3):
        v = hp.reorder(rms[i] ** 2, r2n=True).reshape(-1, n)
        good = v > 0
        ng = good.sum(axis=1)
        s = np.where(good, v, 0.0).sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[i] = hp.reorder(np.sqrt(np.where(ng > 0, s / np.maximum(ng, 1) ** 2, np.inf)),
                                n2r=True)
    return out


def real_noise(nside=64):
    """{ladder band: (label, (3, npix) rms in uK_RJ)} with differencing assemblies combined by
    inverse variance."""
    out = {}
    for label, band, centre, files in REAL:
        w = sum(1.0 / white_rms(fn, nside) ** 2 for fn in files)
        out[band] = (label, cmb_sed(centre) / np.sqrt(w))
    return out


def fisher_sd(M, dM, w, region, nreg, profiled):
    """M (npix, B, K) amplitude columns, dM (npix, B, I) index columns a_c dsed_c/dbeta_c, w (npix,
    B) inverse noise variances. Returns the (nreg, I) Fisher sd of the region indices. Profiled:
    per-pixel Schur complement onto the indices, summed over the region, inverted (the other
    index marginalized). Known: 1 / sqrt(sum F_ii) per index, the convention of `fisher_sd` in
    pysm_intensity_sim.py."""
    Fbb = np.einsum("pbi,pb,pbj->pij", dM, w, dM)
    if profiled:
        Faa = np.einsum("pbi,pb,pbj->pij", M, w, M)
        Fab = np.einsum("pbi,pb,pbj->pij", M, w, dM)
        S = Fbb - np.einsum("pji,pjk->pik", Fab, np.linalg.solve(Faa, Fab))
    else:
        S = Fbb
    Sr = np.zeros((nreg,) + S.shape[1:])
    np.add.at(Sr, region, S)
    if not profiled:
        with np.errstate(divide="ignore"):
            return 1.0 / np.sqrt(np.diagonal(Sr, axis1=1, axis2=2))
    out = np.full(Sr.shape[:2], np.inf)
    ok = np.linalg.det(Sr) > 0
    out[ok] = np.sqrt(np.diagonal(np.linalg.inv(Sr[ok]), axis1=1, axis2=2))
    return out


def schur_scalar(M, d, w):
    """Per-pixel profiled Fisher of one index: d^T W d - d^T W M (M^T W M)^{-1} M^T W d, with M
    (npix, B, K), d (npix, B), w (npix, B)."""
    Faa = np.einsum("pbi,pb,pbj->pij", M, w, M)
    Fab = np.einsum("pbi,pb,pb->pi", M, w, d)
    return np.einsum("pb,pb,pb->p", d, w, d) - np.einsum("pi,pi->p", Fab,
                                                         np.linalg.solve(Faa, Fab[..., None])[..., 0])


def regions(nb, nside=64):
    if nb == 0:
        return np.zeros(hp.nside2npix(nside), int), 1
    return hp.ud_grade(np.arange(12 * nb * nb).astype(float), nside, order_in="RING",
                       order_out="RING").astype(int), 12 * nb * nb


def qmed(x, order):
    """Medians of x in the four quartiles of `order` (ascending: Q1 the smallest `order`)."""
    idx = np.argsort(order)
    return [float(np.median(x[q])) for q in np.array_split(idx, 4)]


def fmt(v):
    return " ".join(f"{x:8.4f}" for x in v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", default="/Users/trygvels/data/e023_pysm_n1/truth.npz")
    ap.add_argument("--out", default="nextgen/experiments/E031/")
    ap.add_argument("--low-nu", type=float, default=11.0, help="low-frequency channel of set C, GHz")
    ap.add_argument("--low-rms", type=float, default=45.0,
                    help="white Q/U noise of the set C channel per nside-64 pixel, uK_RJ")
    ap.add_argument("--haslam-rms", default="/Users/trygvels/data/beyondplanck/"
                    "haslam_rms_n0512_60arcmin_uK_1pct_0.8K.fits")
    ap.add_argument("--low-dec-min", type=float, default=-30.0,
                    help="southern declination limit of the set C channel's footprint, deg")
    ap.add_argument("--tag", default="fisher_report")
    a = ap.parse_args()
    import pysm3
    import pysm3.units as u
    os.makedirs(a.out, exist_ok=True)
    R = []
    t = np.load(a.truth, allow_pickle=True)
    nside, npix = int(t["nside"]), t["a_sync"].size
    side = np.degrees(np.sqrt(hp.nside2pixarea(nside)))
    lat = hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]
    reg16, n16 = regions(16)

    # ---- Table 1: per-band white noise at nside 64 ------------------------------------------
    rn = real_noise(nside)
    R.append(f"E031 report. nside {nside}, pixel side {side:.4f} deg. All noise values are white-noise rms per "
             f"nside-{nside} pixel in uK_RJ.")
    R.append("")
    R.append("Table 1. Real white noise (median over the sky) against the ladder, intensity and Q/U.")
    R.append("band               real I   real QU   ladder x1   ratio I   ladder x0.03  ratio I   "
             "x1 / real QU")
    for i, f in enumerate(FREQS):
        if f not in rn:
            continue
        lab, r = rn[f]
        mi, mq = np.median(r[0]), np.median(np.sqrt(0.5 * (r[1] ** 2 + r[2] ** 2)))
        s1 = RMS16[i] * nside / 16.0
        R.append(f"{lab:18s} {mi:7.3f} {mq:8.3f} {s1:9.1f} {s1 / mi:10.1f} {0.03 * s1:11.2f} "
                 f"{0.03 * s1 / mi:9.2f} {s1 / mq:10.1f}")
    chk = {}
    for f in (100.0, 353.0):
        chk[f] = np.median(rn[f][1][0]) / cmb_sed(f) * side
    R.append(f"HFI check against Planck 2018 I Table 4 (uK_CMB deg, temperature): 100 GHz local "
             f"{chk[100.0]:.3f} against 1.29; 353 GHz local (NPIPE A split) {chk[353.0]:.3f} against 2.56 "
             f"(full mission; 3.6 if the split has 1.4 times its noise)")
    use_table = abs(chk[100.0] / 1.29 - 1) < 0.3 and (abs(chk[353.0] / 2.56 - 1) < 0.3
                                                    or abs(chk[353.0] / 3.6 - 1) < 0.3)
    R.append(f"143 and 217 GHz from the table: {'used' if use_table else 'NOT used (check failed)'}"
             f" (spatial pattern of the 100 GHz rms map)")
    sig_real = np.full((npix, FREQS.size), np.inf)
    for i, f in enumerate(FREQS):
        if f in rn:
            sig_real[:, i] = rn[f][1][0]
        elif use_table:
            pat = rn[100.0][1][0] / np.median(rn[100.0][1][0])
            sig_real[:, i] = PLANCK18_T[f] / side * cmb_sed(f) * pat
            s1 = RMS16[i] * nside / 16.0
            R.append(f"HFI {int(f)} (table)    {np.median(sig_real[:, i]):7.3f}      n/a {s1:9.1f} "
                     f"{s1 / np.median(sig_real[:, i]):10.1f} {0.03 * s1:11.2f} "
                     f"{0.03 * s1 / np.median(sig_real[:, i]):9.2f}")
    R.append("")

    # ---- Intensity: reproduce the stored Fisher values, then the noise cases ------------------
    beta_s, beta_d, T_fixed = t["beta_s"], t["beta_d"], float(t["T_fixed"])
    a_s, a_d = t["a_sync"], t["a_dust"]
    sky_fa = pysm3.Sky(nside=nside, preset_strings=["f1", "a1"])

    def known_sed(comp, nus):
        """(npix, len(nus)) SED of a PySM component normalized at 30 GHz, its 30 GHz amplitude, and
        the number of pixels with zero emission at 30 GHz (those get the median SED of the others;
        a frequency outside the component's tabulated range gets zero emission)."""
        em = []
        for f in list(nus) + [30.0]:
            try:
                em.append(comp.get_emission(f * u.GHz)[0].to(
                    u.uK_RJ, equivalencies=u.cmb_equivalencies(f * u.GHz)).value)
            except Exception:
                em.append(np.zeros(npix))
        em = np.stack(em, axis=1)
        ref, ok = em[:, -1], em[:, -1] > 0
        sed = np.empty_like(em[:, :-1])
        sed[ok] = em[ok, :-1] / ref[ok, None]
        sed[~ok] = np.median(sed[ok], axis=0)
        return sed, ref, int((~ok).sum())

    def intensity_model(nus):
        """Mixing columns at the truth: M3 (cmb, sync, dust), M5 (M3, free-free, AME), and dM
        (a_s dsed_s/dbeta_s, a_d dsed_d/dbeta_d), each (npix, len(nus), K)."""
        nus = np.asarray(nus, dtype=float)
        sed_s = (nus[None, :] / 30.0) ** beta_s[:, None]
        sed_d = dust_sed(nus, beta_d, np.full(npix, T_fixed), 545.0)
        M3 = np.stack([np.broadcast_to(cmb_sed(nus), sed_s.shape), sed_s, sed_d], axis=2)
        dM = np.stack([a_s[:, None] * sed_s * np.log(nus / 30.0)[None, :],
                       a_d[:, None] * sed_d * np.log(nus / 545.0)[None, :]], axis=2)
        M5 = np.concatenate([M3] + [known_sed(c, nus)[0][:, :, None] for c in sky_fa.components],
                            axis=2)
        return M3, M5, dM

    M3, M5, dM = intensity_model(FREQS)
    sig1 = np.broadcast_to(RMS16 * nside / 16.0, (npix, FREQS.size))
    w1 = 1.0 / sig1 ** 2
    k1 = fisher_sd(M3, dM, w1, reg16, n16, False)
    p1 = fisher_sd(M3, dM, w1, reg16, n16, True)
    R.append("Check: stored truth.npz Fisher sd at noise scale 1 reproduced, max relative difference: "
             f"sync known {np.max(np.abs(k1[:, 0] / t['fisher_sd_sync_16'] - 1)):.1e}, sync profiled "
             f"{np.max(np.abs(p1[:, 0] / t['fisher_sd_sync_16_profiled'] - 1)):.1e}, dust known "
             f"{np.max(np.abs(k1[:, 1] / t['fisher_sd_dust_16'] - 1)):.1e}, dust profiled "
             f"{np.max(np.abs(p1[:, 1] / t['fisher_sd_dust_16_profiled'] - 1)):.1e}")
    order_s, order_d = k1[:, 0], k1[:, 1]      # quartiles by the scale 1 amplitude-known sd
    R.append("")
    for comp, name in zip(sky_fa.components, ("free-free (f1)", "AME (a1, two populations summed)")):
        sed, ref, nzero = known_sed(comp, [23.0, 0.408])
        R.append(f"PySM {name}: median amplitude at 30 GHz {np.median(ref):.2f} uK_RJ (sync "
                 f"{np.median(a_s):.2f}); pixels with zero emission at 30 GHz {nzero}; median SED at 0.408 "
                 f"GHz relative to 30 GHz {np.median(sed[:, 1]):.3g}")
    R.append(f"Model I5 has {M5.shape[2]} amplitude columns per pixel (cmb, sync, dust, free-free, AME), "
             "free-free and AME SEDs known per pixel")
    cases = [("x1", M3, M5, dM, w1), ("x0.03", M3, M5, dM, w1 / 0.03 ** 2)]
    if np.all(np.isfinite(sig_real)):
        cases.append(("real", M3, M5, dM, 1.0 / sig_real ** 2))
        # post hoc: the 408 MHz map, effective rms taken per nside-64 pixel (not averaged down)
        sig408 = hp.ud_grade(np.array(hp.read_map(a.haslam_rms, field=0), dtype=float), nside)
        s408 = a_s * (0.408 / 30.0) ** beta_s
        R.append(f"Haslam 408 MHz: effective rms per nside-{nside} pixel median {np.median(sig408) / 1e6:.3f} K; "
                 f"synchrotron truth extrapolated to 408 MHz median {np.median(s408) / 1e6:.2f} K; median S/N "
                 f"{np.median(s408 / sig408):.1f}")
        nusH = np.concatenate([FREQS, [0.408]])
        M3h, M5h, dMh = intensity_model(nusH)
        wH = np.concatenate([1.0 / sig_real ** 2, 1.0 / sig408[:, None] ** 2], axis=1)
        cases.append(("real+408", M3h, M5h, dMh, wH))
    else:
        R.append("real-noise intensity case skipped: bands without noise")
    R.append("")
    R.append("Table 2. Intensity: Fisher sd of the region indices at nside_beta 16 (3072 regions). Quartile "
             "medians, Q1 the regions with the smallest scale 1 amplitude-known sd of that index. The "
             "real+408 rows are post hoc (models I3H and I5H).")
    R.append("noise     model index  kind      Q1       Q2       Q3       Q4      median   frac<0.1 frac<truth sd")
    keep = {}
    for cname, Ma, Mb, dMc, w in cases:
        for mname, M in (("I3", Ma), ("I5", Mb)):
            k = fisher_sd(M, dMc, w, reg16, n16, False)
            p = fisher_sd(M, dMc, w, reg16, n16, True)
            keep[(cname, mname)] = (k, p)
            for j, (iname, order, tsd) in enumerate((("sync", order_s, beta_s.std()),
                                                     ("dust", order_d, beta_d.std()))):
                for kind, v in (("known", k[:, j]), ("profiled", p[:, j])):
                    if mname == "I5" and kind == "known":
                        continue
                    R.append(f"{cname:9s} {mname:5s} {iname:5s}  {kind:8s} {fmt(qmed(v, order))} "
                             f"{np.median(v):8.4f}   {np.mean(v < 0.1):6.3f}   {np.mean(v < tsd):6.3f}")
    R.append(f"truth sd of the PySM index maps: sync {beta_s.std():.4f}, dust {beta_d.std():.4f}")
    R.append("")
    R.append("Table 2b. Intensity, profiled sd of beta_s, median over regions against region size.")
    R.append("noise     model  nside_beta:  whole sky      1        2        4        8       16")
    for cname, Ma, Mb, dMc, w in cases:
        for mname, M in (("I3", Ma), ("I5", Mb)):
            row = []
            for nb in (0, 1, 2, 4, 8, 16):
                rg, nr = regions(nb)
                row.append(np.median(fisher_sd(M, dMc, w, rg, nr, True)[:, 0]))
            R.append(f"{cname:9s} {mname:5s}              " + " ".join(f"{x:8.5f}" for x in row))
    R.append("")
    if ("real", "I5") in keep:
        # which amplitude degeneracy costs the synchrotron index: drop one known-SED column at a time
        R.append("Table 2c (post hoc). Real noise, profiled sd of beta_s per nside-16 region, median, by model:")
        for label, cols in (("I3 (cmb, sync, dust)", [0, 1, 2]), ("I3 + free-free", [0, 1, 2, 3]),
                            ("I3 + AME", [0, 1, 2, 4]), ("I5 (both)", [0, 1, 2, 3, 4])):
            v = fisher_sd(M5[:, :, cols], dM, 1.0 / sig_real ** 2, reg16, n16, True)[:, 0]
            vh = fisher_sd(M5h[:, :, cols], dMh, wH, reg16, n16, True)[:, 0]
            R.append(f"  {label:22s} without 408 MHz {np.median(v):8.4f}   with 408 MHz {np.median(vh):8.4f}")
        R.append("")
    snr = lambda sig, f, i: np.median(a_s * (f / 30.0) ** beta_s / sig[:, i])
    R.append("Median synchrotron S/N per nside-64 pixel, intensity: 23 GHz ladder x1 "
             f"{snr(sig1, 23.0, 0):.2f}, real {snr(sig_real, 23.0, 0):.1f}; 30 GHz ladder x1 "
             f"{snr(sig1, 30.0, 1):.2f}, real {snr(sig_real, 30.0, 1):.1f}")
    R.append("")

    # ---- Polarization --------------------------------------------------------------------------
    sky_p = pysm3.Sky(nside=nside, preset_strings=["s1", "d1"])
    sync, dust = None, None
    for comp in sky_p.components:
        if hasattr(comp, "pl_index"):
            sync = comp
        else:
            dust = comp
    qu_s = sync.get_emission(30.0 * u.GHz)[1:3].to(
        u.uK_RJ, equivalencies=u.cmb_equivalencies(30.0 * u.GHz)).value
    qu_d = dust.get_emission(353.0 * u.GHz)[1:3].to(
        u.uK_RJ, equivalencies=u.cmb_equivalencies(353.0 * u.GHz)).value
    bs = np.asarray(getattr(sync.pl_index, "value", sync.pl_index), dtype=float)
    bd = np.asarray(getattr(dust.mbb_index, "value", dust.mbb_index), dtype=float)
    Td = np.asarray(dust.mbb_temperature.to(u.K).value, dtype=float)
    P_s = np.hypot(*qu_s)
    R.append(f"Polarization truth (PySM s1, d1): P_sync at 30 GHz median {np.median(P_s):.2f} uK_RJ "
             f"(|b| > 30: {np.median(P_s[np.abs(lat) > 30]):.2f}), at 23 GHz median "
             f"{np.median(P_s * (23 / 30.0) ** bs):.2f}; P_dust at 353 GHz median {np.median(np.hypot(*qu_d)):.2f}")
    # amplitude cross-check with the Cosmoglobe DR1 synchrotron posterior mean (30 GHz, 60 arcmin)
    try:
        dr1 = np.array(hp.read_map("/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits",
                                   field=(0, 1, 2)), dtype=float)
        dr1 = hp.ud_grade(dr1, nside)
        fw = np.radians(5.0)
        dr1s = hp.smoothing(dr1, fwhm=fw, pol=True)
        pys = hp.smoothing(np.array([np.zeros(npix), qu_s[0], qu_s[1]]), fwhm=fw, pol=True)
        ratio = np.hypot(dr1s[1], dr1s[2]) / np.hypot(pys[1], pys[2])
        hi = np.abs(lat) > 30
        R.append("Cross-check, DR1 synch P over PySM s1 P at 30 GHz, both smoothed to 5 deg: median "
                 f"{np.median(ratio):.2f} (16/84 pct {np.percentile(ratio, 16):.2f}/"
                 f"{np.percentile(ratio, 84):.2f}); |b| > 30: {np.median(ratio[hi]):.2f}")
    except Exception as e:                                    # the cross-check is optional
        R.append(f"DR1 cross-check failed: {e!r}")

    def pol_set(bands, low=None, lfi_var=1.0, cmb=True):
        """Per-pixel profiled and known Fisher of beta_s summed over Q and U for a band set."""
        nus = [b for b in bands]
        sig = [np.sqrt(0.5 * (rn[b][1][1] ** 2 + rn[b][1][2] ** 2)) * (np.sqrt(lfi_var) if b in LFI
                                                                         else 1.0) for b in bands]
        if low is not None:
            nus.append(low[0])
            sig.append(np.full(npix, low[1]))
        nus = np.array(nus)
        w = 1.0 / np.stack(sig, axis=1) ** 2
        fs = (nus[None, :] / 30.0) ** bs[:, None]
        fd = dust_sed(nus, bd, Td, 353.0)
        cols = [fs, fd] + ([np.broadcast_to(cmb_sed(nus), fs.shape)] if cmb else [])
        M = np.stack(cols, axis=2)
        Sp, Sk = np.zeros(npix), np.zeros(npix)
        for s in range(2):
            d = qu_s[s][:, None] * fs * np.log(nus / 30.0)[None, :]
            Sp += schur_scalar(M, d, w)
            Sk += np.einsum("pb,pb,pb->p", d, w, d)
        return Sp, Sk

    A = [30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0]
    B = [23.0] + A
    low = (a.low_nu, a.low_rms)
    sets = [("A  BP XV set (no K)", pol_set(A)), ("B  A + K band", pol_set(B)),
            (f"C  B + {a.low_nu:g} GHz at {a.low_rms:g}", pol_set(B, low)),
            ("B  no CMB column", pol_set(B, cmb=False)),
            ("B  LFI variance x 2.65", pol_set(B, lfi_var=2.65)),
            ("A  LFI variance x 2.65", pol_set(A, lfi_var=2.65))]
    rsum = lambda x, rg, nr: np.bincount(rg, x, minlength=nr)
    orderB = 1.0 / np.sqrt(rsum(sets[1][1][0], reg16, n16))
    R.append("")
    R.append("Table 3. Polarization: Fisher sd of beta_s per nside-16 region (Q and U share the index; "
             "sync and dust amplitudes profiled, dust SED known; CMB column with a flat prior unless stated). "
             "Quartiles by the set B profiled sd.")
    R.append("set                         kind      Q1       Q2       Q3       Q4      median   frac<0.1 "
             "frac<0.3  whole sky")
    for name, (Sp, Sk) in sets:
        for kind, S in (("profiled", Sp), ("known", Sk)):
            with np.errstate(divide="ignore"):
                v = 1.0 / np.sqrt(rsum(S, reg16, n16))
            R.append(f"{name:27s} {kind:8s} {fmt(qmed(v, orderB))} {np.median(v):8.4f}   "
                     f"{np.mean(v < 0.1):6.3f}   {np.mean(v < 0.3):6.3f}  {1 / np.sqrt(S.sum()):8.5f}")
    # post hoc: set C restricted to the low-frequency survey's footprint (regions fully north of
    # the declination limit); the map coordinates are Galactic
    th, ph = hp.pix2ang(nside, np.arange(npix))
    dec = 90.0 - np.degrees(hp.Rotator(coord=["G", "C"])(th, ph)[0])
    inside = np.bincount(reg16, (dec > a.low_dec_min).astype(float), minlength=n16) == np.bincount(
        reg16, minlength=n16)
    for name, (Sp, Sk) in sets[1:3]:
        v = 1.0 / np.sqrt(rsum(Sp, reg16, n16))[inside]
        R.append(f"{name:27s} profiled, regions fully at dec > {a.low_dec_min:g} deg ({inside.mean():.2f} of the "
                 f"sky): median {np.median(v):.4f}, frac<0.1 {np.mean(v < 0.1):.3f}, frac<0.3 {np.mean(v < 0.3):.3f}")
    R.append("")
    R.append("Table 3b. Polarization, profiled sd of beta_s against region size: median over regions and "
             "the fraction of regions below 0.1.")
    R.append("set                         nside_beta:      1              2              4              8"
             "             16")
    for name, (Sp, Sk) in sets[:3] + sets[4:5]:
        row = []
        for nb in (1, 2, 4, 8, 16):
            rg, nr = regions(nb)
            v = 1.0 / np.sqrt(rsum(Sp, rg, nr))
            row.append(f"{np.median(v):7.4f} ({np.mean(v < 0.1):4.2f})")
        R.append(f"{name:27s}             " + " ".join(row))
    R.append("")
    nq = {b: np.sqrt(0.5 * (rn[b][1][1] ** 2 + rn[b][1][2] ** 2)) for b in (23.0, 30.0)}
    R.append("Median polarized synchrotron S/N per nside-64 pixel (P over the Q/U rms): 23 GHz "
             f"{np.median(P_s * (23 / 30.0) ** bs / nq[23.0]):.2f}, 30 GHz {np.median(P_s / nq[30.0]):.2f}; "
             f"at |b| > 30: 23 GHz {np.median((P_s * (23 / 30.0) ** bs / nq[23.0])[np.abs(lat) > 30]):.2f}, "
             f"30 GHz {np.median((P_s / nq[30.0])[np.abs(lat) > 30]):.2f}")
    kB = 1.0 / np.sqrt(rsum(sets[1][1][0], reg16, n16))
    R.append(f"Prediction 6: ladder x1 sync median sd, known {np.median(k1[:, 0]):.3f}, profiled "
             f"{np.median(p1[:, 0]):.3f}; polarization set B profiled median {np.median(kB):.3f}; "
             + (f"real-noise intensity I3 profiled median {np.median(keep[('real', 'I3')][1][:, 0]):.4f}"
                if ('real', 'I3') in keep else ""))
    out = os.path.join(a.out, a.tag + ".txt")
    with open(out, "w") as f:
        f.write("\n".join(R) + "\n")
    print("\n".join(R))
    print(f"written {out}")


if __name__ == "__main__":
    main()
