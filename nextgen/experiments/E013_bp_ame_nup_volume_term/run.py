"""E013 / Q14: per-pixel AME peak-frequency posteriors on BeyondPlanck data under four lnL variants.

Follows the construction of Commander3's local spectral-index sampler at smoothing scale 3
(nside 16, 60 arcmin FWHM): every component amplitude and every band map is smoothed to 60 arcmin
and degraded to nside 16, and the mixing matrix is then evaluated with the nside-16 values of the
spectral parameters. The non-AME components of the BeyondPlanck posterior-mean sky (BP_mean_v2.h5)
are subtracted from each band in the nu_p frequency range, and per pixel the four one-dimensional
posteriors of E013_bp_ame_nup_volume_term.md are evaluated on a grid of nu_p: chisq (amplitude fixed
at the BP mean), ridge, marginal (ridge - 1/2 ln A) and code-marginal (ridge + 1/2 ln A, what
Commander3 BP10 evaluated), each times the BP prior N(21, 3) GHz on [10, 70] GHz.

Approximations (stated in the record): delta bandpasses at the BP nominal frequencies; the public
BP7 smoothed rms maps instead of the unreleased BP10 ones; Gaussian band beams; per-band monopole
and dipole fitted to the residual at |b| > 30 deg instead of the sampled BP values; radio sources
and the relativistic quadrupole ignored.

    python -I nextgen/experiments/E013_bp_ame_nup_volume_term/run.py --out <dir>
"""
import argparse
import os

import h5py
import healpy as hp
import numpy as np

BP = "/Users/trygvels/data/beyondplanck"
H_K = 0.04799243073366221  # h/k in K per GHz
T_CMB = 2.7255
NSIDE_LO, FWHM_LO = 16, 60.0
NSIDE_WORK = 256                       # resolution at which smoothing is applied
NU_MIN, NU_MAX = 15.0, 100.0           # COMP_NU_P_NU_MIN/MAX07
PRIOR_MEAN, PRIOR_RMS = 21.0, 3.0      # COMP_PRIOR_GAUSS_NU_P_*07
BOUNDS = (10.0, 70.0)                  # COMP_PRIOR_UNI_NU_P_*07
TEMPLATE = "/Users/trygvels/data/templates/spdust2_cnm_pysm2_emissivity.txt"

# label, nominal GHz, native fwhm arcmin, map file, scale to uK, smoothed rms file, rms scale
BANDS = [
    ("030", 28.4, 32.3, f"{BP}/BP_030_IQU_n0512_v2.fits", 1.0,
     f"{BP}/BP7_030_rms_n0016_60arcmin_uK_true.fits", 1.0),
    ("044", 44.1, 27.0, f"{BP}/BP_044_IQU_n0512_v2.fits", 1.0,
     f"{BP}/BP7_044_rms_n0016_60arcmin_uK_true.fits", 1.0),
    ("070", 70.1, 13.2, f"{BP}/BP_070_IQU_n1024_v2.fits", 1.0,
     f"{BP}/BP7_070_rms_n0016_60arcmin_uK_true.fits", 1.0),
    ("Ka", 33.0, 39.6, f"{BP}/WMAP9_Ka1_pol_n512_mK_dipole.fits", 1000.0,
     f"{BP}/WMAP_Ka1_rms_n0016_60arcmin_mK_true.fits", 1000.0),
    ("Q1", 40.6, 30.6, f"{BP}/WMAP9_Q1_pol_n512_uK_dipole.fits", 1.0,
     f"{BP}/WMAP_Q1_rms_n0016_60arcmin_uK_true.fits", 1.0),
    ("Q2", 40.6, 30.6, f"{BP}/WMAP9_Q2_pol_n512_uK_dipole.fits", 1.0,
     f"{BP}/WMAP_Q2_rms_n0016_60arcmin_uK_true.fits", 1.0),
    ("V1", 60.8, 21.0, f"{BP}/WMAP9_V1_pol_n512_uK_dipole.fits", 1.0,
     f"{BP}/WMAP_V1_rms_n0016_60arcmin_uK_true.fits", 1.0),
    ("V2", 60.8, 21.0, f"{BP}/WMAP9_V2_pol_n512_uK_dipole.fits", 1.0,
     f"{BP}/WMAP_V2_rms_n0016_60arcmin_uK_true.fits", 1.0),
]


def cmb_to_rj(nu):
    x = H_K * nu / T_CMB
    return x**2 * np.exp(x) / np.expm1(x) ** 2


def sed_synch(nu, beta, nu_ref=0.408):  # BP intensity reference is 408 MHz (BP_mean parameters/synch/nu_ref)
    return (nu / nu_ref) ** beta


def sed_dust(nu, beta, T, nu_ref=857.0):
    return (nu / nu_ref) ** (beta + 1.0) * np.expm1(H_K * nu_ref / T) / np.expm1(H_K * nu / T)


def sed_ff(nu, Te, nu_ref=40.0):
    """Commander3 free-free SED in RJ units relative to nu_ref (Draine 2011 Gaunt factor)."""
    def g(nu_):
        return np.log(np.exp(5.960 - np.sqrt(3.0) / np.pi * np.log(nu_ * (Te / 1e4) ** -1.5)) + np.e)
    return g(nu) / g(nu_ref) * (nu_ref / nu) ** 2 * np.exp(-H_K * (nu - nu_ref) / Te)


class SpDust:
    """SpDust2 CNM template (peak 30 GHz) shifted in log-frequency to peak at nu_p, RJ units,
    normalized at nu_ref, as in Commander3's comm_spindust_comp_mod."""

    def __init__(self, path, nu_ref=22.0, nu_p0=30.0):
        f, j = np.loadtxt(path).T
        self.lf, self.lj = np.log(f), np.log(j)
        self.nu_ref, self.nu_p0 = nu_ref, nu_p0

    def _j(self, nu):
        return np.exp(np.interp(np.log(nu), self.lf, self.lj, left=-np.inf, right=-np.inf))

    def sed(self, nu, nu_p):
        nu, nu_p = np.broadcast_arrays(np.asarray(nu, float), np.asarray(nu_p, float))
        s = self.nu_p0 / nu_p
        return self._j(nu * s) / self._j(self.nu_ref * s) * (self.nu_ref / nu) ** 2


def c3_alm(h5, comp):
    """Commander3 real-alm layout (index l^2 + l + m, m from -l to l) to healpy complex alms."""
    g = h5["000000"][comp]
    lmax = int(g["amp_lmax"][()])
    flat = np.asarray(g["amp_alm"][0], dtype=np.float64)
    assert flat.size == (lmax + 1) ** 2
    l, m = hp.Alm.getlm(lmax)
    alm = flat[l**2 + l + m].astype(np.complex128)
    pos = m > 0
    alm[pos] = (flat[(l**2 + l + m)[pos]] + 1j * flat[(l**2 + l - m)[pos]]) / np.sqrt(2.0)
    return alm, lmax


def comp_lo(h5, comp):
    """Component amplitude smoothed to FWHM_LO and degraded to NSIDE_LO (the stored alms carry no
    beam; the component's reference beam is irrelevant because we apply the full 60 arcmin)."""
    alm, lmax = c3_alm(h5, comp)
    bl = hp.gauss_beam(np.radians(FWHM_LO / 60.0), lmax=lmax)
    alm = hp.almxfl(alm, bl)
    return hp.ud_grade(hp.alm2map(alm, NSIDE_WORK, lmax=lmax), NSIDE_LO)


def band_lo(path, native_fwhm, scale):
    """Band map in uK (file units times scale), beam-deconvolved and re-smoothed to FWHM_LO, at
    NSIDE_LO. The deconvolution is regularized by the 60 arcmin beam (ratio of Gaussian beams)."""
    m = hp.read_map(path, field=0) * scale
    m = hp.ud_grade(m, NSIDE_WORK)
    lmax = 3 * NSIDE_WORK - 1
    alm = hp.map2alm(m, lmax=lmax, iter=1)
    ratio = hp.gauss_beam(np.radians(FWHM_LO / 60.0), lmax) / hp.gauss_beam(
        np.radians(native_fwhm / 60.0), lmax)
    return hp.ud_grade(hp.alm2map(hp.almxfl(alm, ratio), NSIDE_WORK, lmax=lmax), NSIDE_LO)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)) + "/out")
    ap.add_argument("--ngrid", type=int, default=1200)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    h5 = h5py.File(f"{BP}/BP_mean_v2.h5", "r")
    par = h5["000000"]
    npix = hp.nside2npix(NSIDE_LO)
    lat = 90.0 - np.degrees(hp.pix2ang(NSIDE_LO, np.arange(npix))[0])

    amp = {c: comp_lo(h5, c) for c in ("cmb", "synch", "dust", "ff", "ame")}
    theta = {k: hp.ud_grade(np.asarray(par[v][0], dtype=np.float64), NSIDE_LO) for k, v in
             [("beta_s", "synch/beta_map"), ("beta_d", "dust/beta_map"), ("T_d", "dust/T_map"),
              ("Te", "ff/Te_map"), ("nu_p", "ame/nu_p_map")]}
    spd = SpDust(TEMPLATE)
    print("BP nu_p map: median %.2f GHz, range %.2f to %.2f" % (
        np.median(theta["nu_p"]), theta["nu_p"].min(), theta["nu_p"].max()))

    bands = [b for b in BANDS if NU_MIN <= b[1] <= NU_MAX]
    nb = len(bands)
    resid = np.zeros((nb, npix))
    sigma = np.zeros((nb, npix))
    f_bp = np.zeros((nb, npix))
    high = np.abs(lat) > 30.0
    vecs = np.array(hp.pix2vec(NSIDE_LO, np.arange(npix)))
    for i, (lab, nu, fw, mfile, mscale, rfile, rscale) in enumerate(bands):
        conv = cmb_to_rj(nu)
        d = band_lo(mfile, fw, mscale) * conv                       # uK_RJ
        model = (amp["cmb"] * conv + amp["synch"] * sed_synch(nu, theta["beta_s"])
                 + amp["dust"] * sed_dust(nu, theta["beta_d"], theta["T_d"])
                 + amp["ff"] * sed_ff(nu, theta["Te"]))
        f_bp[i] = spd.sed(nu, theta["nu_p"])
        ame_model = amp["ame"] * f_bp[i]
        r = d - model
        mono, dip = hp.fit_dipole(np.where(high, r - ame_model, hp.UNSEEN), bad=hp.UNSEEN)
        resid[i] = r - mono - dip @ vecs
        rms = hp.read_map(rfile, field=0) * rscale * conv
        sigma[i] = rms if hp.get_nside(rms) == NSIDE_LO else hp.ud_grade(rms, NSIDE_LO)
        chi_ame = np.median(np.abs(resid[i] - ame_model) / sigma[i])
        print(f"band {lab:3s} {nu:5.1f} GHz: md removed ({mono:8.2f}, |d|={np.linalg.norm(dip):8.2f}) "
              f"uK_RJ; median |resid| {np.median(np.abs(resid[i])):8.2f}, median sigma "
              f"{np.median(sigma[i]):7.3f}, median |resid - BP AME|/sigma {chi_ame:6.2f}")

    mask = hp.ud_grade(hp.read_map(f"{BP}/mask_AME_nu_p_BP8_10deg.fits"), NSIDE_LO) > 0.5
    nus = np.array([b[1] for b in bands])
    grid = np.linspace(BOUNDS[0], BOUNDS[1], a.ngrid)
    F = spd.sed(nus[None, :], grid[:, None])                       # (ngrid, nb)
    logprior = -0.5 * ((grid - PRIOR_MEAN) / PRIOR_RMS) ** 2

    names = ["chisq", "ridge", "marginal", "code_marginal"]
    mean = {n: np.full(npix, np.nan) for n in names}
    sd = {n: np.full(npix, np.nan) for n in names}
    at_bound = {n: np.zeros(npix, bool) for n in names}
    pix = np.where(mask)[0]
    for p in pix:
        r, s = resid[:, p], sigma[:, p]
        w = 1.0 / s**2
        A = (F**2 * w).sum(1)
        b = (F * (r * w)).sum(1)
        chi2min = (r**2 * w).sum() - b**2 / A
        lp = {"chisq": -0.5 * ((r - amp["ame"][p] * F) ** 2 * w).sum(1),
              "ridge": -0.5 * chi2min,
              "marginal": -0.5 * chi2min - 0.5 * np.log(A),
              "code_marginal": -0.5 * chi2min + 0.5 * np.log(A)}
        for n in names:
            lg = lp[n] + logprior
            pdf = np.exp(lg - lg.max())
            pdf /= pdf.sum()
            mean[n][p] = (grid * pdf).sum()
            sd[n][p] = np.sqrt(((grid - mean[n][p]) ** 2 * pdf).sum())
            at_bound[n][p] = np.argmax(pdf) in (0, a.ngrid - 1)

    np.savez(os.path.join(a.out, "e013_results.npz"), grid=grid, mask=mask, nu_p_bp=theta["nu_p"],
             a_bp=amp["ame"], resid=resid, sigma=sigma, nus=nus, f_bp=f_bp,
             **{f"mean_{n}": mean[n] for n in names}, **{f"sd_{n}": sd[n] for n in names},
             **{f"bound_{n}": at_bound[n] for n in names})

    lines = [f"E013: {mask.sum()} pixels in the BP nu_p mask at nside {NSIDE_LO}, {nb} bands "
             f"({', '.join(f'{b[1]:g}' for b in bands)} GHz)", ""]
    lines.append(f"{'variant':14s} {'mean nu_p':>9s} {'med sd':>7s} {'shift vs ridge [ridge sd]':>34s} "
                 f"{'vs BP nu_p map [own sd]':>26s} {'MAP at bound':>12s}")
    for n in names:
        dm = (mean[n][pix] - mean["ridge"][pix]) / sd["ridge"][pix]
        db = (mean[n][pix] - theta["nu_p"][pix]) / sd[n][pix]
        lines.append(f"{n:14s} {np.nanmean(mean[n][pix]):9.2f} {np.nanmedian(sd[n][pix]):7.2f} "
                     f"med {np.nanmedian(dm):+6.2f} (p16 {np.nanpercentile(dm, 16):+6.2f}, p84 "
                     f"{np.nanpercentile(dm, 84):+6.2f}) med {np.nanmedian(db):+6.2f} "
                     f"({np.nanmean(np.abs(db) > 2) * 100:3.0f}% >2sd) {at_bound[n][pix].mean() * 100:11.1f}%")
    snr = (amp["ame"][pix] * f_bp[:, pix] / sigma[:, pix]).max(0)
    lines += ["", "By AME peak S/N (max over bands of a_BP f / sigma): median shift from ridge in "
              "ridge sd, and median ridge sd"]
    for lo_, hi_ in [(0, 3), (3, 10), (10, 30), (30, 100), (100, 1e9)]:
        sel = (snr >= lo_) & (snr < hi_)
        if sel.sum() == 0:
            continue
        cm = np.nanmedian(((mean["code_marginal"] - mean["ridge"]) / sd["ridge"])[pix][sel])
        mg = np.nanmedian(((mean["marginal"] - mean["ridge"]) / sd["ridge"])[pix][sel])
        cs = np.nanmedian(((mean["chisq"] - mean["ridge"]) / sd["ridge"])[pix][sel])
        lines.append(f"  S/N [{lo_:g}, {hi_:g}): {sel.sum():5d} px  code_marginal {cm:+.2f}  "
                     f"marginal {mg:+.2f}  chisq {cs:+.2f}  ridge sd {np.nanmedian(sd['ridge'][pix][sel]):.2f} GHz")
    # Full-sky single parameter: BP sampled one nu_p for the whole mask (the BP nu_p map is
    # constant, 25.60 GHz), so the production target was the SUM over mask pixels of the per-pixel
    # lnL. Evaluate the four summed variants on the same grid.
    W = 1.0 / sigma[:, pix] ** 2                                   # (nb, npixmask)
    A_all = np.einsum("gb,bp->gp", F**2, W)                        # (ngrid, npix)
    b_all = np.einsum("gb,bp->gp", F, resid[:, pix] * W)
    chi2min_all = (resid[:, pix] ** 2 * W).sum(0)[None, :] - b_all**2 / A_all
    chisq_all = -0.5 * (((resid[:, pix][None] - amp["ame"][pix][None, None, :] * F[:, :, None]) ** 2
                         * W[None]).sum(1))
    tot = {"chisq": chisq_all.sum(1), "ridge": -0.5 * chi2min_all.sum(1),
           "marginal": -0.5 * chi2min_all.sum(1) - 0.5 * np.log(A_all).sum(1),
           "code_marginal": -0.5 * chi2min_all.sum(1) + 0.5 * np.log(A_all).sum(1)}
    lines += ["", f"Full-sky single nu_p over the {len(pix)} mask pixels (what BP10 sampled; BP posterior "
              f"mean {theta['nu_p'][pix].mean():.2f} GHz):"]
    for n in names:
        lg = tot[n] + logprior
        pdf = np.exp(lg - lg.max()); pdf /= pdf.sum()
        mu = (grid * pdf).sum(); s_ = np.sqrt(((grid - mu) ** 2 * pdf).sum())
        lines.append(f"  {n:14s} mean {mu:6.2f} GHz  sd {s_:5.2f}  MAP {grid[np.argmax(pdf)]:6.2f}"
                     f"{'  (at bound)' if np.argmax(pdf) in (0, a.ngrid - 1) else ''}")
    # refine: the summed posterior is far narrower than the coarse grid spacing
    g0 = grid[np.argmax(tot["ridge"] + logprior)]
    fine = np.linspace(g0 - 0.5, g0 + 0.5, 4001)
    Ff = spd.sed(nus[None, :], fine[:, None])
    Af = np.einsum("gb,bp->gp", Ff**2, W)
    bf = np.einsum("gb,bp->gp", Ff, resid[:, pix] * W)
    c2f = (resid[:, pix] ** 2 * W).sum(0)[None, :] - bf**2 / Af
    csf = -0.5 * (((resid[:, pix][None] - amp["ame"][pix][None, None, :] * Ff[:, :, None]) ** 2
                   * W[None]).sum(1)).sum(1)
    totf = {"chisq": csf, "ridge": -0.5 * c2f.sum(1),
            "marginal": -0.5 * c2f.sum(1) - 0.5 * np.log(Af).sum(1),
            "code_marginal": -0.5 * c2f.sum(1) + 0.5 * np.log(Af).sum(1)}
    lpf = -0.5 * ((fine - PRIOR_MEAN) / PRIOR_RMS) ** 2
    res_full = {}
    for n in names:
        lg = totf[n] + lpf
        pdf = np.exp(lg - lg.max()); pdf /= pdf.sum()
        mu = (fine * pdf).sum(); s_ = np.sqrt(((fine - mu) ** 2 * pdf).sum())
        res_full[n] = (mu, s_)
    lines.append("  fine grid (4001 points over +-0.5 GHz around the ridge peak):")
    for n in names:
        mu, s_ = res_full[n]
        lines.append(f"    {n:14s} mean {mu:7.3f} GHz  sd {s_:6.4f}  shift vs ridge "
                     f"{(mu - res_full['ridge'][0]) / res_full['ridge'][1]:+7.2f} ridge sd  "
                     f"(vs BP 25.60: {(mu - 25.6003) / s_:+8.1f} own sd)")
    lines.append(f"  volume term range over grid: 1/2 sum ln A varies by "
                 f"{0.5 * (np.log(A_all).sum(1).max() - np.log(A_all).sum(1).min()):.1f} in log density; "
                 f"ridge chi2 range {(-2 * tot['ridge']).max() - (-2 * tot['ridge']).min():.1f}")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(a.out, "summary.txt"), "w").write(txt + "\n")


if __name__ == "__main__":
    main()
