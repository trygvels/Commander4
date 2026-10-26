"""Pass-4 Q1 inputs (E014): one zero-level convention for all 20 bands, set by the Cosmoglobe DR1
sky model.

For every band b the offset is

    offset_b = mean over |b| > 30 deg of (map_b - model_b),

where model_b is the DR1 intensity sky model evaluated at the band centre with a delta bandpass,
in the band's unit: CMB + synchrotron + free-free + AME + thermal dust from the released DR1
component maps (CG_*_n1024_v1.fits). The output map is map_b minus the Planck 2018 solar dipole
(where the input still carries it) minus offset_b, so that all 20 maps are consistent with one
sky model up to the residual dipole and the point sources, which are reported.

DR1 SEDs (Commander3 as run for DR1, CG_param_c0001_v1.txt and param_CG_v09.txt; the released
maps hold the parameter values):
  synchrotron  power law, nu_ref 0.408 GHz, beta = I_BETA_MEAN (single value -3.1536);
  free-free    (nu_ref/nu)^2 g(nu)/g(nu_ref) with the BP1 Gaunt factor, T_e = I_TE_MEAN (7000 K),
               nu_ref 40 GHz;
  AME          Commander3 type "exponential" (components/ame/ame_CG_exp.defaults):
               exp(beta (nu/nu_ref - 1)), nu_ref 22 GHz, zero below 10 GHz, beta is the value
               stored in the column named I_NU_P_MEAN (-3.5677, constant; the column name is a
               mislabel of the Commander3 output, DR1 did not use a spdust template);
  thermal dust MBB (nu/nu_ref)^(beta+1) expm1(h nu_ref/kT)/expm1(h nu/kT), nu_ref 70 GHz (the
               FITS header says 545 GHz; finding F13), beta = I_BETA_MEAN (1.4999), T = I_T_MEAN;
  CMB          I_MEAN in uK_CMB, which includes the fitted solar dipole (3363 uK toward
               (264.3, 48.3)).
All amplitudes are uK_RJ at nu_ref; the sum is converted to the band unit with Commander4's
`rj_to_band_unit_factor` (pysm3 CMB equivalencies at band centre), exactly as Commander4 converts a
band. The SED formulas are the same as the Commander4 classes in sky/diffuse_components.py
(ThermalDust, Synchrotron, FreeFree); they are written out here because those classes need a
full component configuration and Commander4 has no exponential AME. Not modelled: the DR1 radio
point sources and the relativistic CMB quadrupole (both small in a |b| > 30 mean).

Delta bandpasses at the band centres of nextgen/params/q1/param_q1_n64_ridge_pr4tuned.yml (printed).
Bandpass integration and colour corrections, which DR1 applied, are not: the foreground means
carry a few percent error at the low-frequency bands.

    source .venv/bin/activate
    python -I nextgen/scripts/dr1_model_zero_levels.py --out ~/data/q1_inputs_pass4
"""
import argparse
import os
import shutil
import sys

import astropy.io.fits as fits
import healpy as hp
import numpy as np

from commander4.units import rj_to_band_unit_factor

DR1 = "/Users/trygvels/data/cosmoglobe_DR1"
OLD_DR1 = "/Users/trygvels/data/q1_inputs"               # pass-2/3 WMAP/LFI inputs
OLD_PR4 = "/Users/trygvels/data/planck_pr4/n512_pass3"  # pass-3 PR4 HFI inputs
BP = "/Users/trygvels/data/beyondplanck"                 # Haslam map used by the q1 parameter files
HASLAM_REMAZEILLES = "/Users/trygvels/data/haslam/haslam408_dsds_Remazeilles2014.fits"
CHISQ_MASK = "/Users/trygvels/data/q1_inputs/mask_common_dx12_n0064_chisq.fits"

DIPOLE_UK, DIPOLE_LB = 3362.08, (264.021, 48.253)       # Planck 2018 solar dipole, uK_CMB
HFI_MJY_PER_KCMB = {"545": 58.04, "857": 2.27}           # band-integrated (Planck 2013 IX), as in
                                                         # prep_pr4_maps.py (MONOSUB there is uK_CMB)
H_K = 6.62607015e-34 / 1.380649e-23 * 1e9                # h/k in K per GHz
NSIDE = 512
BCUT = 30.0

# name, input file (native DR1 or pass-3), output file, band-centre GHz, file unit
DR1_BANDS = [
    ("WMAPK", "CG_023-WMAP_K_IQU_n0512_v1.fits", 22.8, "mK_CMB"),
    ("WMAPKa", "CG_030-WMAP_Ka_IQU_n0512_v1.fits", 33.0, "mK_CMB"),
    ("WMAPQ1", "CG_040-WMAP_Q1_IQU_n0512_v1.fits", 40.6, "mK_CMB"),
    ("WMAPQ2", "CG_040-WMAP_Q2_IQU_n0512_v1.fits", 40.6, "mK_CMB"),
    ("WMAPV1", "CG_060-WMAP_V1_IQU_n0512_v1.fits", 60.8, "mK_CMB"),
    ("WMAPV2", "CG_060-WMAP_V2_IQU_n0512_v1.fits", 60.8, "mK_CMB"),
    ("WMAPW1", "CG_090-WMAP_W1_IQU_n0512_v1.fits", 93.5, "mK_CMB"),
    ("WMAPW2", "CG_090-WMAP_W2_IQU_n0512_v1.fits", 93.5, "mK_CMB"),
    ("WMAPW3", "CG_090-WMAP_W3_IQU_n0512_v1.fits", 93.5, "mK_CMB"),
    ("WMAPW4", "CG_090-WMAP_W4_IQU_n0512_v1.fits", 93.5, "mK_CMB"),
    ("Planck30GHz", "CG_030_IQU_n0512_v1.fits", 28.4, "uK_CMB"),
    ("Planck44GHz", "CG_044_IQU_n0512_v1.fits", 44.1, "uK_CMB"),
    ("Planck70GHz", "CG_070_IQU_n1024_v1.fits", 70.1, "uK_CMB"),
]
HFI_BANDS = [("Planck100GHz", "100", 100.0), ("Planck143GHz", "143", 143.0),
             ("Planck217GHz", "217", 217.0), ("Planck353GHz", "353", 353.0),
             ("Planck545GHz", "545", 545.0), ("Planck857GHz", "857", 857.0)]
HASLAM_MAP = "haslam_n0512_60arcmin_uK.fits"
HASLAM_RMS = "haslam_rms_n0512_60arcmin_uK_1pct_0.8K.fits"


def col(hdu, name):
    v = np.asarray(hdu.data[name], dtype=np.float64).ravel()
    if str(hdu.header.get("ORDERING", "RING")).upper().startswith("NEST"):
        v = hp.reorder(v, n2r=True)
    return v


def solar_dipole(nside):
    vec = hp.ang2vec(*DIPOLE_LB, lonlat=True)
    return DIPOLE_UK * (np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))) * vec[:, None]).sum(0)


def highlat(nside, bcut=BCUT):
    lat = 90.0 - np.degrees(hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))[0])
    return np.abs(lat) > bcut


def gaunt(nu, T_e):
    """BP1 eq. 18 Gaunt factor, as in Commander3 comm_freefree_comp_mod and Commander4 FreeFree."""
    return np.log(np.exp(5.960 - np.sqrt(3.0) / np.pi * np.log(nu * (T_e / 1e4)**-1.5)) + np.e)


class DR1Model:
    """DR1 intensity components at nside 1024 (uK_RJ amplitudes plus per-pixel SED parameters)."""

    def __init__(self):
        hc = fits.open(f"{DR1}/CG_cmb_IQU_n1024_v1.fits")[1]
        hs = fits.open(f"{DR1}/CG_synch_IQU_n1024_v1.fits")[1]
        hf = fits.open(f"{DR1}/CG_freefree_I_n1024_v1.fits")[1]
        ha = fits.open(f"{DR1}/CG_ame_I_n1024_v1.fits")[1]
        hd = fits.open(f"{DR1}/CG_dust_IQU_n1024_v1.fits")[1]
        self.cmb = col(hc, "I_MEAN")                           # uK_CMB, with the DR1 solar dipole
        self.s_amp, self.s_beta = col(hs, "I_MEAN"), col(hs, "I_BETA_MEAN")
        self.s_ref = float(hs.header["NU_REF_T"].split()[0])   # 0.408 GHz
        self.f_amp, self.f_Te = col(hf, "I_MEAN"), col(hf, "I_TE_MEAN")
        self.f_ref = float(hf.header["NU_REF_T"].split()[0])   # 40 GHz
        self.a_amp, self.a_beta = col(ha, "I_MEAN"), col(ha, "I_NU_P_MEAN")  # exponential index
        self.a_ref = float(ha.header["NU_REF_T"].split()[0])   # 22 GHz
        self.d_amp, self.d_beta, self.d_T = col(hd, "I_MEAN"), col(hd, "I_BETA_MEAN"), col(hd, "I_T_MEAN")
        self.d_ref = 70.0  # F13: COMP_NU_REF_T03 = "70 857" in the DR1 parameter file, header says 545
        self.nside = hp.get_nside(self.cmb)
        self.cmb_nodip = self.cmb - solar_dipole(self.nside)
        print(f"DR1 model loaded at nside {self.nside}: synch nu_ref {self.s_ref} GHz, beta "
              f"{np.unique(self.s_beta)}; ff nu_ref {self.f_ref} GHz, T_e {np.unique(self.f_Te)}; "
              f"AME exponential nu_ref {self.a_ref} GHz, beta {np.unique(self.a_beta)}; dust nu_ref "
              f"{self.d_ref} GHz (header {hd.header['NU_REF_T']}), beta {np.unique(self.d_beta)}, "
              f"T median {np.median(self.d_T):.2f} K")
        mono, dip = hp.fit_dipole(self.cmb, gal_cut=BCUT)
        print(f"DR1 CMB map: |b|>{BCUT:g} dipole {np.linalg.norm(dip):.1f} uK toward "
              f"(l,b)=({hp.vec2ang(dip, lonlat=True)[0][0]:.2f},{hp.vec2ang(dip, lonlat=True)[1][0]:.2f}), "
              f"monopole {mono:.2f} uK")

    def foregrounds_rj(self, nu):
        """Sum of the four foregrounds at nu GHz in uK_RJ, nside 1024."""
        synch = self.s_amp * (nu / self.s_ref)**self.s_beta
        ff = self.f_amp * (self.f_ref / nu)**2 * gaunt(nu, self.f_Te) / gaunt(self.f_ref, self.f_Te)
        ame = self.a_amp * np.exp(self.a_beta * (nu / self.a_ref - 1.0)) if nu >= 10.0 else 0.0
        x, x0 = H_K * nu / self.d_T, H_K * self.d_ref / self.d_T
        dust = self.d_amp * (nu / self.d_ref)**(self.d_beta + 1.0) * np.expm1(x0) / np.expm1(x)
        return synch + ff + ame + dust

    def evaluate(self, nu, unit, nside, with_dipole):
        """(model, foregrounds, cmb) at nu GHz in `unit` (Commander4 conversion at band centre),
        degraded to nside by pixel averaging."""
        f = rj_to_band_unit_factor(float(nu), unit)
        f_cmb = f / rj_to_band_unit_factor(float(nu), "uK_CMB")   # uK_CMB -> unit
        fg = hp.ud_grade(self.foregrounds_rj(nu) * f, nside)
        cmb = hp.ud_grade((self.cmb if with_dipole else self.cmb_nodip) * f_cmb, nside)
        return fg + cmb, fg, cmb


def fit_dipole_hl(m, good):
    """Monopole and dipole of m over the good pixels (least squares); returns (mono, amp, l, b)."""
    nside = hp.get_nside(m)
    vec = np.array(hp.pix2vec(nside, np.arange(m.size)))
    A = np.vstack([np.ones(good.sum()), vec[:, good]]).T
    c = np.linalg.lstsq(A, m[good], rcond=None)[0]
    amp = np.linalg.norm(c[1:])
    l, b = hp.vec2ang(c[1:], lonlat=True)
    return c[0], amp, float(l[0]), float(b[0])


def verify(model, outdir):
    """DR1 model against the maps at 30, 70 GHz (DR1 native) and 353 GHz (PR4 pass 3, offset
    restored): slope, intercept and correlation at |b| > 30 after 2 deg smoothing, for the full
    map (dipole removed from both) and for map minus model CMB against the foreground model."""
    print("\n=== Verification of the model evaluation (2 deg smoothing, |b| > 30) ===")
    good = highlat(NSIDE)
    cases = [("30 GHz DR1", f"{DR1}/CG_030_IQU_n0512_v1.fits", "I_MEAN", 28.4, 1.0, True),
             ("70 GHz DR1", f"{DR1}/CG_070_IQU_n1024_v1.fits", "I_MEAN", 70.1, 1.0, True),
             ("353 GHz PR4", f"{OLD_PR4}/pr4_353_IQU_n0512.fits", "I_STOKES", 353.0, 1.0, False)]
    lines = []
    for label, fn, c, nu, scale, dip in cases:
        h = fits.open(fn)[1]
        m = hp.ud_grade(col(h, c) * scale, NSIDE)
        if not dip:
            m = m + float(h.header.get("MONOSUB", 0.0))      # back to the native zero level
        mod, fg, cmb = model.evaluate(nu, "uK_CMB", NSIDE, with_dipole=dip)
        if dip:
            m, mod, cmb = m - solar_dipole(NSIDE), mod - solar_dipole(NSIDE), cmb - solar_dipole(NSIDE)
        sm = [hp.smoothing(x, fwhm=np.radians(2.0))
              for x in (m, mod, fg, cmb)]
        for what, y, x in (("map vs model", sm[0], sm[1]), ("map-cmb vs fg", sm[0] - sm[3], sm[2])):
            A = np.vstack([x[good], np.ones(good.sum())]).T
            slope, icpt = np.linalg.lstsq(A, y[good], rcond=None)[0]
            corr = np.corrcoef(x[good], y[good])[0, 1]
            rms_res = np.std(y[good] - slope * x[good] - icpt)
            line = (f"{label:12s} {what:14s} slope {slope:7.4f} intercept {icpt:9.2f} uK corr {corr:.5f}"
                    f" resid rms {rms_res:8.2f} uK (model rms {np.std(x[good]):.1f})")
            print(line); lines.append(line)
    with open(os.path.join(outdir, "verification.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def band_offsets(m, good, model_nodip, label):
    """Offset (mean of map - model over good pixels) and the dipole of the residual after it."""
    r = m - model_nodip
    off = float(np.mean(r[good]))
    mono, amp, l, b = fit_dipole_hl(r - off, good)
    return off, (amp, l, b, mono)


def process_all(model, outdir, mask512, haslam_dipole=False):
    hl = highlat(NSIDE)
    rows = []
    print(f"\nBand centres (GHz) used with delta bandpasses: 0.408, " +
          ", ".join(f"{nu:g}" for _, _, nu, _ in DR1_BANDS) + ", " +
          ", ".join(f"{nu:g}" for _, _, nu in HFI_BANDS))

    def record(name, nu, unit, old_off, off, off_mask, dipres, model_mean, dipsub, note=""):
        rows.append((name, nu, unit, old_off, off, off_mask, dipres, model_mean, dipsub, note))
        amp, l, b, mono = dipres
        print(f"{name:13s} {nu:7.3f} GHz  old {old_off:12.5g}  new {off:12.5g}  diff {off - old_off:10.4g}"
              f"  [chisq-mask {off_mask:10.4g}]  resid dipole {amp:8.4g} toward ({l:6.1f},{b:5.1f})"
              f"  model mean {model_mean:10.4g} {unit}  {note}")

    # --- Haslam 408 MHz (BP file used by the q1 parameter files; uK, no solar dipole) ---
    hh = fits.open(f"{BP}/{HASLAM_MAP}")
    has = col(hh[1], "TEMPERATURE")
    hr = col(fits.open(f"{BP}/{HASLAM_RMS}")[1], "TEMPERATURE")
    mod, fg, cmb = model.evaluate(0.408, "uK_CMB", NSIDE, with_dipole=False)
    off, dipres = band_offsets(has, hl, mod, "Haslam")
    off_mask = float(np.mean((has - mod)[mask512]))
    rem = col(fits.open(HASLAM_REMAZEILLES)[1], "TEMPERATURE") * 1e6
    note = f"(BP file minus Remazeilles2014 at |b|>30: {np.mean((has - rem)[hl]):.4g} uK)"
    record("Haslam408MHz", 0.408, "uK_CMB", 0.0, off, off_mask, dipres, float(np.mean(mod[hl])), 0.0, note)
    has_out = has - off
    hdip = [("HDIPSUB", 0.0, "uK fitted non-solar dipole removed (amp)")]
    if haslam_dipole:
        # The DR1 md component holds a 3.4 K dipole toward (l, b) = (12, -14) for Haslam (chain
        # sample 000001: 3.2e6, 0.7e6, -0.8e6 uK), which this residual reproduces; opt-in removal.
        amp, l, b = dipres[0], dipres[1], dipres[2]
        vec = hp.ang2vec(l, b, lonlat=True)
        has_out = has_out - amp * (np.array(hp.pix2vec(NSIDE, np.arange(has.size))) * vec[:, None]).sum(0)
        hdip = [("HDIPSUB", amp, "uK fitted non-solar dipole removed (amp)"),
                ("HDIPL", l, "deg"), ("HDIPB", b, "deg")]
    out = np.array([has_out, np.zeros_like(has), np.zeros_like(has)])
    hp.write_map(os.path.join(outdir, HASLAM_MAP), out, overwrite=True, coord="G", dtype=np.float64,
                 column_names=["TEMPERATURE", "Q-POLARISATION", "U-POLARISATION"],
                 column_units=["uK_CMB"] * 3,
                 extra_header=[("DIPSUB", 0.0, "uK_CMB solar dipole subtracted (none)"),
                               ("MONOSUB", off, "uK_CMB zero level subtracted"),
                               ("ZLMETHOD", "dr1model", "mean(map - DR1 model) at |b|>30"),
                               ("SRC", f"{BP}/{HASLAM_MAP}")] + hdip)
    shutil.copyfile(f"{BP}/{HASLAM_RMS}", os.path.join(outdir, HASLAM_RMS))

    # --- 13 DR1 WMAP and LFI bands, native files ---
    for name, fn, nu, unit in DR1_BANDS:
        h = fits.open(f"{DR1}/{fn}")[1]
        sky = np.array([col(h, k) for k in ("I_MEAN", "Q_MEAN", "U_MEAN")])
        nside = hp.get_nside(sky[0])
        good = highlat(nside)
        to_uk = 1000.0 if unit == "mK_CMB" else 1.0
        I_uk = sky[0] * to_uk
        mod_dip, _, _ = model.evaluate(nu, "uK_CMB", nside, with_dipole=True)
        dip_vec = solar_dipole(nside)
        mod_nodip = mod_dip - dip_vec
        # does the native map carry the solar dipole? compare the residual dipole both ways
        _, a_with, _, _ = fit_dipole_hl(I_uk - mod_dip, good)
        _, a_without, _, _ = fit_dipole_hl(I_uk - mod_nodip, good)
        carries = a_with < a_without
        if carries:
            I_uk = I_uk - dip_vec
        off, dipres = band_offsets(I_uk, good, mod_nodip, name)
        msk = mask512 if nside == NSIDE else hp.ud_grade(mask512.astype(float), nside) > 0.5
        off_mask = float(np.mean((I_uk - mod_nodip)[msk]))
        base = fn.replace("_v1.fits", "_nodip_v1.fits")
        old = fits.open(f"{OLD_DR1}/{base}")[1].header
        old_off = float(old["MONOSUB"]) * to_uk
        note = (f"(dipole {'carried' if carries else 'absent'}: resid {a_with:.1f} with / "
                f"{a_without:.1f} uK without model dipole)")
        model_mean = float(np.mean(mod_nodip[good]))
        record(name, nu, "uK_CMB", old_off, off, off_mask, dipres, model_mean,
               DIPOLE_UK if carries else 0.0, note)
        sky[0] = (I_uk - off) / to_uk
        hp.write_map(os.path.join(outdir, base), sky, overwrite=True, coord="G", dtype=np.float64,
                     column_names=["I_MEAN", "Q_MEAN", "U_MEAN"], column_units=[unit] * 3,
                     extra_header=[("DIPSUB", DIPOLE_UK if carries else 0.0, "uK CMB solar dipole subtracted"),
                                   ("MONOSUB", off / to_uk, f"{unit} zero level subtracted"),
                                   ("ZLMETHOD", "dr1model", "mean(map - DR1 model) at |b|>30"),
                                   ("SRC", fn)])
        rms_name = base.replace("_IQU_", "_IQU_rms_")
        shutil.copyfile(f"{OLD_DR1}/{rms_name}", os.path.join(outdir, rms_name))

    # --- six PR4 HFI bands from pass 3 (dipole removed, "dr1" offset in MONOSUB, uK_CMB) ---
    for name, f, nu in HFI_BANDS:
        h = fits.open(f"{OLD_PR4}/pr4_{f}_IQU_n0512.fits")[1]
        hr_ = fits.open(f"{OLD_PR4}/pr4_{f}_IQU_rms_n0512.fits")[1]
        sky = np.array([col(h, k) for k in ("I_STOKES", "Q_STOKES", "U_STOKES")])
        rms_i = col(hr_, "I_RMS")
        unit = str(h.columns["I_STOKES"].unit)
        obs = rms_i > 0
        old_off_uk = float(h.header["MONOSUB"])            # uK_CMB (prep_pr4_maps convention)
        if unit == "MJy/sr":
            old_off = old_off_uk * HFI_MJY_PER_KCMB[f] * 1e-6   # what was subtracted, in MJy/sr
        else:
            old_off = old_off_uk
        I = sky[0] + old_off * obs                           # back to the dipole-free native level
        mod, _, _ = model.evaluate(nu, unit, NSIDE, with_dipole=False)
        good = hl & obs
        off, dipres = band_offsets(I, good, mod, name)
        off_mask = float(np.mean((I - mod)[mask512 & obs]))
        model_mean = float(np.mean(mod[good]))
        note = f"(unobserved I pixels {np.sum(~obs)})" if np.sum(~obs) else ""
        record(name, nu, unit, old_off, off, off_mask, dipres, model_mean, DIPOLE_UK, note)
        sky[0] = (I - off) * obs
        hp.write_map(os.path.join(outdir, f"pr4_{f}_IQU_n0512.fits"), sky, overwrite=True, coord="G",
                     dtype=np.float64, column_names=["I_STOKES", "Q_STOKES", "U_STOKES"],
                     column_units=[unit] * 3,
                     extra_header=[("SRC", h.header.get("SRC", "")),
                                   ("DIPSUB", DIPOLE_UK, "uK_CMB solar dipole removed"),
                                   ("MONOSUB", off, f"{unit} zero level subtracted"),
                                   ("MONOSUBK", off / (HFI_MJY_PER_KCMB[f] * 1e-6) if unit == "MJy/sr" else off,
                                    "same in uK_CMB (band-integrated factor for 545/857)"),
                                   ("ZLMETHOD", "dr1model", "mean(map - DR1 model) at |b|>30")])
        shutil.copyfile(f"{OLD_PR4}/pr4_{f}_IQU_rms_n0512.fits",
                        os.path.join(outdir, f"pr4_{f}_IQU_rms_n0512.fits"))
    return rows


def write_table(rows, outdir):
    hdr = (f"{'band':13s} {'GHz':>8s} {'unit':7s} {'old_offset':>12s} {'new_offset':>12s} {'diff':>10s} "
           f"{'off_chisqmask':>13s} {'resid_dip':>10s} {'dip_l':>6s} {'dip_b':>6s} {'model_mean':>11s} "
           f"{'DIPSUB':>8s}  note")
    lines = ["# Pass-4 zero levels: offset = mean over |b| > 30 deg of (map - DR1 sky model), map units",
             "# old_offset: MONOSUB of the pass-2/3 input in the same unit (WMAP mK converted to uK_CMB;",
             "# 545/857 uK_CMB converted to MJy/sr with 58.04 and 2.27 MJy/sr per K_CMB)",
             "# off_chisqmask: the same mean over the dx12 chisq mask (fsky 0.78) instead of |b| > 30",
             "# resid_dip: dipole fitted at |b| > 30 to (map - model - new_offset), with direction (l, b)",
             "# model_mean: mean of the DR1 model (CMB without dipole plus foregrounds) at |b| > 30",
             hdr]
    for name, nu, unit, old, off, offm, (amp, l, b, mono), mm, dipsub, note in rows:
        lines.append(f"{name:13s} {nu:8.3f} {unit:7s} {old:12.5g} {off:12.5g} {off - old:10.4g} "
                     f"{offm:13.5g} {amp:10.4g} {l:6.1f} {b:6.1f} {mm:11.5g} {dipsub:8.2f}  {note}")
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(outdir, "zero_levels.txt"), "w") as fh:
        fh.write(txt)
    print("\n" + txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/Users/trygvels/data/q1_inputs_pass4")
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--remove-haslam-dipole", action="store_true",
                    help="also remove the fitted non-solar dipole from Haslam (DR1 absorbed a 3.4 K "
                         "dipole in its Haslam md component); default leaves it in the map")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    model = DR1Model()
    hm = fits.open(CHISQ_MASK)[1]
    mask512 = hp.ud_grade(col(hm, hm.columns.names[0]), NSIDE) > 0.5
    print(f"chisq mask fsky at nside {NSIDE}: {mask512.mean():.4f}")
    if not a.no_verify:
        verify(model, a.out)
    rows = process_all(model, a.out, mask512, a.remove_haslam_dipole)
    write_table(rows, a.out)
    for extra in ("mask_b20_n0064.fits", "mask_common_dx12_n0064_chisq.fits"):
        if os.path.exists(f"{OLD_DR1}/{extra}"):
            shutil.copyfile(f"{OLD_DR1}/{extra}", os.path.join(a.out, extra))


if __name__ == "__main__":
    main()
