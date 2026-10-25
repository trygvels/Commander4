"""Convert Planck PR4 (NPIPE) or PR3 (R3.01) single-frequency HFI files, or the BeyondPlanck-input NPIPE
A-split maps hosted at UiO (npipe6v20A_<freq>_map_uK.fits plus an rms file), into Commander4 maps-only
inputs.

Reads I/Q/U and the II/QQ/UU covariance columns from HFI_SkyMap_<freq>_2048_R4.00_full.fits
(NESTED, K_CMB), degrades to --nside with plain pixel averaging (the variance of the mean is the
mean variance divided by the number of merged pixels, i.e. hp.ud_grade(var, power=2)), and writes
two RING-ordered FITS files per band in uK_CMB: <out>/pr4_<freq>_IQU_n<nside>_uK.fits with columns
I_STOKES, Q_STOKES, U_STOKES and <out>/pr4_<freq>_IQU_rms_n<nside>_uK.fits with I_RMS, Q_RMS, U_RMS
(545 and 857 GHz are intensity only; their Q/U columns are zero with zero rms, which Commander4
treats as unobserved). Bad pixels (BAD_DATA sentinel or non-positive variance) get rms 0.

    python -I nextgen/scripts/prep_pr4_maps.py --indir ~/data/planck_pr4 --out ~/data/planck_pr4/n512 --nside 512
"""
import argparse
import os

import astropy.io.fits as fits
import healpy as hp
import numpy as np

FREQS = ["100", "143", "217", "353", "545", "857"]

DIPOLE_UK, DIPOLE_LB = 3362.08, (264.021, 48.253)      # Planck 2018 solar dipole, uK_CMB
HFI_MJY_PER_KCMB = {"545": 58.04, "857": 2.27}            # band-integrated, Planck 2013 IX
DR1 = "/Users/trygvels/data/cosmoglobe_DR1"
H_K = 6.62607015e-34 / 1.380649e-23 * 1e9                 # h/k in K per GHz
T_CMB = 2.7255


def _cmb_to_rj(nu):
    x = H_K * nu / T_CMB
    return x**2 * np.exp(x) / np.expm1(x)**2


def solar_dipole(nside: int) -> np.ndarray:
    vec = hp.ang2vec(*DIPOLE_LB, lonlat=True)
    return DIPOLE_UK * (hp.pix2vec(nside, np.arange(hp.nside2npix(nside))) * vec[:, None]).sum(0)


def dr1_dust_template_uk_cmb(nu: float, nside: int, fwhm_deg: float):
    """(dust, cmb) from Cosmoglobe DR1 at nu GHz, degraded to nside and smoothed to fwhm_deg, uK_CMB.
    Dust: per-pixel modified blackbody with the DR1 beta and T maps, amplitude nu_ref 545 GHz, uK_RJ.
    CMB: DR1 CMB with its dipole fitted and removed. Used only to set zero levels, dust as a
    template with a free scale (its absolute scale under the stated units is 7 to 10x below the NPIPE
    maps, see E014; the zero-level fit does not depend on it)."""
    hd = fits.open(f"{DR1}/CG_dust_IQU_n1024_v1.fits")[1]
    hc = fits.open(f"{DR1}/CG_cmb_IQU_n1024_v1.fits")[1]

    def col(h, k):
        v = np.asarray(h.data[k], float).ravel()
        return hp.reorder(v, n2r=True) if str(h.header.get("ORDERING", "RING")).upper().startswith("NEST") else v
    a, beta, T = (col(hd, k) for k in ("I_MEAN", "I_BETA_MEAN", "I_T_MEAN"))
    x0, x = H_K * 545.0 / T, H_K * nu / T
    dust = hp.smoothing(hp.ud_grade(a * (nu / 545.0)**(beta + 1.0) * np.expm1(x0) / np.expm1(x)
                                    / _cmb_to_rj(nu), nside), fwhm=np.radians(fwhm_deg))
    cmb = hp.smoothing(hp.ud_grade(col(hc, "I_MEAN"), nside), fwhm=np.radians(fwhm_deg))
    cmb = hp.remove_dipole(cmb, fitval=False, verbose=False)
    return dust, cmb


def zero_level(sky_i: np.ndarray, nu: float, method: str, fwhm_deg: float = 2.0, bcut: float = 30.0):
    """Offset to subtract from the I map (uK_CMB). "highlat": mean over |b| > bcut (as for the DR1
    WMAP/LFI maps in prep_q1_inputs.py). "dr1": offset of a least-squares fit of (map - DR1 CMB) against the DR1
    dust template over |b| > bcut after smoothing both to fwhm_deg (zero level where the dust template
    vanishes, the HI-correlation idea with the DR1 dust model as tracer); the slope is printed."""
    nside = hp.get_nside(sky_i)
    lat = 90.0 - np.degrees(hp.pix2ang(nside, np.arange(sky_i.size))[0])
    high = np.abs(lat) > bcut
    if method == "highlat":
        return float(np.mean(sky_i[high])), None
    sm = hp.smoothing(sky_i, fwhm=np.radians(fwhm_deg))
    dust, cmb = dr1_dust_template_uk_cmb(nu, nside, fwhm_deg)
    y = (sm - cmb)[high]
    A = np.vstack([dust[high], np.ones(y.size)]).T
    slope, off = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(off), float(slope)



def _col(hdu, names):
    for n in names:
        if n in hdu.columns.names:
            return np.asarray(hdu.data[n], dtype=np.float64).ravel()
    return None


def convert(path: str, nside: int, outdir: str, freq: str, tag_prefix: str = "pr4",
            zero_method: str = "") -> None:
    hdu = fits.open(path)[1]
    unit_in = str(hdu.columns["I_STOKES"].unit or "K_CMB").strip()
    # PR3 545/857 are in MJy/sr; keep that unit (Commander4 converts with pysm3 units), otherwise
    # convert K_CMB to uK_CMB.
    if (unit_in.upper().startswith("K_CMB") or unit_in.upper() == "KCMB") and freq in HFI_MJY_PER_KCMB:
        # PR4 545/857 GHz are delivered in K_CMB via the band-integrated factors (Planck 2013 IX);
        # undo that and write MJy/sr, which Commander4 converts monochromatically (E014, item 3).
        scale, unit_out = HFI_MJY_PER_KCMB[freq], "MJy/sr"
    elif unit_in.upper().startswith("K_CMB") or unit_in.upper() == "KCMB":
        scale, unit_out = 1e6, "uK_CMB"
    elif "MJY" in unit_in.upper():
        scale, unit_out = 1.0, "MJy/sr"
    else:
        raise ValueError(f"unknown unit {unit_in} in {path}")
    assert hdu.header["ORDERING"] == "NESTED", hdu.header["ORDERING"]
    bad = hdu.header.get("BAD_DATA", -1.6375e30)
    pols = ["I", "Q", "U"]
    sky, rms = [], []
    for p in pols:
        m = _col(hdu, [f"{p}_STOKES"])
        v = _col(hdu, [f"{p}{p}_COV"])
        if m is None:
            sky.append(np.zeros(hp.nside2npix(nside)))
            rms.append(np.zeros(hp.nside2npix(nside)))
            continue
        m = hp.reorder(m, n2r=True)
        v = hp.reorder(v, n2r=True)
        good = (m > bad * 0.5) & np.isfinite(m) & np.isfinite(v) & (v > 0)
        m = np.where(good, m, 0.0)
        v = np.where(good, v, 0.0)
        w = good.astype(np.float64)
        n_in = (hp.nside2npix(hp.get_nside(m)) // hp.nside2npix(nside))
        frac = hp.ud_grade(w, nside)                      # fraction of good sub-pixels
        m_lo = np.divide(hp.ud_grade(m, nside), frac, out=np.zeros_like(frac), where=frac > 0)
        # mean over the n_good good sub-pixels: var = sum(var_i) / n_good^2
        var_sum = hp.ud_grade(v, nside) * n_in
        n_good = frac * n_in
        var_lo = np.divide(var_sum, n_good**2, out=np.zeros_like(frac), where=n_good > 0)
        sky.append(m_lo * scale)
        rms.append(np.sqrt(var_lo) * scale)
    hdr_extra = []
    if zero_method:
        off, slope = zero_level(sky[0], float(freq), zero_method)
        sky[0] = sky[0] - off * (rms[0] > 0)
        hdr_extra.append(("MONOSUB", off, f"{unit_out} zero level ({zero_method}) subtracted"))
        print(f"{freq}: zero level ({zero_method}) {off:.4g} {unit_out} subtracted"
              + (f"; gain vs DR1 model: {slope:.3f}" if slope is not None else ""))
    tag = f"{freq}_IQU_n{nside:04d}"
    hp.write_map(os.path.join(outdir, f"{tag_prefix}_{tag}.fits"), np.array(sky), overwrite=True,
                 column_names=["I_STOKES", "Q_STOKES", "U_STOKES"], column_units=[unit_out] * 3,
                 coord="G", dtype=np.float64, extra_header=[("SRC", os.path.basename(path))] + hdr_extra)
    hp.write_map(os.path.join(outdir, f"{tag_prefix}_{freq}_IQU_rms_n{nside:04d}.fits"), np.array(rms),
                 overwrite=True, column_names=["I_RMS", "Q_RMS", "U_RMS"],
                 column_units=[unit_out] * 3, coord="G", dtype=np.float64)
    print(f"{freq}: unit {unit_out}, median I rms {np.median(rms[0][rms[0] > 0]):.3g}, "
          f"unobserved I pixels {np.sum(rms[0] == 0)}")


def convert_npipe_uio(indir: str, nside: int, outdir: str, freq: str, zero_method: str = "dr1") -> None:
    """BeyondPlanck inputs: npipe6v20A_<freq>_map_uK.fits (I or IQU, nside 2048) and the matching rms
    file (npipe6v20A_<freq>_rms_uK_mcscaled.fits for 353, *_hrscaled.fits for 857). Units and ordering
    are taken from the headers; degrade as in `convert` (variance of the mean over good sub-pixels)."""
    mp = os.path.join(indir, f"npipe6v20A_{freq}_map_uK.fits")
    rp = os.path.join(indir, f"npipe6v20A_{freq}_rms_uK_{'mcscaled' if freq == '353' else 'hrscaled'}.fits")
    hm, hr = fits.open(mp)[1], fits.open(rp)[1]
    unit_in = str(hm.columns[0].unit or "uK_CMB").strip()
    # The UiO files are in uK_CMB (file name "_uK", BeyondPlanck parameter files, and the pixel rms of
    # about 30 uK at 353 GHz) although their headers say K_CMB; trust the file name.
    scale, unit_out = 1.0, "uK_CMB"
    nested = hm.header.get("ORDERING", "RING").upper().startswith("NEST")
    nested_r = hr.header.get("ORDERING", "RING").upper().startswith("NEST")
    ncol = len(hm.columns)
    npix_out = hp.nside2npix(nside)
    sky, rms = [], []
    for i in range(3):
        if i >= ncol:
            sky.append(np.zeros(npix_out)); rms.append(np.zeros(npix_out)); continue
        m = np.asarray(hm.data[hm.columns[i].name], dtype=np.float64).ravel()
        r = np.asarray(hr.data[hr.columns[i].name], dtype=np.float64).ravel()
        if nested: m = hp.reorder(m, n2r=True)
        if nested_r: r = hp.reorder(r, n2r=True)
        v = r**2
        good = np.isfinite(m) & np.isfinite(v) & (v > 0) & (m > -1e20)
        m = np.where(good, m, 0.0); v = np.where(good, v, 0.0)
        w = good.astype(np.float64)
        n_in = hp.nside2npix(hp.get_nside(m)) // npix_out
        frac = hp.ud_grade(w, nside)
        m_lo = np.divide(hp.ud_grade(m, nside), frac, out=np.zeros_like(frac), where=frac > 0)
        var_sum = hp.ud_grade(v, nside) * n_in
        n_good = frac * n_in
        var_lo = np.divide(var_sum, n_good**2, out=np.zeros_like(frac), where=n_good > 0)
        sky.append(m_lo * scale); rms.append(np.sqrt(var_lo) * scale)
    # NPIPE frequency maps keep the solar dipole; remove it, then set the zero level.
    sky[0] = sky[0] - solar_dipole(nside) * (rms[0] > 0)
    off, slope = zero_level(sky[0], float(freq), zero_method)
    sky[0] = sky[0] - off * (rms[0] > 0)
    print(f"{freq}: dipole removed; zero level ({zero_method}) {off:.4g} uK_CMB subtracted"
          + (f"; gain vs DR1 model at |b|>30: {slope:.3f}" if slope is not None else ""))
    # BeyondPlanck converted the 545 and 857 GHz MJy/sr maps to K_CMB with the band-integrated HFI
    # factors (Planck 2013 IX: 58.04 and 2.27 MJy/sr per K_CMB). Commander4 converts a K_CMB band
    # back to brightness with the monochromatic dB/dT at the band centre (1.44 MJy/sr per K_CMB at
    # 857 GHz), which would make 857 GHz 1.6x too faint (E014). Undo BP's conversion and write
    # MJy/sr, which Commander4 converts to uK_RJ monochromatically as for the PR3 files.
    if freq in HFI_MJY_PER_KCMB:
        f = HFI_MJY_PER_KCMB[freq] * 1e-6
        sky = [m * f for m in sky]; rms = [r * f for r in rms]; off *= f
        unit_out = "MJy/sr"
    tag = f"{freq}_IQU_n{nside:04d}"
    hp.write_map(os.path.join(outdir, f"npipe_{tag}.fits"), np.array(sky), overwrite=True,
                 column_names=["I_STOKES", "Q_STOKES", "U_STOKES"], column_units=[unit_out] * 3,
                 coord="G", dtype=np.float64,
                 extra_header=[("SRC", os.path.basename(mp)), ("DIPSUB", DIPOLE_UK, "uK_CMB solar dipole removed"),
                               ("MONOSUB", off, f"{unit_out} zero level ({zero_method}) subtracted")])
    hp.write_map(os.path.join(outdir, f"npipe_{freq}_IQU_rms_n{nside:04d}.fits"), np.array(rms),
                 overwrite=True, column_names=["I_RMS", "Q_RMS", "U_RMS"], column_units=[unit_out] * 3,
                 coord="G", dtype=np.float64)
    print(f"{freq} (npipe A, UiO): unit in {unit_in!r} -> {unit_out}, ncol {ncol}, nested {nested}, "
          f"median I rms {np.median(rms[0][rms[0] > 0]):.3g}, unobserved I pixels {np.sum(rms[0] == 0)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nside", type=int, default=512)
    ap.add_argument("--freqs", nargs="*", default=FREQS)
    ap.add_argument("--release", default="pr4", choices=["pr4", "pr3", "npipe_uio"])
    ap.add_argument("--zero-level", default=None, choices=["", "highlat", "dr1"],
                    help="zero-level method; default dr1 for npipe_uio, none for pr3/pr4")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for f in a.freqs:
        if a.release == "npipe_uio":
            convert_npipe_uio(a.indir, a.nside, a.out, f, a.zero_level or "dr1")
            continue
        if a.release == "pr4":
            name = f"HFI_SkyMap_{f}_2048_R4.00_full.fits"
        else:
            name = f"HFI_SkyMap_{'353-psb' if f == '353' else f}_2048_R3.01_full.fits"
        path = os.path.join(a.indir, name)
        if not os.path.exists(path):
            print(f"missing {path}")
            continue
        convert(path, a.nside, a.out, f, a.release, a.zero_level or "")


if __name__ == "__main__":
    main()
