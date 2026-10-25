"""Convert Planck PR4 (NPIPE) or PR3 (R3.01) single-frequency HFI files into Commander4 maps-only inputs.

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


def _col(hdu, names):
    for n in names:
        if n in hdu.columns.names:
            return np.asarray(hdu.data[n], dtype=np.float64).ravel()
    return None


def convert(path: str, nside: int, outdir: str, freq: str, tag_prefix: str = "pr4") -> None:
    hdu = fits.open(path)[1]
    unit_in = str(hdu.columns["I_STOKES"].unit or "K_CMB").strip()
    # PR3 545/857 are in MJy/sr; keep that unit (Commander4 converts with pysm3 units), otherwise
    # convert K_CMB to uK_CMB.
    if unit_in.upper().startswith("K_CMB") or unit_in.upper() == "KCMB":
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
    tag = f"{freq}_IQU_n{nside:04d}"
    hp.write_map(os.path.join(outdir, f"{tag_prefix}_{tag}.fits"), np.array(sky), overwrite=True,
                 column_names=["I_STOKES", "Q_STOKES", "U_STOKES"], column_units=[unit_out] * 3,
                 coord="G", dtype=np.float64, extra_header=[("SRC", os.path.basename(path))])
    hp.write_map(os.path.join(outdir, f"{tag_prefix}_{freq}_IQU_rms_n{nside:04d}.fits"), np.array(rms),
                 overwrite=True, column_names=["I_RMS", "Q_RMS", "U_RMS"],
                 column_units=[unit_out] * 3, coord="G", dtype=np.float64)
    print(f"{freq}: unit {unit_out}, median I rms {np.median(rms[0][rms[0] > 0]):.3g}, "
          f"unobserved I pixels {np.sum(rms[0] == 0)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--nside", type=int, default=512)
    ap.add_argument("--freqs", nargs="*", default=FREQS)
    ap.add_argument("--release", default="pr4", choices=["pr4", "pr3"])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for f in a.freqs:
        if a.release == "pr4":
            name = f"HFI_SkyMap_{f}_2048_R4.00_full.fits"
        else:
            name = f"HFI_SkyMap_{'353-psb' if f == '353' else f}_2048_R3.01_full.fits"
        path = os.path.join(a.indir, name)
        if not os.path.exists(path):
            print(f"missing {path}")
            continue
        convert(path, a.nside, a.out, f, a.release)


if __name__ == "__main__":
    main()
