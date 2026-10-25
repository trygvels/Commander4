"""Prepare the Cosmoglobe DR1 and Haslam maps for the Q1 maps-only runs (E014).

For each DR1 frequency map: read I/Q/U means and rms, subtract the Planck 2018 CMB solar dipole
(3362.08 uK towards l = 264.021 deg, b = 48.253 deg; Planck 2018 I) from I, then subtract a
monopole fitted to I at |b| > 30 deg (an approximation of Commander3's per-band monopole sampling;
Commander4 has no per-band monopole component). Q/U are left unchanged. Units are preserved and
written to the header (WMAP bands: mK_CMB; LFI bands: uK_CMB). Haslam is copied with the same
monopole treatment disabled (its zero level is part of the BeyondPlanck rms model).

    python -I nextgen/scripts/prep_q1_inputs.py --out ~/data/q1_inputs
"""
import argparse
import glob
import os

import healpy as hp
import numpy as np

DR1 = "/Users/trygvels/data/cosmoglobe_DR1"
DIPOLE_UK = 3362.08
DIPOLE_LB = (264.021, 48.253)


def dipole_map(nside: int, amp_uK: float) -> np.ndarray:
    vec = hp.ang2vec(*DIPOLE_LB, lonlat=True)
    pix_vec = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside))))
    return amp_uK * (vec @ pix_vec)


def process(fn: str, outdir: str, unit_scale_to_uK: float) -> None:
    cols = hp.read_map(fn, field=(0, 1, 2, 3, 4, 5))
    sky, rms = np.array(cols[:3], dtype=np.float64), np.array(cols[3:], dtype=np.float64)
    nside = hp.get_nside(sky[0])
    sky[0] -= dipole_map(nside, DIPOLE_UK) / unit_scale_to_uK
    lat = np.degrees(hp.pix2ang(nside, np.arange(sky.shape[1]))[0])
    high = np.abs(lat - 90.0) > 30.0
    mono = np.median(sky[0][high])
    sky[0] -= mono
    base = os.path.basename(fn).replace("_v1.fits", "_nodip_v1.fits")
    unit = "uK_CMB" if unit_scale_to_uK == 1.0 else "mK_CMB"
    hp.write_map(os.path.join(outdir, base), sky, overwrite=True, coord="G", dtype=np.float64,
                 column_names=["I_MEAN", "Q_MEAN", "U_MEAN"], column_units=[unit] * 3,
                 extra_header=[("DIPSUB", DIPOLE_UK, "uK CMB solar dipole subtracted"),
                               ("MONOSUB", mono, f"{unit} high-latitude monopole subtracted")])
    hp.write_map(os.path.join(outdir, base.replace("_IQU_", "_IQU_rms_")), rms, overwrite=True,
                 coord="G", dtype=np.float64, column_names=["I_RMS", "Q_RMS", "U_RMS"],
                 column_units=[unit] * 3)
    print(f"{base}: unit {unit}, monopole removed {mono:.4g} {unit}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="/Users/trygvels/data/q1_inputs")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    for fn in sorted(glob.glob(f"{DR1}/CG_0*_IQU_n*_v1.fits")):
        is_wmap = "WMAP" in fn
        process(fn, a.out, 1000.0 if is_wmap else 1.0)


if __name__ == "__main__":
    main()
