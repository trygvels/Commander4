"""Sky regions that are what they are named, and a processing mask, for the low-resolution polarization fits.

The four regions used in E022 and E040 ("BP XV surrogate", `data.bpxv_regions`) were reconstructed from the released
BeyondPlanck index map by nearest-plateau assignment. Their "High Latitude" region contains the inner Galactic plane
at 25 to 70 degrees from the centre (55 percent of its squared amplitude lies within 10 degrees of the plane), and no
mask was applied, so Tau A sat in the "Plane" (second referee report on Paper A1, 2026-11-11). The regions here are
defined by geometry:

    0  High Latitude   |b| > 30 deg, outside the Spur
    1  Spur            the northern piece (b > 20 deg, -45 < l < 65 deg) of the pixels that the released BeyondPlanck
                       index map assigns to its Spur value (the North Polar Spur and its surroundings)
    2  Intermediate    10 <= |b| <= 30 deg, outside the Spur
    3  Plane           |b| < 10 deg

and a processing mask removes (i) the 5 percent of pixels with the largest polarized synchrotron amplitude in the
Cosmoglobe DR1 model at 30 GHz smoothed to 5 degrees, and (ii) a disc of 7.5 degrees around Tau A
(l = 184.56, b = -5.78). Svalheim et al. (2023) mask the strongest 5 percent of the signal and strong sources such
as Tau A for their index sampling; this follows that description, not their file, which we do not have.
"""
import healpy as hp
import numpy as np

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
SURROGATE = "/Users/trygvels/data/bpxv_surrogate/comps.npz"
NAMES = ["High Latitude", "Spur", "Intermediate", "Plane"]
TAU_A = (184.56, -5.78)


def dr1_polarized_amplitude(nside=16, fwhm_deg=5.0):
    """Polarized amplitude of the DR1 synchrotron model at 30 GHz (uK_RJ), smoothed, per pixel in RING order."""
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(fwhm_deg), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.hypot(qu[0], qu[1])


def named_regions(nside=16, bright_fraction=0.05, tau_a_radius_deg=7.5, mask=True):
    """Returns (region, keep, block) per (Q, U) element in NESTED order at `nside` (16 or 8).

    region: 0 to 3 as in NAMES. keep: False for masked elements. block: the nside-2 pixel of each element.
    At nside 8 the nside-16 definition is degraded: a pixel takes the majority region of its unmasked children
    and is masked if half or more of its children are."""
    ns0 = 16
    npix = hp.nside2npix(ns0)
    lon, lat = hp.pix2ang(ns0, np.arange(npix), lonlat=True)          # RING
    l180 = np.where(lon > 180.0, lon - 360.0, lon)
    region = np.full(npix, 0)
    region[np.abs(lat) <= 30.0] = 2
    region[np.abs(lat) < 10.0] = 3
    sur = np.load(SURROGATE, allow_pickle=True)["region"].astype(int)
    region[(sur == 1) & (lat > 20.0) & (l180 > -45.0) & (l180 < 65.0)] = 1
    keep = np.ones(npix, dtype=bool)
    if mask:
        P = dr1_polarized_amplitude(ns0) if bright_fraction > 0 else None
        if bright_fraction > 0:
            keep[P >= np.quantile(P, 1.0 - bright_fraction)] = False
        keep[hp.query_disc(ns0, hp.ang2vec(*TAU_A, lonlat=True), np.radians(tau_a_radius_deg))] = False
    region, keep = hp.reorder(region, r2n=True), hp.reorder(keep.astype(float), r2n=True) > 0.5
    if nside != ns0:
        k = (ns0 // nside) ** 2
        rr, kk = region.reshape(-1, k), keep.reshape(-1, k)
        reg_c = np.zeros(rr.shape[0], dtype=int)
        for i in range(rr.shape[0]):
            vals = rr[i][kk[i]] if kk[i].any() else rr[i]
            reg_c[i] = np.bincount(vals, minlength=4).argmax()
        region, keep = reg_c, kk.sum(1) > k / 2
    n1 = region.size
    block = np.arange(n1) // (nside // 2) ** 2
    return np.tile(region, 2), np.tile(keep, 2), np.tile(block, 2)
