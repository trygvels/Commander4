"""Referee check of lowres/regions.py against the description in Sect. 6.1 of Paper A1 (draft 4). Small: reads the
DR1 map once. Counts per region, mask composition, the Spur definition, nside-8 counts."""
import os, sys, warnings
import numpy as np, healpy as hp
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from lowres.regions import NAMES, named_regions, dr1_polarized_amplitude, TAU_A, SURROGATE
z = np.load("nextgen/experiments/E044/lcdm_all.npz")
for ns in (16, 8):
    region, keep, block = named_regions(ns)
    n = region.size // 2
    print(f"nside {ns}: pixels {n}, kept elements {keep.sum()} of {keep.size}; per region (elements) {np.bincount(region[keep], minlength=4)}; pixels per region kept {np.bincount(region[:n][keep[:n]], minlength=4)}; "
          f"masked pixels {(~keep[:n]).sum()}; patches with kept elements {[len(np.unique(block[keep & (region == r)])) for r in range(4)]}")
    if ns == 16:
        assert np.array_equal(keep, z["keep"]) and np.array_equal(region[keep], z["region"])
        print("  identical to the keep and region arrays stored in lcdm_all.npz")
        lon, lat = hp.pix2ang(16, np.arange(n), nest=True, lonlat=True)
        l180 = np.where(lon > 180, lon - 360, lon)
        r, k = region[:n], keep[:n]
        for j in range(4):
            m = r == j
            print(f"  {NAMES[j]:14s}: all pixels {m.sum():5d}, kept {(m & k).sum():5d}; |b| min {np.abs(lat[m]).min():5.1f} max {np.abs(lat[m]).max():5.1f}; b min {lat[m].min():6.1f} max {lat[m].max():6.1f}; l range {l180[m].min():7.1f} {l180[m].max():7.1f}")
        P = hp.reorder(dr1_polarized_amplitude(16), r2n=True)
        thr = np.quantile(P, 0.95)
        bright = P >= thr
        disc = np.zeros(n, bool); disc[hp.ring2nest(16, hp.query_disc(16, hp.ang2vec(*TAU_A, lonlat=True), np.radians(7.5)))] = True
        print(f"  threshold {thr:.1f} uK_RJ; bright pixels {bright.sum()} ({(bright & (np.abs(lat) < 10)).sum()} within 10 deg of the plane; by region of the unmasked definition {np.bincount(r[bright], minlength=4)}); "
              f"Tau A disc {disc.sum()} pixels, {(disc & ~bright).sum()} not already bright; total masked {(bright | disc).sum()}")
        print(f"  share of |b| > 30 pixels that are in High Latitude and kept: {((r == 0) & k).sum()} of {(np.abs(lat) > 30).sum()}; Spur pixels with |b| > 30: {((r == 1) & (np.abs(lat) > 30)).sum()} of {(r == 1).sum()}; masked with |b| > 30: {(~k & (np.abs(lat) > 30)).sum()}")
        sur = hp.reorder(np.load(SURROGATE, allow_pickle=True)["region"].astype(int), r2n=True)
        print(f"  surrogate Spur pixels {(sur == 1).sum()}, of which in the box (b > 20, -45 < l < 65) {((sur == 1) & (lat > 20) & (l180 > -45) & (l180 < 65)).sum()}")
        # maximum kept polarized amplitude per region, in units of threshold
        for j in range(4):
            m = (r == j) & k
            print(f"    kept {NAMES[j]:14s}: DR1 P median {np.median(P[m]):5.1f}, max {P[m].max():5.1f} uK_RJ at 30 GHz")
