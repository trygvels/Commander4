"""Referee check (draft 3): where the patches that dominate the block errors lie on the sky (Galactic l, b of the
nside-2 NESTED pixel centres as used by the saved `block` array), and which known compact polarized sources fall
in them. Read only."""
import numpy as np, os
import healpy as hp
E40 = "experiments/E040"
m = np.load(os.path.join(E40, "block_v3_maps_lcdm.npz"), allow_pickle=True)
block, region = m["block"], m["region"]
n = block.size // 2
nside = hp.npix2nside(n)
print("elements", block.size, "nside", nside, " block ids", block.min(), block.max(), " Q and U share blocks:", np.array_equal(block[:n], block[n:]))
# is the saved block array the nested nside-2 parent of nested or ring nside-16 pixels?
pix = np.arange(n)
par_nest = pix // 64                                   # if the maps are in NESTED order
par_ring = hp.ring2nest(nside, pix) // 64              # if the maps are in RING order (parent in nested nside 2)
th, ph = hp.pix2ang(nside, pix)                        # RING
par_ring_ring = hp.ang2pix(2, th, ph)                  # ring-ordered nside-2 parent
print("block == nested parent of nested pixels:", np.array_equal(block[:n], par_nest),
      "; == nested parent of ring pixels:", np.array_equal(block[:n], par_ring),
      "; == ring parent of ring pixels:", np.array_equal(block[:n], par_ring_ring))
order_nest = np.array_equal(block[:n], par_nest)
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
src = {"Tau A (Crab)": (184.56, -5.78), "Cen A": (309.52, 19.42), "Cyg A": (76.19, 5.76), "Cas A": (111.73, -2.13), "Galactic Center": (0.0, 0.0),
       "Fan region centre": (140.0, 5.0), "3C 273": (289.95, 64.36), "Vir A": (283.78, 74.49), "For A": (240.16, -56.69), "Orion A": (209.0, -19.4), "W51": (49.2, -0.7)}
for name, (l, b) in src.items():
    p = hp.ang2pix(2, l, b, nest=True, lonlat=True) if (order_nest or np.array_equal(block[:n], par_ring)) else hp.ang2pix(2, l, b, lonlat=True)
    print(f"  {name:20s} l {l:6.1f} b {b:+6.1f} -> block {p}")
G, D = m["G"], m["D"]
Dp = np.clip(D, 0, None); h = Dp / Dp.sum(0)
for r, blks in ((0, (18, 15, 17, 22, 12)), (1, (0, 16, 19, 2)), (2, (18, 19, 16, 17)), (3, (24, 21, 4, 26))):
    for bk in blks:
        nestflag = order_nest or np.array_equal(block[:n], par_ring)
        l, b = hp.pix2ang(2, bk, nest=nestflag, lonlat=True)
        nel = ((block == bk) & (region == r)).sum()
        allreg = [int(((block == bk) & (region == k)).sum()) for k in range(4)]
        print(f"{NAMES[r]:16s} block {bk:2d} centre l {l:6.1f} b {b:+6.1f}  elements of this region {nel:3d}; elements by region {allreg}; h {h[bk, r]:.3f}; z {G[bk, r] / np.sqrt(Dp[bk, r]) if Dp[bk, r] > 0 else float('nan'):+.2f}")
