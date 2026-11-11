"""Referee check (draft 3): a coarse picture of the four regions as used in the fits (from the saved `region`
array of the central fit, NESTED nside 16), in Galactic coordinates, and of the fitted polarized amplitude.
H High Latitude, S Spur, C Galactic Center, P Plane. Read only."""
import numpy as np, healpy as hp
c = np.load("experiments/E040/s4_bpxv_n16_wm_lcdm.npz", allow_pickle=True)
reg, a = c["region"], c["a_s"]
n = reg.size // 2
P = np.hypot(a[:n], a[n:])
ch = "HSCP"
ls = np.arange(180, -180.1, -7.5)
print("l from +180 (left) through 0 (centre) to -180 (right), step 7.5 deg; rows b from +82.5 to -82.5 step 7.5")
for b in np.arange(82.5, -83, -7.5):
    row = "".join(ch[reg[hp.ang2pix(16, (l % 360), b, nest=True, lonlat=True)]] for l in ls)
    print(f"{b:+6.1f} {row}")
l, b = hp.pix2ang(16, np.arange(n), nest=True, lonlat=True)
for r in range(4):
    s = reg[:n] == r
    print(f"{ch[r]}: pixels {s.sum()}, |b| min {np.abs(b[s]).min():.1f} median {np.median(np.abs(b[s])):.1f} max {np.abs(b[s]).max():.1f}; pixels with |b| < 5: {(np.abs(b[s]) < 5).sum()}, with |b| < 10: {(np.abs(b[s]) < 10).sum()}; share of sum P^2 at |b| < 10: {(P[s][np.abs(b[s]) < 10] ** 2).sum() / (P[s] ** 2).sum():.2f}, at |b| < 20: {(P[s][np.abs(b[s]) < 20] ** 2).sum() / (P[s] ** 2).sum():.2f}, at |b| > 30: {(P[s][np.abs(b[s]) > 30] ** 2).sum() / (P[s] ** 2).sum():.2f}")
s = (reg[:n] == 0) & (np.abs(b) < 10)
print("High Latitude pixels with |b| < 10: longitudes", np.sort(np.round(l[s], 1)))
