"""Referee check (draft 3): is the Plane patch that carries most of the Plane block error (nside-2 nested pixel 24)
dominated by one compact source? Uses the saved fitted synchrotron amplitudes and residuals of the central fit
(s4_bpxv_n16_wm_lcdm.npz). Read only."""
import numpy as np, os
import healpy as hp
np.set_printoptions(precision=2, suppress=True, linewidth=220)
c = np.load("experiments/E040/s4_bpxv_n16_wm_lcdm.npz", allow_pickle=True)
m = np.load("experiments/E040/block_v3_maps_lcdm.npz", allow_pickle=True)
a, reg, res, blk = c["a_s"], c["region"], c["resid"], m["block"]
n = a.size // 2
P = np.hypot(a[:n], a[n:])                      # polarized amplitude at 23 GHz, uK (brightness)
l, b = hp.pix2ang(16, np.arange(n), nest=True, lonlat=True)
print("regions of Q and U elements identical:", np.array_equal(reg[:n], reg[n:]))
for name, (ls, bs) in {"Tau A": (184.56, -5.78), "Cen A": (309.52, 19.42), "Cyg A": (76.19, 5.76), "Cas A": (111.73, -2.13), "3C 273": (289.95, 64.36), "Vir A": (283.78, 74.49), "3C 279": (305.10, 57.06)}.items():
    p = hp.ang2pix(16, ls, bs, nest=True, lonlat=True)
    same_reg = reg[:n] == reg[p]
    rank = (P[same_reg] > P[p]).sum() + 1
    print(f"{name:7s} pixel {p:4d} block {blk[p]:2d} region {reg[p]}  P {P[p]:7.1f} uK (rank {rank} of {same_reg.sum()} in its region; region median {np.median(P[same_reg]):.1f}, max {P[same_reg].max():.1f})"
          f"  residual Q,U by band (uK): " + "  ".join(f"[{res[k, p]:+.1f},{res[k, p + n]:+.1f}]" for k in range(res.shape[0])))
# Plane block 24: the pixels, their amplitudes and residuals
sel = (blk[:n] == 24) & (reg[:n] == 3)
idx = np.flatnonzero(sel)
order = idx[np.argsort(-P[idx])]
print("\nPlane elements in block 24:", sel.sum(), "pixels; listing the six brightest: pixel, l, b, P(23 GHz), residual |r| in Ka, Q, V (quadrature of Q and U)")
for p in order[:6]:
    print(f"   pix {p:4d} l {l[p]:6.1f} b {b[p]:+5.1f} P {P[p]:7.1f}  resid K {np.hypot(res[0,p],res[0,p+n]):.2f} Ka {np.hypot(res[1,p],res[1,p+n]):.2f} Q {np.hypot(res[2,p],res[2,p+n]):.2f} V {np.hypot(res[3,p],res[3,p+n]):.2f}")
print("   share of sum P^2 of the Plane part of block 24 in its brightest pixel:", (P[order[0]] ** 2) / (P[idx] ** 2).sum())
print("   typical residual rms over all Plane elements by band:", np.sqrt((res[:, np.r_[np.flatnonzero(reg[:n]==3), np.flatnonzero(reg[:n]==3)+n]] ** 2).mean(1)))
# crude weight of each pixel in the region index: a^2 (information is roughly amplitude^2 x lever arm / noise)
for r, nm in enumerate(["High Latitude", "Spur", "Galactic Center", "Plane"]):
    s = reg[:n] == r
    w = P[s] ** 2
    ws = np.sort(w)[::-1]
    print(f"{nm:16s} pixels {s.sum():4d}; effective number of pixels by P^2: {w.sum() ** 2 / (w ** 2).sum():6.1f}; share of sum P^2 in the brightest 1, 5, 10 pixels: {ws[0] / w.sum():.2f} {ws[:5].sum() / w.sum():.2f} {ws[:10].sum() / w.sum():.2f}; |b| of the 5 brightest: {np.abs(b[s][np.argsort(-w)[:5]])}, l: {l[s][np.argsort(-w)[:5]]}")
