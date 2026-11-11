"""Referee check (draft 3): first-order delete-one-patch shifts of each region index, G_B / (D - D_B), from the saved
block scores and informations of the central fit, with the sky position of each patch; and the share of each region's
information by Galactic latitude of the patch centre. Read only; not a refit."""
import numpy as np, healpy as hp
np.set_printoptions(precision=3, suppress=True, linewidth=200)
m = np.load("experiments/E040/block_v3_maps_lcdm.npz", allow_pickle=True)
c = np.load("experiments/E040/s4_bpxv_n16_wm_lcdm.npz", allow_pickle=True)
G, D, beta, sd = m["G"], m["D"], m["beta"], m["sd_model"]
Dp = np.clip(D, 0, None); h = Dp / Dp.sum(0)
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
lb = np.array([hp.pix2ang(2, b, nest=True, lonlat=True) for b in range(48)])
reg, blk = m["region"], m["block"]; n = reg.size // 2
l16, b16 = hp.pix2ang(16, np.arange(n), nest=True, lonlat=True)
for r in range(4):
    tot = Dp[:, r].sum()
    shift = G[:, r] / (tot - Dp[:, r])
    order = np.argsort(-h[:, r])
    print(f"\n{NAMES[r]}: index {beta[r]:.3f} +- {sd[r]:.3f} (sandwich). Patches by share of information:")
    cum = 0
    for bk in [x for x in order[:8] if h[x, r] > 0]:
        cum += h[bk, r]
        s = (blk[:n] == bk) & (reg[:n] == r)
        print(f"  block {bk:2d} centre (l,b)=({lb[bk,0]:5.1f},{lb[bk,1]:+5.1f}) share {h[bk,r]:.3f} (cum {cum:.2f}); pixels of the region in it {s.sum():3d}, their |b| range {np.abs(b16[s]).min():4.1f} to {np.abs(b16[s]).max():4.1f}, median {np.median(np.abs(b16[s])):4.1f}; "
              f"deleting it moves the index by {shift[bk]:+.3f} ({shift[bk]/sd[r]:+.1f} sandwich sd)")
    eq = np.abs(lb[:, 1]) < 1.0
    print(f"  share of the region's information in the eight patches centred on b = 0: {h[eq, r].sum():.2f}; in patches centred at |b| = 19.5: {h[np.abs(np.abs(lb[:,1]) - 19.47) < 1, r].sum():.2f}; at |b| >= 41.8: {h[np.abs(lb[:,1]) > 40, r].sum():.2f}")
    big = np.argsort(-np.abs(shift))[:4]
    print("  largest delete-one shifts:", ", ".join(f"block {bk} ({lb[bk,0]:.0f},{lb[bk,1]:+.0f}): {shift[bk]:+.3f}" for bk in big))
