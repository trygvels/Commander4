"""Referee check (draft 3): how much of each region's block (CR3) variance comes from the block scores of the
OTHER regions through the off-diagonal elements of J_c^-1. Read only; uses the saved G, D, J, T of the map fit."""
import numpy as np, os
np.set_printoptions(precision=3, suppress=True, linewidth=220)
E40 = "experiments/E040"
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
for btag, ctag in (("block_v3_maps_lcdm.npz", "s4_bpxv_n16_wm_lcdm.npz"), ("block_v3_maps_lcdm_floor011.npz", None)):
    m = np.load(os.path.join(E40, btag), allow_pickle=True)
    if ctag is None:
        continue
    c = np.load(os.path.join(E40, ctag), allow_pickle=True)
    G, D = m["G"], m["D"]
    Jc = c["J"] - c["T"]; Ji = np.linalg.inv(Jc)
    Dp = np.clip(D, 0, None); h = Dp / Dp.sum(0)
    sdm = c["sd_sand"]
    print("J_c^-1 as correlations:\n", Ji / np.sqrt(np.outer(np.diag(Ji), np.diag(Ji))))
    print("sqrt(diag J_c^-1)", np.sqrt(np.diag(Ji)), " 1/sqrt(diag J_c)", 1 / np.sqrt(np.diag(Jc)), " sandwich", sdm)
    for name, f in (("CR0", np.ones_like(h)), ("CR3", 1 / (1 - h))):
        Gk = G * f
        M = Gk.T @ Gk
        print(f"\n{name}: middle matrix M = sum_B Gk_B Gk_B^T, as M_rs / sqrt(Jc_rr Jc_ss) (identity-like under the model):\n", M / np.sqrt(np.outer(np.diag(Jc), np.diag(Jc))))
        for r in range(4):
            C = np.outer(Ji[r], Ji[r]) * M          # C[s,t] = Ji[r,s] M[s,t] Ji[t,r]
            tot = C.sum()
            own = C[r, r]
            # variance if the other regions' scores had exactly their model covariance (J_c) and no cross term with region r beyond the model:
            M0 = Jc.copy(); M0[r, r] = M[r, r]
            v0 = Ji[r] @ M0 @ Ji[r]
            print(f"  {NAMES[r]:16s} block sd {np.sqrt(tot):.4f} | own-score term only {np.sqrt(own):.4f} ({own / tot:.2f} of the variance) | other regions' scores replaced by their model covariance: {np.sqrt(v0):.4f} | ratio to sandwich: full {np.sqrt(tot) / sdm[r]:.2f}, own {np.sqrt(own) / sdm[r]:.2f}, model-others {np.sqrt(v0) / sdm[r]:.2f}")
            print("       contributions C[s,t]/total:\n", C / tot)
    # where do the blocks of the dominant cross terms lie: for High Latitude list blocks by their contribution to the full variance
    Gk = G / (1 - h)
    for r in (0, 1):
        per_block = (Gk @ Ji[r]) ** 2          # (sum_s Ji[r,s] Gk_{B,s})^2: contribution of block B to the full CR3 variance of region r
        order = np.argsort(-per_block)
        print(f"\n{NAMES[r]}: blocks by contribution to the FULL CR3 variance (block, share, and which region's score in that block drives it)")
        for b in order[:6]:
            parts = Ji[r] * Gk[b]
            print(f"   block {b:2d} share {per_block[b] / per_block.sum():.2f}  terms Ji[r,s]*Gk[B,s] by region s: {parts}  h_B by region {h[b]}")
    # exceedance of the own-score ratio in the nside-8 null simulations
    z = np.load(os.path.join(E40, "sim_n8_lcdm_blocknull.npz"), allow_pickle=True)
    ratio3 = z["CR3"] / z["s_l"]; ratio0 = z["CR0"] / z["s_l"]
    for r in range(4):
        own3 = np.sqrt((G[:, r] / (1 - h[:, r])) ** 2).sum() * 0  # placeholder
    Gk3 = G / (1 - h)
    own3 = np.sqrt(np.diag(Ji) ** 2 * (Gk3 ** 2).sum(0)) / sdm
    own0 = np.sqrt(np.diag(Ji) ** 2 * (G ** 2).sum(0)) / sdm
    print("\nown-score ratios on the maps: CR3", own3, " CR0", own0)
    print("null simulations (nside 8, LCDM prior, full-matrix form) above the own-score map ratio: CR3", (ratio3 > own3[None]).sum(0), " CR0", (ratio0 > own0[None]).sum(0), "of", ratio3.shape[0])
    print("null simulations above the full map ratio: CR3", (ratio3 > (m['CR3'] / sdm)[None]).sum(0), " CR0", (ratio0 > (m['CR0'] / sdm)[None]).sum(0))
    z2 = np.load(os.path.join(E40, "block_v2_sim_spread0_full.npz"), allow_pickle=True)
    r3 = z2["CR3"] / z2["sd_model"]; r0 = z2["CR0"] / z2["sd_model"]
    print("white-prior null (nside 8) above the own-score map ratio: CR3", (r3 > own3[None]).sum(0), " CR0", (r0 > own0[None]).sum(0), "of", r3.shape[0])
