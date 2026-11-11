"""Referee check (draft 3): Eq. (block) recomputed from the saved block scores G and informations D of the map
fit, the anatomy of the block errors (which patches carry them), and the null ratios of the nside-8 simulations.
Read only; no fits."""
import numpy as np, os
np.set_printoptions(precision=4, suppress=True, linewidth=220)
E40 = "experiments/E040"
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
m = np.load(os.path.join(E40, "block_v3_maps_lcdm.npz"), allow_pickle=True)
c = np.load(os.path.join(E40, "s4_bpxv_n16_wm_lcdm.npz"), allow_pickle=True)
G, D = m["G"], m["D"]
J, T = c["J"], c["T"]; Jc = J - T; Ji = np.linalg.inv(Jc)
print("column sums of G (should vanish at the optimum):", G.sum(0), " largest |G_B|:", np.abs(G).max(0))
print("sum_B D_B / diag(J_c):", D.sum(0) / np.diag(Jc))
Dp = np.clip(D, 0, None); h = Dp / Dp.sum(0)
neff = Dp.sum(0) ** 2 / (Dp ** 2).sum(0)
print("n_eff", neff, " max h", h.max(0), " blocks with D>0:", (D > 0).sum(0), " with D != 0:", (D != 0).sum(0), " D<0:", (D < 0).sum(0))
for name, f in (("CR0", np.ones_like(h)), ("CR2", 1 / np.sqrt(1 - h)), ("CR3", 1 / (1 - h))):
    Gk = G * f
    sd = np.sqrt(np.diag(Ji @ (Gk.T @ Gk) @ Ji))
    sd_diag = np.sqrt((Gk ** 2).sum(0)) / np.diag(Jc)          # scalar form of Eq. (block), region by region
    sd_obs = np.sqrt((Gk ** 2).sum(0)) / D.sum(0)              # bread = observed Hessian
    print(f"{name}: full-matrix {sd}  stored {m[name]}  scalar per region {sd_diag}  bread sum D {sd_obs}")
sdm = c["sd_sand"]
print("CR3 / sandwich", m["CR3"] / sdm)
print("expected inflation of CR3 under the model with independent blocks, sqrt(sum h/(1-h)):", np.sqrt((h / (1 - h)).sum(0)))
print("expected deflation of CR0, sqrt(1 - sum h^2):", np.sqrt(1 - (h ** 2).sum(0)))
print()
# anatomy: contribution of each block to the CR3 variance (scalar form), sorted
for r in range(4):
    contrib = (G[:, r] / (1 - h[:, r])) ** 2
    tot = contrib.sum()
    order = np.argsort(-contrib)
    z = G[:, r] / np.sqrt(np.where(Dp[:, r] > 0, Dp[:, r], np.nan))    # block score in units of its model sd (if Var G_B = D_B)
    print(f"{NAMES[r]}: blocks by share of the CR3 variance (block, share of variance, leverage h, G_B/sqrt(D_B), one-step index offset -G/D)")
    for b in order[:6]:
        print(f"    block {b:2d}  share {contrib[b]/tot:5.2f}  h {h[b, r]:.3f}  z {z[b]:+6.2f}  one-step {(-G[b, r]/D[b, r]) if D[b, r] > 0 else float('nan'):+.3f}")
    # leave the top block out of the variance sum (not a refit): what the block error would be without its term
    for k in (1, 2):
        keep = np.ones(len(contrib), bool); keep[order[:k]] = False
        print(f"    CR3 without the top {k} term(s) of the sum: {np.sqrt(contrib[keep].sum()) / np.diag(Jc)[r]:.4f} (full {np.sqrt(tot) / np.diag(Jc)[r]:.4f}, sandwich {sdm[r]:.4f})")
    # chi-square of block scores under the model if blocks were independent with Var(G_B) = D_B (1 - h_B)
    ok = Dp[:, r] > 0
    chi = (G[ok, r] ** 2 / (Dp[ok, r] * (1 - h[ok, r]))).sum()
    print(f"    sum_B G_B^2 / (D_B (1 - h_B)) = {chi:.1f} for {ok.sum()} blocks; median of the terms {np.median(G[ok, r] ** 2 / (Dp[ok, r] * (1 - h[ok, r]))):.2f}; terms above 9: {(G[ok, r] ** 2 / (Dp[ok, r] * (1 - h[ok, r])) > 9).sum()}")
print()
# null simulations
for tag, key_s in (("sim_n8_lcdm_blocknull.npz", "s_l"), ("block_v2_sim_spread0_full.npz", "sd_model")):
    z = np.load(os.path.join(E40, tag), allow_pickle=True)
    ratio = z["CR3"] / z[key_s]
    ne = z["neff"]
    print(f"== {tag}: {ratio.shape[0]} realizations")
    print("   CR3/sandwich median", np.median(ratio, 0), " 84th", np.percentile(ratio, 84, 0), " 97.5th", np.percentile(ratio, 97.5, 0), " max", ratio.max(0))
    print("   CR0/sandwich median", np.median(z["CR0"] / z[key_s], 0), " max", (z["CR0"] / z[key_s]).max(0))
    print("   n_eff in the simulations: median", np.median(ne, 0), " min", ne.min(0), " max", ne.max(0), "   maps:", neff)
    if "levmax" in z.files:
        print("   largest leverage: median", np.median(z["levmax"], 0), " max", z["levmax"].max(0), "  maps:", h.max(0))
        print("   blocks with negative information per realization (mean):", z["nneg"].mean(0), " cross share mean", z["cross"].mean())
    print("   simulations with ratio above the map value", (ratio > (m["CR3"] / sdm)[None, :]).sum(0), "of", ratio.shape[0])
    print("   log-ratio: mean", np.log(ratio).mean(0), " sd", np.log(ratio).std(0, ddof=1), " map value in sd of log-ratio:", (np.log(m["CR3"] / sdm) - np.log(ratio).mean(0)) / np.log(ratio).std(0, ddof=1))
    for r in range(4):
        cc = np.corrcoef(ne[:, r], ratio[:, r])[0, 1]
        print(f"     {NAMES[r]}: correlation of ratio with n_eff {cc:+.2f}")
    est = z["b_l"] if "b_l" in z.files else z["est"]
    zt = (est - z["truth"]) 
    print("   rms z with CR3", np.sqrt(((zt / z["CR3"]) ** 2).mean(0)), " with sandwich", np.sqrt(((zt / z[key_s]) ** 2).mean(0)), " with CR0", np.sqrt(((zt / z["CR0"]) ** 2).mean(0)))
    print("   scatter of estimate", zt.std(0, ddof=1), " median sandwich", np.median(z[key_s], 0), " median CR3", np.median(z["CR3"], 0), " mean z", (zt / z[key_s]).mean(0), "sd z", (zt / z[key_s]).std(0, ddof=1))
