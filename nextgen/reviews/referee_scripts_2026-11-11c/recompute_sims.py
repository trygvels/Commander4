"""Referee recomputation of the simulation numbers of Sect. 6.3 (nside-8 nulls, zero-mean prior, chi-square null). Reads only."""
import numpy as np
np.set_printoptions(precision=3, suppress=True, linewidth=170)
E = "nextgen/experiments/E044/"
L = lambda t: np.load(E + t + ".npz", allow_pickle=True)
z = L("lcdm_all")
mapr = z["blk_CR3_own"] / z["sd_sand"]
sets = [L("sim_n8_blocknull"), L("sim_n8_blocknull_b"), L("sim_n8_blocknull_c")]
print("truth", [s["truth"] for s in sets], "first", [int(s["first"]) for s in sets], "n", [s["b"].shape[0] for s in sets])
for nm, grp in (("a", sets[:1]), ("b", sets[1:2]), ("a+b", sets[:2]), ("c", sets[2:]), ("all", sets)):
    zz = np.concatenate([(q["b"] - q["truth"]) / q["s"] for q in grp])
    q_ = np.concatenate([q["CR3_own"] / q["s"] for q in grp])
    m = zz.shape[0]
    print(f"set {nm:4s} n={m}: mean z {zz.mean(0)} (se {1/np.sqrt(m):.3f}); sd z {zz.std(0, ddof=1)} (se {1/np.sqrt(2*(m-1)):.3f}); cover1 {np.mean(np.abs(zz)<1, axis=0)}; cover2 {np.mean(np.abs(zz)<2, axis=0)}")
    print(f"          ratio median {np.nanmedian(q_, axis=0)} p84 {np.nanpercentile(q_, 84, axis=0)} p97.5 {np.nanpercentile(q_, 97.5, axis=0)} max {np.nanmax(q_, axis=0)} above map {[(q_[:, j] > mapr[j]).sum() for j in range(4)]} nan {np.isnan(q_).sum(0)}")
    print(f"          mean sandwich sd {np.concatenate([q['s'] for q in grp]).mean(0)}; sd of estimates {np.concatenate([q['b'] - q['truth'] for q in grp]).std(0, ddof=1)}; mean b-truth {np.concatenate([q['b'] - q['truth'] for q in grp]).mean(0)}")
allb = np.concatenate([q["b"] for q in sets])
print("duplicates between sets (same realization twice)?", len(np.unique(np.round(allb[:, 0], 10))), "unique of", allb.shape[0])
print("map ratios", mapr)
zm = L("sim_n8_zm")
dz = (zm["b"] - zm["truth"]) / zm["s"]; dp = (zm["prof"] - zm["truth"]) / zm["s"]; pr = (zm["b"] - zm["prof"]) / zm["s"]
m = dz.shape[0]
print("== zero-mean sim: n", m, "truth", zm["truth"], "first", zm["first"])
print("zm mean z", dz.mean(0), "+-", dz.std(0, ddof=1) / np.sqrt(m), " pooled", dz.mean(), "+-", dz.std(ddof=1) / np.sqrt(dz.size))
print("profile mean z in the same sims", dp.mean(0), "+-", dp.std(0, ddof=1) / np.sqrt(m), " sd z", dp.std(0, ddof=1))
print("paired zm - prof in sd", pr.mean(0), "scatter", pr.std(0, ddof=1), " in index", (zm["b"] - zm["prof"]).mean(0), "zm-truth index", (zm["b"] - zm["truth"]).mean(0))
print("mean sandwich sd in zm sims", zm["s"].mean(0), " (map fit:", z["sd_sand"], ")")
# are the zm sims the same realizations as blocknull set a?
a = sets[0]
print("zm sims identical to first 30 of null set a (profile optima)?", np.allclose(zm["prof"], a["b"][:m]), np.abs(zm["prof"] - a["b"][:m]).max())
c2 = L("sim_n8_chi2null")
print("== chi2 null: n", c2["ratio"].size, "mean", c2["ratio"].mean(), "se", c2["ratio"].std(ddof=1) / np.sqrt(c2["ratio"].size), "sd", c2["ratio"].std(ddof=1), "theory", np.sqrt(2.0 / int(c2["dof"])), "dof", int(c2["dof"]), "nper", c2["nper"])
print("   creg mean", c2["creg"].mean(0), "penalty mean", c2["penalty"].mean(), " penalty per kept element", c2["penalty"].mean() / c2["nper"].sum(), " (map: 55.5 / 5808 =", 55.458 / 5808, ")")
print("   min/max ratio", c2["ratio"].min(), c2["ratio"].max())
