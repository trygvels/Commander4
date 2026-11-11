"""Referee check (draft 3): recompute the numbers of the abstract, Sect. 6 and the conclusions from the .npz files.
Read only; no fits."""
import numpy as np, os
E40 = "experiments/E040"; E42 = "experiments/E042"
np.set_printoptions(precision=4, suppress=True, linewidth=200)
def L(f): return np.load(os.path.join(E40, f), allow_pickle=True)
R = {"lcdm": L("s4_bpxv_n16_wm_lcdm.npz"), "white": L("bpxv_n16_wm.npz"), "flat": L("s4_bpxv_n16_wm_cmbflat.npz")}
ln = np.log(30/23)
for k, z in R.items():
    b, s, f23, f30, pm = z["beta"], z["sd_sand"], z["flat23"], z["flat30"], z["pred_marg"]
    J, T = z["J"], z["T"]
    Jc = J - T
    predJc = np.linalg.solve(Jc, z["nper"] * ln)
    print(f"=== CMB prior {k}: chi2 {float(z['chi2']):.1f} dof {int(z['dof'])} ratio {float(z['chi2'])/int(z['dof']):.4f}")
    print(" profile      ", b, " sd_sand", s, " sd_naive", z["sd_naive"])
    print(" flat23       ", f23, " flat30", f30, " zm", z["zm"])
    print(" flat23-prof  ", f23 - b, " in sd", (f23 - b) / s)
    print(" flat30-prof  ", f30 - b, " in sd", (f30 - b) / s)
    print(" flat30-flat23", f30 - f23, " in sd", (f30 - f23) / s)
    print(" pred H       ", pm, " measured/pred", (f30 - f23) / pm, " percent off", 100 * (1 - (f30 - f23) / pm))
    print(" pred J_c     ", predJc, " measured/pred", (f30 - f23) / predJc)
    print(" zm-prof in sd", (z["zm"] - b) / s)
    # pivot to first order, two ways: measured slope and predicted slope
    piv_meas = 23 * np.exp(-(f23 - b) / (f30 - f23) * ln)
    piv_pred = 23 * np.exp(-(f23 - b) / pm * ln)
    print(" pivot (measured slope)", piv_meas, " pivot (predicted slope)", piv_pred)
    print(" T/J_c diag   ", np.diag(T) / np.diag(Jc), " T/J", np.diag(T)/np.diag(J))
    # is the stored cov_sand Jc^-1 J Jc^-1 ?
    cs = np.linalg.inv(Jc) @ J @ np.linalg.inv(Jc)
    print(" sd from Jc^-1 J Jc^-1", np.sqrt(np.diag(cs)), " stored cov_sand sd", np.sqrt(np.diag(z["cov_sand"])))
    cc = z["cov_sand"]; d = np.sqrt(np.diag(cc)); print(" correlation of the region indices (sandwich):\n", cc / np.outer(d, d))
b0, s0 = R["lcdm"]["beta"], R["lcdm"]["sd_sand"]
print("\n=== moves relative to the LCDM baseline")
print(" white - lcdm", R["white"]["beta"] - b0, " in lcdm sd", (R["white"]["beta"] - b0) / s0)
print(" flat  - lcdm", R["flat"]["beta"] - b0, " in lcdm sd", (R["flat"]["beta"] - b0) / s0)
print(" error ratio flat/lcdm", R["flat"]["sd_sand"] / s0, " flat/white", R["flat"]["sd_sand"] / R["white"]["sd_sand"], " lcdm/white", s0 / R["white"]["sd_sand"])
nest_flat = np.sqrt(R["flat"]["sd_sand"] ** 2 - s0 ** 2)
print(" nested sd flat vs lcdm at nside 16", nest_flat, " move in it", (R["flat"]["beta"] - b0) / nest_flat)
print(" move in the nside-8 paired scatter (0.031,0.037,0.021,0.034)", (R["flat"]["beta"] - b0) / np.array([0.031, 0.037, 0.021, 0.034]))
print(" white-lcdm in nside-8 paired scatter (0.0073,0.0093,0.0028,0.0071)", (R["white"]["beta"] - b0) / np.array([0.0073, 0.0093, 0.0028, 0.0071]))
w0 = L("bpxv_n16.npz")
print(" no window matching (white): unmatched - matched", w0["beta"] - R["white"]["beta"])
print("\n=== arms with the LCDM prior")
for tag in ("l_bpxv_n16_KKa_wm", "l_bpxv_n16_bd148_wm", "l_bpxv_n16_bd158_wm", "l_bpxv_n16_nodrift_wm", "l_bpxv_n16_bp353full_wm", "l_bpxv_n16_delta_wm", "l_bpxv_n16_start27_wm"):
    z = L(tag + ".npz"); b, s = z["beta"], z["sd_sand"]
    extra = ""
    if "KKa" in tag:
        extra = f" nested sd {np.sqrt(s**2 - s0**2)} move/nested {(b - b0) / np.sqrt(s**2 - s0**2)}"
    print(f" {tag:26s} {b} sd {s} move {b - b0} in sd {(b - b0) / s0}{extra}")
print("\n=== low-pass with the LCDM prior")
lp = {}
for Lc in (30, 20, 12):
    z = np.load(os.path.join(E42, f"lowpass{Lc}_lcdm.npz"), allow_pickle=True)
    four, pair = z["K_Ka_Q_V_353"], z["K_Ka_353"]
    lp[Lc] = (four[:4], four[4:], pair[:4], pair[4:])
    nest = np.sqrt(four[4:] ** 2 - s0 ** 2)
    print(f" L={Lc}: four {four[:4]} sd {four[4:]}  minus all {four[:4] - b0}  /nested {(four[:4] - b0) / nest}")
    print(f"        K,Ka {pair[:4]} sd {pair[4:]}  K,Ka - four {pair[:4] - four[:4]}  /nested {(pair[:4] - four[:4]) / np.sqrt(pair[4:]**2 - four[4:]**2)}")
for a, b in ((30, 20), (20, 12), (30, 12)):
    d = lp[b][0] - lp[a][0]; n = np.sqrt(lp[b][1] ** 2 - lp[a][1] ** 2)
    print(f" cuts {a} vs {b}: diff {d} nested sd {n} ratio {d / n}")
# white-prior low-pass for the 'same within 0.025' statement
for Lc in (30, 20, 12):
    f = os.path.join(E42, f"lowpass{Lc}.npz")
    z = np.load(f, allow_pickle=True)
    print(f" white-prior lowpass{Lc} keys {z.files}")
    k = [x for x in z.files if "Q_V" in x and "30" not in x]
    four = z["K_Ka_Q_V_353"]
    dw = four[:4] - R["white"]["beta"]
    dl = lp[Lc][0] - b0
    print(f"   white: four {four[:4]} minus all-scales white {dw}; lcdm shift {dl}; difference of shifts {dl - dw}")
print("\n=== K with 30 GHz")
for tag in ("l_K30_wm", "l_K30_wm_white", "l_K30_wm_a2", "s2_K30_wm", "bpxv_n16_KKa_wm"):
    z = L(tag + ".npz"); print(f" {tag:18s} {z['beta']} sd {z['sd_sand']} chi2/dof {float(z['chi2']) / int(z['dof']):.3f}")
z = L("s3_bpset_cmbflat_wm.npz"); print(" bpset flat:", z["beta"], z["sd_sand"], z["flat23"], z["flat30"], float(z["chi2"]) / int(z["dof"]))
z = L("s2_noK_lfi_wm.npz"); print(" bpset white (profile only):", z["beta"], z["sd_sand"], [k for k in z.files])
