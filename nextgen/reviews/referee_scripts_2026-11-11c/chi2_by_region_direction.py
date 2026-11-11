"""Referee: the paper splits the total chi-square into 'two directions per element that the three amplitudes do not span'
(the flat-CMB fit) and 'one direction per element along the CMB spectrum' (LCDM-prior fit minus flat-CMB fit), for the
whole sky only (1.008 +- 0.019 against 0.956 +- 0.013). The same split by region, from the stored chi2_region of the
two fits (sum_b r_R^T (N^-1 r)_R / (n_R n_bands)). Reads only. Approximations: the region parts of a dense quadratic
form are treated as chi-square with one (or two) degrees of freedom per element; the prior penalty (55 in total)
is spread evenly; the two fits have slightly different indices."""
import numpy as np
E = "nextgen/experiments/E044/"
z, f, k = (np.load(E + t + ".npz", allow_pickle=True) for t in ("lcdm_all", "flat", "KKa"))
names = ["High Latitude", "Spur", "Intermediate", "Plane"]
n = z["nper"].astype(float)
pen = float(z["chi2"]) - float(np.sum(z["chi2_region"] * n * 5))
resL, resF = z["chi2_region"] * n * 5, f["chi2_region"] * n * 5
d = resL - resF
print(f"total: LCDM residual part {resL.sum():.1f} + penalty {pen:.1f} = {float(z['chi2']):.1f}; flat {resF.sum():.1f} = {float(f['chi2']):.1f}")
print(f"{'region':14s} {'n':>5s} | two other directions: chi2/2n, z vs 1 | CMB direction: chi2, expected n - penalty share, ratio, z | z if the noise variance is 0.956 of nominal")
for r in range(4):
    e_other, se_other = resF[r] / (2 * n[r]), np.sqrt(2 / (2 * n[r]))
    exp_c = n[r] * (1 - pen / n.sum())
    zc = (d[r] - exp_c) / np.sqrt(2 * n[r])
    zc_lo = (d[r] - n[r]) / np.sqrt(2 * n[r])                 # if none of the penalty belonged to the region
    zc2 = (d[r] - 0.956 * exp_c) / (0.956 * np.sqrt(2 * n[r]))
    print(f"{names[r]:14s} {int(n[r]):5d} | {e_other:.3f} +- {se_other:.3f}  z {(e_other - 1) / se_other:+.2f} | {d[r]:7.1f}  {exp_c:7.1f}  {d[r] / exp_c:.3f} +- {np.sqrt(2 / n[r]):.3f}  z {zc:+.2f} (at least {zc_lo:+.2f}) | {zc2:+.2f}")
print("sum of the CMB-direction parts", d.sum(), "+ penalty", pen, "=", d.sum() + pen, "for", int(n.sum()), "->", (d.sum() + pen) / n.sum())
# total chi2 of each region in the LCDM fit if the noise variance were overstated uniformly by the flat-fit factor
exp_u = 0.956 * (float(z["dof"]) - pen) / (n.sum() * 5)
se = np.sqrt(2.0 * 3 * n) / (n * 5)
print("LCDM-fit chi2 per element and band", np.round(z["chi2_region"], 3), "; expectation if all noise variance is 0.956 of nominal", round(exp_u, 3), "; z", np.round((z["chi2_region"] - exp_u) / (0.956 * se), 2))
# K and Ka only: one degree of freedom per element, along the CMB spectrum as K and Ka see it
penk = float(k["chi2"]) - float(np.sum(k["chi2_region"] * n * 3))
rk = k["chi2_region"] * n * 3
print(f"K, Ka, 353 fit: chi2 {float(k['chi2']):.1f} for {int(k['dof'])} ({float(k['chi2']) / int(k['dof']):.4f}, {(float(k['chi2']) - int(k['dof'])) / np.sqrt(2 * int(k['dof'])):+.2f} sd); penalty {penk:.1f}")
for r in range(4):
    exp_c = n[r] * (1 - penk / n.sum())
    print(f"   {names[r]:14s} {rk[r]:7.1f} expected {exp_c:7.1f} ratio {rk[r] / exp_c:.3f} +- {np.sqrt(2 / n[r]):.3f}  z {(rk[r] - exp_c) / np.sqrt(2 * n[r]):+.2f}")
