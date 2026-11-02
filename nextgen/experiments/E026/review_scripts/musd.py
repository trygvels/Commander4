import sys; sys.argv = [sys.argv[0]]
exec(open("marg.py").read().split('print("== faint-level optimum')[0])
c = np.load("chain_summary.npz"); bm = c["bm_mean"]
r = np.load("/Users/trygvels/code/Commander4/nextgen/experiments/E023/scan_n1_refs.npz"); ref_s = r["bs_fixT"]
amu_d = np.mean([state_of(f)["ThermalDust"]["amu"] for f in files[::20]], 0); bd_pm = np.mean([state_of(f)["ThermalDust"]["beta"] for f in files[::20]], 0)
ls_s, ls_d = np.log(0.0342), np.log(0.0219); fm = ~bright; faint = fm.astype(float); base = np.where(bright, bm, ref_s)
grid = np.linspace(-0.30, 0.30, 25); i0 = 12
def marg_sd(bs_reg, sd_mu, d):
    M = mixing(bs_reg, bd_pm)[reg]
    P = np.zeros((M.shape[0], 3)); P[:, 0] = 1e-8; P[:, 1] = 1.0 / (np.exp(ls_s) * g_s)**2; P[:, 2] = 1.0 / (np.exp(ls_d) * g_d)**2
    Wt = inv_var.T; A = np.einsum("pbi,pb,pbj->pij", M, Wt, M); Ap = A + np.einsum("pi,ij->pij", P, np.eye(3)); Ainv = np.linalg.inv(Ap)
    m0 = np.zeros((M.shape[0], 3)); m0[:, 2] = amu_d[reg] * g_d
    b0 = np.einsum("pbi,pb,pb->pi", M, Wt, d.T) + P * m0
    v = np.zeros((M.shape[0], 3)); v[:, 1] = P[:, 1] * g_s
    Q = np.bincount(reg, np.einsum("pi,pij,pj->p", v, Ainv, v) - P[:, 1] * g_s**2, 3072) - 1 / sd_mu**2
    L = np.bincount(reg, np.einsum("pi,pij,pj->p", v, Ainv, b0), 3072) + 1 / sd_mu**2
    base_ = 0.5 * np.einsum("pi,pij,pj->p", b0, Ainv, b0).sum() - 0.5 * np.linalg.slogdet(Ap)[1].sum() + 0.5 * np.sum(np.log(P[:, 1:])) - 0.5 * np.sum(P[:, 2] * m0[:, 2]**2)
    mu_hat = -L / Q
    return base_ + np.sum(L * mu_hat + 0.5 * Q * mu_hat**2) - 0.5 * np.sum(np.log(-Q)), -0.5 * np.sum(np.log(-Q[fm]))
rng = np.random.default_rng(5); sims = [sky + rng.standard_normal(sky.shape) * t["rms"][:, None] for _ in range(3)]
for sd_mu in (0.5, 0.2, 0.1, 0.03):
    vals = np.array([marg_sd(base + dd * faint, sd_mu, maps) for dd in grid]); lp, oc = vals[:, 0], vals[:, 1]
    m, s = fit(grid, lp); osl = (oc[i0 + 2] - oc[i0 - 2]) / (grid[i0 + 2] - grid[i0 - 2])
    sim_opt = [fit(grid, np.array([marg_sd(base + dd * faint, sd_mu, ds)[0] for dd in grid]))[0] for ds in sims]
    print(f"mu_r hyperprior sd {sd_mu:.2f}: faint-level optimum {m:+.3f} (sd {s:.3f}) -> {ref_s[fm].mean()+m:+.3f}; Occam slope (faint) {osl:+.1f} nats/unit beta; simulated realizations {np.round(sim_opt,3)}")
