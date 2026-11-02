import sys; sys.argv = [sys.argv[0]]
exec(open("marg.py").read().split('print("== faint-level optimum')[0])
c = np.load("chain_summary.npz"); bm = c["bm_mean"]
r = np.load("/Users/trygvels/code/Commander4/nextgen/experiments/E023/scan_n1_refs.npz"); ref_s = r["bs_fixT"]
amu_d = np.mean([state_of(f)["ThermalDust"]["amu"] for f in files[::20]], 0); bd_pm = np.mean([state_of(f)["ThermalDust"]["beta"] for f in files[::20]], 0)
ls_s, ls_d = np.log(0.0342), np.log(0.0219); fm = ~bright; faint = fm.astype(float)
grid = np.linspace(-0.30, 0.30, 25); i0 = 12

# (1) brute-force check of the analytic mu_r integral on one faint region and one bright region
def mu_terms_full(bs_reg, bd_reg, ls_s_, mu_d, ls_d_, d):
    M = mixing(bs_reg, bd_reg)[reg]
    P = np.zeros((M.shape[0], 3)); P[:, 0] = 1e-8; P[:, 1] = 1.0 / (np.exp(ls_s_) * g_s)**2; P[:, 2] = 1.0 / (np.exp(ls_d_) * g_d)**2
    Wt = inv_var.T; A = np.einsum("pbi,pb,pbj->pij", M, Wt, M); Ap = A + np.einsum("pi,ij->pij", P, np.eye(3)); Ainv = np.linalg.inv(Ap)
    m0 = np.zeros((M.shape[0], 3)); m0[:, 2] = mu_d[reg] * g_d
    b0 = np.einsum("pbi,pb,pb->pi", M, Wt, d.T) + P * m0
    v = np.zeros((M.shape[0], 3)); v[:, 1] = P[:, 1] * g_s
    Q = np.bincount(reg, np.einsum("pi,pij,pj->p", v, Ainv, v) - P[:, 1] * g_s**2, 3072) - 1 / 0.25
    L = np.bincount(reg, np.einsum("pi,pij,pj->p", v, Ainv, b0), 3072) + 1 / 0.25
    return Q, L
Q, L = mu_terms_full(bm, bd_pm, ls_s, amu_d, ls_d, maps)
for rr in (int(np.argsort(-W)[5]), int(np.argsort(-W)[2000])):
    mus = np.linspace(-2, 4, 61); vals = []
    for m in mus:
        mu_map = np.ones(3072); mu_map[rr] = m
        vals.append(data_term(bm, bd_pm, mu_map, ls_s, amu_d, ls_d) - 0.5 * ((m - 1) / 0.5)**2)
    vals = np.array(vals); A2 = np.polyfit(mus, vals, 2)
    print(f"region {rr} (W rank {np.where(np.argsort(-W)==rr)[0][0]}): brute-force quadratic coeff {2*A2[0]:+.4f} vs analytic Q {Q[rr]:+.4f}; mode {-A2[1]/(2*A2[0]):+.3f} vs analytic {-L[rr]/Q[rr]:+.3f}; cubic residual {np.abs(vals - np.polyval(A2, mus)).max():.2e}")

# (2) Occam term slope along the faint level (bright at chain mean), and its region split; data curvature from the scan sd
def occam(bs):
    Qf, _ = mu_terms_full(bs, bd_pm, ls_s, amu_d, ls_d, maps); return -0.5 * np.log(-Qf)
base = np.where(bright, bm, ref_s)
o_m, o_0, o_p = occam(base - 0.1 * faint), occam(base), occam(base + 0.1 * faint)
slope = (o_p - o_m) / 0.2
print(f"Occam term -1/2 log(-Q_r): slope along the faint level {slope[fm].sum():+.1f} nats per unit beta (faint regions), {slope[bright].sum():+.1f} (bright, not shifted here would be 0; shown as slope wrt their own beta)")
o_all_m, o_all_p = occam(base - 0.1), occam(base + 0.1)
print(f"Occam slope along a whole-sky shift: total {((o_all_p - o_all_m)/0.2).sum():+.1f}; bright {((o_all_p - o_all_m)/0.2)[bright].sum():+.1f}; faint {((o_all_p - o_all_m)/0.2)[fm].sum():+.1f}")
print(f"predicted faint-level shift from the Occam term alone: +{slope[fm].sum()} / (1/0.065^2 = {1/0.065**2:.0f}) = {slope[fm].sum()*0.065**2:+.3f}; predicted scalar-index shift: {((o_all_p - o_all_m)/0.2).sum()*0.025**2:+.3f}")
# faint-level scan with the Occam term removed (profile with the hyperprior, correctly charged)
def scan_mode(bs0, d, mode):
    out = []
    for dd in grid:
        bs = bs0 + dd * faint; Qf, Lf = mu_terms_full(bs, bd_pm, ls_s, amu_d, ls_d, d); mu_hat = -Lf / Qf
        val = data_term(bs, bd_pm, mu_hat, ls_s, amu_d, ls_d, d=d) - 0.5 * np.sum(((mu_hat - 1) / 0.5)**2)
        if mode == "marginal": val += -0.5 * np.sum(np.log(-Qf))
        out.append(val)
    return np.array(out)
for mode in ("profile", "marginal"):
    lp = scan_mode(base, maps, mode); m, s = fit(grid, lp); print(f"faint level, realization, mu_r {mode:8s} (hyperprior charged): optimum {m:+.3f} (sd {s:.3f}) -> {ref_s[fm].mean()+m:+.3f}")
# (3) noise realizations
rng = np.random.default_rng(5); opts = []
for k in range(6):
    dsim = sky + rng.standard_normal(sky.shape) * t["rms"][:, None]
    lp = scan_mode(base, dsim, "marginal"); opts.append(fit(grid, lp)[0])
opts = np.array(opts); print(f"faint level, simulated noise realizations, mu_r marginal: optima {np.round(opts,3)}; mean {opts.mean():+.3f} sd {opts.std():.3f} -> level {ref_s[fm].mean()+opts.mean():+.3f}")
lp = scan_mode(base, sky, "marginal"); print(f"faint level, noiseless data plugged in, mu_r marginal: optimum {fit(grid, lp)[0]:+.3f}")
lp = scan_mode(base, sky, "profile"); print(f"faint level, noiseless data plugged in, mu_r profile: optimum {fit(grid, lp)[0]:+.3f}")
print(f"mean 23 GHz noise in the faint set (se units): {(maps-sky)[0][fm[reg]].mean()/(280/np.sqrt(fm[reg].sum())):+.2f}; bright set: {(maps-sky)[0][bright[reg]].mean()/(280/np.sqrt(bright[reg].sum())):+.2f}")
# (4) whole-sky shift at the posterior-mean field with mu_r integrated, and at the reference
for bs0, lab in ((bm, "posterior-mean field"), (ref_s, "reference field")):
    lp = np.array([marg_terms(bs0 + d, bd_pm, ls_s, amu_d, ls_d, maps) for d in grid]); m, s = fit(grid, lp)
    print(f"whole-sky shift, {lab}, mu_r integrated, realization: optimum {m:+.3f} (sd {s:.3f}, gain {lp.max()-lp[i0]:.1f}) -> Fisher-weighted level {np.sum(W*bs0)/W.sum()+m:+.4f}")
# (5) global multiplier (one mu for all regions): profile of the quadratic form
def glob(bs0, sel):
    out = []
    for d in grid:
        bs = bs0 + d * sel; Qf, Lf = mu_terms_full(bs, bd_pm, ls_s, amu_d, ls_d, maps)
        Qg = (Qf + 1 / 0.25).sum() - 1 / 0.25; Lg = (Lf - 1 / 0.25).sum() + 1 / 0.25; mu_hat = -Lg / Qg
        out.append(data_term(bs, bd_pm, np.full(3072, mu_hat), ls_s, amu_d, ls_d) - 0.5 * ((mu_hat - 1) / 0.5)**2 - 0.5 * np.log(-Qg))
    return np.array(out)
for bs0, lab in ((bm, "posterior-mean field"), (ref_s, "reference field")):
    lp = glob(bs0, np.ones(3072)); m, s = fit(grid, lp); print(f"whole-sky shift, {lab}, ONE global multiplier integrated: optimum {m:+.3f} (sd {s:.3f}, gain {lp.max()-lp[i0]:.1f})")
lp = glob(base, faint); m, s = fit(grid, lp); print(f"faint level, ONE global multiplier integrated, faint at ref: optimum {m:+.3f} (sd {s:.3f}) -> {ref_s[fm].mean()+m:+.3f}")
