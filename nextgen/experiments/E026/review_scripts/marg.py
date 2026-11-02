import sys; sys.argv = [sys.argv[0]]
exec(open("decomp.py").read().split("# ---- verification")[0])
bs_t, a_s, nu, rms = t["beta_s"], t["a_sync"], t["nu"], t["rms"]; x_ = np.log(nu / 30.0)
w_known = ((a_s[:, None] * (nu / 30.0)[None, :] ** bs_t[:, None] * x_[None, :])**2 / rms[None, :]**2).sum(1)
W = np.bincount(reg, w_known, 3072); order = np.argsort(-W); bright = np.zeros(3072, bool); bright[order[:768]] = True
faint = (~bright).astype(float)
def marg_terms(bs_reg, bd_reg, ls_s, mu_d, ls_d, d):
    """Exact log of the target integrated over the sync multipliers mu_r (Gaussian in mu_r), at
    fixed beta fields, s, dust multipliers; data d (12, npix)."""
    M = mixing(bs_reg, bd_reg)[reg]
    P = np.zeros((M.shape[0], 3)); P[:, 0] = 1e-8; P[:, 1] = 1.0 / (np.exp(ls_s) * g_s)**2; P[:, 2] = 1.0 / (np.exp(ls_d) * g_d)**2
    Wt = inv_var.T; A = np.einsum("pbi,pb,pbj->pij", M, Wt, M); Ap = A + np.einsum("pi,ij->pij", P, np.eye(3)); Ainv = np.linalg.inv(Ap)
    m0 = np.zeros((M.shape[0], 3)); m0[:, 2] = mu_d[reg] * g_d
    b0 = np.einsum("pbi,pb,pb->pi", M, Wt, d.T) + P * m0
    v = np.zeros((M.shape[0], 3)); v[:, 1] = P[:, 1] * g_s
    Q = np.bincount(reg, np.einsum("pi,pij,pj->p", v, Ainv, v) - P[:, 1] * g_s**2, 3072) - 1 / 0.25
    L = np.bincount(reg, np.einsum("pi,pij,pj->p", v, Ainv, b0), 3072) + 1 / 0.25
    base = 0.5 * np.einsum("pi,pij,pj->p", b0, Ainv, b0).sum() - 0.5 * np.linalg.slogdet(Ap)[1].sum() + 0.5 * np.sum(np.log(P[:, 1:])) - 0.5 * np.sum(P[:, 2] * m0[:, 2]**2)
    # F(mu) = base + L mu + 1/2 Q mu^2 - 1/2 (mu-1)^2/0.25 absorbed in L,Q above; integral over mu: F(mu_hat) - 1/2 log(-Q)
    mu_hat = -L / Q
    return base + np.sum(L * mu_hat + 0.5 * Q * mu_hat**2) - 0.5 * np.sum(np.log(-Q)) - 0.5 * 3072 * 1 / 0.25 * 0  # constant dropped
def fit(grid, lp):
    i = lp.argmax(); s = slice(max(0, i - 3), min(len(grid), i + 4)); A = np.polyfit(grid[s], lp[s], 2); return -A[1] / (2 * A[0]), (np.sqrt(-0.5 / A[0]) if A[0] < 0 else np.nan)

print("== faint-level optimum with mu_r integrated analytically (realization data), bright regions and dust at the draw")
grid = np.linspace(-0.30, 0.30, 25); i0 = 12
for fn in [files[i] for i in (0, 116, 348, 464)]:
    st = state_of(fn); ss, sd_ = st["Synchrotron"], st["ThermalDust"]
    lp = np.array([marg_terms(ss["beta"] + d * faint, sd_["beta"], ss["als"], sd_["amu"], sd_["als"], maps) for d in grid]); m, s = fit(grid, lp)
    lp2 = np.array([marg_terms(ss["beta"] + d, sd_["beta"], ss["als"], sd_["amu"], sd_["als"], maps) for d in grid]); m2, s2 = fit(grid, lp2)
    print(f"{fn.split('/')[-1]}: faint mean {ss['beta'][~bright].mean():+.4f}: faint-shift optimum {m:+.3f} (sd {s:.3f}, gain {lp.max()-lp[i0]:.1f}) -> {ss['beta'][~bright].mean()+m:+.3f}; whole-sky shift optimum {m2:+.3f} (sd {s2:.3f}, gain {lp2.max()-lp2[i0]:.1f})")

print("\n== scalar index: exact 1-d marginal (mu_r integrated, beta prior N(-3.1, 0.3)) at the scalar run's dust state")
sf = sorted(glob.glob(f"{D}/chains_exact_syncscalar_mumap_nc/chains_compsep/chain01_iter*.h5"))[500:]
bd = []; amud = []; alsd = []; als = []; beta_sc = []
for fn in sf[::10]:
    with h5py.File(fn) as f:
        g = f["mcmc/sample_index_fields/state/"]; bd.append(g["ThermalDust/index_field/beta_map"][()]); amud.append(g["ThermalDust/amp_prior/mu"][()]); alsd.append(float(g["ThermalDust/amp_prior/log_s"][()])); als.append(float(g["Synchrotron/amp_prior/log_s"][()]))
for fn in sf:
    with h5py.File(fn) as f: beta_sc.append(float(np.asarray(f["comps/sync/sed/beta"]).ravel()[0]))
beta_sc = np.array(beta_sc)
bd_pm, amud_pm, lsd, lss = np.mean(bd, 0), np.mean(amud, 0), float(np.mean(alsd)), float(np.mean(als))
grid = np.linspace(-3.08, -2.82, 27)
def marginal_curve(d, bd_reg, amud_reg, lsd_, lss_):
    lp = np.array([marg_terms(np.full(3072, b), bd_reg, lss_, amud_reg, lsd_, d) - 0.5 * ((b + 3.1) / 0.3)**2 for b in grid])
    p = np.exp(lp - lp.max()); p /= np.trapezoid(p, grid); mean = np.trapezoid(grid * p, grid); sd = np.sqrt(np.trapezoid((grid - mean)**2 * p, grid)); return mean, sd, lp
mean, sd, lp = marginal_curve(maps, bd_pm, amud_pm, lsd, lss)
print(f"realization data, dust at scalar-run posterior mean: marginal mean {mean:+.4f} sd {sd:.4f}; chain {beta_sc.mean():+.4f} +- {beta_sc.std():.4f} (MC se 0.003); edge drops {lp[0]-lp.max():.1f} {lp[-1]-lp.max():.1f}")
# at a few individual dust draws, to see the conditional's dependence on the dust state
for k in (0, 20, 40, 60):
    mean_k, sd_k, _ = marginal_curve(maps, bd[k], amud[k], alsd[k], als[k]); print(f"  dust state of draw {k*10}: marginal mean {mean_k:+.4f} sd {sd_k:.4f}")
print("\n== pseudo-true scalar index: noise realizations of the exact marginal's mean (dust at scalar-run posterior mean; and at the truth with dust mu 1)")
rng = np.random.default_rng(1)
bd_truth16 = hp.ud_grade(t["beta_d"], 16)
means = []; means_t = []
for k in range(8):
    dsim = sky + rng.standard_normal(sky.shape) * t["rms"][:, None]
    means.append(marginal_curve(dsim, bd_pm, amud_pm, lsd, lss)[0]); means_t.append(marginal_curve(dsim, bd_truth16, np.ones(3072), lsd, lss)[0])
means = np.array(means); means_t = np.array(means_t)
print(f"dust at scalar-run state: pseudo-true {means.mean():+.4f} (scatter over realizations {means.std():.4f}, se {means.std()/np.sqrt(8):.4f}); values {np.round(means,3)}")
print(f"dust at truth (region means), dust mu 1: pseudo-true {means_t.mean():+.4f} (scatter {means_t.std():.4f}); values {np.round(means_t,3)}")
m0, _, _ = marginal_curve(sky, bd_pm, amud_pm, lsd, lss); m0t, _, _ = marginal_curve(sky, bd_truth16, np.ones(3072), lsd, lss)
print(f"noiseless data plugged in (no trace correction): {m0:+.4f} (scalar-run dust), {m0t:+.4f} (truth dust)")
