import sys; sys.argv = [sys.argv[0]]
exec(open("marg.py").read().split('print("== faint-level optimum')[0])
rng = np.random.default_rng(3)
c = np.load("chain_summary.npz"); bm = c["bm_mean"]
r = np.load("/Users/trygvels/code/Commander4/nextgen/experiments/E023/scan_n1_refs.npz"); ref_s = r["bs_fixT"]
st = state_of(files[0]); ss, sd_ = st["Synchrotron"], st["ThermalDust"]
amu_d = np.mean([state_of(f)["ThermalDust"]["amu"] for f in files[::20]], 0); bd_pm = np.mean([state_of(f)["ThermalDust"]["beta"] for f in files[::20]], 0)
ls_s, ls_d = np.log(0.0342), np.log(0.0219)
grid = np.linspace(-0.30, 0.30, 25); i0 = 12
fm = ~bright
def wlevel(b): return float(np.sum(W[fm] * b[fm]) / W[fm].sum())
print(f"faint set: unweighted / Fisher-weighted level: chain mean {bm[fm].mean():+.4f} / {wlevel(bm):+.4f}; reference {ref_s[fm].mean():+.4f} / {wlevel(ref_s):+.4f}; draw0 {ss['beta'][fm].mean():+.4f} / {wlevel(ss['beta']):+.4f}")
def block(bs0, bd0, mud, label, expected_one=False):
    if expected_one:
        lp = np.array([data_term(bs0 + d * faint, bd0, np.ones(3072), ls_s, mud, ls_d, expected=True) for d in grid])
    else:
        lp = np.array([marg_terms(bs0 + d * faint, bd0, ls_s, mud, ls_d, maps) for d in grid])
    m, s = fit(grid, lp)
    print(f"  {label:70s}: optimum {m:+.3f} (sd {s:.3f}, gain {lp.max()-lp[i0]:.1f}) -> faint level {bs0[fm].mean()+m:+.3f}")
dev_chain = np.where(fm, bm - bm[fm].mean(), 0.0)                     # posterior-mean deviation pattern in faint regions
dev_draw = np.where(fm, ss["beta"] - ss["beta"][fm].mean(), 0.0)
shuf = dev_draw.copy(); shuf[fm] = rng.permutation(shuf[fm])
iid = np.where(fm, rng.standard_normal(3072) * 0.15, 0.0)
print("realization data, sync multipliers integrated, dust at chain state, bright regions at the chain mean:")
B = np.where(bright, bm, 0.0)
for bs0, label in ((bm, "posterior-mean field"), (B + np.where(fm, ref_s, 0), "faint = reference pattern"),
                   (B + np.where(fm, ref_s[fm].mean(), 0) + dev_chain, "faint = ref level + chain posterior-mean deviation pattern"),
                   (B + np.where(fm, ref_s[fm].mean(), 0) + dev_draw, "faint = ref level + draw-0 deviation pattern (sd %.3f)" % dev_draw[fm].std()),
                   (B + np.where(fm, ref_s[fm].mean(), 0) + shuf, "faint = ref level + draw-0 pattern shuffled among faint regions"),
                   (B + np.where(fm, ref_s[fm].mean(), 0) + iid, "faint = ref level + iid N(0, 0.15) pattern"),
                   (B + np.where(fm, ss["beta"][fm].mean(), 0), "faint = constant at draw-0 faint mean"),
                   (np.where(bright, ss["beta"], 0) + np.where(fm, ss["beta"], 0), "draw 0 itself")):
    block(bs0, bd_pm, amu_d, label)
print("noiseless data, multipliers = 1 (pure signal-model response to the pattern):")
for bs0, label in ((B + np.where(fm, ref_s, 0), "faint = reference pattern"), (B + np.where(fm, ref_s[fm].mean(), 0) + dev_draw, "faint = ref level + draw-0 deviation pattern"),
                   (B + np.where(fm, ref_s[fm].mean(), 0) + iid, "faint = ref level + iid N(0, 0.15) pattern"), (B + np.where(fm, ref_s[fm].mean(), 0) + 2 * iid, "faint = ref level + iid N(0, 0.30) pattern")):
    block(bs0, bd_pm, amu_d, label, expected_one=True)
