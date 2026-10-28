"""Which model/sim mismatch produces the faint-region dust pull: CMB prior width, or the template
construction (g = a/r in the sim versus a = g r in the model's story)? Resimulate in numpy."""
import numpy as np, healpy as hp, sys
src = open(sys.argv[1]).read(); sys.argv = [sys.argv[0], "gmrf16"]; exec(src.split("bs64, bd64 = ")[0])   # reuse data, SEDs and marg() from the scan script
bs64, bd64 = hp.ud_grade(bs16, 64), hp.ud_grade(bd16, 64)
ts64, td64 = hp.ud_grade(hp.ud_grade(t["beta_s"], 16), 64), hp.ud_grade(hp.ud_grade(t["beta_d"], 16), 64)
fsd16 = t["fisher_sd_dust_16_profiled"]; qd = np.argsort(fsd16)[-768:]
seld = np.isin(hp.ud_grade(np.arange(3072.0), 64).astype(int), qd)
cmb = t["cmb"]; print("true CMB rms all / faint cap (uK_CMB):", cmb.std().round(2), cmb[seld].std().round(2))
a_s, a_d = t["a_sync"], t["a_dust"]
print("a_s/g_s mean, sd:", (a_s/gs).mean().round(4), (a_s/gs).std().round(4), " a_d/g_d:", (a_d/gd).mean().round(4), (a_d/gd).std().round(4))
deltas = np.arange(-0.10, 0.251, 0.025)
def scan(label, dat, state, cmb_sd=100.0):
    global d; d0 = d; d = dat
    bs, bd, mus, ss, mud, sd = state; out = []
    for dl in deltas: out.append(marg(bs, bd + dl*seld, mus, ss, mud, sd, seld, cmb_sd=cmb_sd))
    d = d0; out = np.array(out); tot = out.sum(1); i0 = np.argmin(abs(deltas)); tot -= tot[i0]
    # quadratic fit of the peak
    j = np.argmax(tot); j = min(max(j, 1), len(deltas)-2); c = np.polyfit(deltas[j-1:j+2], tot[j-1:j+2], 2); pk = -c[1]/(2*c[0])
    print(f"{label:70s} argmax {deltas[np.argmax(tot)]:+.3f} (quadratic peak {pk:+.3f}), gain {tot.max():8.1f}; chi2 part at +0.125 {out[np.argmin(abs(deltas-0.125)),0]-out[i0,0]:8.1f}, logdet part {out[np.argmin(abs(deltas-0.125)),1]-out[i0,1]:8.1f}")
truth = (ts64, td64, 0.9990, 0.2986, 0.9976, 0.2988)
scan("real sim data, cmb_sd 100 (the model)", d, truth, 100.0)
scan("real sim data, cmb_sd 84.1 (true per-pixel rms)", d, truth, 84.1)
scan("real sim data, cmb_sd 1e4 (near flat)", d, truth, 1e4)
scan("real sim data, cmb_sd 100, s_d 0.2988 but dust template = true a_d (r = 1)", d, truth, 100.0) if False else None
rng = np.random.default_rng(11)
Ms, Md = sync_sed(ts64), dust_sed(td64)
def simulate(a_s_, a_d_, seed):
    r_ = np.random.default_rng(seed)
    return (cmb_mix[None, :]*cmb[:, None] + Ms*a_s_[:, None] + Md*a_d_[:, None] + r_.standard_normal(d.shape)*rms)
scan("resim: same a (g = a/r story), fresh noise", simulate(a_s, a_d, 1), truth, 100.0)
scan("resim: same a, fresh noise, seed 2", simulate(a_s, a_d, 2), truth, 100.0)
rs, rd = 1 + 0.3*rng.standard_normal(len(a_s)), 1 + 0.3*rng.standard_normal(len(a_d))
dd = simulate(gs*rs, gd*rd, 3)
tr2 = (ts64, td64, rs.mean(), rs.std(), rd.mean(), rd.std())
scan("resim: a = g r (model story, templates g unchanged), fresh noise", dd, tr2, 100.0)
scan("resim: a = g r, cmb_sd 84.1", dd, tr2, 84.1)
# no CMB at all
dd2 = dd - cmb_mix[None, :]*cmb[:, None]
scan("resim: a = g r, CMB removed from data, cmb_sd 100", dd2, tr2, 100.0)
scan("resim: a = g r, CMB removed, cmb_sd 1 (CMB off)", dd2, tr2, 1.0)
# hyperparameter scans on the model-story data
allp = np.ones(len(d), bool)
for lab, dat, st in (("model-story data", dd, tr2),):
    d0 = d; d = dat
    grid = np.arange(0.95, 1.031, 0.005)
    v = np.array([sum(marg(st[0], st[1], mu, st[3], st[4], st[5], allp)) for mu in grid]); print(f"mu_s scan, {lab}: argmax {grid[np.argmax(v)]:.3f} (realized {rs.mean():.4f})")
    v = np.array([sum(marg(st[0], st[1], st[2], st[3], mu, st[5], allp)) for mu in grid]); print(f"mu_d scan, {lab}: argmax {grid[np.argmax(v)]:.3f} (realized {rd.mean():.4f})")
    d = d0
