"""PIT along two more directions of the full Commander4 target: log sigma_s, and the faint-level
direction (shift of the faint regions' beta with the bright regions fixed, z and mu adjusted)."""
import os; os.environ["OMP_NUM_THREADS"] = "1"
import sys; sys.argv = [sys.argv[0]]
exec(open("decomp.py").read().split("# ---- verification")[0])
from commander4.compsep import index_field as ifd, index_field_sampler as ifs
from scipy.stats import kstest
cmb.amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=1e4)
sync.amp_prior = ifd.AmpPrior("hier", template=g_s, mu_nside=16); dust.amp_prior = ifd.AmpPrior("hier", template=g_d, mu_nside=16)
comps = [cmb, sync, dust]; specs, views = ifs.discover_index_specs(comps, None, centred=False); names = ["Synchrotron", "ThermalDust"]
layout = ifd.StateLayout(specs, names, ifs.hier_component_dims(comps, names))
like = ifd.RegionIndexLikelihood.from_arrays(freqs, maps, inv_var, comps, layout, NSIDE, nthreads=2)
fld = {s.comp_name: s.field for s in specs}; F = fld["Synchrotron"]
bs_t, a_s, nu, rms = t["beta_s"], t["a_sync"], t["nu"], t["rms"]; x_ = np.log(nu / 30.0)
w_known = ((a_s[:, None] * (nu / 30.0)[None, :] ** bs_t[:, None] * x_[None, :])**2 / rms[None, :]**2).sum(1)
W = np.bincount(reg, w_known, 3072); order = np.argsort(-W); bright = np.zeros(3072, bool); bright[order[:768]] = True
faint_ind = (~bright).astype(float)
pick = files[::35][:20]
def pack(st, zs, mus, lss, us):
    bl = {"Synchrotron": F.pack(zs, mus, lss, us), "ThermalDust": fld["ThermalDust"].pack(st["ThermalDust"]["z"], st["ThermalDust"]["mu"], st["ThermalDust"]["ls"], st["ThermalDust"]["u"])}
    return layout.pack(bl, {c: (st[c]["amu"], st[c]["als"]) for c in names})
def pit(grid, lp):
    p = np.exp(lp - lp.max()); cum = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(grid))]); return float(np.interp(0.0, grid, cum / cum[-1])), lp[0] - lp.max(), lp[-1] - lp.max()
res = {"log_sigma": [], "faint_level": []}
for i, fn in enumerate(pick):
    st = state_of(fn); s = st["Synchrotron"]; z0, mu0, ls0, u0 = s["z"], s["mu"], s["ls"], s["u"]
    sigma = np.exp(ls0); rho = float(F.transform.rho(u0)); w, _ = F.field.weights(rho); V = F.basis.V
    # (i) log sigma: grid +-0.9 around the draw (sigma posterior sd 0.047 / 0.15 ~ 0.3 in log)
    g1 = np.linspace(-0.9, 0.9, 13); lp1 = np.array([like.logp_and_grad(pack(st, z0, mu0, ls0 + d, u0))[0] for d in g1])
    c1, e1a, e1b = pit(g1, lp1); res["log_sigma"].append(c1)
    # (ii) faint-level: beta' = beta + d * 1_faint; mu' = mu + d * mean(1_faint); z' = z + d * V^T (1_faint - mean) / (sigma w)
    pat = faint_ind - faint_ind.mean(); dz = (V.T @ pat) / (sigma * w)
    g2 = np.linspace(-0.25, 0.25, 13); lp2 = np.array([like.logp_and_grad(pack(st, z0 + d * dz, mu0 + d * faint_ind.mean(), ls0, u0))[0] for d in g2])
    c2, e2a, e2b = pit(g2, lp2); res["faint_level"].append(c2)
    # check the mapping reproduces the intended beta
    b_chk = F.beta(F.pack(z0 + 0.1 * dz, mu0 + 0.1 * faint_ind.mean(), ls0, u0)) - F.beta(F.pack(z0, mu0, ls0, u0))
    print(f"draw {i:2d}: log sigma PIT {c1:.3f} (edges {e1a:.1f} {e1b:.1f}); faint-level PIT {c2:.3f} (edges {e2a:.1f} {e2b:.1f}; range {lp2.max()-lp2.min():.1f}); mapping check {np.abs(b_chk - 0.1*faint_ind).max():.1e}", flush=True)
for k, v in res.items():
    v = np.array(v); ks = kstest(v, "uniform"); print(f"{k}: PIT mean {v.mean():.3f} +- {1/np.sqrt(12*len(v)):.3f}, KS D {ks.statistic:.3f} p {ks.pvalue:.3f}, values {np.round(np.sort(v),2)}")
