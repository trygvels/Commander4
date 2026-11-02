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
pick = files[::35][:20]
def pack(st, zs, mus, lss, us):
    bl = {"Synchrotron": F.pack(zs, mus, lss, us), "ThermalDust": fld["ThermalDust"].pack(st["ThermalDust"]["z"], st["ThermalDust"]["mu"], st["ThermalDust"]["ls"], st["ThermalDust"]["u"])}
    return layout.pack(bl, {c: (st[c]["amu"], st[c]["als"]) for c in names})
def pit(grid, lp):
    p = np.exp(lp - lp.max()); cum = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(grid))]); return float(np.interp(0.0, grid, cum / cum[-1]))
dirs = {"z_0": 0, "z_50": 50, "z_1000": 1000, "u": None}
res = {k: [] for k in dirs}
for i, fn in enumerate(pick):
    st = state_of(fn); s = st["Synchrotron"]; z0, mu0, ls0, u0 = s["z"], s["mu"], s["ls"], s["u"]
    line = []
    for name, k in dirs.items():
        def f(d):
            if k is None: return like.logp_and_grad(pack(st, z0, mu0, ls0, u0 + d))[0]
            z = z0.copy(); z[k] += d; return like.logp_and_grad(pack(st, z, mu0, ls0, u0))[0]
        g = np.linspace(-4, 4, 9); lp = np.array([f(d) for d in g])          # coarse
        A = np.polyfit(g, lp, 2); sd = np.sqrt(-0.5 / A[0]) if A[0] < 0 else 1.5; ctr = np.clip(-A[1] / (2 * A[0]) if A[0] < 0 else 0.0, -3, 3)
        hw = min(4.5 * sd, 6.0); g2 = np.linspace(ctr - hw, ctr + hw, 13)
        if not (g2[0] < 0 < g2[-1]): g2 = np.linspace(min(g2[0], -0.5 * hw), max(g2[-1], 0.5 * hw), 13)
        lp2 = np.array([f(d) for d in g2]); c = pit(g2, lp2); res[name].append(c)
        line.append(f"{name} PIT {c:.3f} (sd {sd:.2f}, edges {lp2[0]-lp2.max():.1f} {lp2[-1]-lp2.max():.1f})")
    print(f"draw {i:2d}: " + "; ".join(line), flush=True)
for k, v in res.items():
    v = np.array(v); ks = kstest(v, "uniform"); print(f"{k}: PIT mean {v.mean():.3f} +- 0.065, KS D {ks.statistic:.3f} p {ks.pvalue:.3f}")
