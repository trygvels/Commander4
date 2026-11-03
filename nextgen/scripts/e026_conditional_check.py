#!/usr/bin/env python
"""E026: Gibbs-consistency check of the index-field chain along the field-mean directions.

For K post-burn draws, the full target (Commander4 RegionIndexLikelihood with the non-centred
GMRF fields and the hier priors with mu maps) is evaluated along the field mean mu_c of one
component with every other coordinate held at the draw's value. The conditional CDF at the
draw's own mu_c must be uniform on (0, 1) over draws if the chain samples its target (the draw
is, conditionally on the rest, a draw from this 1-d density). Reports the CDF values, a KS
statistic, and the conditional optimum relative to the draw.
usage: e026_conditional_check.py --dir DIR --tag TAG --sync-template F --dust-template F [--k 20]
Negative controls (VALIDATION, D025): --index-prior-rms, --noise-scale or a wrong template corrupt
the target; the check must then fail.
"""
import argparse, glob, sys, time
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--sync-template", required=True); ap.add_argument("--dust-template", required=True)
ap.add_argument("--mu-nside", type=int, default=16); ap.add_argument("--cmb-sd", type=float, default=1e4)
ap.add_argument("--k", type=int, default=20); ap.add_argument("--burn", type=int, default=500); ap.add_argument("--nside-beta", type=int, default=16)
ap.add_argument("--half-width", type=float, nargs=2, default=[0.10, 0.015], help="scan half-widths (sync, dust)")
ap.add_argument("--points", type=int, default=11)
ap.add_argument("--mu-sd", type=float, default=0.5, help="sd of the N(1, sd^2) hyperprior on the template multipliers (E029: 0.1)")
ap.add_argument("--index-prior-rms", type=float, default=0.3, help="negative control: rms of the Gaussian prior on the field mean (chains use 0.3)")
ap.add_argument("--noise-scale", type=float, default=1.0, help="negative control: multiply the noise rms (chains use 1.0)")
ap.add_argument("--t-map", default=None, help="fixed per-region dust T map (FITS, E028 stage A); the target uses it instead of the scalar T_fixed")
a = ap.parse_args(); sys.argv = [sys.argv[0]]
import numpy as np, healpy as hp, h5py
from scipy.stats import kstest
from pixell.bunch import Bunch
from commander4.compsep import index_field as ifd
from commander4.compsep import index_field_sampler as ifs
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust

D = a.dir; NSIDE, NB = 64, a.nside_beta
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; freqs = np.array(bands, float)
maps = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands])
inv_var = 1.0 / (a.noise_scale * np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands]))**2
t = np.load(f"{D}/truth.npz"); T_fixed = float(t["T_fixed"])

def comp(cls, name, short, beta0, **sed):
    p = Bunch(polarization="I", shortname=short, spatially_varying_MM=False, Cl_prior_amplitude=None,
              lmax=191, units="uK_RJ", sample_spectral_index=True, spectral_index_nside=NB,
              spectral_index_prior=Bunch(type="gaussian", mean=beta0, rms=a.index_prior_rms),
              spectral_index_field=Bunch(type="gmrf", basis="graph", alpha=2.0, sigma_scale=0.3,
                                         rho_prior="loguniform"), beta=beta0, **sed)
    object.__setattr__(p, "_name", name)
    return cls(p, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True, eval_pol="I", comp_name=name)

cmb = comp(CMB, "CMB", "cmb", 0.0, nu_ref=1e-3); cmb.comp_params.sample_spectral_index = False
sync = comp(Synchrotron, "Synchrotron", "sync", -3.1, nu_ref=30.0)
dust = comp(ThermalDust, "ThermalDust", "dust", 1.6, T=T_fixed, nu_ref=545.0)
if a.t_map:
    tm = hp.read_map(a.t_map, field=0); dust.set_T_map(tm, hp.npix2nside(tm.size)); print(f"fixed T map {a.t_map}: mean {tm.mean():.3f} K")
cmb.amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=a.cmb_sd)
sync.amp_prior = ifd.AmpPrior("hier", template=np.maximum(hp.read_map(a.sync_template, field=0), 1e-3), mu_nside=a.mu_nside, mu_prior=(1.0, a.mu_sd))
dust.amp_prior = ifd.AmpPrior("hier", template=np.maximum(hp.read_map(a.dust_template, field=0), 1e-3), mu_nside=a.mu_nside, mu_prior=(1.0, a.mu_sd))
comps = [cmb, sync, dust]
specs, views = ifs.discover_index_specs(comps, None, centred=False)
names = ["Synchrotron", "ThermalDust"]
layout = ifd.StateLayout(specs, names, ifs.hier_component_dims(comps, names))
like = ifd.RegionIndexLikelihood.from_arrays(freqs, maps, inv_var, comps, layout, NSIDE, nthreads=8)
fields = {s.comp_name: s.field for s in specs}

files = sorted(glob.glob(f"{D}/chains_{a.tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
pick = files[:: max(1, len(files) // a.k)][: a.k]
print(f"{len(files)} post-burn draws, {len(pick)} used (every {max(1, len(files) // a.k)}th)")

def state_of(fn):
    st = {}
    with h5py.File(fn) as f:
        g = f["mcmc/sample_index_fields/state"]
        for c in names:
            fi = g[c]["index_field"]
            st[c] = dict(z=fi["z"][()], mu=float(fi["mu"][()]), ls=float(fi["log_sigma"][()]), u=float(fi["u"][()]),
                         amu=np.atleast_1d(g[c]["amp_prior/mu"][()]), als=float(g[c]["amp_prior/log_s"][()]))
    return st

def pack(st, shift):
    bl = {c: fields[c].pack(st[c]["z"], st[c]["mu"] + shift.get(c, 0.0), st[c]["ls"], st[c]["u"]) for c in names}
    return layout.pack(bl, {c: (st[c]["amu"], st[c]["als"]) for c in names})

t0 = time.time(); x0 = pack(state_of(pick[0]), {}); lp0, _ = like.logp_and_grad(x0)
print(f"one evaluation {time.time() - t0:.1f} s; first draw log target {lp0:.2f}")
cdf = {c: [] for c in names}; opt = {c: [] for c in names}; csd = {c: [] for c in names}
for i, fn in enumerate(pick):
    st = state_of(fn)
    for c, hw in zip(names, a.half_width):
        grid = np.linspace(-hw, hw, a.points)
        lp = np.array([like.logp_and_grad(pack(st, {c: d}))[0] for d in grid])
        # quadratic fit for the conditional optimum and sd; CDF at 0 by trapezoid on the grid
        A = np.polyfit(grid, lp, 2); mode = -A[1] / (2 * A[0]); sd = np.sqrt(-0.5 / A[0]) if A[0] < 0 else np.nan
        p = np.exp(lp - lp.max()); cum = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(grid))])
        F0 = np.interp(0.0, grid, cum / cum[-1])
        cdf[c].append(F0); opt[c].append(mode); csd[c].append(sd)
        print(f"draw {i:2d} {fn.split('/')[-1]} {c:12s}: mu {st[c]['mu']:+.4f}; conditional optimum at {mode:+.4f} (sd {sd:.4f}, "
              f"{mode / sd:+.1f} sd from the draw); CDF at the draw {F0:.3f}; log target range {lp.max() - lp.min():.1f}; "
              f"edge check lp(-hw)-max {lp[0] - lp.max():.1f} lp(+hw)-max {lp[-1] - lp.max():.1f}", flush=True)
for c in names:
    F = np.array(cdf[c]); ks = kstest(F, "uniform")
    print(f"{c}: CDF values mean {F.mean():.3f} (uniform: 0.5 +- {1 / np.sqrt(12 * len(F)):.3f}), fraction > 0.5 {np.mean(F > 0.5):.2f}, "
          f"KS D {ks.statistic:.3f} p {ks.pvalue:.3f}; conditional optimum minus draw: mean {np.mean(opt[c]):+.4f} "
          f"sd over draws {np.std(opt[c]):.4f}, in conditional sds mean {np.mean(np.array(opt[c]) / np.array(csd[c])):+.2f}; "
          f"conditional sd median {np.nanmedian(csd[c]):.4f}")
print("DONE")
