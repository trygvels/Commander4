#!/usr/bin/env python
"""E026: does the block-shift optimum of the sync field depend on the field's roughness?

The exact-target scan at the chain's posterior-mean state (a smooth field: the coefficients z
are shrunk by averaging) preferred a whole-sky sync shift of -0.12 to -0.18, while the
Gibbs-consistency check at typical draws (field sd 0.15 around the mean) finds the conditional
optimum at the draw. Here, for a few draws, the non-centred coefficients z are scaled by a
factor (1 = the draw, 0 = a constant field at mu) and the block-shift optimum along mu is
located by a quadratic fit to an 11-point scan, with the draw's multipliers and dust state;
then the same at the posterior-mean state. Reuses the Commander4 target of
e026_conditional_check.py.
usage: e026_roughness_check.py --dir DIR --tag TAG --sync-template F --dust-template F
"""
import argparse, glob, sys
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--sync-template", required=True); ap.add_argument("--dust-template", required=True)
ap.add_argument("--mu-nside", type=int, default=16); ap.add_argument("--cmb-sd", type=float, default=1e4)
ap.add_argument("--k", type=int, default=6); ap.add_argument("--burn", type=int, default=500)
ap.add_argument("--scales", type=float, nargs="+", default=[1.0, 0.5, 0.25, 0.0])
a = ap.parse_args(); sys.argv = [sys.argv[0]]
import numpy as np, healpy as hp, h5py
from pixell.bunch import Bunch
from commander4.compsep import index_field as ifd
from commander4.compsep import index_field_sampler as ifs
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust

D = a.dir; NSIDE, NB = 64, 16
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; freqs = np.array(bands, float)
maps = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands])
inv_var = 1.0 / np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands])**2
t = np.load(f"{D}/truth.npz"); T_fixed = float(t["T_fixed"])

def comp(cls, name, short, beta0, **sed):
    p = Bunch(polarization="I", shortname=short, spatially_varying_MM=False, Cl_prior_amplitude=None,
              lmax=191, units="uK_RJ", sample_spectral_index=True, spectral_index_nside=NB,
              spectral_index_prior=Bunch(type="gaussian", mean=beta0, rms=0.3),
              spectral_index_field=Bunch(type="gmrf", basis="graph", alpha=2.0, sigma_scale=0.3,
                                         rho_prior="loguniform"), beta=beta0, **sed)
    object.__setattr__(p, "_name", name)
    return cls(p, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True, eval_pol="I", comp_name=name)

cmb = comp(CMB, "CMB", "cmb", 0.0, nu_ref=1e-3); cmb.comp_params.sample_spectral_index = False
sync = comp(Synchrotron, "Synchrotron", "sync", -3.1, nu_ref=30.0)
dust = comp(ThermalDust, "ThermalDust", "dust", 1.6, T=T_fixed, nu_ref=545.0)
cmb.amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=a.cmb_sd)
sync.amp_prior = ifd.AmpPrior("hier", template=np.maximum(hp.read_map(a.sync_template, field=0), 1e-3), mu_nside=a.mu_nside)
dust.amp_prior = ifd.AmpPrior("hier", template=np.maximum(hp.read_map(a.dust_template, field=0), 1e-3), mu_nside=a.mu_nside)
comps = [cmb, sync, dust]
specs, views = ifs.discover_index_specs(comps, None, centred=False)
names = ["Synchrotron", "ThermalDust"]
layout = ifd.StateLayout(specs, names, ifs.hier_component_dims(comps, names))
like = ifd.RegionIndexLikelihood.from_arrays(freqs, maps, inv_var, comps, layout, NSIDE, nthreads=8)
fields = {s.comp_name: s.field for s in specs}
files = sorted(glob.glob(f"{D}/chains_{a.tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
pick = files[:: max(1, len(files) // a.k)][: a.k]

def state_of(fn):
    st = {}
    with h5py.File(fn) as f:
        g = f["mcmc/sample_index_fields/state"]
        for c in names:
            fi = g[c]["index_field"]
            st[c] = dict(z=fi["z"][()], mu=float(fi["mu"][()]), ls=float(fi["log_sigma"][()]), u=float(fi["u"][()]),
                         amu=np.atleast_1d(g[c]["amp_prior/mu"][()]), als=float(g[c]["amp_prior/log_s"][()]),
                         beta=fi["beta_map"][()])
    return st

def pack(st, shift=0.0, zscale=1.0):
    bl = {}
    for c in names:
        s = st[c]; z = s["z"] * (zscale if c == "Synchrotron" else 1.0)
        bl[c] = fields[c].pack(z, s["mu"] + (shift if c == "Synchrotron" else 0.0), s["ls"], s["u"])
    return layout.pack(bl, {c: (st[c]["amu"], st[c]["als"]) for c in names})

def optimum(st, zscale):
    grid = np.linspace(-0.25, 0.10, 15)
    lp = np.array([like.logp_and_grad(pack(st, d, zscale))[0] for d in grid])
    i = lp.argmax(); sel = slice(max(0, i - 3), min(len(grid), i + 4))
    A = np.polyfit(grid[sel], lp[sel], 2)
    return -A[1] / (2 * A[0]), np.sqrt(-0.5 / A[0]) if A[0] < 0 else np.nan, lp.max() - lp[np.argmin(np.abs(grid))]

for i, fn in enumerate(pick):
    st = state_of(fn)
    sd_field = np.std(st["Synchrotron"]["beta"])
    out = []
    for zs in a.scales:
        m, sd, gain = optimum(st, zs); out.append(f"z x {zs:.2f}: optimum {m:+.3f} (sd {sd:.3f}, gain over 0 {gain:.1f})")
    print(f"draw {i} {fn.split('/')[-1]}: mu {st['Synchrotron']['mu']:+.4f}, field sd {sd_field:.3f}, sigma {np.exp(st['Synchrotron']['ls']):.3f}; "
          + "; ".join(out), flush=True)
# posterior-mean state: mean of z, mu, log sigma, u, multipliers over all post-burn draws
acc = None
for fn in files:
    st = state_of(fn)
    if acc is None:
        acc = {c: {k: np.zeros_like(np.asarray(v, float)) for k, v in st[c].items()} for c in names}
    for c in names:
        for k, v in st[c].items():
            acc[c][k] += np.asarray(v, float) / len(files)
pm = {c: {k: (float(v) if np.ndim(v) == 0 else v) for k, v in acc[c].items()} for c in names}
for zs in a.scales:
    m, sd, gain = optimum(pm, zs)
    print(f"posterior-mean state (z, mu, hyperparameters, multipliers averaged): z x {zs:.2f}: optimum {m:+.3f} (sd {sd:.3f}, gain over 0 {gain:.1f}); "
          f"field sd at this z scale {np.std(fields['Synchrotron'].unpack(fields['Synchrotron'].pack(pm['Synchrotron']['z'] * zs, pm['Synchrotron']['mu'], pm['Synchrotron']['ls'], pm['Synchrotron']['u']))[0]) if hasattr(fields['Synchrotron'], 'unpack') else float('nan'):.3f}")
print("DONE")
