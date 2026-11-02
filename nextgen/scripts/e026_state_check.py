#!/usr/bin/env python
"""E026: the full target's log density (Commander4 RegionIndexLikelihood, non-centred GMRF fields,
hier priors with mu maps) at the chain's final state against (a) the same state with the index
fields moved to the noise-free region-fit reference, (b) the chain state with the sync field mean
shifted by a given amount, (c) the chain state with the dust mu map set to 1. Says whether the
chain's offsets are where the target puts its mass."""
import argparse, glob, sys
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--refs", required=True); ap.add_argument("--sync-template", required=True); ap.add_argument("--dust-template", required=True)
ap.add_argument("--mu-nside", type=int, default=16); ap.add_argument("--cmb-sd", type=float, default=1e4)
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
like = ifd.RegionIndexLikelihood.from_arrays(freqs, maps, inv_var, comps, layout, NSIDE, nthreads=4)
fields = {s.comp_name: s.field for s in specs}

files = sorted(glob.glob(f"{D}/chains_{a.tag}/chains_compsep/chain01_iter*.h5"))
post = files[500:]
acc = {c: {"beta": [], "mu": [], "ls": [], "fmu": [], "fls": [], "fu": []} for c in names}
for fn in post:
    with h5py.File(fn) as f:
        st = f["mcmc/sample_index_fields/state"]
        for c in names:
            g = st[c]["index_field"]
            acc[c]["beta"].append(g["beta_map"][()]); acc[c]["fmu"].append(float(g["mu"][()]))
            acc[c]["fls"].append(float(g["log_sigma"][()])); acc[c]["fu"].append(float(g["u"][()]))
            acc[c]["mu"].append(np.atleast_1d(st[c]["amp_prior/mu"][()])); acc[c]["ls"].append(float(st[c]["amp_prior/log_s"][()]))
mean = {c: {k: np.mean(np.array(v), axis=0) for k, v in acc[c].items()} for c in names}
beta_chain = {c: mean[c]["beta"] for c in names}
hyp = {c: (mean[c]["mu"] if a.mu_nside else float(mean[c]["mu"][0]), float(mean[c]["ls"])) for c in names}
print(f"posterior means over {len(post)} draws: field means {mean['Synchrotron']['fmu']:.4f} {mean['ThermalDust']['fmu']:.4f}, "
      f"sigma {np.exp(mean['Synchrotron']['fls']):.4f} {np.exp(mean['ThermalDust']['fls']):.4f}, beta-map means "
      f"{beta_chain['Synchrotron'].mean():.4f} {beta_chain['ThermalDust'].mean():.4f}, mu_d map mean {np.mean(mean['ThermalDust']['mu']):.4f}")

def with_betas(betas, hyp_use=hyp):
    """Posterior-mean field hyperparameters (log sigma, u) and amplitude hyperparameters; the
    field mean set to the mean of `betas` and the coefficients solved for the rest."""
    bl = {}
    for c in names:
        fld = fields[c]; ls, u = mean[c]["fls"], mean[c]["fu"]
        sigma = np.exp(ls); rho = float(fld.transform.rho(u))
        mu_f = float(np.mean(betas[c])); resid = betas[c] - mu_f
        V = fld.basis.V; w = fld.field.weights(rho)[0]
        z_new = (V.T @ resid) / (sigma * w)
        bl[c] = fld.pack(z_new, mu_f, ls, u)
    return layout.pack(bl, hyp_use)

x_chain = with_betas(beta_chain)
lp_chain, _ = like.logp_and_grad(x_chain)
print(f"posterior-mean state: log target {lp_chain:.2f}")
refs = np.load(a.refs)
ref = {"Synchrotron": refs["bs_fixT"], "ThermalDust": refs["bd_fixT"]}
for label, betas in (("region-fit reference", ref),
                     ("reference, dust only (sync at chain)", {"Synchrotron": beta_chain["Synchrotron"], "ThermalDust": ref["ThermalDust"]}),
                     ("reference, sync only (dust at chain)", {"Synchrotron": ref["Synchrotron"], "ThermalDust": beta_chain["ThermalDust"]})):
    x = with_betas(betas); lp, _ = like.logp_and_grad(x)
    print(f"{label:42s}: log target {lp:.2f}  (chain minus this: {lp_chain - lp:+.2f})  field means {np.mean(betas['Synchrotron']):.4f} {np.mean(betas['ThermalDust']):.4f}")
for shift in (-0.13, -0.065, +0.065):
    b = {"Synchrotron": beta_chain["Synchrotron"] + shift, "ThermalDust": beta_chain["ThermalDust"]}
    x = with_betas(b); lp, _ = like.logp_and_grad(x)
    print(f"chain with sync field shifted by {shift:+.3f}: log target {lp:.2f} (chain minus this {lp_chain - lp:+.2f})")
for shift in (-0.04, -0.02, +0.02):
    b = {"Synchrotron": beta_chain["Synchrotron"], "ThermalDust": beta_chain["ThermalDust"] + shift}
    x = with_betas(b); lp, _ = like.logp_and_grad(x)
    print(f"chain with dust field shifted by {shift:+.3f}: log target {lp:.2f} (chain minus this {lp_chain - lp:+.2f})")
if a.mu_nside:
    hyp1 = dict(hyp); hyp1["ThermalDust"] = (np.ones_like(hyp["ThermalDust"][0]), hyp["ThermalDust"][1])
    x = with_betas(beta_chain, hyp1); lp, _ = like.logp_and_grad(x)
    print(f"chain with dust mu map set to 1: log target {lp:.2f} (chain minus this {lp_chain - lp:+.2f}); chain mu_d map mean {np.mean(hyp['ThermalDust'][0]):.4f}")
    b = {"Synchrotron": beta_chain["Synchrotron"], "ThermalDust": beta_chain["ThermalDust"] - 0.04}
    x = with_betas(b, hyp1); lp, _ = like.logp_and_grad(x)
    print(f"chain with dust mu map 1 and dust field -0.04: log target {lp:.2f} (chain minus this {lp_chain - lp:+.2f})")
