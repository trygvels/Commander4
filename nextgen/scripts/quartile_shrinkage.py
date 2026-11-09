#!/usr/bin/env python
"""Independent-region shrinkage prediction (D027) for an index-field chain, per S/N quartile.

pm_pred_r = (ref_r / fs_r^2 + mu / sigma^2) / (1 / fs_r^2 + 1 / sigma^2), with ref_r the region-fit
reference, fs_r the region Fisher sd at the true amplitudes (E027 refs), mu and sigma the
posterior-mean field mean and marginal sd of the chain (spatial correlation ignored). A quartile
whose chain-minus-reference offset matches pred-minus-reference is prior-dominated, not biased
by the likelihood.
usage: quartile_shrinkage.py --dir DIR --tag TAG --refs REFS.npz --component ThermalDust
       [--ref-key bd_trueT] [--q-key qd16] [--fisher-key fisher_sd_dust] [--burn 500]"""
import argparse, glob, numpy as np, h5py
ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--refs", required=True); ap.add_argument("--component", default="ThermalDust")
ap.add_argument("--ref-key"); ap.add_argument("--q-key"); ap.add_argument("--fisher-key")
ap.add_argument("--burn", type=int, default=500)
a = ap.parse_args()
dust = a.component == "ThermalDust"
ref_key = a.ref_key or ("bd_trueT" if dust else "bs_trueT")
q_key = a.q_key or ("qd16" if dust else "qs16")
f_key = a.fisher_key or ("fisher_sd_dust" if dust else "fisher_sd_sync")
r = np.load(a.refs); ref, q, fs = r[ref_key], r[q_key].astype(int), r[f_key]
files = sorted(glob.glob(f"{a.dir}/chains_{a.tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
bm, mu, ls = [], [], []
for fn in files:
    with h5py.File(fn) as f:
        g = f[f"mcmc/sample_index_fields/state/{a.component}"]
        bm.append(g["index_field/beta_map"][()]); mu.append(float(g["index_field/mu"][()]))
        ls.append(float(g["index_field/log_sigma"][()]))
bm = np.array(bm); pm = bm.mean(0); sd = bm.std(0)
mu_m = float(np.mean(mu)); sigma = float(np.mean(np.exp(ls)))
pred = (ref / fs**2 + mu_m / sigma**2) / (1 / fs**2 + 1 / sigma**2)
w = (1 / fs**2) / (1 / fs**2 + 1 / sigma**2)
print(f"{a.tag} {a.component}: {len(files)} post-burn draws; mu {mu_m:+.4f}, sigma {sigma:.4f}; "
      f"{ref_key} mean {ref.mean():+.4f}; per-quartile {ref_key} means "
      + " ".join(f"{ref[q == k].mean():+.4f}" for k in range(4)))
print("quartile  chain-ref  pred-ref  mu-ref  median w  median fs  median post sd  rms(chain-pred)")
for k in range(4):
    m = q == k
    print(f"  Q{k+1}    {np.mean(pm[m] - ref[m]):+.4f}   {np.mean(pred[m] - ref[m]):+.4f}  {mu_m - ref[m].mean():+.4f}"
          f"   {np.median(w[m]):.2f}     {np.median(fs[m]):.3f}     {np.median(sd[m]):.3f}        "
          f"{np.sqrt(np.mean((pm[m] - pred[m])**2)):.3f}")
print(f"all sky  chain-ref {np.mean(pm - ref):+.4f}  pred-ref {np.mean(pred - ref):+.4f}  "
      f"rms(chain-pred) {np.sqrt(np.mean((pm - pred)**2)):.3f}")
