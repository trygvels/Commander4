#!/usr/bin/env python
"""Dust index field of a chain against the true-T (bd_trueT) and fixed-T (bd_fixT) region-fit
references per dust S/N quartile (E027 labels), plus the dust mu-map mean.
usage: e028_dust_vs_refs.py --dir DIR --refs REFS.npz --tags TAG [TAG ...] [--burn 500]"""
import argparse, glob, numpy as np, h5py
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--refs", required=True)
ap.add_argument("--tags", nargs="+", required=True); ap.add_argument("--burn", type=int, default=500)
a = ap.parse_args(); r = np.load(a.refs); q = r["qd16"] if "qd16" in r.files else r["qd"]
for tag in a.tags:
    files = sorted(glob.glob(f"{a.dir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
    bd, mu = [], []
    for fn in files:
        with h5py.File(fn) as f:
            g = f["mcmc/sample_index_fields/state/ThermalDust"]
            bd.append(g["index_field/beta_map"][()]); mu.append(g["amp_prior/mu"][()])
    bd = np.array(bd); mu = np.array(mu); pm = bd.mean(0); sd = bd.std(0)
    print(f"{tag}: {len(files)} post-burn draws; dust mu map mean {mu.mean():.4f} +- {mu.mean(1).std():.4f}; field mean {pm.mean():+.4f}")
    for name in ("bd_trueT", "bd_fixT"):
        ref = r[name]; z = (pm - ref) / sd
        print(f"  against {name}: " + "; ".join(f"Q{k+1} z median {np.median(z[q==k]):+.2f} mean(pm-ref) {np.mean((pm-ref)[q==k]):+.3f} rmse {np.sqrt(np.mean((pm-ref)[q==k]**2)):.3f}" for k in range(4)))
