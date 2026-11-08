"""Hyperparameter traces and ESS of an index-field chain: hyper_ess.py <dir> <tag> [burn]."""
import glob, sys, numpy as np, h5py
D, tag = sys.argv[1], sys.argv[2]; burn = int(sys.argv[3]) if len(sys.argv) > 3 else 500
files = sorted(glob.glob(f"{D}/chains_{tag}/chains_compsep/chain01_iter*.h5"))
post = files[burn:]
def ess(x):
    x = x - x.mean(); n = len(x); c = np.correlate(x, x, "full")[n-1:] / np.sum(x*x); s = 0.0
    for k in range(0, n-1, 2):
        p = c[k] + c[k+1]
        if p <= 0: break
        s += p
    return n / (2*s - 1) if s > 0.5 else n
keys = []
with h5py.File(post[0]) as f:
    st = f["mcmc/sample_index_fields/state"]
    for c in st:
        for sub in ("index_field/mu", "index_field/log_sigma", "index_field/u", "amp_prior/mu", "amp_prior/log_s"):
            if sub in st[c]: keys.append(f"{c}/{sub}")
vals = {k: [] for k in keys}; grads = []
with h5py.File(post[0]) as f:
    maps = {k: np.ndim(f["mcmc/sample_index_fields/state"][k][()]) > 0 for k in keys}
for fn in post:
    with h5py.File(fn) as f:
        st = f["mcmc/sample_index_fields/state"]
        for k in keys:
            v = np.asarray(st[k][()], dtype=float)
            vals[k].append(float(v) if v.ndim == 0 else float(v.mean()))
print(f"{tag}: {len(files)} files, {len(post)} post-burn")
for k in keys:
    x = np.array(vals[k]); print(f"  {k:40s} mean {x.mean():+.5f} sd {x.std():.5f} ESS {ess(x):7.1f}" + ("  (map mean)" if maps.get(k) else ""))
# region index ESS summary
with h5py.File(post[0]) as f:
    comps = [c for c in f["comps"] if "index_field" in f["comps"][c]]
for c in comps:
    b = np.array([h5py.File(fn)[f"comps/{c}/index_field/beta_map"][()] for fn in post])
    e = np.array([ess(b[:, i]) for i in range(b.shape[1])])
    print(f"  {c} beta_map: {b.shape[1]} regions, ESS min {e.min():.0f} 5% {np.percentile(e,5):.0f} median {np.median(e):.0f}")
