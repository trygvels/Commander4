"""E026: where does the +0.13 sync aggregate of gmrf16_exact_mumap_nc come from?

Compares (1) the single-sync-index run (exact_syncscalar_mumap_nc) with its pseudo-true value
(the Fisher-weighted truth mean), (2) the posterior shape of the field run's mean (skewness),
(3) the field run's data-dominated regions (top-N by region Fisher weight, amplitudes known)
against the noise-free region-fit reference, and (4) the Fisher-weighted level of the chain's
field against the Fisher-weighted truth and reference, per template-S/N quartile.
usage: e026_sync_level_check.py --dir DIR --refs nextgen/experiments/E023/scan_n1_refs.npz [--burn 500]
"""
import argparse, glob
import numpy as np, h5py, healpy as hp
from scipy.stats import skew

S = "mcmc/sample_index_fields/state/Synchrotron/"


def load(path, keys, burn):
    files = sorted(glob.glob(f"{path}/chains_compsep/*.h5"))[burn:]
    out = {k: [] for k in keys}
    for fn in files:
        with h5py.File(fn) as f:
            for k in keys:
                out[k].append(np.asarray(f[k]))
    return {k: np.array(v) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--refs", required=True)
    ap.add_argument("--burn", type=int, default=500)
    ap.add_argument("--field-run", default="chains_gmrf16_exact_mumap_nc")
    ap.add_argument("--scalar-run", default="chains_exact_syncscalar_mumap_nc")
    a = ap.parse_args()
    t = np.load(f"{a.dir}/truth.npz"); r = np.load(a.refs)
    bs, amp, nu, rms = t["beta_s"], t["a_sync"], t["nu"], t["rms"]
    x = np.log(nu / 30.0)
    w = ((amp[:, None] * (nu / 30.0)[None, :] ** bs[:, None] * x[None, :]) ** 2 / rms[None, :] ** 2).sum(1)
    nside_b = hp.npix2nside(r["bs_fixT"].size)
    reg = hp.ud_grade(np.arange(r["bs_fixT"].size, dtype=float), int(t["nside"])).astype(int)
    W = np.bincount(reg, w, r["bs_fixT"].size); fsd = 1 / np.sqrt(W)
    ref, q, truth = r["bs_fixT"], r["qs16"], hp.ud_grade(np.asarray(bs, float), nside_b)
    wm = lambda wt, b: float((wt * b).sum() / wt.sum())
    print(f"truth: unweighted mean {bs.mean():+.4f}; Fisher-weighted (amplitudes known) {wm(w, bs):+.4f}; "
          f"per quartile Fisher-weighted " + " ".join(f"Q{k+1} {wm(W[q == k], truth[q == k]):+.4f}" for k in range(4))
          + "; Fisher-weight share per quartile " + " ".join(f"{W[q == k].sum() / W.sum():.3f}" for k in range(4)))
    sc = load(f"{a.dir}/{a.scalar_run}", ["comps/sync/sed/beta"], a.burn)["comps/sync/sed/beta"].ravel()
    print(f"scalar run: beta_s mean {sc.mean():+.4f} median {np.median(sc):+.4f} sd {sc.std():.4f} skew {skew(sc):+.2f} "
          f"(n {sc.size}); minus Fisher-weighted truth {sc.mean() - wm(w, bs):+.4f} ({(sc.mean() - wm(w, bs)) / sc.std():+.1f} sd)")
    d = load(f"{a.dir}/{a.field_run}", [S + "index_field/mu", S + "index_field/beta_map"], a.burn)
    mu, bm = d[S + "index_field/mu"], d[S + "index_field/beta_map"]
    pm, sd = bm.mean(0), bm.std(0)
    print(f"field run: mu mean {mu.mean():+.4f} median {np.median(mu):+.4f} sd {mu.std():.4f} skew {skew(mu):+.2f}; "
          f"field mean over regions {bm.mean(1).mean():+.4f}; per-region skew of beta_r mean {skew(bm, axis=0).mean():+.3f}")
    print(f"field run: Fisher-weighted level of the posterior-mean field {wm(W, pm):+.4f}; reference {wm(W, ref):+.4f}; "
          f"truth {wm(W, truth):+.4f}")
    for k in range(4):
        m = q == k
        print(f"  Q{k+1}: posterior mean minus reference {np.mean(pm[m] - ref[m]):+.3f} (median {np.median(pm[m] - ref[m]):+.3f}); "
              f"quartile mean minus mu {np.mean(bm[:, m].mean(1) - mu):+.4f}; region sd median {np.median(sd[m]):.3f}; "
              f"region Fisher sd median {np.median(fsd[m]):.2f}")
    order = np.argsort(-W)
    for n in (10, 30, 100, 300, 768):
        s = order[:n]; z = (pm[s] - ref[s]) / sd[s]
        print(f"  top {n:3d} regions by Fisher weight (share {W[s].sum() / W.sum():.2f}; Fisher sd median {np.median(fsd[s]):.3f}, "
              f"posterior sd median {np.median(sd[s]):.3f}): mean(post - ref) {np.mean(pm[s] - ref[s]):+.3f}, z median "
              f"{np.median(z):+.2f}, |z| < 1 {np.mean(np.abs(z) < 1):.2f}, |z| < 2 {np.mean(np.abs(z) < 2):.2f}")


if __name__ == "__main__":
    main()
