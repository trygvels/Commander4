"""E032: four-chain analysis of index-field runs that wrote both in-flight chains.

For each run tag (a Commander4 output directory chains_<tag>) both chain01 and chain02 files are
read. Reports, against the pre-registered items of the E032 record:

  - per chain: hyperparameter means, sds and effective sample sizes;
  - split R-hat over the 2 x nchain half-chains (Gelman et al. 2013, rank-free version) for the
    hyperparameters and for every region index, pairwise and between-start differences of the
    field means in units of the pooled posterior sd;
  - field level, quartile offsets and calibration against the region-fit reference (mapT);
  - GMRF and amplitude scatters, pixel RMSE against the truth with the constant and region-mean
    baselines, and the largest region excursion.

    python nextgen/scripts/e032_analyze.py --dir D --tags TAG_A TAG_B --refs REFS.npz [--burn 500]

`--links` also creates chains_<tag>_c2/chains_compsep with chain02 files linked under chain01 names,
so that the single-chain check scripts (conditional check, E033 split, multiplier maps) can be run
on the second chain of each run unchanged.
"""
import argparse
import glob
import os
import sys

import h5py
import healpy as hp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e020_analyze import region_map  # noqa: E402

STATE = "mcmc/sample_index_fields/state"
SCALARS = ("index_field/mu", "index_field/log_sigma", "index_field/u", "amp_prior/log_s")
COMPS = (("Synchrotron", "bs", "beta_s", "qs16"), ("ThermalDust", "bd", "beta_d", "qd16"))


def ess(x):
    x = x - x.mean(); n = len(x)
    if not np.any(x):
        return float(n)
    c = np.correlate(x, x, "full")[n - 1:] / np.sum(x * x); s = 0.0
    for k in range(0, n - 1, 2):
        p = c[k] + c[k + 1]
        if p <= 0:
            break
        s += p
    return n / (2 * s - 1) if s > 0.5 else float(n)


def split_rhat(chains):
    """Split R-hat of a list of equal-length 1-d traces (each is cut in two halves)."""
    n = min(len(c) for c in chains) // 2
    seq = np.array([h for c in chains for h in (c[:n], c[n:2 * n])])
    W = seq.var(axis=1, ddof=1).mean()
    B = n * seq.mean(axis=1).var(ddof=1)
    if W <= 0:
        return np.nan
    return float(np.sqrt(((n - 1) / n * W + B / n) / W))


def split_rhat_cols(chains):
    """Split R-hat per column for a list of (ndraw, ncol) arrays."""
    n = min(c.shape[0] for c in chains) // 2
    seq = np.array([h for c in chains for h in (c[:n], c[n:2 * n])])          # (m, n, ncol)
    W = seq.var(axis=1, ddof=1).mean(axis=0)
    B = n * seq.mean(axis=1).var(axis=0, ddof=1)
    return np.sqrt(((n - 1) / n * W + B / n) / np.where(W > 0, W, np.nan))


def read_chain(d, tag, chain, burn):
    files = sorted(glob.glob(f"{d}/chains_{tag}/chains_compsep/chain{chain:02d}_iter*.h5"))
    post = files[burn:]
    if not post:
        return None
    out = {"n": len(files), "npost": len(post), "grad": []}
    for comp, *_ in COMPS:
        for k in SCALARS:
            out[f"{comp}/{k}"] = []
        out[f"{comp}/beta_map"] = []; out[f"{comp}/amp_mu"] = []
    for fn in post:
        with h5py.File(fn, "r") as f:
            st = f[STATE]
            for comp, *_ in COMPS:
                for k in SCALARS:
                    if k in st[comp]:
                        out[f"{comp}/{k}"].append(float(np.asarray(st[comp][k][()]).mean()))
                out[f"{comp}/beta_map"].append(np.asarray(st[comp]["index_field/beta_map"][()], float))
                if "amp_prior/mu" in st[comp]:
                    out[f"{comp}/amp_mu"].append(np.asarray(st[comp]["amp_prior/mu"][()], float))
    for k in list(out):
        if isinstance(out[k], list):
            out[k] = np.array(out[k])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--tags", nargs="+", required=True)
    ap.add_argument("--refs", required=True); ap.add_argument("--burn", type=int, default=500)
    ap.add_argument("--ref-key", default="mapT"); ap.add_argument("--links", action="store_true")
    a = ap.parse_args()
    truth = np.load(f"{a.dir}/truth.npz", allow_pickle=True)
    nside = int(truth["nside"]); reg = region_map(16, nside); refs = np.load(a.refs)

    if a.links:
        for tag in a.tags:
            src = sorted(glob.glob(f"{a.dir}/chains_{tag}/chains_compsep/chain02_iter*.h5"))
            dst = f"{a.dir}/chains_{tag}_c2/chains_compsep"; os.makedirs(dst, exist_ok=True)
            for fn in src:
                ln = os.path.join(dst, os.path.basename(fn).replace("chain02_", "chain01_"))
                if not os.path.lexists(ln):
                    os.symlink(fn, ln)
            print(f"# linked {len(src)} chain02 files of {tag} as chains_{tag}_c2")

    chains, labels = [], []
    for it, tag in enumerate(a.tags):
        for c in (1, 2):
            r = read_chain(a.dir, tag, c, a.burn)
            if r is None:
                print(f"# {tag} chain {c}: no post-burn files"); continue
            chains.append(r); labels.append(f"{'ABCDEFGH'[it]}{c}")
            print(f"# chain {labels[-1]}: {tag} chain{c:02d}, {r['n']} files, {r['npost']} post-burn")
    if not chains:
        return
    start_of = [lab[0] for lab in labels]

    print("\n## 1. Hyperparameters per chain: mean, sd, ESS; then split R-hat over all half-chains, the largest pairwise")
    print("##    difference of chain means and the between-start difference, both in units of the pooled posterior sd")
    for comp, *_ in COMPS:
        for k in SCALARS:
            key = f"{comp}/{k}"
            if any(len(c[key]) == 0 for c in chains):
                continue
            tr = [c[key] for c in chains]
            pooled_sd = np.sqrt(np.mean([t.var(ddof=1) for t in tr]))
            means = np.array([t.mean() for t in tr])
            pair = (means.max() - means.min()) / pooled_sd if pooled_sd > 0 else np.nan
            starts = sorted(set(start_of))
            bs = [np.mean([m for m, s in zip(means, start_of) if s == st]) for st in starts]
            between = (bs[0] - bs[-1]) / pooled_sd if len(starts) > 1 and pooled_sd > 0 else np.nan
            per = "  ".join(f"{lab} {t.mean():+.4f} {t.std(ddof=1):.4f} {ess(t):5.0f}" for lab, t in zip(labels, tr))
            print(f"{key:38s} {per}")
            print(f"{'':38s} split R-hat {split_rhat(tr):.3f}  max pairwise {pair:.2f} sd  between starts {between:+.2f} sd"
                  f"  pooled mean {np.mean(means):+.4f} sd {pooled_sd:.4f}")

    print("\n## 2. Region indices: split R-hat over all half-chains, per component")
    for comp, *_ in COMPS:
        rh = split_rhat_cols([c[f"{comp}/beta_map"] for c in chains])
        e = np.array([[ess(c[f"{comp}/beta_map"][:, i]) for i in range(0, rh.size, 16)] for c in chains])
        print(f"{comp:12s} R-hat median {np.nanmedian(rh):.3f}  95 pct {np.nanpercentile(rh, 95):.3f}  max {np.nanmax(rh):.3f}  "
              f"fraction below 1.05: {np.mean(rh < 1.05):.3f}  below 1.10: {np.mean(rh < 1.10):.3f}  "
              f"(ESS per chain on every 16th region: median {np.median(e):.0f}, min {e.min():.0f})")

    print("\n## 3. Field level, quartiles and calibration against the region-fit reference; pooled over chains, then per chain")
    for comp, rk, tk, qk in COMPS:
        ref = refs[f"{rk}_{a.ref_key}"]; q = refs[qk]; tmap = np.asarray(truth[tk], float)
        tmean = tmap.mean()
        pooled = np.concatenate([c[f"{comp}/beta_map"] for c in chains], axis=0)
        sets = [("pooled", pooled)] + [(lab, c[f"{comp}/beta_map"]) for lab, c in zip(labels, chains)]
        const_rmse = np.sqrt(np.mean((tmap - tmean) ** 2))
        floor = np.sqrt(np.mean((tmap - (np.bincount(reg, weights=tmap) / np.bincount(reg))[reg]) ** 2))
        print(f"# {comp}: truth mean {tmean:+.4f}, reference ({a.ref_key}) mean {ref.mean():+.4f}; baselines: constant at the "
              f"truth mean RMSE {const_rmse:.4f}, truth region means {floor:.4f}")
        for lab, b in sets:
            pm, psd = b.mean(axis=0), b.std(axis=0, ddof=1)
            z = (pm - ref) / psd
            qo = " ".join(f"{(pm - ref)[q == k].mean():+.4f}" for k in range(4))
            cal = " ".join(f"{np.mean(np.abs(z[q == k]) < 1):.2f}" for k in range(4))
            zm = " ".join(f"{np.median(z[q == k]):+.2f}" for k in range(4))
            rmse = np.sqrt(np.mean((pm[reg] - tmap) ** 2))
            print(f"  {lab:7s} field mean {pm.mean():+.4f}  minus ref {pm.mean() - ref.mean():+.4f}  minus truth {pm.mean() - tmean:+.4f}"
                  f" | quartile offsets Q1..Q4 {qo} | |z|<1 {cal} | median z {zm} | pixel RMSE {rmse:.4f}"
                  f" | max |pm - ref| {np.abs(pm - ref).max():.3f}")

    print("\n## 4. Scatters per chain: GMRF sigma = exp(log_sigma), amplitude scatter log s")
    for comp, *_ in COMPS:
        line = []
        for lab, c in zip(labels, chains):
            ls = c[f"{comp}/index_field/log_sigma"]; s = c[f"{comp}/amp_prior/log_s"]
            line.append(f"{lab} sigma {np.exp(ls).mean():.4f}+-{np.exp(ls).std():.4f} log_s "
                        f"{(s.mean() if len(s) else np.nan):+.3f}")
        print(f"{comp:12s} " + "  ".join(line))

    print("\n## 5. Template multiplier maps: mean over regions and chains agreement (correlation of posterior means between chains)")
    for comp, *_ in COMPS:
        ms = [c[f"{comp}/amp_mu"].mean(axis=0) for c in chains if c[f"{comp}/amp_mu"].ndim == 2]
        if len(ms) > 1:
            cc = [np.corrcoef(ms[0], m)[0, 1] for m in ms[1:]]
            print(f"{comp:12s} map means " + " ".join(f"{m.mean():.4f}" for m in ms) + "  corr with first chain " +
                  " ".join(f"{x:.3f}" for x in cc))


if __name__ == "__main__":
    main()
