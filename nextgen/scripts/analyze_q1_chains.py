"""E014 analysis: compare spectral-index chains from Commander4 maps-only runs by likelihood mode.

Reads chains_compsep/chain0*_iter*.h5 under <chain_root>/<mode>/, extracts the MCMC group's beta
values and acceptance rates per Gibbs iteration, and reports, after burn-in, mean, sd, ESS (Geyer
initial positive sequence) and the shift of each mode from `ridge` in ridge sd. Both chains (the
two pipelined Commander4 chains) are pooled after a per-chain Gelman-Rubin check.

    python -I nextgen/scripts/analyze_q1_chains.py --root ~/data/c4_chains/q1_n64 --burn 400
"""
import argparse
import glob
import json
import os

import h5py
import numpy as np

GROUP = "mcmc/sample_spectral_indices"


def ess_geyer(x):
    x = np.asarray(x, float) - np.mean(x)
    n = x.size
    if n < 10 or np.allclose(x, 0):
        return float(n)
    f = np.fft.rfft(x, 2 * n)
    acov = np.fft.irfft(f * np.conj(f))[:n] / n
    rho = acov / acov[0]
    # initial positive sequence on pair sums
    pairs = rho[1:n - 1:2] + rho[2:n:2]
    tau = -1.0 + 2.0 * rho[0]
    for g in pairs:
        if g <= 0:
            break
        tau += 2.0 * g
    return n / max(tau, 1e-12)


def read_mode(d):
    files = sorted(glob.glob(os.path.join(d, "chains_compsep", "chain*_iter*.h5")))
    out = {}
    for fn in files:
        base = os.path.basename(fn)
        chain = int(base[5:7])
        it = int(base.split("iter")[1].split(".")[0])
        with h5py.File(fn, "r") as f:
            if GROUP not in f:
                continue
            g = f[GROUP]
            row = {k: float(g["params"][k][()]) for k in g["params"]}
            row["accept"] = float(g["accept_rate"][()])
            row["chi2"] = float(f["chi2/reduced"][()])
        out.setdefault(chain, {})[it] = row
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--modes", nargs="*", default=["ridge", "conditional", "marginal", "resolve_chisq"])
    ap.add_argument("--burn", type=int, default=400)
    a = ap.parse_args()
    summary = {}
    for mode in a.modes:
        d = os.path.join(a.root, mode)
        if not os.path.isdir(d):
            print(f"{mode}: no directory")
            continue
        chains = read_mode(d)
        if not chains:
            print(f"{mode}: no chain files")
            continue
        params = [k for k in next(iter(next(iter(chains.values())).values())) if k not in ("accept", "chi2")]
        res = {"n_iter": {c: max(v) for c, v in chains.items()}}
        for p in params:
            per_chain = []
            for c, v in chains.items():
                its = sorted(i for i in v if i > a.burn)
                per_chain.append(np.array([v[i][p] for i in its]))
            pooled = np.concatenate(per_chain)
            means = [x.mean() for x in per_chain]
            vars_ = [x.var(ddof=1) for x in per_chain]
            rhat = np.sqrt((np.var(means, ddof=1) * len(pooled) / len(per_chain) + np.mean(vars_)) / np.mean(vars_)) \
                if len(per_chain) > 1 and np.mean(vars_) > 0 else float("nan")
            ess = sum(ess_geyer(x) for x in per_chain)
            res[p] = {"mean": float(pooled.mean()), "sd": float(pooled.std(ddof=1)), "ess": float(ess),
                      "se": float(pooled.std(ddof=1) / np.sqrt(max(ess, 1))), "rhat": float(rhat),
                      "n": int(pooled.size)}
        acc = np.concatenate([[v[i]["accept"] for i in sorted(v) if i > a.burn] for v in chains.values()])
        chi2 = np.concatenate([[v[i]["chi2"] for i in sorted(v) if i > a.burn] for v in chains.values()])
        res["accept"] = float(np.mean(acc))
        res["red_chi2"] = float(np.mean(chi2))
        summary[mode] = res
    if "ridge" in summary:
        ref = summary["ridge"]
        print(f"{'mode':14s} {'param':6s} {'mean':>10s} {'sd':>8s} {'ESS':>7s} {'Rhat':>6s} "
              f"{'shift/ridge sd':>15s} {'sd ratio':>9s} {'accept':>7s} {'red.chi2':>9s}")
        for mode, res in summary.items():
            for p in [k for k in res if isinstance(res[k], dict) and "mean" in res[k]]:
                r = res[p]
                shift = (r["mean"] - ref[p]["mean"]) / ref[p]["sd"]
                se_shift = np.hypot(r["se"], ref[p]["se"]) / ref[p]["sd"]
                print(f"{mode:14s} {p:6s} {r['mean']:10.4f} {r['sd']:8.4f} {r['ess']:7.0f} {r['rhat']:6.3f} "
                      f"{shift:+8.2f} +- {se_shift:4.2f} {r['sd'] / ref[p]['sd']:9.2f} {res['accept']:7.3f} "
                      f"{res['red_chi2']:9.1f}")
    out = os.path.join(a.root, "summary_q1.json")
    json.dump(summary, open(out, "w"), indent=1)
    print("written", out)


if __name__ == "__main__":
    main()
