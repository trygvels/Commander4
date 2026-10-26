"""Bias and coverage curves from the coverage.py JSONs (E012 v2 Reference 3, no pass/fail).

For every recorded quantity q with truth q*, posterior mean m, sd s and central-interval
quantiles (q_2.5, q_16, q_84, q_97.5):
    standardized bias   z = (m - q*) / s,
    coverage_68         1[q_16  <= q* <= q_84],
    coverage_95         1[q_2.5 <= q* <= q_97.5].
A calibrated posterior has mean z = 0, sd z = 1, coverage 0.68 and 0.95 (up to binomial error
sqrt(c (1 - c) / n), reported). Index pixels are pooled over replicates and binned by the per-pixel
synchrotron S/N proxy snr_index of coverage.py (log10 scale); bin edges are quantiles of the pooled
S/N of the whole directory (--n-bins) or given by --edges, so that the modes and noise scales share
the same bins. Hyperparameters (mu, log sigma, u) are summarized per (noise scale, mode) over
replicates. Also reported: the median posterior sd (width), the median absolute error |m - q*|,
ESS, divergences and zero-density evaluations per cell.

Writes coverage_summary.txt and coverage_summary.json into --dir and prints the text table.

Run from nextgen/prototypes as
    python -m indexfield.coverage_analyze --dir indexfield/checks/out/coverage/<tag>
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

from indexfield import coverage

HYPER = ("mu", "log_sigma", "u")


def load(dir_: str) -> list[dict]:
    recs = []
    for fn in sorted(glob.glob(os.path.join(dir_, "rep*_noise*_*.json"))):
        with open(fn) as f:
            r = json.load(f)
        for k in ("truth", "mean", "sd", "ess", "snr_index", "snr_pix_mean"):
            r[k] = np.asarray(r[k], dtype=float)
        r["quantiles"] = np.asarray(r["quantiles"], dtype=float)      # (4, nq)
        recs.append(r)
    if not recs:
        raise SystemExit(f"no rep*_noise*_*.json in {dir_}")
    return recs


def calibration(truth, mean, sd, q):
    """z, inside-68, inside-95 arrays for stacked quantities (q: (4, n) quantiles)."""
    z = (mean - truth) / sd
    in68 = (q[1] <= truth) & (truth <= q[2])
    in95 = (q[0] <= truth) & (truth <= q[3])
    return z, in68, in95


def summarize(z, in68, in95, sd, abs_err, ess=None, rep=None):
    """Pooled statistics of z. `rep` (replicate index per entry) adds the replicate-clustered
    standard error of mean z: pixels of one replicate share the field hyperparameters, the truth
    and the noise realization, so their z are correlated and the pooled SE understates (E016
    residual check: pooled SE 0.04, clustered SE 0.15 for 12 replicates x 48 pixels)."""
    n = int(z.size)
    if n == 0:
        return dict(n=0)
    c68, c95 = float(in68.mean()), float(in95.mean())
    out = dict(n=n, bias_mean=float(z.mean()), bias_se=float(z.std(ddof=1) / np.sqrt(n)) if n > 1
               else float("nan"), bias_sd=float(z.std(ddof=1)) if n > 1 else float("nan"),
               bias_rms=float(np.sqrt(np.mean(z ** 2))),
               cov68=c68, cov68_se=float(np.sqrt(c68 * (1 - c68) / n)),
               cov95=c95, cov95_se=float(np.sqrt(c95 * (1 - c95) / n)),
               sd_median=float(np.median(sd)), abs_err_median=float(np.median(abs_err)))
    if rep is not None:
        reps = np.unique(rep)
        means = np.array([z[rep == r].mean() for r in reps])
        out["n_rep"] = int(reps.size)
        out["bias_se_rep"] = float(means.std(ddof=1) / np.sqrt(reps.size)) if reps.size > 1 \
            else float("nan")
        out["bias_mean_rep"] = means.tolist()
    if ess is not None:
        out["ess_median"] = float(np.median(ess))
        out["ess_min"] = float(np.min(ess))
    return out


def fmt_row(label, s):
    if s.get("n", 0) == 0:
        return f"  {label:<22s} {'(no data)':>9s}"
    rep_se = f" [{s['bias_se_rep']:4.2f}]" if "bias_se_rep" in s else " " * 7
    return (f"  {label:<22s} {s['n']:>6d}  {s['bias_mean']:+6.2f} +- {s['bias_se']:4.2f}"
            f"{rep_se}  {s['bias_sd']:5.2f}  {s['cov68']:5.2f} +- {s['cov68_se']:4.2f}  "
            f"{s['cov95']:5.2f} +- {s['cov95_se']:4.2f}  {s['sd_median']:8.4f}  "
            f"{s['abs_err_median']:8.4f}" + (f"  {s['ess_median']:7.0f}" if "ess_median" in s
                                            else ""))


def bin_label(lo, hi):
    """'beta S/N [a, b)' with open ends for infinite edges (edges in log10)."""
    if np.isfinite(lo) and np.isfinite(hi):
        return f"beta S/N [{10 ** lo:.3g}, {10 ** hi:.3g})"
    if np.isfinite(hi):
        return f"beta S/N < {10 ** hi:.3g}"
    return f"beta S/N >= {10 ** lo:.3g}"


HEADER = (f"  {'bin / quantity':<22s} {'n':>6s}  {'mean z [rep SE]':>21s}  {'sd z':>5s}  "
          f"{'cov68':>13s}  {'cov95':>13s}  {'med sd':>8s}  {'med |err|':>8s}  {'med ESS':>7s}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True, help="coverage.py output directory")
    ap.add_argument("--n-bins", type=int, default=5, help="quantile S/N bins (log10 snr_index)")
    ap.add_argument("--edges", type=float, nargs="+", default=None,
                    help="explicit log10 S/N bin edges instead of quantiles")
    ap.add_argument("--snr", default="snr_index", choices=["snr_index", "snr_pix_mean"])
    args = ap.parse_args(argv)
    recs = load(args.dir)
    cfg_path = os.path.join(args.dir, "config.json")
    config = json.load(open(cfg_path)) if os.path.exists(cfg_path) else {}
    noises = sorted({r["noise_scale"] for r in recs})
    modes = [m for m in coverage.MODES if any(r["mode"] == m for r in recs)]
    n_beta = recs[0]["n_beta"]
    lsnr_all = np.log10(np.concatenate([r[args.snr] for r in recs]))
    if args.edges is not None:
        edges = np.asarray(args.edges, dtype=float)
    else:
        edges = np.quantile(lsnr_all, np.linspace(0, 1, args.n_bins + 1))
        edges[0], edges[-1] = -np.inf, np.inf
    lines = [f"coverage study: {args.dir}", f"{len(recs)} cells; noise scales {noises}; modes "
             f"{modes}; n_beta {n_beta}; S/N variable {args.snr}; log10 S/N bin edges "
             f"{np.round(edges, 3).tolist()}", "config: " + json.dumps(
                 {k: config[k] for k in config if k not in ('out_dir',)}, default=str)]
    summary = dict(dir=args.dir, config=config, snr_variable=args.snr, edges=edges.tolist(),
                   noise_scales=noises, modes=modes, n_beta=n_beta, cells={})
    for ns in noises:
        for mode in modes:
            rr = [r for r in recs if r["noise_scale"] == ns and r["mode"] == mode]
            if not rr:
                continue
            key = f"noise{ns:g}_{mode}"
            cell = dict(noise_scale=ns, mode=mode, n_rep=len(rr),
                        reps=sorted(r["rep"] for r in rr),
                        n_divergent=int(sum(r["n_divergent"] for r in rr)),
                        n_cells_with_divergences=int(sum(r["n_divergent"] > 0 for r in rr)),
                        n_loglike_errors_sampling=int(sum(r["n_loglike_errors_sampling"]
                                                          for r in rr)),
                        ebfmi_min=float(min(r["ebfmi"] for r in rr)),
                        seconds_mean=float(np.mean([r["seconds_total"] for r in rr])),
                        ess_min_median=float(np.median([r["ess"].min() for r in rr])),
                        step_size_median=float(np.median([r["step_size"] for r in rr])),
                        hyper={}, beta_bins=[])
            lines += ["", f"noise x{ns:g}, mode {mode}: {len(rr)} replicates, divergences "
                      f"{cell['n_divergent']} (in {cell['n_cells_with_divergences']} chains), "
                      f"zero-density evaluations {cell['n_loglike_errors_sampling']}, min E-BFMI "
                      f"{cell['ebfmi_min']:.2f}, median min ESS {cell['ess_min_median']:.0f}, "
                      f"{cell['seconds_mean']:.0f} s per chain", HEADER]
            for j, name in enumerate(HYPER):
                tr = np.array([r["truth"][j] for r in rr])
                mn = np.array([r["mean"][j] for r in rr])
                sd = np.array([r["sd"][j] for r in rr])
                q = np.column_stack([r["quantiles"][:, j] for r in rr])
                ess = np.array([r["ess"][j] for r in rr])
                z, in68, in95 = calibration(tr, mn, sd, q)
                s = summarize(z, in68, in95, sd, np.abs(mn - tr), ess)
                s["truth"] = tr.tolist()
                s["mean"] = mn.tolist()
                s["sd"] = sd.tolist()
                cell["hyper"][name] = s
                lines.append(fmt_row(name, s))
            # beta pixels pooled over replicates, binned by S/N
            tr = np.concatenate([r["truth"][3:] for r in rr])
            mn = np.concatenate([r["mean"][3:] for r in rr])
            sd = np.concatenate([r["sd"][3:] for r in rr])
            q = np.concatenate([r["quantiles"][:, 3:] for r in rr], axis=1)
            ess = np.concatenate([r["ess"][3:] for r in rr])
            lsnr = np.log10(np.concatenate([r[args.snr] for r in rr]))
            rep = np.concatenate([np.full(r["truth"].size - 3, r["rep"]) for r in rr])
            z, in68, in95 = calibration(tr, mn, sd, q)
            s_all = summarize(z, in68, in95, sd, np.abs(mn - tr), ess, rep)
            cell["beta_all"] = s_all
            lines.append(fmt_row("beta (all pixels)", s_all))
            for b in range(len(edges) - 1):
                sel = (lsnr >= edges[b]) & (lsnr < edges[b + 1])
                s = summarize(z[sel], in68[sel], in95[sel], sd[sel], np.abs(mn - tr)[sel], ess[sel],
                              rep[sel])
                s.update(lo=float(edges[b]), hi=float(edges[b + 1]),
                         snr_median=float(10 ** np.median(lsnr[sel])) if sel.any() else None)
                cell["beta_bins"].append(s)
                lines.append(fmt_row(bin_label(edges[b], edges[b + 1]), s))
            summary["cells"][key] = cell
    lines += ["", "columns: n = pooled count; mean z = mean (post. mean - truth) / post. sd with "
              "its pooled SE and, for beta rows, the replicate-clustered SE in brackets (the one "
              "to use: pixels of a replicate share truth and noise); sd z = sd of z (1 if "
              "calibrated); cov68, cov95 = fraction of truths inside "
              "the central 68 / 95 percent intervals with binomial SE; med sd = median posterior "
              "sd; med |err| = median |post. mean - truth|; med ESS = median Geyer ESS. Curves, "
              "not a test (E012 principle 4)."]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(args.dir, "coverage_summary.txt"), "w") as f:
        f.write(text + "\n")
    with open(os.path.join(args.dir, "coverage_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(f"\nwrote {os.path.join(args.dir, 'coverage_summary.txt')} and .json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
