"""Bias and coverage table for the scalar-index study (E017) from scalar_coverage.py JSONs.

Per (noise scale, mode): mean z = (post. mean - truth) / post. sd over replicates with its SE
(replicates are independent, so this is the right SE), sd z, 68 and 95 percent coverage with
binomial SE, median posterior sd, median |error|, the full-sky S/N (median over replicates) and,
for the hier modes, ESS and sampler diagnostics plus the amplitude hyperparameter means. Also the
mean shift between modes on the same data in units of the w0 sd (the E014 "chain minus grid"
style comparison), since all modes see the same replicate.

Run from nextgen/prototypes as
    python -m indexfield.scalar_coverage_analyze --dir indexfield/checks/out/scalar/<tag>
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

from indexfield import scalar_coverage as sc
from indexfield.coverage_analyze import summarize, calibration


def load(dir_):
    recs = []
    for fn in sorted(glob.glob(os.path.join(dir_, "rep*_noise*_*.json"))):
        with open(fn) as f:
            recs.append(json.load(f))
    if not recs:
        raise SystemExit(f"no rep*_noise*_*.json in {dir_}")
    return recs


HEADER = (f"  {'mode':<13s} {'n':>3s}  {'S/N':>5s}  {'mean z':>13s}  {'sd z':>5s}  {'cov68':>13s}  "
          f"{'cov95':>13s}  {'med sd':>8s}  {'med |err|':>9s}  {'shift/w0 sd':>12s}  {'ESS':>5s}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    args = ap.parse_args(argv)
    recs = load(args.dir)
    noises = sorted({r["noise_scale"] for r in recs})
    modes = [m for m in sc.MODES if any(r["mode"] == m for r in recs)]
    lines = [f"scalar-index coverage study: {args.dir}", f"{len(recs)} cells; noise scales "
             f"{noises}; modes {modes}"]
    summary = dict(dir=args.dir, noise_scales=noises, modes=modes, cells={})
    for ns in noises:
        lines += ["", f"noise x{ns:g}", HEADER]
        w0 = {r["rep"]: r for r in recs if r["noise_scale"] == ns and r["mode"] == "w0"}
        for mode in modes:
            rr = sorted([r for r in recs if r["noise_scale"] == ns and r["mode"] == mode],
                        key=lambda r: r["rep"])
            if not rr:
                continue
            tr = np.array([r["truth"][0] for r in rr])
            mn = np.array([r["mean"][0] for r in rr])
            sd = np.array([r["sd"][0] for r in rr])
            q = np.column_stack([np.asarray(r["quantiles"])[:, 0] for r in rr])
            z, in68, in95 = calibration(tr, mn, sd, q)
            ess = np.array([r["ess"][0] for r in rr]) if "ess" in rr[0] else None
            s = summarize(z, in68, in95, sd, np.abs(mn - tr), ess)
            s["snr_median"] = float(np.median([r["snr_full_sky"] for r in rr]))
            s["reps"] = [r["rep"] for r in rr]
            common = [r for r in rr if r["rep"] in w0]
            if common and mode != "w0":
                sh = np.array([(r["mean"][0] - w0[r["rep"]]["mean"][0]) / w0[r["rep"]]["sd"][0]
                               for r in common])
                s["shift_w0_mean"] = float(sh.mean())
                s["shift_w0_se"] = float(sh.std(ddof=1) / np.sqrt(sh.size)) if sh.size > 1 \
                    else float("nan")
                shift = f"{s['shift_w0_mean']:+6.2f} +- {s['shift_w0_se']:4.2f}"
            else:
                shift = " " * 14
            if mode.startswith("hier"):
                amp = np.array([r["mean"][1:] for r in rr])
                amp_sd = np.array([r["sd"][1:] for r in rr])
                s["amp_mean"] = amp.mean(axis=0).tolist()
                s["amp_sd_median"] = np.median(amp_sd, axis=0).tolist()
                s["n_divergent"] = int(sum(r["n_divergent"] for r in rr))
                s["bad_reps"] = [r["rep"] for r in rr if r["n_divergent"] > 0 or r["ebfmi"] < 0.3]
                s["ebfmi_min"] = float(min(r["ebfmi"] for r in rr))
                s["seconds_mean"] = float(np.mean([r["seconds_total"] for r in rr]))
            summary["cells"][f"noise{ns:g}_{mode}"] = s
            lines.append(f"  {mode:<13s} {s['n']:>3d}  {s['snr_median']:5.0f}  "
                         f"{s['bias_mean']:+6.2f} +- {s['bias_se']:4.2f}  {s['bias_sd']:5.2f}  "
                         f"{s['cov68']:5.2f} +- "
                         f"{s['cov68_se']:4.2f}  {s['cov95']:5.2f} +- {s['cov95_se']:4.2f}  "
                         f"{s['sd_median']:8.2e}  {s['abs_err_median']:9.2e}  {shift}"
                         + (f"  {s['ess_median']:5.0f}" if "ess_median" in s else ""))
        for mode in modes:
            s = summary["cells"].get(f"noise{ns:g}_{mode}")
            if s and "amp_mean" in s:
                a, w = s["amp_mean"], s["amp_sd_median"]
                lines.append(f"    {mode}: (mu_sync, mu_dust, log s_sync, log s_dust) mean "
                             f"over replicates ({a[0]:.3f}, {a[1]:.3f}, {a[2]:.2f}, {a[3]:.2f}), "
                             f"median post. sd ({w[0]:.3f}, {w[1]:.3f}, {w[2]:.2f}, "
                             f"{w[3]:.2f}); generator "
                             f"(1.027, 1.027, {np.log(sc.GEN_S):.2f}, {np.log(sc.GEN_S):.2f}); "
                             f"divergences {s['n_divergent']}, min E-BFMI {s['ebfmi_min']:.2f}, "
                             f"{s['seconds_mean']:.0f} s per chain"
                             + (f"; CHAINS WITH DIVERGENCES OR E-BFMI < 0.3 (rows above include "
                                f"them): reps {s['bad_reps']}" if s["bad_reps"] else ""))
    lines += ["", "columns: n = replicates; S/N = median full-sky S/N proxy; mean z = mean "
              "(post. mean - truth) / post. sd with its SE over replicates; sd z = sd of z (1 if "
              "calibrated); cov68, cov95 = fraction of truths inside the central intervals with "
              "binomial SE; med sd = median posterior sd; med |err| = median |post. mean - "
              "truth|; shift/w0 sd = mean over replicates of (mode mean - w0 mean) / w0 sd on the "
              "same data; ESS = median Geyer ESS of beta (hier modes; grid modes are exact)."]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(args.dir, "scalar_summary.txt"), "w") as f:
        f.write(text + "\n")
    with open(os.path.join(args.dir, "scalar_summary.json"), "w") as f:
        json.dump(summary, f, indent=1, default=float)
    print(f"\nwrote {os.path.join(args.dir, 'scalar_summary.txt')} and .json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
