"""E040: collect the runs under nextgen/experiments/E040 into tables.

    python nextgen/scripts/e040_summarize.py [--dir nextgen/experiments/E040] [--ref bpxv_n16_wm]

Prints one row per run with the profile index and its debiased sandwich error per region, the shift from
the reference run in units of the combined error, and the chi-square per degree of freedom.
"""
import argparse
import glob
import os
import re

import numpy as np


def header_info(txt):
    out = {}
    with open(txt) as f:
        for line in f:
            if line.startswith("# profile fit"):
                m = re.search(r"chi2_min ([\d.]+); nominal dof (\d+).*ratio ([\d.]+)", line)
                if m:
                    out["chi2"], out["dof"], out["ratio"] = float(m.group(1)), int(m.group(2)), float(m.group(3))
                out["converged"] = "converged True" in line
            elif line.startswith("# chi2 per band"):
                out["perband"] = line.split(":", 1)[1].strip()
            elif line.startswith("# chi2 by region"):
                out["perregion"] = line.split(":", 1)[1].strip()
            elif line.startswith("# T / J_c"):
                out["tj"] = line.split(":", 1)[1].strip()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="nextgen/experiments/E040")
    ap.add_argument("--ref", default="bpxv_n16_wm")
    ap.add_argument("--match", default="bpxv")
    a = ap.parse_args()
    runs = {}
    for fn in sorted(glob.glob(os.path.join(a.dir, "*.npz"))):
        tag = os.path.basename(fn)[:-4]
        txt = fn[:-4] + ".txt"
        if not os.path.exists(txt):
            continue
        z = np.load(fn, allow_pickle=True)
        runs[tag] = (z, header_info(txt))
    if a.ref not in runs:
        print("reference run", a.ref, "not found; runs:", list(runs))
        return
    zr = runs[a.ref][0]
    br, sr = zr["beta"], zr["sd_sand"]
    print(f"reference {a.ref}: " + "  ".join(f"{b:+.4f} +- {s:.4f}" for b, s in zip(br, sr)))
    print("| run | " + " | ".join(f"region {i}" for i in range(br.size)) + " | shift from reference in combined sd | chi2 / dof | converged |")
    print("|---|" + "---|" * (br.size + 3))
    for tag, (z, h) in runs.items():
        if a.match not in tag or z["beta"].size != br.size:
            continue
        b, s = z["beta"], z["sd_sand"]
        comb = np.sqrt(s ** 2 + sr ** 2)
        print(f"| {tag} | " + " | ".join(f"{x:+.4f} +- {e:.4f}" for x, e in zip(b, s)) + " | "
              + " ".join(f"{(x - y) / c:+.1f}" for x, y, c in zip(b, br, comb)) + f" | {h.get('ratio', float('nan')):.4f} | {h.get('converged')} |")
    print()
    for tag, (z, h) in runs.items():
        if "flat23" in z.files:
            print(f"## {tag}: estimators")
            nper = z["nper"]
            J, T = z["J"], z["T"]
            Jc = J - T
            pred = np.linalg.solve(Jc, nper * np.log(30.0 / 23.0))
            print("| region | n | profile | sd curvature | sd sandwich | flat23 - profile | flat30 - profile | flat30 - flat23 | predicted | ratio | zero-mean - profile | in sandwich sd |")
            print("|---|---|---|---|---|---|---|---|---|---|---|---|")
            for r in range(z["beta"].size):
                d = z["flat30"][r] - z["flat23"][r]
                zm = z["zm"][r] - z["beta"][r] if "zm" in z.files else float("nan")
                print(f"| {r} | {nper[r]} | {z['beta'][r]:+.4f} | {z['sd_naive'][r]:.4f} | {z['sd_sand'][r]:.4f} | {z['flat23'][r] - z['beta'][r]:+.4f} | "
                      f"{z['flat30'][r] - z['beta'][r]:+.4f} | {d:+.4f} | {pred[r]:+.4f} | {d / pred[r]:.3f} | {zm:+.4f} | {zm / z['sd_sand'][r]:+.2f} |")
            print("chi2 per band:", h.get("perband"))
            print("chi2 by region:", h.get("perregion"))
            print("T / J_c:", h.get("tj"))
            print()


if __name__ == "__main__":
    main()
