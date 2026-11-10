"""E039: summary of the 48-region simulation from one or more output files of e039_lowres_validation.py
(the run was split after a restart). Uses the truth table of the header and the per-realization optima.

    python nextgen/scripts/e039_merge_hp2.py nextgen/experiments/E039/sim_n8_hp2_part1.txt [part2.txt ...]
"""
import sys

import numpy as np


def parse(fn):
    truth, P, chi = {}, {}, {}
    for line in open(fn):
        w = line.split()
        if line.startswith("#") and len(w) == 7 and w[1].isdigit():
            truth[int(w[1])] = (int(w[2]), float(w[3]), float(w[4]), float(w[5]))
        elif line.startswith("real "):
            i = w.index("chi2/dof")
            P[int(w[1])] = np.array([float(x) for x in w[3:i]])
            chi[int(w[1])] = float(w[i + 1])
    return truth, P, chi


def main():
    truth, P, chi = {}, {}, {}
    for fn in sys.argv[1:]:
        t, p, c = parse(fn)
        truth.update(t); P.update(p); chi.update(c)
    nreg = len(truth)
    bt = np.array([truth[r][1] for r in range(nreg)])
    sdn = np.array([truth[r][2] for r in range(nreg)])
    sds = np.array([truth[r][3] for r in range(nreg)])
    keys = sorted(P)
    X = np.array([P[k] for k in keys])
    zs, zn = (X - bt) / sds, (X - bt) / sdn
    t = (sds / sdn) ** 2 - 1.0
    print(f"# {len(keys)} realizations ({keys[0]} to {keys[-1]}), {nreg} regions; sandwich sd at the truth: min {sds.min():.3f}, median {np.median(sds):.3f}, "
          f"max {sds.max():.3f}; T/J at the truth: min {t.min():.3f}, median {np.median(t):.3f}, max {t.max():.3f}")
    edge = np.mean((X < -4.59) | (X > -1.61))
    print(f"# fraction of region fits on a bound: {edge:.4f}; chi2/dof mean {np.mean([chi[k] for k in keys]):.4f}, sd {np.std([chi[k] for k in keys], ddof=1):.4f}")
    print("class by sandwich sd | regions | mean z +- se | sd z sandwich | sd z naive | predicted naive sqrt(1+T/J) | coverage 68 sandwich | naive")
    for lo, hi in ((0, 0.05), (0.05, 0.1), (0.1, 0.25), (0.25, 9)):
        sel = (sds >= lo) & (sds < hi)
        if not sel.any():
            continue
        a, b = zs[:, sel], zn[:, sel]
        # standard error of the mean with realizations as clusters
        se = np.std(a.mean(axis=1), ddof=1) / np.sqrt(a.shape[0])
        print(f"{lo:g} to {hi:g} | {int(sel.sum())} | {a.mean():+.3f} +- {se:.3f} | {a.std(ddof=1):.3f} | {b.std(ddof=1):.3f} | "
              f"{np.sqrt(np.mean(1 + t[sel])):.3f} | {np.mean(np.abs(a) < 1):.3f} | {np.mean(np.abs(b) < 1):.3f}")
    ok = sds < 0.25
    a = zs[:, ok]
    se = np.std(a.mean(axis=1), ddof=1) / np.sqrt(a.shape[0])
    print(f"all {int(ok.sum())} regions with sandwich sd < 0.25: mean z {a.mean():+.3f} +- {se:.3f}, sd z sandwich {a.std(ddof=1):.3f}, "
          f"naive {zn[:, ok].std(ddof=1):.3f}, coverage {np.mean(np.abs(a) < 1):.3f}; largest abs mean z of a region {np.abs(a.mean(axis=0)).max():.3f} "
          f"(se of one region {1 / np.sqrt(a.shape[0]):.3f})")


if __name__ == "__main__":
    main()
