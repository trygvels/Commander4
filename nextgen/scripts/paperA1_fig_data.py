"""Paper A1, data figure: the polarized synchrotron index per region, by estimator and by band set (E040).

    python nextgen/scripts/paperA1_fig_data.py --dir nextgen/experiments/E040 --out nextgen/papers/paperA1/figures/fig_data.pdf
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REGIONS = ["High Latitude", "Spur", "Galactic Center", "Plane"]
BANDSETS = [("bpxv_n16_wm", "K Ka Q V"), ("bpxv_n16_KKa_wm", "K Ka"), ("bpxv_n16_KQV_wm", "K Q V"),
            ("bpxv_n16_noK_wm", "Ka Q V"), ("s2_K30_wm", "K 30"), ("s2_30Ka_wm", "30 Ka"),
            ("s2_noK_lfi_wm", "30 44 70 Ka Q V"), ("s2_all_wm", "all eight")]
# published values quoted in the text, with the region they are closest to: (label, value, error)
PUBLISHED = {"Spur": [("Svalheim+ 2023", -3.17, 0.06)], "Plane": [("Svalheim+ 2023", -3.03, 0.07), ("Dunkley+ 2009, low lat.", -3.00, 0.04)],
             "High Latitude": [("Dunkley+ 2009, high lat.", -3.08, 0.06)], "Galactic Center": []}


def load(d, tag):
    fn = os.path.join(d, tag + ".npz")
    return np.load(fn, allow_pickle=True) if os.path.exists(fn) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="nextgen/experiments/E040")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    base = load(a.dir, "bpxv_n16_wm")
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 5.2), sharey=True)
    rows = []
    for r, ax in enumerate(axes):
        y = 0
        ticks, labels = [], []
        # estimators on the baseline band set
        for key, lab, mk, col in (("beta", "profile, sandwich error", "o", "k"), ("flat23", "flat prior at 23 GHz", "<", "C0"),
                                  ("flat30", "flat prior at 30 GHz", ">", "C3"), ("zm", "zero-mean prior", "x", "C2")):
            if base is not None and key in base.files:
                if key == "beta":
                    ax.errorbar(base[key][r], y, xerr=base["sd_sand"][r], fmt=mk, color=col, ms=5, capsize=2)
                else:
                    ax.plot(base[key][r], y, mk, color=col, ms=6)
            ticks.append(y); labels.append(lab); y -= 1
        y -= 0.6
        for tag, lab in BANDSETS:
            z = load(a.dir, tag)
            if z is not None:
                ok = bool(z["ok"][r]) if "ok" in z.files else True
                ax.errorbar(z["beta"][r], y, xerr=z["sd_sand"][r] if ok else 0.0, fmt="s", color="C4" if "30" in lab or "all" in lab else "0.3",
                            ms=4.5, capsize=2)
                rows.append((REGIONS[r], lab, float(z["beta"][r]), float(z["sd_sand"][r])))
            ticks.append(y); labels.append(lab); y -= 1
        y -= 0.6
        for lab, v, e in PUBLISHED[REGIONS[r]]:
            ax.errorbar(v, y, xerr=e, fmt="D", color="C1", ms=4.5, capsize=2)
            ax.text(v, y - 0.45, lab, fontsize=6.5, ha="center", va="top", color="C1")
            y -= 1.2
        if base is not None:
            ax.axvline(base["beta"][r], color="0.8", lw=0.8, zorder=0)
        ax.set_title(REGIONS[r], fontsize=10)
        ax.set_xlabel(r"$\beta_{\rm s}$")
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(y, 0.8)
        ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    fig.savefig(a.out)
    for row in rows:
        print(f"{row[0]:16s} {row[1]:18s} {row[2]:+.4f} +- {row[3]:.4f}")


if __name__ == "__main__":
    main()
