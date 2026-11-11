"""Paper A1, data figure: the polarized synchrotron index per region, by estimator and by band set
(E040, E042).

    python nextgen/scripts/paperA1_fig_data.py --out nextgen/papers/paperA1/figures/fig_data.pdf
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REGIONS = ["High Latitude", "Spur", "Galactic Center", "Plane"]
E40, E42 = "nextgen/experiments/E040", "nextgen/experiments/E042"
PUBLISHED = {"Spur": [("Svalheim+ 2023", -3.17, 0.06)], "Plane": [("Svalheim+ 2023", -3.03, 0.07), ("Dunkley+ 2009, $|b|<20^\\circ$", -3.00, 0.04)],
             "High Latitude": [("Dunkley+ 2009, $|b|>20^\\circ$", -3.08, 0.06)], "Galactic Center": []}


def npz(path):
    return np.load(path, allow_pickle=True) if os.path.exists(path) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    base = npz(f"{E40}/bpxv_n16_wm.npz")
    rows = [("head", "same maps, four ways to remove the amplitudes")]
    if base is not None:
        rows += [("pt", "profile likelihood", base["beta"], base["sd_sand"], "k", "o"),
                 ("pt", "flat prior, reference 23 GHz", base["flat23"], None, "C0", "<"),
                 ("pt", "flat prior, reference 30 GHz", base["flat30"], None, "C3", ">"),
                 ("pt", "zero-mean prior, free scale", base["zm"], None, "C2", "x")]
    rows.append(("head", "band sets, all scales"))
    for tag, lab, col in (("bpxv_n16_wm", "K Ka Q V", "0.25"), ("bpxv_n16_KKa_wm", "K Ka", "0.25"), ("s2_K30_wm", "K 30 (PR4)", "C4"),
                          ("s2_noK_lfi_wm", "30 44 70 Ka Q V", "C4")):
        z = npz(f"{E40}/{tag}.npz")
        if z is not None:
            rows.append(("pt", lab, z["beta"], z["sd_sand"], col, "s"))
    z = npz(f"{E42}/pr3_arm.npz")
    if z is not None and "K_30PR3_353" in z.files:
        v = z["K_30PR3_353"]
        rows.append(("pt", "K 30 (PR3)", v[:4], v[4:], "C5", "s"))
    rows.append(("head", "band sets, index from multipoles 2 to 20"))
    z = npz(f"{E42}/lowpass20.npz")
    if z is not None:
        for key, lab, col in (("K_Ka_Q_V_353", "K Ka Q V", "0.25"), ("K_Ka_353", "K Ka", "0.25"), ("K_30_353", "K 30 (PR4)", "C4"),
                              ("30_Ka_353", "30 (PR4) Ka", "C4")):
            if key in z.files:
                rows.append(("pt", lab, z[key][:4], z[key][4:], col, "D"))
    rows.append(("head", "published"))
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 6.0), sharey=True)
    for r, ax in enumerate(axes):
        y, ticks, labels = 0.0, [], []
        for row in rows:
            if row[0] == "head":
                y -= 0.5
                if r == 0:
                    ax.text(-0.02, y, row[1], transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8, style="italic")
                y -= 1.0
                continue
            _, lab, val, err, col, mk = row
            if err is not None and np.isfinite(err[r]):
                ax.errorbar(val[r], y, xerr=err[r], fmt=mk, color=col, ms=4.5, capsize=2, lw=1)
            else:
                ax.plot(val[r], y, mk, color=col, ms=6)
            ticks.append(y); labels.append(lab); y -= 1.0
        for lab, v, e in PUBLISHED[REGIONS[r]]:
            ax.errorbar(v, y, xerr=e, fmt="^", color="C1", ms=4.5, capsize=2, lw=1)
            ax.text(v, y - 0.42, lab, fontsize=6.5, ha="center", va="top", color="C1")
            y -= 1.3
        if base is not None:
            ax.axvline(base["beta"][r], color="0.85", lw=0.8, zorder=0)
        ax.set_title(REGIONS[r], fontsize=10)
        ax.set_xlabel(r"$\beta_{\rm s}$")
        ax.set_xlim(-3.75, -2.2)
        ax.set_xticks([-3.5, -3.0, -2.5])
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_ylim(y, 0.2)
        ax.tick_params(axis="x", labelsize=8)
    fig.subplots_adjust(left=0.27, right=0.99, top=0.95, bottom=0.08, wspace=0.06)
    fig.savefig(a.out)
    print("rows:", [r[1] for r in rows])


if __name__ == "__main__":
    main()
