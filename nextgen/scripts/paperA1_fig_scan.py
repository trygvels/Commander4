"""Paper A1: the flat-prior marginal along the index of the High Latitude region, for two priors on the CMB.

Reads the scans of `e044_flat_scan.py` (the other three indices held at the profile optimum) and draws, for the LCDM
prior and for a flat prior on the CMB, the profile log-likelihood and the flat-prior marginal log-likelihood with the
synchrotron amplitude defined at 23 and at 30 GHz, each relative to its value at the profile optimum. Two panels, one
shared vertical axis.

    python nextgen/scripts/paperA1_fig_scan.py --lcdm nextgen/experiments/E044/scan_lcdm.npz \
        --flat nextgen/experiments/E044/scan_flat.npz --out nextgen/papers/paperA1/figures/fig_scan.pdf
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

INK, MUTED, GRID = "#0b0b0b", "#52514e", "#dddcd6"
C23, C30 = "#2a78d6", "#eb6834"          # categorical slots 1 and 2 (validated: CVD dE 24.7, contrast >= 3:1)
SERIES = [(2, "profile likelihood", INK, "-", "o"), (3, "flat prior, amplitude at 23 GHz", C23, "-", "s"),
          (4, "flat prior, amplitude at 30 GHz", C30, "-", "D")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lcdm", required=True)
    ap.add_argument("--flat", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ymin", type=float, default=-400.0)
    ap.add_argument("--ymax", type=float, default=0.0, help="0: from the data")
    a = ap.parse_args()
    panels = [(r"$\Lambda$CDM prior on the CMB", np.load(a.lcdm)), ("flat prior on the CMB", np.load(a.flat))]
    top = a.ymax
    if top <= 0:
        top = 1.12 * max(float(np.max(z["rows"][z["rows"][:, 0] == 0][:, 3:5])) for _, z in panels)
        top = max(top, 60.0)
    fig, axs = plt.subplots(1, 2, figsize=(6.9, 3.3), sharey=True)
    for ax, (title, z) in zip(axs, panels):
        rows = z["rows"][z["rows"][:, 0] == 0]
        x = rows[:, 1]
        b0 = float(z["profile_opt"][int(z["region"])])
        ax.axhline(0.0, color=MUTED, lw=0.6, zorder=1)
        ax.axvline(b0, color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=1)
        for col, label, color, ls, mk in SERIES:
            y = rows[:, col]
            ax.plot(x, y, ls, color=color, lw=1.6, marker=mk, ms=3.4, markeredgecolor="white", markeredgewidth=0.5, label=label,
                    zorder=3, clip_on=True)
            k = int(np.argmax(y))
            if col > 2 and a.ymin < y[k] < top:              # the largest value on the grid: an interior maximum, or the edge
                ax.plot([x[k]], [y[k]], marker=mk, ms=6.5, color=color, markeredgecolor=INK, markeredgewidth=0.9, zorder=4)
        ax.set_xlim(-4.7, -1.5)
        ax.set_ylim(a.ymin, top)
        ax.set_title(title, fontsize=9, color=INK, loc="left")
        ax.set_xlabel("index of the High Latitude region", fontsize=8.5, color=INK)
        ax.tick_params(labelsize=8, colors=MUTED, length=3)
        ax.grid(axis="y", color=GRID, lw=0.5, zorder=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(MUTED)
            ax.spines[sp].set_linewidth(0.6)
        ax.text(b0, a.ymin + 0.03 * (top - a.ymin), " profile optimum", fontsize=7.5, color=MUTED, ha="left", va="bottom")
    axs[0].set_ylabel(r"$\ln\mathcal{L}$ minus its value at the profile optimum", fontsize=8.5, color=INK)
    axs[0].legend(loc="upper left", fontsize=7.5, frameon=False, handlelength=1.6, labelcolor=INK, borderaxespad=0.2)
    fig.tight_layout(w_pad=1.6)
    fig.savefig(a.out)
    print("wrote", a.out, "ymax", round(top, 1))


if __name__ == "__main__":
    main()
