"""Paper A1: the picture behind the block errors. For every nside-2 patch of sky that carries information on the
index of a region, the score of the patch in units of the error it has under the model,

    z_B = -G_B / sqrt(D_B)        (positive: the patch alone would move to a flatter index),

against the share h_B of the region's information that the patch carries. Under the model z_B scatters by about
one. The squared ratio of the block error (own scores, CR3) to the model error is the information-weighted mean
of z_B^2 / (1 - h_B)^2, so the picture shows what that ratio is made of, for every patch. Small multiples, one
panel per region, one shared vertical axis.

(Until the third referee read the vertical axis was the one-step index of the patch with its error, and patches
with an error above 0.12 were not drawn, which left 4 of the 38 High Latitude patches.)

    python nextgen/scripts/paperA1_fig_blocks.py --data nextgen/experiments/E044/lcdm_all.npz \
        --out nextgen/papers/paperA1/figures/fig_blocks.pdf
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NAMES = ["High Latitude", "Spur", "Intermediate", "Plane"]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#dddcd6"
POINT = "#2a78d6"                       # one series, one hue
BAND = "#e6e5df"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ylim", type=float, default=4.6)
    a = ap.parse_args()
    z = np.load(a.data)
    if "blk_G" in z.files:                                   # E044 output: block errors from each region's own scores
        G, D, beta, sd, cr3, neff = z["blk_G"], z["blk_D"], z["beta"], z["sd_sand"], z["blk_CR3_own"], z["blk_neff"]
    else:
        G, D, beta, sd, cr3, neff = z["G"], z["D"], z["beta"], z["sd_model"], z["CR3"], z["neff"]
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 5.4), sharey=True, sharex=True)
    axes = axs.ravel()
    for r, ax in enumerate(axes):
        d, g = D[:, r], G[:, r]
        pos = d > 0
        share = d[pos] / d[pos].sum()
        zs = -g[pos] / np.sqrt(d[pos])
        ax.axhspan(-1.0, 1.0, color=BAND, lw=0, zorder=0)
        ax.axhline(0.0, color=MUTED, lw=0.6, zorder=1)
        for y in (-2.0, 2.0):
            ax.axhline(y, color=MUTED, lw=0.5, ls=(0, (2, 2)), zorder=1)
        inside = np.abs(zs) <= a.ylim
        ax.plot(100 * share[inside], zs[inside], "o", ms=4.6, color=POINT, markeredgecolor="white", markeredgewidth=0.6, zorder=3, clip_on=False)
        for xx, zz in zip(100 * share[~inside], zs[~inside]):  # a value outside the frame: an arrow at the edge, with the value
            ax.annotate("", xy=(xx, np.sign(zz) * a.ylim * 0.99), xytext=(xx, np.sign(zz) * a.ylim * 0.86),
                        arrowprops=dict(arrowstyle="-|>", color=POINT, lw=1.0))
            ax.text(xx * 1.12, np.sign(zz) * a.ylim * 0.90, f"{zz:+.1f}".replace("-", "−"), fontsize=7, color=INK, va="center")
        ax.set_xscale("log")
        ax.set_xlim(0.02, 60)
        ax.set_ylim(-a.ylim, a.ylim)
        ax.set_xticks([0.1, 1, 10])
        ax.set_xticklabels(["0.1", "1", "10"])
        ax.set_yticks([-4, -2, 0, 2, 4])
        ax.set_title(f"{NAMES[r]}   {beta[r]:.2f}".replace("-", "−"), fontsize=9, color=INK, loc="left", pad=13)
        ratio = cr3[r] / sd[r]
        note = f"{int(pos.sum())} patches, effective number {neff[r]:.1f}; block error {ratio:.2f} times the model error"
        ax.text(0.0, 1.02, note, transform=ax.transAxes, fontsize=7.2, color=MUTED, ha="left", va="bottom")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
            ax.spines[s].set_linewidth(0.6)
        ax.tick_params(colors=MUTED, labelsize=7.5, length=3)
    for ax in axs[1]:
        ax.set_xlabel("share of the region's information (percent)", fontsize=8.5, color=INK)
    fig.supylabel("score of the patch in units of its error under the model", fontsize=8.5, color=INK)
    axes[0].text(0.025, 1.0, "one standard deviation\nunder the model", fontsize=7.2, color=MUTED, ha="left", va="bottom")
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote", a.out)
    for r in range(4):
        d, g = D[:, r], G[:, r]
        pos = d > 0
        h = d[pos] / d[pos].sum()
        zz = -g[pos] / np.sqrt(d[pos])
        print(f"{NAMES[r]:16s} patches with information {int(pos.sum())}; rms z {np.sqrt(np.mean(zz ** 2)):.2f}; largest |z| {np.abs(zz).max():.1f}; "
              f"sqrt of the weighted mean of z^2 / (1 - h)^2 = {np.sqrt(np.sum(h * zz ** 2 / (1 - h) ** 2)):.2f} (block over model error {cr3[r] / sd[r]:.2f})")


if __name__ == "__main__":
    main()
