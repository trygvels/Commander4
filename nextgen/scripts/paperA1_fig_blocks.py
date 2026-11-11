"""Paper A1: the picture behind the block errors. For every nside-2 patch of sky, the index that the patch
alone would move to in one Newton step from the region's fitted index, delta_B = -G_B / D_B, with the error
1 / sqrt(D_B) it would have under the model, against the share of the region's information that the patch
carries. Small multiples, one panel per region, one shared vertical axis.

    python nextgen/scripts/paperA1_fig_blocks.py --data nextgen/experiments/E040/block_v3_maps_lcdm.npz \
        --out nextgen/papers/paperA1/figures/fig_blocks.pdf
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

NAMES = ["High Latitude", "Spur", "Intermediate", "Plane"]
INK, MUTED = "#0b0b0b", "#52514e"
POINT = "#2a78d6"                       # one series, one hue
BAND_MODEL, BAND_BLOCK = "#b9b8b2", "#e6e5df"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-err", type=float, default=0.12, help="patches with a larger one-step error are not drawn")
    ap.add_argument("--ylim", type=float, default=0.42)
    a = ap.parse_args()
    z = np.load(a.data)
    if "blk_G" in z.files:                                   # E044 output: block errors from each region's own scores
        G, D, beta, sd, cr3, neff = z["blk_G"], z["blk_D"], z["beta"], z["sd_sand"], z["blk_CR3_own"], z["blk_neff"]
    else:
        G, D, beta, sd, cr3, neff = z["G"], z["D"], z["beta"], z["sd_model"], z["CR3"], z["neff"]
    fig, axs = plt.subplots(2, 2, figsize=(6.6, 5.6), sharey=True, sharex=True)
    axes = axs.ravel()
    for r, ax in enumerate(axes):
        d, g = D[:, r], G[:, r]
        pos = d > 0
        share = np.where(pos, d, 0.0) / d[pos].sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            delta, err = -g / d, 1.0 / np.sqrt(d)
        sel = pos & (err < a.max_err)
        x = 100 * share[sel]
        if np.isfinite(cr3[r]) and neff[r] >= 4.0:
            ax.axhspan(-cr3[r], cr3[r], color=BAND_BLOCK, lw=0, zorder=0)
        ax.axhspan(-sd[r], sd[r], color=BAND_MODEL, lw=0, zorder=1)
        ax.axhline(0.0, color=MUTED, lw=0.6, zorder=2)
        ax.errorbar(x, delta[sel], yerr=err[sel], fmt="o", ms=4.0, color=POINT, ecolor=POINT, elinewidth=1.0, capsize=0,
                    markeredgecolor="white", markeredgewidth=0.6, zorder=3, clip_on=True)
        out_of_range = np.abs(delta[sel]) > a.ylim
        for xx, dd in zip(x[out_of_range], delta[sel][out_of_range]):
            ax.annotate("", xy=(xx, np.sign(dd) * a.ylim * 0.98), xytext=(xx, np.sign(dd) * a.ylim * 0.86),
                        arrowprops=dict(arrowstyle="-|>", color=POINT, lw=1.0))
        ax.set_xscale("log")
        ax.set_xlim(0.5, 60)
        ax.set_ylim(-a.ylim, a.ylim)
        ax.set_xticks([1, 3, 10, 30])
        ax.set_xticklabels(["1", "3", "10", "30"])
        ax.set_title(f"{NAMES[r]}   {beta[r]:.2f}".replace("-", "\u2212"), fontsize=9, color=INK, loc="left", pad=13)
        ax.text(0.0, 1.02, f"{int(sel.sum())} of {int(pos.sum())} patches drawn, effective number {neff[r]:.1f}", transform=ax.transAxes,
                fontsize=7.5, color=MUTED, ha="left", va="bottom")
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelsize=7.5)
    for ax in axs[1]:
        ax.set_xlabel("share of the region's information (percent)", fontsize=8.5, color=INK)
    fig.supylabel("index of the patch minus index of the region", fontsize=8.5, color=INK)
    # direct labels for the two bands, once, in the first panel
    axes[0].annotate("error if the region\nhas one index", xy=(45, sd[0]), xytext=(6.0, 0.33), fontsize=7.5, color=INK,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6), va="center")
    axes[0].annotate("block error", xy=(45, -cr3[0] * 0.92), xytext=(10.0, -0.30), fontsize=7.5, color=INK,
                     arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6), va="center")
    fig.tight_layout()
    fig.savefig(a.out)
    print("wrote", a.out)
    for r in range(4):
        d, g = D[:, r], G[:, r]
        pos = d > 0
        zz = (-g[pos] / d[pos]) * np.sqrt(d[pos])
        print(f"{NAMES[r]:16s} patches with information {int(pos.sum())}; sum of (delta / error)^2 = {np.sum(zz ** 2):.1f}; "
              f"largest shares {np.sort(d[pos] / d[pos].sum())[::-1][:3].round(2).tolist()}")


if __name__ == "__main__":
    main()
