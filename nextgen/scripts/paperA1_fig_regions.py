"""Paper A1: the regions and the processing mask of the real-data demonstration (`lowres.regions.named_regions`).

    python -I nextgen/scripts/paperA1_fig_regions.py --out nextgen/papers/paperA1/figures/fig_regions.pdf
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import healpy as hp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres.regions import NAMES, TAU_A, named_regions  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]          # validated categorical slots 1 to 4, fixed order
MASKED = "#8f8e88"
INK = "#0b0b0b"
LABEL_AT = {0: (120.0, 55.0), 1: (25.0, 48.0), 2: (-110.0, 20.0), 3: (100.0, 0.0)}       # (l, b) of the direct labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    region, keep, _ = named_regions(16)
    npix = region.size // 2
    m = np.where(keep[:npix], region[:npix], 4).astype(float)
    fig = plt.figure(figsize=(6.6, 3.7))
    hp.mollview(m, nest=True, fig=fig.number, cmap=ListedColormap(COLORS + [MASKED]), min=-0.5, max=4.5, cbar=False, title="", notext=True,
                badcolor="white", bgcolor="white")
    hp.graticule(dpar=30, dmer=60, color="white", alpha=0.35, lw=0.4, verbose=False)
    for k, (lon, lat) in LABEL_AT.items():
        n_el = int(((region == k) & keep).sum())
        hp.projtext(lon, lat, f"{NAMES[k]}\n{n_el} elements", lonlat=True, fontsize=7.5, color=INK, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85))
    hp.projtext(TAU_A[0] - 12.0, TAU_A[1] - 14.0, "Tau A", lonlat=True, fontsize=7, color=INK, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    hp.projtext(-35.0, -9.0, "masked", lonlat=True, fontsize=7, color=INK, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    fig.savefig(a.out, bbox_inches="tight", pad_inches=0.05)
    print("wrote", a.out, "; kept elements per region", [int(((region == k) & keep).sum()) for k in range(4)], "; masked pixels", int((~keep[:npix]).sum()))


if __name__ == "__main__":
    main()
