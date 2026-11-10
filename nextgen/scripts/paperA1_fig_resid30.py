"""Paper A1: maps of observed minus predicted for PR4 30 GHz and for the Ka control (E042 g).

    python nextgen/scripts/paperA1_fig_resid30.py --dir nextgen/experiments/E042 --out nextgen/papers/paperA1/figures/fig_resid30.pdf
"""
import argparse

import healpy as hp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="nextgen/experiments/E042")
    ap.add_argument("--out", required=True)
    ap.add_argument("--vmax", type=float, default=6.0)
    a = ap.parse_args()
    z30 = np.load(f"{a.dir}/predict30.npz")
    zka = np.load(f"{a.dir}/predictKa.npz")
    n = z30["D"].size // 2
    fig = plt.figure(figsize=(11.0, 6.4))
    panels = [(z30["target"][:n], "30 GHz, observed $Q$", 30.0), (z30["D"][:n], "30 GHz, observed $-$ predicted, $Q$", a.vmax),
              (zka["D"][:n], "Ka, observed $-$ predicted, $Q$", a.vmax),
              (z30["target"][n:], "30 GHz, observed $U$", 30.0), (z30["D"][n:], "30 GHz, observed $-$ predicted, $U$", a.vmax),
              (zka["D"][n:], "Ka, observed $-$ predicted, $U$", a.vmax)]
    for i, (m, title, vm) in enumerate(panels):
        hp.mollview(hp.reorder(m, n2r=True), fig=fig.number, sub=(2, 3, i + 1), title=title, min=-vm, max=vm, cmap="RdBu_r",
                    unit=r"$\mu$K", cbar=True, notext=True, margins=(0.005, 0.02, 0.005, 0.04))
    fig.savefig(a.out, dpi=150)
    for name, z in (("30 GHz", z30), ("Ka", zka)):
        D = z["D"]
        print(f"{name}: rms of D in Q {np.sqrt(np.mean(D[:n] ** 2)):.2f}, in U {np.sqrt(np.mean(D[n:] ** 2)):.2f} uK; "
              f"rms at abs(b) > 30: {np.sqrt(np.mean(D[np.tile(np.abs(hp.pix2ang(16, np.arange(n), nest=True, lonlat=True)[1]) > 30, 2)] ** 2)):.2f} uK")


if __name__ == "__main__":
    main()
