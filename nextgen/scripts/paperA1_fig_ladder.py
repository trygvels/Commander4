"""Paper A1, coverage figure from the E041 ladder.

    python nextgen/scripts/paperA1_fig_ladder.py --data nextgen/experiments/E041/ladder.npz \
        --out nextgen/papers/paperA1/figures/fig_ladder.pdf
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    z = np.load(a.data, allow_pickle=True)
    cols, R = list(z["cols"]), z["rows"]
    c = lambda k: R[:, cols.index(k)]  # noqa: E731
    meas = c("sd_sand_truth") < 0.25
    t = c("t_truth")
    fig, ax = plt.subplots(1, 3, figsize=(12.0, 3.7))
    tt = np.logspace(-4, 0.5, 100)
    ax[0].semilogx(tt, np.sqrt(1 + tt), "k-", lw=1.0, label=r"$\sqrt{1 + T/J}$")
    ax[0].semilogx(tt, (1 + tt) ** 1.5 / np.sqrt(1 + 2 * tt), "k:", lw=1.0, label=r"$(1+t)^{3/2}/(1+2t)^{1/2}$")
    ax[0].axhline(1.0, color="0.6", lw=0.8)
    for key, mk, col, lab in (("p_sdz_curv", "o", "C0", "curvature width"), ("p_sdz_plain", "^", "C1", "plug-in sandwich"),
                              ("p_sdz_deb", "s", "C3", "debiased sandwich")):
        ax[0].plot(t[meas], c(key)[meas], mk, color=col, ms=4.5, label=lab)
    ax[0].set_xlabel(r"$T/J$ at the truth")
    ax[0].set_ylabel(r"sd of $(\hat\beta - \beta)/\sigma$")
    ax[0].set_ylim(0.85, 2.2)
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].axhline(0.683, color="0.6", lw=0.8)
    for key, mk, col in (("p_cov_curv", "o", "C0"), ("p_cov_plain", "^", "C1"), ("p_cov_deb", "s", "C3")):
        ax[1].semilogx(t[meas], c(key)[meas], mk, color=col, ms=4.5)
    ax[1].semilogx(t[meas], c("z_cov")[meas], "x", color="C2", ms=5, label="zero-mean prior")
    ax[1].fill_between([1e-4, 3], 0.683 - 2 * 0.0233, 0.683 + 2 * 0.0233, color="0.9", zorder=0)
    ax[1].set_xlim(1e-4, 3)
    ax[1].set_xlabel(r"$T/J$ at the truth")
    ax[1].set_ylabel(r"coverage of $\pm 1\sigma$")
    ax[1].legend(frameon=False, fontsize=8)
    for key, pk, mk, col, lab in (("f23_shift", "f23_pred", "v", "C0", "reference 23 GHz"), ("f30_shift", "f30_pred", "^", "C3", "reference 30 GHz")):
        ax[2].loglog(np.abs(c(pk)), np.abs(c(key)), mk, color=col, ms=4.5, label=lab)
    xx = np.logspace(-4, 2, 10)
    ax[2].loglog(xx, xx, "k-", lw=0.8)
    ax[2].axhline(1.7, color="0.6", lw=0.8, ls="--")
    ax[2].set_xlim(1e-4, 100)
    ax[2].set_ylim(1e-4, 5)
    ax[2].set_xlabel(r"predicted shift, $|\sum_p \ln(\nu_{\rm ref}/\bar\nu_p)|/J$")
    ax[2].set_ylabel("measured shift of the flat-prior optimum")
    ax[2].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out)


if __name__ == "__main__":
    main()
