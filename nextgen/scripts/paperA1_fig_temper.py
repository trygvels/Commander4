"""Paper A1, Fig. 1: temperature of the tempered profile posterior that the re-solve move converges to.

Left: T against k for one index move per amplitude draw: T = 1/2 + 1/(4 c_k), its leading term, and a
toy chain. Right: T_n / T_1 for n moves per draw: the record-process law and the toy.

    python nextgen/scripts/paperA1_fig_temper.py --out nextgen/papers/paperA1/figures/fig_temper.pdf \
        --data nextgen/experiments/E035/fig_temper_toy.npz
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import quad  # noqa: E402
from scipy.special import gammaln  # noqa: E402
from scipy.stats import norm  # noqa: E402


def T_exact(k):
    lc = -(k + 1) * np.log(2.0) + gammaln(k + 1) - 2 * gammaln(k / 2 + 1)
    return 0.5 + 0.25 * np.exp(-lc)


def R_n(n):
    nphi = sum(j * quad(lambda t: norm.pdf(t) ** 2 * norm.sf(t) ** (j - 1), -12, 12)[0] for j in range(1, n + 1))
    hn = sum(1.0 / i for i in range(1, n + 2)) - 1.0
    return hn / (np.sqrt(np.pi) * nphi)


def toy(k, nstep, rng, step_frac=0.15, nchain=4000, burn=400, keep=800):
    """Composite chain for chi2_min = beta^2: refresh X ~ chi2_k, then nstep moves. Returns Var(beta)
    and its standard error from the spread of per-chain second moments."""
    T1 = T_exact(k)
    width = np.sqrt(T1 * (R_n(nstep) if nstep > 1 else 1.0))
    step = step_frac * width
    beta = width * rng.standard_normal(nchain)
    acc2 = np.zeros(nchain)
    for it in range(burn + keep):
        X = rng.chisquare(k, nchain)
        for _ in range(nstep):
            bp = beta + step * rng.standard_normal(nchain)
            Xp = rng.chisquare(k, nchain)
            ok = np.log(rng.random(nchain)) < -0.5 * (bp ** 2 - beta ** 2 + Xp - X)
            beta = np.where(ok, bp, beta)
            X = np.where(ok, Xp, X)
        if it >= burn:
            acc2 += beta ** 2
    m2 = acc2 / keep
    return float(m2.mean()), float(m2.std() / np.sqrt(nchain))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", required=True)
    a = ap.parse_args()
    ks = [4, 10, 30, 100, 300, 1000, 10000, 100000]
    ns = [1, 2, 3, 5, 10, 20, 50, 100]
    if os.path.exists(a.data):
        z = np.load(a.data)
        v1, e1, vn, en = z["v1"], z["e1"], z["vn"], z["en"]
    else:
        rng = np.random.default_rng(3511)
        r1 = [toy(k, 1, rng) for k in ks]
        v1, e1 = np.array([r[0] for r in r1]), np.array([r[1] for r in r1])
        vn, en = np.zeros((2, len(ns))), np.zeros((2, len(ns)))
        for i, k in enumerate((10000, 100000)):
            for j, n in enumerate(ns):
                vn[i, j], en[i, j] = toy(k, n, rng, burn=300, keep=500 if n < 50 else 300)
        np.savez(a.data, ks=ks, ns=ns, v1=v1, e1=e1, vn=vn, en=en)
    kk = np.logspace(0.3, 5.3, 200)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.6))
    ax[0].loglog(kk, T_exact(kk), "k-", lw=1.2, label=r"$T = \frac{1}{2} + 1/(4c_k)$")
    ax[0].loglog(kk, np.sqrt(np.pi * kk / 8), "k--", lw=0.9, label=r"$\sqrt{\pi k/8}$")
    ax[0].errorbar(ks, v1, yerr=3 * e1, fmt="o", ms=4.5, color="C3", label="toy chain")
    ax[0].set_xlabel("amplitudes re-solved per proposal, $k$")
    ax[0].set_ylabel(r"variance / profile posterior variance")
    ax[0].legend(frameon=False, fontsize=9)
    nn = np.arange(1, 101)
    ax[1].semilogx(nn, [R_n(int(n)) for n in nn], "k-", lw=1.2, label="record-process law")
    for i, (k, c, mk) in enumerate(((10000, "C0", "s"), (100000, "C3", "o"))):
        ax[1].errorbar(ns, vn[i] / T_exact(k), yerr=3 * en[i] / T_exact(k), fmt=mk, ms=4.5, color=c,
                       label=f"toy, $k = 10^{int(np.log10(k))}$")
    ax[1].set_xlabel("index moves per amplitude draw, $n$")
    ax[1].set_ylabel("$T_n / T_1$")
    ax[1].set_ylim(0.3, 1.05)
    ax[1].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(a.out)
    print("k      T_exact   toy       se       (toy-T)/se")
    for k, v, e in zip(ks, v1, e1):
        print(f"{k:6d} {T_exact(k):9.3f} {v:9.3f} {e:8.3f} {(v - T_exact(k)) / e:+6.1f}")
    print("n   R_n     toy(1e4)/T1  toy(1e5)/T1")
    for j, n in enumerate(ns):
        print(f"{n:3d} {R_n(n):.4f}  {vn[0, j] / T_exact(10000):.4f} +- {en[0, j] / T_exact(10000):.4f}   "
              f"{vn[1, j] / T_exact(100000):.4f} +- {en[1, j] / T_exact(100000):.4f}")


if __name__ == "__main__":
    main()
