"""E035: stationary law of the legacy re-solve index move (toy).

chi2(beta, a) = chi2_min(beta) + X. The legacy move proposes beta' and fresh amplitudes, so
X' ~ chi2_k, and accepts on exp(-(chi2' - chi2)/2). The Gibbs sweep refreshes X ~ chi2_k once per
iteration. Ridge posterior N(0, I) in `dim` dimensions (chi2_min = |beta|^2). Reports the sd of
beta in units of the ridge sd against the prediction 0.7916 k^(1/4), and two exact controls.
"""
import numpy as np

COEF = (2.0 * np.pi) ** 0.25 / 2.0


def run(k, numstep, tau, refresh, nchain=4000, burn=100, keep=300, dim=1, seed=0):
    rng = np.random.default_rng([seed, k, numstep, int(1000 * tau), dim, {"chi2": 0, "none": 1, "tilted": 2}[refresh]])
    beta = np.zeros((nchain, dim))
    X = rng.chisquare(k, nchain) if refresh != "tilted" else rng.gamma(k / 2.0, 1.0, nchain)
    out = np.empty((keep, nchain, dim))
    xs = np.empty((keep, nchain))
    nacc = 0
    for it in range(burn + keep):
        if refresh == "chi2":
            X = rng.chisquare(k, nchain)
        elif refresh == "tilted":
            X = rng.gamma(k / 2.0, 1.0, nchain)
        for _ in range(numstep):
            prop = beta + tau * rng.standard_normal((nchain, dim))
            Xp = rng.chisquare(k, nchain)
            loga = -0.5 * ((prop ** 2).sum(1) - (beta ** 2).sum(1)) - 0.5 * (Xp - X)
            acc = np.log(rng.random(nchain)) < loga
            beta[acc] = prop[acc]
            X[acc] = Xp[acc]
            if it >= burn:
                nacc += acc.sum()
        if it >= burn:
            out[it - burn] = beta
            xs[it - burn] = X
    h = keep // 2
    b = out[..., 0]
    sd, sd1, sd2 = b.std(), b[:h].std(), b[h:].std()
    kurt = ((b - b.mean()) ** 4).mean() / sd ** 4 - 3.0
    # Monte Carlo error of sd from the scatter of per-chain variances
    v = b.var(0)
    se = 0.5 * v.std() / np.sqrt(nchain) / sd
    return sd, se, sd1, sd2, kurt, nacc / (keep * numstep * nchain), xs.mean()


if __name__ == "__main__":
    print("# P1: one index step per refresh, proposal 0.3 predicted widths, refresh X ~ chi2_k")
    print(f"{'k':>8s} {'pred':>8s} {'sd':>8s} {'+-':>6s} {'sd/pred':>8s} {'half1':>8s} {'half2':>8s} {'kurt':>6s} {'acc':>6s}")
    for k in (10, 100, 1000, 10000, 100000, 190825):
        pred = COEF * k ** 0.25
        sd, se, s1, s2, ku, acc, _ = run(k, 1, 0.3 * pred, "chi2", burn=300, keep=600)
        print(f"{k:8d} {pred:8.3f} {sd:8.3f} {se:6.3f} {sd / pred:8.3f} {s1:8.3f} {s2:8.3f} {ku:6.2f} {acc:6.3f}", flush=True)
    print("# P2: no refresh (pure pseudo-marginal chain), proposal 2.4 ridge sd; exact for ridge, mean X = k/2")
    print(f"{'k':>8s} {'sd':>8s} {'+-':>6s} {'half1':>8s} {'half2':>8s} {'acc':>7s} {'mean X':>8s} {'k/2':>6s}")
    for k in (4, 8, 16):
        sd, se, s1, s2, ku, acc, mx = run(k, 20, 2.4, "none", nchain=4000, burn=400, keep=1200)
        print(f"{k:8d} {sd:8.3f} {se:6.3f} {s1:8.3f} {s2:8.3f} {acc:7.4f} {mx:8.3f} {k / 2:6.1f}", flush=True)
    print("# P3: refresh from the extended target's conditional X ~ Gamma(k/2, 1); exact for ridge")
    for k in (4, 8, 16):
        sd, se, s1, s2, ku, acc, mx = run(k, 20, 2.4, "tilted", nchain=4000, burn=200, keep=600)
        print(f"{k:8d} {sd:8.3f} {se:6.3f} {s1:8.3f} {s2:8.3f} {acc:7.4f} {mx:8.3f} {k / 2:6.1f}", flush=True)
    print("# P4: dependence on proposal width (in predicted widths) and steps per refresh, refresh X ~ chi2_k")
    print(f"{'k':>8s} {'dim':>4s} {'nstep':>6s} {'tau/W':>6s} {'pred':>8s} {'sd':>8s} {'+-':>6s} {'sd/pred':>8s} {'half1':>8s} {'half2':>8s} {'kurt':>6s} {'acc':>6s}")
    for k in (1000, 190825):
        pred = COEF * k ** 0.25
        for dim in (1, 2):
            for numstep in (1, 5, 20):
                for f in (0.3, 1.0, 2.4):
                    if dim == 2 and f == 0.3:
                        continue
                    sd, se, s1, s2, ku, acc, _ = run(k, numstep, f * pred, "chi2", dim=dim,
                                                     burn=300 if f < 1 else 100, keep=300)
                    print(f"{k:8d} {dim:4d} {numstep:6d} {f:6.1f} {pred:8.3f} {sd:8.3f} {se:6.3f} {sd / pred:8.3f} "
                          f"{s1:8.3f} {s2:8.3f} {ku:6.2f} {acc:6.3f}", flush=True)
