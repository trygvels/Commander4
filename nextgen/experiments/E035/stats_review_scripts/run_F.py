"""Part F: E014-like settings in my own 2D toy. chi2_min = bs^2 + bd^2 (ridge sd 1 each), numstep 20,
recorded proposal sds 5.6 (beta_s) and 93 (beta_d) ridge sd, joint proposal, k = 190825.
The acceptance-noise variance is scaled by f (X ~ chi2_{f k}): f = 1 is the exact-draw toy.
Measured in E014: width ratios 17.1 (beta_s, ESS 74) and 16.5 (beta_d, ESS 1184), acceptance 0.073."""
import sys, time
from math import sqrt
import numpy as np
from sr_lib import T_small

def sim2(keff, n, sig, nch, nit, burn, rng):
    sig = np.asarray(sig)
    beta = np.zeros((nch, 2))
    beta[nch // 2:] = 30.0 * rng.standard_normal((nch - nch // 2, 2))   # over-dispersed half
    s2 = np.zeros((nch, 2)); nacc = 0; first = 0
    for it in range(nit):
        X = rng.chisquare(keff, nch)
        for j in range(n):
            bp = beta + sig * rng.standard_normal((nch, 2))
            Xp = rng.chisquare(keff, nch)
            la = -0.5 * ((bp ** 2).sum(1) - (beta ** 2).sum(1) + Xp - X)
            acc = np.log(rng.random(nch)) < la
            beta[acc] = bp[acc]; X[acc] = Xp[acc]
            if it >= burn:
                nacc += int(acc.sum())
        if it >= burn:
            s2 += beta ** 2
    m = nit - burn
    v = s2 / m
    h = nch // 2
    sd = np.sqrt(v.mean(0)); se = 0.5 * v.std(0, ddof=1) / sqrt(nch) / sd
    return sd, se, np.sqrt(v[:h].mean(0)), np.sqrt(v[h:].mean(0)), nacc / (m * n * nch)

k = 190825
rng = np.random.default_rng(int(sys.argv[1]))
t0 = time.time()
print(f"# one-move small-step law at f=1: sqrt(T) = {sqrt(T_small(k)):.3f}")
print(f"{'f':>5} {'sd_s':>8} {'+-':>6} {'sd_d':>8} {'+-':>6} {'acc':>7}   under-start sd (s,d)   over-start sd (s,d)")
for f in [float(x) for x in sys.argv[2].split(",")]:
    sd, se, lo, hi, acc = sim2(f * k, 20, (5.6, 93.0), 2500, 650, 250, rng)
    print(f"{f:5.2f} {sd[0]:8.3f} {se[0]:6.3f} {sd[1]:8.3f} {se[1]:6.3f} {acc:7.4f}   {lo[0]:.3f} {lo[1]:.3f}   {hi[0]:.3f} {hi[1]:.3f}", flush=True)
print("time", round(time.time() - t0, 1))
