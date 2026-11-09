"""Part D: (i) legacy move alone (no redraw): prediction var(beta) = 1 (profile), E[X] = k/2.
(ii) redraw at HALF covariance (X ~ chi2_k / 2) + legacy moves: prediction var(beta) = 1 exactly.
(iii) negative control: redraw at full covariance: T(k) > 1.
Acceptance bound at stationarity for (i),(ii) with beta frozen: 2 I_{1/3}(k/2, k/2)."""
import time
from math import sqrt
import numpy as np
from scipy.special import betainc
from sr_lib import *
rng = np.random.default_rng(71)
t0 = time.time()
sig = 1.5
print(f"{'k':>4} {'mode':>22} {'var':>8} {'se':>7} {'(var-1)/se':>10} {'under':>7} {'over':>7} {'E[X]/k':>7} {'acc':>9} {'2I_1/3':>9} {'T(k)':>7}")
for k in [1, 4, 16, 64]:
    bound = 2 * betainc(k / 2, k / 2, 1 / 3)
    tau = 4.0 / bound
    burn = int(15 * tau) + 200
    nit = burn + int(40 * tau) + 2000
    nch = int(min(8000, max(1500, 6e7 / nit)))
    for mode, redraw, n in [("legacy alone", None, 1), ("half-cov redraw, n=1", 0.5, 1),
                            ("half-cov redraw, n=5", 0.5, 5), ("full redraw (control)", 1.0, 1)]:
        nit_m = nit if n == 1 else max(burn // 3 + 600, nit // 4)
        burn_m = burn if n == 1 else burn // 3
        if redraw == 1.0:
            nit_m, burn_m = 1500, 300
        r = simulate(k, n, sig, nch, nit_m, burn_m, rng, redraw=redraw, init_sd=1.0)
        print(f"{k:4d} {mode:>22} {r['var']:8.4f} {r['se']:7.4f} {(r['var']-1)/r['se']:10.2f} {r['var_under']:7.4f} "
              f"{r['var_over']:7.4f} {r['meanX']/k:7.4f} {r['acc']:9.2e} {bound:9.2e} {T_small(k):7.3f}", flush=True)
print("acceptance bound 2 I_{1/3}(k/2,k/2) at larger k:")
for k in [100, 300, 1000, 3000, 10000]:
    print(f"  k={k:6d}  {2*betainc(k/2, k/2, 1/3):.3e}   (8/9)^(k/2) = {(8/9)**(k/2):.3e}")
print("time", round(time.time() - t0, 1))
