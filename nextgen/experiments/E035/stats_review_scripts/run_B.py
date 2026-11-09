"""Part B: n = 1, step size from small to much larger than the stationary width.
eps = sig / sqrt(T(k)). Reports var/T from the exact-kernel grid and from simulation."""
import sys, time
from math import sqrt
import numpy as np
from sr_lib import *
ks = [int(float(s)) for s in sys.argv[1].split(",")]
rng = np.random.default_rng(int(sys.argv[2]))
t0 = time.time()
print(f"{'k':>7} {'eps':>6} {'T_pred':>9} {'grid/T':>8} {'sim/T':>8} {'se/T':>7} {'kurtG':>6} {'kurt':>6} {'accG':>7} {'acc':>7}")
for k in ks:
    T = T_small(k)
    tab = alpha_table(k)
    for eps in [0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]:
        sig = eps * sqrt(T)
        g2, gk, gacc, _ = stationary_grid(k, sig, nb=1401, Lfac=8.0, tab=tab)
        tau = max(1.0, 4.0 / (eps * eps * 0.5), 2.0 / max(gacc, 1e-3))
        burn = int(20 * tau) + 50
        nit = burn + 1500
        r = simulate(k, 1, sig, 6000, nit, burn, rng, init_sd=sqrt(T))
        print(f"{k:7d} {eps:6.2f} {T:9.3f} {g2/T:8.4f} {r['var']/T:8.4f} {r['se']/T:7.4f} {gk:6.3f} {r['kurt']:6.3f} "
              f"{gacc:7.4f} {r['acc']:7.4f}", flush=True)
print("time", round(time.time() - t0, 1))
