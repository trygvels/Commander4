"""Part D2: legacy move alone, repeat with new seeds; k=64 started from X ~ chi2_k (far) and from
X ~ chi2_k/2 (the predicted stationary law), to separate invariance from convergence."""
import time
from math import sqrt
import numpy as np
from sr_lib import *
rng = np.random.default_rng(72)
t0 = time.time()
print(f"{'k':>4} {'start X':>10} {'nit':>7} {'var':>8} {'se':>7} {'(var-1)/se':>10} {'under':>7} {'over':>7} {'E[X]/k':>7} {'acc':>9}")
for k, x0, nch, nit, burn in [(1, 1.0, 8000, 6000, 1000), (2, 1.0, 8000, 6000, 1000), (8, 1.0, 8000, 6000, 1000),
                              (64, 1.0, 1500, 30000, 10000), (64, 0.5, 1500, 30000, 10000)]:
    r = simulate(k, 1, 1.5, nch, nit, burn, rng, redraw=None, init_sd=1.0, x0_scale=x0)
    print(f"{k:4d} {('chi2_k' if x0 == 1.0 else 'chi2_k/2'):>10} {nit:7d} {r['var']:8.4f} {r['se']:7.4f} {(r['var']-1)/r['se']:10.2f} "
          f"{r['var_under']:7.4f} {r['var_over']:7.4f} {r['meanX']/k:7.4f} {r['acc']:9.2e}", flush=True)
print("time", round(time.time() - t0, 1))
