"""Grid (exact kernel) stationary variance for n=1 as the step shrinks: does it converge to T(k)?"""
import time
from math import sqrt
from sr_lib import *
t0 = time.time()
print(f"{'k':>7} {'T_pred':>9} " + " ".join(f"eps={e:<5}" for e in [0.4, 0.2, 0.1, 0.05]))
for k in [10, 100, 1000, 100000]:
    tab = alpha_table(k)
    T = T_small(k)
    row = []
    for e in [0.4, 0.2, 0.1, 0.05]:
        nb = 1201 if e >= 0.1 else 2401
        g2, gk, gacc, _ = stationary_grid(k, e * sqrt(T), nb=nb, tab=tab)
        row.append(g2 / T)
    print(f"{k:7d} {T:9.3f} " + " ".join(f"{x:9.5f}" for x in row), flush=True)
print("time", round(time.time() - t0, 1))
