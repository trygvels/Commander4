"""Part E: table of Var(beta)/T1(k) for n moves per redraw and step sig = eps*sqrt(T1(k))."""
import sys, time
from math import sqrt
import numpy as np
from sr_lib import *
k = int(float(sys.argv[1]))
ns = [int(s) for s in sys.argv[2].split(",")]
epss = [float(s) for s in sys.argv[3].split(",")]
rng = np.random.default_rng(int(sys.argv[4]))
budget = float(sys.argv[5])
t0 = time.time()
T1 = T_small(k)
print(f"# k={k}, T1 = {T1:.3f}; entries: Var/T1 (se) [acceptance]; small-step prediction R_n in last column")
print(f"{'n':>4} " + " ".join(f"{'eps='+str(e):>24}" for e in epss) + f" {'R_n pred':>9}")
for n in ns:
    Na, Np = record_sums(n)
    row = []
    for eps in epss:
        sig = eps * sqrt(T1)
        tau = max(8.0, 1.0 / (Na * eps * eps))
        burn = int(20 * tau) + 20
        nit = burn + int(60 * tau) + 100
        nch = int(min(10000, max(600, budget / (n * nit))))
        r = simulate(k, n, sig, nch, nit, burn, rng, init_sd=sqrt(T1))
        row.append(f"{r['var']/T1:7.4f} ({r['se']/T1:.4f}) [{r['acc']:.3f}]")
    print(f"{n:4d} " + " ".join(f"{s:>24}" for s in row) + f" {var_scaled_n(n)/var_scaled_n(1):9.4f}", flush=True)
print("time", round(time.time() - t0, 1))
