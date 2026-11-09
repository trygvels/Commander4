"""Part A: n = 1 (redraw then ONE legacy move), small steps. var(beta) vs T(k) = 1/2 + 1/(4 c_k)."""
import sys, time
from math import sqrt
import numpy as np
from sr_lib import *
ks = [int(float(s)) for s in sys.argv[1].split(",")]
eps = float(sys.argv[2])
nch, nit, burn = int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
rng = np.random.default_rng(int(sys.argv[6]))
t0 = time.time()
print(f"# n=1, sig = {eps} * sqrt(T(k)); nch={nch} nit={nit} burn={burn}; profile variance = 1")
print(f"{'k':>7} {'T_pred':>9} {'lead':>9} {'grid':>9} {'sim var':>9} {'se':>7} {'sim/T':>7} {'(sim-T)/se':>10} "
      f"{'sim/lead':>8} {'under':>8} {'over':>8} {'kurt':>6} {'kurtG':>6} {'acc':>6} {'accG':>6}")
for k in ks:
    T = T_small(k)
    sig = eps * sqrt(T)
    g2, gk, gacc, _ = stationary_grid(k, sig)
    r = simulate(k, 1, sig, nch, nit, burn, rng, init_sd=sqrt(T))
    print(f"{k:7d} {T:9.3f} {T_lead(k):9.3f} {g2:9.3f} {r['var']:9.3f} {r['se']:7.3f} {r['var']/T:7.4f} "
          f"{(r['var']-T)/r['se']:10.2f} {r['var']/T_lead(k):8.4f} {r['var_under']:8.3f} {r['var_over']:8.3f} "
          f"{r['kurt']:6.3f} {gk:6.3f} {r['acc']:6.4f} {gacc:6.4f}", flush=True)
print("time", round(time.time() - t0, 1))
