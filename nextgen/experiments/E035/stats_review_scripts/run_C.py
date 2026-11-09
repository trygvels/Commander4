"""Part C: n legacy moves per amplitude redraw. Small-step prediction (record process):
Var(beta) = sqrt(k) * N_acc(n) / (2 sqrt2 N_phi(n)), acceptance fraction N_acc(n)/n."""
import sys, time
from math import sqrt
import numpy as np
from sr_lib import *
k = int(float(sys.argv[1]))
eps = float(sys.argv[2])
ns = [int(s) for s in sys.argv[3].split(",")]
rng = np.random.default_rng(int(sys.argv[4]))
budget = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0e8
t0 = time.time()
T1 = T_small(k)
print(f"# k={k}, sig = {eps} * sqrt(V_n), V_n = leading-order prediction; T(k) n=1 exact = {T1:.3f}, sqrt(pi k/8) = {T_lead(k):.3f}")
print(f"{'n':>5} {'V_lead':>9} {'sim var':>9} {'se':>7} {'sim/V_lead':>10} {'sim/T1':>8} {'pred R_n':>8} {'kurt':>6} "
      f"{'acc':>7} {'acc_pred':>8} {'under':>8} {'over':>8} {'nch':>6} {'nit':>6}")
for n in ns:
    Na, Np = record_sums(n)
    V = sqrt(k) * Na / (2 * sqrt(2) * Np)
    sig = eps * sqrt(V)
    tau = max(8.0, 1.0 / (Na * eps * eps))
    burn = int(20 * tau) + 20
    nit = burn + int(60 * tau) + 100
    nch = int(min(10000, max(1000, budget / (n * nit))))
    r = simulate(k, n, sig, nch, nit, burn, rng, init_sd=sqrt(V))
    print(f"{n:5d} {V:9.3f} {r['var']:9.3f} {r['se']:7.3f} {r['var']/V:10.4f} {r['var']/T1:8.4f} "
          f"{var_scaled_n(n)/var_scaled_n(1):8.4f} {r['kurt']:6.3f} {r['acc']:7.4f} {Na/n:8.4f} "
          f"{r['var_under']:8.3f} {r['var_over']:8.3f} {nch:6d} {nit:6d}", flush=True)
print("time", round(time.time() - t0, 1))
