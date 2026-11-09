import time
from math import sqrt, pi
import numpy as np
from sr_lib import *
t0 = time.time()
print("check closed form c_k against FFT quadrature, and T(k)")
print(f"{'k':>8} {'c_k closed':>12} {'c_k FFT':>12} {'alpha(0) FFT':>12} {'1/2+c_k':>10} {'T=1/2+1/(4c)':>13} {'sqrt(pi k/8)':>13}")
for k in [3, 10, 30, 100, 1000, 10000, 100000, 200000]:
    z, alp, G2 = alpha_table(k)
    i0 = np.argmin(np.abs(z))
    print(f"{k:8d} {c_k(k):12.6e} {G2[i0]:12.6e} {alp[i0]:12.6f} {0.5+c_k(k):10.6f} {T_small(k):13.4f} {T_lead(k):13.4f}")
print("records: n, N_acc, N_phi, Var/sqrt(k), ratio to n=1, mean accept fraction")
v1 = var_scaled_n(1)
for n in [1, 2, 3, 5, 10, 20, 30, 50, 100, 300, 1000]:
    Na, Np = record_sums(n)
    print(f"{n:5d} {Na:8.4f} {Np:8.4f} {var_scaled_n(n):8.4f} {var_scaled_n(n)/v1:8.4f} {Na/n:8.4f}")
print("sqrt(pi/8) =", sqrt(pi/8))
print("time", time.time() - t0)
