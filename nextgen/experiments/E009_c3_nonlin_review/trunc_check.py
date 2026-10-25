"""Check: MH with a resample-until-inside-bounds Gaussian random-walk proposal (as in
sample_specind_alm lines 526-548) and plain min(1, pi(y)/pi(x)) acceptance. Claim: the
stationary density is pi(x) Z(x), Z(x) = P(x + s*xi inside bounds), not pi(x).
Toy: pi = N(0,1) truncated to [-1, 2], step s. Compare histogram with pi and with pi*Z."""
import numpy as np
from math import erf, sqrt
rng = np.random.default_rng(1)
lo, hi = -1.0, 2.0
def logpi(x): return -0.5*x*x
def Phi(z): return 0.5*(1+erf(z/sqrt(2)))
def Z(x, s): return Phi((hi-x)/s) - Phi((lo-x)/s)

def run(s, n=2_000_000):
    x = 0.3; out = np.empty(n)
    for i in range(n):
        q = 0
        while True:
            q += 1
            y = x + s*rng.standard_normal()
            if lo < y < hi: break
            if q >= 1000: y = x; break
        if np.log(rng.random()) < logpi(y) - logpi(x): x = y
        out[i] = x
    return out

edges = np.linspace(lo, hi, 31); c = 0.5*(edges[1:]+edges[:-1])
for s in (1.5, 0.3):
    xs = run(s)
    h, _ = np.histogram(xs, edges, density=True)
    p = np.exp(logpi(c)); p /= np.trapezoid(p, c)
    pz = np.exp(logpi(c))*np.array([Z(v, s) for v in c]); pz /= np.trapezoid(pz, c)
    print(f"s={s}: max |hist-pi|/pi = {np.max(np.abs(h-p)/p):.3f}, "
          f"max |hist-pi*Z|/(pi*Z) = {np.max(np.abs(h-pz)/pz):.3f}, "
          f"Z at bound {Z(lo+1e-9,s):.3f}")
# BP10 numbers: prior sigma 0.1 around -3.3, bounds [-4.5,-1.5], step ~ <0.1
print("BP10-like Z at -3.3 for s=0.1:", Z(-3.3, 0.1), " at -4.4:", Z(-4.4, 0.1))
