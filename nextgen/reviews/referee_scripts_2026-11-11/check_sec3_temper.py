"""Referee check of Sect. 3 of paper A1 (tempering law of the re-solve move).

Checks, all independent of the authors' scripts:
  (a) c_k = E[exp(-Z/2); Z>0] for Z = X'-X, X,X' iid chi2_k, by quadrature against the closed form
      2^-(k+1) C(k, k/2); T = 1/2 + 1/(4 c_k) against the table; large-k form.
  (b) R_n = (H_{n+1}-1) / (sqrt(pi) N_phi(n)) by quadrature; R_2 = 5/6; asymptotic exponent.
  (c) acceptance of the half-covariance repair, 2 I_{1/3}(k/2,k/2), and its decay rate.
  (d) a small toy chain (k = 1000, one move per redraw, small step) for Var(beta)/T.
  (e) the small-step value of sqrt(R_20 T_1) at k = 190825 (what Eq. Rn gives for the E014 settings).
Runs in a few seconds on one thread.
"""
import numpy as np
from scipy import special, stats, integrate

print("(a) c_k and T(k)")
for k in [1, 2, 10, 100, 1000, 10_000, 100_000, 200_000, 190_825]:
    lg = special.gammaln(k + 1) - 2 * special.gammaln(k / 2 + 1) - (k + 1) * np.log(2)
    ck = np.exp(lg)
    T = 0.5 + 1 / (4 * ck)
    lead = np.sqrt(np.pi * k / 8)
    # quadrature for small k only (direct double integral reduces to 1D)
    if k <= 1000:
        # c_k = int f(x') e^{-x'/2} [int_0^{x'} f(x) e^{x/2} dx] dx'
        f = stats.chi2(k)
        inner = lambda xp: integrate.quad(lambda x: np.exp(f.logpdf(x) + x / 2 - xp / 2), 0, xp)[0]
        lo, hi = max(0, k - 12 * np.sqrt(2 * k) - 20), k + 12 * np.sqrt(2 * k) + 40
        cq = integrate.quad(lambda xp: f.pdf(xp) * inner(xp), lo, hi, limit=200)[0]
    else:
        cq = np.nan
    print(f"  k={k:7d} c_k={ck:.6e} quad={cq:.6e} T={T:9.3f} lead={lead:9.3f} T-lead-0.5={T-lead-0.5:+.4f}"
          f" acc=0.5+c_k={0.5+ck:.3f} sqrtT/k^0.25={np.sqrt(T)/k**0.25:.4f}")

print("(b) R_n")
t = np.linspace(-12, 12, 480_001)
phi = stats.norm.pdf(t)
sf = stats.norm.sf(t)
def Nphi(n):
    s = 0.0
    for j in range(1, n + 1):
        s += j * np.trapezoid(phi**2 * sf ** (j - 1), t)
    return s
H = lambda n: np.sum(1.0 / np.arange(1, n + 1))
R = {}
for n in [1, 2, 3, 5, 10, 20, 50, 100, 1000]:
    R[n] = (H(n + 1) - 1) / (np.sqrt(np.pi) * Nphi(n))
    print(f"  n={n:5d} R_n={R[n]:.4f}  sqrt(R_n)={np.sqrt(R[n]):.4f}")
print(f"  R_2 - 5/6 = {R[2]-5/6:+.2e}")
print(f"  local exponent of R_n in ln n between n=100 and 1000: "
      f"{np.log(R[1000]/R[100])/np.log(np.log(1000)/np.log(100)):+.3f} (paper: width as (ln n)^-1/4, so R as -1/2)")

print("(c) half-covariance repair")
for k in [4, 8, 16, 64, 100, 1000]:
    a = 2 * special.betainc(k / 2, k / 2, 1 / 3)
    print(f"  k={k:5d} 2 I_1/3(k/2,k/2) = {a:.3e}")
k1, k2 = 2000, 4000
r = (np.log(2 * special.betainc(k2 / 2, k2 / 2, 1 / 3)) - np.log(2 * special.betainc(k1 / 2, k1 / 2, 1 / 3))) / (k2 - k1)
print(f"  decay rate d ln(acc)/dk = {r:.5f}; 0.5 ln(8/9) = {0.5*np.log(8/9):.5f}")

print("(d) toy chain, chi2_min = beta^2, k = 1000, one move per redraw, step 0.1 sqrt(T)")
rng = np.random.default_rng(20261111)
k = 1000
ck = np.exp(special.gammaln(k + 1) - 2 * special.gammaln(k / 2 + 1) - (k + 1) * np.log(2))
T = 0.5 + 1 / (4 * ck)
nch, nburn, nkeep = 4000, 3000, 6000
beta = rng.normal(0, np.sqrt(T), nch)  # start at the predicted law; chi2_min=beta^2 -> profile var 1
s1 = np.zeros(nch); s2 = np.zeros(nch); acc = 0
step = 0.1 * np.sqrt(T)
for it in range(nburn + nkeep):
    X = rng.chisquare(k, nch)
    bp = beta + step * rng.normal(size=nch)
    Xp = rng.chisquare(k, nch)
    logr = -0.5 * (bp**2 + Xp - beta**2 - X)
    a = np.log(rng.random(nch)) < logr
    beta = np.where(a, bp, beta)
    if it >= nburn:
        s1 += beta; s2 += beta**2; acc += a.mean()
var_ch = s2 / nkeep - (s1 / nkeep) ** 2
tot = s2.sum() / (nch * nkeep) - (s1.sum() / (nch * nkeep)) ** 2
se = np.std(s2 / nkeep) / np.sqrt(nch)
print(f"  T={T:.3f}  Var(beta)/profile var = {tot:.3f} +- {se:.3f}  acceptance={acc/nkeep:.3f} (0.5+c_k={0.5+ck:.3f})")

print("(e) what Eq. (Rn) gives for the real-map settings (k=190825, n=20, small steps)")
k = 190825
T1 = 0.5 + np.sqrt(np.pi * k / 8)
print(f"  T_1={T1:.1f}, sqrt(T_1)={np.sqrt(T1):.2f}, sqrt(R_20 T_1)={np.sqrt(R[20]*T1):.2f}; "
      f"small-step acceptance per proposal (H_21-1)/20={(H(21)-1)/20:.3f}")
print(f"  with the measured noise factor 2.44 in variance: sqrt(R_20 (0.5+sqrt(pi k 2.44/8)))="
      f"{np.sqrt(R[20]*(0.5+np.sqrt(np.pi*k*2.44/8))):.2f}")
