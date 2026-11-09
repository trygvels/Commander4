"""Measure the acceptance noise of the real E014 legacy chain from its own log.
Step 1 of every iteration: current amplitudes are the fresh draw of the amplitude group, proposed
amplitudes are a fresh draw, so dloglike = -(D + X' - X)/2 with X, X' independent fresh draws and
D = chi2_min(b') - chi2_min(b). Regress dloglike on a polynomial difference P(b') - P(b); the residual
variance is Var(X' - X)/4 = Var(chi2 at a fresh draw)/2, which is k = 190825 if X ~ chi2_k."""
import numpy as np
a = np.load("e014_steps.npz")["a"]
ch, it, st, dl, acc, cd, cs, pd_, ps = a.T
K = 190825
RS, RD = 0.00036, 0.000043          # exact ridge sds from E014 (beta_s, beta_d)
burn = 400
print("continuity of the state between iterations, and stationary widths from the log:")
for c in (1, 2):
    m = ch == c
    last = m & (st == 20); first = m & (st == 1)
    fd = np.where(acc[last] > 0, pd_[last], cd[last]); fs = np.where(acc[last] > 0, ps[last], cs[last])
    ok = np.allclose(fd[:-1], cd[first][1:]) and np.allclose(fs[:-1], cs[first][1:])
    keep = it[last] > burn
    d, s = fd[keep], fs[keep]
    def ess(x):
        x = x - x.mean(); n = len(x); v = x.var()
        tau = 1.0
        for lag in range(1, n // 3):
            r = np.dot(x[:-lag], x[lag:]) / ((n - lag) * v)
            if r < 0.05: break
            tau += 2 * r
        return n / tau
    print(f"  chain {c}: continuous={ok}; beta_d {d.mean():.5f} sd {d.std(ddof=1):.6f} = {d.std(ddof=1)/RD:.2f} ridge sd (ESS {ess(d):.0f});"
          f" beta_s {s.mean():.5f} sd {s.std(ddof=1):.5f} = {s.std(ddof=1)/RS:.2f} ridge sd (ESS {ess(s):.0f});"
          f" acc {acc[m & (it > burn)].mean():.4f}")
post = it > burn
print("acceptance by step index (both chains, post burn-in), with 1/(j+1) for reference:")
print("  " + " ".join(f"{j}:{acc[post & (st == j)].mean():.3f}" for j in (1, 2, 3, 4, 5, 10, 15, 20)))

def feats(d, s, order):
    cols = []
    for p in range(1, order + 1):
        for q in range(p + 1):
            cols.append(d ** (p - q) * s ** q)
    return np.array(cols).T

sel = post & (st == 1) & np.isfinite(dl)
md, ms = cd[sel].mean(), cs[sel].mean()
# scale to units of the legacy width so the design matrix is well conditioned
ud0, us0 = (cd[sel] - md) / (16.5 * RD), (cs[sel] - ms) / (17.1 * RS)
ud1, us1 = (pd_[sel] - md) / (16.5 * RD), (ps[sel] - ms) / (17.1 * RS)
y = dl[sel]
print(f"step-1 proposals used: {sel.sum()};  raw sd of dloglike {y.std():.1f};  toy premise sqrt(k) = {np.sqrt(K):.1f}")
print(f"{'poly order':>10} {'npar':>5} {'resid sd':>9} {'+-':>6} {'f = var/k':>10} {'+-':>6} {'kurtosis':>9} {'mean resid':>10}")
for order in (2, 3, 4, 5):
    Xd = feats(ud1, us1, order) - feats(ud0, us0, order)
    coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    r = y - Xd @ coef
    n, p = Xd.shape
    v = (r ** 2).sum() / (n - p)
    kurt = (r ** 4).mean() / (r ** 2).mean() ** 2
    sev = v * np.sqrt((kurt - 1) / n)
    print(f"{order:10d} {p:5d} {np.sqrt(v):9.1f} {0.5*sev/np.sqrt(v):6.1f} {v/K:10.3f} {sev/K:6.3f} {kurt:9.2f} {r.mean():10.1f}")
    if order == 3:
        r3, X3 = r, Xd
# robustness: by chain, by size of the proposed jump in beta_d, robust scale
c1 = ch[sel] == 1
jump = np.abs(pd_[sel] - cd[sel]) / RD
print("cubic fit residual, subsets:")
for name, mk in [("chain 1", c1), ("chain 2", ~c1), ("|jump_d| < 60 ridge sd", jump < 60), ("|jump_d| > 60 ridge sd", jump >= 60),
                 ("first half", it[sel] < 1200), ("second half", it[sel] >= 1200)]:
    rr = r3[mk]
    mad = 1.4826 * np.median(np.abs(rr - np.median(rr)))
    print(f"  {name:>24}: n={mk.sum():5d} sd {rr.std():7.1f}  f={rr.var()/K:.3f}  robust sd (MAD) {mad:7.1f}  f_MAD={mad**2/K:.3f}")
# Gaussian-noise check of the general law: T = 1/2 + sqrt(2 pi) sd(Z)/8 with sd(Z) = 2 resid sd
sdZ = 2 * r3.std()
print(f"sd(Z) = {sdZ:.0f} (toy premise {2*np.sqrt(K):.0f});  one-move small-step sqrt(T) from the measured noise = "
      f"{np.sqrt(0.5 + np.sqrt(2*np.pi)*sdZ/8):.2f} (toy premise {np.sqrt(0.5 + np.sqrt(2*np.pi)*2*np.sqrt(K)/8):.2f})")
