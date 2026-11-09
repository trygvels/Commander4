import numpy as np, glob, sys
fs = sorted(glob.glob(sys.argv[1] + "/rt_single_*.npz"))
Xp = np.concatenate([np.load(f)["Xp"] for f in fs]); Xr = np.concatenate([np.load(f)["Xr"] for f in fs]); Xq = np.concatenate([np.load(f)["Xq"] for f in fs])
z = np.load(fs[0]); B, k = float(z["B"]), int(z["k"]); n = len(Xp)
def rep(name, X):
    v = X.var(ddof=1); kurt = ((X - X.mean()) ** 4).mean() / X.var() ** 2
    print(f"{name:>44}: mean {X.mean()/k:6.3f} k; sd {np.sqrt(v):7.1f}; var/(2k) = {v/(2*k):.3f} +- {v/(2*k)*np.sqrt((kurt-1)/n):.3f}; kurt {kurt:.2f}")
print(f"pooled {n} draws from {len(fs)} runs; B = {B:.0f} = {B/k:.3f} k")
rep("pixel-space draw", Xp); rep("after round trip (as global_loglike)", Xr); rep("quadratic part (P delta)^T A (P delta)", Xq)
cross = Xr - Xq - B
rep("cross term -2 (P delta)^T A (Q a_hat)", cross)
print(f"mean of cross {cross.mean():.1f} +- {cross.std()/np.sqrt(n):.1f}; corr(quad, cross) = {np.corrcoef(Xq, cross)[0,1]:.3f}; var(cross)/(4B) = {cross.var(ddof=1)/(4*B):.3f}")
