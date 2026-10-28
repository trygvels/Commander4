"""E022: diagnostics of the released BeyondPlanck chain's per-region polarized beta_s
(read by bp_chain_beta_s.py). For a two-block Gibbs scheme var(delta beta) / 2 = E[var_cond]
exactly, r^2 = var(beta) / E[var_cond], lag-1 = 1 - 1/r^2; drift is tested by half-chain
means in units of the chain's own Monte Carlo error. Usage: bp_chain_analyze.py file.npz [...]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__)); from bpxv_gibbs_coupling import tau_ips
NAMES = ("region 0 (prior-only)", "region 1 (Spur)", "region 2 (prior-only)", "region 3 (Plane)")
BURN = 10
for fn in sys.argv[1:]:
    z = np.load(fn); P = z["beta_pixreg_val"][:, 1, :].astype(float); n_all = P.shape[0]
    print(f"== {os.path.basename(fn)}: {n_all} samples, burn {BURN}")
    for r in range(4):
        x = P[BURN:, r]; n = x.size; xc = x - x.mean(); lag1 = np.sum(xc[1:] * xc[:-1]) / np.sum(xc**2)
        tau, err = tau_ips(x); var = x.var(ddof=1); dx = np.diff(x); e_vc = 0.5 * np.var(dx)
        r2 = var / e_vc; tau_ar1 = (1 + lag1) / (1 - lag1) if lag1 < 1 else np.inf
        h = n // 2; m1, m2 = x[:h].mean(), x[h:].mean(); mc_half = np.sqrt(var * tau / h)  # MC error of one half mean
        acf = [np.sum(xc[k:] * xc[:-k]) / np.sum(xc**2) for k in (5, 10, 20, 50, 100, 200)]
        t = np.arange(n) - (n - 1) / 2; slope = np.sum(t * xc) / np.sum(t**2)  # per sample
        slope_err = np.sqrt(var * tau / np.sum(t**2))  # white-noise formula with n_eff = n / tau
        sd_detr = np.std(xc - slope * t, ddof=2)
        print(f"    trend {1000*slope:+.4f} per 1000 samples (+- {1000*slope_err:.4f}, {slope/slope_err:+.1f} sigma); detrended sd {sd_detr:.4f}")
        print(f"  {NAMES[r]:22s} mean {x.mean():.4f} sd {np.sqrt(var):.4f} | sd_cond from increments {np.sqrt(e_vc):.4f} -> r {np.sqrt(r2):.2f}, lag-1 identity {1-1/r2:.3f} vs measured {lag1:.3f} | tau IPS {tau:.1f}+-{err:.1f}, AR(1) from lag-1 {tau_ar1:.1f} | ESS {n/tau:.0f} | halves {m1:.4f} / {m2:.4f}, diff {(m2-m1):+.4f} = {(m2-m1)/(mc_half*np.sqrt(2)):+.1f} MC err | acf at lags 5,10,20,50,100,200: {' '.join(f'{a:.2f}' for a in acf)}")
