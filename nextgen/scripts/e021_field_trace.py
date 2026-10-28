"""Mixing of the dust index field in an E021 chain: per-region traces, acf, offsets by S/N."""
import glob, sys, numpy as np, healpy as hp, h5py
D = "/Users/trygvels/data/e021_consistent"; tag = sys.argv[1]; burn = int(sys.argv[2])
files = sorted(glob.glob(f"{D}/chains_{tag}/chains_compsep/chain01_iter*.h5"))
bs, bd, hyp = [], [], []
for fn in files:
    with h5py.File(fn) as f:
        bs.append(f["comps/sync/index_field/beta_map"][()]); bd.append(f["comps/dust/index_field/beta_map"][()])
        st = f["mcmc/sample_index_fields/state"]
        hyp.append([st["ThermalDust/index_field/mu"][()], st["ThermalDust/index_field/log_sigma"][()], st["ThermalDust/index_field/u"][()],
                    st["Synchrotron/index_field/mu"][()], st["Synchrotron/index_field/log_sigma"][()], st["Synchrotron/index_field/u"][()]])
bs, bd, hyp = np.array(bs), np.array(bd), np.array(hyp)
t = np.load(f"{D}/truth.npz"); print(list(t.keys()))
td16 = hp.ud_grade(t["beta_d"], 16); ts16 = hp.ud_grade(t["beta_s"], 16)
fsd = t["fisher_sd_dust_16_profiled"] if "fisher_sd_dust_16_profiled" in t else None
def acf(x, k): x = x - x.mean(); return float(np.sum(x[:-k]*x[k:])/np.sum(x*x))
def ess(x):
    x = x - x.mean(); n = len(x); c = np.correlate(x, x, "full")[n-1:]/np.sum(x*x)
    s = 0.0
    for k in range(0, n-1, 2):
        p = c[k] + c[k+1]
        if p <= 0: break
        s += p
    return n/(2*s-1) if s > 0.5 else n
post = slice(burn, None)
for name, b, tr in (("dust", bd, td16), ("sync", bs, ts16)):
    m = b[post].mean(0); sd = b[post].std(0); err = m - tr
    e = np.array([ess(b[post, r]) for r in range(b.shape[1])])
    q = np.argsort(fsd if fsd is not None else sd)  # low to high Fisher sd = high to low S/N
    print(f"== {name}: n_post {b[post].shape[0]}; ESS per region: min {e.min():.1f} 5% {np.percentile(e,5):.1f} median {np.median(e):.1f}")
    print(f"   posterior-mean error (m - truth): mean {err.mean():+.4f} rms {np.sqrt((err**2).mean()):.4f}; region sd median {np.median(sd):.4f}")
    for lab, idx in (("highest S/N quartile", q[:768]), ("lowest S/N quartile", q[-768:])):
        print(f"   {lab}: truth mean {tr[idx].mean():.4f} sd {tr[idx].std():.4f}; post mean {m[idx].mean():.4f}; err mean {err[idx].mean():+.4f} rms {np.sqrt((err[idx]**2).mean()):.4f}; sd median {np.median(sd[idx]):.4f}; ESS median {np.median(e[idx]):.1f}")
    worst = np.argsort(e)[:3]
    for r in worst:
        x = b[:, r]
        print(f"   region {r}: truth {tr[r]:.3f} start {x[0]:.3f} post mean {m[r]:.3f} sd {sd[r]:.3f} ESS {e[r]:.1f}; acf lag 1/10/50/100 {acf(x[post],1):.2f}/{acf(x[post],10):.2f}/{acf(x[post],50):.2f}/{acf(x[post],100):.2f}; trace every 100: {np.round(x[::100],3)}")
    # field-mean trajectory
    fm = b.mean(1); print(f"   field mean every 100: {np.round(fm[::100],4)}; truth mean {tr.mean():.4f}")
H = ["mu_d","logsig_d","u_d","mu_s","logsig_s","u_s"]
for i, h in enumerate(H):
    x = hyp[:, i]; print(f"{h:9s} post mean {x[post].mean():+.4f} sd {x[post].std():.4f} ESS {ess(x[post]):.1f} acf1 {acf(x[post],1):.2f} acf10 {acf(x[post],10):.2f}; every 100: {np.round(x[::100],3)}")
