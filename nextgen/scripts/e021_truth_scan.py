"""Exact collapsed hier marginal (numpy, 12x12 per pixel) at the truth state of a consistent sky:
scan of a common shift of the faint-quartile dust and sync index and of the amplitude hyperparameters."""
import sys, numpy as np, healpy as hp
D = sys.argv[1]
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; nu = np.array(bands, float)
d = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands]).T
rms = np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands]).T
t = np.load(f"{D}/truth.npz")
gs = np.maximum(hp.read_map(f"{D}/template_sync_30GHz.fits", field=0), 1e-3)
gd = np.maximum(hp.read_map(f"{D}/template_dust_545GHz.fits", field=0), 1e-3)
HK = 4.799243e-11; T = float(t["T_fixed"])
x = HK*nu*1e9/2.7255; cmb_mix = x**2*np.exp(x)/np.expm1(x)**2
def dust_sed(b, ref=545.0):
    xx, xr = HK*nu*1e9/T, HK*ref*1e9/T
    return (nu/ref)**(b[:, None] + 1.0) * np.expm1(xr)/np.expm1(xx)
def sync_sed(b): return (nu/30.0)**b[:, None]
def marg(bs, bd, mus, ss, mud, sd, sel, cmb_sd=100.0):
    n = sel.sum()
    M = np.stack([np.broadcast_to(cmb_mix, (n, 12)), sync_sed(bs[sel]), dust_sed(bd[sel])], -1)
    m = np.stack([np.zeros(n), mus*gs[sel], mud*gd[sel]], -1)
    Tv = np.stack([np.full(n, cmb_sd**2), (ss*gs[sel])**2, (sd*gd[sel])**2], -1)
    C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + np.einsum("ni,ij->nij", rms[sel]**2, np.eye(12))
    r = d[sel] - np.einsum("nij,nj->ni", M, m)
    sol = np.linalg.solve(C, r[..., None])[..., 0]
    return -0.5*np.einsum("ni,ni->", r, sol), -0.5*np.linalg.slogdet(C)[1].sum()
ts64, td64 = t["beta_s"], t["beta_d"]
rs, rd = t["a_sync"]/gs, t["a_dust"]/gd
truth = (ts64, td64, rs.mean(), rs.std(), rd.mean(), rd.std())
print("realized a/g: sync", rs.mean().round(5), rs.std().round(5), "dust", rd.mean().round(5), rd.std().round(5))
reg = hp.ud_grade(np.arange(3072.0), 64).astype(int)
deltas = np.arange(-0.10, 0.251, 0.025)
for comp, key in (("d", "fisher_sd_dust_16_profiled"), ("s", "fisher_sd_sync_16_profiled")):
    q = np.argsort(t[key])[-768:]; sel = np.isin(reg, q)
    out = []
    for dl in deltas:
        out.append(marg(ts64 + (dl*sel if comp == "s" else 0), td64 + (dl*sel if comp == "d" else 0), *truth[2:], sel))
    out = np.array(out); tot = out.sum(1); i0 = np.argmin(abs(deltas)); tot -= tot[i0]
    j = min(max(np.argmax(tot), 1), len(deltas)-2); c = np.polyfit(deltas[j-1:j+2], tot[j-1:j+2], 2)
    curv = -tot[j-1:j+2] @ np.array([1, -2, 1]) / 0.025**2
    print(f"{'dust' if comp=='d' else 'sync'} faint quartile ({sel.sum()} px) common shift at the truth: argmax {deltas[np.argmax(tot)]:+.3f}, quadratic peak {-c[1]/(2*c[0]):+.4f}, sd of the common shift {1/np.sqrt(curv):.4f}; total at +0.125 {tot[np.argmin(abs(deltas-0.125))]:.1f} (chi2 part {out[np.argmin(abs(deltas-0.125)),0]-out[i0,0]:.1f}, logdet part {out[np.argmin(abs(deltas-0.125)),1]-out[i0,1]:.1f})")
allp = np.ones(len(d), bool); grid = np.arange(0.97, 1.031, 0.005)
v = np.array([sum(marg(ts64, td64, mu, truth[3], truth[4], truth[5], allp)) for mu in grid]); j = np.argmax(v); c = np.polyfit(grid[j-1:j+2], v[j-1:j+2], 2)
print(f"mu_s scan at the truth: peak {-c[1]/(2*c[0]):.4f} (realized {rs.mean():.4f}); curvature sd {np.sqrt(-1/(2*c[0])):.4f}")
v = np.array([sum(marg(ts64, td64, truth[2], truth[3], mu, truth[5], allp)) for mu in grid]); j = np.argmax(v); c = np.polyfit(grid[j-1:j+2], v[j-1:j+2], 2)
print(f"mu_d scan at the truth: peak {-c[1]/(2*c[0]):.4f} (realized {rd.mean():.4f}); curvature sd {np.sqrt(-1/(2*c[0])):.4f}")
grid = np.arange(-1.30, -1.10, 0.01)
v = np.array([sum(marg(ts64, td64, truth[2], np.exp(ls), truth[4], truth[5], allp)) for ls in grid]); j = np.argmax(v); c = np.polyfit(grid[j-1:j+2], v[j-1:j+2], 2)
print(f"log_s_s scan at the truth: peak {-c[1]/(2*c[0]):.4f} (realized {np.log(rs.std()):.4f})")
v = np.array([sum(marg(ts64, td64, truth[2], truth[3], truth[4], np.exp(ls), allp)) for ls in grid]); j = np.argmax(v); c = np.polyfit(grid[j-1:j+2], v[j-1:j+2], 2)
print(f"log_s_d scan at the truth: peak {-c[1]/(2*c[0]):.4f} (realized {np.log(rd.std()):.4f})")
