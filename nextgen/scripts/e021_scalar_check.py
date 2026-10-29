"""Scalar (single-index) run of E021 attempt 2: chain posterior of beta_s, beta_d, the RMSE of the
truth about them, and an exact collapsed-marginal scan (all pixels, 12x12 per pixel) in each scalar
index at the chain's amplitude hyperparameters: does the chain sit at the mode of its target, and
how far is that mode from the linearized Fisher-weighted sky mean (-2.9878, 1.5973)?"""
import glob, numpy as np, healpy as hp, h5py
D = "/Users/trygvels/data/e021_consistent_v2"; burn = 500
files = sorted(glob.glob(f"{D}/chains_scalar/chains_compsep/chain01_iter*.h5"))[burn:]
st = lambda fn, c, k: h5py.File(fn)[f"mcmc/sample_index_fields/state/{c}/{k}"][()]
B = np.array([[st(fn, "Synchrotron", "beta"), st(fn, "ThermalDust", "beta"), st(fn, "Synchrotron", "amp_prior/mu"), st(fn, "Synchrotron", "amp_prior/log_s"),
               st(fn, "ThermalDust", "amp_prior/mu"), st(fn, "ThermalDust", "amp_prior/log_s")] for fn in files]).squeeze()
m, sd = B.mean(0), B.std(0)
print(f"scalar chain: beta_s {m[0]:.5f} +- {sd[0]:.5f}, beta_d {m[1]:.5f} +- {sd[1]:.5f}; amp mu_s {m[2]:.5f} log_s_s {m[3]:.4f} mu_d {m[4]:.5f} log_s_d {m[5]:.4f}")
t = np.load(f"{D}/truth.npz"); ts, td = t["beta_s"], t["beta_d"]
for lab, bs, bd in (("chain", m[0], m[1]), ("Fisher-weighted mean", -2.9878, 1.5973), ("unweighted mean", ts.mean(), td.mean())):
    print(f"  RMSE of truth about {lab}: sync {np.sqrt(np.mean((ts-bs)**2)):.4f} dust {np.sqrt(np.mean((td-bd)**2)):.4f}")
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; nu = np.array(bands, float)
d = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands]).T
rms = np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands]).T
gs = np.maximum(hp.read_map(f"{D}/template_sync_30GHz.fits", field=0), 1e-3)
gd = np.maximum(hp.read_map(f"{D}/template_dust_545GHz.fits", field=0), 1e-3)
cmb_mix = np.array([h5py.File(files[0])[f"comps/cmb/mixing/Band{b:03d}"][()] for b in bands])
HK = 4.799243e-11
def dust_sed(b, T=20.9287, ref=545.0):
    x, xr = HK*nu*1e9/T, HK*ref*1e9/T
    return (nu/ref)**(b + 1.0) * (np.exp(xr) - 1.0)/(np.exp(x) - 1.0)
def marg(bs, bd, mus, ss, mud, sdd, cmb_sd=100.0):
    n = len(d)
    M = np.stack([np.broadcast_to(cmb_mix, (n, 12)), np.broadcast_to((nu/30.0)**bs, (n, 12)), np.broadcast_to(dust_sed(bd), (n, 12))], -1)
    mm = np.stack([np.zeros(n), mus*gs, mud*gd], -1)
    Tv = np.stack([np.full(n, cmb_sd**2), (ss*gs)**2, (sdd*gd)**2], -1)
    C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + np.einsum("ni,ij->nij", rms**2, np.eye(12))
    r = d - np.einsum("nij,nj->ni", M, mm)
    sol = np.linalg.solve(C, r[..., None])[..., 0]
    return -0.5*np.einsum("ni,ni->", r, sol) - 0.5*np.linalg.slogdet(C)[1].sum()
hy = (m[2], np.exp(m[3]), m[4], np.exp(m[5]))
for comp, c0, step in (("sync", m[0], 0.002), ("dust", m[1], 0.001)):
    grid = c0 + step*np.arange(-6, 7); v = []
    for g in grid:
        v.append(marg(g if comp == "sync" else m[0], g if comp == "dust" else m[1], *hy))
    v = np.array(v); i = np.argmax(v)
    # quadratic peak
    a, b, c = np.polyfit(grid[max(i-2,0):i+3], v[max(i-2,0):i+3], 2); pk = -b/(2*a); curv_sd = np.sqrt(-1/(2*a)) if a < 0 else np.nan
    print(f"exact-marginal scan in beta_{comp[0]} at the chain hyperparameters: grid argmax {grid[i]:.4f}, quadratic peak {pk:.5f} (curvature sd {curv_sd:.5f}); chain {c0:.5f} +- {sd[0 if comp=='sync' else 1]:.5f}; logp rel to max:", np.round(v - v.max(), 1))
# and at the realized amplitude targets instead of the chain hyperparameters
hy2 = (1.00309, 0.29888, 0.99986, 0.29897)
for comp, c0, step in (("sync", m[0], 0.004), ("dust", m[1], 0.002)):
    grid = c0 + step*np.arange(-6, 7); v = np.array([marg(g if comp == "sync" else m[0], g if comp == "dust" else m[1], *hy2) for g in grid]); i = np.argmax(v)
    print(f"  same at the realized amplitude targets: argmax beta_{comp[0]} {grid[i]:.4f}; logp rel to max:", np.round(v - v.max(), 1))
