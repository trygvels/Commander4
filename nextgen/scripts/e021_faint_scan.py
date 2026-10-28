"""Exact collapsed hier marginal (numpy, 12x12 per pixel) scanned in the faint-region dust/sync index
and in the amplitude hyperparameters, at the gmrf16 posterior-mean state and at the truth."""
import glob, sys, numpy as np, healpy as hp, h5py
D = "/Users/trygvels/data/e021_consistent"; tag = sys.argv[1] if len(sys.argv) > 1 else "gmrf16"; burn = 500
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; nu = np.array(bands, float)
d = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands]).T      # (npix, 12)
rms = np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands]).T
t = np.load(f"{D}/truth.npz")
gs = np.maximum(hp.read_map(f"{D}/template_sync_30GHz.fits", field=0), 1e-3)
gd = np.maximum(hp.read_map(f"{D}/template_dust_545GHz.fits", field=0), 1e-3)
files = sorted(glob.glob(f"{D}/chains_{tag}/chains_compsep/chain01_iter*.h5"))[burn:]
with h5py.File(files[0]) as f: cmb_mix = np.array([f[f"comps/cmb/mixing/Band{b:03d}"][()] for b in bands])
bs16 = np.mean([h5py.File(fn)["comps/sync/index_field/beta_map"][()] for fn in files], 0)
bd16 = np.mean([h5py.File(fn)["comps/dust/index_field/beta_map"][()] for fn in files], 0)
hy = np.mean([[h5py.File(fn)[f"mcmc/sample_index_fields/state/{c}/amp_prior/{k}"][()] for c in ("Synchrotron", "ThermalDust") for k in ("mu", "log_s")] for fn in files], 0)
print("chain amp hyper mean (mu_s, log_s_s, mu_d, log_s_d):", np.round(hy, 4), " realized: 0.9990", round(np.log(0.2986), 4), "0.9976", round(np.log(0.2988), 4))
HK = 4.799243e-11
def dust_sed(b, T=20.9287, ref=545.0):
    x, xr = HK*nu*1e9/T, HK*ref*1e9/T
    return (nu/ref)**(b[:, None] + 1.0) * (np.exp(xr) - 1.0)/(np.exp(x) - 1.0)
def sync_sed(b): return (nu/30.0)**b[:, None]
def marg(bs, bd, mus, ss, mud, sd, sel, cmb_sd=100.0):
    """Sum over pixels sel of log N(d; M m, N + M T M^T) split into (chi2 part, logdet part)."""
    M = np.stack([np.broadcast_to(cmb_mix, (sel.sum(), 12)), sync_sed(bs[sel]), dust_sed(bd[sel])], -1)  # (n,12,3)
    m = np.stack([np.zeros(sel.sum()), mus*gs[sel], mud*gd[sel]], -1)
    Tv = np.stack([np.full(sel.sum(), cmb_sd**2), (ss*gs[sel])**2, (sd*gd[sel])**2], -1)
    C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + np.einsum("ni,ij->nij", rms[sel]**2, np.eye(12))
    r = d[sel] - np.einsum("nij,nj->ni", M, m)
    sol = np.linalg.solve(C, r[..., None])[..., 0]
    chi = -0.5*np.einsum("ni,ni->", r, sol); ld = -0.5*np.linalg.slogdet(C)[1].sum()
    return chi, ld
bs64, bd64 = hp.ud_grade(bs16, 64), hp.ud_grade(bd16, 64)
ts64, td64 = hp.ud_grade(hp.ud_grade(t["beta_s"], 16), 64), hp.ud_grade(hp.ud_grade(t["beta_d"], 16), 64)
fsd16 = t["fisher_sd_dust_16_profiled"]; fss16 = t["fisher_sd_sync_16_profiled"]
qd = np.argsort(fsd16)[-768:]; qs = np.argsort(fss16)[-768:]     # lowest S/N quartile regions
seld = np.isin(hp.ud_grade(np.arange(3072.0), 64).astype(int), qd); sels = np.isin(hp.ud_grade(np.arange(3072.0), 64).astype(int), qs)
print(f"faint dust quartile: {seld.sum()} pixels; truth beta_d mean {td64[seld].mean():.4f}; chain mean {bd64[seld].mean():.4f}")
deltas = np.arange(-0.15, 0.301, 0.025)
def scan(label, comp, state, sel, **kw):
    bs, bd, mus, ss, mud, sd = state
    out = []
    for dl in deltas:
        b2 = (bs + dl*sel) if comp == "s" else bs; d2 = (bd + dl*sel) if comp == "d" else bd
        out.append(marg(b2, d2, mus, ss, mud, sd, sel, **kw))
    out = np.array(out); tot = out.sum(1); i0 = np.argmin(abs(deltas)); 
    print(f"{label}: delta  chi2part  logdetpart  total (rel to delta 0); argmax total at delta {deltas[np.argmax(tot)]:+.3f}")
    for dl, (c, l), tt in zip(deltas, out - out[i0], tot - tot[i0]): print(f"   {dl:+.3f} {c:12.1f} {l:12.1f} {tt:12.1f}")
chain = (bs64, bd64, hy[0], np.exp(hy[1]), hy[2], np.exp(hy[3]))
truth = (ts64, td64, 0.9990, 0.2986, 0.9976, 0.2988)
scan("DUST faint quartile shift, chain state", "d", chain, seld)
scan("DUST faint quartile shift, truth state", "d", truth, seld)
scan("DUST faint quartile shift, truth state, s_d = 0.03 (template-pinned)", "d", (ts64, td64, 0.999, 0.2986, 0.9976, 0.03), seld)
scan("DUST faint quartile shift, truth state, s_d = 3 (near flat)", "d", (ts64, td64, 0.999, 0.2986, 0.9976, 3.0), seld)
scan("SYNC faint quartile shift, truth state", "s", truth, sels)
# amplitude hyperparameter mu_s at fixed everything else, all pixels
allp = np.ones(len(d), bool)
for lab, st in (("chain betas", chain), ("truth betas", truth)):
    grid = np.arange(0.95, 1.011, 0.005); v = []
    for mu in grid: v.append(sum(marg(st[0], st[1], mu, st[3], st[4], st[5], allp)))
    v = np.array(v); print(f"mu_s scan at {lab}: argmax {grid[np.argmax(v)]:.3f}; logp rel to max:", np.round(v - v.max(), 1))
    grid = np.arange(0.95, 1.011, 0.005); v = []
    for mu in grid: v.append(sum(marg(st[0], st[1], st[2], st[3], mu, st[5], allp)))
    v = np.array(v); print(f"mu_d scan at {lab}: argmax {grid[np.argmax(v)]:.3f}; logp rel to max:", np.round(v - v.max(), 1))
