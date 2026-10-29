"""Is the iid16 dust/sync bias a property of the iid target? (1) errors by S/N quartile for
gmrf16 and iid16; (2) the linear-Gaussian iid prediction: mu_hat = precision-weighted truth mean
with w_r = 1/(sigma^2 + fsd_r^2), regional shrinkage m_r = (t_r/fsd_r^2 + mu/sigma^2)/(1/fsd_r^2
+ 1/sigma^2), its RMSE and quartile errors; (3) exact-marginal + iid-prior block-shift scan of the
two middle dust quartiles and the lowest sync quartile at the iid16 chain state."""
import glob, sys, numpy as np, healpy as hp, h5py
D = "/Users/trygvels/data/e021_consistent_v2"; burn = 500
t = np.load(f"{D}/truth.npz")
td16 = hp.ud_grade(t["beta_d"], 16); ts16 = hp.ud_grade(t["beta_s"], 16)
fsd16 = t["fisher_sd_dust_16_profiled"]; fss16 = t["fisher_sd_sync_16_profiled"]
def load(tag):
    files = sorted(glob.glob(f"{D}/chains_{tag}/chains_compsep/chain01_iter*.h5"))[burn:]
    bs = np.array([h5py.File(fn)["comps/sync/index_field/beta_map"][()] for fn in files])
    bd = np.array([h5py.File(fn)["comps/dust/index_field/beta_map"][()] for fn in files])
    st = lambda fn, c, k: h5py.File(fn)[f"mcmc/sample_index_fields/state/{c}/index_field/{k}"][()]
    hy = np.mean([[st(fn, "Synchrotron", "mu"), np.exp(st(fn, "Synchrotron", "log_sigma")),
                   st(fn, "ThermalDust", "mu"), np.exp(st(fn, "ThermalDust", "log_sigma"))] for fn in files], 0)
    ap = np.mean([[h5py.File(fn)[f"mcmc/sample_index_fields/state/{c}/amp_prior/{k}"][()] for c in ("Synchrotron", "ThermalDust") for k in ("mu", "log_s")] for fn in files], 0)
    return bs.mean(0), bd.mean(0), hy, ap
def quart(label, m, tr, fsd):
    q = np.argsort(fsd); err = m - tr; out = []
    for i in range(4):
        idx = q[i*768:(i+1)*768]
        out.append(f"Q{i+1}(S/N high->low) truth {tr[idx].mean():.4f} err {err[idx].mean():+.4f} rms {np.sqrt((err[idx]**2).mean()):.4f}")
    print(f"{label}: overall err mean {err.mean():+.4f} rms {np.sqrt((err**2).mean()):.4f}\n   " + "\n   ".join(out))
res = {}
for tag in ("gmrf16", "iid16"):
    bs, bd, hy, ap = load(tag); res[tag] = (bs, bd, hy, ap)
    print(f"== {tag}: sync mu {hy[0]:.4f} sigma {hy[1]:.4f}; dust mu {hy[2]:.4f} sigma {hy[3]:.4f}")
    quart(f"{tag} dust", bd, td16, fsd16); quart(f"{tag} sync", bs, ts16, fss16)
print("truth means: dust", round(td16.mean(), 4), "sync", round(ts16.mean(), 4))
# (2) linear-Gaussian iid prediction, iterating mu to its fixed point (empirical Bayes at fixed sigma)
for name, tr, fsd, sig in (("dust", td16, fsd16, res["iid16"][2][3]), ("sync", ts16, fss16, res["iid16"][2][1])):
    w = 1.0/(sig**2 + fsd**2); mu_hat = np.sum(w*tr)/np.sum(w)
    m_pred = (tr/fsd**2 + mu_hat/sig**2)/(1/fsd**2 + 1/sig**2)
    print(f"iid prediction {name}: sigma {sig:.4f}; precision-weighted truth mean {mu_hat:.4f} (unweighted {tr.mean():.4f}); "
          f"predicted shrinkage RMSE {np.sqrt(((m_pred-tr)**2).mean()):.4f}")
    quart(f"   predicted {name}", m_pred, tr, fsd)
    # same with the posterior-mean field of iid16 in place of the truth: does the chain's mu match its own weighted field mean?
    m_chain = res["iid16"][1] if name == "dust" else res["iid16"][0]
    print(f"   chain field: weighted mean {np.sum(w*m_chain)/np.sum(w):.4f}, unweighted {m_chain.mean():.4f}, chain mu {res['iid16'][2][2 if name=='dust' else 0]:.4f}")
# (3) exact marginal + iid prior block shift at the iid16 chain state
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; nu = np.array(bands, float)
d = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands]).T
rms = np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands]).T
gs = np.maximum(hp.read_map(f"{D}/template_sync_30GHz.fits", field=0), 1e-3)
gd = np.maximum(hp.read_map(f"{D}/template_dust_545GHz.fits", field=0), 1e-3)
fn0 = sorted(glob.glob(f"{D}/chains_iid16/chains_compsep/chain01_iter*.h5"))[burn]
cmb_mix = np.array([h5py.File(fn0)[f"comps/cmb/mixing/Band{b:03d}"][()] for b in bands])
HK = 4.799243e-11
def dust_sed(b, T=20.9287, ref=545.0):
    x, xr = HK*nu*1e9/T, HK*ref*1e9/T
    return (nu/ref)**(b[:, None] + 1.0) * (np.exp(xr) - 1.0)/(np.exp(x) - 1.0)
def sync_sed(b): return (nu/30.0)**b[:, None]
def marg(bs, bd, mus, ss, mud, sd, sel, cmb_sd=100.0):
    M = np.stack([np.broadcast_to(cmb_mix, (sel.sum(), 12)), sync_sed(bs[sel]), dust_sed(bd[sel])], -1)
    m = np.stack([np.zeros(sel.sum()), mus*gs[sel], mud*gd[sel]], -1)
    Tv = np.stack([np.full(sel.sum(), cmb_sd**2), (ss*gs[sel])**2, (sd*gd[sel])**2], -1)
    C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + np.einsum("ni,ij->nij", rms[sel]**2, np.eye(12))
    r = d[sel] - np.einsum("nij,nj->ni", M, m)
    sol = np.linalg.solve(C, r[..., None])[..., 0]
    return -0.5*np.einsum("ni,ni->", r, sol) - 0.5*np.linalg.slogdet(C)[1].sum()
bs16, bd16, hy, ap = res["iid16"]
bs64, bd64 = hp.ud_grade(bs16, 64), hp.ud_grade(bd16, 64)
reg64 = hp.ud_grade(np.arange(3072.0), 64).astype(int)
deltas = np.arange(-0.06, 0.061, 0.01)
def scan(label, comp, regs):
    sel = np.isin(reg64, regs); tr = (td16 if comp == "d" else ts16)[regs].mean()
    mu, sig = (hy[2], hy[3]) if comp == "d" else (hy[0], hy[1]); b16 = bd16 if comp == "d" else bs16
    like, prior = [], []
    for dl in deltas:
        b2s = bs64 + dl*sel if comp == "s" else bs64; b2d = bd64 + dl*sel if comp == "d" else bd64
        like.append(marg(b2s, b2d, ap[0], np.exp(ap[1]), ap[2], np.exp(ap[3]), sel))
        prior.append(-0.5*np.sum((b16[regs] + dl - mu)**2)/sig**2)
    like, prior = np.array(like), np.array(prior); tot = like + prior; i0 = np.argmin(abs(deltas))
    print(f"{label}: regions {len(regs)}, truth mean {tr:.4f}, chain mean {b16[regs].mean():.4f}; "
          f"argmax like {deltas[np.argmax(like)]:+.3f}, argmax like+prior {deltas[np.argmax(tot)]:+.3f}")
    for dl, l, p, tt in zip(deltas, like - like[i0], prior - prior[i0], tot - tot[i0]): print(f"   {dl:+.3f} like {l:10.1f} prior {p:10.1f} total {tt:10.1f}")
qd = np.argsort(fsd16); qs = np.argsort(fss16)
scan("DUST Q2 (2nd S/N quartile) block shift, iid16 chain state", "d", qd[768:1536])
scan("DUST Q3 block shift, iid16 chain state", "d", qd[1536:2304])
scan("SYNC Q4 (lowest S/N) block shift, iid16 chain state", "s", qs[2304:])
