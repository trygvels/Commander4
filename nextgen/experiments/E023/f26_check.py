"""Independent check of F26: expected log-likelihood of the collapsed per-pixel marginal per
quartile versus a block index shift, for several CMB prior sds and CMB realizations."""
import numpy as np, healpy as hp, sys
D = "/Users/trygvels/data/e023_pysm_n1"
R = "/Users/trygvels/code/Commander4/nextgen/experiments/E023/scan_n1_refs.npz"
H_K = 4.799243073366221e-11; T_CMB = 2.7255
F = np.array([23., 30., 33., 41., 44., 61., 70., 94., 100., 143., 217., 353.])
def cmb_mix(nu):
    x = H_K*nu*1e9/T_CMB; return x**2*np.exp(x)/np.expm1(x)**2
def sync(bs): return (F[None]/30.)**bs[:, None]
def dust(bd, T):
    x = H_K*F[None]*1e9/T[:, None]; xr = H_K*545e9/T[:, None]
    return (F[None]/545.)**(bd[:, None]+1.)*np.expm1(xr)/np.expm1(x)
t = np.load(f"{D}/truth.npz"); refs = np.load(R)
sky = np.asarray(t["sky"], float); cmb = np.asarray(t["cmb"], float)
rms = np.array([hp.read_map(f"{D}/band{int(f):03d}_rms.fits", field=0, verbose=False) for f in F]).T
a_s, a_d, T = t["a_sync"], t["a_dust"], t["T_d"]
reg = hp.ud_grade(np.arange(3072, dtype=float), 64).astype(int)
bs0, bd0 = refs["bs_trueT"][reg], refs["bd_trueT"][reg]
qs, qd = refs["qs16"][reg], refs["qd16"][reg]
cm = cmb_mix(F)
THETA = {100.: (0.9944, 0.0316, 1.0021, 0.0227), 1e4: (0.9910, 0.0484, 1.0021, 0.0206)}

def ell(sel, s, rmsm, tau, bs, bd, th):
    """Expected log-lik parts over pixels sel: chi2, logdet, trace, plus the CMB-only quadratic
    piece c^2 m^T C^-1 m (cq) and the linear cross piece (lin)."""
    mus, ss, mud, sd = th
    n = sel.sum()
    M = np.stack([np.broadcast_to(cm, (n, 12)), sync(bs[sel]), dust(bd[sel], T[sel])], -1)
    m = np.stack([np.zeros(n), mus*a_s[sel], mud*a_d[sel]], -1)
    Tv = np.stack([np.full(n, tau**2), (ss*a_s[sel])**2, (sd*a_d[sel])**2], -1)
    N = rmsm[sel]**2
    C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + N[:, None, :]*np.eye(12)[None]
    r = s[sel] - np.einsum("nij,nj->ni", M, m)
    L = np.linalg.cholesky(C)
    ld = -np.sum(np.log(np.einsum("nii->ni", L)))
    Ci = np.linalg.inv(C)
    chi = -0.5*np.einsum("ni,nij,nj->", r, Ci, r)
    tr = -0.5*np.einsum("nii,ni->", Ci, N)
    return chi, ld, tr

def scan(label, sel, s, rmsm, tau, comp, deltas, th):
    out = []
    for dl in deltas:
        bs = bs0 + (dl*sel if comp == "s" else 0); bd = bd0 + (dl*sel if comp == "d" else 0)
        out.append(ell(sel, s, rmsm, tau, bs, bd, th))
    out = np.array(out); out -= out[np.argmin(np.abs(deltas))]
    tot = out.sum(1); i = np.argmax(tot)
    # parabolic argmax on the three points around i (or edge)
    if 0 < i < len(deltas)-1:
        y0, y1, y2 = tot[i-1], tot[i], tot[i+1]; h = deltas[1]-deltas[0]
        pk = deltas[i] + 0.5*h*(y0-y2)/(y0-2*y1+y2); sdv = np.sqrt(-h**2/(y0-2*y1+y2))
    else:
        pk, sdv = deltas[i], np.nan
    print(f"{label:34s} tau {tau:7g}: argmax {pk:+.3f} (sd {sdv:.4f}) gain {tot[i]:7.2f}; "
          f"chi2 at -0.1/+0.1 {out[0,0]:8.2f}/{out[-1,0]:8.2f}; logdet {out[0,1]:+.2f}/{out[-1,1]:+.2f}; "
          f"trace {out[0,2]:+.2f}/{out[-1,2]:+.2f}", flush=True)
    return tot

deltas = np.arange(-0.10, 0.1001, 0.02)
variants = {"sky": sky, "cmb flipped": sky - 2*cmb[:, None]*cm[None],
            "cmb rolled": sky - cmb[:, None]*cm[None] + np.roll(cmb, 24576)[:, None]*cm[None],
            "no cmb": sky - cmb[:, None]*cm[None]}
# check that sky = cmb*cm + fg: residual of sky - cmb*cm against a_s*sync + a_d*dust
fg = a_s[:, None]*sync(t["beta_s"]) + a_d[:, None]*dust(t["beta_d"], T)
print("sky - cmb*cm - fg rms per band:", np.round(np.sqrt(np.mean((sky - cmb[:, None]*cm[None] - fg)**2, 0)), 3))
scans = [("d_Q2", qd == 1, "d"), ("d_Q3", qd == 2, "d"), ("d_Q4", qd == 3, "d"), ("s_Q4", qs == 3, "s"), ("d_Q1", qd == 0, "d")]
res = {}
for label, sel, comp in scans:
    for vname, s in variants.items():
        for tau in (100., 84., 300., 1e4):
            if vname != "sky" and tau not in (100.,):
                continue
            th = THETA[1e4] if tau == 1e4 else THETA[100.]
            res[(label, vname, tau)] = scan(f"{label} {vname}", sel, s, rms, tau, comp, deltas, th)
    # symmetric / antisymmetric split under the CMB sign flip at tau 100
    a, b = res[(label, "sky", 100.)], res[(label, "cmb flipped", 100.)]
    print(f"   {label} tau 100: symmetric part (quadratic in cmb) at -0.1/+0.1: {(a[0]+b[0])/2:+.2f}/{(a[-1]+b[-1])/2:+.2f}; "
          f"antisymmetric (linear in cmb) at -0.1/+0.1: {(a[0]-b[0])/2:+.2f}/{(a[-1]-b[-1])/2:+.2f}", flush=True)
# noise scaling: 0.03 x rms, dust Q2, tau 100, sky and flipped
print("--- noise scale 0.03")
for vname in ("sky", "cmb flipped"):
    scan(f"d_Q2 {vname} noise 0.03", qd == 1, variants[vname], 0.03*rms, 100., "d", np.arange(-0.004, 0.0041, 0.001), THETA[100.])
# tau 1e4 with the tau-100 hyperparameters (sensitivity to the fixed theta)
scan("d_Q2 sky theta100", qd == 1, sky, rms, 1e4, "d", deltas, THETA[100.])
