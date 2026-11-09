"""Direct test of the open item of E035 on the E014 inputs (read-only, no MPI, no MH code).
At the exact ridge mean of (beta_s, beta_d): draw per-pixel amplitudes a = a_hat + delta,
delta_p ~ N(0, A_p^-1) on the full sky (as solve_compsep_perpix does), then evaluate the chi2 inside
the index mask (a) at the pixel-space draw and (b) after the storage round trip of each component map
(map2alm_healpix niter=3, lmax=191, then alm_to_map), as global_loglike sees it. X = chi2 - chi2_min.
Exact draws give X ~ chi2_k, k = 5 x 38165 = 190825: mean k, variance 2k."""
import os, sys, time
import numpy as np
import healpy as hp
from pixell import curvedsky
from commander4.parameters.parse import load_params
from commander4.sky.comp_list import CompList
from commander4.math_utils.sht import alm_to_map

pf, cache, ndraw, seed, prec = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
t0 = time.time()
params, _, _ = load_params(pf)
names = [n for n in params.compsep.bands if params.compsep.bands[n].enabled]
d, ivar, freqs = [], [], []
for n in names:
    z = np.load(os.path.join(cache, f"band_{n}.npz"))
    d.append(z["d"]); ivar.append(z["ivar"]); freqs.append(float(z["freq"]))
d, ivar, freqs = np.array(d), np.array(ivar), np.array(freqs)
keep = np.load(os.path.join(cache, "mask.npy"))
comps = list(CompList.init_from_params(params.components, params).split_for_eval_pol("I"))
for c in comps:
    if type(c).__name__ == "Synchrotron": c.beta = -3.34244
    if type(c).__name__ == "ThermalDust": c.beta = 1.36539
M = np.stack([c.get_sed(freqs) for c in comps], axis=1).astype(np.float64)
nb, npix = d.shape; nc = M.shape[1]; nside = 64; lmax = 191
k = nc * int(keep.sum())
print(f"{nb} bands, {nc} comps {[type(c).__name__ for c in comps]}, kept {int(keep.sum())}, k = {k}, precision {prec}")
A = np.einsum("bi,bp,bj->pij", M, ivar, M)
bvec = np.einsum("bi,bp->pi", M, ivar * d)
ahat = np.linalg.solve(A, bvec[..., None])[..., 0]          # (npix, nc)
L = np.linalg.cholesky(A)
Lt = np.ascontiguousarray(L.transpose(0, 2, 1))
cond = np.linalg.cond(A[keep])
print(f"condition number of A_p inside the mask: median {np.median(cond):.2e}, 99th pct {np.percentile(cond, 99):.2e}")
rdt, cdt = (np.float32, np.complex64) if prec == "single" else (np.float64, np.complex128)

def chi2_mask(a):              # a: (npix, nc) amplitudes in pixel space
    r = d - (M @ a.T)          # (nb, npix)
    return float(np.sum((ivar * r * r)[:, keep]))

def roundtrip(a):
    out = np.empty_like(a)
    for ic in range(nc):
        alm = curvedsky.map2alm_healpix(np.array([a[:, ic]], dtype=rdt), niter=3, spin=0, lmax=lmax)
        alm = alm.astype(cdt, copy=False)
        alm = hp.smoothalm(alm, 0.0, inplace=False) if prec == "double" else alm
        out[:, ic] = alm_to_map(alm.astype(cdt, copy=False), nside, lmax, spin=0)[0]
    return out

def chi2_mask_as_run(a_rt):    # mimic get_sky_at_nu (float32 accumulation) + local_loglike
    tot = 0.0
    for ib in range(nb):
        sky = np.zeros(npix, dtype=np.float32 if prec == "single" else np.float64)
        for ic in range(nc):
            sky += (a_rt[:, ic].astype(rdt) * M[ib, ic]).astype(sky.dtype)
        res = d[ib, keep] - sky[keep]
        w = res * np.sqrt(ivar[ib, keep])
        tot += float(np.sum(w ** 2))
    return tot

c0 = chi2_mask(ahat)
ahat_rt = roundtrip(ahat)
B = chi2_mask_as_run(ahat_rt) - c0
print(f"chi2_min inside mask {c0:.1f} ({c0/(nb*keep.sum()):.3f} per band-pixel); round trip of the mean adds B = {B:.1f} = {B/k:.4f} k")
rng = np.random.default_rng(seed)
Xp, Xr, Xq = [], [], []
for i in range(ndraw):
    eta = rng.standard_normal((npix, nc))
    delta = np.linalg.solve(Lt, eta[..., None])[..., 0]
    a = ahat + delta
    Xp.append(chi2_mask(a) - c0)
    Xr.append(chi2_mask_as_run(roundtrip(a)) - c0)
    # fluctuation alone through the round trip: (P delta)^T A (P delta) inside the mask
    pd_ = roundtrip(delta)
    Xq.append(float(np.einsum("pi,pij,pj->", pd_[keep], A[keep], pd_[keep])))
    if i == 4:
        print(f"  5 draws in {time.time()-t0:.1f}s", flush=True)
Xp, Xr, Xq = map(np.array, (Xp, Xr, Xq))
def rep(name, X):
    n = len(X); v = X.var(ddof=1)
    print(f"{name:>46}: mean {X.mean():11.1f} = {X.mean()/k:6.3f} k;  sd {np.sqrt(v):7.1f};  var/(2k) = {v/(2*k):6.3f} +- {v/(2*k)*np.sqrt(2/(n-1)):.3f}")
rep("pixel-space draw (exact, expect 1 k, 1)", Xp)
rep("after alm round trip, as global_loglike sees", Xr)
rep("fluctuation alone after round trip", Xq)
print(f"corr(pixel, round trip) = {np.corrcoef(Xp, Xr)[0,1]:.3f}; ndraw {ndraw}; time {time.time()-t0:.1f}s")
np.savez(os.path.join(cache, f"rt_{prec}_{seed}.npz"), Xp=Xp, Xr=Xr, Xq=Xq, B=B, c0=c0, k=k)
