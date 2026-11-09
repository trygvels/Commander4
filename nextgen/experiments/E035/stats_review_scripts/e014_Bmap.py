"""Where does the round-trip offset B live? Per-pixel contribution (Q a_hat)^T A (Q a_hat) inside the mask
by Galactic latitude, by component, and its sensitivity to the number of map2alm iterations."""
import os, sys
import numpy as np, healpy as hp
from pixell import curvedsky
from commander4.parameters.parse import load_params
from commander4.sky.comp_list import CompList
from commander4.math_utils.sht import alm_to_map
pf, cache = sys.argv[1], sys.argv[2]
params, _, _ = load_params(pf)
names = [n for n in params.compsep.bands if params.compsep.bands[n].enabled]
d, ivar, freqs = [], [], []
for n in names:
    z = np.load(os.path.join(cache, f"band_{n}.npz")); d.append(z["d"]); ivar.append(z["ivar"]); freqs.append(float(z["freq"]))
d, ivar, freqs = np.array(d), np.array(ivar), np.array(freqs)
keep = np.load(os.path.join(cache, "mask.npy"))
comps = list(CompList.init_from_params(params.components, params).split_for_eval_pol("I"))
for c in comps:
    if type(c).__name__ == "Synchrotron": c.beta = -3.34244
    if type(c).__name__ == "ThermalDust": c.beta = 1.36539
cn = [type(c).__name__ for c in comps]
M = np.stack([c.get_sed(freqs) for c in comps], axis=1).astype(np.float64)
A = np.einsum("bi,bp,bj->pij", M, ivar, M); b = np.einsum("bi,bp->pi", M, ivar * d)
ahat = np.linalg.solve(A, b[..., None])[..., 0]; nc = len(comps); k = nc * int(keep.sum())
def rt(a, niter):
    out = np.empty_like(a)
    for ic in range(nc):
        alm = curvedsky.map2alm_healpix(np.array([a[:, ic]]), niter=niter, spin=0, lmax=191)
        out[:, ic] = alm_to_map(alm, 64, 191, spin=0)[0]
    return out
lat = 90.0 - np.degrees(hp.pix2ang(64, np.arange(12 * 64 ** 2))[0])
for niter in (0, 3, 10, 30):
    q = ahat - rt(ahat, niter)
    bp = np.einsum("pi,pij,pj->p", q, A, q)
    line = f"niter {niter:2d}: B = {bp[keep].sum():10.0f} = {bp[keep].sum()/k:.3f} k"
    if niter == 3:
        bp3, q3 = bp, q
    print(line)
print("niter 3, B per kept pixel by |b| (deg):")
for lo, hi in [(0, 10), (10, 20), (20, 30), (30, 45), (45, 60), (60, 90)]:
    m = keep & (np.abs(lat) >= lo) & (np.abs(lat) < hi)
    if m.sum():
        print(f"  {lo:2d}-{hi:2d}: {m.sum():6d} px, mean {bp3[m].mean():8.2f}, median {np.median(bp3[m]):7.2f}, share of B {bp3[m].sum()/bp3[keep].sum():.3f}")
srt = np.sort(bp3[keep])[::-1]
print(f"share of B in the top 1% / 10% of kept pixels: {srt[:len(srt)//100].sum()/srt.sum():.3f} / {srt[:len(srt)//10].sum()/srt.sum():.3f}")
print("rms of the round-trip residual of each mean amplitude map inside the mask, relative to its marginal sd sqrt((A^-1)_cc):")
Ainv = np.linalg.inv(A[keep])
for ic in range(nc):
    print(f"  {cn[ic]:>13}: {np.sqrt(np.mean(q3[keep, ic] ** 2 / Ainv[:, ic, ic])):.3f}")
