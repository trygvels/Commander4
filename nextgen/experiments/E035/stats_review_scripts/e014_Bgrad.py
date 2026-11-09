"""Does the round-trip offset B(beta) = chi2(P a_hat(beta)) - chi2_min(beta) tilt the profile the legacy
move sees? Finite differences of chi2_min and B around the exact ridge mean, double precision.
Predicted shift of the legacy mean relative to the ridge mean: -(1/2) C_ridge grad B."""
import os, sys, time
import numpy as np
from pixell import curvedsky
from commander4.parameters.parse import load_params
from commander4.sky.comp_list import CompList
from commander4.math_utils.sht import alm_to_map
pf, cache = sys.argv[1], sys.argv[2]
params, _, _ = load_params(pf)
names = [n for n in params.compsep.bands if params.compsep.bands[n].enabled]
d, ivar, freqs = [], [], []
for n in names:
    z = np.load(os.path.join(cache, f"band_{n}.npz"))
    d.append(z["d"]); ivar.append(z["ivar"]); freqs.append(float(z["freq"]))
d, ivar, freqs = np.array(d), np.array(ivar), np.array(freqs)
keep = np.load(os.path.join(cache, "mask.npy"))
comps = list(CompList.init_from_params(params.components, params).split_for_eval_pol("I"))
nc = len(comps); k = nc * int(keep.sum())
def both(bs, bd, prec=np.float64):
    for c in comps:
        if type(c).__name__ == "Synchrotron": c.beta = bs
        if type(c).__name__ == "ThermalDust": c.beta = bd
    M = np.stack([c.get_sed(freqs) for c in comps], axis=1).astype(np.float64)
    A = np.einsum("bi,bp,bj->pij", M, ivar, M)
    b = np.einsum("bi,bp->pi", M, ivar * d)
    ahat = np.linalg.solve(A, b[..., None])[..., 0]
    def chi2(a):
        r = d - M @ a.T
        return float(np.sum((ivar * r * r)[:, keep]))
    art = np.empty_like(ahat)
    cdt = np.complex128 if prec == np.float64 else np.complex64
    for ic in range(nc):
        alm = curvedsky.map2alm_healpix(np.array([ahat[:, ic]], dtype=prec), niter=3, spin=0, lmax=191)
        art[:, ic] = alm_to_map(alm.astype(cdt, copy=False), 64, 191, spin=0)[0]
    c0 = chi2(ahat)
    return c0, chi2(art) - c0
bs0, bd0, ss, sd = -3.34244, 1.36539, 0.00036, 0.000043
h = 5.0
t0 = time.time()
pts = {(i, j): both(bs0 + i * h * ss, bd0 + j * h * sd) for i in (-1, 0, 1) for j in (-1, 0, 1)}
c = {kk: v[0] for kk, v in pts.items()}; B = {kk: v[1] for kk, v in pts.items()}
# derivatives in units of ridge sd (u_s, u_d)
def grad_hess(f):
    g = np.array([(f[1, 0] - f[-1, 0]) / (2 * h), (f[0, 1] - f[0, -1]) / (2 * h)])
    H = np.array([[(f[1, 0] - 2 * f[0, 0] + f[-1, 0]) / h ** 2, (f[1, 1] - f[1, -1] - f[-1, 1] + f[-1, -1]) / (4 * h * h)],
                  [0.0, (f[0, 1] - 2 * f[0, 0] + f[0, -1]) / h ** 2]])
    H[1, 0] = H[0, 1]
    return g, H
gc, Hc = grad_hess(c); gB, HB = grad_hess(B)
C = np.linalg.inv(Hc / 2.0)                 # ridge covariance in units of the quoted ridge sds
print(f"chi2_min {c[0,0]:.1f}; B = {B[0,0]:.1f} = {B[0,0]/k:.4f} k (double precision)")
print(f"ridge covariance from the chi2_min Hessian, in quoted ridge sd units: sd_s {np.sqrt(C[0,0]):.3f}, sd_d {np.sqrt(C[1,1]):.3f}, "
      f"corr {C[0,1]/np.sqrt(C[0,0]*C[1,1]):.3f}; grad chi2_min at the quoted mean (per sd): {gc[0]:.2f}, {gc[1]:.2f}")
print(f"grad B per ridge sd: d/ds {gB[0]:.2f}, d/dd {gB[1]:.2f};  Hessian of B / Hessian of chi2_min (diag): "
      f"{HB[0,0]/Hc[0,0]:.4f}, {HB[1,1]/Hc[1,1]:.4f}")
shift = -np.linalg.solve(Hc + HB, gB)
print(f"predicted legacy-minus-ridge mean shift from B alone: beta_s {shift[0]:+.2f} ridge sd, beta_d {shift[1]:+.2f} ridge sd")
print("measured (E014 log, burn 400): beta_s -2.7 (chain 1), -0.9 (chain 2); beta_d +3.0 (chain 1), +3.0 (chain 2)")
bsingle = both(bs0, bd0, prec=np.float32)[1]
print(f"B in single precision round trip: {bsingle:.1f}; time {time.time()-t0:.1f}s")
