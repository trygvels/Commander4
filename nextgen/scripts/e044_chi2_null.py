"""E044: the distribution of chi2_min and of the chi-square per element and band of each region under the model.

Reads realizations drawn by `e040j_sim_n16_stages.py draw` (regions, mask, covariances and prior of the map fit), does
the LCDM-prior profile fit of each and prints what `e044_fit.py` prints for the maps: chi2_min over its nominal
degrees of freedom (bands * n - 2 n - regions; the CMB amplitude is held by its prior, whose penalty is inside
chi2_min) and sum_b r_R^T (N^-1 r)_R / (n_R n_bands) per region.

    python -I nextgen/scripts/e044_chi2_null.py --dir ~/data/c4_tools/sim_n8_named --nreal 40 --out FILE.npz
"""
import argparse
import gc
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e040j_sim_n16_stages as S  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import RegionLikelihood  # noqa: E402
from lowres.mask import restrict_prior  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--nreal", type=int, default=40)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    a.mode = "blocknull"
    t0 = time.time()
    z, ns, region, block, Bs, ninv, low, keep = S.load(a)
    assert not keep.all()
    kw = dict(cmb_prior_sd=None, cmb_prior_prec=restrict_prior(np.load(os.path.join(a.dir, "cmb_cov.npy")), keep))
    n, nreg, nb = region.size, 4, len(Bs)
    nper = np.bincount(region, minlength=nreg)
    dof = nb * n - 2 * n - nreg
    print(f"# E044 chi-square under the model, nside {ns}: {n} elements kept, nominal dof {dof}, sd of the ratio {np.sqrt(2.0 / dof):.4f}", flush=True)
    R, C, P = [], [], []
    for r_ in range(min(a.nreal, z["data"].shape[0])):
        L = RegionLikelihood(S.bands_for(z, Bs, ninv, low, r_, keep), region, **kw)
        bh, res = fit_profile(L, S.BETA_TRUE)
        sol = sandwich_plugin(L, bh)[4]
        rb = L.residuals(sol)
        creg = np.zeros(nreg)
        for b, r in zip(L.bands, rb):
            creg += np.bincount(region, weights=r * b.apply(r), minlength=nreg)
        chi2 = 2.0 * res.fun
        R.append(chi2 / dof); C.append(creg / (nper * nb)); P.append(chi2 - creg.sum())
        print(f"real {r_} chi2_min/dof {R[-1]:.4f} | per element and band " + " ".join(f"{x:.4f}" for x in C[-1]) + f" | prior penalty {P[-1]:.1f} | conv {res.success} [{time.time() - t0:.0f} s]", flush=True)
        del L, sol, rb
        gc.collect()
        if a.out:
            np.savez(a.out, ratio=np.array(R), creg=np.array(C), penalty=np.array(P), dof=dof, nper=nper)
    R, C, P = np.array(R), np.array(C), np.array(P)
    print(f"## {len(R)} realizations: chi2_min/dof mean {R.mean():.4f} +- {R.std(ddof=1) / np.sqrt(len(R)):.4f}, sd {R.std(ddof=1):.4f} (theory {np.sqrt(2.0 / dof):.4f}); "
          f"prior penalty mean {P.mean():.1f}")
    print("## per element and band, mean: " + " ".join(f"{x:.4f}" for x in C.mean(0)) + "; sd: " + " ".join(f"{x:.4f}" for x in C.std(0, ddof=1)))


if __name__ == "__main__":
    main()
