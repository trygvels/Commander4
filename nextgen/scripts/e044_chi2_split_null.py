"""E044 H14: the sampling distribution of the chi-square split by region and direction, under the model.

For each simulated realization (drawn by `e040j_sim_n16_stages.py draw` with the regions, mask, covariances and prior
of the map fit): the LCDM-prior profile fit and the flat-CMB-prior profile fit, and for each the residual chi-square
by region, sum_b r_R^T (N^-1 r)_R, as `e044_fit.py` computes it for the maps. Under the model the difference of the two
totals is chi-square with one degree of freedom per element (the direction of a CMB spectrum); the script measures how
the per-region parts of that difference and of the flat-fit residual scatter.

    python -I nextgen/scripts/e044_chi2_split_null.py --dir ~/data/c4_tools/sim_n8_named --nreal 40 --out FILE.npz
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


def by_region(L, beta, region, nreg):
    sol = sandwich_plugin(L, beta)[4]
    rb = L.residuals(sol)
    c = np.zeros(nreg)
    for b, r in zip(L.bands, rb):
        c += np.bincount(region, weights=r * b.apply(r), minlength=nreg)
    return c


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
    prec = restrict_prior(np.load(os.path.join(a.dir, "cmb_cov.npy")), keep)
    n, nreg, nb = region.size, 4, len(Bs)
    nper = np.bincount(region, minlength=nreg).astype(float)
    print(f"# E044 chi-square split under the model, nside {ns}: elements per region {nper.astype(int)}; expected sd of the difference over n_R "
          + " ".join(f"{x:.3f}" for x in np.sqrt(2.0 / nper)) + "; of the flat-fit chi-square over 2 n_R " + " ".join(f"{x:.3f}" for x in np.sqrt(1.0 / nper)), flush=True)
    CL, CF, TL, TF = [], [], [], []
    for r_ in range(min(a.nreal, z["data"].shape[0])):
        bands = S.bands_for(z, Bs, ninv, low, r_, keep)
        L = RegionLikelihood(bands, region, cmb_prior_sd=None, cmb_prior_prec=prec)
        bl, rl = fit_profile(L, S.BETA_TRUE)
        cl = by_region(L, bl, region, nreg)
        del L
        gc.collect()
        L = RegionLikelihood(bands, region, cmb_prior_sd=None)
        bf, rf = fit_profile(L, S.BETA_TRUE)
        cf = by_region(L, bf, region, nreg)
        del L
        gc.collect()
        CL.append(cl); CF.append(cf); TL.append(2.0 * rl.fun); TF.append(2.0 * rf.fun)
        pen = TL[-1] - cl.sum()
        d = (cl - cf) / (nper - pen * nper / nper.sum())
        print(f"real {r_} difference over expectation " + " ".join(f"{x:.3f}" for x in d) + " | flat over 2 n_R " + " ".join(f"{x:.3f}" for x in cf / (2 * nper))
              + f" | totals {TL[-1]:.1f} {TF[-1]:.1f} | conv {rl.success} {rf.success} [{time.time() - t0:.0f} s]", flush=True)
        if a.out:
            np.savez(a.out, cl=np.array(CL), cf=np.array(CF), tl=np.array(TL), tf=np.array(TF), nper=nper)
    CL, CF, TL, TF = np.array(CL), np.array(CF), np.array(TL), np.array(TF)
    pen = (TL - CL.sum(1))[:, None]
    D = (CL - CF) / (nper - pen * nper / nper.sum())
    F = CF / (2 * nper)
    m = D.shape[0]
    print(f"## {m} realizations. Difference (LCDM minus flat) by region over its expectation: mean " + " ".join(f"{x:.3f}" for x in D.mean(0))
          + "; sd " + " ".join(f"{x:.3f}" for x in D.std(0, ddof=1)) + "; expected sd " + " ".join(f"{x:.3f}" for x in np.sqrt(2.0 / nper)))
    print("## flat-fit chi-square by region over 2 n_R: mean " + " ".join(f"{x:.3f}" for x in F.mean(0)) + "; sd " + " ".join(f"{x:.3f}" for x in F.std(0, ddof=1))
          + "; expected sd " + " ".join(f"{x:.3f}" for x in np.sqrt(1.0 / nper)))
    print("## correlation of the two, per region: " + " ".join(f"{np.corrcoef(D[:, j], F[:, j])[0, 1]:+.2f}" for j in range(nreg))
          + f"; total difference over n: mean {((TL - TF) / n).mean():.4f}, sd {((TL - TF) / n).std(ddof=1):.4f} (expected {np.sqrt(2.0 / n):.4f})")


if __name__ == "__main__":
    main()
