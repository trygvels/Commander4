#!/usr/bin/env python
"""E026 analysis: index-field chains on the realistic-noise PySM sky against the noise-free
region-fit reference, per template-S/N quartile, plus the hier hyperparameters (scalar mu or mu
map against the realized template error) and the ESS of every block.

usage: e026_analyze.py --dir DIR --refs E023/scan_n1_refs.npz --tags TAG [TAG ...] [--burn 500]
       [--cs-scatter 0.3 --cs-fwhm 10 --cs-seed 11]   (the cs template's error field, for the mu map)
"""
import argparse, glob
import numpy as np, healpy as hp, h5py


def ess(x):
    x = x - x.mean(); n = len(x)
    if n < 4 or np.all(x == 0):
        return float(n)
    c = np.correlate(x, x, "full")[n - 1:] / np.sum(x * x); s = 0.0
    for k in range(0, n - 1, 2):
        p = c[k] + c[k + 1]
        if p <= 0:
            break
        s += p
    return n / (2 * s - 1) if s > 0.5 else float(n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--refs", required=True)
    ap.add_argument("--tags", nargs="+", required=True); ap.add_argument("--burn", type=int, default=500)
    ap.add_argument("--cs-scatter", type=float, default=0.3); ap.add_argument("--cs-fwhm", type=float, default=10.0)
    ap.add_argument("--cs-seed", type=int, default=11)
    a = ap.parse_args()
    t = np.load(f"{a.dir}/truth.npz"); nside = int(t["nside"]); npix = hp.nside2npix(nside)
    refs = np.load(a.refs)
    ref = {"sync": refs["bs_fixT"], "dust": refs["bd_fixT"]}
    quart = {"sync": refs["qs16"], "dust": refs["qd16"]}
    truth16 = {"sync": hp.ud_grade(np.asarray(t["beta_s"], float), 16),
               "dust": hp.ud_grade(np.asarray(t["beta_eff3"], float), 16)}
    # the cs template's error field (same construction as the FITS writer): exp(-f) per region
    lmax = 3 * nside - 1; sig_b = np.radians(a.cs_fwhm) / np.sqrt(8 * np.log(2)); ell = np.arange(lmax + 1)
    cl = np.exp(-ell * (ell + 1) * sig_b**2); np.random.seed(a.cs_seed); expf = {}
    for k in ("sync", "dust"):
        f = hp.synfast(cl, nside, lmax=lmax, new=True); f = f / f.std() * a.cs_scatter
        expf[k] = hp.ud_grade(np.exp(-f), 16)       # region mean of exp(-f) (the realized a/g)
    comp_of = {"sync": "Synchrotron", "dust": "ThermalDust"}
    for tag in a.tags:
        files = sorted(glob.glob(f"{a.dir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))
        post = files[a.burn:]
        print(f"\n== {tag}: {len(files)} files, {len(post)} post-burn")
        if len(post) < 10:
            print("   too few post-burn draws"); continue
        beta = {k: [] for k in ("sync", "dust")}; mu = {k: [] for k in ("sync", "dust")}
        ls = {k: [] for k in ("sync", "dust")}; grads = []
        for fn in post:
            with h5py.File(fn) as f:
                for k in ("sync", "dust"):
                    beta[k].append(f[f"comps/{k}/index_field/beta_map"][()])
                    g = f[f"mcmc/sample_index_fields/state/{comp_of[k]}/amp_prior"]
                    mu[k].append(np.atleast_1d(g["mu"][()])); ls[k].append(float(g["log_s"][()]))
                if "mcmc/sample_index_fields/nuts/n_leapfrog" in f:
                    grads.append(float(f["mcmc/sample_index_fields/nuts/n_leapfrog"][()]))
        if grads:
            print(f"   gradients per iteration (post-burn mean) {np.mean(grads):.1f}")
        for k in ("sync", "dust"):
            b = np.array(beta[k]); m = b.mean(0); sd = b.std(0); e = np.array([ess(b[:, r]) for r in range(b.shape[1])])
            z = (m - ref[k]) / sd; zt = (m - truth16[k]) / sd
            print(f"   {k}: region ESS min {e.min():.0f} 5pct {np.percentile(e, 5):.0f} median {np.median(e):.0f}; "
                  f"max |mean - ref| {np.abs(m - ref[k]).max():.3f}")
            for q in range(4):
                sel = quart[k] == q
                print(f"      Q{q+1}: rmse vs ref {np.sqrt(np.mean((m - ref[k])[sel]**2)):.4f} (vs region-mean truth "
                      f"{np.sqrt(np.mean((m - truth16[k])[sel]**2)):.4f}); z vs ref median {np.median(z[sel]):+.2f} "
                      f"MAD {np.median(np.abs(z[sel] - np.median(z[sel]))):.2f} |z|<1 {np.mean(np.abs(z[sel]) < 1):.2f} "
                      f"|z|<2 {np.mean(np.abs(z[sel]) < 2):.2f}; sd median {np.median(sd[sel]):.4f}; z vs truth median {np.median(zt[sel]):+.2f}")
            M = np.array(mu[k]); L = np.exp(np.array(ls[k]))
            print(f"      s: mean {L.mean():.4f} sd {L.std():.4f} ESS {ess(L):.0f}")
            if M.shape[1] == 1:
                print(f"      mu (scalar): mean {M[:, 0].mean():.4f} sd {M[:, 0].std():.4f} ESS {ess(M[:, 0]):.0f}")
            else:
                mm = M.mean(0); ms = M.std(0); me = np.array([ess(M[:, r]) for r in range(M.shape[1])])
                nside_mu = hp.npix2nside(M.shape[1]); target = hp.ud_grade(expf[k], nside_mu) if nside_mu != 16 else expf[k]
                zmu = (mm - target) / ms
                print(f"      mu map ({M.shape[1]}): ESS min {me.min():.0f} median {np.median(me):.0f}; mean {mm.mean():.4f}; "
                      f"against the realized exp(-f) region means: correlation {np.corrcoef(mm, target)[0, 1]:.3f}, "
                      f"rms diff {np.sqrt(np.mean((mm - target)**2)):.4f}, z median {np.median(zmu):+.2f} MAD "
                      f"{np.median(np.abs(zmu - np.median(zmu))):.2f}")
                for q in range(4):
                    sel = quart[k] == q if nside_mu == 16 else np.ones(M.shape[1], bool)
                    print(f"         Q{q+1}: corr {np.corrcoef(mm[sel], target[sel])[0, 1]:.3f} rms diff "
                          f"{np.sqrt(np.mean((mm - target)[sel]**2)):.4f} sd median {np.median(ms[sel]):.4f}")


if __name__ == "__main__":
    main()
