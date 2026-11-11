"""Addition 1, claim (c): profile optimum and debiased sandwich with the LCDM prior of `lcdm_cmb_precision`, nside 4.

Arms (same synchrotron, dust and noise draws):
  A  CMB drawn from the returned covariance, fitted with the returned precision                 (the claim)
  B  the same data fitted with a white prior of the same rms per element                        (wrong shape)
  C  CMB drawn from the returned covariance plus a white part of 0.44 of its variance, the share that a
     pixel-averaged LCDM sky has above the band limit at nside 16 (t1b), fitted with the returned precision
  D  the data of C fitted with the returned covariance plus that white part as the prior          (repaired prior)
"""
import sys
import time

import numpy as np

sys.path.insert(0, "nextgen")
sys.path.insert(0, "nextgen/reviews/stats_review_scripts_2026-11-11b")
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import RegionLikelihood  # noqa: E402
from lowres.window import lcdm_cmb_precision  # noqa: E402
import toy_n4 as tn  # noqa: E402

nreal = int(sys.argv[1]) if len(sys.argv) > 1 else 200
toy = tn.Toy()
prec, cov = lcdm_cmb_precision(toy.WM)
Lc = np.linalg.cholesky(cov)
rms = np.sqrt(np.mean(np.diag(cov)))
extra_var = 0.44 * np.mean(np.diag(cov))
prec_fix = np.linalg.inv(cov + extra_var * np.eye(tn.N))
print(f"# t1d: nside 4, {tn.N} elements, regions (high, mid, plane) with {np.bincount(toy.region).tolist()} elements, truth {tn.BETA_TRUE}; "
      f"LCDM prior rms {rms:.3f} uK per element; noise sd {[round(float(np.sqrt(np.mean(np.diag(c)))), 3) for c in toy.covs]}", flush=True)
sig = toy.signal()
t0 = time.time()
res = {k: [] for k in "ABCD"}
for r in range(nreal):
    rng = np.random.default_rng([2026, 1, r])
    a_c = Lc @ rng.standard_normal(tn.N)
    nz = toy.noise(rng)
    a_x = np.sqrt(extra_var) * rng.standard_normal(tn.N)
    d1 = [s + a_c + e for s, e in zip(sig, nz)]
    d2 = [d + a_x for d in d1]
    for key, dv, kw in (("A", d1, dict(cmb_prior_prec=prec)), ("B", d1, dict(cmb_prior_sd=rms)),
                        ("C", d2, dict(cmb_prior_prec=prec)), ("D", d2, dict(cmb_prior_prec=prec_fix))):
        L = RegionLikelihood(toy.bands(dv), toy.region, **kw)
        bh, o = fit_profile(L, tn.BETA_TRUE)
        _, cs, Jp, T, _, ok = sandwich_plugin(L, bh)
        res[key].append(np.concatenate([bh, np.sqrt(np.diag(cs)), np.diag(T) / np.diag(Jp - T), [float(o.success)]]))
    if r in (0, 4) or r % 50 == 49:
        print(f"  realization {r}: A {res['A'][-1][:3].round(4)} +- {res['A'][-1][3:6].round(4)} [{time.time() - t0:.0f} s]", flush=True)
names = {"A": "A  CMB from the prior covariance, LCDM prior", "B": "B  same data, white prior of the same rms",
         "C": "C  CMB with 44 percent extra white variance, LCDM prior as coded", "D": "D  data of C, prior with the white part added"}
for key in "ABCD":
    a = np.array(res[key])
    est, sd = a[:, :3], a[:, 3:6]
    zm, zs, cov1, sem, ses = tn.zstats(est, sd, tn.BETA_TRUE)
    print(f"## {names[key]}: {len(a)} realizations, converged {int(a[:, -1].sum())}, median T/J_c " + " ".join(f"{x:.3f}" for x in np.median(a[:, 6:9], axis=0)))
    print("   mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (est - tn.BETA_TRUE).mean(0)) + "  +- " + " ".join(f"{x:.4f}" for x in est.std(0, ddof=1) / np.sqrt(len(a))))
    print("   mean z                " + " ".join(f"{x:+.3f}" for x in zm) + f"  (se {sem:.3f})")
    print("   sd of z               " + " ".join(f"{x:.3f}" for x in zs) + f"  (se {ses:.3f})")
    print("   coverage of +-1 sd    " + " ".join(f"{x:.3f}" for x in cov1) + "  (0.683 nominal)")
    print("   median quoted error   " + " ".join(f"{x:.4f}" for x in np.median(sd, axis=0)) + "; scatter of the estimate " + " ".join(f"{x:.4f}" for x in est.std(0, ddof=1)))
np.savez("nextgen/reviews/stats_review_scripts_2026-11-11b/out_t1d.npz", **{k: np.array(v) for k, v in res.items()})
print(f"# done, {time.time() - t0:.0f} s")
