"""Addition 3 on simulations: the low-pass estimator at nside 4 (384 elements per band, basis to lmax 11).

Band 23 GHz (the anchor) and 353 GHz keep all modes; 33, 41 and 61 GHz enter through multipoles 2 to L (L = 5: 64 of
280 basis columns), built exactly as in e039c_lowpass_sim.py. CMB drawn from the LCDM covariance of
`lcdm_cmb_precision`, the same covariance as the prior. For each realization: all-scales fit and low-pass fit.

Then, on one realization, what a contamination does to each estimator (shift of the index in units of its error):
  (1) in a low-passed band, made of basis columns above the cut;      (2) in a low-passed band, outside span(Y);
  (3) in the anchor band, basis columns above the cut;                (4) in the anchor band, outside span(Y);
  (5) as (1) with the other projector, pinv(Y_k).
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

nreal = int(sys.argv[1]) if len(sys.argv) > 1 else 300
LCUT = int(sys.argv[2]) if len(sys.argv) > 2 else 5
toy = tn.Toy()
prec, cov = lcdm_cmb_precision(toy.WM)
Lc = np.linalg.cholesky(cov)
Y, Yp, ell = toy.WM.Y, toy.WM.Yp, toy.WM.ell
keep = ell <= LCUT
Yk, Ypk = Y[:, keep], Yp[keep]
LOW = (1, 2, 3)


def lowpass_ninv(Lmap):
    out = []
    for i in range(len(tn.FREQS)):
        if i in LOW:
            Cc = Lmap @ toy.covs[i] @ Lmap.T
            nl = Lmap.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Lmap
            out.append(0.5 * (nl + nl.T))
        else:
            out.append(toy.ninv[i])
    return out


ninv_lp = lowpass_ninv(Ypk)


def fit(dv, lowpass, Lmap=Ypk, ninvs=None, proj=None):
    if lowpass:
        proj = Yk @ Ypk if proj is None else proj
        dd = [proj @ d if i in LOW else d for i, d in enumerate(dv)]
        L = RegionLikelihood(toy.bands(dd, ninv_lp if ninvs is None else ninvs), toy.region, cmb_prior_prec=prec)
    else:
        L = RegionLikelihood(toy.bands(dv), toy.region, cmb_prior_prec=prec)
    bh, o = fit_profile(L, tn.BETA_TRUE)
    _, cs, Jp, T, _, ok = sandwich_plugin(L, bh)
    return bh, np.sqrt(np.diag(cs)), np.diag(T) / np.diag(Jp - T), o.success


print(f"# t3c: nside 4, anchor 23 GHz and 353 GHz at all scales, 33, 41, 61 GHz through multipoles 2 to {LCUT} ({int(keep.sum())} of {Y.shape[1]} columns); "
      f"truth {tn.BETA_TRUE}; {nreal} realizations", flush=True)
sig = toy.signal()
t0 = time.time()
rows = []
for r in range(nreal):
    rng = np.random.default_rng([2026, 3, r])
    a_c = Lc @ rng.standard_normal(tn.N)
    dv = [s + a_c + e for s, e in zip(sig, toy.noise(rng))]
    ra, rl = fit(dv, False), fit(dv, True)
    rows.append(np.concatenate([ra[0], ra[1], rl[0], rl[1], rl[2], [float(ra[3]), float(rl[3])]]))
    if r % 100 == 0:
        print(f"  realization {r}: all scales {ra[0].round(4)} +- {ra[1].round(4)} | low-pass {rl[0].round(4)} +- {rl[1].round(4)} [{time.time() - t0:.0f} s]", flush=True)
a = np.array(rows)
for name, est, sd in (("all scales", a[:, 0:3], a[:, 3:6]), (f"index from l <= {LCUT}", a[:, 6:9], a[:, 9:12])):
    zm, zs, c1, sem, ses = tn.zstats(est, sd, tn.BETA_TRUE)
    print(f"## {name}")
    print("   mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (est - tn.BETA_TRUE).mean(0)) + "  +- " + " ".join(f"{x:.4f}" for x in est.std(0, ddof=1) / np.sqrt(len(a))))
    print("   mean z                " + " ".join(f"{x:+.3f}" for x in zm) + f"  (se {sem:.3f})")
    print("   sd of z               " + " ".join(f"{x:.3f}" for x in zs) + f"  (se {ses:.3f})")
    print("   coverage of +-1 sd    " + " ".join(f"{x:.3f}" for x in c1) + "  (0.683 nominal)")
    print("   median quoted error   " + " ".join(f"{x:.4f}" for x in np.median(sd, axis=0)) + "; scatter of the estimate " + " ".join(f"{x:.4f}" for x in est.std(0, ddof=1)))
print("   low-pass: median T/J_c " + " ".join(f"{x:.4f}" for x in np.median(a[:, 12:15], axis=0)) + f"; converged all-scales {int(a[:, 15].sum())}, low-pass {int(a[:, 16].sum())} of {len(a)}")
dlt = a[:, 6:9] - a[:, 0:3]
pred = np.sqrt(np.clip(np.median(a[:, 9:12], axis=0) ** 2 - np.median(a[:, 3:6], axis=0) ** 2, 0, None))
print("## low-pass minus all scales, paired: mean " + " ".join(f"{x:+.4f}" for x in dlt.mean(0)) + " +- " + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1) / np.sqrt(len(a)))
      + "; sd " + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1)) + "; (sd_low^2 - sd_all^2)^(1/2) from the quoted errors " + " ".join(f"{x:.4f}" for x in pred))

# contamination tests on one realization
rng = np.random.default_rng([2026, 3, 10 ** 6])
a_c = Lc @ rng.standard_normal(tn.N)
dv = [s + a_c + e for s, e in zip(sig, toy.noise(rng))]
base_a, base_l = fit(dv, False), fit(dv, True)
Pout = np.eye(tn.N) - Y @ Yp
hi_cols = Y[:, ~keep]
c_hi = hi_cols @ rng.standard_normal(hi_cols.shape[1])
c_hi *= 0.3 / np.sqrt(np.mean(c_hi ** 2))                                   # 0.3 uK rms, 1.3 times the noise per element
c_out = Pout @ rng.standard_normal(tn.N)
c_out *= 0.3 / np.sqrt(np.mean(c_out ** 2))
alt = np.linalg.pinv(Yk)
ninv_alt = lowpass_ninv(alt)
base_alt = fit(dv, True, ninvs=ninv_alt, proj=Yk @ alt)
print("## contamination of 0.3 uK rms added to one band: shift of the index / its quoted error (high, mid, plane)")
for label, band, c in (("(1) 33 GHz (low-passed), basis columns above the cut", 1, c_hi), ("(2) 33 GHz (low-passed), outside span(Y)", 1, c_out),
                       ("(3) 23 GHz (anchor), basis columns above the cut", 0, c_hi), ("(4) 23 GHz (anchor), outside span(Y)", 0, c_out)):
    dc = [d + (c if i == band else 0.0) for i, d in enumerate(dv)]
    sa, sl_ = fit(dc, False), fit(dc, True)
    line = f"   {label:55s} all scales " + " ".join(f"{x:+8.3f}" for x in (sa[0] - base_a[0]) / base_a[1]) + " | low-pass " + " ".join(f"{x:+.2e}" for x in (sl_[0] - base_l[0]) / base_l[1])
    if band == 1 and c is c_hi:
        s2 = fit(dc, True, ninvs=ninv_alt, proj=Yk @ alt)
        line += " | (5) low-pass with pinv(Y_k) " + " ".join(f"{x:+.3f}" for x in (s2[0] - base_alt[0]) / base_alt[1])
    print(line)
print(f"   estimates on the clean data: all scales {base_a[0].round(4)}, low-pass {base_l[0].round(4)}, low-pass with pinv(Y_k) {base_alt[0].round(4)} +- {base_alt[1].round(4)}")
np.savez(f"nextgen/reviews/stats_review_scripts_2026-11-11b/out_t3c_L{LCUT}.npz", rows=a)
print(f"# done, {time.time() - t0:.0f} s")
