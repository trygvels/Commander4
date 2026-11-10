"""Check 3: coverage of intervals for the profile optimum.

For each noise realization: beta_hat = fit_profile; then
  truth-sandwich   J^-1 (J+T) J^-1 with J, T from the reviewer's explicit matrices at the true beta and amplitudes
  truth-naive      J^-1 (same J)
  plugin           fit.sandwich_plugin(L, beta_hat): J_c^-1 J(a_hat) J_c^-1, J_c = J(a_hat) - T
  plugin-curv      J_c^-1 (first return value of sandwich_plugin)
  plain            fit.sandwich(L, beta_hat): J(a_hat)^-1 (J(a_hat)+T) J(a_hat)^-1 (the biased plug-in)
Coverage = fraction of |beta_hat - beta_true| < k sd for k = 1, 2 (nominal 0.683, 0.954), pooled over regions.

usage: check_03_coverage.py <config> <nrep>
"""
import sys
import time

import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info
from lowres.fit import fit_profile, sandwich, sandwich_plugin

np.set_printoptions(precision=4, linewidth=170, suppress=True)
CONFIGS = {
    # name: (problem kwargs, likelihood kwargs, reference kwargs, cmb mode, cmb sd of the truth)
    "S1": (dict(n=320, nreg=2, amp=0.5, noise_slope=2.5), dict(use_dust=False, use_cmb=False),
           dict(use_dust=False, use_cmb=False), "none", 0.0),
    "S2": (dict(n=320, nreg=2, amp=0.35, noise_slope=2.5), dict(use_dust=False, use_cmb=False),
           dict(use_dust=False, use_cmb=False), "none", 0.0),
    "S3": (dict(n=1280, nreg=2, amp=0.35, noise_slope=2.5), dict(use_dust=False, use_cmb=False),
           dict(use_dust=False, use_cmb=False), "none", 0.0),
    "S4": (dict(n=96, nreg=2, amp=0.7, noise_slope=2.5), dict(use_dust=False, use_cmb=False),
           dict(use_dust=False, use_cmb=False), "none", 0.0),
    "F1": (dict(n=384, nreg=2, amp=0.6, noise_slope=2.5), dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "redraw", 0.3),
    "F2": (dict(n=384, nreg=2, amp=0.6, noise_slope=2.5), dict(cmb_prior_sd=1.0), dict(cmb_sd=1.0), "fixed", 1.0),
    "F3": (dict(n=384, nreg=2, amp=0.6, noise_slope=2.5), dict(cmb_prior_sd=1.0), dict(cmb_sd=1.0), "redraw", 1.0),
    "F4": (dict(n=384, nreg=2, amp=0.6, noise_slope=2.5), dict(), dict(), "fixed", 1.0),
    # WMAP/LFI-like toy of check_02b: 8 bands, noise flat in frequency, CMB prior sd 0.3
    #   F5: the CMB in the data has sd 0.6 (prior variance too small by 4); F6: no CMB in the data; F7: matched (control)
    "F5": (dict(n=384, nreg=2, amp=2.5, amp_d=3.0, noise353=0.5, nus=(23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0)),
           dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "redraw", 0.6),
    "F6": (dict(n=384, nreg=2, amp=2.5, amp_d=3.0, noise353=0.5, nus=(23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0)),
           dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "fixed", 0.0),
    "F7": (dict(n=384, nreg=2, amp=2.5, amp_d=3.0, noise353=0.5, nus=(23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0)),
           dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "redraw", 0.3),
}
name, nrep = sys.argv[1], int(sys.argv[2])
pk, lk, rk, mode, sc = CONFIGS[name]
prob = Problem(dense=(0, 2, 5), seed=7, noise=1.0, amp_c=1.0, **pk)
n, nreg, truth = prob.n, prob.nreg, prob.beta_true
zero = np.zeros(n)
a_c_fixed = prob.a_c * (sc / prob.a_c.std()) if sc > 0 else zero
info = ref_info(prob, truth, prob.a_s, **rk)
J, T = info["J"], info["T"]
Ji = np.linalg.inv(J)
sd_sw, sd_nv = np.sqrt(np.diag(Ji @ (J + T) @ Ji)), np.sqrt(np.diag(Ji))
print(f"config {name}: n={n} nreg={nreg} likelihood {lk} cmb truth '{mode}' sd {sc}  nrep={nrep}")
print("J diag", np.diag(J), " T/J", np.diag(T) / np.diag(J), " sd sandwich", sd_sw, " sd naive", sd_nv)

rng = np.random.default_rng(99)
rows, nbound, nfail = [], 0, 0
t0 = time.time()
sync_only = (lk.get("use_dust") is False)
for it in range(nrep):
    if mode == "none":
        mean = prob.mean(a_d=zero, a_c=zero)
    elif mode == "fixed":
        mean = prob.mean(a_c=a_c_fixed)
    else:
        mean = prob.mean(a_c=sc * rng.standard_normal(n))
    L = RegionLikelihood(prob.banddata(mean + prob.noise(rng)), prob.region, **lk)
    bh, res = fit_profile(L, truth)
    nbound += int(np.any((bh < -4.6 + 1e-6) | (bh > -1.6 - 1e-6)))
    nfail += int(not res.success)
    cc, cs, Jp, Th, _, ok = sandwich_plugin(L, bh)
    Jpi = np.linalg.inv(Jp)
    plain = np.sqrt(np.abs(np.diag(Jpi @ (Jp + Th) @ Jpi)))
    rows.append(np.concatenate([bh - truth, np.sqrt(np.diag(cs)), np.sqrt(np.diag(cc)), plain, ok.astype(float),
                                np.diag(Jp), np.diag(Th)]))
rows = np.array(rows)
err, sd_pl, sd_cc, sd_plain, ok, Jp_d, T_d = [rows[:, i * nreg:(i + 1) * nreg] for i in range(7)]
print(f"elapsed {time.time() - t0:.0f} s; fits at a bound: {nbound}; optimizer not converged: {nfail}")
print("mean(beta_hat - truth)       ", err.mean(0), "+-", err.std(0, ddof=1) / np.sqrt(nrep),
      "  in units of sandwich sd:", err.mean(0) / sd_sw)
print("empirical sd / truth-sandwich", err.std(0, ddof=1) / sd_sw, "+-", 1 / np.sqrt(2 * nrep), " / truth-naive", err.std(0, ddof=1) / sd_nv)
print("rms error / truth-sandwich   ", np.sqrt((err ** 2).mean(0)) / sd_sw)
print("mean J(a_hat) at beta_hat / (J+T)", Jp_d.mean(0) / np.diag(J + T), "  mean T(beta_hat)/T", T_d.mean(0) / np.diag(T))
print("fraction with all regions ok ", ok.all(1).mean(), "  per region", ok.mean(0))


def cover(z, sel=None):
    z = z if sel is None else z[sel]
    m = z.size
    c1, c2 = np.mean(np.abs(z) < 1), np.mean(np.abs(z) < 2)
    return (f"68.3%: {c1:.3f} +- {np.sqrt(c1 * (1 - c1) / m):.3f}   95.4%: {c2:.3f} +- {np.sqrt(c2 * (1 - c2) / m):.3f}   "
            f"sd(z) {z.std():.3f} +- {1 / np.sqrt(2 * m):.3f}")


okm = ok.astype(bool)
np.savez(f"/Users/trygvels/code/Commander4/nextgen/experiments/E039/stats_review_scripts/out_03_rows_{name}.npz",
         err=err, sd_plugin=sd_pl, sd_curv=sd_cc, sd_plain=sd_plain, ok=okm, Jhat=Jp_d, That=T_d, J=J, T=T, sd_sw=sd_sw, sd_nv=sd_nv)
print("coverage (pooled over regions and realizations)")
print("  truth-sandwich         ", cover(err / sd_sw))
print("  truth-naive J^-1       ", cover(err / sd_nv))
print("  plugin (ok only)       ", cover(err / sd_pl, okm))
miss1 = np.mean(np.where(okm, np.abs(err / sd_pl) < 1, False))
miss2 = np.mean(np.where(okm, np.abs(err / sd_pl) < 2, False))
print(f"  plugin, not-ok counted as a miss:  68.3%: {miss1:.3f}   95.4%: {miss2:.3f}")
print("  plugin-curv J_c^-1 (ok)", cover(err / sd_cc, okm))
print("  plain plug-in          ", cover(err / sd_plain))
print("median plugin sd / truth-sandwich sd (ok only)", np.nanmedian(np.where(okm, sd_pl / sd_sw, np.nan), axis=0),
      "  16th, 84th percentile", np.nanpercentile(np.where(okm, sd_pl / sd_sw, np.nan), [16, 84], axis=0).ravel())
print("corr(|err|, plugin sd) per region", [round(float(np.corrcoef(np.abs(err[okm[:, r], r]), sd_pl[okm[:, r], r])[0, 1]), 3) for r in range(nreg)])
print("corr(err, plugin sd) per region  ", [round(float(np.corrcoef(err[okm[:, r], r], sd_pl[okm[:, r], r])[0, 1]), 3) for r in range(nreg)])
