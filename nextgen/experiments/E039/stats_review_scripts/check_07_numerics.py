"""Check 7: numerical fragility.

(a) chi2 = d^T N^-1 d - rhs^T A^-1 rhs (the form in `solve`) against the QR reference and against the residual
    form sum_b r_b^T N_b^-1 r_b + prior penalty, as the signal-to-noise and the conditioning of A grow.
(b) analytic profile gradient against finite differences of the QR reference in the same configurations.
(c) singular A (an element with no data in any band; flat priors): what happens.
(d) information() with the synchrotron prior switched on, against explicit matrices.
(e) an inverse covariance that is not exactly symmetric.
(f) region labels with a gap (an empty region).
(g) sandwich_plugin when J_c is barely positive or negative.
"""
import numpy as np

from rv_common import Problem, RegionLikelihood, BandData, fd_grad, ref_info, ref_solve
from lowres.fit import fit_profile, sandwich, sandwich_plugin

np.set_printoptions(precision=5, linewidth=170)
beta = np.array([-3.0, -3.2, -2.9])

print("== (a), (b) cancellation and conditioning ==")
print(f"{'case':34s} {'cond(A)':>9s} {'dNd':>10s} {'chi2 (QR)':>12s} {'code - QR':>11s} {'resid - QR':>11s} {'grad rel err':>13s}")
for label, pk in (("amp 8, bands as in the tests", dict(amp=8.0)),
                  ("amp 1e3", dict(amp=1e3)), ("amp 1e5", dict(amp=1e5)), ("amp 1e7", dict(amp=1e7)),
                  ("amp 1e3, bands 22.8 23.0 23.2 33 41 353", dict(amp=1e3, nus=(22.8, 23.0, 23.2, 33.0, 41.0, 353.0))),
                  ("amp 1e3, noise 353 = 1e-4", dict(amp=1e3, noise353=1e-4)),
                  ("amp 1e5, dust 1e5, noise 353 = 1e-3", dict(amp=1e5, amp_d=1e5, noise353=1e-3))):
    p = Problem(n=48, nreg=3, dense=(0, 2, 5), seed=1, **pk)
    d = p.mean() + p.noise(np.random.default_rng(5))
    L = RegionLikelihood(p.banddata(d), p.region)
    sol = L.solve(beta)
    ref = ref_solve(p, d, beta)
    res = L.residuals(sol)
    chi_res = sum(r @ b.apply(r) for r, b in zip(res, L.bands))
    A = L.assemble(beta)[0]
    _, g = L.profile(beta, grad=True)
    gfd = fd_grad(lambda b: 0.5 * ref_solve(p, d, b)["chi2"], beta, h=1e-6)
    print(f"{label:34s} {np.linalg.cond(A):9.1e} {L.dNd:10.2e} {ref['chi2']:12.5f} {sol['chi2'] - ref['chi2']:+11.1e} "
          f"{chi_res - ref['chi2']:+11.1e} {np.max(np.abs(g - gfd)) / np.max(np.abs(g)):13.1e}")
print("(code - QR: error of the value the optimizer sees; a profile difference of 0.5 is one sigma.)")

print("\n== (c) singular A ==")
p = Problem(n=24, nreg=2, dense=(), seed=2)
d = p.mean() + p.noise(np.random.default_rng(1))
bands = p.banddata(d)
for b in bands:
    b.ninv[5] = 0.0                                    # element 5 unobserved in every band (reviewer's own arrays)
for kw in ({}, {"cmb_prior_sd": 0.3}):
    try:
        v = RegionLikelihood(bands, p.region, **kw).profile(p.beta_true)
        print(f"unobserved element, {kw}: profile returns {v}")
    except Exception as e:  # noqa: BLE001
        print(f"unobserved element, {kw}: raises {type(e).__name__}: {str(e)[:70]}")
bands = p.banddata(d)
for b in bands[2:]:
    b.ninv[5] = 0.0                                    # element 5 observed in two bands only, three free amplitudes
try:
    v = RegionLikelihood(bands, p.region).profile(p.beta_true)
    print("element seen in 2 bands, 3 flat amplitudes: profile returns", v)
except Exception as e:  # noqa: BLE001
    print(f"element seen in 2 bands, 3 flat amplitudes: raises {type(e).__name__}: {str(e)[:70]}")
bands = p.banddata(d)
for b in bands[2:]:
    b.ninv[5] = 1e-14                                  # nearly unobserved: no exception expected
L = RegionLikelihood(bands, p.region)
A = L.assemble(p.beta_true)[0]
print(f"element nearly unobserved in 4 of 6 bands (weight 1e-14): cond(A) {np.linalg.cond(A):.1e}, profile {L.profile(p.beta_true):.6f}, "
      f"fitted |a_s| at that element {abs(L.solve(p.beta_true)['x'][5]):.2e}; no warning is given")

print("\n== (d) information() with a synchrotron prior ==")
p = Problem(n=36, nreg=3, dense=(0, 2), seed=3)
d = p.mean() + p.noise(np.random.default_rng(2))
lt = np.log(np.array([8.0, 15.0, 4.0]))
L = RegionLikelihood(p.banddata(d), p.region, cmb_prior_sd=0.3, sync_prior=True)
Jc, Tc, sol = L.information(beta, lt)
info = ref_info(p, beta, sol["x"][:p.n], cmb_sd=0.3, log_tau=lt)
print("max rel diff J", np.max(np.abs(Jc - info["J"])) / np.max(np.abs(info["J"])), " T", np.max(np.abs(Tc - info["T"])) / np.max(np.abs(info["T"])))

print("\n== (e) inverse covariance not exactly symmetric ==")
p = Problem(n=36, nreg=3, dense=(0, 2), seed=3)
bands = p.banddata(d)
L0 = RegionLikelihood(bands, p.region)
for eps in (1e-8, 1e-4, 1e-2):
    K = np.random.default_rng(0).standard_normal((p.n, p.n))
    K = eps * (K - K.T) * np.abs(bands[0].ninv).max()
    b2 = [BandData(b.name, b.d, (b.ninv + K if i == 0 else b.ninv), b.sync, b.dust) for i, b in enumerate(bands)]
    L1 = RegionLikelihood(b2, p.region)
    print(f"antisymmetric part {eps:.0e} of max|N^-1| in band 0: profile changes by {L1.profile(beta) - L0.profile(beta):+.3e} "
          f"(a symmetric quadratic form would not change); no error raised")

print("\n== (f) region labels with a gap ==")
reg = p.region.copy()
reg[reg == 1] = 0                                    # region 1 is now empty, labels are {0, 2}
L = RegionLikelihood(bands, reg)
bh, res = fit_profile(L, np.array([-3.0, -2.2, -3.0]))
print("fit with an empty region: beta_hat", bh, "(the empty region stays at its start value, success =", res.success, ")")
try:
    out = sandwich(L, bh)
    print("sandwich: returns sd", np.sqrt(np.diag(out[1])))
except Exception as e:  # noqa: BLE001
    print(f"sandwich: raises {type(e).__name__}: {e}")
cc, cs, Jp, Th, _, ok = sandwich_plugin(L, bh)
print("sandwich_plugin: ok", ok, " sd", np.sqrt(np.diag(cs)))

print("\n== (g) sandwich_plugin at very low signal-to-noise ==")
p = Problem(n=96, nreg=3, dense=(0, 2), seed=5, amp=0.25, noise_slope=2.5)
info = ref_info(p, p.beta_true, p.a_s, use_dust=False, use_cmb=False)
print("T/J at the truth", np.diag(info["T"]) / np.diag(info["J"]))
rng = np.random.default_rng(3)
zero = np.zeros(p.n)
nok, sds, jc_rel = [], [], []
for _ in range(400):
    L = RegionLikelihood(p.banddata(p.mean(a_d=zero, a_c=zero) + p.noise(rng)), p.region, use_dust=False, use_cmb=False)
    bh, _ = fit_profile(L, p.beta_true)
    cc, cs, Jp, Th, _, ok = sandwich_plugin(L, bh)
    nok.append(ok)
    sds.append(np.sqrt(np.diag(cs)))
    jc_rel.append(np.diag(Jp - Th) / np.diag(info["J"]))
nok, sds, jc_rel = np.array(nok), np.array(sds), np.array(jc_rel)
Ji = np.linalg.inv(info["J"])
sd_true = np.sqrt(np.diag(Ji @ (info["J"] + info["T"]) @ Ji))
print("fraction ok per region", nok.mean(0), "; J_c / J: mean", jc_rel.mean(0), "sd", jc_rel.std(0))
fin = np.where(nok, sds, np.nan)
print("plug-in sd / truth sandwich sd among ok: 5th, 50th, 95th percentile per region")
print(np.nanpercentile(fin / sd_true, [5, 50, 95], axis=0))
print("fraction of ok cases whose plug-in sd is above 3 times, and below one third of, the truth sandwich sd:",
      np.nanmean((fin / sd_true)[nok] > 3), np.nanmean((fin / sd_true)[nok] < 1 / 3))
