"""Check 4b: the linear response of the flat-marginal optimum to the reference frequency, d beta_hat / d ln(nu_ref).

Check 4 showed that for 23 -> 30 GHz the shift is far outside the linear range at low signal-to-noise (n_r L >> J).
Here the first-order statement itself is tested: nu_ref is moved by exp(+-0.002) around a base value, so the shift is
tiny, and the measured slope is compared with
    J^-1 n               the claim, J at the true amplitudes
    H_marg^-1 n          exact first order: H_marg = Hessian of the marginal objective at its optimum (finite differences)
    (J + T + S)^-1 n     reviewer's expected Hessian: Hess (1/2) ln det A = T + S, S_rr' = tr(B M_rr') - tr(B M_r B M_r')
The base nu_ref is the pivot at which d (1/2) ln det A / d beta = 0 at the true beta, so that the marginal optimum is
near the truth and J at the truth is the relevant J. Averages over noise realizations, with standard errors.
Also: finite-difference Hessian of (marginal - profile) against T + S.
"""
import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info, fd_hess, sync_ref
from lowres.fit import fit_marginal

np.set_printoptions(precision=4, linewidth=170, suppress=True)
EPS = 0.002
NREAL = 40


def S_matrix(p, info, ref, nreg):
    B, dMs, n = info["B"], info["dM"], p.n
    S = np.zeros((nreg, nreg))
    for r in range(nreg):
        M2 = np.zeros_like(dMs[r])
        for b, nu in enumerate(p.nus):
            s, ds = sync_ref(nu, p.beta_true, ref)
            M2[b * n:(b + 1) * n, :n] = np.diag((ds * np.log(nu / ref))[p.region] * (p.region == r))
        S[r, r] += np.trace(B @ M2)
        for q in range(nreg):
            S[r, q] -= np.trace(B @ dMs[r] @ B @ dMs[q])
    return S


import sys
SO, RSO = dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False)
C3, RC3 = dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3)
if len(sys.argv) > 1 and sys.argv[1] == "lowsn":
    # amplitudes chosen so that T/J at the pivot is about 0.25, 0.5 and 1
    CASES = [(320, 0.27, SO, RSO, "sync only"), (320, 0.19, SO, RSO, "sync only"), (320, 0.136, SO, RSO, "sync only"),
             (384, 0.46, C3, RC3, "3 comp, cmb prior"), (384, 0.32, C3, RC3, "3 comp, cmb prior"), (384, 0.23, C3, RC3, "3 comp, cmb prior")]
else:
    CASES = [(320, 2.0, SO, RSO, "sync only"), (320, 0.5, SO, RSO, "sync only"), (320, 0.35, SO, RSO, "sync only"),
             (384, 3.0, C3, RC3, "3 comp, cmb prior"), (384, 0.8, C3, RC3, "3 comp, cmb prior"), (384, 0.6, C3, RC3, "3 comp, cmb prior")]
for (n, amp, kw, rk, label) in CASES:
    p = Problem(n=n, nreg=2, dense=(0, 2, 5), seed=7, noise=1.0, amp=amp, amp_c=0.3, noise_slope=2.5)
    sync_only = kw.get("use_dust") is False
    zero = np.zeros(n)
    npr = np.bincount(p.region, minlength=2).astype(float)
    # pivot: tr(B M_r)(ref) = 0  <=>  ln(ref_piv / ref0) = tr(B M_r) / n_r ; use the mean over regions
    i0 = ref_info(p, p.beta_true, p.a_s, nu_ref=23.0, **rk)
    lp = np.mean([np.trace(i0["B"] @ i0["dM"][r]) / npr[r] for r in range(2)])
    ref = 23.0 * np.exp(lp)
    info = ref_info(p, p.beta_true, p.a_s, nu_ref=ref, **rk)
    J, T = info["J"], info["T"]
    S = S_matrix(p, info, ref, 2)
    mean = p.mean(a_d=zero, a_c=zero, nu_ref=ref) if sync_only else p.mean(nu_ref=ref)
    # note: p.mean uses amplitudes defined at `ref`, consistent with info above
    rng = np.random.default_rng(31)
    slopes, slopes_H, bias = [], [], []
    for it in range(NREAL):
        dd = mean + p.noise(rng)
        opt = {}
        for sgn in (-1, 0, 1):
            L = RegionLikelihood(p.banddata(dd, nu_ref=ref * np.exp(sgn * EPS)), p.region, **kw)
            b, _, res = fit_marginal(L, p.beta_true, tol=1e-14)
            opt[sgn] = (b, L)
        slopes.append((opt[1][0] - opt[-1][0]) / (2 * EPS))
        bias.append(opt[0][0] - p.beta_true)
        if it < 10:
            Hm = fd_hess(opt[0][1].marginal, opt[0][0], h=2e-3)
            slopes_H.append(np.linalg.solve(Hm, npr))
    slopes, slopes_H, bias = np.array(slopes), np.array(slopes_H), np.array(bias)
    # Hessian of (1/2) ln det A at the true beta: finite differences of (marginal - profile) on any data set
    L = opt[0][1]
    Hld = fd_hess(lambda b: L.marginal(b) - L.profile(b), p.beta_true, h=2e-3)
    print(f"\n{label}: n={n}, pivot nu_ref {ref:.2f} GHz, T/J {np.diag(T) / np.diag(J)}, (T+S)/J {np.diag(T + S) / np.diag(J)}")
    print("   Hess (1/2) ln det A, finite differences diag", np.diag(Hld), "  T + S diag", np.diag(T + S), "  rel diff", np.diag(Hld) / np.diag(T + S) - 1)
    print("   marginal optimum - truth at the pivot (mean)", bias.mean(0), "+-", bias.std(0, ddof=1) / np.sqrt(NREAL),
          "  fits at a bound:", int(np.sum(np.any(np.abs(np.abs(bias + p.beta_true + 3.1) - 1.5) < 1e-6, axis=1))))
    m, se = slopes.mean(0), slopes.std(0, ddof=1) / np.sqrt(NREAL)
    print("   measured d beta_hat / d ln nu_ref (mean)     ", m, "+-", se)
    print("   claim J^-1 n                                 ", np.linalg.solve(J, npr), "  claim / measured", np.linalg.solve(J, npr) / m)
    print("   (J + T + S)^-1 n                             ", np.linalg.solve(J + T + S, npr), "  ratio to measured", np.linalg.solve(J + T + S, npr) / m)
    print("   H_marg^-1 n, first 10 realizations: ratio to the measured slope of the same realization, mean and range",
          (slopes_H / slopes[:10]).mean(0), (slopes_H / slopes[:10]).min(), (slopes_H / slopes[:10]).max())
