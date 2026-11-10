"""Check 4: reference-frequency identity and the first-order shift of the flat-prior marginal optimum.

(a) profile(ref') = profile(ref); marginal(ref') - marginal(ref) = -sum_r n_r beta_r ln(ref'/ref), flat sync prior,
    with and without a CMB prior, dense noise, unequal region sizes.
(b) The same difference with the zero-mean synchrotron prior switched on (the identity is not claimed there).
(c) Shift of the flat-marginal optimum between ref and ref' against
      J^-1 n L            (the claim; J at the true amplitudes, and J at the fitted amplitudes)
      H_marg^-1 n L       (reviewer: H_marg = Hessian of the marginal objective at its optimum)
      (J + T + S)^-1 n L  (reviewer's expectation: S_rr' = tr(B M_rr') - tr(B M_r B M_r'))
    at three signal-to-noise levels.
"""
import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info, fd_hess, sync_ref
from lowres.fit import fit_marginal, fit_profile

np.set_printoptions(precision=5, linewidth=170, suppress=True)
REF, REF2 = 23.0, 30.0
LR = np.log(REF2 / REF)

print("== (a), (b) identity ==")
prob = Problem(n=60, nreg=3, dense=(0, 2, 5), seed=4)
d = prob.mean() + prob.noise(np.random.default_rng(8))
beta = np.array([-3.0, -3.25, -2.85])
nper = np.bincount(prob.region, minlength=prob.nreg)
print("region sizes", nper)
for name, kw in (("flat everything", {}), ("cmb prior 0.3", {"cmb_prior_sd": 0.3}), ("sync only", {"use_dust": False, "use_cmb": False})):
    L1 = RegionLikelihood(prob.banddata(d, nu_ref=REF), prob.region, **kw)
    L2 = RegionLikelihood(prob.banddata(d, nu_ref=REF2), prob.region, **kw)
    dm = L2.marginal(beta) - L1.marginal(beta)
    pred = -np.sum(nper * beta) * LR
    print(f"{name:18s} profile diff {L2.profile(beta) - L1.profile(beta):+.2e}   marginal diff {dm:+.10f}  claim {pred:+.10f}  "
          f"rel.err {abs(dm - pred) / abs(pred):.1e}")
lt = np.log(np.array([8.0, 15.0, 4.0]))
L1 = RegionLikelihood(prob.banddata(d, nu_ref=REF), prob.region, cmb_prior_sd=0.3, sync_prior=True)
L2 = RegionLikelihood(prob.banddata(d, nu_ref=REF2), prob.region, cmb_prior_sd=0.3, sync_prior=True)
print(f"sync prior on, same log_tau: marginal diff {L2.marginal(beta, lt) - L1.marginal(beta, lt):+.6f}  (flat-prior claim would be {-np.sum(nper * beta) * LR:+.6f})")
lt2 = lt + beta * LR        # tau scales as the amplitude does: tau' = tau (ref'/ref)^beta
print(f"sync prior on, log_tau' = log_tau + beta ln(ref'/ref): marginal diff {L2.marginal(beta, lt2) - L1.marginal(beta, lt):+.2e}  (exactly invariant if tau is rescaled)")

print("\n== (c) first-order shift of the flat-marginal optimum, ref 23 -> 30 GHz ==")
print("columns: region values of the shift beta_hat(30) - beta_hat(23)")
for (n, amp, slope, kw, rk, label) in [(320, 2.0, 2.5, dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False), "sync only, high S/N"),
                                       (320, 0.5, 2.5, dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False), "sync only, T/J~0.46"),
                                       (320, 0.35, 2.5, dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False), "sync only, T/J~0.93"),
                                       (384, 3.0, 2.5, dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "3 comp, cmb prior, high S/N"),
                                       (384, 0.8, 2.5, dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "3 comp, cmb prior, T/J~0.26"),
                                       (384, 0.6, 2.5, dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), "3 comp, cmb prior, T/J~0.46")]:
    p = Problem(n=n, nreg=2, dense=(0, 2, 5), seed=7, noise=1.0, amp=amp, amp_c=0.3, noise_slope=slope)
    sync_only = kw.get("use_dust") is False
    zero = np.zeros(n)
    mean = p.mean(a_d=zero, a_c=zero) if sync_only else p.mean()
    dd = mean + p.noise(np.random.default_rng(21))
    L1 = RegionLikelihood(p.banddata(dd, nu_ref=REF), p.region, **kw)
    L2 = RegionLikelihood(p.banddata(dd, nu_ref=REF2), p.region, **kw)
    b1, _, r1 = fit_marginal(L1, p.beta_true)
    b2, _, r2 = fit_marginal(L2, p.beta_true)
    bp, _ = fit_profile(L1, p.beta_true)
    npr = np.bincount(p.region, minlength=2).astype(float)
    info = ref_info(p, p.beta_true, p.a_s, **rk)
    J, T = info["J"], info["T"]
    Jh, Th, _ = L1.information(b1)
    Hm = fd_hess(L1.marginal, b1, h=2e-3)
    Hp = fd_hess(L1.profile, b1, h=2e-3)
    # S = tr(B M_rr') - tr(B M_r B M_r') from explicit matrices at the true beta
    B, dMs = info["B"], info["dM"]
    ncomp = 1 if sync_only else 3
    S = np.zeros((2, 2))
    for r in range(2):
        # d2M/dbeta_r^2: only the diagonal r = r' is nonzero
        M2 = np.zeros_like(dMs[r])
        for b, nu in enumerate(p.nus):
            s, ds = sync_ref(nu, p.beta_true, REF)
            M2[b * n:(b + 1) * n, :n] = np.diag((ds * np.log(nu / REF))[p.region] * (p.region == r))
        S[r, r] += np.trace(B @ M2)
        for q in range(2):
            S[r, q] -= np.trace(B @ dMs[r] @ B @ dMs[q])
    print(f"\n{label}:  n={n}  T/J {np.diag(T) / np.diag(J)}  profile opt {bp}  marginal opt(23) {b1}  marginal opt(30) {b2}")
    print("   actual shift                      ", b2 - b1)
    print("   claim  J^-1 n L (J at truth)      ", np.linalg.solve(J, npr * LR), "  ratio claim/actual", np.linalg.solve(J, npr * LR) / (b2 - b1))
    print("   claim  J(a_hat)^-1 n L (plug-in)  ", np.linalg.solve(Jh, npr * LR), "  ratio", np.linalg.solve(Jh, npr * LR) / (b2 - b1))
    print("   claim  (J(a_hat)-T)^-1 n L        ", np.linalg.solve(Jh - Th, npr * LR), "  ratio", np.linalg.solve(Jh - Th, npr * LR) / (b2 - b1))
    print("   H_marg^-1 n L (numerical Hessian) ", np.linalg.solve(Hm, npr * LR), "  ratio", np.linalg.solve(Hm, npr * LR) / (b2 - b1))
    print("   (J+T+S)^-1 n L (expected Hessian) ", np.linalg.solve(J + T + S, npr * LR), "  ratio", np.linalg.solve(J + T + S, npr * LR) / (b2 - b1))
    print("   diag: J", np.diag(J), " T", np.diag(T), " S", np.diag(S), " H_marg", np.diag(Hm), " H_profile", np.diag(Hp))
