"""Addition 1, claim (b): RegionLikelihood with a dense CMB prior precision against the reviewer's reference
(CMB integrated into a covariance that is dense across bands), finite differences, and exact score moments.
Small generic case: 40 elements, 5 bands, 3 regions, dense noise in three bands, random dense prior covariance."""
import sys

import numpy as np

sys.path.insert(0, "nextgen")
sys.path.insert(0, "nextgen/reviews/stats_review_scripts_2026-11-11b")
from lowres.fit import sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
import rv2_common as rv  # noqa: E402

rng = np.random.default_rng(11)
n, nreg = 40, 3
region = np.sort(rng.integers(0, nreg, n))
region[:nreg] = np.arange(nreg)
freqs = rv.FREQS
covs = []
for i, f in enumerate(freqs):
    covs.append(rv.random_spd(n, rng, cond=20.0, scale=0.05) if i in (0, 1, 3) else np.diag(rng.uniform(0.03, 0.2, n)))
Cc = rv.random_spd(n, rng, cond=300.0, scale=0.002)                           # prior covariance, far from white
beta0 = np.array([-3.1, -2.9, -3.3])
a_s = rng.standard_normal(n) * np.array([2.0, 4.0, 1.0])[region]
a_d = rng.standard_normal(n) * 3.0
ref = rv.Reference(freqs, covs, region, cmb_cov=Cc)


def draw(rng):
    a_c = np.linalg.cholesky(Cc) @ rng.standard_normal(n)
    out = []
    for i, f in enumerate(freqs):
        s, _ = rv.sync_fun(f)(beta0)
        out.append(s[region] * a_s + rv.dust_coef(f) * a_d + a_c + np.linalg.cholesky(covs[i]) @ rng.standard_normal(n))
    return out


def code_like(dv):
    bands = []
    for i, f in enumerate(freqs):
        ninv = np.linalg.inv(covs[i])
        ninv = 0.5 * (ninv + ninv.T) if i in (0, 1, 3) else np.diag(ninv).copy()
        bands.append(BandData(f"b{i}", dv[i], ninv, rv.sync_fun(f), rv.dust_coef(f)))
    return RegionLikelihood(bands, region, cmb_prior_prec=np.linalg.inv(Cc))


dv = draw(rng)
L = code_like(dv)
print("# t1c: dense CMB prior precision in RegionLikelihood against the reference (40 elements, 5 bands, 3 regions)")
betas = [beta0, beta0 + np.array([0.2, -0.3, 0.1]), np.array([-2.5, -3.6, -2.0])]
pc = np.array([L.profile(b) for b in betas])
pr = np.array([ref.profile(b, dv) for b in betas])
mc = np.array([L.marginal(b) for b in betas])
mr = np.array([ref.marginal(b, dv) for b in betas])
print("profile, code - reference at three beta:", " ".join(f"{x:+.2e}" for x in pc - pr), f"(values {pc[0]:.4f} {pc[1]:.4f} {pc[2]:.4f})")
print("marginal, code - reference (must be one constant: ln det of the CMB block of A):", " ".join(f"{x:+.10f}" for x in mc - mr),
      f"; spread {np.ptp(mc - mr):.2e}")
for b in betas[:2]:
    f, g = L.profile(b, grad=True)
    gm = L.marginal(b, grad=True)[1]
    print(f"gradient at {b}: profile max |analytic - finite difference| {np.abs(g - rv.fd_grad(L.profile, b)).max():.2e} (|g| {np.abs(g).max():.2e}); "
          f"marginal {np.abs(gm - rv.fd_grad(L.marginal, b)).max():.2e} (|g| {np.abs(gm).max():.2e})")

# information at the true amplitudes: code J and T evaluated at the fitted amplitudes of noiseless data
dv0 = []
for i, f in enumerate(freqs):
    s, _ = rv.sync_fun(f)(beta0)
    dv0.append(s[region] * a_s + rv.dust_coef(f) * a_d)                        # no CMB, no noise: a_hat_s = a_s exactly
L0 = code_like(dv0)
Jc_, Tc_, sol0 = L0.information(beta0)
Jr, Tr = ref.info(beta0, a_s)
print(f"noiseless data: max |a_hat_s - a_s| {np.abs(sol0['x'][:n] - a_s).max():.2e}; J code/reference - 1: {np.abs(Jc_ / Jr - 1).max():.2e}; "
      f"T code/reference - 1: {np.abs(Tc_ / Tr - 1).max():.2e}")
print("   J diag", np.diag(Jr).round(2), " T diag", np.diag(Tr).round(3), " T/J", (np.diag(Tr) / np.diag(Jr)).round(3))
# exact moments of the score under the model with the CMB drawn from the prior
mean, cov = ref.score_moments(beta0, a_s, a_d)
print(f"exact score moments (CMB drawn from the prior): |mean| max {np.abs(mean).max():.2e}; "
      f"max |Cov(g) - (J + T)| / max(J + T) {np.abs(cov - (Jr + Tr)).max() / np.abs(Jr + Tr).max():.2e}")
H = rv.fd_hess(lambda b: ref.expected_objective(b, beta0, a_s, a_d), beta0)
print(f"Hessian of the expected profile objective at the truth against J: max |H - J| / max J {np.abs(H - Jr).max() / np.abs(Jr).max():.2e}")
# plug-in: E[J(a_hat)] = J + T
M, Dm, At, Ati, Wt = ref.parts(beta0)
print(f"E[J(a_hat)] - J - T = tr(X Cov(a_hat_s)) - T with Cov(a_hat_s) = [At^-1]_ss: identity by construction; "
      f"check [A^-1]_ss of the code against [At^-1]_ss of the reference: {np.abs(L.solve(beta0, need_inverse=True)['Ainv'][:n, :n] - Ati[:n, :n]).max():.2e}")
# Monte Carlo of the code's score at the truth, as an end-to-end check of the exact moments
nmc = 4000
g = np.zeros((nmc, nreg))
jp = np.zeros((nmc, nreg))
for k in range(nmc):
    Lk = code_like(draw(rng))
    sol = Lk.solve(beta0)
    g[k] = Lk._grad_chi2_half(sol)
print(f"Monte Carlo of the code's score, {nmc} draws: mean / se " + " ".join(f"{x:+.2f}" for x in g.mean(0) / (g.std(0) / np.sqrt(nmc)))
      + "; variance / (J + T) diag " + " ".join(f"{x:.3f}" for x in g.var(0) / np.diag(Jr + Tr)) + f" (se {np.sqrt(2 / nmc):.3f})")
# what the score does when the CMB is NOT drawn from the prior (negative control): true covariance = 4 x prior
ref4 = rv.Reference(freqs, covs, region, cmb_cov=4 * Cc)
Mx, Dmx, _, _, _ = ref.parts(beta0)
G = Ati[:n] @ Mx.T @ ref.Si
m_ = Mx @ np.concatenate([a_s, a_d])
bias = []
for r in range(nreg):
    B = Wt @ (Dmx * (region == r)[None, :]) @ G
    bias.append(-(np.trace(B @ ref4.Sigma) + m_ @ B @ m_))
bias = np.array(bias)
print("negative control, true CMB covariance 4 times the prior: exact mean score / sd(score) " + " ".join(f"{x:+.2f}" for x in bias / np.sqrt(np.diag(Jr + Tr)))
      + "; first-order bias of the index in units of the sandwich sd " + " ".join(f"{x:+.2f}" for x in -np.linalg.solve(Jr, bias) / np.sqrt(np.diag(np.linalg.inv(Jr) @ (Jr + Tr) @ np.linalg.inv(Jr)))))
