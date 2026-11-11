"""Referee re-derivations at toy size (n = 96 elements, 6 bands; largest matrix 576 x 576; one process; under a minute).

1. Masking: Schur complement of the inverse covariance equals the inverse of the kept block; prior restriction.
2. The reference-frequency identity: ln L_flat(beta; 30) - ln L_flat(beta; 23) = sum_r n_r beta_r ln(30/23), with a CMB prior.
3. Degrees of freedom: E[chi2_min] with the prior penalty inside = bands*n - 2n - regions when the CMB is drawn from
   the prior; flat CMB prior: bands*n - 3n - regions; the difference of the two nested minima has n.
4. Eq. (mismatch): my own dense formulas against the paper's and against 1500 simulated fits.
   My derivation: Sigma_w = N + U S_w U^T, P = Sigma_w^-1 - Sigma_w^-1 M (M^T Sigma_w^-1 M)^-1 M^T Sigma_w^-1 (M: synchrotron
   and dust columns), g_r = (1/2) d^T P_r d = -(M_r a)^T P eps + (1/2) eps^T P_r eps, eps ~ N(0, Sigma_t). Hence
     E[g_r]      = (1/2) tr(P_r U dS U^T)
     Cov(g_r,g_s)= (M_r a)^T P Sigma_t P (M_s a) + (1/2) tr(P_r Sigma_t P_s Sigma_t)
                 = J + q^T dS q   +   T + tr(P_r Sigma_w P_s U dS U^T) + (1/2) tr(P_r U dS U^T P_s U dS U^T)
     E[Hess]     = (1/2) a^T M^T P_rs M a + (1/2) tr(P_rs Sigma_t) = J + (1/2) tr(P_rs U dS U^T)
   The paper keeps J + T + q^T dS q and J: it drops two terms that are first order in dS (and of the order of T).
Never uses cholesky with lower=True on C-ordered arrays.
"""
import os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "lowres", "tests"))
from lowres.fit import fit_profile, fit_marginal, sandwich_plugin
from lowres.likelihood import BandData, RegionLikelihood
from lowres.mask import marginalize, restrict_prior
from lowres.mismatch import cmb_prior_mismatch
from test_block_prior_lowpass import setup, draw

t0 = time.time()
rng = np.random.default_rng(2026)
# ---- 1. masking
n = 60
X = rng.standard_normal((n, n)); N = X @ X.T / n + np.eye(n)
keep = np.ones(n, bool); keep[rng.choice(n, 9, replace=False)] = False
_, nk = marginalize(np.zeros(n), np.linalg.inv(N), keep)
print("1. Schur complement vs inverse of kept block: max abs diff", np.abs(nk - np.linalg.inv(N[np.ix_(keep, keep)])).max(),
      "; restrict_prior:", np.abs(restrict_prior(N, keep) - np.linalg.inv(N[np.ix_(keep, keep)])).max())
# chi-square of kept data under the marginal equals the full chi-square minimized over the masked data
d = rng.standard_normal(n); Ni = np.linalg.inv(N); m = ~keep
dm = -np.linalg.solve(Ni[np.ix_(m, m)], Ni[np.ix_(m, keep)] @ d[keep]); full = d.copy(); full[m] = dm
print("   kept chi2", d[keep] @ nk @ d[keep], " = full chi2 minimized over the masked elements", full @ Ni @ full)

# ---- 2. identity with a dense CMB prior, three regions, a dense band
n, nreg = 96, 3
region, a_s, a_d, specs = setup(n, nreg, seed=5, noise=1.0, amp=6.0, dense=(0, 2))
Xc = rng.standard_normal((n, 7)); S = 0.05 * Xc @ Xc.T / 7 + 0.03 * np.eye(n)
truth = np.array([-3.1, -3.0, -2.9])
a_c = np.linalg.cholesky(S) @ rng.standard_normal(n)
bands23 = draw(specs, a_s, a_d, a_c, truth[region], np.random.default_rng(1))
def with_ref(bands, nu):
    return [BandData(b.name, b.d, b.ninv, (lambda be, B=B: B.power_law(be, nu)), b.dust) for b, (B, _, _, _) in zip(bands, specs)]
for kw, lab in ((dict(cmb_prior_prec=np.linalg.inv(S)), "dense CMB prior"), (dict(), "flat CMB prior")):
    L23, L30 = RegionLikelihood(with_ref(bands23, 23.0), region, **kw), RegionLikelihood(with_ref(bands23, 30.0), region, **kw)
    nper = np.bincount(region, minlength=nreg)
    worst = 0.0
    b0 = np.array([-3.0, -3.2, -2.7])
    for _ in range(5):
        b1 = b0 + rng.uniform(-0.8, 0.8, nreg)
        lhs = -(L30.marginal(b1) - L30.marginal(b0)) + (L23.marginal(b1) - L23.marginal(b0))
        rhs = np.sum(nper * (b1 - b0)) * np.log(30 / 23)
        worst = max(worst, abs(lhs - rhs))
        assert abs(L30.profile(b1) - L23.profile(b1)) < 1e-7 * abs(L23.profile(b1))
    print(f"2. identity, {lab}: max |lhs - rhs| over 5 random pairs {worst:.2e} (rhs of order {abs(rhs):.1f}); profile objective independent of the reference")

# ---- 3. degrees of freedom
n, nreg = 96, 2
region, a_s, a_d, specs = setup(n, nreg, seed=13, noise=1.0, amp=6.0, dense=(0,))
truth = np.array([-3.1, -3.0])
Xc = rng.standard_normal((n, 9)); S = 0.3 * Xc @ Xc.T / 9 + 0.05 * np.eye(n)     # CMB prior with strong modes, comparable to the noise
root = np.linalg.cholesky(S); Sinv = np.linalg.inv(S)
nb = len(specs)
r3 = np.random.default_rng(77)
cl, cf = [], []
for _ in range(400):
    bd = draw(specs, a_s, a_d, root @ r3.standard_normal(n), truth[region], r3)
    Lp = RegionLikelihood(bd, region, cmb_prior_prec=Sinv); Lf = RegionLikelihood(bd, region)
    cl.append(2 * fit_profile(Lp, truth)[1].fun); cf.append(2 * fit_profile(Lf, truth)[1].fun)
cl, cf = np.array(cl), np.array(cf)
se = lambda x: x.std(ddof=1) / np.sqrt(x.size)
print(f"3. chi2_min with the prior penalty: mean {cl.mean():.2f} +- {se(cl):.2f}, expected bands*n - 2n - regions = {nb * n - 2 * n - nreg}; variance {cl.var(ddof=1):.0f} (2 dof = {2 * (nb * n - 2 * n - nreg)})")
print(f"   flat CMB prior: mean {cf.mean():.2f} +- {se(cf):.2f}, expected {nb * n - 3 * n - nreg}; difference of the nested minima: mean {(cl - cf).mean():.2f} +- {se(cl - cf):.2f}, expected n = {n}; "
      f"variance {(cl - cf).var(ddof=1):.0f} (2n = {2 * n}); correlation of the difference with the flat chi2: {np.corrcoef(cl - cf, cf)[0, 1]:+.3f}")

# ---- 4. mismatch (configuration of the repository's unit test, 1500 fits)
n, sd = 96, 0.3
region, a_s, a_d, specs = setup(n, 2, seed=13, noise=1.0, amp=6.0, dense=(0,))
r4 = np.random.default_rng(14)
X5 = r4.standard_normal((n, 5))
S_true = 0.6 * X5 @ X5.T / 5 + 0.02 * np.eye(n)
root = np.linalg.cholesky(S_true)
est, quoted = [], []
for i in range(1500):
    L = RegionLikelihood(draw(specs, a_s, a_d, root @ r4.standard_normal(n), truth[region], r4), region, cmb_prior_sd=sd)
    bh = fit_profile(L, truth)[0]
    est.append(bh)
    if i < 200:
        quoted.append(np.sqrt(np.diag(sandwich_plugin(L, bh)[1])))
est = np.array(est)
clean = [BandData(B.name, B.power_law(truth[region], 23.0)[0] * a_s + dust * a_d, ninv, (lambda b, B=B: B.power_law(b, 23.0)), dust) for B, dust, ninv, _ in specs]
Lc = RegionLikelihood(clean, region, cmb_prior_sd=sd)
dS = S_true - sd ** 2 * np.eye(n)
out = cmb_prior_mismatch(Lc, truth, 1.0 / sd ** 2, dS)
# the paper's formula at the truth, as the unit test evaluates it (J at true amplitudes; note J_c = J - T is used inside)
print("4. 1500 fits: scatter", est.std(0, ddof=1), "+-", est.std(0, ddof=1) / np.sqrt(2 * 1499), "; bias", est.mean(0) - truth, "+-", est.std(0, ddof=1) / np.sqrt(1500))
print("   lowres.mismatch (as in the unit test): sd_true", np.sqrt(np.diag(out["cov_true"])), "sd_matched", np.sqrt(np.diag(out["cov_matched"])), "bias", out["bias"])
print("   mean quoted (debiased sandwich) error of the white fit, 200 fits:", np.mean(quoted, axis=0))
# my own dense formulas
nb = len(specs)
def Mmat(beta):
    s = np.vstack([np.diag(B.power_law(beta[region], 23.0)[0]) for B, _, _, _ in specs])       # (nb n, n)
    dcol = np.vstack([dust * np.eye(n) for _, dust, _, _ in specs])
    return np.hstack([s, dcol])
U = np.vstack([np.eye(n)] * nb)
Nfull = np.zeros((nb * n, nb * n))
for i, (B, dust, ninv, rt) in enumerate(specs):
    Nfull[i * n:(i + 1) * n, i * n:(i + 1) * n] = np.linalg.inv(ninv) if ninv.ndim == 2 else np.diag(1.0 / ninv)
Sw = Nfull + sd ** 2 * U @ U.T
St = Nfull + U @ S_true @ U.T
Swi = np.linalg.inv(Sw)
def Pm(beta):
    M = Mmat(beta)
    return Swi - Swi @ M @ np.linalg.solve(M.T @ Swi @ M, M.T @ Swi)
h = 1e-4
E = np.eye(2)
P0 = Pm(truth)
Pr = [(Pm(truth + h * E[r]) - Pm(truth - h * E[r])) / (2 * h) for r in range(2)]
Prs = [[(Pm(truth + h * E[r] + h * E[s]) - Pm(truth + h * E[r] - h * E[s]) - Pm(truth - h * E[r] + h * E[s]) + Pm(truth - h * E[r] - h * E[s])) / (4 * h * h) for s in range(2)] for r in range(2)]
a = np.concatenate([a_s, a_d])
M0 = Mmat(truth)
Mr = [(Mmat(truth + h * E[r]) - Mmat(truth - h * E[r])) / (2 * h) @ a for r in range(2)]
dSig = U @ dS @ U.T
J = np.array([[Mr[r] @ P0 @ Mr[s] for s in range(2)] for r in range(2)])
T = np.array([[0.5 * np.trace(Pr[r] @ Sw @ Pr[s] @ Sw) for s in range(2)] for r in range(2)])
q = np.array([[Mr[r] @ P0 @ dSig @ P0 @ Mr[s] for s in range(2)] for r in range(2)])
Tt = np.array([[0.5 * np.trace(Pr[r] @ St @ Pr[s] @ St) for s in range(2)] for r in range(2)])
Eg = np.array([0.5 * np.trace(Pr[r] @ dSig) for r in range(2)])
Ht = np.array([[0.5 * (M0 @ a) @ Prs[r][s] @ (M0 @ a) + 0.5 * np.trace(Prs[r][s] @ St) for s in range(2)] for r in range(2)])
Hw = np.array([[0.5 * (M0 @ a) @ Prs[r][s] @ (M0 @ a) + 0.5 * np.trace(Prs[r][s] @ Sw) for s in range(2)] for r in range(2)])
Ji = np.linalg.inv(J)
print("   my J diag", np.diag(J), "T diag", np.diag(T), "T/J", np.diag(T) / np.diag(J), "; q^T dS q diag", np.diag(q), "; T_true - T diag", np.diag(Tt - T), "; E[H]_matched - J", np.diag(Hw - J), "; E[H]_true - J", np.diag(Ht - J))
print("   paper's formula with exact J, T:  sd", np.sqrt(np.diag(Ji @ (J + T + q) @ Ji)), " bias", -Ji @ Eg, " matched sd", np.sqrt(np.diag(Ji @ (J + T) @ Ji)))
Hi = np.linalg.inv(Ht)
print("   with the two dropped first-order terms (T under Sigma_true, expected curvature under Sigma_true): sd", np.sqrt(np.diag(Hi @ (J + q + Tt) @ Hi)), " bias", -Hi @ Eg)
print("   keeping T_true but J for the curvature: sd", np.sqrt(np.diag(Ji @ (J + q + Tt) @ Ji)))
print(f"   [{time.time() - t0:.0f} s]")
