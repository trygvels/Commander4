"""Addition 3: algebra of the low-pass construction of e042e_lowpass.py / e039c_lowpass_sim.py, and the validity of
`RegionLikelihood.information` with a singular N'^-1. Exact checks against the reviewer's reference in which the
low-passed bands enter as coefficient vectors c_b = Y_k^+ d_b with covariance Y_k^+ N_b Y_k^+T.
nside 2 for the exact reference (96 elements per band), nside 4 for the projector properties."""
import sys

import healpy as hp
import numpy as np

sys.path.insert(0, "nextgen")
sys.path.insert(0, "nextgen/reviews/stats_review_scripts_2026-11-11b")
from lowres.fit import sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.window import WindowMatcher  # noqa: E402
import rv2_common as rv  # noqa: E402

rng = np.random.default_rng(21)
print("# t3a (1): projector properties, nside 4, basis to lmax 11, kept multipoles 2 to 5")
WM = WindowMatcher(4, target_arcmin=60.0)
Y, Yp, ell = WM.Y, WM.Yp, WM.ell
n = Y.shape[0]
keep = ell <= 5
Yk, Ypk = Y[:, keep], Yp[keep]
Pout = np.eye(n) - Y @ Yp                                        # projector on the complement of span(Y)
print(f"   Y: {Y.shape}, kept columns {int(keep.sum())}; condition number of Y^T Y {np.linalg.cond(Y.T @ Y):.2f}")
print(f"   |Y+ Y - 1| {np.abs(Yp @ Y - np.eye(Y.shape[1])).max():.1e};  |Y_k+ Y_k - 1| {np.abs(Ypk @ Yk - np.eye(Yk.shape[1])).max():.1e};  "
      f"|Y_k+ Y_high| {np.abs(Ypk @ Y[:, ~keep]).max():.1e};  |Y_k+ (1 - Y Y+)| {np.abs(Ypk @ Pout).max():.1e}")
alt = np.linalg.pinv(Yk)                                          # the other candidate: least squares on the kept columns only
print(f"   alternative projector pinv(Y_k): |pinv(Y_k) Y_high| {np.abs(alt @ Y[:, ~keep]).max():.3f} (rms {np.sqrt(np.mean((alt @ Y[:, ~keep]) ** 2)):.3f}); "
      f"|pinv(Y_k) (1 - Y Y+)| {np.abs(alt @ Pout).max():.1e}")
# covariance and the script's N'^-1
A_ = rng.standard_normal((n, n))
Ncov = 0.04 * (np.eye(n) + 0.5 * (Y * (ell / 2.0) ** -2.0) @ Y.T) + 0.002 * A_ @ A_.T / n
ninv = np.linalg.inv(Ncov)
Cc = Ypk @ np.linalg.inv(ninv) @ Ypk.T                           # as in the scripts
nlp = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
nlp = 0.5 * (nlp + nlp.T)
w = np.linalg.eigvalsh(nlp)
print(f"   N'^-1: rank {int(np.sum(w > 1e-9 * w.max()))} of {n} (kept modes {Yk.shape[1]}); |N'^-1 N N'^-1 - N'^-1| / |N'^-1| {np.abs(nlp @ Ncov @ nlp - nlp).max() / np.abs(nlp).max():.1e}; "
      f"|N'^-1 Y_k Y_k+ - N'^-1| / |N'^-1| {np.abs(nlp @ Yk @ Ypk - nlp).max() / np.abs(nlp).max():.1e}")
d, m = rng.standard_normal(n), rng.standard_normal(n)
dlp = Yk @ (Ypk @ d)
c, cm = Ypk @ d, Ypk @ m
q_code = (dlp - m) @ nlp @ (dlp - m)
q_ref = (c - cm) @ np.linalg.solve(Cc, c - cm)
print(f"   quadratic form with (d', N'^-1) against the Gaussian form of c = Y_k+ d with Cov(c): {q_code:.8f} against {q_ref:.8f} (difference {q_code - q_ref:.1e})")
print(f"   N'^-1 annihilates: basis columns above the cut {np.abs(nlp @ Y[:, ~keep]).max() / np.abs(nlp).max():.1e}; the complement of span(Y) {np.abs(nlp @ Pout).max() / np.abs(nlp).max():.1e}")

print("# t3a (2): exact reference at nside 2 (96 elements, basis to lmax 5, 64 columns; kept multipoles 2 and 3, 24 columns)")
WM2 = WindowMatcher(2, target_arcmin=60.0)
Y, Yp, ell = WM2.Y, WM2.Yp, WM2.ell
n = Y.shape[0]
keep = ell <= 3
Yk, Ypk = Y[:, keep], Yp[keep]
npix = n // 2
lat = np.abs(hp.pix2ang(2, np.arange(npix), nest=True, lonlat=True)[1])
region = np.tile((lat < 30).astype(int), 2)
freqs = rv.FREQS
covs = [0.02 * (np.eye(n) + 0.5 * (Y * (ell / 2.0) ** -2.0) @ Y.T) + 0.002 * (lambda a: a @ a.T / n)(rng.standard_normal((n, n))) for _ in range(4)] + [np.diag(rng.uniform(0.01, 0.05, n))]
Ccmb = 0.14 ** 2 * (0.2 * np.eye(n) + 0.8 * (Y * (ell / 2.0) ** -3.0) @ Y.T)
Ccmb = 0.5 * (Ccmb + Ccmb.T)
beta0 = np.array([-3.1, -2.9])
a_s = rng.standard_normal(n) * np.array([2.0, 5.0])[region]
a_d = rng.standard_normal(n) * 3.0
maps = [None, Ypk, Ypk, Ypk, None]                                 # 23 GHz and 353 GHz keep all modes
ref = rv.Reference(freqs, covs, region, cmb_cov=Ccmb, maps=maps)


def code_like(dv, cmb_prior=True):
    bands = []
    for i, f in enumerate(freqs):
        if maps[i] is None:
            ni = np.linalg.inv(covs[i])
            bands.append(BandData(f"b{i}", dv[i], 0.5 * (ni + ni.T), rv.sync_fun(f), rv.dust_coef(f)))
        else:
            Cc = Ypk @ covs[i] @ Ypk.T
            nl = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
            bands.append(BandData(f"b{i}", Yk @ (Ypk @ dv[i]), 0.5 * (nl + nl.T), rv.sync_fun(f), rv.dust_coef(f)))
    kw = dict(cmb_prior_prec=np.linalg.inv(Ccmb)) if cmb_prior else {}
    return RegionLikelihood(bands, region, **kw)


def draw(rng, noise=True):
    a_c = np.linalg.cholesky(Ccmb) @ rng.standard_normal(n) if noise else np.zeros(n)
    return [rv.sync_fun(f)(beta0)[0][region] * a_s + rv.dust_coef(f) * a_d + a_c
            + (np.linalg.cholesky(covs[i]) @ rng.standard_normal(n) if noise else 0.0) for i, f in enumerate(freqs)]


dv = draw(rng)
L = code_like(dv)
betas = [beta0, beta0 + np.array([0.15, -0.2]), np.array([-2.6, -3.5])]
print("   profile, code - reference:", " ".join(f"{L.profile(b) - ref.profile(b, dv):+.2e}" for b in betas), f"(value {L.profile(beta0):.4f})")
mm = [L.marginal(b) - ref.marginal(b, dv) for b in betas]
print(f"   marginal, code - reference, spread over three beta (a constant is allowed): {np.ptp(mm):.2e}")
g = L.profile(betas[1], grad=True)[1]
print(f"   profile gradient against finite differences: {np.abs(g - rv.fd_grad(L.profile, betas[1])).max():.2e} (|g| {np.abs(g).max():.2e})")
L0 = code_like(draw(rng, noise=False))
Jc_, Tc_, sol0 = L0.information(beta0)
Jr, Tr = ref.info(beta0, a_s)
print(f"   noiseless data: |a_hat_s - a_s| {np.abs(sol0['x'][:n] - a_s).max():.1e}; J code / reference - 1 {np.abs(Jc_ / Jr - 1).max():.1e}; T code / reference - 1 {np.abs(Tc_ / Tr - 1).max():.1e}; "
      f"T/J {(np.diag(Tr) / np.diag(Jr)).round(4)}")
mean, cov = ref.score_moments(beta0, a_s, a_d)
print(f"   exact moments of the score of the coefficient likelihood: |mean| {np.abs(mean).max():.1e}; |Cov - (J + T)| / |J + T| {np.abs(cov - (Jr + Tr)).max() / np.abs(Jr + Tr).max():.1e}")
H = rv.fd_hess(lambda b: ref.expected_objective(b, beta0, a_s, a_d), beta0)
print(f"   Hessian of the expected objective at the truth against J: {np.abs(H - Jr).max() / np.abs(Jr).max():.1e}")
nmc = 3000
gs = np.zeros((nmc, 2))
for k in range(nmc):
    Lk = code_like(draw(rng))
    gs[k] = Lk._grad_chi2_half(Lk.solve(beta0))
print(f"   Monte Carlo of the code's score at the truth, {nmc} draws: mean / se " + " ".join(f"{x:+.2f}" for x in gs.mean(0) / (gs.std(0) / np.sqrt(nmc)))
      + "; variance / (J + T) " + " ".join(f"{x:.3f}" for x in gs.var(0) / np.diag(Jr + Tr)) + f" (se {np.sqrt(2 / nmc):.3f})")
# the information the high multipoles of the two full bands still carry about the index
ref_all = rv.Reference(freqs, covs, region, cmb_cov=Ccmb)
Ja, Ta = ref_all.info(beta0, a_s)
ref_2 = rv.Reference(freqs[[0, 4]], [covs[0], covs[4]], region, cmb_cov=Ccmb)
J2, T2 = ref_2.info(beta0, a_s)
print(f"   information J (diag): all scales {np.diag(Ja).round(1)}; low-pass arm {np.diag(Jr).round(1)}; the two full bands alone (23 and 353 GHz, all scales) {np.diag(J2).round(3)}"
      f" = {(np.diag(J2) / np.diag(Jr)).round(4)} of the low-pass arm")
# determinacy: without the CMB prior the normal matrix is singular (three amplitudes per element, two full bands)
A, *_ = code_like(dv, cmb_prior=False).assemble(beta0)
w = np.linalg.eigvalsh(A)
print(f"   no CMB prior: eigenvalues of A below 1e-10 of the largest: {int(np.sum(w < 1e-10 * w.max()))} of {A.shape[0]} "
      f"(n - min(regions, low-passed bands) x kept modes = {n - min(2, 3) * int(keep.sum())}: one free direction per element, less what the low-passed bands see); with the prior: smallest / largest {np.linalg.eigvalsh(L.assemble(beta0)[0]).min() / w.max():.1e}")
