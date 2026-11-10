# Statistical review of nextgen/lowres: region-index likelihood, fit helpers, spectra, LFI noise model

- **Written:** 2026-11-10 21:51 CET
- **Reviewer:** stats-reviewer subagent, independent of the author's derivations
- **Code reviewed:** working tree read 2026-11-10 21:02 to 21:05 CET; identical to HEAD 35f03b4 for the five files
  under review (md5: likelihood.py b6caee73, fit.py 3ce194f7, spectra.py c8a0a180, lfi_noise.py b62cf8bf,
  tests/test_likelihood.py 6f3ba88b). `data.py` was read once and changed afterwards (21:12); it is not reviewed.
- **Order of work:** derivations first (`nextgen/experiments/E039/stats_review_scripts/derivations_before_reading.txt`,
  time-stamped 21:06), then numerical checks, then E039, E040, E041 and the Paper A1 draft were read (section 6).
- **Scripts and outputs:** `nextgen/experiments/E039/stats_review_scripts/` (`check_NN_*.py`, `out_NN_*.txt`,
  reference implementation `rv_common.py`). One or two threads throughout.
- **Existing tests:** 6 passed in 105 s (`out_00_existing_tests.txt`).

## Verdict

1. `profile`, `marginal` (both synchrotron priors), all gradients, and all of `lfi_noise.py`: correct for the stated targets (brute force 1e-13, finite differences 6e-8).
2. `information`, `sandwich`, `sandwich_plugin`: correct with caveats. Exact with flat priors; with the CMB prior only over CMB skies drawn from that prior, otherwise biased (+0.7 sd measured at 192 elements when the true rms is twice sigma_c; grows as sqrt(n)).
3. Reference frequency: identity exact (7e-14). The first-order shift J^-1 n L is the high signal-to-noise limit only; the linear response is (J + T + S)^-1 n L, and J^-1 is high by 1.5, 2.0, 3.0 at T/J = 0.25, 0.5, 1.0.
4. Coverage: sandwich at the truth and debiased plug-in are nominal for T/J 0.27 to 0.95 (J^-1, J_c^-1 and the plain plug-in fail); the sandwich sd is 3 percent low; the plug-in is unreliable above T/J of 1.5 even with `ok` True.
5. `spectra.py`: consistent, `np.gradient` acceptable (at most 2e-4); `wmap_band` code and docstring disagree on detector weighting (0.034 in the index). **Overall: correct with caveats.**

Here L = ln(ref'/ref), and T + S is the Hessian of (1/2) ln det A (section 1, item 6).

## 1. Derivations

Notation. Stack the bands: d = M(beta) x + n, x = (a_s, a_d, a_c), M_b = [diag(s_b(beta_r(p))), f_b I, I],
W = N^-1 (block diagonal over bands), Lam = prior precision = blockdiag(Lam_s, 0, Lam_c) with
Lam_c = I / sigma_c^2 and Lam_s = diag(tau_r(p)^-2) when the synchrotron prior is on. A = M^T W M + Lam,
rhs = M^T W d, xhat = A^-1 rhs = B d with B = A^-1 M^T W, P = W - W M A^-1 M^T W. M_r = dM/dbeta_r: only the
synchrotron block, only elements of region r. n_r = number of elements in region r.

**Items 1 and 2.** The joint negative log density of (d, x) given beta, tau is

    Q(x) + sum_p ln tau_r(p) + const,   Q = (1/2)(d - Mx)^T W (d - Mx) + (1/2) x^T Lam x
                                          = (1/2) d^T P d + (1/2)(x - xhat)^T A (x - xhat).

So min_x Q = (1/2) d^T P d = (1/2)(d^T W d - rhs^T A^-1 rhs) = chi2_min / 2 (item 1), and
-ln int exp(-Q) dx = chi2_min / 2 + (1/2) ln det A up to (2 pi)^(dim/2) (item 2). With the synchrotron prior the
normalization (2 pi tau^2)^(-1/2) per element adds sum_p ln tau_r(p) = sum_r n_r ln tau_r. The CMB prior
normalization does not depend on beta or tau. Identity used for the brute-force check: with every prior proper,
ln det(N + M Sig M^T) = ln det N + ln det Sig + ln det A, so the code value equals -ln N(d; 0, N + M Sig M^T)
minus (1/2) ln det N, minus n ln sigma_c, minus the 2 pi terms.

**Item 3.** Envelope theorem: d(chi2_min/2)/dbeta_r = dQ/dbeta_r at xhat = -(M_r xhat)^T W (d - M xhat)
= -sum_b sum_{p in r} ds_b(p) a_s(p) [N_b^-1 r_b](p).
d((1/2) ln det A)/dbeta_r = (1/2) tr(A^-1 dA/dbeta_r) = tr(A^-1 M_r^T W M)
= sum_b sum_{p in r} ds_b(p) [N_b^-1 M_b (A^-1)_{.,s}]_pp. The prior precision does not depend on beta.
d/dlog tau_r = sum_{p in r} [ -a_s(p)^2 / tau_r^2 - (A^-1)_pp / tau_r^2 + 1 ].

**Item 4, flat priors (Lam = 0).** Then P M = 0, B M = I, P N P = P, B N B^T = A^-1, B N P = 0.
dP/dbeta_r = P_r = -P M_r B - B^T M_r^T P, so the score is g_r = (1/2) d^T P_r d = -xhat^T M_r^T P d.
At the true beta, P d = P n and xhat = a + B n:

    g_r = -(M_r a)^T P n - n^T B^T M_r^T P n.

Mean: -tr(M_r^T P N B^T) = -tr(M_r^T P M A^-1) = 0. Covariance: the linear part gives
(M_r a)^T P N P (M_r' a) = J_rr'. The quadratic part n^T G_r n with G_r = B^T M_r^T P has covariance
tr(G_r N G_r'^T N) + tr(G_r N G_r' N); the first is tr(M_r^T P M_r' B N B^T) = tr(A^-1 M_r^T P M_r') = T_rr', the
second contains P N B^T = P M A^-1 = 0. Linear and quadratic parts are uncorrelated (odd moment). So
Cov(g) = J + T. Second derivative: P_rr' = -(X + X^T) with
X = P_r' M_r B + P M_rr' B + P M_r (A^-1 M_r'^T P - B M_r' B); taking E[(1/2) d^T P_rr' d] the signal part gives J
and the trace parts give -T + T = 0. So E[Hessian of chi2_min/2] = J, and to first order
Cov(beta_hat) = J^-1 (J + T) J^-1. All three statements are exact for any fixed true amplitudes. The covariance
statement is first order: it needs the Hessian to be close to its mean.

**Item 4 with a Gaussian CMB prior.** Minimizing over a_c first shows that P, and the (s, d) block of A^-1,
are identical to those of a reduced model with flat (a_s, a_d) and noise N~ = N + C Lam_c^-1 C^T (C = the CMB
columns). The code's J and T (built from the full A, with Lam_c inside) are therefore exactly the J and T of
that reduced model. Consequences:

- If the CMB sky is random with covariance sigma_c^2 I and the expectation is over noise and CMB, every
  flat-prior statement holds with the code's J and T: unbiased score, covariance J + T, expected Hessian J.
- For one fixed CMB sky (expectation over noise only), with x the true amplitudes, P M = B^T Lam is not zero and

      E[g_r | x] = -x^T M_r^T B^T Lam x + x^T Lam A^-1 M_r^T B^T Lam x - t_r,    t_r = tr(Lam A^-1 M_r^T B^T).

  The first term is linear in a_c (zero mean over skies, it is the part of J + T that is CMB sample variance).
  The score covariance over noise alone is smaller than J + T, and B N B^T = A^-1 - A^-1 Lam A^-1, so
  E[J(a_hat)] is not J + T either.
- If the CMB is random with covariance kappa^2 sigma_c^2 I (prior variance wrong by kappa^2):

      E[g_r] = (kappa^2 - 1) t_r,      first-order bias of beta_hat = -(kappa^2 - 1) J^-1 t.

  kappa = 0 is "no CMB in the data". t_r is a sum over elements of one sign, so the bias in units of the sd
  grows as sqrt(n_r). t_r vanishes when the prior is much tighter or much wider than the noise on the CMB and
  is largest when they are comparable. It does not depend on the synchrotron amplitudes, so bias/sd falls with
  signal-to-noise. With the arrays that `information` already forms,
  t_r = sigma_c^-2 sum_{p in r} sum_k (C A^-1)[p, cmb k] (A^-1)[cmb k, p] (verified, `check_02f`).

**Item 5.** At the true beta with flat priors, a_hat_s = a_s + e_s, Cov(e) = A^-1, so
E[J(a_hat)]_rr' = J_rr' + tr(X_rr' (A^-1)_ss) = J + T exactly (X = D^T P D, D the stack of diag(ds_b)).
J_c = J(a_hat) - T is unbiased for J. J_c^-1 J(a_hat) J_c^-1 is a ratio estimator: consistent, not unbiased, and
noisy when sd(J_c) is comparable to J. Same conditions on the CMB prior as item 4.

**Item 6.** s_b(beta; ref') = s_b(beta; ref) (ref/ref')^beta, so M' = M K with K = blockdiag(diag(kappa_p), I, I),
kappa_p = (ref/ref')^beta_r(p). With a flat synchrotron prior A' = K A K (a CMB prior commutes with K) and
rhs' = K rhs: chi2_min is unchanged and (1/2) ln det A' = (1/2) ln det A - sum_r n_r beta_r ln(ref'/ref). This
matches the claim. With the zero-mean prior the objective is invariant only if tau_r is rescaled with the
amplitude, tau' = tau (ref'/ref)^beta.

First-order shift. The marginal gradient changes by -n L, L = ln(ref'/ref), so the optimum moves by
H_marg^-1 n L with H_marg the Hessian of the **marginal** objective, not of the profile. In expectation

    E[H_marg] = J + Hess (1/2) ln det A,    Hess_rr' (1/2) ln det A = T_rr' + S_rr',
    S_rr' = tr(B M_rr') - tr(B M_r B M_r').

For one component and diagonal noise S = T (both equal sum_p V_p), so E[H_marg] = J + 2T. The claimed
J^-1 n L is therefore the T/J -> 0 limit, and overstates the linear response by 1 + (T + S)/J. The same
correction applies to "flat optimum minus profile optimum = J^-1 g": the Newton step from the profile optimum
is (J + T + S)^-1 g.

**Item 7.** (r I + Y D Y^T)^-1 = (1/r)[I - Y (r D^-1 + G)^-1 Y^T], G = Y^T Y, and
det = r^n det(D / r) det(r D^-1 + G). For a 2 by 2 SPD block S: sqrt(S) = (S + s I)/t with s = sqrt(det S),
t = sqrt(tr S + 2 s); inverse square root = adj(sqrt(S)) / s. Both match the code.

## 2. Code check

Matches my derivation, line by line:

- `likelihood.py:61-77` constant blocks and CMB prior; `:96-110` assembly including the tau^-2 diagonal;
  `:114-117` Cholesky, log determinant, chi2_min; `:138-144` profile gradient; `:146-158` log-determinant
  gradient (the transpose in `b.ninv.T` is right); `:165-182` marginal value, prior normalization and tau
  gradient; `:185-202` J and T (T = sum over region pairs of (A^-1)_ss times X, elementwise).
- `fit.py:31-38`, `:53-72`: J^-1 (J + T) J^-1 and J_c^-1 J(a_hat) J_c^-1 as stated.
- `lfi_noise.py:44-51`, `:88-99`, `:106-122`: as derived.
- `spectra.py:25-35`, `:56-68`: g, modified blackbody in brightness temperature, band average, power law and its
  derivative.

Discrepancies and weak points (ranked in section 5):

| where | what |
|---|---|
| `likelihood.py:22-26` (docstring), `fit.py:44-47` | "The profile score is unbiased", "expectation J + T": stated without the condition that a prior-constrained component is a draw from its prior. |
| `fit.py:31-35` | "Exact when those amplitudes are the truth (noiseless data)": with a CMB prior and CMB in the noiseless vector the fitted a_s is not the truth (J off by -4 and +15 percent in `check_02` case C). Exact if the clean vector has no CMB, as in `e039_lowres_validation.py:92`. |
| `fit.py:55` | `ok = diag(Jc) > 0` is the only guard on the plug-in error. |
| `fit.py:6-10`, `:13-28` | no flag for a fit that ends on a bound; `res.success` is False for harmless line-search stops. |
| `likelihood.py:117` | chi2_min by subtraction of two large numbers. |
| `likelihood.py:38-49`, `:114` | no symmetry check on a dense N^-1; Cholesky reads the lower triangle, the right-hand side uses the full matrix. |
| `likelihood.py:109`, `:160-161`, `:185-189`, `fit.py:36`, `:53` | `profile`, and `information`/`sandwich`/`sandwich_plugin` without `log_tau`, raise TypeError when `sync_prior=True` (loud, not silent). |
| `likelihood.py:75` | `cmb_prior_sd=0` means a flat prior (truthiness test), not "no CMB". |
| `spectra.py:47-49` | comment says trapezoid weights; `np.gradient` gives the full first and last interval to the end points (trapezoid: half). |
| `spectra.py:12-14` | "by less than 1e-3 for the WMAP K to V bands": V band is 2.5e-3. |
| `spectra.py:84-85` vs `:92-93` | docstring: detectors weighted by area; code: each detector divided by its area, so equal weights. |
| `spectra.py:70-74` | 4001-point grid; not unique for a band that contains the minimum of nu^beta / g (286 GHz for beta = -3.1). |
| `spectra.py:21` | h/k is the CODATA 2014 value; 3e-7 from 2018. |

## 3. Numerical checks

All checks use my own stacked-matrix reference (`rv_common.py`: explicit M, N, priors; QR of the whitened
least-squares problem; spectra typed independently) against the code under review.

### 3.1 Targets and gradients (`check_01`, n = 48, 6 bands, 3 dense, 3 unequal shuffled regions)

| check | result |
|---|---|
| profile, marginal against QR reference, 6 prior configurations | relative difference at most 3.3e-14; amplitudes agree to 7e-13 |
| marginal, all priors proper, against -ln N(d; 0, N + M Sig M^T) with constants | 8438.39923968 both; difference 1.8e-12 |
| flat priors as the wide-prior limit (prior sd 1e4), difference between two beta | +3.044220 against +3.044222; +7.213284 against +7.213283; -3.918353 against -3.918354 |
| direct two-dimensional quadrature over the amplitudes (2 elements, dense noise), flat and zero-mean prior | difference 1.5e-14 |
| gradients against central differences, 7 configurations, beta and log tau | at most 5.7e-8 relative |
| `information` with the synchrotron prior against explicit matrices (`check_07` d) | 7e-15 |

### 3.2 Score, its covariance, expected Hessian, E[J(a_hat)] at the true beta (`check_02`, `02c`, `02d`, `02e`)

n = 40 in two regions, 6 bands, 3 dense. T/J = 0.46, 0.57 (case A) and 0.41, 0.50 (others). Two numbers per cell
are the two regions. sigma_c = 2.0, comparable to the noise (1 to 2.5 per band), the regime where a prior bites.

| case | realizations | score mean / score sd | score var / (J+T) | E[Hessian] / J | E[J(a_hat)] / (J+T) |
|---|---|---|---|---|---|
| A flat priors, fixed truth | 80 000 | 0.000 +- 0.004, 0.002 +- 0.004 | 0.998 +- 0.005, 0.995 +- 0.005 | 1.003 +- 0.004, 1.001 +- 0.004 | 0.9995 +- 0.0009, 1.0002 +- 0.0008 |
| B CMB prior, CMB redrawn from it | 400 000 | 0.000 +- 0.002, 0.000 +- 0.002 | 0.999 +- 0.002, 0.999 +- 0.002 | 0.994 +- 0.003, 1.001 +- 0.003 | 0.9999 +- 0.0004, 1.0004 +- 0.0003 |
| C CMB prior, one fixed sky, rms = sigma_c | 20 000 | +0.313 +- 0.006, -0.363 +- 0.007 | 0.794 +- 0.008, 0.896 +- 0.009 | 1.25, 0.79 | 0.938 +- 0.002, 1.065 +- 0.002 |
| D CMB prior, CMB redrawn with 2 sigma_c | 20 000 | -0.98 +- 0.01, -1.25 +- 0.01 | 1.84 +- 0.02, 1.85 +- 0.02 | -0.07, -0.32 | 1.104 +- 0.002, 1.109 +- 0.002 |
| E CMB prior, no CMB in the data | 20 000 | +0.321 +- 0.006, +0.415 +- 0.006 | 0.785 +- 0.008, 0.786 +- 0.008 | 1.36, 1.42 | 0.968 +- 0.002, 0.963 +- 0.001 |

Closed form for the conditional score mean: case C +0.304, -0.362; case E +0.326, +0.414 (score sd units).
Case D equals -3 times case E, as (kappa^2 - 1) t predicts. T from the code equals T from explicit matrices to
1e-13. Hessian by central differences of the analytic gradient on a quarter of the realizations.

Audit of my own Monte Carlo errors: two of about 40 standardized deviations in single runs were beyond 3 sigma
(case B, seeds 777 and 4242). A fourth run of 200 000 with batch means shows independent samples (2000 batches of
100: ratio of between-batch to within-batch error 0.99 to 1.02 +- 0.016) and agrees with the exact values. The
pooled rows above are the result; I read the two outliers as chance.

### 3.3 Size of the CMB-prior bias (`check_02b`, closed form; `check_03` F5 to F7, Monte Carlo)

Toy at roughly WMAP/LFI depth: 8 bands (23 to 70 GHz and 353), noise 1 to 2.5 per element in every band below
100 GHz, 0.5 at 353, sigma_c = 0.3, synchrotron amplitudes of median 3 and rms 7, T/J 0.05 to 0.12. First-order bias of
beta_hat in units of the sandwich sd:

| elements per region | no CMB in the data | true CMB rms 0.6 (twice sigma_c) | rms over skies of the conditional mean, matched fixed sky |
|---|---|---|---|
| 48 | -0.17 | +0.50 | 0.26 |
| 192 | -0.28 | +0.83 | 0.26 |
| 768 | -0.46 | +1.39 | 0.28 |

Dependence on sigma_c at 192 elements (no CMB in the data): -0.003 (0.03), -0.04 (0.1), -0.28 (0.3), -0.83 (1.0),
-0.30 (3.0), -0.03 (10). Monte Carlo confirmation with fits (384 elements in 2 regions, 800 realizations):

| configuration | predicted bias / sd | measured bias / sd | coverage 68.3 / 95.4, sandwich at truth | same, plug-in |
|---|---|---|---|---|
| F7 prior 0.3, truth 0.3 (control) | 0, 0 | -0.05 +- 0.04, -0.05 +- 0.04 | 0.674 +- 0.012 / 0.946 +- 0.006 | 0.673 / 0.947 |
| F5 prior 0.3, truth 0.6 | +0.74, +0.82 | +0.67 +- 0.04, +0.78 +- 0.04 | 0.534 +- 0.012 / 0.876 +- 0.008 | 0.526 / 0.848 |
| F6 prior 0.3, no CMB | -0.25, -0.27 | -0.34 +- 0.03, -0.31 +- 0.03 | 0.684 / 0.949 | 0.689 / 0.954 |

### 3.4 Coverage of the profile optimum (`check_03`, `03b`)

Sync only (S) or three components (F), dense noise in 3 of 6 bands, 2 regions. T/J of 0.2 to 1 with a usable sd
needs hundreds of elements per region, so these runs are larger than the n of 20 to 100 asked for; S4 is the
small one. "truth" = J and T from explicit matrices at the true beta and amplitudes (equal to
`sandwich(L_clean, truth)` when the clean vector has no prior-constrained component). Cells: coverage of
+-1 sd / +-2 sd, nominal 0.683 / 0.954. Errors: 0.007 / 0.004 (S1, S2, S4), 0.015 / 0.007 (S3), 0.012 / 0.006 (F).

| config | elements per region | T/J | sandwich sd | bias / sd | empirical sd / sandwich sd | sandwich at truth | J^-1 at truth | `sandwich_plugin` | J_c^-1 | plain plug-in |
|---|---|---|---|---|---|---|---|---|---|---|
| S1 | 160 | 0.46 | 0.21, 0.19 | -0.02, -0.05 (+-0.02) | 1.02, 1.00 (+-0.016) | 0.679 / 0.951 | 0.587 / 0.900 | 0.681 / 0.956 | 0.587 / 0.905 | 0.562 / 0.886 |
| S2 | 160 | 0.95, 0.92 | 0.35, 0.32 | -0.04, -0.07 (+-0.02) | 1.07, 1.03 (+-0.016) | 0.672 / 0.944 | 0.513 / 0.844 | 0.686 / 0.955 | 0.517 / 0.854 | 0.462 / 0.797 |
| S3 | 640 | 0.67, 0.63 | 0.11, 0.10 | -0.03, -0.02 (+-0.05) | 1.05, 1.04 (+-0.03) | 0.664 / 0.948 | 0.545 / 0.867 | 0.661 / 0.954 | 0.548 / 0.860 | 0.519 / 0.831 |
| S4 | 48 | 0.27, 0.51 | 0.30, 0.42 | -0.04, -0.01 (+-0.02) | 1.00, 1.07 (+-0.016) | 0.676 / 0.945 | 0.604 / 0.902 | 0.687 / 0.938 | 0.604 / 0.896 | 0.587 / 0.883 |
| F1 CMB prior 0.3, redrawn | 192 | 0.46 | 0.30, 0.28 | -0.09, +0.01 (+-0.04) | 1.03, 1.02 (+-0.025) | 0.682 / 0.951 | 0.588 / 0.896 | 0.674 / 0.946 | 0.594 / 0.894 | 0.568 / 0.876 |
| F3 CMB prior 1.0, redrawn | 192 | 0.48 | 0.37, 0.34 | -0.05, +0.02 (+-0.04) | 1.04, 1.03 (+-0.025) | 0.676 / 0.949 | 0.569 / 0.890 | 0.670 / 0.952 | 0.570 / 0.892 | 0.547 / 0.869 |
| F4 flat CMB, one fixed sky | 192 | 0.48 | 0.38, 0.36 | -0.04, -0.06 (+-0.04) | 1.02, 0.97 (+-0.025) | 0.693 / 0.956 | 0.600 / 0.902 | 0.699 / 0.963 | 0.604 / 0.906 | 0.570 / 0.889 |
| F2 CMB prior 1.0, one fixed sky of rms 1.0 | 192 | 0.48 | 0.37, 0.34 | -0.19, -0.21 (+-0.03) | 0.95, 0.92 (+-0.025) | 0.707 / 0.964 | 0.622 / 0.911 | 0.722 / 0.971 | 0.632 / 0.913 | 0.602 / 0.900 |

Reading:

- Sandwich at the truth and the debiased plug-in are calibrated; J^-1, J_c^-1 and the plain plug-in are not
  (the negative controls fail as they should). `ok` was True in every realization of these runs.
- Empirical sd over first-order sandwich sd, pooled over the S and F rows with the right model:
  1.029 +- 0.005. The first-order formula is slightly narrow; the excess grows with T/J (1.01 at 0.46, 1.05 at
  0.93 for 160 elements).
- F2 against F3: same likelihood, same prior. Over skies drawn from the prior (F3) all is nominal. For one
  fixed sky (F2) the optimum is displaced by 0.2 sd (6 sigma of the Monte Carlo) and the scatter is 5 to 8
  percent smaller than the sandwich says. With a flat CMB prior (F4) a fixed sky is unbiased.
- One-sided rates of the plug-in +-2 sd interval (nominal 0.0228 each side): truth below the interval (optimum
  too flat) 0.025 (S1), 0.030 (S2), 0.046 (S4), 0.036 (F1); truth above 0.020, 0.015, 0.016, 0.018. Errors
  0.0024 (S) and 0.0037 (F). The plug-in sd is anticorrelated with the error (correlation -0.4 to -0.8): a flat
  fluctuation comes with a smaller quoted error. Two-sided coverage is fine, one-sided is not at 48 elements.
- Mean J(a_hat) at beta_hat over (J + T) at the truth: 0.998 to 1.006 (S1 to S3), 1.02 to 1.05 (S4, F1).
- Start values (`check_08`): 450 fits from -3.1, -2.0 and -4.2 end at the same optimum (difference 0.000 sd).
  The 3 to 7 of 800 fits with `success=False` are line-search stops at the optimum (`check_09`: remaining
  Newton step 0.000 sd).

Very low signal-to-noise (`check_07` g, 32 elements per region, T/J = 1.6 to 2.9, 400 realizations): `ok` True in
96 to 97.5 percent; among those the plug-in sd over the true sd has 5th, 50th, 95th percentiles 0.34 to 0.47, 0.83 to 0.96, 3.5 to 4.2;
8 percent are above 3 and 3 percent below one third. J_c / J has mean 1.3 to 1.8 and sd 1.2 to 2.1.

### 3.5 Reference frequency (`check_04`, `04b`, `04c`)

Identity (n = 60, regions of 12, 24, 24 elements, 23 to 30 GHz, three prior configurations): profile changes by
at most 3.6e-12; marginal changes by +48.4642574297 against the claimed +48.4642574297 (relative 7.5e-14). With
the zero-mean prior and the same tau the change is -250.69, as expected; with tau rescaled it is 1e-13.

Linear response, measured as (optimum at ref e^0.002 minus optimum at ref e^-0.002) / 0.004 around the pivot
reference (where the volume gradient vanishes), mean over 40 noise realizations, two regions:

| model | T/J | (T+S)/J | claim J^-1 n over measured | (J+T+S)^-1 n over measured | error of measured |
|---|---|---|---|---|---|
| sync only, n = 320 | 0.005 | 0.009 | 1.007, 1.009 | 0.998, 1.000 | 0.2 % |
| | 0.074 | 0.148 | 1.136, 1.142 | 0.989, 0.997 | 1.1 % |
| | 0.15 | 0.30 | 1.284, 1.290 | 0.986, 0.996 | 1.4 % |
| | 0.25 | 0.50 | 1.483, 1.488 | 0.984, 0.994 | 1.6 % |
| | 0.51 | 1.02 | 1.991, 1.989 | 0.982, 0.993 | 1.8 % |
| | 1.00 | 1.98 | 2.952, 2.939 | 0.983, 0.994 | 1.8 % |
| 3 components, CMB prior 0.3, n = 384 | 0.006 | 0.004 | 1.000, 1.002 | 0.996, 0.999 | 0.4 % |
| | 0.082 | 0.051 | 1.038, 1.044 | 0.989, 0.995 | 1.5 % |
| | 0.146 | 0.090 | 1.077, 1.082 | 0.989, 0.994 | 1.9 % |
| | 0.25 | 0.153 | 1.144, 1.145 | 0.993, 0.995 | 2.2 % |
| | 0.51 | 0.317 | 1.334, 1.313 | 1.015, 1.001 | 2.5 % |
| | 0.99 | 0.61 | 1.704, 1.618 | 1.060, 1.009 | 3 % |

Per realization, H_marg^-1 n with the finite-difference Hessian of the marginal objective reproduces the measured
slope to 0.9999 to 1.0004. The finite-difference Hessian of (marginal minus profile) equals T + S to four
digits. In the three-component case S is negative (-0.38 T), so the correction is smaller than for one
component, but it has the same sign.

Finite shift 23 to 30 GHz on one data set (`check_04`): claim over actual 1.03, 1.00 (sync only, T/J = 0.03) and
0.99, 1.00 (three components, T/J = 0.02). At T/J of 0.26 and above the actual shift is 1.1 to 1.9 in the index
and the 23 GHz optimum sits on or near the lower bound -4.6: no first-order formula applies there.

### 3.6 LFI noise model (`check_05`)

| check | result |
|---|---|
| `neg2loglike` against z^T R^-1 z + ln det R, dense, random non-orthogonal Y, (n, modes) = (60, 20), (60, 59), (40, 40), 5 parameter sets | absolute difference at most 4e-13; 7.6e-8 on 2.4e4 at (ln r, ln A) = (-6, 6) |
| same with `spin2_basis` at nside 2 and 4 | at most 2.3e-13; max abs(G - I) 0.17 and 0.14, eigenvalues of G 0.60 to 1.36 |
| `sqrt_blocks` against `scipy.linalg.sqrtm`, 2000 blocks including nearly singular | 6e-14 relative; inverse times root minus identity 5e-14 |
| `lfi_ninv` against inv(W^1/2 R W^1/2) with dense square roots (loader replaced in memory) | 1.6e-15 relative |
| `whitened_half_difference` against dense | 5.6e-16 |
| mean gradient of `neg2loglike` at the truth, 4000 draws | z = -0.19, -0.55, +0.47; without the log determinant (negative control) -70 |

### 3.7 Spectra (`check_06`, run with `python -I`)

- Delta band: power law, its derivative and the modified blackbody equal f(nu0)/g(nu0) to 1e-8 to 1.4e-6; the
  whole difference is the h/k constant. CMB column 1 to 16 digits. Derivative against finite differences
  1e-9. `power_law` returns shape (nbeta,).
- `np.gradient` as quadrature. Interior weights equal trapezoid weights on any grid; the two end weights are
  twice the trapezoid ones. Smooth bandpass, jittered non-uniform grid (3:1 step spread), error of the
  band-averaged power law against a converged integral:

| bandpass | points | np.gradient | trapezoid |
|---|---|---|---|
| response 3e-4 of peak at the ends | 20, 60, 200 | +7.0e-4, -5.1e-6, +9.0e-6 | +6.3e-4, -2.0e-5, +7.9e-7 |
| truncated at 5 % of peak | 20, 60, 200 | +3.9e-4, +2.1e-4, +9.5e-5 | -4.0e-4, -1.4e-5, -2.2e-6 |
| truncated at 30 % of peak | 20, 60, 200 | +4.3e-3, -2.4e-4, -2.2e-4 | -1.7e-4, -4.4e-6, +1.8e-7 |

  Acceptable when the tabulated response goes to zero at both ends; first order in the step otherwise. On the
  real files (code weights against trapezoid, beta = -3.1): WMAP K to V at most 6e-7; LFI 30 2e-8; LFI 44
  -8.4e-5 (upper end at 1.6 % of peak); LFI 70 +1.8e-4 (ends at 9.4 and 7.9 % of peak), mean log frequency
  changes by -9.5e-5; HFI 353 4e-8 (dust 4e-7).
- `effective_frequency`: delta-band response at the returned frequency equals the band average to 2e-5 to
  1.3e-4 (grid step). Monotonic over every band in use. A flat band from 220 to 420 GHz has two roots.
- Convention claim in the module docstring, sum(w f / g) / sum(w) against sum(w f) / sum(w g), beta = -3.1:
  K -4.7e-4, Ka -6.4e-4, Q -9.3e-4, V -2.5e-3.
- `wmap_band` weighting. Unit area per detector (the code) against the raw sum (area weighting, the docstring):
  s(-3.1) differs by -1.9e-2 (K), -5e-4 (Ka), +2.4e-3 (Q), +1.5e-2 (V); regressed on ln nu this is 0.034 in the
  index. Detector areas in the files span 1.29 to 2.84 (K) and 2.16 to 6.38 (V).
- `hfi_band`: dropping samples below the cut and zeroing them in place give the same integrals to 1e-11. The
  response between 420 and 583 GHz that survives the 1e-5 cut carries 2.0 % (PSB) to 2.2 % of the dust integral.

### 3.8 Numerical fragility (`check_07`)

chi2_min as computed (`likelihood.py:117`) against the QR reference, and the residual form
sum_b r_b^T N_b^-1 r_b for comparison (n = 48):

| case | cond(A) | d^T N^-1 d | error of the code value | error of the residual form | gradient error (relative) |
|---|---|---|---|---|---|
| amplitudes as in the tests | 5e4 | 5e4 | 3e-12 | 0 | 1.5e-9 |
| amplitudes 1e3 | 5e4 | 9e7 | 3e-8 | 8e-11 | 4e-10 |
| amplitudes 1e5 | 5e4 | 9e11 | 5e-4 | 1e-7 | 7e-10 |
| amplitudes 1e7 | 5e4 | 9e15 | 0.11 | 0.012 | 1e-10 |
| amplitudes 1e3, 353 GHz noise 1e-4 | 2e11 | 2e11 | 8e-4 | 1e-10 | 8e-8 |
| amplitudes and dust 1e5, 353 GHz noise 1e-3 | 2e9 | 7e19 | 4.7e3 | 2e-5 | 1e-7 |

The gradient is computed from residuals and stays accurate; the value loses eps times d^T N^-1 d times a
conditioning factor. Safe by a wide margin in every configuration I ran. Not verified on the real matrices.

Other: an element with zero weight in every band, or seen in fewer bands than it has free amplitudes, raises
LinAlgError (loud). With weight 1e-14 in four of six bands cond(A) is 4e15 and nothing is flagged. An
antisymmetric part of 1e-4 of max abs(N^-1) in one band changes the profile by 1.6 with no error (the loaders in
`data.py` symmetrize, so this is latent). An empty region: `sandwich` raises "Singular matrix",
`sandwich_plugin` flags it.

## 4. Tests: what they do and do not establish

Good: closed forms for one component with diagonal noise (independent reference); finite-difference gradients;
Monte Carlo of the sandwich and the plug-in with a visible negative control in the plug-in test (plain plug-in
sd of z above 1.08); the reference-frequency identity to 1e-9.

Gaps:

- No independent value check of `marginal` with three components and dense noise (only gradients). Done here.
- The CMB prior is never exercised where it matters: truth rms 0.2, prior 0.3, noise 3 to 12 per band
  (`test_likelihood.py:16-36`, `:114`, `:135`); the CMB is fixed, not redrawn, and at that depth it is
  invisible. A test with sigma_c comparable to the noise, in the three variants B, C, D of section 3.2, would
  pin the condition down.
- `test_likelihood.py:128` (`Sw > Ji`) is true by construction (T is positive semi-definite): no control.
- `test_likelihood.py:129`: the bias tolerance has an additive 0.02 on top of 4 standard errors, about 0.4 sd of
  the estimator here; a bias of a third of an sd passes.
- Every fit starts at the truth (`:119`, `:146`). I found no start dependence, but the tests would not see one.
- Tolerances of 0.12 on sd ratios with 300 fits are 3 standard errors (fine). `rtol=0.05` on the mean of
  J(a_hat) has no stated Monte Carlo error.
- Not tested: `information` with the synchrotron prior; one-sided coverage; T/J above 0.5; the plug-in when
  `ok` is False; the linear response to the reference frequency; `lfi_noise.py` and `spectra.py` have no tests.

## 5. Defects, ranked

**Wrong result (a stated claim or a number can be wrong)**

1. First-order reference-frequency shift (claim 6; also E039 V4, E040 R2, Paper A1 `sec_flat.tex` "delta beta =
   J^-1 g" and the coupled form). The curvature must be that of the marginal, J + T + S in expectation. J^-1
   overstates the shift by 1 + (T+S)/J: 14 % at T/J = 0.074, 48 % at 0.25, a factor 2 at 0.5 and 3 at 1.0 for
   one component; 4 %, 14 %, 33 %, 66 % at T/J = 0.08, 0.25, 0.5, 1.0 with three components and a CMB prior.
   Failure scenario: V4 or R2 with a 25 % tolerance fails in any region with T/J above about 0.13 (one
   component; 0.4 in my three-component toy) for a reason that is not a bug. Where such a test passes at
   low signal-to-noise, nonlinearity may be offsetting the missing factor (a conjecture, not tested here; E041
   L4 reports 4 % at 23 GHz and 15 % at 30 GHz). No code change needed: T + S is the Hessian of
   `marginal(beta) - profile(beta)`.
2. Unbiased score, J + T and E[J(a_hat)] = J + T with a Gaussian CMB prior (claims 4 and 5;
   `likelihood.py:22-26`, `fit.py:44-47`; the E039 derivation uses P M = 0, which is false with a prior). True
   only over CMB skies drawn from N(0, sigma_c^2 I). Failure scenario: real CMB Q/U rms per element 0.6 with
   sigma_c = 0.3 gives +0.7 to +0.8 sd bias at 192 elements per region and 68 % coverage of 0.53 (measured);
   no CMB at that depth gives -0.3 sd; both scale as sqrt(n_r). E039 redraws the CMB with the prior sd, so E039
   cannot see this. The flat-CMB arm of E040 (R7) is the unbiased reference; (kappa^2 - 1) J^-1 t predicts the
   difference between arms and t_r is one trace from arrays already formed.
3. `wmap_band` (`spectra.py:84-93`): docstring and code disagree on detector weighting; 0.034 in the index
   between the two, two to three times the E039 errors of the bright regions. I did not determine which is right.

**Fragile**

4. `sandwich_plugin` guard (`fit.py:55`): `ok` is True in 96 % of cases at T/J 1.6 to 2.9 with 32 elements, where
   the sd is off by more than 3 either way in 11 %. A guard on the significance of J_c (for example
   T/J_c below 1, or J_c above 3 sd(J_c)) is needed before small faint regions are quoted.
5. First-order sandwich sd low by 3 % on average, 7 % at T/J near 0.9 with 160 elements. Plug-in +-2 sd
   intervals are one-sided: at 48 elements the truth is below the interval 4.6 % of the time (nominal 2.3 %).
6. `fit_profile`, `fit_marginal` (`fit.py:6-28`): a fit on a bound is returned without a flag (2 of 2000 at
   T/J 0.93, 7 of 2000 at 48 elements; the flat marginal at 23 GHz sat on -4.6 in two of my six cases), and
   `res.success` False does not mean a bad fit. Callers do not check either.
7. chi2_min by subtraction (`likelihood.py:117`): absolute error 8e-4 at cond(A) = 2e11, 4.7e3 at
   d^T N^-1 d = 7e19. The residual form costs one extra pass and is good to 2e-5 in the same cases.
8. No symmetry or conditioning check on the inputs (`likelihood.py:38-49`, `:114`). Nearly unobserved elements
   pass silently with cond(A) = 4e15.
9. `sandwich` on noiseless data (`fit.py:31-38`) returns J at shrunk amplitudes when the clean vector contains a
   prior-constrained component (-4 and +15 % in J in my case C).
10. `data.py:load_lfi` (changed after my read, not reviewed): the cache file name does not depend on the fitted
    noise parameters, so a refit of the E038 model would silently reuse the old matrix.

**Cosmetic**

11. `spectra.py:47`: the weights are not trapezoid at the two ends; at most 1.8e-4 on the real bandpasses (LFI 70).
12. `spectra.py:12-14`: the two bandpass conventions differ by 2.5e-3 at V band, not below 1e-3.
13. `likelihood.py:109`, `:160`: `profile`, `information`, `sandwich` unusable with `sync_prior=True` (TypeError).
14. `likelihood.py:75`: `cmb_prior_sd=0` silently means a flat prior.
15. `fit.py:62-71`: after dropping regions the covariance of the rest is the conditional one (others held
    fixed), and rows and columns of dropped regions are inf, including cross terms.
16. `spectra.py:70-74`: resolution 1e-4; not unique across the minimum of nu^beta / g. `spectra.py:21`: h/k 3e-7.
17. `lfi_noise.py:96`: general solve on a triangular factor (correct, slower than a triangular solve).

## 6. Disagreements with the records read afterwards

- **E039 record, derivation paragraph.** "P M = 0 ... so the mean of the score is zero" is stated for a model
  whose priors include N(0, 0.2^2) on the CMB. With that prior P M = B^T Lam. The result is right for the E039
  simulation because `simulate()` redraws the CMB from the prior for every realization. The record should say
  so. Paper A1 `app.tex` does state the condition ("provided the true amplitudes of that component are a draw
  from the prior. If they are not, the score is biased by a term proportional to the mismatch"); it gives no
  size. `sec_profile.tex` ("The amplitudes are fixed unknowns, no distribution is assumed for them") holds for
  the flat-prior components only.
- **E039 V4, E040 R2, `sec_flat.tex`.** See defect 1. In the E039 nside-8 output so far (T/J below 0.01) the
  shift of region 0 is +0.069 against n L / J = 0.069: consistent with both formulas, as it must be at that
  depth. The statement in the Table caption of `sec_flat.tex` that the prediction "degrades where the shift
  exceeds the posterior width" misses the second reason: even an infinitesimal shift is smaller than J^-1 g by
  1 + (T+S)/J. For one component the boxed shift becomes ln(ref/pivot) n / (J + 2T), which saturates at low
  signal-to-noise instead of diverging.
- **E040 amendment 2 and `sec_profile.tex`.** E[J(a_hat)] = J + T and the debiased sandwich: confirmed, with the
  CMB-prior condition and the limits in defects 4 and 5.
- No disagreement on items 1, 2, 3, 7.

## 7. Not checked

- Nothing was run on the real WMAP or Planck maps and covariances: cond(A), T/J, t_r and the plug-in guard for
  the real configuration are open. Only the bandpass files were read.
- Instrument conventions: which WMAP detector weighting is right, the LFI and HFI bandpass definitions and
  units, the 1e-5 cut of the HFI transmission (2 % of the dust integral sits above 420 GHz), whether the
  ring-half difference contains all the noise of the full-mission map (noise common to the two halves cancels),
  and whether R fitted on the half difference transfers to W of the full map.
- `spin2_basis` beyond G and the Woodbury algebra (E/B convention, pixel window: the basis is sampled at pixel
  centres while the data are pixel averages).
- `data.py` (loaders, `degrade`, region maps) and `window.py`.
- Statistical behaviour of the zero-mean-prior optimum (bias, coverage); only its value and gradients.
- Coverage with more than 2 regions, with strongly unequal noise patterns, or with a non-white CMB.
- Whether the profile objective can be multimodal on real data (none seen in 450 simulated fits).
- No MCMC or sampler is part of this module; detailed balance, adaptation and RNG-stream checks do not apply.
