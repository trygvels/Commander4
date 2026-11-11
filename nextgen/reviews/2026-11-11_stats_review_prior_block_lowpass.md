# Statistical review: LCDM CMB prior, block-robust errors, low-pass estimator (nextgen/lowres and E039, E040, E042 scripts)

- **Written:** 2026-11-11 06:25 CET
- **Reviewer:** stats-reviewer subagent, independent of the author's derivations (only the claimed targets were given).
- **Code reviewed:** HEAD 2c92919, files unmodified in the working tree when read (05:41 CET). md5: `lowres/window.py` 509e30b7,
  `lowres/likelihood.py` d3ae3381, `lowres/fit.py` b89388b2, `lowres/lfi_noise.py` b62cf8bf,
  `scripts/e040f_block_errors_v2.py` 7b8fad6d, `scripts/e042e_lowpass.py` 66754e73, `scripts/e039c_lowpass_sim.py` 664aa2cd,
  `scripts/e039e_lowpass_sim_n16.py` 17b9ce57. The paper text (`papers/paperA1/sec_data.tex`) was read at 05:35 and has been
  edited by others since; quotes refer to that reading.
- **Order of work:** targets derived from the stated model before the code was opened, but written down only here (no
  time-stamped derivation file this time). Then code, then numerical checks, then the E040, E042 records and the paper.
- **Scripts and outputs:** `nextgen/reviews/stats_review_scripts_2026-11-11b/` (`tNx_*.py`, `out_*.txt`, reference
  implementation `rv2_common.py`, toy `toy_n4.py`). Working grids are nside 4 or smaller, or generic 96 to 192 element
  layouts; nside 32 appears only as a quadrature grid inside nside-4 pixels. Two threads, one process at a time, longest
  test 5.5 minutes, under 0.5 GB. Nothing under `lowres`, `scripts` or `papers` was edited.
- **Existing tests:** `nextgen/lowres/tests/test_likelihood.py`, 6 passed in 86 s (the LAPACK canary deselected by me to
  stay small). None of them exercises any of the three additions.

## Verdicts

| Claim | Verdict | One line |
|---|---|---|
| 1a prior covariance | correct with conditions | Exact for a field band-limited at 3 nside - 1 (1e-14 against an analytic two-point function). It is not the covariance of the pixel-averaged sky: at nside 16 that has 44 percent more variance, nearly white, 126 times the floor. |
| 1b objectives, gradients, J, T, sandwich with a dense prior | correct | Brute force 3e-13, finite differences 5e-8, exact score moments 1e-15. |
| 1c unbiased and calibrated when the CMB is drawn from the prior | correct | 500 realizations: mean z -0.01, +0.03, +0.05 (se 0.045), sd z 1.01, 0.97, 0.97 (se 0.032). |
| 2a sum of element scores is the gradient, zero at the optimum | correct | 1e-8 to 7e-8 against finite differences; 1e-5 to 1e-8 of the largest block score at the optimum. |
| 2b sum of D_B is the profile information | correct | Equal to the finite-difference Hessian to 1e-15; over J_c 0.99 to 1.013. |
| 2c CR0 low, CR3 jackknife, rms of z | correct with conditions | Exact for independent blocks. Fails (conservatively in every case I found) when the covariances couple blocks. The t reference with nu = n_eff is rough. |
| 2d CR3 tracks the rms about the region mean | correct with conditions | Holds for offsets independent between blocks (ratio 1.05 to 1.14). 0.60 to 0.70 when offsets are correlated over two blocks. The estimand is a population mean. |
| 3 low-pass likelihood, unbiased, calibrated, J and T with singular N'^-1 | correct | Identities to 1e-13 to 1e-15; 400 realizations: mean z -0.05, -0.06, -0.02 (se 0.05), sd z 0.95, 1.00, 1.05 (se 0.035). |
| 3 "amplitudes stay determined because one band keeps all modes" | correct with conditions | Needs a proper CMB prior as well; with a flat one A is singular. |

**Overall: correct with caveats.** No error was found in the algebra or the code of any of the three additions. The
caveats are about what the block error and the prior are claimed to be.

## Defects

### Code defects (all low severity)

1. `scripts/e040f_block_errors_v2.py:51,61`: `ok` from `sandwich_plugin` is discarded and `Ji = np.linalg.inv(Jp - T)` is
   taken over all regions. A region with J_c <= 0 gets finite, meaningless block errors while the sandwich returns inf.
   Failure scenario: a faint region in a subset fit; the table then shows a block error for a region with no measurement.
2. `scripts/e040f_block_errors_v2.py:57-60`: only dG_{B,r}/dbeta_r is used, negative D_B are set to zero silently, and
   neither their number nor the cross-region share is printed. In a model-true case with strong large-scale noise modes
   every realization had a negative D_B, the cross-region share of |dG_{B,r}/dbeta_s| was 0.49, and the code's CR3 was
   0.94 and 0.87 of the full-Jacobian form (`out_t2b_modes_tau0.txt`). With diagonal or mildly dense covariances the
   share is 0.000 and 0.005 and the two forms agree to 1e-3. On the maps this is not monitored.
3. `scripts/e040f_block_errors_v2.py:59`: the clip at 0.9. No effect for h <= 0.9. Above it CR3 is capped and CR2 turns
   anti-conservative: with one block at h = 0.95 (0.99), E[CR2] = 0.52 (0.11) and E[CR3] = 4.8 (1.0) times the true variance,
   where unclipped CR3 gives 19 (99) (`out_t2a.txt`). A block error should not be returned when max h > 0.9. Not
   reached on the maps (largest h 0.48).
4. `scripts/e040f_block_errors_v2.py:9` (docstring): CR2 is unbiased when Var(G_B) = D_B (1 - h_B), that is when the
   information equality holds block by block and blocks are independent, not "for equal-variance scores".
5. `lowres/likelihood.py:228-232`: `t_cmb` (the sensitivity to a wrong prior amplitude) is returned only for the white
   prior; with `cmb_prior_prec` nothing is returned although the docstring formula (line 35) is general.
   `lowres/likelihood.py:101-103`: no check of shape, symmetry or positive definiteness of the precision, and a precision
   passed with `use_cmb=False` is ignored silently.
6. `scripts/e042e_lowpass.py:52`, `e039c_lowpass_sim.py:63-65`, `e039e_lowpass_sim_n16.py:64-65`: the construction is
   exact only if Y_k^+ Y_k = 1, which needs Y of full column rank. No assertion. Condition number of Y^T Y: 1.8 at nside 2,
   2.3 at nside 4; not checked by me at nside 8 or 16 (outside my size limit). One line to add:
   `assert np.abs(Ypk @ Yk - np.eye(Yk.shape[1])).max() < 1e-8`.
7. `scripts/e042e_lowpass.py:79`: prints "low-passed copy" for 353 GHz, for which none is made (cosmetic).
   `scripts/e042e_lowpass.py:91-93`: a fit with `converged False` (K, Ka, 353 at l <= 12 in `E042/lowpass12.txt`) is stored
   in the npz like the others.
8. `lowres/window.py:59`: CAMB is run to lmax 300, which gives the lensing BB 7 percent low at l <= 47 (EE within 1e-3).
   Irrelevant here (BB is 0.013 uK rms per element), stated for completeness.

### Claims worded too strongly

1. **"The LCDM covariance of E and B modes on the grid"** (`window.py:52-55`, paper Sect. data, model paragraph). It is
   the covariance of multipoles 2 to 3 nside - 1 through the pixel window. The maps are pixel averages of skies that are
   not band-limited. Variance per element of pixel-averaged LCDM Q or U, 60 arcmin window (`out_t1b.txt`):

   | nside | in band | above the band limit | share | times the floor variance (1e-4 uK^2) |
   |---|---|---|---|---|
   | 16 | 0.0288 uK^2 (0.170 uK) | 0.0126 uK^2 (0.112 uK) | 0.44 | 126 |
   | 8 | 0.0209 (0.145) | 0.0026 (0.051) | 0.13 | 26 |
   | 4 | 0.0190 (0.138) | 0.0004 (0.020) | 0.02 | 4 |

   The in-band values reproduce the 0.17 uK of the paper and the 0.145 uK of `E040/cmbprior_sim_n8.txt`. The part above the
   band limit is close to white on the grid (30 to 60 percent of it inside the span of the basis, the rest outside, where
   the prior allows 0.01 uK). So the pixel-averaged LCDM sky has 0.20 uK rms per element, not 0.17, and the floor is not a
   description of what lies outside the basis. Size of the effect on the index (exact first-order bias, `out_t1e.txt`):
   +0.07 to +0.18 quoted sd with errors 11 percent too small when the CMB rms is 0.6 of the noise per element (Monte Carlo,
   `out_t1d.txt` arm C: mean z +0.16, +0.16, +0.14 +- 0.045, sd z 1.10, 1.11, 1.06; repaired by adding the white part to
   the prior, arm D: mean z -0.01, +0.04, +0.05, sd z 0.99, 0.99, 0.96); +0.003 sd at 384 elements when the ratio is 0.08 as
   on the WMAP maps, so of order 0.01 sd at 6144 elements if it grows as sqrt(n). Negligible for this paper, not for
   low-noise data. Every simulation in the repository draws the CMB from the same truncated covariance, so none can see it.
   Fix: `floor_uk` of 0.11 at nside 16 (0.05 at nside 8) as a sensitivity arm, or the exact pixel-average covariance.
   Second, smaller: the isotropic pixel window is 10 to 50 percent off in variance for elements within 24 degrees of a
   pole at nside 4 (48 of 384 elements, `out_t1f.txt`), because Q and U are averaged across a rotating frame there; the
   mean over the sky is right to 3 percent.
2. **"The error of 'the index of this region' is three to eight times the error of the estimator under the model"**
   (paper, block-error paragraph), and claim 2d. The block error is the error of the estimate as a measurement of the
   mean of a population of exchangeable, independent blocks. Three qualifications, each measured:
   - The estimate stays within 1.05 sandwich sd of the realized information-weighted mean index of the region for a
     spread of 0.1 (rms 0.0207, 0.0297 against sandwich 0.0196, 0.0284; `out_t2b_white_tau0.1.txt`). For that estimand the
     sandwich error is the right one.
   - If the index offsets are correlated between neighbouring blocks (correlation length two blocks), CR3 is 0.60 and
     0.70 of the true rms about the region mean and +-1 CR3 sd covers 0.44 and 0.45 (`out_t2b_white_tau0.1_corr2.txt`).
     This is the anti-conservative direction, and a smooth index field on 30 degree blocks is the plausible case. The
     simulation in the script (line 98) draws independent offsets only.
   - The variance estimate has fewer degrees of freedom than n_eff once the spread dominates: Satterthwaite 1.8 to 2.2
     at n_eff 4.3, 3.6 to 4.7 at n_eff 8, 4.1 to 5.3 at n_eff 9 (`out_t2a.txt`); +-2 CR3 sd covered 0.86 to 0.93.
3. **"Under the model that ratio has a median near one"** and the use of the map ratios as evidence that the regions are
   not one index each. True for independent blocks. With a sky coherent over many blocks and eight large-scale noise modes
   of 50 times the white variance (the generalized least squares fit gives them no weight), model true, the block errors
   are 1.44 to 1.51 (CR0), 1.56 to 1.75 (CR2) and 1.74 to 2.13 (CR3) times the true scatter, while the sandwich is
   calibrated (rms z 1.00) (`out_t2b_modes_tau0.txt`). Exact reason: the block scores are negatively correlated, sum_B
   Var(G_B) / Var(sum_B G_B) = 1.20 and 1.68 (`out_t2b_exact_block_covariance.txt`). The recorded nside-8 simulation with
   the WMAP matrices shows the same sign at High Latitude (CR3 rms z 0.82 with the model true, 0.59 with a spread of 0.1,
   where CR3 is 0.058 against an rms of 0.037). The null ratio therefore has to come from simulations with the covariances
   and prior of the fit: nside 16, LCDM prior. That exists only at nside 8 with the white prior. The Plane and Galactic
   Center ratios (6 to 9) are safe against any null I produced; High Latitude and Spur (CR0 over sandwich 2.4 and 2.0 with
   the LCDM prior, 2.7 and 1.6 with the white prior) are 1.2 to 1.9 times my worst null median (1.4, with a 16 to 84
   percent spread of the kind recorded in E040) and are not safe until that simulation exists.
4. **"1.1 to 1.3 is expected for a variance estimate with 4 to 8 degrees of freedom"** (paper; E040 record "CR3 is the
   closest to that"). For independent blocks and the model true, CR3 gives rms z 1.07 to 1.18 at n_eff 4 to 10, not
   sqrt(n_eff / (n_eff - 2)) = 1.13 to 1.39, because CR3 is biased high by sum h / (1 - h) (1.14 at n_eff 8, 1.35 at 4.3)
   and that offsets the tails; +-1 CR3 sd covers 0.68 to 0.70 and +-2 sd 0.92 to 0.94. The sqrt(nu / (nu - 2)) reference
   at nu = n_eff describes CR2 (rms z 1.31 against 1.37 at n_eff 4.3, 1.14 against 1.15 at 8.0), not CR3. With 60 realizations the rms of z is known to 0.1 to 0.2, so the recorded 0.82,
   1.20, 1.02 cannot tell these apart; "reproduces the scatter" holds at that precision only.
5. **"CR3 is the delete-one-block jackknife to first order."** True for deleting the block's score terms. It equals
   deleting the block's data only when the covariances do not couple blocks: CR3 over the true delete-data jackknife
   (refit without the block) is 1.000 (16 to 84 percent: 0.99 to 1.01) with diagonal noise and a white prior, 0.97 and
   0.93 (0.82 to 1.07) with noise correlated over two blocks and a dense prior, 1.5 to 1.7 with the strong modes.
6. **"The per-pixel amplitudes stay determined because one band keeps all modes"** (claim 3, `e042e_lowpass.py:9-10`). Two
   full bands and three amplitudes per element: the CMB prior is what closes it. Without it A has 48 null directions of
   288 at nside 2 with two regions (`out_t3a.txt`). All three scripts use a proper prior, so nothing is wrong in the runs.
   "The pixel scale does not inform them" is exact for the low-passed bands and approximate for the anchor (below).

## 1. CMB prior from the LCDM spectra

**Derivation.** healpy's `alm2map` is real-linear in (Re a_lm, Im a_lm). For an isotropic Gaussian field Re and Im of
a_lm (m > 0) have variance C_l / 2 and a_l0 has C_l, so Cov = sum over real modes of C_l v v^T with v the map of a unit
real coefficient times 1 (m = 0) or 1/sqrt 2 (m > 0). `spin2_basis` returns v sqrt(Omega_pix) with exactly these factors,
hence Cov = Y diag(C_l / Omega) Y^T. Beam and pixel window multiply a_lm, so they enter squared. With a prior precision
Lam on the CMB block, minimizing over a_c is integrating it out: chi2_min equals the generalized least squares chi2 of
(a_s, a_d) under Sigma = N + U Lam^-1 U^T (U stacks the CMB column; Sigma is dense across bands), and ln det A differs
from ln det(M_sd^T Sigma^-1 M_sd) by ln det(U^T N^-1 U + Lam), constant in beta. So all flat-prior results carry over
with N replaced by Sigma: score mean 0, Cov(score) = J + T, E[Hessian] = J, E[J(a_hat)] = J + T, with W, A^-1 computed
from the A that includes Lam. This holds over CMB skies drawn from Lam^-1 only. The code does exactly this.

**Tests.**
- `t1a` (analytic two-point function with Wigner d functions and explicit frame rotations, pixel positions from
  `pix2vec(nest=True)`, no harmonic transform): max |code - analytic| / max = 1.5e-15 (nside 2), 1.3e-14 (nside 4), all
  four blocks. Negative controls on the same statistic: E and B swapped 0.33 to 0.81, U sign flipped 0.57 to 0.64, RING
  order 0.60, missing sqrt 2 0.50. `lcdm_cmb_precision` minus its floor against the analytic form with the CAMB spectra:
  3e-15; controls: EE and BB columns swapped 0.50, window not squared 0.04 to 0.07, temperature pixel window 0.02.
  prec times cov is the identity to 4e-12; eigenvalues of cov from 1e-4 to 0.49.
- `t1c` (40 elements, 5 bands, 3 regions, dense noise in three bands, random dense prior with condition number 300):
  profile code minus reference 3e-13; marginal minus reference constant to 3e-13; gradients against finite differences
  5e-8; J and T against the reference 2e-14 and 1e-15; exact score mean 5e-15 and |Cov - (J + T)| 9e-16; Hessian of the
  expected objective against J 4e-7; Monte Carlo of the code's score (4000): mean within 1.1 se, variance / (J + T) 0.98
  to 0.99 (se 0.022). Negative control: true CMB covariance four times the prior gives a first-order bias of +0.22 to
  +0.84 sandwich sd.
- `t1d` (nside 4, 384 elements, 3 regions, dense red noise, 500 realizations, code's fit and plug-in sandwich):
  A, CMB from the returned covariance: mean z -0.013, +0.030, +0.046 (se 0.045); sd z 1.009, 0.970, 0.969 (se 0.032);
  coverage 0.664, 0.702, 0.730; median error 0.0565, 0.0289, 0.0122 against scatter 0.0569, 0.0280, 0.0118.
  B, white prior of the same rms on the same data: sd z 0.93, 0.86, 0.93 (errors too large here; E040 found them too
  small at nside 8: the sign depends on the configuration, the miscalibration does not). C and D as in item 1 above.

## 2. Block-robust errors

**Derivation.** Envelope theorem: d(chi2_min / 2)/dbeta_r = -sum_b sum_{p in r} (N_b^-1 r_b)_p ds_b(p) a_s(p), so
sum_p g_p is the gradient (2a), with or without a prior that does not depend on beta. Each g_p has mean zero at the truth
(M^T W = 0 as a matrix identity). sum_B D_B is the derivative of the gradient with the amplitudes re-solved, the observed
profile Hessian; E[H] = J(a) = E[J_c] (2b). For independent blocks and one index, with h_B = D_B / D and
Var(G_B) = D_B (1 - h_B): E[CR0] = (1 - sum h^2) / D, E[CR2] = 1 / D, E[CR3] = sum h / (1 - h) / D, and deleting block B
moves the estimate by (G_B / D) / (1 - h_B), so the sum of squared moves is CR3 with bread D = sum D_B, which is
G / (G - 1) times the jackknife variance. With a random offset of variance tau^2 per block the estimate has variance
tau^2 sum h^2 + 1 / D about the population mean. With dense N or a dense prior, Cov(G_B, G_B') = (D a)_B^T W_BB' (D a)_B'
is not zero and E[D_B] = Cov(G_B, sum G) is not Var(G_B): none of the three expectations above is then guaranteed.

**Tests.** `t2a` (abstract, 4e5 draws), `t2b` (the script's own `block_scores` and `block_errors` on a 192-element,
24-block, two-region likelihood, n_eff 10 and 4, one block straddling the boundary).

Model true, sqrt(mean variance) over true rms, regions 0 and 1 (`out_t2b_*_tau0.txt`; 2000, 2000, 1500 realizations):

| arm | sum_B Var(G_B) / Var(sum G), exact | sandwich | CR0 | CR2 | CR3 | CR3 rms z | CR3 cover +-1, +-2 | CR3 / delete-data jackknife |
|---|---|---|---|---|---|---|---|---|
| white (independent blocks) | 1.000, 1.000 | 1.00, 1.01 | 0.95, 0.89 | 1.00, 1.03 | 1.05, 1.20 | 1.07, 1.17 | 0.68, 0.70; 0.94, 0.92 | 1.000, 0.999 |
| dense noise and prior | 1.03, 0.97 | 1.01, 1.00 | 0.98, 0.86 | 1.03, 1.00 | 1.09, 1.19 | 1.04, 1.26 | 0.69, 0.68; 0.94, 0.90 | 0.97, 0.93 |
| strong large-scale modes | 1.20, 1.68 | 1.00, 1.01 | 1.44, 1.51 | 1.56, 1.75 | 1.74, 2.13 | 0.64, 0.61 | 0.89, 0.89; 0.99, 0.99 | 1.51, 1.68 |

Expected for independent blocks: CR0 sqrt(1 - 1/n_eff) = 0.95, 0.87. `t2a` confirms the three expectations to three
digits for eight share patterns and CR3 = G / (G - 1) times the jackknife exactly.

Spread of 0.1 between blocks, rms about the region truth (`out_t2b_*_tau0.1*.txt`):

| arm | sandwich / rms | CR0 / rms | CR3 / rms (sqrt mean var) | CR3 median / rms | CR3 cover +-1 | +-2 |
|---|---|---|---|---|---|---|
| white, independent offsets | 0.53, 0.50 | 0.93, 0.81 | 1.05, 1.14 | 0.99, 0.93 | 0.67, 0.64 | 0.93, 0.89 |
| dense, independent offsets | 0.39, 0.33 | 0.97, 0.73 | 1.10, 1.09 | 1.02, 0.85 | 0.70, 0.61 | 0.93, 0.86 |
| strong modes, independent offsets | 0.59, 0.54 | 1.50, 1.42 | 1.91, 2.20 | 1.78, 1.76 | 0.91, 0.91 | 1.00, 1.00 |
| white, offsets correlated over two blocks | 0.34, 0.37 | 0.53, 0.50 | 0.60, 0.70 | 0.55, 0.58 | 0.44, 0.45 | 0.73, 0.72 |

Mean of the estimate minus the region truth with a spread of 0.1: -0.0021 +- 0.0008 and -0.0006 +- 0.0013 (a second-order
shift of 0.06 rms at most).

**Answers to the questions put.**
- 2c is right for independent blocks. It fails with covariances that couple blocks; in every case I built the failure
  is conservative (block error too large), up to a factor 2. 2d is right for independent offsets and fails
  anti-conservatively for offsets correlated between blocks.
- Cross-region terms: negligible unless the covariance couples regions strongly (share 0.000, 0.005, 0.49 in the three
  arms; code over full-Jacobian CR3: 1.000, 1.000, 0.94 and 0.87).
- Clipping: see code defect 3.
- Bread: for the jackknife identity and for validity under misfit the consistent choice is the observed Hessian
  sum_B D_B (it is what h is made of, and J_c is a model quantity). Numerically it makes no difference here: the two
  agree to 1 or 2 percent in every arm (2.5 to 97.5 percent of realizations within 0.981 to 1.019) and on the maps (1.00
  to 1.01), and CR0 with either
  agrees to three digits. Keep J_c, but stop if |sum D / J_c - 1| exceeds a few percent instead of only printing it.

## 3. Low-pass estimator

**Derivation.** For any matrix L of full row rank, c = L d is Gaussian with mean L M a and covariance L N L^T, and its
log likelihood is -(d - M a)^T N'^-1 (d - M a) / 2 with N'^-1 = L^T (L N L^T)^-1 L. Using d' = Y_k L d in place of d
changes nothing iff L Y_k = 1. With L the rows of Y^+ = (Y^T Y)^-1 Y^T at the kept columns, L Y_k = 1, L Y_high = 0 and
L v = 0 for v orthogonal to span(Y): c is blind to basis modes above the cut and to everything outside the basis, which
is what "multipoles 2 to L only" requires. The other candidate, pinv(Y_k), also gives a valid likelihood but is not blind
to the modes above the cut on a grid where the columns are not orthogonal. The stack (full bands, coefficient vectors) is
an ordinary Gaussian linear model, so the flat-prior results apply with N^-1 replaced by N'^-1, provided the code uses
N^-1 only through products (it does: `scaled`, `apply`) and N'^-1 N N'^-1 = N'^-1 (true). T = tr(A^-1 M_r^T P M_r')
needs no inverse of N'^-1.

**Tests.**
- `t3a`, nside 4, cut at 5 (64 of 280 columns): |Y^+ Y - 1| 5e-15, |Y_k^+ Y_k - 1| 5e-15, |Y_k^+ Y_high| 7e-16,
  |Y_k^+ (1 - Y Y^+)| 7e-16; rank of N'^-1 64; |N'^-1 N N'^-1 - N'^-1| 2e-15; quadratic form with (d', N'^-1) against the
  Gaussian form of c: difference 4.5e-13 in 2778. pinv(Y_k) on the columns above the cut: max 0.05, rms 0.003.
- `t3a`, nside 2, exact reference with coefficient vectors (23 and 353 GHz full, three bands at l <= 3, dense noise,
  dense CMB prior, two regions): profile code minus reference 2e-11; marginal constant to 2e-11; gradient against finite
  differences 9e-7 of 36; J and T against the reference 4e-15 and 7e-16; exact score mean 2e-13, |Cov - (J + T)| 2e-15;
  Hessian of the expected objective against J 6e-7; Monte Carlo of the code's score (3000): mean within 0.6 se, variance
  / (J + T) 1.02, 0.99 (se 0.026). Information of the two full bands alone: 0.000 of the low-pass arm (two bands, two
  free amplitudes per element).
- `t3c`, nside 4, 400 realizations, LCDM CMB and prior, cut at 5: all scales mean z +0.02, -0.04, -0.02, sd z 0.99, 0.92,
  1.01; low-pass mean z -0.05, -0.06, -0.02 (se 0.05), sd z 0.95, 1.00, 1.05 (se 0.035), coverage 0.71, 0.67, 0.67; median
  error 0.159, 0.111, 0.044 against scatter 0.154, 0.112, 0.046. Mean of the low-pass estimate minus truth -0.019 +- 0.008
  and -0.013 +- 0.006 where the error is 0.16 and 0.11: a second-order shift of 0.12 sd that scales with the variance and
  is irrelevant at errors of 0.01 to 0.03. Paired difference low-pass minus all scales: sd 0.147, 0.111, 0.044 against
  (sd_low^2 - sd_all^2)^(1/2) = 0.149, 0.107, 0.042, so the rule the paper uses for nested fits holds here too.
- Contamination of 0.3 uK rms (1.3 noise sd per element) in one band, shift in units of the quoted error: in a
  low-passed band, from basis columns above the cut or from outside the span: all scales up to 2.2, low-pass 1e-14;
  in the anchor band: all scales up to 1.5, low-pass 0.04 (the region boundaries couple the anchor's small scales to the
  large-scale prediction of the other bands); with pinv(Y_k) instead: 0.04.

## Test gaps in the repository

- No unit test for `lcdm_cmb_precision`, for `cmb_prior_prec`, for `block_scores` or `block_errors`, or for the low-pass
  construction. The functions in `t1a` (analytic covariance), `rv2_common.Reference` (Woodbury reference, coefficient
  maps) and the identities of `t3a (1)` are small enough to lift into `lowres/tests`.
- The simulations validate the prior against itself (CMB drawn from the truncated covariance) and the block error
  against independent offsets and, for the null ratio, against nside 8 with a white prior. Missing: a null distribution
  of block over sandwich at the configuration of the map fit (LCDM prior; nside 16 or at least 8); offsets coherent over
  neighbouring blocks; a CMB with power above the band limit.
- Tolerances: the rms of z from 60 realizations with 4 to 8 effective blocks is uncertain by 0.1 to 0.2; it cannot
  separate CR2 from CR3 or 1.07 from 1.3.

## Not checked

- Nothing at nside 8 or 16, nothing on the WMAP or Planck maps or matrices: the condition number of Y^T Y at nside 16,
  whether any D_B is negative on the maps, the cross-region share on the maps, the null ratio with the real matrices.
- The scaling of the truncation bias to 6144 elements is an extrapolation (sqrt(n) from 384 elements, white toy sky).
- The effective window above l = 4 nside at nside 16 is taken from the exact pixel average at nside 4 at the same
  l / nside (checked against the tabulated window where it exists: 0.41 against 0.43 at l = 47, 0.19 against 0.19 at 64).
- `t1f` needs the file that `t1b` writes (`pixavg_modes_n4.npz`, 26 MB, deleted after use). The first two `t2b` outputs
  were produced before the exact-covariance line and the correlation option were added to the script; the numbers do
  not depend on either.
