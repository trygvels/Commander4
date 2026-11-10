# E041: coverage ladder for one region index, from one data set's own error estimates

- **Date started / finished:** 2026-11-10 21:03 (pre-registration) / 2026-11-10 21:09
- **Branch and commit:** `nextgen` at c331e93 plus the uncommitted `fit.sandwich_plugin`
- **Open question addressed:** Q2, Q20; Paper A1 Sect. 5 and its coverage figure
- **Validation level:** L1 and L2 (exact one-dimensional targets, known truth, many realizations)

## Question

E036 mapped the estimators on one realistic sky at three noise levels with grid posteriors, and
E039 tests the theory of the sandwich covariance with J and T evaluated at the truth. Neither
tests what is done on real data: an interval built from one data set's fitted amplitudes. Over a
ladder of per-amplitude S/N and of the number of amplitudes that share the index, how do the
intervals of (i) the profile with its curvature width, (ii) the profile with the plain plug-in
sandwich, (iii) the profile with the debiased plug-in sandwich, (iv) the flat-prior marginal at
references of 23 GHz, 30 GHz and the pivot, and (v) the zero-mean prior with a maximized scale
cover the truth, and does the flat-prior shift follow its formula across the ladder?

What changes: this is the coverage figure of Paper A1 and the test of the error recipe used in
E040. If (iii) fails where the region is measured, the real-data errors are not usable as coded.

## Target and priors

One power-law component, delta bands at 22.8, 28.4, 33.0, 40.6, 44.1 GHz in brightness
temperature, one index for n amplitudes, flat on [-4.6, -1.6], truth -2.9 (off the centre of the
range). Noise diagonal, sd per band 1.93, 1.15, 1.95, 1.77, 1.38 uK times a ladder factor s, times
a fixed log-normal factor per amplitude and band (log-sd 0.2), so that the lever arm and the pivot
differ between amplitudes. Amplitudes fixed per cell: 8.4 uK times a log-normal factor (log-sd
1.2) with random sign, heavy-tailed like the real polarized sky and drawn from none of the priors.

Objectives per amplitude p, with w = 1 / sigma^2, f_b = (nu_b / nu_ref)^beta, A_p = sum w f^2,
b_p = sum w f d: profile -(1/2) sum b_p^2 / A_p; flat-prior marginal adds (1/2) sum ln A_p;
zero-mean prior a_p ~ N(0, tau^2) with tau maximized jointly: -(1/2) sum b_p^2 tau^2 / (1 + tau^2
A_p) + (1/2) sum ln(1 + tau^2 A_p).

Errors. Profile: (i) observed curvature of the objective at the optimum; (ii) plain plug-in
sandwich J_p^-1 (J_p + T) J_p^-1 with J_p = sum u_hat V, u_hat = b_p^2 / A_p; (iii) debiased
plug-in J_c = J_p - T, variance J_p / J_c^2. Flat marginal and zero-mean prior: curvature of their
own objective (for the zero-mean prior, the beta element of the inverse 2 by 2 Hessian).

## Setup

n = 32, 128, 512, 2048, 8192; s = 0.5, 1, 2, 4, 8, 16, 32, 64 (median S/N per amplitude at 22.8 GHz
from 8.7 to 0.07); 400 noise realizations per cell; one thread. Script `nextgen/scripts/e041_coverage_ladder.py`.

## Reference

The truth is known. Formulas to be tested, all evaluated at the truth: sd of z with the
curvature width = sqrt(1 + T / J); flat-prior shift = sum_p ln(nu_ref / nu_p) / J with nu_p the
weighted geometric mean frequency of amplitude p.

**Predictions, written before the first run.** "Measured cells" are those with a sandwich sd at
the truth below 0.25.

- L1. Profile optimum: abs(mean z) below 0.15 with the debiased sandwich width in every measured
  cell.
- L2. Profile with the curvature width: sd of z equals sqrt(1 + T / J) within 0.07 in every
  measured cell.
- L3. Profile with the debiased plug-in sandwich: sd of z between 0.93 and 1.07 and 68 percent
  coverage between 0.65 and 0.72 in every measured cell. With the plain plug-in sandwich the sd of
  z exceeds 1.07 wherever T / J is above 0.15.
- L4. Flat-prior marginal: the mean shift of the optimum equals the formula within 15 percent,
  or within 0.005, wherever the predicted shift is below 0.3 in absolute value, at references of
  23 and 30 GHz; at the pivot reference the shift is below a tenth of the 23 GHz shift.
- L5. Zero-mean prior with a maximized scale: abs(mean z) below 0.2, sd of z between 0.9 and 1.15,
  coverage between 0.62 and 0.74 in every measured cell with s of 4 or less. No prediction for s of
  8 and 16.
- L6 (algebra, to be confirmed numerically in a side cell with the per-amplitude noise factor
  switched off): with a noise pattern common to all amplitudes the zero-mean-prior optimum
  coincides with the profile optimum, and its curvature variance is the profile's times
  q / (q - n), q = sum b_p^2 / A, which is the plug-in form of 1 + T / J. So in that case the
  zero-mean prior with a free scale is the profile with its sandwich width.

**Design amendment, 2026-11-10 21:04, before the run.** A smoke test of two cells (n = 32, s = 1 and 8, 40
realizations, not kept) showed T / J of 2e-4 and 1e-2: with heavy-tailed amplitudes a few bright
ones carry the information, and the ladder as first written never reaches the regime where the
sandwich matters. T / J of order 0.3 to 1 with a measured index needs many faint amplitudes, so
n = 8192 and s = 32, 64 are added. The predictions are unchanged.

## Command

```bash
cd /Users/trygvels/code/Commander4
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python nextgen/scripts/e041_coverage_ladder.py \
  --n 32 128 512 2048 8192 --s 0.5 1 2 4 8 16 32 64 --out nextgen/experiments/E041/ladder.npz > nextgen/experiments/E041/ladder.txt
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python nextgen/scripts/e041_coverage_ladder.py --common-noise \
  --n 128 --out nextgen/experiments/E041/ladder_common.npz > nextgen/experiments/E041/ladder_common.txt
```

## Result

Run 21:04 to 21:07 on one thread (`E041/ladder.txt`, `ladder.npz`; `ladder_common.txt` for the side
cell). 40 cells, 400 realizations each. 30 cells are measured (sandwich sd at the truth below
0.25). T / J at the truth runs from 2e-4 to 13; among the measured cells it reaches 1.56.

Profile optimum with three widths, a selection of cells (all 40 are in the output file):

| n | s | sandwich sd at truth | T / J | mean z | sd of z: curvature | sqrt(1 + T/J) | plain plug-in | debiased | coverage: curvature | plain | debiased |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 128 | 1 | 0.032 | 0.003 | -0.02 | 1.044 | 1.002 | 1.044 | 1.042 | 0.665 | 0.665 | 0.665 |
| 128 | 4 | 0.131 | 0.051 | -0.01 | 0.989 | 1.025 | 0.991 | 0.966 | 0.690 | 0.690 | 0.703 |
| 512 | 16 | 0.151 | 0.240 | -0.02 | 1.076 | 1.113 | 1.092 | 0.962 | 0.675 | 0.662 | 0.713 |
| 2048 | 8 | 0.049 | 0.112 | 0.00 | 1.092 | 1.055 | 1.097 | 1.035 | 0.613 | 0.605 | 0.650 |
| 2048 | 16 | 0.113 | 0.449 | +0.03 | 1.216 | 1.204 | 1.278 | 1.009 | 0.557 | 0.542 | 0.660 |
| 8192 | 16 | 0.051 | 0.389 | +0.04 | 1.167 | 1.178 | 1.216 | 0.991 | 0.585 | 0.570 | 0.682 |
| 8192 | 32 | 0.139 | 1.555 | +0.01 | 1.588 | 1.598 | 2.013 | 0.997 | 0.485 | 0.383 | 0.693 |
| 8192 | 64 (not measured) | 0.468 | 6.218 | +0.09 | 2.608 | 2.687 | 5.129 | 0.968 | 0.310 | 0.150 | 0.693 |

- **L1 pass.** Mean z with the debiased width is within 0.10 in all 30 measured cells (standard
  error 0.05).
- **L2 pass.** With the curvature width the sd of z equals sqrt(1 + T / J) within 0.054 in all 30
  measured cells (its own error is 0.035). Outside them it keeps following the formula where the
  bounds do not truncate (2.61 against 2.69 at T / J = 6.2).
- **L3 pass on the width, one cell out on the coverage, and one cell out by 0.001 on the
  control.** Debiased plug-in sandwich: sd of z between 0.946 and 1.052 in all 30 measured cells;
  coverage between 0.650 and 0.718 in 29, and 0.635 in one (n = 2048, s = 4, where the curvature
  width gives 0.632: the realizations of that cell, not the estimator; 2.1 binomial standard
  errors from 0.683). The plain plug-in sandwich is too narrow wherever T / J is above 0.15, by a
  factor 1.07 to 5.1 in sd of z, except in one cell where it is 1.069 against the 1.07 I wrote
  (n = 128, s = 8, T / J = 0.20). It is narrower than the curvature width, so it is worse than no
  correction.
- **L4 pass.** Flat-prior marginal: wherever the predicted shift is below 0.3 the measured shift
  is within 4 percent of the formula at a 23 GHz reference and within 15 percent at 30 GHz (worst
  ratio 0.855 at a predicted +0.288). At the pivot reference the shift is below 1e-3 in every
  measured cell. Beyond 0.3 the shift saturates at the edge of the allowed range. In units of the
  flat-prior posterior width the shift reaches 33 and 36 at n = 8192, with zero coverage from
  n = 512 and s = 8 on.
- **L5 fails in one cell, and the cells without a prediction show why.** Zero-mean prior with a
  maximized scale, measured cells with s of 4 or less: within the tolerances in 19 of 20; mean z
  of -0.31 at n = 2048, s = 4. At higher noise the same n gives -0.42, -0.58 and -0.53 (s = 8, 16,
  32), n = 32 gives +0.11 to +0.27, n = 512 gives +0.11 to +0.24, and n = 8192 stays within 0.15.
  The sign and size belong to the amplitude draw of each n, not to the noise level.
- **L6 pass.** With a noise pattern common to all amplitudes (eight cells) the zero-mean-prior
  optimum equals the profile optimum to the optimizer tolerance (rms difference 0.00 sd in every
  cell) and its sd of z equals that of the debiased plug-in sandwich to 0.008 (for example 1.024
  against 1.023 at n = 2048, s = 16, where the curvature width gives 1.248).

## Negative control

The flat-prior marginal fails as it must (L4). The plain plug-in sandwich fails as predicted
(L3).

## Conclusion

1. **The error recipe used for the real maps is validated in one dimension:** the profile optimum
   with the debiased plug-in sandwich, J_c = J(a_hat) - T and variance J(a_hat) / J_c^2, has sd of z
   within 0.055 of 1 and nominal coverage wherever the index is measured, for T / J up to 1.6,
   from one data set alone.
2. **The curvature of the profile likelihood is too narrow by sqrt(1 + T / J)**, confirmed
   over four decades of T / J. This is the published observation that the profile likelihood is
   overconfident at low S/N (F12, F18), with its size.
3. **The plain plug-in sandwich is wrong in the harmful direction.** It was what `fit.sandwich`
   would have reported on real data before the amendment of 2026-11-10.
4. **The zero-mean prior with a free scale is the profile with its sandwich width when the noise
   pattern is common to the amplitudes, and is biased when it is not.** With heavy-tailed
   amplitudes and a noise pattern that varies between amplitudes, the Gaussian prior is
   misspecified and the few bright amplitudes, with whatever noise pattern they happen to have,
   set a bias of up to 0.6 posterior sd whose sign depends on the sky. E036 found the same
   prior calibrated within 0.11 except for -0.30 in one cell; that cell was this effect. The
   profile score is unbiased for any fixed amplitudes, so the profile with the debiased sandwich
   is the estimator of record for E040, and the zero-mean prior is reported beside it.
5. **The flat-prior shift follows its formula** to 4 percent (23 GHz) and 15 percent (30 GHz) up
   to shifts of 0.3, and vanishes at the pivot reference.

What this does not establish: several coupled regions, dense noise, or more than one component
(E039 and its plug-in addendum); non-Gaussian noise; the zero-mean prior with a scale map.

Changes: `fit.sandwich_plugin` is the error estimate of E040 (amendment in that record); F44's
sentence that the zero-mean prior with a free scale is calibrated is qualified; Paper A1 Sect. 5.


## Follow-ups
