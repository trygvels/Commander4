# E036: amplitude priors for signed Q/U amplitudes at real WMAP/LFI depth, coverage of a region index

- **Date started / finished:** 2026-11-09 23:07 / 2026-11-10 00:39
- **Branch and commit:** `nextgen` at 1eaac33
- **Open question addressed:** Q2 (calibrated priors) in the regime E031 identified as the one that
  matters, Q3; the design assumption of the proposed polarized track
  (`../reviews/2026-11-09_program_review.md`)
- **Validation level:** L4 in miniature (a realistic amplitude field that is not drawn from any of
  the priors, known index, many noise realizations), exact grid posteriors, no sampler

## Question

Every amplitude-prior result so far (F12, F15, F16, F18) is for intensity, positive templates and
a prior whose mean is tied to a template. Polarized amplitudes are signed, the per-pixel S/N at
real depth is a few, and the frequency lever arm is short. For one spectral index per region, at
real WMAP/LFI Q/U white-noise depth and on a realistic polarized synchrotron amplitude field,
which target gives an unbiased and calibrated index posterior: the flat-prior marginal, the
profile likelihood, or a zero-mean Gaussian amplitude prior with an inferred scale? And does the
zero-mean prior need a template for its scale map?

What changes with the answer: if the zero-mean prior with a free scale per region is calibrated
without a template, it becomes the amplitude prior of the polarized track and the template
question of F21 to F38 does not arise there. If it is not, the track has no validated amplitude
prior and that is the first thing to solve.

## Target and priors

Per pixel p and Stokes parameter, bands b: d_pb = a_p (nu_b / nu_ref)^beta + n_pb, with n_pb
Gaussian, independent, variance sigma_pb^2 known. One index beta per region, flat on the grid
[-4.6, -1.6]. Write w = 1 / sigma^2, f_b = (nu_b / nu_ref)^beta, A_p = sum_b w f^2,
b_p = sum_b w f d, S_p = sum_b w d^2.

- (a) Flat-prior marginal: log L = sum_p [-(S_p - b_p^2 / A_p) / 2 - (1/2) log A_p]. Depends on
  nu_ref through log A_p: moving the reference from nu_1 to nu_2 adds beta ln(nu_1 / nu_2) per
  amplitude. Evaluated at nu_ref = 23, 30 and 60 GHz.
- (b) Profile: log L = sum_p -(S_p - b_p^2 / A_p) / 2. Independent of nu_ref.
- (c) Zero-mean prior, one scale per region: a_p ~ N(0, tau_r^2), prior flat in log tau_r on
  [0.05, 2000] uK_RJ: log L = sum_p [-(S_p - b_p^2 tau^2 / (1 + tau^2 A_p)) / 2 - (1/2) log(1 +
  tau^2 A_p)], marginalized over tau_r on a grid. Independent of nu_ref up to the grid range.
- (c1) The same with a scale map: tau_p = s_r g_p, g_p a template of the polarized intensity
  carrying a log-normal error (rms 0.3, FWHM 10 deg), s_r flat in log on [0.03, 30].
- (c2) The same with the oracle scale map g_p = abs(a_p), which no analysis has.

## Setup

- Truth: the Cosmoglobe DR1 polarized synchrotron posterior mean (Q, U at 30 GHz, uK_RJ),
  smoothed to 5 deg and averaged to nside 16, scaled to the bands with one index beta = -3.1 over
  the sky. It is a data-derived field with the real sky's morphology and heavy tails, and is drawn
  from none of the priors.
- Bands: WMAP K, LFI 30, WMAP Ka, WMAP Q, LFI 44 at 22.8, 28.4, 33.0, 40.6, 44.1 GHz, with the
  real Q/U white-noise rms per nside-16 pixel from the DR1 rms maps (E031: medians 1.96, 1.15,
  1.99, 1.85, 1.40 uK_RJ). Synchrotron only: no CMB, no dust, white noise. Those come in the next
  step and can only make things harder.
- Regions: nside 2 (48 regions, 128 amplitudes each) and nside 4 (192 regions, 32 each). 100
  noise realizations, seeds fixed. Grid step 0.01 in beta.
- Regions are classed by their Fisher sd at the truth with profiled amplitudes, which does not
  depend on the noise realization: below 0.1, 0.1 to 0.25, above 0.25.
- One thread, numpy.

## Reference

The truth is known (beta = -3.1 everywhere). For each estimator, region class and region size:
mean and sd of z = (posterior mean - truth) / posterior sd over realizations and regions, the
mean shift in index units, and the coverage of the central 68 and 95 percent intervals.

**Predictions, written before the first run:**

- P1. (c), zero-mean prior with a free scale per region: mean z within +-0.15 and sd of z between
  0.9 and 1.15, coverage of the 68 percent interval between 63 and 73 percent, in the two
  data-dominated classes at nside 2. At nside 4 (32 amplitudes) I allow +-0.25. Reason: with a
  free scale and noise that is nearly uniform inside a region, the Gaussian likelihood depends on
  the data through the band covariance only, whose expectation the model reproduces for any
  amplitude distribution, and normal-theory inference on loadings is known to be robust to
  non-normal factors when the errors are normal (from memory: Anderson and Amemiya 1988, Browne
  and Shapiro 1988; to be verified).
- P2. (a), flat-prior marginal: biased, with a sign set by the reference. Steeper than the truth
  at nu_ref = 23 (below the information-weighted band), flatter at nu_ref = 60. Mean z beyond
  +-0.5 in the 0.1 to 0.25 class at nside 2, and the three references differ from each other by
  more than one posterior sd there.
- P3. (b), profile: mean z within +-0.3, but intervals too narrow where the per-pixel S/N is low:
  sd of z above 1.1 and 68 percent coverage below 62 percent in the 0.1 to 0.25 class.
- P4. (c1), scale map with a 30 percent error: no worse than (c) by more than 0.15 in mean z, and
  no better. If so the template adds nothing and is dropped.
- P5. (c2), oracle scale map: calibrated within the P1 tolerances, with a posterior sd at most 10
  percent below that of (c).
- Decision rule: if P1 holds, (c) is adopted for the polarized track, subject to the same test
  with CMB, dust and correlated noise. If P1 fails by more than twice its tolerance, the cause is
  found before anything is built on it.

**Follow-up cells, predictions written 2026-11-09 23:18 after the real-depth run was read and before
these were run.** The real-depth run showed the profile calibrated (P3 failed), so the test did
not separate the profile from the zero-mean prior. A known property explains it: for linear
amplitudes with known noise the profile score is exactly unbiased, and its curvature differs from
the score variance by a term of relative order 1 / (S/N)^2 per amplitude, which is 5 percent at
the median S/N of 4.45. The same sky with the noise multiplied by 2 and by 4 (per-element S/N of
2.2 and 1.1 at K band; 4 corresponds to nside-64 pixels at real depth) should separate them.

- F1, noise times 2: profile sd of z between 1.1 and 1.3 and 68 percent coverage below 0.64 in
  the classes with Fisher sd above 0.1; zero-mean prior sd of z within 0.9 to 1.15.
- F2, noise times 4: profile sd of z above 1.3 and coverage below 0.60 in every class above 0.1;
  zero-mean prior within 0.9 to 1.15 with coverage 0.63 to 0.73 and mean z within +-0.2.
- F3, both: the flat-prior marginal shifts scale as the noise variance (4 and 16 times the
  real-depth shifts in index units) until the grid edge, following
  shift = ln(nu_ref / nu_pivot) sigma_1^2 with sigma_1 the index sd of a single amplitude.

**Second stage, predictions written 2026-11-09 23:21 before any of its cells was run: a free CMB and a
free dust amplitude in every pixel and Stokes parameter.** E022 found that the amplitude-index
coupling of the BP XV polarized configuration is carried by the free CMB amplitude (F22). The
model gains two amplitudes per element with known SEDs (CMB; dust as a modified blackbody with
index 1.55 and 19.6 K) and flat priors, and the bands V, 70 and 353 GHz. Amplitudes whose SED does
not depend on the index contribute no index-dependent volume term, so they are projected out
exactly and the three estimators keep their form with a projected weight matrix. Dust truth: the
DR1 polarized dust map (referenced at 70 GHz, F13), 5 deg, nside 16. CMB truth: zero, because at
this resolution its signal is far below the noise and only its free amplitude matters.

- G1. The region Fisher sd grows by a factor 1.2 to 2.5 against the synchrotron-only model
  (median 0.060 at nside 2).
- G2. The profile stays calibrated at real depth: mean z within +-0.15, sd of z at most 1.15 and
  68 percent coverage at least 0.62 in the classes below 0.25.
- G3. The zero-mean synchrotron prior with flat CMB and dust priors stays within the P1
  tolerances.
- G4. The flat-prior marginal follows the same formula with the larger single-amplitude
  variance, so its shifts grow by the square of the G1 factor.
- G5. With a Gaussian CMB prior of 0.2 uK per element, the size of the signal, the Fisher sd
  returns to within 10 percent of the synchrotron-only value (F22: ratio 1.02 to 1.04).
- The code path changed to a general weight matrix for this stage. A rerun of the real-depth
  cell with the new code must reproduce the first table to the last printed digit.

**Third stage, predictions written 2026-11-09 23:44 before its cells were run.** Two points from the
independent challenge (`../reviews/2026-11-09_way_forward_challenge.md` Sect. E).

(a) *The truth sat at the centre of the grid.* With the truth at -3.1 in a flat prior on
[-4.6, -1.6], reversion to the prior reads as zero bias. Cells with the truth at -2.8 and at
-3.4, synchrotron only, real depth.

- H1. In the classes below 0.25 the profile and the zero-mean prior are unchanged: shift within
  +-0.02 and sd of z within 0.9 to 1.15.
- H2. In the class above 0.25 at nside 4, where 14 to 17 percent of the posteriors touch the grid
  edge, the shift moves toward the grid centre: between -0.05 and -0.12 for a truth of -2.8, and
  between -0.05 and +0.05 for -3.4 (the -0.05 seen at -3.1 is a skew toward steep indices that
  the reversion adds to or works against).

(b) *One scale shared by all index regions on a scale map*, the arm the challenge showed to be
unsound in a two-zone toy and the closest polarized relative of the template prior of F21 to F38.
The prior is a ~ N(0, (s g_p)^2) with one s for the sky, stated at a reference of 23 GHz, at the
information pivot, or at 30 GHz, with the scale map moved between references with the true
index. Scale maps: the true local rms of the amplitudes (10 deg), and the same 30 percent high at
abs(b) > 30 deg. 50 realizations, real depth, truth -3.1, shifts per latitude zone for regions
with Fisher sd below 0.25.

- H3. With the true scale map the shared-scale prior is unbiased at every reference: zone shifts
  within +-0.02.
- H4. With the map 30 percent high at high latitude the high-latitude zone shifts with a sign set
  by the reference: negative at 23 GHz, positive at 30 GHz, within +-0.02 at the pivot, with a
  size between 0.03 and 0.3 at the two off-pivot references.
- H5. The low-latitude zone moves by less than 0.02 in every arm.
- If H4 holds, the displaced faint sky with a reference-dependent sign of F38 and E030 chain 3
  has a reproducible relative in an exact grid experiment, and "never one scale across index
  regions on a template" is a tested rule and not only algebra.

## Command

```bash
cd /Users/trygvels/code/Commander4
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -I nextgen/scripts/e036_qu_amplitude_priors.py \
  > nextgen/experiments/E036/coverage_real_depth.txt
for s in 2 4; do OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -I \
  nextgen/scripts/e036_qu_amplitude_priors.py --noise-scale $s > nextgen/experiments/E036/coverage_noise_x$s.txt; done
# second stage (script extended with --nuisance)
P="OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -I nextgen/scripts/e036_qu_amplitude_priors.py"
eval $P > nextgen/experiments/E036/coverage_real_depth_regression.txt
eval $P --nuisance cmb_dust > nextgen/experiments/E036/coverage_cmb_dust_flat.txt
eval $P --nuisance cmb_dust --cmb-prior-sd 0.2 > nextgen/experiments/E036/coverage_cmb_dust_cmbprior.txt
# third stage (script extended with --global-arm)
eval $P --beta -2.8 > nextgen/experiments/E036/coverage_truth_m28.txt
eval $P --beta -3.4 > nextgen/experiments/E036/coverage_truth_m34.txt
eval $P --global-arm --nreal 50 > nextgen/experiments/E036/coverage_global_scale.txt
```

## Result

**Real depth** (`E036/coverage_real_depth.txt`, 100 noise realizations, finished 23:17). Median
Q/U noise per nside-16 pixel 1.93, 1.15, 1.95, 1.77, 1.38 uK_RJ; truth amplitude at 23 GHz median
8.4 uK_RJ per element (16 and 84 percent points 2.1 and 24.7); median S/N per element at K band
4.45. Shift is the mean of posterior mean minus truth in index units; z in posterior sd.

Regions at nside 2 (128 amplitudes each):

| class by Fisher sd | estimator | shift | mean z | sd of z | 68 percent coverage | 95 percent |
|---|---|---|---|---|---|---|
| below 0.03 (9 regions) | flat marginal, ref 23 / 30 / 60 GHz | -0.005 / +0.007 / +0.037 | -0.26 / +0.35 / +1.94 | 1.02 / 1.03 / 1.18 | 0.69 / 0.66 / 0.21 | 0.94 / 0.95 / 0.51 |
| | profile | 0.000 | -0.01 | 1.02 | 0.69 | 0.96 |
| | zero-mean, free scale | -0.001 | -0.05 | 1.02 | 0.69 | 0.96 |
| 0.03 to 0.1 (28) | flat marginal, ref 23 / 30 / 60 | -0.061 / +0.062 / +0.341 | -0.94 / +1.01 / +6.06 | 1.00 / 1.03 / 2.12 | 0.49 / 0.48 / 0.00 | 0.84 / 0.83 / 0.01 |
| | profile | -0.002 | -0.03 | 0.97 | 0.71 | 0.96 |
| | zero-mean, free scale | -0.002 | -0.06 | 0.98 | 0.70 | 0.96 |
| 0.1 to 0.25 (11) | flat marginal, ref 23 / 30 / 60 | -0.401 / +0.341 / +1.32 | -2.37 / +2.53 / +33.9 | 1.17 / 1.09 / 23 | 0.11 / 0.08 / 0.00 | 0.36 / 0.30 / 0.00 |
| | profile | -0.019 | -0.07 | 1.05 | 0.67 | 0.94 |
| | zero-mean, free scale | +0.002 | +0.06 | 1.10 | 0.64 | 0.93 |
| | zero-mean, scale map with 30 percent error | -0.018 | -0.06 | 1.02 | 0.68 | 0.94 |
| | zero-mean, oracle scale map | -0.023 | -0.10 | 1.01 | 0.68 | 0.95 |

Regions at nside 4 (32 amplitudes each): the profile and the three zero-mean variants agree to
0.001 in shift in every class, with mean z between -0.01 and -0.11, sd of z 0.97 to 1.08 and 68
percent coverage 0.67 to 0.72, up to the class above 0.25 (39 regions, median Fisher sd 0.40).
The flat-prior marginal shifts there are -0.019 / +0.020 / +0.116 (class 0.03 to 0.1), -0.113 /
+0.098 / +0.551 (0.1 to 0.25) and -0.62 / +0.52 at the grid edge (above 0.25).

Verdicts on the predictions:

- **P1 pass.** The zero-mean prior with a free scale per region: mean z -0.05, -0.06, +0.06, sd
  of z 1.02, 0.98, 1.10, coverage 0.69, 0.70, 0.64 at nside 2; at nside 4 mean z within 0.11.
- **P2 pass, and stronger than predicted.** The flat-prior marginal is biased with a sign set by
  the reference: steeper than the truth at 23 GHz, flatter at 30 and 60 GHz. Moving the reference
  from 23 to 30 GHz moves the posterior mean by 0.012, 0.12 and 0.74 in the three classes. The
  shifts follow shift = ln(nu_ref / nu_pivot) sigma_1^2, with sigma_1 the index sd of a single
  amplitude and nu_pivot the information-weighted geometric mean frequency (26.4 GHz here): for
  the middle class sigma_1^2 = 0.059^2 x 128 = 0.446 and ln(30 / 23) x 0.446 = 0.118, against the
  measured 0.122. The shift does not average down with the number of pixels that share the index.
- **P3 fail.** I predicted the profile's intervals too narrow in the 0.1 to 0.25 class. They are
  calibrated: sd of z 1.05 and coverage 0.67 at nside 2, 1.01 and 0.67 at nside 4, 1.08 and 0.67
  above 0.25. At a per-element S/N of 4.45 the profile needs no help.
- **P4 pass on "no worse", fail on "no better" in one class.** The scale map with a 30 percent
  error equals the free scale in mean z everywhere. In the 0.1 to 0.25 class at nside 2, where one
  scale must cover a 30 degree region, it is slightly better calibrated (sd of z 1.02 against
  1.10, coverage 0.68 against 0.64). At nside 4 the two are identical.
- **P5 pass.** The oracle scale map is calibrated and its posterior sd is within 2 percent of the
  free scale's: a scale map buys nothing in precision here.

**Noise times 2** (`E036/coverage_noise_x2.txt`, finished 23:26, run with the original script;
median S/N per element at K band 2.22). Classes are by Fisher sd, so a class holds different
regions than at real depth.

| regions | class | profile: shift, mean z, sd of z, cov68 | zero-mean, free scale | zero-mean, scale map | flat marginal shift, ref 23 / 30 |
|---|---|---|---|---|---|
| nside 2 | 0.03 to 0.1 (17) | -0.003, -0.03, 0.97, 0.70 | -0.012, -0.17, 0.98, 0.69 | -0.003, -0.03, 0.97, 0.70 | -0.063 / +0.071 |
| nside 2 | 0.1 to 0.25 (21) | -0.009, -0.04, 1.00, 0.68 | -0.002, 0.00, 1.04, 0.67 | -0.010, -0.04, 0.98, 0.70 | -0.378 / +0.333 |
| nside 2 | above 0.25 (8) | -0.073, -0.17, 1.27, 0.60 | -0.044, -0.04, 1.21, 0.62 | -0.068, -0.14, 1.08, 0.67 | grid edge |
| nside 4 | 0.03 to 0.1 (24) | -0.002, 0.00, 1.00, 0.67 | -0.002, -0.01, 1.00, 0.67 | -0.002, -0.01, 1.00, 0.67 | -0.023 / +0.028 |
| nside 4 | 0.1 to 0.25 (68) | -0.014, -0.04, 1.00, 0.68 | -0.018, -0.06, 0.99, 0.69 | -0.016, -0.06, 0.98, 0.69 | -0.142 / +0.124 |
| nside 4 | above 0.25 (98) | -0.047, -0.09, 1.09, 0.67 | -0.049, -0.09, 0.94, 0.74 | -0.048, -0.09, 0.94, 0.74 | grid edge |

- **F1 fails in three of four cells.** I predicted the profile's sd of z between 1.1 and 1.3 with
  coverage below 0.64 in the classes above 0.1. It is 1.00 and 0.68 in the 0.1 to 0.25 classes at
  both region sizes, 1.09 and 0.67 above 0.25 at nside 4, and reaches 1.27 and 0.60 only above
  0.25 at nside 2. The profile holds at a per-element S/N of 2.2.
- The zero-mean prior stays within its band except above 0.25 at nside 2 (sd of z 1.21, eight
  regions), where the scale map version is better (1.08).
- **F3 was posed wrongly and holds in its corrected form.** At a fixed Fisher sd the flat-prior
  shift is the same as at real depth, because the formula can be written as ln(nu_ref / nu_pivot)
  times the number of amplitudes times the region variance. For 128 amplitudes it predicts -0.072
  at a Fisher sd of 0.064 and -0.34 at 0.140, against the measured -0.063 and -0.378.

**Regression check of the extended script** (23:34 to 23:43): the real-depth table is reproduced
to the last printed digit (`E036/regression_diff.txt` is empty). The first attempt at the noise
times 4 cell crashed on a variable name clash that the extension introduced
(`E036/coverage_noise_x4_crash.err`); it was rerun with the corrected script.

**Noise times 4** (`E036/coverage_noise_x4.txt`, finished 23:52; median S/N per element at K band
1.11, the regime of nside-64 pixels at real depth).

| regions | class | profile: shift, mean z, sd of z, cov68 | zero-mean, free scale | zero-mean, scale map | flat marginal shift, ref 23 / 30 |
|---|---|---|---|---|---|
| nside 2 | 0.03 to 0.1 (8) | -0.003, -0.03, 1.01, 0.68 | -0.014, -0.18, 1.01, 0.68 | -0.005, -0.06, 1.01, 0.69 | -0.073 / +0.089 |
| nside 2 | 0.1 to 0.25 (18) | -0.016, -0.05, 1.01, 0.68 | -0.061, -0.30, 1.03, 0.64 | -0.026, -0.09, 0.97, 0.70 | -0.516 / +0.438 |
| nside 2 | above 0.25 (22) | -0.057, -0.11, 1.40, 0.58 | -0.024, 0.00, 0.99, 0.73 | -0.048, -0.09, 0.95, 0.74 | grid edge |
| nside 4 | 0.1 to 0.25 (35) | -0.018, -0.05, 1.02, 0.66 | -0.022, -0.07, 1.00, 0.67 | -0.021, -0.06, 1.00, 0.67 | -0.148 / +0.145 |
| nside 4 | above 0.25 (152) | -0.042, -0.07, 1.08, 0.69 | -0.049, -0.09, 0.84, 0.80 | -0.048, -0.09, 0.84, 0.80 | grid edge |

- **F2 fails for the profile in three of four cells, in the profile's favour.** I predicted sd of
  z above 1.3 and coverage below 0.60 in every class above 0.1. That happens only above 0.25 at
  nside 2 (1.40 and 0.58). Elsewhere the profile stays at 1.01 to 1.08 and 0.66 to 0.69, at a
  per-element S/N of 1.1.
- **F2 fails for the zero-mean prior with a free scale in one cell.** In the 0.1 to 0.25 class at
  nside 2 its mean z is -0.30 (shift -0.061), outside the +-0.2 I allowed. With a scale map, even
  one carrying a 30 percent error, the shift is -0.026. One scale for a 30 degree region whose
  amplitudes span an order of magnitude, under noise that is not uniform inside the region, is
  where the free scale is exposed; the challenge report lists this as its residual risk (ii).
  Above 0.25 at nside 4 the zero-mean priors over-cover (0.80), with 47 percent of the posteriors
  touching the grid edge.
- The flat-prior marginal is 2.5 to 8 posterior sd from the truth in every class above 0.1.

**Second stage, free CMB and dust amplitudes with flat priors** (`E036/coverage_cmb_dust_flat.txt`,
finished 2026-11-10 00:02; eight bands, 353 GHz noise 0.19 uK_RJ; truth dust amplitude median 2.15
uK_RJ at 353 GHz).

| regions | class | profile: shift, mean z, sd of z, cov68 | zero-mean, free scale | flat marginal shift, ref 23 / 30 |
|---|---|---|---|---|
| nside 2 | below 0.03 (3) | 0.000, 0.00, 0.91, 0.72 | -0.001, 0.00, 0.91, 0.72 | +0.001 / +0.016 |
| nside 2 | 0.03 to 0.1 (22) | -0.003, -0.04, 1.01, 0.68 | -0.006, -0.07, 1.01, 0.68 | +0.006 / +0.167 |
| nside 2 | 0.1 to 0.25 (17) | +0.003, +0.05, 1.02, 0.66 | +0.011, +0.08, 1.00, 0.67 | +0.058 / +0.859 |
| nside 2 | above 0.25 (6) | -0.007, 0.00, 1.08, 0.62 | -0.025, -0.06, 1.00, 0.68 | grid edge |
| nside 4 | 0.03 to 0.1 (34) | -0.003, -0.03, 1.00, 0.69 | -0.003, -0.03, 0.99, 0.69 | +0.002 / +0.049 |
| nside 4 | 0.1 to 0.25 (78) | -0.006, -0.01, 1.02, 0.68 | -0.009, -0.03, 1.01, 0.68 | +0.007 / +0.296 |
| nside 4 | above 0.25 (78) | -0.015, -0.02, 1.05, 0.68 | -0.018, -0.03, 0.94, 0.74 | grid edge |

- **G1 pass.** The region Fisher sd grows by 1.63 (median 0.098 against 0.060 at nside 2; 0.213
  against 0.130 at nside 4).
- **G2 pass.** The profile stays calibrated with two more free amplitudes per element: mean z
  within 0.05, sd of z 0.91 to 1.08, coverage 0.62 to 0.72 in every class.
- **G3 pass.** The zero-mean synchrotron prior with flat CMB and dust priors: mean z within 0.08,
  sd of z 0.91 to 1.01, coverage 0.67 to 0.74.
- **G4 fails as posed and holds in a corrected form.** I assumed the pivot would stay at 26.4
  GHz. With the CMB and dust projected out, the index information sits at the lowest band and
  the pivot moves to about 22 to 23 GHz. So a flat prior stated at 23 GHz is now nearly unbiased
  (+0.006 and +0.058) while one stated at 30 GHz is biased flat by +0.17 at an index error of
  0.064 and by +0.86 at 0.13. The difference between the two references, 0.161 in the middle
  class, is the formula's ln(30 / 23) times the single-amplitude variance (0.139 from the class
  median). The shift is still ln(nu_ref / nu_pivot) times that variance; the pivot depends on
  which other amplitudes are free. E022 measured +0.27 (Spur) and +0.13 (Plane) between the
  flat-prior marginal and the profile on real data at a 30 GHz reference with a free CMB
  amplitude, the same sign and size.

**Second stage with a Gaussian CMB prior of 0.2 uK per element**
(`E036/coverage_cmb_dust_cmbprior.txt`, finished 00:10).

- **G5 pass.** The region Fisher sd returns to the synchrotron-only value: median 0.059 against
  0.060 at nside 2 and 0.126 against 0.130 at nside 4. With a flat CMB prior it was 0.098 and
  0.213. A free CMB amplitude with a flat prior costs a factor 1.63 in the index error and one
  with a prior the size of the signal costs nothing.
- The flat-prior marginal returns to its synchrotron-only pattern (-0.060 at 23 GHz and +0.063
  at 30 GHz in the 0.03 to 0.1 class at nside 2), so the pivot is back at 26.4 GHz. The pivot of
  the flat-prior shift is set by the CMB prior as much as by the bands.
- The profile and the zero-mean prior stay calibrated: mean z within 0.16, sd of z 0.96 to 1.10,
  coverage 0.65 to 0.71 in every class.

**Third stage (a), off-centre truths** (`E036/coverage_truth_m28.txt`, `coverage_truth_m34.txt`,
finished 00:21 and 00:26; synchrotron only, real depth, 100 realizations each).

| regions | class | profile shift at truth -2.8 / -3.1 / -3.4 | zero-mean, free scale | profile sd of z |
|---|---|---|---|---|
| nside 2 | 0.03 to 0.1 | -0.002 / -0.002 / -0.003 | +0.001 / -0.002 / -0.006 | 0.97 / 0.97 / 0.97 |
| nside 2 | 0.1 to 0.25 | -0.021 / -0.019 / -0.018 | -0.014 / +0.002 / -0.020 | 1.07 / 1.05 / 1.05 |
| nside 4 | 0.03 to 0.1 | -0.002 / -0.002 / -0.002 | -0.002 / -0.002 / -0.003 | 0.99 / 0.99 / 0.99 |
| nside 4 | 0.1 to 0.25 | -0.012 / -0.012 / -0.011 | -0.013 / -0.013 / -0.013 | 1.01 / 1.01 / 1.01 |
| nside 4 | above 0.25 | -0.089 / -0.051 / -0.010 | -0.101 / -0.053 / -0.001 | 1.08 / 1.08 / 1.07 |

- **H1 pass.** In the classes below 0.25 nothing moves with the truth: shifts agree to 0.003 for
  the profile, and the one value outside +-0.02 (-0.021) is the same at all three truths.
- **H2 pass.** In the class above 0.25 the shift follows the distance of the truth from the grid
  centre, -0.089, -0.051, -0.010, a reversion of about 13 percent on top of a skew of -0.05. So
  the calibration of the data-dominated classes does not come from the truth sitting at the
  centre, and the numbers of the faintest class are prior-dependent, as the challenge said.
- The flat-prior marginal behaves the same at all three truths. For the 24 high-latitude regions
  at nside 2 with Fisher sd below 0.25 it is -0.22 / +0.15 at references of 23 / 30 GHz for a
  truth of -2.8 and -0.18 / +0.19 for -3.4, against -0.011 and -0.010 for the profile.

**Third stage (b), the shared-scale arm** (`E036/coverage_global_scale.txt`, finished 00:38; 50
realizations, real depth, truth -3.1; information pivot 26.03 GHz). Mean shift of the regions
with Fisher sd below 0.25, by latitude zone; errors over realizations are 0.002 to 0.003.

| regions, zone | profile | free scale per region | shared scale, true map: ref 23 / pivot / 30 | shared scale, map 30 percent high at abs(b) > 30: ref 23 / pivot / 30 |
|---|---|---|---|---|
| nside 2, abs(b) > 30 (24) | -0.010 | +0.002 | -0.004 / -0.009 / -0.013 | -0.030 / -0.009 / +0.017 |
| nside 2, abs(b) <= 30 (24) | -0.001 | -0.003 | -0.006 / -0.002 / +0.002 | +0.006 / -0.003 / -0.012 |
| nside 4, abs(b) > 30 (78) | -0.007 | -0.007 | -0.001 / -0.007 / -0.013 | -0.013 / -0.007 / +0.001 |
| nside 4, abs(b) <= 30 (75) | -0.007 | -0.009 | -0.013 / -0.009 / -0.002 | -0.001 / -0.009 / -0.018 |

- **H3 pass.** With the true scale map the shared-scale prior stays within 0.013 of the truth at
  every reference.
- **H4: the sign and the null hold, the size is below what I predicted.** With the map 30 percent
  high at high latitude, the high-latitude zone moves by -0.026 at a 23 GHz reference and by
  +0.029 at 30 GHz relative to the true-map arm, and by 0.000 at the pivot (nside 2). At nside 4
  the same differences are -0.012, 0.000 and +0.014. I predicted 0.03 to 0.3, so three of the
  four sizes miss the range from below. The challenge report's formula gives about 0.045 for a
  prior-to-noise variance ratio of 30, which is where these data-informative regions sit; the
  larger shifts it quotes are for fainter regions, which this zone table excludes.
- **H5 pass.** The low-latitude zone moves the other way by less than 0.02.
- The effect is coherent over a zone: each zone mean is 8 standard errors from its true-map
  value. For a latitude trend of the index it is a bias of the order of the pooled error.

## Negative control

The flat-prior marginal is the control and fails as it must. At real depth its mean z at
references of 23 and 30 GHz is 0.9 to 2.5 for the 30 degree regions above a Fisher sd of 0.03 and
0.3 to 1.9 for the 15 degree regions, and at a 60 GHz reference it is 6 to 34.

## Conclusion

1. **The amplitude prior is not the problem of the polarized track, provided it is not flat.**
   At real WMAP/LFI Q/U white-noise depth, on a realistic amplitude field, the profile likelihood
   and a zero-mean Gaussian prior with one free scale per index region both give calibrated
   region indices: mean z within 0.11, 68 percent coverage 0.64 to 0.72, down to region index
   errors of 0.4. Neither needs a template. This holds with free CMB and dust amplitudes, with
   the truth moved off the grid centre, and down to four times the real noise except in the
   faintest 30 degree regions.
2. **The flat-prior marginal is shifted by ln(nu_ref / nu_pivot) times the index variance of a
   single amplitude.** The shift does not average down with the number of amplitudes sharing the
   index, flips sign between references on either side of the pivot, and the pivot depends on
   the bands and on the CMB prior: 26.4 GHz with a CMB prior or no CMB, the lowest band with a
   flat CMB prior. Sizes at real depth: 0.06 to 0.4 for synchrotron only, +0.17 to +0.86 at a 30
   GHz reference with a flat CMB prior.
3. **The CMB prior matters, the synchrotron prior hardly.** A free CMB amplitude under a flat
   prior costs a factor 1.63 in the index error; with a prior the size of the signal it costs
   nothing.
4. **A scale shared across index regions on a scale map is biased by a wrong map with a
   reference-dependent sign**, as the challenge report derived, by 0.01 to 0.03 here for a 30
   percent error. It is the polarized relative of the template prior of F21 to F38, and it is
   small at this depth for regions the data inform.
5. Five of my predictions failed. Four of them underestimated the profile likelihood (P3, F1,
   F2) or overestimated an effect (H4), and one assumed a pivot that moved (G4). None of the
   failures favours the prior I set out to justify over the simpler profile.

What this does not establish: anything with correlated noise or a dense covariance, a spatially
varying truth index, or real data. The sandwich factor for the profile's width at low S/N was not
checked region by region. The truth amplitudes are a posterior mean of real data.

Changes: F44 added; Q20 answered for white noise; D029 item 5; Paper A outline thesis item 2.

## Follow-ups

- The same grid machinery on the E022 surrogate with dense WMAP covariances and K band, at
  references of 23 and 30 GHz, with and without a CMB prior. It would test whether the pivot
  explains why the BeyondPlanck chains show no flat-prior pull (Q17): for the BP XV band set
  without K band the pivot is near 30 GHz, the reference BeyondPlanck used.
- The sandwich factor computed per region and overlaid on these tables and on E012 and E017.
- A truth with a latitude step in the index, and a data-built scale map.
