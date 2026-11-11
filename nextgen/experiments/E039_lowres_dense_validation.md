# E039: region indices of polarized synchrotron with the released WMAP noise covariances, validation on simulations

- **Date started / finished:** 2026-11-10 20:32 / running
- **Branch and commit:** `nextgen` at 747ed50
- **Open question addressed:** Q21, Q20; Gates 1 and 3 of the challenge report, reduced to regions
- **Validation level:** L0 and L1 (unit tests against closed forms), then L4 in the sense of a
  realistic amplitude field that is drawn from no prior, with the real noise covariance

## Question

Paper A1's data section estimates one synchrotron index per sky region from WMAP K, Ka, Q, V and
Planck 353 GHz at nside 16, with the dense WMAP inverse covariances as released. Before any real
index is read: does the dense collapsed likelihood of `nextgen/lowres/` return unbiased region
indices with calibrated errors when the noise is drawn from those covariances, and do the three
targets behave as E036 found with white noise (profile calibrated once its width carries the
sandwich factor, flat-prior marginal shifted by the reference frequency, zero-mean prior equal to
the profile)?

What changes with the answer: a pass clears the estimator for the real maps (E040). A fail in
bias or coverage is found and fixed first.

## Target and priors

Per band b, in thermodynamic units, elements p (Q then U): d_b = s_b(beta_r(p)) a_s + f_b a_d + a_c
+ n_b, n_b ~ N(0, N_b). s_b is the band-averaged power law at nu_ref = 23 GHz (WMAP nine-year
bandpass profiles; delta band at 353 GHz), f_b the band-averaged modified blackbody (index 1.53,
19.6 K, reference 353 GHz), the CMB column is 1. N_b^-1: the released WMAP res 4 matrices; for 353
GHz the PR4 per-pixel Q/U covariance combined over sub-pixels. Priors: flat on a_d; N(0, 0.2^2
uK^2) per element on a_c; a_s flat (profile and flat-prior marginal) or N(0, tau_r^2) with tau_r
maximized jointly (zero-mean prior). beta_r flat on [-4.6, -1.6]. Objectives and the sandwich
covariance J^-1 (J + T) J^-1 are in `nextgen/lowres/likelihood.py`; derivation of T below.

Derivation of the score covariance. With P = N^-1 - N^-1 M A^-1 M^T N^-1 and H = M A^-1 M^T N^-1,
the profile is chi2_min = d^T P d and its score is s_r = -(1/2) d^T P_r d with P_r = dP / dbeta_r.
For d = M a + n: P M = 0 gives P_r M a = -P M_r a, so the mean of the score is zero and its
signal part has variance (M_r a)^T P N P (M_r' a) = (M_r a)^T P (M_r' a) = J_rr'. The noise part is
(1/2) tr(P_r N P_r' N) = (1/2) tr(H_r H_r'), and with H_r = (1 - H) M_r G + M A^-1 M_r^T N^-1 (1 - H),
G = A^-1 M^T N^-1, the two squared terms vanish (G (1 - H) = 0 and (1 - H) M = 0) and the cross
terms give T_rr' = tr(A^-1 M_r^T P M_r'). The expected curvature of chi2_min / 2 is J. For one
amplitude per element and diagonal noise this is J = sum u V and T = sum V, the form checked in
the E036 addendum.

## Setup

- Unit tests (`nextgen/lowres/tests/test_likelihood.py`, five tests, pass at 747ed50 plus the
  module): analytic gradients of the three objectives against finite differences in four
  configurations; equality with the per-element closed forms for diagonal noise; J and T against
  sum u V and sum V; the sandwich covariance against the scatter of 300 fits with dense noise,
  three components and low S/N; the reference-frequency identity of the flat-prior marginal.
- Simulation: truth synchrotron amplitudes from the Cosmoglobe DR1 polarized synchrotron map (5
  deg, scaled from 30 to 23 GHz with each region's true index), dust amplitudes from the PR4 353
  GHz map at the working resolution, CMB drawn N(0, 0.2^2). Noise drawn from the real covariances.
  Bands K, Ka, Q, V, 353.
- Region sets and true indices, chosen away from the centre of the allowed range: the four
  BeyondPlanck XV regions with -3.15 (High Latitude), -3.25 (Spur), -2.85 (Galactic Center), -2.95
  (Plane); and the 48 nside-2 pixels with -3.1 + 0.12 g, g a fixed standard normal draw.
- Bulk at nside 8 (data and covariances degraded as the real nside-8 analysis would be), a spot
  check at nside 16. Profile fits for every realization; flat-prior marginal at references of 23
  and 30 GHz and the zero-mean prior for a subset.

## Reference

The truth is known. z = (estimate - truth) / sd per region, with the naive sd from J^-1 and with
the sandwich sd.

**Predictions, written before the first simulation run:**

- V1. Profile, unbiased: mean z (sandwich) within +-0.2 in every BP XV region, and within +-0.3
  for the nside-2 regions whose sandwich sd is below 0.25, with 200 realizations.
- V2. Profile, width: sd of z with the sandwich sd between 0.9 and 1.1; with the naive sd between
  1.0 and 1.5 and never below the sandwich value by more than 0.03.
- V3. Coverage of +-1 sandwich sd between 0.64 and 0.72 over the regions with sd below 0.25.
- V4. Flat-prior marginal: optimum minus profile optimum has opposite signs at references of 23
  and 30 GHz or, if the pivot lies below 23 GHz as with a free CMB in E036, is larger at 30 than
  at 23 by ln(30 / 23) times the number of amplitudes times the squared profile sd, within 25
  percent, in the regions where that product is below 0.3.
- V5. Zero-mean prior with tau maximized: within 0.3 sandwich sd of the profile optimum in the BP
  XV regions.
- V6. The reduced chi-square at the profile optimum is 1.00 +- 0.03 (the model is the truth).

**Addendum planned 2026-11-10 21:32, expectations written before it is run: the errors one data set gives.**
V2 and V3 use J and T at the truth, so they test the theory and not what is done on real maps.
While the simulations ran I found that J evaluated at the fitted amplitudes has expectation J + T
(E041, F46), and `fit.sandwich_plugin` now uses J_c = J(a_hat) - T. The addendum repeats the first
40 realizations of the two nside-8 sets with errors from each realization
(`nextgen/scripts/e039b_plugin_check.py`).

- A1. Four BP XV regions (T / J at the truth below 0.02): the curvature, plain plug-in and
  debiased widths agree within 2 percent, and the sd of z is between 0.88 and 1.12 (160 values).
- A2. 48 regions, those with a sandwich sd at the truth below 0.25: with the debiased sandwich the
  pooled sd of z is between 0.92 and 1.08 and the 68 percent coverage between 0.64 and 0.72; in
  the classes with T / J above 0.1 the curvature width gives an sd of z larger than the debiased
  one by at least 0.03, and the plain plug-in a larger one still.
- A3. Fewer than 5 percent of the measured region fits have a J_c that is not positive.

## Command

```bash
cd /Users/trygvels/code/Commander4
.venv/bin/python -m pytest -q -p no:cacheprovider nextgen/lowres/tests
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions bpxv --nreal 200 --full 40 \
  > nextgen/experiments/E039/sim_n8_bpxv.txt
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions hp2 --nreal 200 --full 0 \
  > nextgen/experiments/E039/sim_n8_hp2.txt
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 16 --regions bpxv --nreal 12 --full 0 \
  > nextgen/experiments/E039/sim_n16_bpxv.txt
```

## Result

**Stage 1, the four BP XV regions, read 2026-11-10 22:03** (`E039/sim_n8_bpxv.txt`, 200 realizations at nside 8,
the first 40 with every estimator, 20:34 to 21:55; `E039/sim_n16_bpxv.txt`, 12 realizations at nside
16, 20:34 to 21:18). The 48-region set is still running and is read separately. Sandwich sd at
the truth: 0.0153, 0.0166, 0.0108, 0.0146 (nside 8) and 0.0145, 0.0158, 0.0104, 0.0137 (nside 16);
T / J at the truth at most 0.006 (nside 8) and 0.018 (nside 16).

| | High Latitude | Spur | Galactic Center | Plane | pooled |
|---|---|---|---|---|---|
| nside 8: mean shift | -0.0003 | 0.0000 | +0.0007 | +0.0002 | |
| mean z (sandwich) | -0.02 +- 0.07 | 0.00 +- 0.07 | +0.07 +- 0.07 | +0.01 +- 0.07 | +0.014 |
| sd of z, sandwich / naive | 1.036 / 1.039 | 0.979 / 0.979 | 0.936 / 0.936 | 1.024 / 1.024 | 0.993 / 0.994 |
| coverage of +-1 sandwich sd | 0.670 | 0.715 | 0.685 | 0.695 | 0.691 |
| flat23 - profile (40 realizations) | -0.0008 | 0.0000 | -0.0002 | -0.0004 | |
| flat30 - flat23, measured | +0.0687 | +0.0270 | +0.0076 | +0.0142 | |
| predicted, scalar ln(30/23) n sd^2 | +0.0692 | +0.0156 | +0.0015 | +0.0085 | |
| predicted, coupled J^-1 (n ln(30/23)) | +0.0723 | +0.0278 | +0.0078 | +0.0145 | |
| flat30 - profile in sandwich sd | +4.45 | +1.62 | +0.69 | +0.94 | |
| zero-mean - profile, mean (rms) in sd | +0.03 (0.03) | 0.00 (0.01) | +0.06 (0.06) | +0.04 (0.04) | |
| nside 16 (12 realizations): mean z | +0.30 +- 0.30 | +0.09 +- 0.25 | -0.54 +- 0.27 | +0.28 +- 0.28 | +0.03 |
| sd of z, sandwich | 1.035 | 0.855 | 0.951 | 0.980 | 0.989 |

Reduced chi-square at the optimum: 1.0006 (sd 0.022 over 200) at nside 8, 0.9977 (sd 0.0091 over 12)
at nside 16.

- **V1 pass** for the four regions: mean z within 0.07 of zero where 0.2 was allowed.
- **V2: the sandwich clause passes** (0.936 to 1.036). **The naive clause fails as I worded it in two
  regions**: I wrote "between 1.0 and 1.5", and the Spur and Galactic Center give 0.979 and 0.936.
  With T / J below 0.006 the two widths are the same number, and the lower bound of 1.0 left no
  room for the sampling error of 0.05. Nothing is learned about the sandwich correction from this
  set; that is what E041 and the addendum are for.
- **V3 pass**: coverage 0.670 to 0.715, pooled 0.691.
- **V4: fails as first worded in three of four regions, passes in the coupled form in all four.**
  The scalar formula holds only for the High Latitude region (+0.0687 against +0.0692). The Spur,
  Galactic Center and Plane move by 1.7, 5 and 1.7 times their scalar predictions, because the
  regions are coupled through the dense covariances and the shared amplitudes' noise, and the
  large region drags its neighbours. The coupled form, adopted after a two-realization smoke test
  and before this run, is within 5 percent in every region (ratios 0.950, 0.971, 0.974, 0.979).
  The pivot of this band set with a CMB prior sits at 23 GHz: the flat prior stated there is
  neutral to 0.05 sd, and stated at 30 GHz it moves the High Latitude index by 4.5 sandwich sd.
  The independent review (below) corrects the formula further: the response uses the Hessian of
  the marginal, J + Hess (1/2) ln det A, not J. Here T / J is below 0.006 and the two agree.
- **V5 pass**: the zero-mean prior is within 0.06 sd of the profile optimum in every region.
- **V6 pass**: 1.0006 and 0.9977.

**Independent review** (`../reviews/2026-11-10_stats_review_lowres.md`, stats-reviewer agent, own
derivations written before reading this record; scripts in `E039/stats_review_scripts/`). Verdict:
correct with caveats. Objectives, gradients and the LFI noise model agree with brute force to
1e-13. What it changed:

1. *The first-order shift of the flat-prior optimum uses the Hessian of the marginal*, J + T + S
   with T + S the Hessian of (1/2) ln det A (2T for one component), not J. I checked this on the
   E041 setup with a reference close to the pivot: the measured shift over g / (J + 2T) is 1.008
   and 1.031 at T / J of 0.03 and 0.11, where over g / J it is 0.955 and 0.842. The ratios near one
   that E041 reported against g / J were in cells with T / J below 0.06. `fit.marginal_hessian`
   added; Paper A1 Sect. 4 corrected.
2. *With a CMB prior the profile score is unbiased over CMB skies drawn from the prior, not for
   one fixed sky.* If the true CMB variance is kappa^2 times the prior variance the optimum is
   biased by -(kappa^2 - 1) J^-1 t to first order, growing as the square root of the number of
   elements in units of the sd. E039 draws the CMB from the prior and cannot see this.
   `information` now returns t; E040 prints the sensitivity per region, and its flat, 0.1 and 0.4
   uK arms measure it on the data.
3. *`wmap_band`: code and docstring disagreed on the detector weighting.* The code was right: with
   each detector at unit thermal response and equal weights the files reproduce the per-radiometer
   effective frequencies of the WMAP team's `effective_freq.pro` to 0.005 GHz (K11 22.055 against
   22.055, K12 22.870 against 22.875, and so on through V22), and the team averages radiometers
   with equal weights. Summing the peak-normalized curves misses the K-band value by 0.15 GHz.
   The same routine applies a nine-year drift correction (0.13, 0.12, 0.11, 0.06 percent for K to
   V) that I had not applied: 0.002 on an index within WMAP, 0.018 between K band and 30 GHz. It
   is now the default, with a no-drift arm in E040.
4. The first-order sandwich sd is 3 percent low on average (up to 7 percent at T / J near 1), and
   the debiased plug-in is unreliable above T / J of 1.5 even when J_c is positive. E040 quotes a
   sandwich error only where T / J_c is below 1.
5. Two-sided coverage of the plug-in interval is nominal; one-sided is not at 48 elements (the
   quoted error is anticorrelated with the error itself).

**The four-region configuration is cleared for the real maps.** The 48-region configuration is
not, until its simulation and the plug-in addendum are read.

### Addendum A1 result (plug-in errors, four regions), read 2026-11-11 02:46

`E039/plugin_n8_bpxv.txt`, 40 realizations, 160 region-realizations, T / J at the truth below 0.02 in all
four regions: sd of z 0.974 (curvature), 0.974 (plain plug-in sandwich), 0.973 (debiased sandwich);
coverage of one sd 0.700 for all three; mean z of the debiased version +0.015 +- 0.077. **A1 holds**: at
this signal-to-noise the three error estimates are the same and all are calibrated. The run is one of
those exposed to the Cholesky overflow of E043 and is rerun there (V4).

### Addendum B, pre-registered 2026-11-11 03:20: the estimator that reads the index from large scales only

Referee item 2 on Paper A1 draft 1: the low-pass arm of E042 was run on the maps without a simulation
test. `nextgen/scripts/e039c_lowpass_sim.py`: the E039 nside-8 simulations (same truth, noise and seeds;
current band model), K band and 353 GHz at all scales, Ka, Q and V through their harmonic coefficients at
multipoles 2 to 12 with their covariance (330 of 1144 modes), profile optimum and debiased sandwich error,
100 realizations, next to the all-scales fit of the same realizations.

- B1. Unbiased: mean z within 2 standard errors of zero in each region (se 0.1), pooled within 0.1.
- B2. Calibrated: sd of z within 0.85 to 1.15 per region, coverage of one sd within 0.58 to 0.78.
- B3. Less precise than the all-scales fit, by the information removed: I expect the error to grow by 10
  to 40 percent (the real-map fits went from 0.012 to 0.018 at High Latitude between all scales and
  multipole 12).
- B4. The difference between the two estimates of one realization has zero mean and a spread consistent
  with sqrt(sd_low^2 - sd_all^2) (the all-scales fit uses a superset of the data and is efficient).
- If B1 or B2 fails the estimator is not used on the maps anywhere in Paper A1 or B until understood.
- What this cannot show: that the estimator removes the pixel-scale bias of the real maps. The simulations
  have none (signal smooth at 5 degrees, no beam mismatch). It shows only that the estimator is unbiased
  and calibrated when the model is true.

### Addendum B result at nside 8 (read 2026-11-11 06:56)

`E039/lowpass_n8_l12.txt`: 100 realizations, sound Cholesky path, canary (no OVERFLOW line). Regions in the
usual order.

| Fit | Mean estimate minus truth | Mean z (se 0.10) | sd of z | Coverage of 1 sd | Median quoted error | Scatter of the estimate |
|---|---|---|---|---|---|---|
| all scales | -0.0018, -0.0002, +0.0005, +0.0010 | -0.11, -0.01, +0.05, +0.08 | 1.05, 1.04, 0.93, 0.98 | 0.65, 0.67, 0.69, 0.76 | 0.0152, 0.0166, 0.0107, 0.0145 | 0.0160, 0.0173, 0.0099, 0.0142 |
| index from multipoles up to 12 | -0.0018, -0.0008, +0.0004, +0.0007 | -0.09, -0.04, +0.04, +0.05 | 0.95, 1.03, 0.95, 0.99 | 0.75, 0.66, 0.72, 0.67 | 0.0184, 0.0190, 0.0123, 0.0154 | 0.0176, 0.0195, 0.0118, 0.0152 |

Low-pass minus all scales, same realization: mean -0.0000, -0.0006, -0.0001, -0.0003; sd 0.0092, 0.0097,
0.0057, 0.0058, against sqrt(sd_low^2 - sd_all^2) = 0.0104, 0.0092, 0.0061, 0.0052.

- B1 holds: mean z within one standard error of zero in every region.
- B2 holds: sd of z 0.95 to 1.03, coverage 0.66 to 0.75.
- B3 holds for three regions: the error grows by 21, 14 and 15 percent; in the Plane by 6 percent (written:
  10 to 40).
- B4 holds: the paired difference has zero mean and a spread of 0.89, 1.05, 0.94 and 1.12 of
  sqrt(sd_low^2 - sd_all^2). So the difference between the low-pass and the all-scales fit of one data set can
  be judged against that number.
- The independent review of the construction (`reviews/2026-11-11_stats_review_prior_block_lowpass.md`)
  confirms it algebraically (identities to 1e-13) and on 400 realizations at nside 4, and measured that
  contamination of a low-passed band outside the kept multipoles moves the estimate by 1e-14 of its error.
- What this gives for the maps (E042): with the LCDM prior the index from multipoles up to 20 is steeper
  than the all-scales index by 0.024, 0.031, 0.055 and 0.084, which is 2.8, 3.3, 8.0 and 10.3 times
  sqrt(sd_low^2 - sd_all^2). Under the model that difference has zero mean, so the scale dependence is a
  departure from the model in every region, not a property of the estimator. Its cause stays open (Q25).
- The spot check at nside 16 with the LCDM prior is part of `e040h_sim_n16.py` (queue 12).

## Negative control

The flat-prior marginal (V4) is the known-biased target.

## Conclusion

Pending.

## Follow-ups
