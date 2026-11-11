# E042: anatomy of the disagreement between WMAP K band and Planck 30 GHz in polarization

- **Date started / finished:** 2026-11-10 22:56 (pre-registration) /
- **Branch and commit:** `nextgen` at c8351bd
- **Open question addressed:** Q21; the step named in advance in E040 stage 2 ("if S1 holds ...")
- **Validation level:** L5 (real data); the estimator is the one validated in E039 and E041

## Question

E040 stage 2, read at 22:55 on 2026-11-10: with K band and PR4 30 GHz alone the profile index of
the High Latitude region is -3.428 +- 0.019, against -3.062 +- 0.016 with K and Ka and -3.094 +-
0.012 with K, Ka, Q and V (Spur -3.180 against -3.074; Galactic Center -3.263 against -3.046;
Plane -3.035 against -2.980). S1 held. In amplitude the High Latitude gap is 8 percent between
maps 0.22 apart in ln nu. Weiland et al. (2022) report the same sign in the Fan region with PR3
(-3.27 against -2.90) and attribute it to large-scale systematics, mostly at 30 GHz; Watts et al.
(2023) attribute the K-band side to poorly measured modes that the released covariance does not
fully describe. Which is it, on these maps, with these covariances: a number that multiplies one
map, an angle, or structure at particular angular scales in one map, and in which one?

What changes: this is Sect. 7.5 of Paper A1 and decides whether 30 GHz can enter a combined index
at all. If the gap lives in a few large-scale modes of one map, that map's covariance is what
needs repair, and the index from the other bands stands. If it is a multiplicative factor, the
calibration or the bandpass is the suspect. If nothing removes it, it is reported as open.

## Target and priors

The profile likelihood of E040 (four BP XV regions, nside 16, window-matched, dense WMAP
covariances, LFI noise model of E038, CMB prior 0.2 uK, dust from 353 GHz), with one of these
changes at a time:

- (a) *Scale cut on 30 GHz*: the E and B modes of the 30 GHz map at multipoles 2 to L are
  marginalized by adding a variance of 1e6 times the white level to those modes in the whitened
  noise model, for L = 3, 5, 8, 12, 20.
- (b) *Scale cut on K band*: the same modes projected out of K band, N^-1 -> N^-1 - N^-1 V (V^T N^-1
  V)^-1 V^T N^-1 with V the harmonics at multipoles 2 to L, same L.
- (c) *Worst K-band modes*: the eigenmodes of the released K-band covariance whose variance is
  above 10 times the median are projected out in the same way.
- (d) *Amplitude factor*: bands K, Ka, 30, 353 with one free factor c multiplying the model at 30
  GHz, per region (four values) and global (one value), profiled together with the indices.
- (e) *Angle*: bands K, Ka, 30, 353 with the 30 GHz map rotated by a free polarization angle psi,
  global.
- (f) *Other maps for 30 GHz*: PR3 (2018) with its released low-resolution noise covariance if it
  can be put in the same form; Cosmoglobe DR1 K and 30 GHz with their white-noise weights.

## Setup

Script `nextgen/scripts/e042_k30_anatomy.py`, new file, importing `nextgen/lowres` read-only. One
dense job at a time behind the E040 queue.

## Reference

No truth. The reference is the K, Ka, 353 index of each region (E040) and the internal consistency
of the three descriptions.

**Expectations, written before any of these fits.**

- X1. The gap is a large-scale effect: with multipoles 2 to 8 of the 30 GHz map marginalized, the
  High Latitude index from K, 30, 353 moves at least half of the way from -3.43 to -3.06, and its
  error grows by less than a factor 2.5.
- X2. The same cut on K band moves it by less than half as much as the cut on 30 GHz.
- X3. Projecting out the worst K-band modes (c) moves the High Latitude index by less than 0.1.
- X4. The amplitude factor is not one number: the per-region factors differ from each other by
  more than 3 combined sd, with the High Latitude factor between 0.88 and 0.95 and the Plane factor
  between 0.97 and 1.01.
- X5. The angle fit prefers abs(psi) below 3 degrees and does not remove the gap (the index from
  K, 30 stays more than 0.2 steeper than from K, Ka at high latitude when psi is fixed at its
  fitted value).
- X6. No prediction for (f).
- If X1 fails (the gap survives the removal of multipoles 2 to 20 at 30 GHz), the large-scale
  reading is wrong and the amplitude and bandpass descriptions come first.

**Added 2026-11-10 23:00, before any E042 fit: (g) prediction residual.** After the pre-registration E040
returned the band set without K band (30, 44, 70, Ka, Q, V, 353): High Latitude -2.994 +- 0.019,
Spur -3.059 +- 0.024, Galactic Center -3.048 +- 0.012, Plane -2.766 +- 0.018. So 30 GHz is faint
against K band (steep K to 30 index) and faint against Ka to V (flat index when 30 GHz is the
lowest band). A direct look: fit K, Ka, Q, V, 353 without any Planck LFI band, predict the 30 GHz
map from the fitted amplitudes and indices, and study the difference D = observed - predicted
(`nextgen/scripts/e042b_predict30.py`): its regression on the predicted synchrotron per region
(a fractional deficit), its chi-square under the 30 GHz noise model, its E and B power by
multipole bin against that model, and its map.

- X7. The regression slope of D on the predicted synchrotron is negative in all four regions
  (30 GHz fainter than predicted), between -0.03 and -0.10 at high latitude and between 0 and
  -0.03 in the Plane.
- X8. After removing the per-region slope the residual has excess power over the noise model at
  multipoles 2 to 6 (ratio above 2) and is consistent with it above multipole 20 (ratio between
  0.8 and 1.3).
- X9. The same residual for Ka predicted from K, Q, V, 353 (the control) has per-region slopes
  within 0.02 of zero.

**(f) made concrete, 2026-11-11 00:30, before any fit with a PR3 or DR1 map.** Three arms.

- *PR3 with a fitted noise model* (`e042c_pr3_arm.py`): the PR3 30 GHz map through the same steps
  as PR4, its noise model fitted to the PR3 ring-half difference.
- *PR3 with the released covariance* (`e042d_pr3_ncvm.py`): the released nside-16 map and
  `offset_covmat` (IRSA, `LFI_NoiseCovMat_030_0016_R3.20.tgz`, fetched 2026-11-11 00:20). Checked
  before any fit: the matrix is symmetric, positive definite with a floor at the 0.02 uK
  regularization, pixel sd 0.85 uK in Q and U; the released map correlates at 0.93 with my
  noise-weighted downgrade of the PR3 map and its transfer function against it is the cosine
  window (1.00 to multipole 10, 0.94 at 20, 0.83 at 24, 0.47 at 32) times about 0.97. The 30 GHz
  information is restricted to multipoles 2 to 24; the PR4 map is cut to the same multipoles for
  a like-for-like comparison.
  - Y1. The PR3 ring-half difference through the same steps has a chi-square per mode between
    0.8 and 1.3 under the released covariance. Outside that, the covariance is rescaled by the
    measured value for the fits and the fact is reported.
- *Cosmoglobe DR1* (`e040c_dr1_arm.py`): the DR1 K, Ka and 30 GHz maps with their white-noise
  weights, which understate the errors.
- Y2. No prediction for the indices of these arms. What they are for: if the K to 30 GHz index at
  high latitude is within 0.1 of the K to Ka index for one of these map versions and not for PR4,
  the gap belongs to the processing of the 30 GHz map (or, for DR1, to the joint reprocessing of
  both), and the paper says which.

## Command

Pending.

## Result

**(g) Prediction residual, 2026-11-10 23:37 to 23:43** (`E042/predict30.txt`, `predictKa.txt`; a first run
crashed after the slopes on a formatting error of mine and is kept as `predict30_first_run_crashed.txt`;
the numbers it printed are identical).

| | 30 GHz (PR4) predicted from K, Ka, Q, V, 353 | Ka predicted from K, Q, V, 353 (control) |
|---|---|---|
| rms of the map, of the prediction, of D (uK) | 16.05, 17.48, 3.74 | 11.54, 11.06, 2.42 |
| chi-square of D per element, target noise only | 4.29 | 1.23 |
| slope of D on the predicted synchrotron: High Latitude | -0.114 +- 0.003 | +0.019 +- 0.006 |
| Spur | -0.024 +- 0.004 | +0.014 +- 0.007 |
| Galactic Center | -0.077 +- 0.002 | +0.037 +- 0.004 |
| Plane | -0.024 +- 0.004 | +0.010 +- 0.006 |
| one slope for the sky | -0.073 +- 0.002 | +0.024 +- 0.003 |
| chi-square per element after removing four slopes | 3.80 | 1.21 |
| power of whitened D over white, multipoles 2 to 3, E / B, slopes removed | 33 / 43 | 3.7 / 2.2 |
| 4 to 6 | 16 / 19 | 3.5 / 2.3 |
| 7 to 10 | 7.3 / 12 | 1.6 / 1.3 |
| 11 to 20 | 5.9 / 5.5 | 1.3 / 1.2 |
| 21 to 29 | 3.9 / 4.1 | 1.1 / 1.1 |
| 30 to 47 | 2.8 / 3.1 | 1.2 / 1.2 |

The errors on the slopes and the power ratios count the target's noise only; the noise of the
prediction is not in them, so the ratios are upper limits on the excess and the control column
is the yardstick.

- **X7: the sign holds in all four regions, the High Latitude size is beyond my range.** PR4 30 GHz
  is fainter than the WMAP-only prediction by 11.4 percent at high latitude (I allowed 3 to 10),
  7.7 percent in the Galactic Center and 2.4 percent in the Spur and the Plane.
- **X8 holds at large scales and fails at small ones as I wrote it.** After the slopes are removed
  the residual has 16 to 43 times the white expectation at multipoles 2 to 6. Above multipole 20
  it is at 2.8 to 4.1, not the 0.8 to 1.3 I wrote: I left the noise of the prediction out of the
  expectation. Against the control (1.1 to 1.2 there) the 30 GHz residual is still 2.5 to 3.5
  times too large at the pixel scale and 5 to 20 times at multipoles 2 to 10.
- **X9 holds in three regions**: the Ka control slopes are +0.010 to +0.019, and +0.037 in the
  Galactic Center. Ka is 1 to 2 percent brighter than K, Q and V predict, the K, Ka against four
  bands difference of E040 seen from the other side.
- B-mode power in the 30 GHz residual is as large as E-mode power at every scale. Synchrotron
  has less B than E; a residual with equal E and B is not sky signal misassigned by a wrong
  index, it is map-level structure.
- So the gap is two things at once: a deficit of 30 GHz that depends on the region (2 to 11
  percent), and large-scale structure in both E and B that no amplitude factor describes.
- Caveat on the slopes: the regressor is the fitted synchrotron, which carries noise, so each
  slope has an errors-in-variables bias of minus the noise share of the regressor's variance,
  about -1.5 percent at high latitude and less elsewhere. The control shares it. The joint fits
  below do not have it.

**(a) and (b) Scale cuts, 2026-11-11 00:16 to 01:01** (`E042/base_cuts.txt`; K, 30, 353; the modes of one
map at multipoles 2 to L projected out; profile index with debiased sandwich error).

| cut | modes | High Latitude | Spur | Galactic Center | Plane |
|---|---|---|---|---|---|
| none | 0 | -3.428 +- 0.019 | -3.180 +- 0.023 | -3.263 +- 0.013 | -3.035 +- 0.020 |
| 30 GHz, L = 3 | 24 | -3.329 +- 0.019 | -3.055 +- 0.025 | -3.202 +- 0.013 | -2.974 +- 0.021 |
| 30 GHz, L = 5 | 64 | -3.318 +- 0.020 | -3.026 +- 0.028 | -3.194 +- 0.013 | -2.913 +- 0.023 |
| 30 GHz, L = 8 | 154 | -3.270 +- 0.022 | -2.976 +- 0.032 | -3.154 +- 0.014 | -2.767 +- 0.026 |
| 30 GHz, L = 12 | 330 | -3.224 +- 0.024 | -2.966 +- 0.037 | -3.103 +- 0.015 | -2.703 +- 0.030 |
| 30 GHz, L = 20 | 874 | -3.129 +- 0.031 | -2.903 +- 0.044 | -3.043 +- 0.018 | -2.683 +- 0.035 |
| K band, L = 3 | 24 | -3.344 +- 0.020 | -3.071 +- 0.026 | -3.195 +- 0.013 | -2.964 +- 0.022 |
| K band, L = 5 | 64 | -3.334 +- 0.020 | -3.053 +- 0.029 | -3.188 +- 0.013 | -2.895 +- 0.024 |
| K band, L = 8 | 154 | -3.291 +- 0.022 | -3.007 +- 0.033 | -3.150 +- 0.014 | -2.734 +- 0.027 |
| K band, L = 12 | 330 | -3.247 +- 0.025 | -3.004 +- 0.038 | -3.099 +- 0.015 | -2.667 +- 0.030 |
| K band, L = 20 | 874 | -3.159 +- 0.032 | -2.960 +- 0.046 | -3.040 +- 0.018 | -2.644 +- 0.035 |
| for reference: K, Ka, 353 | | -3.062 +- 0.016 | -3.074 +- 0.020 | -3.046 +- 0.011 | -2.980 +- 0.016 |
| K, Ka, 30, 353 | | -3.215 +- 0.013 | -3.111 +- 0.016 | -3.130 +- 0.009 | -3.004 +- 0.013 |

The reduced chi-square of the K, 30 fit falls from 2.64 (no cut) to 1.79 (L = 20), and the
30 GHz chi-square per element from 1.85 to 1.19.

- **X1 fails narrowly, and the expectation behind it was too simple.** With multipoles 2 to 8 of
  30 GHz marginalized the High Latitude index moves 43 percent of the way to the K, Ka value
  (-3.270), where I required half; at L = 12 it is 56 percent, at L = 20, 82 percent. The error
  grows by 1.16 (L = 8), inside the 2.5 allowed. But the index does not converge on the K, Ka
  value and stop: in every region it keeps flattening as large scales are removed, and in the
  Spur and the Plane it passes the K, Ka value at L = 3 to 5 and ends 0.17 and 0.30 flatter. **The
  index between K band and 30 GHz depends on angular scale**: steepest at the largest scales,
  flattest at the smallest, by 0.30 (High Latitude), 0.28 (Spur), 0.22 (Galactic Center) and
  0.35 (Plane) between no cut and L = 20. In amplitude the ratio of the 30 GHz map to the K-band
  map rises by 5 to 8 percent from the largest scales to multipoles above 20.
- **X2 was ill-posed.** In a fit with two synchrotron bands a mode informs the index only if both
  bands measure it, so projecting it out of either band removes the same information: the K-band
  cuts reproduce the 30 GHz cuts within 0.03 at every L. The comparison cannot say which map
  carries the scale dependence. The control with Ka (below) can.

**Exploratory, not pre-registered, prompted by the residual map** (`E042/leakage_check.txt`,
`nextgen/papers/paperA1/figures/fig_resid30.pdf`). The map of D shows a band along the inner
Galactic plane and a blob at the Galactic Center. Regression of D on the PR4 30 GHz intensity at
nside 16:

| zone | Q: slope on I, correlation | U: slope on I, correlation |
|---|---|---|
| abs(b) < 5, inner Galaxy (108 pixels) | -0.13 percent, -0.54 | +0.10 percent, +0.63 |
| abs(b) < 5, outer Galaxy (212) | +0.12 percent, +0.67 | +0.05 percent, +0.36 |
| 5 < abs(b) < 20 (768) | +0.03 percent, +0.20 | +0.03 percent, +0.18 |
| abs(b) > 30 (1600) | +0.01 percent, +0.08 | 0.00 percent, 0.00 |

After the per-region amplitude factors are removed the U residual in the plane keeps its
correlation with intensity (+0.58) and has none with the predicted synchrotron (+0.03). So in
the plane the residual follows the intensity at the 0.05 to 0.13 percent level with a sign that
changes between the inner and outer Galaxy in Q: what a leakage of intensity into polarization
does (the intensity reaches 22 mK in the brightest pixel, so 0.1 percent is 22 uK there). Which
map carries it is not decided by this; the Ka control, whose prediction is built mostly on K
band, shows only a small feature at the Galactic Center, which points at 30 GHz. At abs(b) > 30
there is no correlation with intensity, and the residual there anticorrelates with the predicted
synchrotron (-0.37): the high-latitude deficit is a different thing from the plane residual.

**Negative control, 2026-11-11 01:01 to 01:17** (`E042/modes_ctrl.txt`): the same cuts in the K, Ka, 353 fit.

| cut | High Latitude | Spur | Galactic Center | Plane |
|---|---|---|---|---|
| none | -3.062 +- 0.016 | -3.074 +- 0.020 | -3.046 +- 0.011 | -2.980 +- 0.016 |
| Ka, L = 8 | -3.048 +- 0.019 | -2.992 +- 0.033 | -2.992 +- 0.013 | -2.754 +- 0.028 |
| K band, L = 8 | -3.051 +- 0.019 | -2.994 +- 0.033 | -2.993 +- 0.013 | -2.735 +- 0.029 |
| Ka, L = 20 | -2.979 +- 0.030 | -2.912 +- 0.046 | -2.937 +- 0.017 | -2.615 +- 0.037 |
| K band, L = 20 | -2.989 +- 0.030 | -2.916 +- 0.048 | -2.936 +- 0.018 | -2.599 +- 0.039 |

- **The negative control fails.** I required the K, Ka index to move by less than 2 sd of its
  enlarged error. It moves by +0.014, +0.082, +0.054 and +0.226 at L = 8 (0.7, 2.5, 4.1 and 8.2
  sd) and by +0.083, +0.162, +0.109 and +0.365 at L = 20. Within WMAP too, the index flattens as
  large scales are removed.
- So the scale dependence of (a) is not a property of 30 GHz. In amplitude, the ratio to K band at
  multipoles above 20 against all scales rises by 3, 6, 4 and 14 percent for Ka and by 7, 6, 5
  and 8 percent for PR4 30 GHz (the four regions in order): K band has less power at the pixel
  scale than the other bands, relative to what the window model gives it.
- Mechanism, consistent with the last bin of the ratio table (every ratio to K band up by 6 to 10
  percent at multipoles 35 to 47): the nside-16 maps are pixel averages of maps that are not
  band-limited. Structure above multipole 47, of which the bright plane has plenty, aliases into
  the pixel scale with each band's own beam, and K band has the widest. The window operator works
  to multipole 47 and cannot touch this. The E038 (c) simulation did not show it at this size
  because its sky had no compact structure and its fits used all scales.
- **Consequence for (a):** the statement "the index between K band and 30 GHz depends on angular
  scale" is withdrawn as a statement about 30 GHz. Net of the control, only the High Latitude
  region keeps a 30 GHz-specific large-scale effect (about 4 percent in amplitude between all
  scales and multipoles above 20).
- **Consequence for E040:** every fit that uses the pixel scale carries this bias toward a flatter
  index, largest where the sky is bright and structured. It has to be measured from the other
  side, with fits that do not use the pixel scale.

**Low-pass arm, expectations written 2026-11-11 01:17 before it is run** (`nextgen/scripts/e042e_lowpass.py`).
Every band restricted to multipoles 2 to L by keeping its harmonic coefficients with their
covariance; L = 20 first, then 12 and 30. Band sets K, Ka, Q, V; K, Ka; K, 30; 30, Ka (all with
353).

- Z1. At L = 20 the K, Ka, Q, V indices are steeper than the all-scale baseline by 0 to 0.04 at
  high latitude and by 0.02 to 0.12 in the Plane and the Spur.
- Z2. At L = 20 the K, 30 (PR4) index is between -3.7 and -3.3 in the Plane and between -3.6 and
  -3.3 at high latitude: with the pixel scale gone, the plane deficit of PR4 that the ratio table
  shows (0.44 to 0.47 against 0.505) is no longer masked.
- Z3. At L = 20 the 30, Ka index is flatter than -2.6 in the Plane.
- Z4. The indices at L = 12, 20 and 30 for K, Ka, Q, V agree with each other within 2 sd of the
  L = 12 error; if they do, L = 20 becomes the baseline of the paper's tables and the all-scale
  fit is reported as the arm that carries the aliasing bias. If they do not, there is no
  resolution-independent index for these regions at this depth and the paper says so.

**Exploratory, not pre-registered: amplitude ratios by multipole and PR4 against PR3, 2026-11-11
01:02 to 01:03** (`E042/cross_ratio.txt`, `E042/pr3_vs_pr4_direct.txt`). Prompted by the scale
dependence of (a). Two model-light diagnostics on the window-matched nside-16 maps.

(i) The amplitude ratio of two maps per multipole bin from cross spectra with a third map, which
has no noise bias: R(x / K; ref) = sum c_x c_ref / sum c_K c_ref over the E and B coefficients.
For an index of -3.06 the expected ratios are 0.505 (30 GHz over K) and 0.321 (Ka over K).

| multipoles | PR4 30 / K, full sky | PR3 30 / K, full sky | Ka / K, full sky | PR4 30 / K, abs(b) > 20 | PR3 30 / K, abs(b) > 20 | PR4 30 / K, abs(b) < 20 | PR3 30 / K, abs(b) < 20 |
|---|---|---|---|---|---|---|---|
| 2 to 3 | 0.448 | 0.464 | 0.326 | 0.447 | 0.420 | 0.442 | 0.488 |
| 4 to 6 | 0.450 | 0.493 | 0.321 | 0.483 | 0.482 | 0.444 | 0.489 |
| 7 to 10 | 0.446 | 0.490 | 0.322 | 0.486 | 0.474 | 0.440 | 0.485 |
| 11 to 16 | 0.454 | 0.492 | 0.326 | 0.500 | 0.495 | 0.455 | 0.493 |
| 17 to 24 | 0.462 | 0.504 | 0.332 | 0.508 | 0.504 | 0.461 | 0.504 |
| 25 to 34 | 0.470 | 0.507 | 0.331 | 0.505 | 0.492 | 0.469 | 0.508 |
| 35 to 47 | 0.497 | 0.535 | 0.353 | 0.524 | 0.542 | 0.497 | 0.530 |

(ii) Direct regression of the PR4 map on the PR3 map, and of each on K band, by zone (Q and U
together, through the origin, unweighted):

| zone | PR4 on PR3 | PR4 on K | PR3 on K | correlation PR4, PR3 |
|---|---|---|---|---|
| abs(b) < 5, inner Galaxy | 0.854 | 0.440 | 0.503 | 0.978 |
| abs(b) < 5, outer Galaxy | 0.841 (Q 0.92, U 0.44) | 0.499 | 0.498 | 0.872 |
| 5 to 20 | 0.955 | 0.486 | 0.484 | 0.968 |
| 20 to 45 | 0.988 | 0.471 | 0.460 | 0.983 |
| above 45 | 1.034 | 0.426 | 0.383 | 0.969 |
| all sky | 0.886 | 0.456 | 0.490 | 0.966 |

(The slopes on K band at high latitude are attenuated by the noise of K band; the cross-spectrum
ratios of (i) are not.)

What these say, as observations to be checked and not yet as results:

- **Two different things are in the gap, and only the first is secure.** With bootstrap errors over
  the modes of each bin (`E042/cross_ratio_errors.txt`, added 01:04 after the first reading):
  (A) for abs(b) < 20 the PR4 30 GHz over K ratio is 0.440 to 0.469 with errors of 0.005 to 0.009
  at multipoles 4 to 34, 5 to 8 sd below 0.505 in every bin; the PR3 ratio is 0.485 to 0.507, 2 to
  4 percent low up to multipole 16 and consistent above; and PR4 over PR3 is 0.905 to 0.926 with
  errors of 0.011 to 0.018. In the Galactic plane PR4 is 7 to 10 percent fainter than PR3 at all
  scales, and PR3 is within 4 percent of the K-band extrapolation.
  (B) For abs(b) > 20 both 30 GHz maps are consistent with the K-band extrapolation within 2 sd
  in every bin; over multipoles 2 to 10 together PR4 is 0.481 +- 0.008 and PR3 0.474 +- 0.011,
  about 5 percent and 3 sd low. PR4 over PR3 there is 1.01 to 1.03. I first wrote "4 to 17
  percent" for (B) from the point values alone; with errors it is a 3 sd hint, not a result.
- (A) fits the leakage check above: the plane residual follows the intensity at 0.05 to 0.13
  percent, and the polarized fraction in the inner plane at 30 GHz is of order one percent, so a
  difference of 0.1 percent of intensity between the two releases' bandpass-leakage corrections is
  a 10 percent difference in polarization there. WMAP K and Ka side with PR3.
- (B) is the large-scale discrepancy that Planck 2018 II (App. A), Weiland et al. (2022) and
  Watts et al. (2023) describe; this analysis does not say which map carries it.
- The Ka over K ratio is flat at 0.32 to 0.33 to multipole 34: no scale dependence within WMAP.
  In the last bin (35 to 47) every ratio to K band rises by 6 to 10 percent, Ka included, so K
  band has less power at the pixel scale than the Gaussian 52.8 arcmin window gives it. That is a
  limit of the window matching and it bears on every fit at the 0.01 level of E038 (c).
- **X2 can now be answered in the form it should have had**: within WMAP the ratio does not depend
  on scale, so the scale dependence of the K to 30 GHz index belongs to the 30 GHz maps or to the
  pair, not to K band alone.

## Negative control

The scale cut applied to Ka in the K, Ka, 353 fit: the index must not move by more than 2 sd of
its enlarged error.

## Conclusion

Pending.

## Follow-ups
