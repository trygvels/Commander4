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

## Negative control

The scale cut applied to Ka in the K, Ka, 353 fit: the index must not move by more than 2 sd of
its enlarged error.

## Conclusion

Pending.

## Follow-ups
