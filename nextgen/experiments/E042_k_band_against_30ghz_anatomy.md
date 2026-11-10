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

Pending.

## Negative control

The scale cut applied to Ka in the K, Ka, 353 fit: the index must not move by more than 2 sd of
its enlarged error.

## Conclusion

Pending.

## Follow-ups
