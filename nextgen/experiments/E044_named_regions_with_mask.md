# E044: the estimator comparison on regions that are what they are named, with a processing mask

- **Date started / finished:** 2026-11-11 08:45 / open
- **Branch and commit:** `nextgen`; code in `nextgen/lowres/regions.py`, `mask.py`, `mismatch.py`, `block.py`; driver `nextgen/scripts/e044_fit.py`
- **Open question addressed:** Q21, Q24; the real-data section of Paper A1 after the second referee read
- **Validation level:** L4 (real maps), with L1 checks of the new pieces (unit tests) and L2 nulls at the configuration of the fit

## Why this experiment exists

The second independent referee read of Paper A1 (`reviews/2026-11-11_referee_report_paperA1_draft3.md`) found that
the four regions used in E022 and E040 are not what they are named, and I confirmed it (2026-11-11 08:44):

- They are a reconstruction from the released BeyondPlanck index map (nearest plateau; "Galactic Center" was my
  cut at 25 degrees among the prior-only pixels). The region called High Latitude contains the inner Galactic
  plane on both sides of that cut: 140 of its 2194 pixels lie within 10 degrees of the plane and hold 55 percent
  of its squared fitted amplitude; pixels above 30 degrees latitude hold 27 percent.
- No mask was applied. Tau A sits in the Plane region with a fitted polarized amplitude of 175 uK against a
  region median of 38, and by the referee's first-order calculation its patch carries 58 percent of the Plane's
  block variance and moves the Plane index by 0.07.
- The block error of "High Latitude" (0.044) is 65 percent the Galactic Center's scores, through the coupling of
  the two indices; "the High Latitude region is not one index" was not supported.

Every statement of E040, F48 to F54 and Paper A1 draft 3 that names a region is a statement about those
surrogate regions without a mask. The differential comparison of estimators is not affected (every row used the
same regions), but its labels are, and so is everything that interprets misfit.

## Question

On regions defined by geometry, with the brightest emission and Tau A masked: what do the four estimators give on
the same maps, does the reference-frequency identity hold at the predicted size, and how do the profile index,
its error and the misfit diagnostics change?

## Target and priors

As E040: profile likelihood per region index with flat priors on the synchrotron and dust amplitudes per element
and a Gaussian prior on the CMB amplitudes (LCDM covariance on the kept elements, the kept block of
`lcdm_cmb_precision`; white 0.2 uK and flat as arms); flat-prior marginal with the synchrotron amplitude defined
at 23 and at 30 GHz; zero-mean prior with a free scale per region.

## Setup

- Data: WMAP nine-year K, Ka, Q, V (res 4, released dense inverse covariances), Planck PR4 353 GHz; measured
  bandpasses with the drift correction; all bands brought to a 60 arcmin window at nside 16.
- **Regions** (`lowres.regions.named_regions`): High Latitude, |b| > 30 degrees outside the Spur; Spur, the
  northern piece (b > 20 degrees, -45 < l < 65 degrees) of the pixels that the released BeyondPlanck index map
  assigns to its Spur value; Intermediate, 10 to 30 degrees of latitude outside the Spur; Plane, |b| < 10
  degrees.
- **Mask:** the 5 percent of pixels with the largest polarized amplitude in the Cosmoglobe DR1 synchrotron model
  at 30 GHz smoothed to 5 degrees (154 pixels, threshold 28.6 uK_RJ; 141 of them within 10 degrees of the plane),
  and a disc of 7.5 degrees around Tau A (14 pixels). 168 of 3072 pixels. This follows the description of the
  processing mask of Svalheim et al. (2023) as the referee quoted it, not their file, which we do not have.
- Masked elements are removed exactly: the inverse covariance of the kept elements is the Schur complement
  (`lowres.mask.marginalize`), after the window matching. Kept: 5808 of 6144 elements; 2974, 266, 1698 and 870 in
  the four regions; nside-2 patches with kept elements: 38, 6, 39, 24.
- New pieces, each with a unit test (`lowres/tests/test_block_prior_lowpass.py`, 11 tests): the Schur
  complement against the inverse of the kept block; block errors from a region's own scores (`CR3_own`), which do
  not inherit a neighbour's misfit; the first-order bias and covariance of a fit whose CMB prior is not the truth
  (`lowres.mismatch`, against 500 and 1500 simulated fits: scatter 0.113, 0.089 predicted 0.113, 0.091 where the
  matched formula gives 0.090, 0.083; bias +0.016 +- 0.003, +0.021 +- 0.002 predicted +0.021, +0.020).
- Every run on the sound Cholesky path under the malloc canary (E043). One heavy job at a time.

## Expectations, written 2026-11-11 09:00 before any fit of this experiment at nside 16

End-to-end tests of the driver at nside 8 (which has the degrade bias of E040 R8 and is not a result) were run
and seen; nothing at nside 16.

- **H1, the reference frequency.** With the LCDM prior the flat-prior optimum stated at 23 GHz is within one
  sandwich sd of the profile optimum in every region (the pivot is set by the bands, 23.3 GHz). Stated at 30 GHz
  it is flatter in every region. Sizes from n sigma^2 ln(30/23): Plane 0.03 to 0.08, Spur 0.02 to 0.06,
  Intermediate 0.1 to 0.3, High Latitude 0.3 or more (low signal to noise, many amplitudes; beyond first order).
  The first-order prediction with the curvature of the fitted marginal is within 15 percent of the measured move
  where the move is below 0.3 and high above that. The prediction that needs no flat-prior fit (J_c) is worse,
  by up to 25 percent.
- **H2, the profile index.** I have no sharp expectation and will not quote one as a measurement. Ranges: Plane
  -3.15 to -3.0 (the referee's first-order value without Tau A on the old region was -3.04); Spur -3.2 to -3.05;
  Intermediate -3.2 to -3.0; High Latitude -3.3 to -3.0 with an error of 0.03 to 0.06.
- **H3, goodness of fit.** Chi-square per element and band within 10 percent of 0.6 at High Latitude, in the Spur
  and at Intermediate latitude, within 25 percent in the Plane.
- **H4, block errors.** From each region's own scores, ratio to the sandwich error: High Latitude below 2
  (inside the null); Intermediate 1.5 to 3; Plane 2 to 4. The Spur has 6 patches and gets no block error. A
  region is called "more than one index" only if its own-score ratio exceeds every value of a null simulated
  with these regions, this mask, these covariances and this prior, and even then the text names the other
  readings (noise model, dust, beams, leakage).
- **H5, the CMB prior.** `lowres.mismatch` on the white-prior fit under the LCDM covariance: the quoted error is
  too small by 5 to 40 percent, most at High Latitude; the bias is below 0.3 sd. The white-prior optimum differs
  from the LCDM one by less than twice sqrt(var_white,true - var_lcdm) in every region. With a flat CMB prior the
  errors are 2 to 2.5 times larger and the optimum is within two of sqrt(sd_flat^2 - sd_lcdm^2) of the LCDM one.
- **H6, Tau A.** On the old regions with only the Tau A disc masked: the Plane index moves from -2.964 by -0.05
  to -0.10 and its chi-square per element and band falls from 0.770 to below 0.68; the other three regions move
  by less than 0.01.
- **H7, angular scale.** If compact structure made the scale dependence of E042, masking shrinks it. With the
  index from multipoles up to 20: the Plane steepens by 0.01 to 0.04 (it was 0.084), and High Latitude and
  Intermediate are within two nested sd of the all-scales value. If the Plane still steepens by more than 0.06,
  compact structure in the masked pixels was not the cause.
- **H8, arms.** Dust index 1.48 and 1.58: the synchrotron index moves by 0.03 or less in every region, more at
  High Latitude than in the Plane. Delta bandpasses: steeper by 0.05 to 0.08. No drift correction: 0.002 or less.
  K and Ka alone: within 2.5 nested sd of the four-band index in every region.

## Command

```bash
C=~/data/c4_tools/libcanary.dylib; P=".venv/bin/python -I nextgen/scripts/e044_fit.py"; E=nextgen/experiments/E044
DYLD_INSERT_LIBRARIES=$C $P --prior lcdm --estimators profile flat zm --block --out $E/lcdm_all.npz > $E/lcdm_all.txt
DYLD_INSERT_LIBRARIES=$C $P --geometry bpxv --tau-a-only --prior lcdm --block --out $E/bpxv_tauA.npz > $E/bpxv_tauA.txt
DYLD_INSERT_LIBRARIES=$C $P --prior white --mismatch --estimators profile flat --out $E/white.npz > $E/white.txt
DYLD_INSERT_LIBRARIES=$C $P --prior flat --estimators profile flat --out $E/flat.npz > $E/flat.txt
# arms and low-pass: nextgen/scripts/e044_queue.sh
```

## Result

### R1, the central fit (read 2026-11-11 09:43)

`E044/lcdm_all.txt`: LCDM prior, all estimators, block errors; sound path, canary, no OVERFLOW line. chi2_min
16953.8 for 17420 nominal degrees of freedom (0.973). Every fit converged.

| | High Latitude | Spur | Intermediate | Plane |
|---|---|---|---|---|
| elements kept | 2974 | 266 | 1698 | 870 |
| rms fitted synchrotron amplitude at 23 GHz (uK) | 11.3 | 33.4 | 17.5 | 26.3 |
| profile index | -3.1469 | -3.1393 | -3.0911 | -3.0390 |
| sandwich error | 0.0291 | 0.0275 | 0.0194 | 0.0150 |
| block error, joint (CR3) | 0.0480 | 0.0288 | 0.0292 | 0.0193 |
| block error, own scores (CR3) | 0.0475 | 0.0274 | 0.0280 | 0.0192 |
| own-score block error over sandwich | 1.63 | 1.00 | 1.44 | 1.28 |
| effective number of patches; largest leverage | 12.9; 0.22 | 3.5; 0.46 | 13.2; 0.16 | 10.3; 0.18 |
| flat prior at 23 GHz minus profile | -0.0275 | -0.0079 | -0.0112 | -0.0042 |
| flat prior at 30 GHz minus profile | +0.4509 | +0.1459 | +0.1930 | +0.0683 |
| flat 30 GHz minus flat 23 GHz | +0.4783 | +0.1538 | +0.2042 | +0.0725 |
| the same in sandwich errors | 16.4 | 5.6 | 10.5 | 4.8 |
| first order, curvature of the fitted marginal | +0.6382 | +0.1834 | +0.2377 | +0.0780 |
| first order, J_c | +0.6921 | +0.1959 | +0.2505 | +0.0802 |
| zero-mean prior minus profile | -0.1038 | -0.0429 | -0.0299 | -0.0035 |
| chi-square per element and band | 0.577 | 0.591 | 0.572 | 0.614 |
| T / J_c | 0.047 | 0.004 | 0.013 | 0.004 |

Null for the block errors at nside 8 with these regions, this mask and the LCDM prior (`E044/sim_n8_blocknull.txt`, 60
realizations, read 09:34): pooled mean z -0.08 +- 0.065, sd of z 1.05, coverage 0.65. Own-score CR3 over sandwich,
median / 84th / 97.5th percentile / maximum: High Latitude 1.34 / 1.66 / 2.09 / 2.17; Spur 1.32 / 2.09 / 3.29 / 3.58;
Intermediate 1.12 / 1.35 / 1.67 / 1.77; Plane 1.13 / 1.50 / 1.73 / 1.84. Effective number of patches 12.1, 3.1,
11.7, 9.5. One thing I did not expect: the mean z of the Spur is -0.39 +- 0.13 (mean estimate minus truth -0.011);
the other three are within 0.1. A second, independent set of 60 realizations is running to tell a fluctuation from
a bias of this small region.

Against the expectations:

- **H3 holds, better than written.** Chi-square per element and band is 0.577, 0.591, 0.572, 0.614 against 0.6: every
  region within 5 percent. On the old regions without a mask it was 0.58, 0.59, 1.78, 0.77. The misfit of the earlier
  fits was the brightest emission and Tau A.
- **H2 holds.** -3.147 +- 0.029, -3.139 +- 0.028, -3.091 +- 0.019, -3.039 +- 0.015, all inside the ranges written. The
  Plane is at the value the referee's first-order calculation gave for the old Plane without the Tau A patch
  (-3.04). The four indices now steepen with latitude; with these errors High Latitude and Plane differ by
  0.108 +- 0.033.
- **H4 fails on the low side for two regions.** Own-score block error over sandwich: 1.63, 1.00, 1.44, 1.28. I
  wrote 1.5 to 3 for Intermediate and 2 to 4 for the Plane. Against the null, all four are inside: High Latitude
  1.63 against an 84th percentile of 1.66; Intermediate 1.44 against a 97.5th of 1.67; Plane 1.28 against an
  84th of 1.50. **With the mask there is no evidence in the block scores that any region is more than one index.**
  The large block errors of E040 were the masked emission, as the second referee argued.
- **H1.** The flat prior at 23 GHz is within one sandwich error of the profile optimum in every region (0.95, 0.29,
  0.58, 0.28 sd), as written. At 30 GHz it is flatter in every region, by 0.478, 0.154, 0.204 and 0.073 against
  the 23 GHz statement: 16.4, 5.6, 10.5 and 4.8 sandwich errors. Sizes: Plane and Intermediate inside the ranges
  written, High Latitude "0.3 or more" (0.48); the Spur is 0.154 where I wrote 0.02 to 0.06 (I used an error of
  0.02 for a region whose error is 0.028 and forgot the coupling to High Latitude). First order with the
  curvature of the fitted marginal: measured over predicted 0.75, 0.84, 0.86, 0.93; within 15 percent only for
  Intermediate and Plane, 16 percent for the Spur, and 25 percent high at High Latitude where the move is 0.48.
  With J_c: 0.69, 0.79, 0.82, 0.90, worse by 4 to 8 points as expected, and beyond the 25 percent I wrote at
  High Latitude. The first-order prediction is always high.
- **Not expected: the zero-mean prior.** It is steeper than the profile optimum by 0.104, 0.043, 0.030 and 0.004:
  3.6, 1.6, 1.5 and 0.2 sandwich errors. On the old regions it was within 0.8. At the low signal to noise of the
  new High Latitude region (T / J_c 0.047) with the scan-modulated WMAP noise the two estimators part. A
  simulation of the zero-mean prior at this configuration is running (`E044/sim_n8_zm.txt`); its first realization
  has the same sign in all four regions.

### R2, the old regions with only the disc around Tau A masked (read 2026-11-11 09:54)

`E044/bpxv_tauA.txt` (LCDM prior; 28 of 6144 elements removed). Against the unmasked fit of E040 (`s4_bpxv_n16_wm_lcdm`):

| Old region | Index, unmasked | Index, Tau A masked | Move | chi2 per element and band, before and after | Own-score block error over sandwich |
|---|---|---|---|---|---|
| High Latitude | -3.0825 | -3.0884 +- 0.0125 | -0.006 | 0.584, 0.585 | 2.25 |
| Spur | -3.0799 | -3.0986 +- 0.0167 | -0.019 | 0.592, 0.591 | 1.57 |
| Galactic Center | -3.0807 | -3.0823 +- 0.0089 | -0.002 | 1.782, 1.786 | 8.6 |
| Plane | -2.9641 | -3.0483 +- 0.0144 | -0.084 | 0.770, 0.596 | 2.28 |

- **H6 holds for the Plane**: the index moves by -0.084 (written: -0.05 to -0.10; the referee's first-order value was
  -0.074), six of its statistical errors, and the chi-square falls from 0.770 to 0.596 (written: below 0.68), that
  is to the expectation. One source did both. The joint block error of that region falls from 0.100 to 0.032.
- H6 fails in one detail: the old Spur moves by -0.019 where I wrote "less than 0.01". Its strips along the edge
  of the plane reach the masked disc, and it is coupled to the Plane through the noise.
- The old Galactic Center still does not fit (1.79); that is the brightest emission, which the full mask removes.

### Simulations at this configuration, nside 8 (read 2026-11-11 10:07)

**Second null set** (`E044/sim_n8_blocknull_b.txt`, 60 realizations with other seeds): mean z -0.08, -0.02, +0.02, -0.02;
sd of z 1.12, 1.05, 1.02, 0.87. The Spur's -0.39 +- 0.13 of the first set was a fluctuation: pooled over 120
realizations the four mean z are -0.05, -0.20, +0.02, +0.03 (se 0.09) and the sd of z 1.11, 1.02, 1.01, 0.98.
Own-score block error over sandwich, pooled: median 1.33, 1.30, 1.11, 1.12; 97.5th percentile 2.02, 3.11, 1.69,
1.72. Of the 120 simulations 15, 85, 12 and 31 lie above the map values (1.63, 1.00, 1.44, 1.28). All four
regions are inside the null.

**The zero-mean prior** (`E044/sim_n8_zm.txt`, 24 of 30 realizations read): in the same simulated data the
zero-mean optimum is steeper than the profile optimum by 1.05, 0.54, 0.44 and 0.22 sandwich errors (0.030, 0.015,
0.009, 0.004 in index), with a scatter between realizations of 0.02 or less. The profile optimum is unbiased in
these simulations; so the zero-mean prior with a free scale per region is biased at this configuration, by an
amount that does not depend on the noise draw. On the maps at nside 16 the difference is 0.104, 0.043, 0.030,
0.004: the same sign and order, about three times larger in index. E041 had found up to 0.6 sd on its ladder
when the noise pattern differs between amplitudes; here, with three components, dense noise and a faint
High Latitude region, it is larger. This changes one sentence of the paper's abstract and one recommendation:
the zero-mean prior is not a substitute for the profile likelihood on these maps.

## Negative control

The flat-prior marginal at two reference frequencies; and the old regions with only Tau A masked (H6).

## Conclusion

(to be filled)
