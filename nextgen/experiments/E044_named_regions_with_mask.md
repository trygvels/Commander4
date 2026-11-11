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

- **H9, the zero-mean prior at nside 16 (added 2026-11-11 10:13, after R1 and the nside-8 simulation were read,
  before any nside-16 simulation exists).** The nside-8 simulation gives zero-mean minus profile of -1.05, -0.54,
  -0.44, -0.22 sandwich errors; the maps at nside 16 give -3.6, -1.6, -1.5, -0.2. If the difference on the maps
  is the bias of the estimator at this configuration (more and fainter amplitudes at nside 16), two simulated
  realizations at nside 16 give -2.5 to -4.5 at High Latitude, -1.0 to -2.2 in the Spur, -1.0 to -2.0 at
  Intermediate latitude and -0.1 to -0.5 in the Plane, the same in both realizations to 0.1. If they give about
  -1 as at nside 8, then two thirds of the difference on the maps is something the simulated sky (the DR1 model
  smoothed to 5 degrees, a Gaussian CMB, the released noise) does not have, and the paper must say that the
  simulated bias explains only a third of it. Run: `nextgen/scripts/e044_queue_sim16b.sh`, behind the other
  nside-16 simulations.

- **H10, calibration of the sandwich error at High Latitude (added 2026-11-11 10:23, after both null sets were
  read).** Both sets of 60 give a standard deviation of the standardized error of 1.11 and 1.12 at High Latitude
  (pooled 1.104 +- 0.065, scatter of the estimates 0.0321 against a mean sandwich error of 0.0290; coverage 0.62
  +- 0.04). The simulation is drawn from the model of the fit, so the expectation is 1. That is 1.6 standard
  errors in one of four regions, which a fluctuation gives one time in three, but it is the faintest region and
  the one the abstract quotes. A third, independent set of 240 realizations (seed 43) decides: if the sandwich
  error is right its standard deviation of z at High Latitude is 1.00 +- 0.046 and lies below 1.09 with 97.5
  percent probability; if the 10 percent is real it is 1.10 +- 0.05. Rule: above 1.09 in the new set alone, the
  paper says that the sandwich error is about 10 percent low at High Latitude in this configuration and the
  cause is looked for; below, the pooled value of 360 realizations is reported with its error. The other three
  regions are read the same way and none is expected to move.

- **H11, the flat-prior marginal with a free CMB amplitude at High Latitude (added 2026-11-11 10:56, after R5
  was read, before any scan).** In R5 the optimizer ends at the upper bound (-1.6) for High Latitude in both
  conventions, and the curvature of the marginal at the profile optimum is negative there. A scan along the
  High Latitude index (`e044_flat_scan.py`, step 0.1 from -3.6 to -1.6) tests whether that is the shape of the
  function or an artefact of a coupled optimizer started where the curvature is negative. Expected with the
  flat CMB prior: no interior maximum above the profile optimum; ln L_flat rises monotonically from -3.04 to
  -1.6 for 23 GHz and for 30 GHz, by more than 10 (in units of ln L) for 30 GHz, in both slices (others at
  the profile optimum; others at the flat-prior optimum). Expected with the LCDM prior: an interior maximum
  within 0.05 of the fitted flat-prior optima (-3.174 at 23 GHz, -2.696 at 30 GHz) in the slice with the
  others at the profile optimum, shifted by the coupling. The direct evaluation at 30 GHz equals the 23 GHz
  value plus sum_r n_r beta_r ln(30/23) to 1e-6 relative. If the flat-CMB scan has an interior maximum, the
  statement "no interior optimum" is withdrawn and the bound hit is reported as an optimizer result only.
  Amended 10:58, before any scan: the grid is the whole allowed range, -4.6 to -1.6 in steps of 0.15, because a
  marginal that is convex at the profile optimum could also rise toward the lower edge. Expected: the value at
  the upper edge is the largest on the grid for both references.

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

**The zero-mean prior** (`E044/sim_n8_zm.txt`, 24 of 30 realizations read at 10:07; all 30 at 10:23, unchanged:
zero-mean minus truth -0.023, -0.025, -0.008, -0.004 in index, pooled mean z -0.57 +- 0.09): in the same simulated
data the zero-mean optimum is steeper than the profile optimum by 1.05, 0.54, 0.44 and 0.22 sandwich errors (0.030,
0.015, 0.009, 0.004 in index), with a scatter between realizations of 0.02 or less. The profile optimum is unbiased in
these simulations; so the zero-mean prior with a free scale per region is biased at this configuration, by an
amount that does not depend on the noise draw. On the maps at nside 16 the difference is 0.104, 0.043, 0.030,
0.004: the same sign and order, about three times larger in index. E041 had found up to 0.6 sd on its ladder
when the noise pattern differs between amplitudes; here, with three components, dense noise and a faint
High Latitude region, it is larger. This changes one sentence of the paper's abstract and one recommendation:
the zero-mean prior is not a substitute for the profile likelihood on these maps.

### R3, the white CMB prior judged under the LCDM covariance (read 2026-11-11 10:21)

`E044/white.txt` (white prior of 0.2 uK per element; profile and flat estimators; `lowres.mismatch` with the LCDM
covariance as the truth; canary, no OVERFLOW line; every fit converged). chi2_min 16932.8 / 17420 = 0.972.

| | High Latitude | Spur | Intermediate | Plane |
|---|---|---|---|---|
| profile index, white prior | -3.1728 | -3.1365 | -3.1035 | -3.0396 |
| error the white-prior fit quotes | 0.0264 | 0.0248 | 0.0187 | 0.0143 |
| the same over the LCDM fit's error | 0.906 | 0.902 | 0.966 | 0.956 |
| first-order true error if the CMB is LCDM | 0.0321 | 0.0301 | 0.0204 | 0.0157 |
| true over quoted | 1.218 | 1.214 | 1.086 | 1.098 |
| true over the LCDM fit's error | 1.103 | 1.096 | 1.050 | 1.049 |
| first-order bias | -0.0023 | -0.0004 | -0.0009 | -0.0003 |
| bias over quoted error | -0.09 | -0.02 | -0.05 | -0.02 |
| white minus LCDM optimum | -0.0258 | +0.0028 | -0.0124 | -0.0005 |
| expected sd of that difference, sqrt(true^2 - LCDM^2) | 0.0136 | 0.0123 | 0.0062 | 0.0047 |
| difference minus bias, in that sd | -1.74 | +0.26 | -1.87 | -0.05 |
| flat 30 GHz minus flat 23 GHz, white prior | +0.4156 | +0.0939 | +0.1769 | +0.0675 |
| the same over first order (marginal curvature; J_c) | 0.79; 0.74 | 0.89; 0.84 | 0.89; 0.85 | 0.93; 0.91 |
| chi-square per element and band | 0.573 | 0.588 | 0.570 | 0.616 |

- **H5 holds for the white prior.** The quoted error is too small by 22, 21, 9 and 10 percent (written: 5 to 40,
  most at High Latitude; the Spur equals it). The bias is 0.09 quoted sd or less (written: below 0.3). The two
  optima differ by 1.74, 0.26, 1.87 and 0.05 of the expected sd of their difference (written: below 2). Two of
  the four are near the limit and have the same sign; the two indices are correlated at 0.2, so this is weaker
  than two independent 1.8 sd events, and I read it as consistent, not as confirmed to better than that.
- The white-prior fit quotes an error that is 3 to 10 percent smaller than the LCDM fit's, and if the sky is
  LCDM its true error is 5 to 10 percent larger than the LCDM fit's. The white prior looks better and is worse.
- Size against the earlier statements: E040 stage 4c (nside 8, other regions, simulation) had "a quarter to a
  third"; relative to the matched prior on the same skies 14 to 27 percent. Here, exact to first order at the
  configuration of the map fit: 9 to 22 percent.
- Flat prior with the white CMB prior: the pivot stays at 23.2 to 23.3 GHz and the move between the two
  reference frequencies is 0.42, 0.09, 0.18, 0.07 (15.7, 3.8, 9.4, 4.7 of its own quoted errors), smaller than
  with the LCDM prior in the Spur (0.154), where the white prior holds the small-scale CMB tighter.
- Limits: the calculation is first order in the mismatch and takes the LCDM covariance of the paper (multipoles
  2 to 47, floor 0.01 uK) as the truth. Its unit test (1500 fits) is at a toy configuration.

### R4, the chi-square of the central fit against its null (read 2026-11-11 10:39)

A check I should have made when R1 was read. chi2_min of the central fit, with the penalty of the CMB prior
inside it (55), is 16953.8 for 17420 degrees of freedom: ratio 0.9732, **2.50 standard deviations low** (sd of the
ratio sqrt(2 / dof) = 0.0107). Per element and band the expectation is 0.598, not 0.6, and with sampling errors
of 0.009, 0.030, 0.012, 0.017 the four regions are at -2.3, -0.2, -2.2 and +0.9 standard deviations. The deficit
is in the two faint regions. Per band (residual chi-square per element): K 0.109, Ka 0.878, Q 0.954, V 0.969.

Null (`nextgen/scripts/e044_chi2_null.py`, `E044/sim_n8_chi2null.txt`, 40 realizations of the first nside-8 set,
regions, mask and prior of the fit): chi2_min / dof = 1.0019 +- 0.0027, scatter 0.0173 against sqrt(2 / dof) =
0.0215; per element and band 0.597, 0.601, 0.602, 0.592. The expectation is right.

- R1 said "H3 holds, better than written". H3 holds as written (every region within 5 percent of 0.6). But
  "the model fits" needs the other half: nothing is misfitted high, and the total is low at 2.5 sd, which a
  fluctuation gives one time in 80.
- Readings: a fluctuation; or released WMAP covariances that overstate the noise at high latitude (1.3
  percent in rms if spread evenly over the bands; Q and V carry most of the residual degrees of freedom). In
  the second case every error of this experiment is too large by about that factor, which is small against
  everything else here. Not rescaled. The white-prior fit has the same ratio (0.972), so the CMB prior is not
  the cause, and the old regions had the same low values where they were faint (0.584, 0.592).
- Not checked: whether the deficit is in particular modes of the covariances (the transmission-imbalance
  modes, the large-scale modes with 65 to 164 times the median variance found in E038). That is a question for
  Paper B.

### R5, a flat prior on the CMB (read 2026-11-11 10:53)

`E044/flat.txt` (one free CMB amplitude per element; profile and flat estimators; canary, no OVERFLOW line;
converged). chi2_min 11101.4 for 11612 degrees of freedom: ratio 0.956, **3.35 standard deviations low** (here
chi2_min has no prior penalty and exactly two degrees of freedom per element). Per element and band 0.385, 0.408,
0.372, 0.385 against 0.4: -2.0, +0.3, -2.9 and -1.1 sampling errors.

| | High Latitude | Spur | Intermediate | Plane |
|---|---|---|---|---|
| profile index, flat CMB prior | -3.0406 | -3.1219 | -3.0167 | -3.0384 |
| sandwich error | 0.0618 | 0.0599 | 0.0444 | 0.0353 |
| the same over the LCDM fit's error | 2.12 | 2.18 | 2.29 | 2.36 |
| flat minus LCDM optimum | +0.1063 | +0.0174 | +0.0744 | +0.0006 |
| in units of sqrt(sd_flat^2 - sd_lcdm^2) | +1.95 | +0.33 | +1.86 | +0.02 |
| T / J_c | 0.083 | 0.008 | 0.027 | 0.009 |
| n sd^2 (index variance carried by one amplitude) | 11.4 | 0.96 | 3.35 | 1.08 |
| flat amplitude priors at 23 GHz | -1.600 (bound) | -2.778 | -2.414 | -2.809 |
| flat amplitude priors at 30 GHz | -1.600 (bound) | -2.577 | -1.737 | -2.519 |
| curvature of the marginal at the profile optimum (diagonal) | -366 | +340 | +646 | +920 |
| J (diagonal) | 316 | 288 | 549 | 832 |

- **H5 holds for the flat prior by its letter**: errors 2.1 to 2.4 times larger (written: 2 to 2.5); the optimum
  within two nested sd in every region (written: within two). But two regions are at +1.95 and +1.86, the same
  two that were at -1.74 and -1.87 for the white prior, with opposite sign. Read together: the more tightly the
  CMB amplitude is held, the steeper the index at High and Intermediate latitude (-3.173, -3.147, -3.041 at
  High Latitude for white, LCDM, flat; -3.104, -3.091, -3.017 at Intermediate), by an amount at the edge of
  what the calculation allows in each step. The three fits share the data, so these are not three tests. I
  do not have a cause and the paper says what was measured.
- **Not expected in this form: with a free CMB amplitude and flat priors on all amplitudes the High Latitude
  optimum is at the upper bound of the allowed range (-1.6) in both conventions.** I had written "the move at
  high latitude" as a number in the abstract template; there is no interior optimum to move. The curvature of
  the flat-prior marginal at the profile optimum is negative for High Latitude (the volume term's curvature is
  negative and larger than J: -366 against J = 316; smallest eigenvalue of the 4 by 4 matrix -376). Sect. 4 of
  the paper says this can happen with further free amplitudes; E040 stage 3 saw it on the old regions for a
  band set without K band. Here it is with K band, on |b| > 30 degrees, masked. The index variance carried by
  one amplitude is 11 at High Latitude, so any offset of the reference from the pivot by a tenth in ln nu is a
  shift of order one. The other three optima are interior but conditional on a neighbour at the bound, and are
  0.2 to 1.3 flatter than the profile optimum. A scan along the High Latitude index is pre-registered as H11.
- The chi-square is low again, more so than with the LCDM prior (R4): 3.35 sd with two residual degrees of
  freedom per element. A fluctuation of that size is one in 1200. This strengthens the second reading of R4:
  the released covariances overstate the noise of the residual of this fit by about 2 percent in rms.
- **The deficit is not uniform over the combinations of bands (noted 11:40).** The LCDM-prior fit and the
  flat-prior fit are nested, and under the model the difference of their chi2_min is chi-square with one degree
  of freedom per element (the fitted CMB amplitude against its noise plus prior variance). It is 16953.8 -
  11101.4 = 5852.4 for 5808: ratio 1.008 +- 0.019. The other two degrees of freedom per element have 0.956 +-
  0.013. The two differ by 0.052 +- 0.023 (2.3 sd). A uniformly overstated noise would put both at 0.956. So
  either the deficit is confined to the directions that the three amplitudes do not span, or it is everywhere
  and the direction of the CMB spectrum carries an excess of about 5 percent of its noise variance beyond the
  LCDM prior. The second would also give an index that depends on how tightly the CMB amplitude is held, which
  is the pattern of R3 and R5. Not pursued here (Paper B; literature in
  `research/2026-11-11_wmap_lowres_covariance_normalization.md`, which finds no published statement that the
  matrices overstate the noise, two uncommented V-band values of 0.96 to 0.98 with the same matrices, and
  model-free tests that need the single-year files).

### R6, the arms (all read by 2026-11-11 12:08; written 12:08)

Every fit converged; no OVERFLOW line in any error file. Baseline: R1 (LCDM prior, all scales, four WMAP bands and
353 GHz). For the arms that use less of the same data under the same model, the change is also given in units of the
standard deviation that the model predicts for it.

| Arm | High Latitude | Spur | Intermediate | Plane | chi2_min / dof |
|---|---|---|---|---|---|
| baseline | -3.1469 | -3.1393 | -3.0911 | -3.0390 | 0.9732 |
| sandwich error | 0.0291 | 0.0275 | 0.0194 | 0.0150 | |
| index from l <= 30: change | -0.0009 | +0.0087 | -0.0204 | +0.0156 |  |
| index from l <= 30: error | 0.0307 | 0.0309 | 0.0248 | 0.0219 | |
| index from l <= 30: change over sqrt(sd_arm^2 - sd_base^2) | -0.10 | +0.62 | -1.31 | +0.97 | |
| index from l <= 20: change | +0.0036 | -0.0283 | -0.0468 | +0.0106 |  |
| index from l <= 20: error | 0.0343 | 0.0402 | 0.0395 | 0.0349 | |
| index from l <= 20: change over sqrt(sd_arm^2 - sd_base^2) | +0.20 | -0.96 | -1.36 | +0.33 | |
| index from l <= 12: change | +0.0263 | -0.1928 | -0.1927 | +0.1412 |  |
| index from l <= 12: error | 0.0624 | 0.2216 | 0.1873 | 0.1086 | |
| index from l <= 12: change over sqrt(sd_arm^2 - sd_base^2) | +0.48 | -0.88 | -1.03 | +1.31 | |
| K and Ka only: change | +0.0124 | -0.0047 | +0.0366 | +0.0248 | 0.9619 |
| K and Ka only: error | 0.0364 | 0.0346 | 0.0247 | 0.0194 | |
| K and Ka only: change over sqrt(sd_arm^2 - sd_base^2) | +0.57 | -0.22 | +2.40 | +2.00 | |
| dust index 1.48: change | -0.0039 | -0.0054 | -0.0060 | -0.0208 | 0.9746 |
| dust index 1.58: change | +0.0036 | +0.0049 | +0.0055 | +0.0188 | 0.9730 |
| no drift correction: change | -0.0011 | -0.0011 | -0.0011 | -0.0014 | 0.9733 |
| 353 GHz, all bolometers: change | +0.0005 | +0.0007 | +0.0008 | +0.0027 | 0.9732 |
| delta bandpasses: change | -0.0500 | -0.0505 | -0.0483 | -0.0464 | 0.9730 |
| LCDM prior, floor 0.11 uK: change | -0.0020 | -0.0005 | -0.0008 | -0.0004 | 0.9712 |
| K and Ka only minus four bands, both from l <= 20 | +0.0177 | -0.0066 | +0.0564 | +0.0153 | |
| the same over its expected sd | +0.72 | -0.22 | +1.91 | +0.52 | |

For comparison, the low-pass arm on the old regions without a mask (E042 `lowpass20_lcdm` against E040
`s4_bpxv_n16_wm_lcdm`): changes of -0.024, -0.031, -0.055, -0.084, which are 2.8, 3.3, 8.1 and 10.4 of the expected
standard deviation.

- **H7 holds in its reading and fails in its number.** I wrote that the Plane would still steepen by 0.01 to 0.04 at
  l <= 20; it changes by +0.011 (+0.33 expected sd), the other sign. High Latitude and Intermediate are within two
  (+0.20, -1.36), as written. With the mask no region depends on the angular scale beyond the noise at l <= 30, 20
  or 12 (largest 1.4). On the old regions it was 3 to 10. The scale dependence of F49 went with the masked emission.
- **H8, dust index: holds in size, fails in pattern.** Changes of 0.004 to 0.021 for 0.05 in the dust index (written:
  0.03 or less); the largest is in the Plane, not at High Latitude as I wrote (the dust is bright in the Plane).
  By-product, not planned: chi2_min is 16977.1, 16953.8 and 16949.4 at 1.48, 1.53 and 1.58, and a parabola through
  the three has its minimum at 1.567 +- 0.016. Planck 2018 XI has 1.53 +- 0.02 for polarization from the HFI bands.
  Consistent at 1.4 sd; at an uncertainty of 0.02 the dust index is worth 0.002 to 0.008 in the synchrotron index.
- **H8, bandpasses: holds.** Delta bands steepen every region by 0.046 to 0.051 (written: 0.05 to 0.08; the Plane
  is just under). No drift correction: 0.0011 to 0.0014 (written: 0.002 or less). All bolometers for the 353 GHz
  bandpass: 0.0005 to 0.0026.
- **H8, K and Ka alone: holds, near its limit.** +0.57, -0.22, +2.39 and +2.01 expected sd (written: within 2.5).
  Read from l <= 20 on both sides: +0.72, -0.21, +1.91, +0.52. The index from K and Ka alone is flatter at
  Intermediate latitude and in the Plane by 0.037 and 0.025; that is the sign of a spectrum that steepens with
  frequency, or of something at Q and V that the model lacks. Two regions at 2 sd in one arm, not a detection.
  Its chi2_min is 0.962 of its expectation (2.1 sd low; one residual degree of freedom per element).
- Floor of the LCDM prior at 0.11 uK: 0.0004 to 0.0020. As on the old regions (0.0006), no effect.
- Over all arms that test the model (flat CMB prior, three low-pass cuts, K and Ka, K and Ka at l <= 20: 24
  numbers that share their data) the largest is 2.39. No arm rejects the model on these regions with this mask.

## Negative control

The flat-prior marginal at two reference frequencies; and the old regions with only Tau A masked (H6).

## Conclusion

(to be filled)
