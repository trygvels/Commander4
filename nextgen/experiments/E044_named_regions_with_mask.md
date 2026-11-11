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

(to be filled)

## Negative control

The flat-prior marginal at two reference frequencies; and the old regions with only Tau A masked (H6).

## Conclusion

(to be filled)
