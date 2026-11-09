# E034: map ingest checks from the pipeline audit (model-error floor, degrade before smooth)

- **Date started / finished:** 2026-11-09 22:44 / 2026-11-09 22:47
- **Branch and commit:** `nextgen` at 52061ad (working tree with E032 and E033 records uncommitted)
- **Open question addressed:** Q1 (real maps), Q7, Q15; audit items B14/D9 and B1/D5 of
  `../reviews/2026-11-09_pipeline_assumptions_audit.md`
- **Validation level:** L0 (deterministic checks of the data path against independent references)

## Question

The audit read the map ingest code and raised two hypotheses that would change how the real-map
results of E014 (F3, F14, F17) are read. Neither was run by the audit.

1. *Floor (B14).* `add_signal_fraction_to_rms` is added per native pixel before the degradation
   and the smoothing, so a model-error floor meant as 5 percent of the signal averages down like
   white noise, about 26 times at E014 settings. If so, the "rms (with floor)" column of the E014
   high-latitude table (3.3, 3.1, 5.1, 13.4, 23 uK_RJ at 30, 70, 100, 353, 857 GHz) is not what
   the chain used, the chain's HFI bands were the tightest bands and not slack ones, and the
   high-latitude reduced chi-squares of that table (0.0 to 0.1 for HFI) are wrong by one to two
   orders of magnitude.
2. *Order (B1).* Maps are degraded with `hp.ud_grade` before they are smoothed to the common
   beam. Power below the coarse pixel scale is aliased by an amount that depends on the band's
   native beam, so bands stop sharing a beam and a sky, at a level the noise model does not know.

What changes with the answer: if either holds at a level comparable to the per-pixel noise that
reaches the likelihood, every real-map number so far carries an ingest error on top of the model
errors already recorded, the fix (smooth at native resolution, then degrade; add the floor after
both) comes before any further real-map run, and F3, F14 and F17 get a caveat on their units of
"posterior sd".

## Target and priors

No sampler. What "correct" means for each check:

- Floor: the variance that reaches the likelihood should be the propagated white-noise variance
  plus the square of the stated model-error fraction times a signal scale. A model error is not
  white, so it must not average down under degradation and smoothing.
- Order: with all bands brought to one beam, the difference between two bands observing the same
  CMB sky must vanish up to noise. Any rms of the difference is unmodelled.

## Setup

- Floor: the E014 tuned parameter file `nextgen/params/q1/param_q1_n64_ridge_npipetuned.yml` (17
  bands, nside 512 inputs, eval_nside 64, common beam 120 arcmin, floor 0.05 on the three HFI
  bands), read through Commander4's own `read_data_map_from_file`, as configured and with the
  floor key removed; plus the degraded rms straight from the files. Median over abs(b) > 30 deg.
  Real maps under `~/data/q1_inputs` and `~/data/beyondplanck` (downloaded data: Python run with -I).
- Order: one Gaussian CMB realization at nside 512, lmax 1535, from a smooth stand-in TT spectrum
  (first peak 5500 uK^2 at l = 220, damping tail), observed with Gaussian beams of 5, 13.2, 32.3,
  52.8 and 60 arcmin and the nside-512 pixel window. Code path: `ud_grade` to nside 64, then
  `hp.smoothing` by sqrt(120^2 - fwhm^2) (what the noise-weighted smoothing reduces to for uniform
  weights). No noise, no foregrounds, no point sources, so the result is a lower bound.
- Single thread, beside the running E032 chain.

## Reference

- Floor: columns computed three ways from the same files (reader as configured, reader without
  the floor, degraded file rms with the nominal floor added afterwards).
- Order: the same alms synthesized directly at the nside-64 pixel centres with the 120 arcmin
  beam, with and without the nside-64 pixel window; and the other order (smooth at nside 512, then
  degrade).

**Expectations, written before the first run:**

- Floor. (a) The column "degraded rms with the nominal floor" reproduces the E014 table column to
  20 percent (3.3, 3.1, 5.1, 13.4, 23 at 30, 70, 100, 353, 857 GHz). (b) The rms the likelihood
  gets is 2.9 to 3.3 times below the degraded rms for bands without a floor, and 13 to 26 times
  below the E014 column for the three floored bands; the audit's point predictions are 1.0, 0.9,
  0.4, 0.5, 0.9 uK_RJ at 30, 70, 100, 353, 857 GHz. (c) The floor that reaches the likelihood is
  the nominal one divided by 20 to 30.
- Order. My hand estimate of the aliased rms in the coarse map before smoothing is 20 to 27 uK_CMB
  for the 5 and 13 arcmin beams, about 12 for 32 arcmin and 5 to 6 for 53 to 60 arcmin. After the
  smoothing I expect the code path to differ from the reference by 1 to 7 uK for the narrow beams
  and by less than 2 uK for the wide ones, uncertain by a factor of three because folded power near
  the band limit is suppressed more than white noise. Between a 13 arcmin and a 60 arcmin band I
  expect a difference of 1 to 5 uK in the code path and below 0.1 uK for smooth-then-degrade.
- Reading rule: a between-band difference above 30 percent of the smaller per-pixel noise rms of
  the pair (floor check, column A) makes the ingest order a first-order error for real-map runs at
  these settings. Below 10 percent it is harmless here.

## Command

```bash
cd /Users/trygvels/code/Commander4
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python -I nextgen/scripts/e034_ingest_checks.py --test order \
  > nextgen/experiments/E034/order_cmb_seed1.txt
.venv/bin/python -I nextgen/scripts/e034_ingest_checks.py --test floor \
  --param nextgen/params/q1/param_q1_n64_ridge_npipetuned.yml > nextgen/experiments/E034/floor_e014.txt
```

## Result

Both runs took 20 s each on one thread (22:44:48 to 22:45:10 and 22:45:21 to 22:45:41). Outputs in
`E034/order_cmb_seed1.txt` and `E034/floor_e014.txt`.

**Floor (real E014 maps, median over abs(b) > 30 deg, uK_RJ).** A is the rms the likelihood got, B
the same without the floor, D the degraded rms with the nominal floor and no smoothing reduction.

| band | E014 table column | D | A (chain) | B (no floor) | A / D | nominal floor | floor in A |
|---|---|---|---|---|---|---|---|
| Haslam 0.408 | 1.3e5 | 1.264e5 | 4.16e4 | 4.16e4 | 0.329 | none | none |
| WMAP K | 5.6 | 5.61 | 1.79 | 1.79 | 0.319 | none | none |
| WMAP Ka, Q, V, W | 4 to 15 | 5.7 to 15.0 | 1.75 to 4.37 | same | 0.29 to 0.31 | none | none |
| Planck 30 | 3.3 | 3.25 | 0.975 | 0.975 | 0.300 | none | none |
| Planck 44 | (in range) | 3.82 | 1.13 | 1.13 | 0.297 | none | none |
| Planck 70 | 3.1 | 3.11 | 0.905 | 0.905 | 0.291 | none | none |
| Planck 100 | 5.1 | 5.75 | 0.375 | 0.313 | 0.065 | 5.65 | 0.205 |
| Planck 353 | 13.4 | 13.44 | 0.494 | 0.083 | 0.037 | 13.44 | 0.487 |
| Planck 857 | 23 | 23.03 | 0.834 | 0.017 | 0.036 | 23.03 | 0.834 |

Expectation (a) holds: D reproduces the E014 column (100 GHz 13 percent off, the rest within 2
percent). Expectation (b) holds: the chain's rms is 3.0 to 3.4 times below the column for bands
without a floor and 15, 27 and 28 times below it for 100, 353 and 857 GHz; the audit's point
predictions (1.0, 0.9, 0.4, 0.5, 0.9) are met to 10 percent. Expectation (c) holds: the floor that
reaches the likelihood is the nominal one divided by 27.6 in all three bands (8 from the
degradation times 3.45 from the smoothing).

**Order (simulated CMB, one realization, uK_CMB; reference map rms 62.7).**

| native beam (arcmin) | aliased rms before smoothing | code minus reference | smooth-then-degrade minus reference |
|---|---|---|---|
| 5.0 | 25.2 | 2.09 | 0.52 |
| 13.2 | 23.5 | 1.90 | 0.52 |
| 32.3 | 16.4 | 1.24 | 0.52 |
| 52.8 | 8.6 | 0.77 | 0.52 |
| 60.0 | 6.4 | 0.68 | 0.52 |

Reference: the 120 arcmin sky with the nside-64 pixel window. Between-band differences in the code
path: 0.79 (13.2 against 32.3 arcmin), 1.40 (13.2 against 52.8), 1.53 (13.2 against 60), 0.64 (32.3
against 52.8), 1.74 (5 against 60). The aliased rms before smoothing matches the hand estimate for
the narrow beams and is 30 to 40 percent above it for 32 and 53 arcmin. The differences after
smoothing are at the lower end of the expected ranges.

## Negative control

Smooth-then-degrade gives between-band differences of 0.0000 uK in every pair, and a common
difference of 0.52 uK from the reference that every band shares and the CMB amplitude absorbs. The
control behaves as it must.

## Conclusion

1. **The model-error floor is confirmed as ineffective (audit B14).** `add_signal_fraction_to_rms`
   adds a white term per native pixel, so at E014 settings it reaches the likelihood 27.6 times
   below its nominal value. The "rms (with floor)" column of the E014 high-latitude table is the
   degraded rms with the nominal floor and without the smoothing reduction. It is not what the
   chains used.
2. **The E014 chains saw the HFI bands as the tightest bands, not as slack ones.** Per-pixel rms
   in the chain: 0.37, 0.49 and 0.83 uK_RJ at 100, 353 and 857 GHz, against 0.9 to 1.1 for LFI and
   1.7 to 4.4 for WMAP. Rescaling the table's high-latitude reduced chi-squares by (D / A)^2 gives
   values of order 13 to 60 in every band (for the bands without a floor the table divided by
   the rms before smoothing, which is 3.0 to 3.4 times the propagated rms the chain used) (K 44, 30 GHz 42, 70 GHz 13, 100 GHz near 20, 353 and 857
   GHz near 27, Haslam 59), where the table shows 0.0 to 6.4. So the model did not fit at high
   latitude at the noise level the chains used, and every "posterior sd" of E014 is the width of a
   fit with a high-latitude reduced chi-square of tens. The rescaling is arithmetic on the table,
   not a new residual computation; the table's residuals come from a plain-average reconstruction.
3. **Degrading before smoothing leaves band-dependent CMB aliasing (audit B1).** For a CMB-only
   sky the bands differ by 0.6 to 1.7 uK_CMB per nside-64 pixel after the common smoothing, and by
   zero in the other order. In uK_RJ the band-dependent part is about 1.1 at 100 GHz against the
   K band (three times the 100 GHz noise of 0.37), 0.7 at 70 GHz against 30 GHz (0.8 times the
   70 GHz noise) and 0.6 between 30 GHz and K (0.65 times the 30 GHz noise). By the reading rule
   fixed in advance (30 percent of the noise) the ingest order is a first-order error for real-map
   runs at these settings. It accounts for a large part of the 100 and 70 GHz misfit of item 2 and
   for little of the K, 30 GHz and Haslam misfit.
4. What stands from E014: the chain agrees with the exact grid target of the model as
   implemented (F2 on real maps), the conditional cannot mix (F14, a statement about the ratio of
   two widths in one error model), and zero levels move the indices (F17, a statement about
   means). What does not stand as written: any number expressed in posterior sd, the reading of
   which bands carry the fit, and the dust index as a measurement (the three HFI bands entered
   with 235 to 760 times the intended weight). The absolute volume-term shift of +0.019 in beta_s
   belongs to this error model: the low-frequency bands that set beta_s carried no floor, so item 1
   does not obviously change it, but it is a number for a model that does not fit and is quoted
   from now on without the "54 sd".

What this does not establish: the size of the aliasing for foregrounds and point sources (the test
is CMB only, so it is a lower bound), and the effect on the indices, which needs the chains rerun
with a corrected ingest.

Changes: F41 added; caveats on F3, F14, F17; the reader fix (smooth at native resolution, then
degrade; floor added after both, as a fraction of the smoothed signal) is queued behind the E032
chains because the working tree is the installed code. Upstream is affected by both items and is
not contacted without Trygve's decision.

## Follow-ups

- Implement the corrected ingest behind a switch, rerun the E014 ridge chain and the grid target,
  and re-read F3 and F17 in a model that fits at high latitude.
- The same order test with a PySM foreground sky and point sources.
- Audit test 9 (pixel covariance of smoothed noise against the stored diagonal, Q7).
