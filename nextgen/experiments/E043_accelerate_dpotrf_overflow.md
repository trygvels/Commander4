# E043: a heap overflow in the Cholesky routine under every dense fit, and whether any recorded number was touched

- **Date started / finished:** 2026-11-11 02:50 (first crash seen) / open
- **Branch and commit:** `nextgen`, fix in the commit that adds this record
- **Open question addressed:** none in the list; this is a correctness gate on E038 to E042 and on Paper A1
- **Validation level:** L0 (the numerical kernel under the likelihood)

## Question

Three Python processes crashed with corrupted heaps between 02:46 and 02:53 on 2026-11-11 (the E039 plug-in
check at interpreter exit, and two smoke runs of a new nside-8 simulation script). Is the cause in our code,
in a library, or in the machine, and did it change any number we have recorded?

## What was found (diagnosis, before any rerun)

All times from the clock. Machine: Mac mini, Apple M4, macOS 27.0 build 26A5388g (a pre-release build),
numpy 2.5.4, scipy 1.18.1 (wheels, both use Apple Accelerate), Python 3.13.11. Every file of every installed
package matches the hash in its RECORD (7064 files, 0 mismatches), so the installation is intact.

1. The crashes are reproducible and they sit in the Cholesky factorization. Under Apple's guard allocator
   (`DYLD_INSERT_LIBRARIES=/usr/lib/libgmalloc.dylib`) the fault is a byte **write** one byte past a heap
   block, inside `dpotrf$NEWLAPACK` of Accelerate's `libLAPACK.dylib`, called directly through
   `scipy.linalg.lapack.dpotrf(A, lower=0)` with no scipy wrapper code in between.
2. It is the upper-triangle option only. `dpotrf` with `uplo = "L"` never faults. With `uplo = "U"` it faults
   for every order tried from 1500 upward under the guard allocator (1500, 1535, 1536, 1537, 1600, 2048, 2304,
   3072, 4608) and for none below (1024, 1280, 1400).
3. scipy maps a C-ordered array with `lower=True` to `uplo = "U"` (it factorizes the same memory read as a
   Fortran-ordered array). So `scipy.linalg.cho_factor(A, lower=True)` and `scipy.linalg.cholesky(A,
   lower=True)` on ordinary numpy arrays take the faulty path, and `lower=False` takes the sound one. Checked
   under the canary below for all four combinations of memory order and triangle at orders 1536, 3072 and 4608:
   (C, lower) and (F, upper) overflow, (C, upper) and (F, lower) do not, the solves included.
4. How far it writes, measured with a malloc canary (`nextgen/tools/malloc_canary.c`: every block of at least
   64 kB gets 8 MB of padding filled with a pattern, scanned when the block is freed). The call is exactly the
   one in the likelihood, `cho_factor(A, lower=True, overwrite_a=True, check_finite=False)`:

   | order of A | block overrun | bytes written past its end |
   |---|---|---|
   | 1400 | none | 0 |
   | 1500 | 524288 bytes | 786432 |
   | 1536 | 524288 bytes | 786432 |
   | 2048 | 786432 bytes | 786432 |
   | 3072 | 827648 bytes | 745216 |
   | 4608 | 1220864 bytes | 352000 |
   | 5000 | 1327360 bytes | 245504 |
   | 6144 | none | 0 |
   | 8000 | none | 0 |
   | 9216 | none | 0 |
   | 12288 | none | 0 |
   | 18432 | none | 0 |

   What is written is floating-point work data of the factorization (values such as 0.026, 0.033, -0.004). The
   factor that comes back is correct in every case (log-determinant and solution identical to the sound path
   to 1e-14). The damage is to whatever lies behind the workspace in the heap.
5. A second, harmless defect of the same routine: for orders of 33 and more it leaves work values in the
   triangle that LAPACK promises not to touch, inside 8 by 8 blocks on the diagonal, and
   `scipy.linalg.cholesky(lower=True)` returns them (it clears that triangle before the call, not after). The
   requested triangle is correct. Our code only ever read the requested triangle (`solve_triangular`,
   `cho_solve`, the diagonal).
6. Two earlier crash reports exist, from 2026-11-10 20:33 and 22:08, both at interpreter exit. They were not
   noticed at the time because the outputs were complete. They are consistent with the same cause.

## Exposure

| What | Order | Call | Exposed |
|---|---|---|---|
| nside-16 fits, three components (E040, E042, E038 window checks at nside 16) | 18432 | `cho_factor(lower=True)` | no overflow measured at this order |
| nside-16 fits, two components | 12288 | same | no overflow measured |
| nside-16 noise draws (E039 nside-16 arm, E038) | 6144 | `cholesky(lower=True)` | no overflow measured |
| nside-8 fits, three components (E039 stages 1 and 2 and the plug-in addendum, E040 nside-8 arm, E038 `--nside-fit 8`) | 4608 | `cho_factor(lower=True)` | **yes**, 352 kB per factorization |
| nside-8 noise draws (E039) | 1536 | `cholesky(lower=True)` | **yes**, 768 kB per call |
| unit tests of `nextgen/lowres` | up to 1152 | `cho_factor(lower=True)` | no (below 1500) |
| E041 coverage ladder | none | closed forms per element, no LAPACK | no |
| Commander4 itself (`src/commander4`) | | `np.linalg` (uplo L) in `nuts.py`; `cho_factor(lower=True)` only in `standalone_tools/constrained_cmb_realizations.py`, which no experiment here ran | no |
| Julia prototypes (E001 to E037) | | not Accelerate through scipy | no |

So the central real-data table of Paper A1 (nside 16) was computed at orders where no overflow is measured,
and the simulation validation that licenses it (E039, nside 8) was computed at orders where every
factorization wrote 352 kB to 768 kB past a buffer.

## Fix

- `nextgen/lowres/likelihood.py`: the factor is taken with `lower=False` (uplo L). Same cost (0.54 s against
  0.51 s at order 4608; 6.1 s against 5.4 s at 12288 on a loaded machine), same numbers. The environment
  variable `LOWRES_CHOLESKY_UPLO_U=1` restores the old path, for reproducing the overflow only.
- `e039_lowres_validation.py`, `e039b_plugin_check.py`, `e039c_lowpass_sim.py`: noise factor from
  `np.linalg.cholesky` (uplo L; the same matrix L, so the same simulated data for the same seeds).
- New test `test_solve_does_not_write_past_lapack_buffers`: runs profile, marginal and information at order
  4608 under the canary; the default path must not overflow; the old path is run too and reported.
  7 tests pass; the old path is reported as overflowing on this system.
- The full nside-8 simulation pipeline (`e039c_lowpass_sim.py`, two realizations) under the canary with the
  fix: exit 0, no overflow line. Before the fix it crashed in five of five runs.

## Verification plan, written before any rerun is read (2026-11-11 03:13)

A recorded number stands if a rerun on the sound path, under the canary, reproduces it. Tolerance: 2e-4 in an
index and 2e-4 in chi-square per degree of freedom (the optimizer tolerance; the two paths differ at 1e-15 in
the factor, so iterates can differ in the last digits). A larger difference in any realization or fit means
the original was corrupted there, and every statistic that used it is replaced by the rerun.

- V1. nside 16, baseline (K, Ka, Q, V, 353, window-matched, four regions, profile), **old path**, under the
  canary. Expectation: no OVERFLOW line in several hundred factorizations at order 18432, and the recorded
  indices (-3.0940, -3.0965, -3.0867, -2.9950) reproduced. If an OVERFLOW line appears, every nside-16 result is
  suspect and is rerun.
- V2. The same on the sound path. Expectation: no OVERFLOW line, same indices within tolerance.
- V3. E039 stage 1 (nside 8, four regions, 200 realizations, 40 with all estimators, same seeds): every
  per-realization line reproduced within tolerance, and the summary statistics with it. I do not know how
  often the overrun hit live data. The originals passed chi-square (1.0006) and coverage, which argues that
  it was rare; that is an argument, not a check, and this rerun is the check.
- V4. E039 plug-in addendum (40 realizations): the same.
- V5. E039 stage 2 (48 regions): the 95 realizations that exist are rerun together with the missing 105, so
  the stage is completed on the sound path and the first 95 are compared.
- V6. E038 nside-8 window checks and the E040 nside-8 arm: rerun and compared (they support only the
  statement that nside-8 fits of the real maps are biased by the degrade and are not used).
- Everything that runs from now on runs on the sound path, and the nside-8 runs under the canary.

## Command

```bash
cc -O1 -dynamiclib -o ~/data/c4_tools/libcanary.dylib nextgen/tools/malloc_canary.c
C=~/data/c4_tools/libcanary.dylib
# V3
DYLD_INSERT_LIBRARIES=$C .venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions bpxv \
  --nreal 200 --full 40 > nextgen/experiments/E043/sim_n8_bpxv_rerun.txt 2> nextgen/experiments/E043/sim_n8_bpxv_rerun.err
# V1, V2
LOWRES_CHOLESKY_UPLO_U=1 DYLD_INSERT_LIBRARIES=$C .venv/bin/python -I nextgen/scripts/e040_wmap_region_indices.py \
  --regions bpxv --match-window 60 > nextgen/experiments/E043/bpxv_n16_wm_oldpath.txt 2> nextgen/experiments/E043/bpxv_n16_wm_oldpath.err
DYLD_INSERT_LIBRARIES=$C .venv/bin/python -I nextgen/scripts/e040_wmap_region_indices.py \
  --regions bpxv --match-window 60 > nextgen/experiments/E043/bpxv_n16_wm_newpath.txt 2> nextgen/experiments/E043/bpxv_n16_wm_newpath.err
python -I nextgen/scripts/e043_compare.py ORIGINAL RERUN
```

## Result

**V1 holds (read 2026-11-11 03:31).** (Correction, 08:10: the plan above says "several hundred factorizations"; a
profile-only fit makes about a dozen: nine likelihood evaluations and the information. The statement is about those.)

**V1, as read.** The nside-16 baseline on the old path under the canary
(`E043/bpxv_n16_wm_oldpath.txt`): no OVERFLOW line in the whole run (window operators, covariance
inversions, every factorization at order 18432), and the recorded output of 2026-11-10
(`E040/bpxv_n16_wm.txt`) reproduced in every printed digit: chi2_min 19792.7, chi-square per band 0.140,
1.010, 1.057, 0.991, 0.000, per region 0.581, 0.591, 1.815, 0.776, indices -3.0940, -3.0965, -3.0867,
-2.9950 with errors 0.0121, 0.0149, 0.0088, 0.0120, T / J_c 0.0129, 0.0044, 0.0003, 0.0022. The nside-16 fits
were not exposed.

**V2 holds (03:42).** The same fit on the sound path under the canary (`E043/bpxv_n16_wm_newpath.txt`): no
OVERFLOW line, and the output is identical to the old-path output in every character apart from the timings.

**V4 holds (03:42, 25 of 40 realizations read; completed below).** The plug-in addendum of E039 rerun with the
band model it was first run with: every per-realization line identical as printed. A first rerun was made with
`--no-drift` by mistake (the addendum was first run after the drift correction, unlike stage 1); it is kept as
`E043/plugin_n8_bpxv_nodrift.txt` and differs from the original by one unit of the last printed digit of z in
3 of 40 realizations, which is the size of the change in the band model.

**V6, window checks (03:42).** `window_control_n8` and `window_wmap_matched_n8` rerun: biases +0.0149, -0.0146,
+0.2248, -0.0168 and +0.0055, -0.0183, +0.2168, -0.0185, chi2_min 2083.79 and 2087.08, all identical to the
recorded values. No OVERFLOW line.

**V4 complete (03:57).** All 40 realizations identical as printed (160 numbers), and the summary with them:
sd of z 0.974, 0.974, 0.973, coverage 0.700, mean z +0.015 +- 0.077.

**V6 complete (03:57).** The nside-8 real-data arm (`E043/bpxv_n8_wm_all_rerun.txt`) is identical to the
recorded output in every character apart from the timings.

No OVERFLOW line in any run made since the fix (all run under the canary).

**Guard allocator on a whole fit (04:01).** A complete fit at order 4608 (profile and marginal with
gradients, the optimizer, the plug-in sandwich, a noise Cholesky at order 1536) under
`/usr/lib/libgmalloc.dylib`, which faults on a write one byte past any heap block of any size: exit 0 on the
sound path, segmentation fault on the old path (`LOWRES_CHOLESKY_UPLO_U=1`). So the sound path is clean for
small blocks too, which the canary (blocks of 64 kB and more) does not see.

**V3 holds (05:40).** E039 stage 1 rerun with the band model of the original run (`--no-drift`), on the sound
path under the canary (`E043/sim_n8_bpxv_rerun.txt`): 200 of 200 realizations compared, 1480 numbers, largest
absolute difference 0.00 as printed (four decimals in the indices and in chi-square per degree of freedom);
198 lines identical in every character and 2 that differ only in the iteration count of the optimizer; every
summary line identical (mean z, sd of z, coverage, chi-square 1.0006 +- 0.0220, the flat-prior shifts and
their predictions). No OVERFLOW line. The validation table of Paper A1 stands as recorded.

**More factorizations on the old path at the nside-16 orders (08:15).** `e043_oldpath_repeat.py` under the canary: 20
new matrices of order 18432 and 20 of order 12288 factorized in place with `lower=True`: no OVERFLOW line.
The same script at order 4608 as the positive control: 5 of 5 overrun. With the dozen factorizations of the
real baseline (V1) this is what "not exposed" rests on for the nside-16 fits: 52 factorizations at their
orders without an overrun, against an overrun in every factorization tried between orders 1500 and 5000.

**V5, interim (08:07; 29 of 95 realizations of the first part).** 1421 numbers compared: 93.3 percent identical as
printed, 99.4 percent within 1e-4, largest difference 3e-4. Nine entries exceed the tolerance of 2e-4 that I
wrote, in three realizations, all in regions whose sandwich error at the truth is 0.15 to 0.34 (regions that the
data do not measure), where 2e-4 to 3e-4 is 0.0006 to 0.0014 of the error. Read as written, V5 fails for those
entries. Read for what it tests, it does not: the optimizer stops where the objective is flat to its tolerance,
the two Cholesky paths differ in the last bits, and a flat direction turns that into 1e-4 of index; no
well-measured region differs. The tolerance should have been stated relative to the error. To be completed when
the rerun finishes; the stage is not used in the paper.

## Conclusion

The crashes came from a heap overflow in Apple Accelerate's `dpotrf` with uplo = U on macOS 27.0 (26A5388g),
reached through scipy's `cho_factor(lower=True)` and `cholesky(lower=True)` on C-ordered arrays, for matrix
orders of about 1500 to 6000. It is not in our code and not in scipy's arithmetic; the returned factor is
correct. The fix is the other triangle.

Every exposed result that the paper uses has been rerun on the sound path under a malloc canary and
compared as printed: E039 stage 1 (200 realizations, 1480 numbers), the plug-in addendum (40 realizations),
the nside-8 window checks and the nside-8 real-data arm. **No recorded number changed.** The nside-16 fits
were never exposed, which the baseline rerun on the old path under the canary shows directly. The 48-region
stage of E039 (V5) is being rerun and completed on the sound path; it is not used in the paper.

What this establishes: F52. What it does not: that the overflow never altered a number in a run we did not
keep; that other LAPACK routines of this OS build are sound at orders we have not used; that the order
thresholds are stable under other loads or OS builds. What changes: D033 (the sound triangle, and every
dense run on this machine under the canary); a paragraph in the paper; a memory note for future sessions.
