# E038: Gate 0 of the polarized track, are the public low-resolution inputs usable as released?

- **Date started / finished:** 2026-11-09 23:48 / running (items a.1, b and b2 done; a.2 and c open)
- **Branch and commit:** `nextgen` at 0413181
- **Open question addressed:** Q21; Gate 0 of `../reviews/2026-11-09_way_forward_challenge.md` Sect. F
- **Validation level:** L0 (input checks, no inference)

## Question

The proposed K band against 30 GHz null test and everything after it assume that the WMAP
nine-year res 4 Q/U maps and inverse covariances are usable as released, that a trustworthy LFI
30 GHz low-resolution covariance exists publicly, and that WMAP and LFI can be put through one
window. Are they? What changes: if (b) fails, LFI enters with an inflated diagonal and "exact" is
dropped for it; if (c) fails, one window for both sets comes before any index.

## Target and priors

None. Deterministic checks.

## Setup

- (a.1) For each WMAP band: symmetry of N^-1, its spectrum and null space, the pixel rms implied
  by its inverse against the white-noise diagonal, how far the covariance is from white, and the
  chi-square of the map. Script `nextgen/scripts/gate0_wmap_ninv_check.py`. K band fetched
  2026-11-09 23:43 from LAMBDA (`data/map/dr5/skymaps/9yr/raw/` and `raw/ninv/`, 151 MB) into
  `~/data/wmap9_lowres_K/`; Ka, Q, V are in `~/data/wmap9_lowres/`.
- (a.2) E022's conditional check with the dense Ka, Q, V matrices reproduced within 5 percent.
- (b) An LFI 30 GHz dense low-resolution covariance: half-split null test with a PTE between 0.01
  and 0.99, diagonal within 20 percent of BeyondPlanck's simulated rms.
- (c) A beta = -3.1 sky pushed through the WMAP res 4 path and the LFI low-resolution path
  returns region indices within 0.01 with the noise off.

## Reference

I wrote no expectations down before a.1 was run. The only prior statement is the challenge
report's Gate 0 text: N^-1 symmetric to 1e-10 and positive semi-definite, with the processing
mask as its null space. The E031 white-noise estimate of the pixel rms (1.96 uK_RJ per nside-16
pixel for K, 1.99 for Ka) is the other number fixed in advance.

## Command

```bash
cd /Users/trygvels/code/Commander4
for b in "K 22.8 ~/data/wmap9_lowres_K" "Ka 33.0 ~/data/wmap9_lowres"; do set -- $b
  OMP_NUM_THREADS=1 .venv/bin/python -I nextgen/scripts/gate0_wmap_ninv_check.py --label $1 --nu $2 \
    --ninv $3/wmap_band_quninv_r4_9yr_$1_v5.fits --map $3/wmap_band_qumap_r4_9yr_$1_v5.fits \
    > nextgen/experiments/E038/wmap_$1_ninv_check.txt; done
```

## Result

**(a.1), K and Ka, 2026-11-09 23:46 to 23:48.**

| | K (22.8 GHz) | Ka (33.0 GHz) |
|---|---|---|
| Symmetry, max abs(N - N^T) / max abs(N) | 3.8e-10 | 1.0e-10 |
| Negative or null eigenvalues | none | none |
| Condition number | 425 | 298 |
| Pixel rms from the inverse, median (16 and 84 percent), uK_RJ | 2.01 (1.67, 2.42) | 2.04 (1.71, 2.45) |
| The same from the diagonal of N^-1 alone | 1.84 | 1.90 |
| Variance of the worst mode over the median mode | 92 | 65 |
| Second, fifth, tenth, fiftieth worst | 72, 25, 18, 4.9 | 39, 24, 16, 4.6 |
| Share of the total variance in the worst 10 and 100 modes | 0.047, 0.114 | 0.037, 0.101 |
| Mean multipole of the three worst modes (Q part) | 2.9, 6.2, 3.9 | 2.3, 3.9, 4.2 |
| Map chi-square for 6144 elements | 1.16e6 | 1.40e5 |

Against the Gate 0 text: symmetry holds to 1e-10 for Ka and to 4e-10 for K, which is the storage
precision. Both matrices are positive definite. There is no null space: the processing mask is
not projected out in these files, unlike what the Gate 0 text assumed. The pixel rms agrees with the E031 white-noise estimate
(1.96 and 1.99) to 3 percent. The noise is far from white on large scales: a few modes at
multipoles 2 to 6 carry 20 to 90 times the variance of the median mode. A white-noise Fisher
estimate of anything that lives at those multipoles, such as a latitude trend of the index or
the reionization signal, is optimistic by that order.

**(a.1), Q and V, 2026-11-10 00:42 to 00:44.** Q (40.6 GHz): symmetry 8.8e-10, positive definite,
condition 649, pixel rms 2.01 uK_RJ (white diagonal 1.81), worst mode 164 times the median, then
82, 43 (fifth), 17 (tenth); worst ten modes 6.3 percent of the variance; map chi-square 4.3e4. V
(60.8 GHz): symmetry 1.4e-10, positive definite, condition 450, pixel rms 2.24 uK_RJ (2.06), worst
mode 101 times the median, then 54, 29, 21; worst ten modes 4.5 percent; map chi-square 1.3e4.
The pixel rms is within 8 percent of the E031 estimates (1.85 and 2.38). The three worst modes of
every band sit at mean multipoles 2 to 6.

(a.2), (c): not started.

**(b) restated for the data in hand, expectations written 2026-11-10 20:44 before the check was run.** No
public dense low-resolution LFI covariance was reachable (DATA.md, holdings table). What is in
hand since 2026-11-10 20:38 is the official PR4 30, 44 and 70 GHz full-mission maps with their
per-pixel white-noise covariance columns and the two ring-half splits. Check: bring the two ring
halves to nside 16 by the inverse-covariance weighted average used for 353 GHz, and compare the
half difference (d1 - d2) / 2 with the covariance (N1 + N2) / 4 that the map files imply. Script
`nextgen/scripts/e038_lfi_ringhalf_check.py`. A ring-half difference sees white noise and
correlated noise faster than half a pointing period. It does not see errors common to both halves
(gain, bandpass leakage, sidelobes), so a pass is necessary and not sufficient.

- B1. Chi-square per element of the half difference between 1.0 and 1.6 for each band: the
  white-noise columns miss the correlated noise that destriping leaves.
- B2. The excess is on large scales: the ratio of the half-difference power to the white-noise
  expectation is above 1.3 at multipoles 2 to 6 and within 0.15 of 1 at multipoles 30 to 47.
- B3. If both hold, LFI enters the region fits with its white covariance multiplied by the measured
  ratio and is called approximate; the WMAP bands keep their released dense matrices.
- B4 (sign and units check): the nside-16 Q and U maps of LFI 30 GHz and WMAP K correlate at above
  0.8 over the elements with abs(K) above 10 uK, with a ratio near (28.1 / 22.4)^-3.1 times the
  unit factors, 0.50 in thermodynamic units.

**(b) result, 2026-11-10 20:44 to 20:46** (`E038/lfi_ringhalf_030.txt`, `lfi_ringhalf_044_070.txt`; nside
16, 200 white-noise simulations per band for the expected spectra).

| | 30 GHz | 44 GHz | 70 GHz |
|---|---|---|---|
| White-noise sd per element from the map file, uK_CMB | 1.136 | 1.334 | 1.173 |
| Chi-square per element of the half difference (sd 0.018) | 1.771 | 1.401 | 1.209 |
| Power over white expectation, multipoles 2 to 3, EE / BB | 3.3 / 6.9 | 3.7 / 2.1 | 3.7 / 2.9 |
| 4 to 6 | 4.0 / 4.6 | 2.6 / 2.3 | 1.2 / 2.4 |
| 7 to 10 | 2.4 / 2.7 | 2.2 / 1.6 | 1.4 / 1.4 |
| 11 to 20 | 2.4 / 2.3 | 1.6 / 1.5 | 1.2 / 1.1 |
| 21 to 29 | 1.9 / 1.8 | 1.6 / 1.3 | 1.4 / 1.3 |
| 30 to 47 | 1.67 / 1.66 | 1.37 / 1.34 | 1.14 / 1.17 |

- **B1 fails at 30 GHz** (1.77, above the 1.6 I allowed) and holds at 44 and 70 GHz (1.40, 1.21).
- **B2 fails.** The large-scale excess is there (2 to 7 at multipoles 2 to 6), and the excess does
  not vanish at the pixel scale: 1.67, 1.36 and 1.15 at multipoles 30 to 47, where I predicted
  within 0.15 of 1. The white-noise columns of the PR4 map files understate the noise of an
  nside-16 pixel average at every scale, by more at lower frequency, the order of the knee
  frequencies of the radiometers.
- **B3 does not apply as written.** One factor on the white covariance cannot describe a ratio
  that runs from 1.7 to 5 across scales. A scale-dependent model is fitted below.
- **B4 holds for the sign and the units and not for the ratio.** K band against LFI 30 GHz on the
  elements with abs(K) above 10 uK: correlation 0.985 (Q, 1693 elements) and 0.970 (U, 1295),
  unweighted slopes 0.453 and 0.470 against 0.502 for an index of -3.1 with band-averaged spectra.
  Read as an index between 22.4 and 28.1 GHz the slopes are -3.56 and -3.40. An unweighted slope
  on the bright elements gives the poorly measured large-scale modes of the K-band map full
  weight, so this is not a measurement; it is the discrepancy the K band against 30 GHz null test
  is for.

**(b2) a scale-dependent LFI noise model, expectations written 2026-11-10 20:47 before it was fitted.** Model:
N = W^(1/2) R W^(1/2), with W the white covariance of the full-mission map file (2 by 2 per pixel)
and R = r I + Y diag(A l^-gamma) Y^T, Y the spin-2 harmonics of E and B to multipole 47 on the
nside-16 grid, three parameters (r, A, gamma) per band fitted by maximum likelihood to the
whitened ring-half difference. Script `nextgen/scripts/e038_lfi_noise_model.py`, code
`nextgen/lowres/lfi_noise.py`.

- C1. r between 1.5 and 1.8 at 30 GHz, 1.2 and 1.5 at 44 GHz, 1.05 and 1.25 at 70 GHz.
- C2. Under the fitted model, simulated from it 200 times, the quantity sum over a zone of
  z (R^-1 z) for the real half difference lies within 3 sd of the simulations in each of six
  latitude zones and in each of the four BeyondPlanck XV regions, and the binned E and B power of
  the table above within 3 sd in every bin.
- C3. If C2 fails in more than two of the 22 cells per band, the isotropic model is not adequate
  and LFI stays out of the region fits until a model with scan-direction structure is built.
- Limits known in advance: the model is isotropic in the whitened frame where the correlated
  noise follows the scan; it is fitted to one realization; the ring-half difference may overstate
  the noise of the full-mission map (baselines are better determined with all the data) and sees
  nothing common to the two halves.

**(b2) result, 2026-11-10 20:51 to 21:06** (`E038/lfi_noise_model.txt`, parameters in
`E038/lfi_noise_model.npz`; 200 simulations from each fitted model).

| | 30 GHz | 44 GHz | 70 GHz |
|---|---|---|---|
| r | 1.464 | 1.280 | 1.169 |
| A at multipole 10 | 1.161 | 0.448 | 0.164 |
| gamma | 1.265 | 1.246 | 1.692 |
| Improvement in -2 ln L over white | 1368 | 439 | 141 |
| Variance over white at multipoles 2, 5, 10, 20, 40 | 10.4, 4.3, 2.6, 2.0, 1.7 | 4.6, 2.3, 1.7, 1.5, 1.4 | 3.7, 1.7, 1.3, 1.2, 1.2 |
| Cells beyond 3 sd, beyond 2 sd, of 22 | 0, 0 | 0, 1 | 0, 2 |
| Largest deviations | GC +1.9, B 30 to 47 -1.7 | Plane -2.3, zone 0 -1.7 | zone 1 -2.9, E 21 to 29 +2.7 |
| Rms deviation over the 22 cells | 0.98 | 1.03 | 1.20 |

- **C1 fails at 30 GHz by 0.04** (r = 1.464 where I wrote 1.5 to 1.8: the power-law term still adds
  0.2 at multipole 40, and the sum, 1.66, is the measured ratio there) and holds at 44 and 70 GHz.
- **C2 holds for all three bands**: no cell beyond 3 sd. The 70 GHz model is the weakest, with two
  cells near the limit (the latitude zone 10 to 20 degrees at -2.9 and E modes at multipoles 21 to
  29 at +2.7) and an rms deviation of 1.20 where 1.0 +- 0.15 is expected.
- C3 is not triggered. **The LFI bands enter the region fits with N^-1 = W^(-1/2) R^-1 W^(-1/2)**
  (`lowres.lfi_noise.lfi_ninv`) and are called approximate: the model is isotropic in the whitened
  frame, is fitted to one ring-half difference, sees nothing common to the two halves, and may
  overstate the full-mission noise. A noise-scale arm (white covariance times r, and the model with
  A doubled) goes with every result that uses LFI.

Item (b) as first worded, a public dense covariance with a half-split PTE between 0.01 and 0.99,
is not met: there is none in reach. What replaces it is an empirical model that passes the
half-split test it was fitted to in 22 sub-tests per band. Whether that is enough is a judgement
that the K band against 30 GHz comparison has to carry as a stated caveat.

**(c) restated, expectations written 2026-11-10 21:14 before the check was run.** The released WMAP res 4
maps are pixel averages of maps with beams of 0.88 (K), 0.66 (Ka), 0.51 (Q) and 0.35 (V) degrees;
the PR4 maps have 0.54, 0.47 and 0.22 degrees (30, 44, 70 GHz) and 0.08 degrees (353 GHz). The
region fits treat every band as seeing the same nside-16 sky. Check, noise off: a synchrotron
sky with index -3.1 everywhere (the DR1 polarized synchrotron map, 60 arcmin, plus a synthetic
small-scale field that continues its E and B spectra as power laws to multipole 1024, modulated
by the local polarized intensity) and a dust sky (the DR1 polarized dust map, 10 arcmin), each
band smoothed with its own beam, averaged to nside 16 and fitted with the real inverse
covariances as weights. Script `nextgen/scripts/e038_window_check.py`.

- D1. With all bands at a common 1 degree beam (the control) every BP XV region returns -3.1
  within 0.001.
- D2. With the true beams the WMAP-only fit (K, Ka, Q, V, 353) is biased by less than 0.01 in
  every BP XV region, with the largest bias in the Galactic Center or the Plane. Expected sign:
  positive. K band has the widest beam and loses the most power in the bright pixels that carry
  the weight, which reads as a flatter spectrum. (The first wording of this item contradicted
  itself on the sign; corrected before the run, the expectation is positive.)
- D3. Adding LFI 30 GHz changes the bias by less than 0.005.
- D4. For the 48 nside-2 regions the largest bias is below 0.03.
- If D2 fails, the bands are brought to a common window before any index is quoted, or the bias
  is quoted as a systematic per region.

**(c) first results, 2026-11-10 21:21 to 21:32** (`E038/window_control.txt`, `window_wmap.txt`,
`window_wmap_lfi30.txt`). Bias of the fitted index, noise off, four BP XV regions (High Latitude,
Spur, Galactic Center, Plane):

| arm | bias |
|---|---|
| all bands at a common 60 arcmin beam (control) | 0.0000, 0.0000, 0.0000, 0.0000 |
| true beams, K, Ka, Q, V, 353 | +0.0072, +0.0048, +0.0111, +0.0064 |
| true beams, adding LFI 30 GHz | +0.0132, +0.0091, +0.0196, +0.0134 |

- **D1 pass.** **D2: the sign is as expected (positive everywhere), and the size fails in one
  region** (Galactic Center +0.011 against the 0.01 allowed). **D3 fails**: adding 30 GHz raises the
  bias by 0.004 to 0.009 where I allowed 0.005, because the K to 30 GHz lever arm is short and a
  given beam mismatch costs more index.
- The synthetic small-scale field is an upper estimate: the B-mode spectrum of the DR1 map over
  multipoles 20 to 80 fits a power law of index -1.2, which is noise in the map and not sky, and
  the extrapolation carries it to multipole 1024. The arm without small scales gives the lower
  estimate.
- These biases are 0.3 to 1.1 of the statistical error of the WMAP-only fit. They cannot be left
  as they are.

**(c2) window matching, expectation written 2026-11-10 21:34 before it was run.** `nextgen/lowres/window.py`
brings each band to a common 60 arcmin window with S = 1 + Y diag(f_l - 1) Y^+ on the nside-16
grid, f_l the ratio of the beam windows to multipole 47, and transforms N^-1 with S^-1 on both
sides. The map is invertible, so no information is lost (chi-square invariance checked to 12
digits); what it cannot touch is structure above multipole 47.

- D5. With window matching the bias falls below 0.003 in every BP XV region for both band sets.
- D6. If D5 holds, window matching is the baseline of E040 and the unmatched fit an arm; the
  residual bias is quoted as a systematic. If the residual is above 0.005, the cause is found
  first.

## Negative control

None for a.1. For (b): the 200 white-noise simulations, whose binned power scatters by the "white
sd" column of the output files.

## Conclusion

Partial. The WMAP K, Ka, Q and V inverse covariances are usable as released for a dense analysis. The
E031 white-noise numbers for the index per region hold at the pixel scale and are optimistic for
the largest scales. The white-noise columns of the PR4 LFI map files understate the noise of an
nside-16 pixel average at every scale, by 1.2 to 1.8 in chi-square and by 3 to 10 in variance at
the lowest multipoles; a three-parameter scale-dependent model per band fitted to the ring-half
difference passes its own 22 sub-tests and is what LFI enters with. K band and LFI 30 GHz agree in
sign and units, and their unweighted ratio on bright elements is 6 to 10 percent below what an
index of -3.1 gives.

## Follow-ups

- a.2 and c.
- A noise model with scan-direction structure for LFI, or the PR3 low-resolution covariance matrices if the archive serves them again.
