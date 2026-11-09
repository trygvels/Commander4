# E031: the noise regime of the simulation ladder against real maps (intensity and polarization)

- **Date started / finished:** 2026-11-09 21:18 / 2026-11-09 21:30 (run times corrected from the session log at 21:37)
- **Branch and commit:** nextgen, c475106
- **Open question addressed:** Q16 and Q18 (which regime the template-prior results F26 to F38
  apply to), Q9 (what real data can constrain), and the 2026-11-09 review
- **Validation level:** L0 (analytic Fisher information at a known truth; no sampler)

## Question
The ladder from E023 on calls noise scale 1 of `pysm_intensity_sim.py` "realistic noise". Is it?
How does the per-pixel white noise of the two ladder levels (scale 1 and 0.03) compare with the
real WMAP, LFI and HFI maps at nside 64, and what is the per-region (nside 16) Fisher sd of the
spectral indices at real noise, in intensity (the ladder's own sky and model) and in polarization
(the Paper B regime)? The answer decides which findings of F26 to F38 describe real data, and what
Paper B can constrain with and without K band and a low-frequency polarized channel.

## Target and priors
No sampler. The quantity is the Fisher information of the spectral indices at the truth for the
Gaussian white-noise per-pixel likelihood d_pb = sum_c a_pc f_c(nu_b; theta_p) + n_pb,
n_pb ~ N(0, sigma_pb^2), with the amplitudes profiled out under flat priors (Schur complement of
the per-pixel Fisher matrix onto the indices, summed over the pixels of an nside-16 region and
inverted), and the same with the amplitudes known (no Schur term). This is the quantity
`fisher_sd` and `fisher_sd_profiled` of `nextgen/scripts/pysm_intensity_sim.py` report in
`truth.npz`, generalized to per-pixel noise. Models:

- **I3** (the ladder's model): CMB, synchrotron power law (reference 30 GHz), one modified
  blackbody (reference 545 GHz, T_d known). Indices beta_s and beta_d per region.
- **I5**: I3 plus free-free and AME amplitude columns with known per-pixel SEDs (PySM f1 and a1),
  to show what the two missing intensity components cost at real noise. No Haslam, no H-alpha.
- **P** (polarization, per Stokes Q and U, one beta_s per region shared by Q and U): synchrotron
  power law, dust with a known SED (beta_d and T_d fixed at the PySM d1 values), CMB either with a
  flat prior or absent (its E and B signal is far below the noise per nside-64 pixel).

## Setup
- Truth, intensity: `/Users/trygvels/data/e023_pysm_n1/truth.npz` (PySM d1 s1 c1, nside 64, 12
  bands 23 to 353 GHz). Free-free and AME amplitudes and SEDs from PySM f1 and a1 at the same bands.
- Truth, polarization: PySM s1 and d1 Q and U at nside 64 (the s1 polarization template is the
  WMAP nine-year K band smoothed to 5 deg with Gaussian small scales). Cross-check of the
  amplitude with the Cosmoglobe DR1 synchrotron Q and U posterior mean smoothed to 5 deg.
- Real noise: white-noise rms maps of Cosmoglobe DR1 (WMAP K, Ka, Q1, Q2, V1, V2, W1 to W4, LFI 30,
  44, 70; `~/data/q1_inputs/CG_*_IQU_rms_*`), Planck PR3 100 GHz
  (`~/data/planck_pr3/n512/pr3_100_IQU_rms_n0512.fits`) and NPIPE A-split 353 GHz
  (`~/data/beyondplanck/npipe/n512/npipe_353_IQU_rms_n0512.fits`; about 1.4 times the
  full-mission noise). Averaged to nside 64 as sigma_64^2 = sum sigma_i^2 / N^2 over the N
  sub-pixels, differencing assemblies combined by inverse variance, converted to uK_RJ. 143 and
  217 GHz are not on disk (truncated downloads, DATA.md); their noise is set from the Planck 2018
  sensitivities (0.55 and 0.78 uK_CMB deg in temperature, Planck 2018 I Table 4, quoted from
  memory) only if the same table reproduces the local 100 and 353 GHz values (1.29 and 2.56
  uK_CMB deg) to 30 percent; otherwise the two bands are dropped from the real-noise case.
- White noise is the optimistic bound: correlated noise raises the LFI variance at nside 16 by 2.5
  to 2.8 (E022, Q7 side measurement) and WMAP polarization has poorly constrained large-scale
  modes. No beams, no smoothing, as in the ladder.
- Polarized band sets: A (the BP XV set) LFI 30, 44, 70, WMAP Ka, Q, V, 353; B is A plus K band;
  C is B plus one low-frequency channel at 11 GHz with a white Q/U noise of 45 uK_RJ per nside-64
  pixel (my conversion of the 35 to 40 uK per 1 deg beam I remember for the QUIJOTE MFI wide
  survey, 2301.05113; unverified, to be replaced by the value in
  `nextgen/research/2026-11-09_science_landscape_and_data.md`), full sky as an upper bound on
  the gain (the survey covers about 70 percent).
- Machine: Mac mini M4, single thread, seconds.

## Reference
The code must first reproduce the stored `fisher_sd_sync_16`, `fisher_sd_sync_16_profiled`,
`fisher_sd_dust_16` and `fisher_sd_dust_16_profiled` of `truth.npz` at noise scale 1 to
1e-6 relative (same formula, uniform noise). The real-noise per-band medians for WMAP and LFI
must reproduce the table computed inline at 21:15 (below). The polarization Fisher is checked
against E022: summed over the whole sky the profiled sd for set A must be within a factor 3 of
0.01 to 0.02 (E022 posterior sds of four large regions are 0.02 to 0.05 with nside-16 data and
BP's simulated rms, which carries 1.6 times the white sd for LFI).

## Pre-registered expectations (written 2026-11-09 21:18, before any Fisher number at real noise was computed)
Known at this point: the nside-64 white-noise medians for WMAP and LFI intensity (uK_RJ): K 5.47,
LFI 30 3.25, Ka 5.56, Q 5.18, LFI 44 3.80, V 6.67, LFI 70 3.10, W 6.76, against the ladder's 280,
240, 240, 220, 240, 240, 280, 160 at scale 1 (ratios 24 to 90) and 8.4, 7.2, 7.2, 6.6, 7.2, 7.2,
8.4, 4.8 at scale 0.03 (ratios 0.7 to 2.7). Q/U medians: K 7.59, LFI 30 4.57, Ka 7.70, Q 7.04,
LFI 44 5.62, V 9.12, LFI 70 4.37, W 9.25.

1. HFI: the local 100 and 353 GHz maps give about 1.1 and 0.3 uK_RJ per nside-64 pixel, ratios of
   about 90 and 1000 to the ladder's scale 1 (96 and 320), and 3 and 30 to scale 0.03.
2. Intensity, model I3, real noise: the amplitude-known Fisher sd of beta_s per nside-16 region has
   quartile medians of 0.015 to 0.10 (the scale 1 medians 1.05, 2.68, 3.41, 3.86 divided by about
   50), the profiled sd is within a factor 3 of that, and at least 75 percent of regions have a
   profiled sd below 0.1. At scale 1 no quartile median is below 1; at scale 0.03 the medians are
   within a factor 3 of the real-noise ones. Dust: real-noise profiled sd below 0.02 in every
   quartile (scale 1 medians 0.09 to 2.3).
3. Intensity, model I5 (free-free and AME amplitudes free, no Haslam): the profiled sd of beta_s
   grows by more than a factor 5 over I3 at real noise, so that the faint half is again above 0.1.
   This is the weakest of the predictions: I have no prior calculation of it.
4. Polarization, set B (with K band), flat CMB prior: median profiled sd of beta_s per nside-16
   region between 0.5 and 3; fewer than 10 percent of regions below 0.1. Set A (no K band): median
   2 to 4 times larger than B. Whole sky as one scalar, set A: 0.005 to 0.06.
5. Set C (an 11 GHz channel at 45 uK_RJ per pixel): median region sd at least a factor 2 below
   set B.
6. The scale 1 intensity ladder is closer to the real polarization regime than to the real
   intensity regime: the median region sd of beta_s at scale 1 (profiled 6.0, amplitudes known 2.7)
   is within a factor 5 of the polarization set B median.

## Additions after the first run (2026-11-09 21:28, written before the rerun)
The first run (21:21) passed both reference checks and produced every number of items 1 to 6. Two
additions and one repair were decided after seeing it; they are marked as post hoc wherever they
are used.

- Repair: 18 of the 49152 pixels of PySM `a1` have zero emission at 30 GHz, so the AME column
  normalized at 30 GHz was NaN there. The regions holding those pixels came out as inf, and so did
  the I5 whole-sky entry of Table 2b. The rerun gives those pixels the median AME SED of the
  others. No I3 or polarization number is touched.
- Model I5H (post hoc): I5 plus the 408 MHz Haslam map, synchrotron a pure power law from 0.408 to
  353 GHz, free-free and AME with their known SEDs. The Haslam uncertainty is the BeyondPlanck
  effective rms map (`haslam_rms_n0512_60arcmin_uK_1pct_0.8K.fits`, median 1.05 K: 1 percent of the
  signal plus a 0.8 K floor) taken as independent between nside-64 pixels and not averaged down
  within them. That is optimistic for everything that is correlated on large scales (zero level,
  calibration, striping: F17). Expectation, written before computing: the profiled sd of beta_s
  per nside-16 region has quartile medians between 0.01 and 0.1, more than a factor 10 below I5
  without Haslam (first run: 0.62, 1.23, 1.50, 2.20). Reason: Haslam fixes the amplitude at 408 MHz
  to about 7 percent per pixel, so the index error is the relative error of the 30 GHz synchrotron
  amplitude in the five-component fit divided by ln(30/0.408) = 4.3.
- Set C inside the QUIJOTE footprint (post hoc): the same set C numbers restricted to regions with
  every pixel at declination above -30 deg, since the first run's set C row assumed full sky.

## Command
```bash
cd /Users/trygvels/code/Commander4; source .venv/bin/activate
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -I nextgen/scripts/e031_noise_regime_fisher.py --out nextgen/experiments/E031/
```

## Result
Runs: 21:21 (kept as `E031/fisher_report_run1.txt`) and the rerun with the additions at 21:28
(`E031/fisher_report.txt`, 12 s). All numbers below are from the rerun; the two runs agree on every
I3 and polarization number.

**Reference checks.** The stored Fisher values of `truth.npz` are reproduced to 8e-15 relative
(criterion 1e-6). The WMAP and LFI intensity medians reproduce the inline table to the digits
given. The whole-sky polarized sd for set A is 0.018, inside 0.01 to 0.02. All three pass.

**Table 1, white noise per nside-64 pixel, uK_RJ, median over the sky.**

| band | real I | real Q/U | ladder x1 | x1 over real I | ladder x0.03 | x0.03 over real I | x1 over real Q/U |
|---|---|---|---|---|---|---|---|
| WMAP K 23 | 5.47 | 7.85 | 280 | 51 | 8.4 | 1.5 | 36 |
| LFI 30 | 3.25 | 4.59 | 240 | 74 | 7.2 | 2.2 | 52 |
| WMAP Ka 33 | 5.56 | 7.94 | 240 | 43 | 7.2 | 1.3 | 30 |
| WMAP Q 41 | 5.18 | 7.42 | 220 | 42 | 6.6 | 1.3 | 30 |
| LFI 44 | 3.80 | 5.61 | 240 | 63 | 7.2 | 1.9 | 43 |
| WMAP V 61 | 6.67 | 9.52 | 240 | 36 | 7.2 | 1.1 | 25 |
| LFI 70 | 3.10 | 4.38 | 280 | 90 | 8.4 | 2.7 | 64 |
| WMAP W 94 | 6.76 | 9.65 | 160 | 24 | 4.8 | 0.7 | 17 |
| HFI 100 (PR3) | 1.09 | 1.66 | 96 | 88 | 2.9 | 2.6 | 58 |
| HFI 143 (table) | 0.36 | | 64 | 176 | 1.9 | 5.3 | |
| HFI 217 (table) | 0.29 | | 96 | 337 | 2.9 | 10 | |
| HFI 353 (NPIPE A) | 0.29 | 0.75 | 320 | 1090 | 9.6 | 33 | 428 |

The Q/U column is the rms of Q and U together; the Q-only values quoted in the pre-registration
are 0 to 5 percent lower. The local 100 GHz map gives 1.284 uK_CMB deg against 1.29 in the table I
remembered from Planck 2018 I, so the 143 and 217 GHz rows use that table with the spatial pattern
of the 100 GHz map. The 353 GHz map is the NPIPE A split: 3.47 uK_CMB deg against 2.56 for the
full mission, which is the 1.4 expected for half the detectors.

**Table 2, intensity, Fisher sd of beta_s and beta_d per nside-16 region (3072 regions), quartile
medians, Q1 the regions with the smallest scale 1 amplitude-known sd.**

| noise | model | index | amplitudes | Q1 | Q2 | Q3 | Q4 | median | below 0.1 |
|---|---|---|---|---|---|---|---|---|---|
| x1 | I3 | sync | known | 1.00 | 2.15 | 3.17 | 4.58 | 2.68 | 0.002 |
| x1 | I3 | sync | profiled | 2.23 | 4.83 | 7.12 | 10.2 | 5.99 | 0.000 |
| x1 | I3 | dust | profiled | 0.54 | 2.41 | 5.71 | 17.1 | 3.88 | 0.021 |
| x0.03 | I3 | sync | known | 0.030 | 0.065 | 0.095 | 0.137 | 0.081 | 0.66 |
| x0.03 | I3 | sync | profiled | 0.067 | 0.145 | 0.214 | 0.307 | 0.180 | 0.23 |
| x0.03 | I3 | dust | profiled | 0.016 | 0.072 | 0.171 | 0.512 | 0.116 | 0.45 |
| real | I3 | sync | known | 0.017 | 0.039 | 0.056 | 0.077 | 0.048 | 0.97 |
| real | I3 | sync | profiled | 0.029 | 0.065 | 0.093 | 0.126 | 0.078 | 0.70 |
| real | I3 | dust | profiled | 0.0016 | 0.0072 | 0.0175 | 0.048 | 0.012 | 1.00 |
| real | I5 | sync | profiled | 0.62 | 1.23 | 1.50 | 2.19 | 1.36 | 0.002 |
| real | I5 | dust | profiled | 0.0028 | 0.011 | 0.026 | 0.069 | 0.018 | 0.96 |
| real + 408 MHz (post hoc) | I3 | sync | profiled | 0.0023 | 0.0042 | 0.0059 | 0.0082 | 0.0051 | 1.00 |
| real + 408 MHz (post hoc) | I5 | sync | profiled | 0.049 | 0.094 | 0.129 | 0.163 | 0.101 | 0.49 |

The truth sd of the PySM index maps is 0.063 (sync) and 0.037 (dust). Median synchrotron S/N per
nside-64 pixel: 0.31 (ladder x1) against 16.0 (real) at 23 GHz, 0.17 against 12.6 at 30 GHz.
Post hoc decomposition at real noise (profiled median, without and with 408 MHz): three components
0.078 and 0.0051; plus free-free 0.70 and 0.028; plus AME 0.12 and 0.020; plus both 1.36 and 0.10.
Against region size the I3 real-noise median runs 0.078 (nside 16), 0.038 (8), 0.018 (4), 0.0082
(2); I5 runs 1.36, 0.64, 0.31, 0.14.

**Table 3, polarization, profiled Fisher sd of beta_s per region (PySM s1 and d1, real Q/U white
noise, Q and U share the index, dust SED known, flat CMB prior unless stated).**

| set | nside 16 Q1 | Q2 | Q3 | Q4 | median 16 | below 0.1 | median 8 | median 4 | median 2 | whole sky |
|---|---|---|---|---|---|---|---|---|---|---|
| A, BP XV set (no K band) | 0.79 | 1.72 | 2.78 | 5.96 | 2.18 | 0.000 | 1.04 | 0.52 | 0.23 | 0.018 |
| B, A plus K band | 0.27 | 0.61 | 0.99 | 2.12 | 0.77 | 0.011 | 0.37 | 0.18 | 0.087 | 0.0063 |
| C, B plus 11 GHz at 45 uK_RJ | 0.060 | 0.135 | 0.219 | 0.466 | 0.171 | 0.26 | 0.083 | 0.041 | 0.019 | 0.0014 |
| B, no CMB column | 0.15 | 0.34 | 0.55 | 1.16 | 0.43 | 0.048 | | | | 0.0034 |
| B, LFI variance x 2.65 | 0.36 | 0.81 | 1.31 | 2.79 | 1.02 | 0.007 | 0.49 | 0.24 | 0.114 | 0.0083 |

Amplitudes known instead of profiled: set A median 0.84, B 0.36, C 0.065, a factor 2.1 to 2.6
below the profiled values. Inside the declination limit of -30 deg (74 percent of the sky, post
hoc): B 0.79, C 0.17. Fractions of regions with a profiled sd below 0.1 for set B: 0.01 (nside 16),
0.06 (8), 0.24 (4), 0.65 (2), 0.92 (1). Median polarized synchrotron S/N per nside-64 pixel: 2.3 at
23 GHz, 1.9 at 30 GHz (1.8 and 1.5 at |b| > 30). The DR1 synchrotron polarization amplitude over
PySM s1 at 30 GHz, both smoothed to 5 deg, has median 0.91 (16 and 84 percent 0.68 and 1.17), so
the s1 amplitude is representative.

**Verdict on the pre-registered expectations.**

1. HFI: pass. 1.09 and 0.29 uK_RJ, ratios 88 and 1090 to scale 1, 2.6 and 33 to scale 0.03.
2. Intensity I3: mixed. Pass: known quartile medians 0.017 to 0.077 at real noise; profiled a
   factor 1.6 to 1.7 above; no scale 1 quartile median below 1 (1.0005 at the limit); scale 0.03
   medians a factor 1.7 (known) and 2.3 (profiled) above real. Fail: 70 percent of regions, not 75,
   have a profiled sd below 0.1; the dust profiled sd is 0.048 in the faintest quartile, not below
   0.02 (0.0016, 0.0072, 0.0175 in the other three).
3. Intensity I5: pass, and far beyond the predicted factor 5. The median grows by 17, to 1.36, and
   every quartile is above 0.6. Free-free does it (0.70 alone); AME alone gives 0.12.
4. Polarization: all four parts pass. Set B median 0.77; 1.1 percent of regions below 0.1; set A
   median 2.8 times larger; whole-sky set A 0.018.
5. Set C: pass, a factor 4.5 below set B (and the same factor inside the footprint).
6. Mixed. The direction holds: the scale 1 ladder (profiled 6.0) is a factor 7.7 from polarization
   set B and a factor 77 from real intensity I3. The criterion fails: 7.7 is not within a factor 5.
   The ladder is a factor 2.7 from set A, the configuration without K band.
7. Post hoc I5H: the gain over I5 passes (a factor 13); the quartile medians are 0.049, 0.094,
   0.129, 0.163, so the upper bound of 0.1 fails in the faint half.

**Limits of these numbers.** A Fisher sd at the truth measures information; where it is of order
1 or larger it is not a posterior width. White noise only: correlated noise costs LFI a variance
factor 2.5 to 2.8 at nside 16 (row "LFI variance x 2.65": 30 percent on the set B median), WMAP
polarization has poorly measured large-scale modes, and no beam or smoothing is applied. The
free-free and AME SEDs are known exactly, which no real analysis has. The 408 MHz uncertainty is
treated as independent between nside-64 pixels although its zero level and calibration errors are
coherent over the sky (F17). The 11 GHz noise is a remembered number. The polarized dust index is
not scoped: there is no local Q/U noise for 143 and 217 GHz.

## Negative control
The x0.03 rows equal the x1 rows times 0.03 to four digits (2.683 to 0.0805, 5.994 to 0.1798), as
a pure rescaling of uniform noise must. The amplitude-known polarization rows with and without the
CMB column are identical (0.3606), as they must be. The stored scale 1 values are reproduced to
machine precision before any noise is changed.

## Conclusion
Recorded 2026-11-09 21:30.

1. **The label "realistic noise" for scale 1 of the PySM ladder is wrong for intensity and is
   withdrawn.** Scale 1 is 24 to 90 times the white noise of the WMAP and LFI intensity maps at the
   same pixel size and 88 to 1090 times that of HFI. It is also 17 to 64 times the WMAP and LFI
   Q/U noise. The level called "0.03" is the one near real WMAP and LFI depth (0.7 to 2.7), and
   still 2.6 to 33 times above HFI. The levels were inherited from the validation harness, where
   they were chosen to make a fast correctness test, and were never compared with a real map.
2. **No real intensity analysis lives in the regime of scale 1.** With three components at real
   noise the synchrotron index of an nside-16 region is measured to 0.03 to 0.13 and the dust index
   to 0.002 to 0.05, against 2.2 to 10 and 0.5 to 17 in the ladder. Real dust in intensity is data
   dominated in every region, so priors on it are irrelevant and model error is everything.
3. **Real intensity synchrotron is limited by the component model, not by noise.** Free-free with
   a known SED alone raises the region sd from 0.078 to 0.70; with AME it is 1.36. The 408 MHz map
   brings it back to 0.10, at an uncertainty that is systematic and coherent. The three-component
   intensity sky therefore tests neither the noise nor the degeneracy structure of a real intensity
   analysis. This agrees with F17 (zero levels set the real-map indices) and points at Q15.
4. **Real polarized synchrotron is where the weak-data machinery applies.** With K band the region
   sd is 0.77 at nside 16, 0.37 at nside 8, 0.18 at nside 4 and 0.087 at nside 2, against a
   scatter of 0.06 in the PySM s1 index map (the published range of polarized index variations is
   not verified in this record; see the science report). That is prior dominated at fine
   regions, intermediate at nside 4 to 8, and data dominated only for a dozen large regions. One
   11 GHz channel at the assumed depth moves every region size by a factor 4.5. The BP XV band set
   without K band is a factor 2.8 weaker than set B, which is the one real configuration close to
   the scale 1 ladder (a factor 2.7).
5. **Consequence for F26 to F38 and D021, D024, D025.** They stand as statements about a regime
   where the index sd of a region is 30 to 100 times the true scatter. They were established for an
   intensity model whose real counterpart is not in that regime, with a positive-template amplitude
   prior that does not exist for Q and U. Their transfer to polarization is a hypothesis to test,
   not a result. The scale 0.03 results (F23, F24, F27) sit at a region sd of 1 to 5 times the true
   scatter, which is the regime of real WMAP and LFI intensity depth and, in that ratio, of
   polarized synchrotron at nside 4 to 8.

## Follow-ups
- Relabel the ladder in the records: "harness noise" for scale 1 and "WMAP/LFI depth" for 0.03
  (proposed; the paper wording needs Trygve).
- Build the polarized Q/U validation configuration at real noise (band sets A, B, C) with a signed
  hierarchical amplitude prior, and repeat the ladder there. This replaces further intensity
  template work.
- Verify the 11 GHz depth and footprint, and add 13 GHz, C-BASS 5 GHz and S-PASS 2.3 GHz rows, from
  `nextgen/research/2026-11-09_science_landscape_and_data.md` when it exists.
- Scope the polarized dust index once Q/U noise for 143 and 217 GHz is local.
- Intensity: any further work goes to per-band offsets and the 408 MHz anchor (Q15), on a
  five-component sky.
