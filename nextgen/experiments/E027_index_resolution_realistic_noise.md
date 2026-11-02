# E027: Index-field resolution at realistic noise with exact templates (gmrf4 and gmrf8 against the E026 gmrf16 rung)

- **Date started / finished:** 2026-11-02 /
- **Branch and commit:** nextgen, 191ea97 (plus the `--ref-npz` option of the scan script)
- **Open question addressed:** Q2 at realistic noise (E021 answered the resolution question at 0.03 times the noise), Q3 / M2.5 (regime boundary for Paper A Sect. 5), Q16 (does the fixed-T absorption of E026 depend on the index resolution)
- **Validation level:** L3 (misspecified simulation with known truth; the sky is not from the model class: indices vary inside regions, T_d varies while the model's T is fixed)

## Question
At realistic noise, with exact templates and the per-region template multiplier field (E026
gmrf16_exact_mumap_nc), does a coarser index field (nside 8, nside 4) recover the sky's indices
better against the truth than nside 16, and does the realistic-noise regime found in E026 (sync
field prior-dominated in every region, aggregate 3 sd high; bright dust +0.06 with the dust mu
map at 1.029) persist, weaken or vanish with region size? The answer places the regime boundary
for Paper A (the resolution the data support at Planck plus WMAP depth) and decides whether the
E026 offsets are properties of the hier model at this depth independent of the index field.

## Target and priors
Identical to E026 gmrf16_exact_mumap_nc (its Target and priors section: CMB amplitude prior
N(0, 1e4^2), hier priors with the exact templates `template_{sync,dust}_exact_*.fits`,
`mu_prior [1, 0.5]` per entry, `log_s_prior [log 0.3, 1]`, `mu_nside: 16` on both hier priors,
GMRF fields with graph-Laplacian basis, alpha 2, sigma_scale 0.3, rho loguniform, field-mean
hyperpriors N(-3.1, 0.3^2) sync and N(1.6, 0.3^2) dust, fixed T 20.9287 K, non-centred
parameterization, dense FD Hessian metric with the analytic mu block, 1200 iterations, warm-up
400, burn 500) with one change per run: `spectral_index_nside` 8 (gmrf8_exact_mumap_nc, 768
regions per component) and 4 (gmrf4_exact_mumap_nc, 192 regions). The template multiplier
field stays at nside 16 in every run so that only the index resolution changes. Parameter
files `param_gmrf8_exact_mumap_nc.yml`, `param_gmrf4_exact_mumap_nc.yml` in the data directory
(diff against the nside-16 file: the two `spectral_index_nside` lines and the output directory).

## Setup
Data `/Users/trygvels/data/e023_pysm_n1` (PySM d1s1c1, nside 64, 12 bands 23 to 353 GHz, white
noise at the E012 levels, seed 2027: the E026 realization). 12 MPI ranks on the Mac mini M4,
runs queued behind the E026 chain-state scan (`run_e027_queue.sh`: gmrf8 then gmrf4). Region
Fisher sds on this sky (`sim_report.txt`, amplitudes known, medians): sync 0.65 (nside 4),
1.33 (8), 2.68 (16) against a sky sd of 0.063; dust 0.096, 0.21, 0.45 against 0.037 (16/84
percentiles 0.013/0.36 at nside 4, 0.084/1.7 at nside 16).

Run log: the queue started gmrf8 at 17:14 on 2026-11-02, was stopped at 17:16 while the E026
chain was tested for sampler bias (the exact-target scan at the chain's state preferred a
sync level 0.12 to 0.18 steeper), and restarted after the Gibbs-consistency check
passed at 17:19 (E026 items 6 and 7; D025). Each rung gets the same check before it is read.

## Reference
1. Pixel-level RMSE and bias against the truth at nside 64 (beta_s, and beta_eff3 for dust at
   the fixed T), as E021 (`e020_analyze.py --tags`), the metric that includes the
   representation error of a coarse field.
2. Noise-free region fits at the run's nside (fixed T, the E023 `region_fit`, started at the
   region mean of beta_s and beta_eff3), per S/N quartile at that nside (`e027_refs.py`,
   written before the analysis): the pseudo-true value of a flat-amplitude region model, as in
   E026. z and RMSE per quartile with `e026_analyze.py --refs`. The quartiles are by the
   region Fisher sd of the index with the true amplitudes, not by the data-built template as
   in E023 and E026 (those carry the noise selection recorded in E026: the faint sync quartile
   sits where the 23 GHz noise is low). At nside 16 the two labellings agree for 51 percent
   (sync) and 66 percent (dust) of the regions, so the nside-16 rung is re-analysed with the
   E027 labels (`E027/refs_nside16.npz`; the fixed-T region fits are identical to E023's) and
   the E026 per-quartile numbers are not compared directly with E027's.
3. The E026 gmrf16_exact_mumap_nc chain as the nside-16 rung (same data, priors and sampler).
4. The Fisher-weighted truth means per quartile for the sync aggregate (E026 check: -2.973 in
   the bright quartile, -2.994 in the three faint ones).

Pre-registered reading:
- (a) Sync: at every rung the region Fisher sd exceeds the sky sd by a factor of 10 or more, so
  the field is expected to stay prior-dominated and near-constant; the quantity to compare is
  the aggregate against the Fisher-weighted truth. If the +0.13 of nside 16 is a property of
  the hier model at this depth (the free multipliers, the fixed-T dust), it should persist at
  nside 8 and 4 within the aggregate's sd; if it comes from the number of weakly constrained
  regions, it should shrink with their number (768, 192). Sharpened after E026 item 7
  (written before any E027 result): the offset is the coupling of the level to the field's
  roughness, and the roughness is set by the prior-dominated field scale (sigma_s 0.15 +- 0.05
  against 0.063). At nside 8 and 4 the region Fisher sds (1.33, 0.65 medians) still exceed the
  sky sd, but fewer regions and more data per region should constrain sigma better; the
  prediction is that the offset of the unweighted mean falls with the posterior sd of sigma,
  and that the roughness check (block-shift optimum at the draw against the smoothed state)
  shrinks accordingly. If sigma stays prior-dominated at nside 4 the fix is an informative
  prior on sigma, not resolution.
- (b) Dust: the bright quartile is data-dominated at every rung; the pixel-level RMSE against
  the truth should fall from nside 16 to 8 in the faint quartiles (less prior shrinkage
  needed) and rise again where the within-region variation of the truth (sd of beta_eff3
  inside a region: computed in `e027_refs.py`) exceeds the region's posterior sd. The rung with
  the smallest pixel-level RMSE is the resolution the data support at this depth.
- (c) If the dust mu map stays at 1.03 and the bright-dust offset at +0.06 at every rung, the
  fixed-T absorption is independent of the index resolution and belongs to the amplitude
  model (D023 caveat, M2.6); if it changes with nside, the index field is part of it.
- (d) Efficiency per rung: hyperparameter ESS, gradients per iteration, wall time, for the
  Paper A table.

## Command
```bash
for t in gmrf8_exact_mumap_nc gmrf4_exact_mumap_nc; do
  mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_$t.yml
done
python nextgen/scripts/e027_refs.py --dir /Users/trygvels/data/e023_pysm_n1 --nside 8 4
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e023_pysm_n1 --burn 500 --tags gmrf8_exact_mumap_nc gmrf4_exact_mumap_nc gmrf16_exact_mumap_nc
python nextgen/scripts/e026_analyze.py --dir /Users/trygvels/data/e023_pysm_n1 --refs nextgen/experiments/E027/refs_nside8.npz --tags gmrf8_exact_mumap_nc
```

## Result
**gmrf8 rung (17:19 to 19:05, 1.75 h, 1200 iterations, 700 post-burn, 0 tracebacks;
`E027/gmrf8_analysis.txt`, `gmrf8_hyper_ess.txt`, `gmrf8_sync_level.txt`,
`e020_analysis_gmrf8_16.txt`; E027 quartile labels for both rungs).** 31 gradients per
iteration (nside 16: 63). Hyperparameter ESS: sync mu 168, log sigma 100, u 38; dust mu 125,
log sigma 29, u 14; region index ESS medians 462 (sync) and 177 (dust), minima 82 and 24.
Fields: sigma_s 0.130 +- 0.043, rho 63 +- 41 deg, mu -2.879 +- 0.040; sigma_d 0.065 +- 0.012,
rho 77 +- 27 deg, mu 1.581 +- 0.014; amplitude hyperparameters mu_s map mean 0.996, mu_d map
mean 1.029 (unchanged from nside 16), s 0.036 / 0.022.
- Pixel-level RMSE against the truth (the metric with the representation error): sync 0.1347
  (nside 16: 0.1467; paired difference +0.0121 +- 0.0012 in favour of nside 8), dust 0.0627
  (nside 16: 0.0597; +0.0022 +- 0.0006 in favour of nside 16); representation floors 0.0235 /
  0.0237 (nside 16: 0.0130 / 0.0187). Per S/N quartile (low to high) sync 0.139, 0.147, 0.144,
  0.106 (nside 16: 0.152, 0.161, 0.156, 0.113); dust 0.039, 0.066, 0.077, 0.062 (0.032, 0.064,
  0.074, 0.061).
- Against the nside-8 region-fit reference (fixed T) per template-S/N quartile: sync z medians
  +0.52, +0.92, +0.99, +0.83, RMSE 0.103, 0.127, 0.133, 0.131, posterior sd medians 0.11 to
  0.14 (region Fisher sds 0.49, 1.06, 1.55, 2.26); dust z +1.44, +1.74, +1.03, +0.10, RMSE 0.069,
  0.087, 0.071, 0.035 (nside 16 with the same labels: 0.072, 0.086, 0.065, 0.031).
- Sync level: the field mean sits +0.12 above the sky mean as at nside 16 (faint quartiles
  +0.11 to +0.12 above the reference; the bright quartile +0.055; the top 10 / 100 regions by
  Fisher weight +0.053 / +0.043 with z medians +0.5 / +0.35); the Fisher-weighted level is
  -2.970 against a truth of -2.999. Pre-registered reading (a): the offset does not shrink with
  the number of index regions (768 against 3072), as the corrected F29 mechanism predicts (the
  Occam push comes from the 3072 multipliers, which are unchanged), and against the
  fewer-regions hypothesis.
- Dust (c): the mu_d map stays at 1.029 and the bright quartile at +1.44 z against the fixed-T
  reference: the fixed-T absorption is independent of the index resolution, as predicted.
- D025 Gibbs-consistency check, K = 100, dust step 0.5 sd (`conditional_check_gmrf8.txt`, 19:10):
  sync PIT mean 0.490, KS p 0.62, optimum minus draw +0.07 sd; dust 0.481, p 0.49, +0.08 sd:
  passes.

**gmrf4 rung (19:05 to 20:15, 1.17 h, 1200 iterations, 700 post-burn, 0 tracebacks;
`E027/gmrf4_analysis.txt`, `gmrf4_hyper_ess.txt`, `gmrf4_sync_level.txt`,
`e020_analysis_gmrf4_8_16.txt`, analysed 20:21).** 30.9 gradients per iteration (nside 8:
31.0, nside 16: 62.9). Hyperparameter ESS: sync mu 269, log sigma 24, u 49; dust mu 53,
log sigma 47, u 64; region index ESS medians 342 (sync) and 141 (dust), minima 37 and 26.
Fields: sigma_s 0.145 +- 0.042, rho 48 +- 34 deg, mu -2.872 +- 0.039; sigma_d 0.060 +- 0.011,
rho 123 +- 32 deg, mu 1.581 +- 0.013; mu_s map mean 0.997, mu_d map mean 1.030, s 0.038 / 0.022
(all three rungs agree on every hyperparameter within one sd).
- Pixel-level RMSE against the truth: sync 0.1430 (nside 8: 0.1347, nside 16: 0.1467; paired
  differences +0.0083 +- 0.0014 for nside 4 minus 8 and +0.0121 +- 0.0012 for 16 minus 8),
  dust 0.0674 (0.0627, 0.0597; +0.0063 +- 0.0014 for 4 minus 16, +0.0022 +- 0.0006 for 8 minus
  16); representation floors 0.0367 / 0.0277 (nside 8: 0.0235 / 0.0237; 16: 0.0130 / 0.0187).
  Per S/N quartile (low to high) sync 0.145, 0.151, 0.152, 0.122; dust 0.046, 0.071, 0.082,
  0.066 (the finer rungs are better in every dust quartile, most in the faintest: 0.046, 0.039,
  0.032 at nside 4, 8, 16).
- GLS CMB residual rms at thinned index draws: 93.45 uK (nside 8: 93.38, nside 16: 93.60; true-SED
  floor 93.18); paired against nside 8: +0.07 +- 0.05 (nside 4) and +0.22 +- 0.05 (nside 16).
- Against the nside-4 region-fit reference (fixed T) per Fisher-sd quartile: sync z medians
  +0.58, +0.87, +0.88, +0.70, RMSE 0.125, 0.139, 0.130, 0.124, posterior sd medians 0.13 to
  0.15 (region Fisher sds 0.22, 0.51, 0.75, 1.08); dust z +1.48, +1.72, +1.17, +0.34, RMSE
  0.058, 0.087, 0.074, 0.048.
- Sync level: the field mean sits +0.13 above the sky mean (nside 8: +0.12, nside 16: +0.13);
  faint quartiles +0.11 to +0.12 above the reference, the bright quartile +0.074, the top 10
  regions by Fisher weight (65 percent of the information at this nside) +0.051 with z median
  +0.5; the Fisher-weighted level of the posterior-mean field is -2.964 against a region-level
  truth of -3.020 and a reference of -3.014. Reading (a): the offset is the same at 192, 768
  and 3072 regions, so it is not a function of the number of weakly constrained index regions;
  the sharpened (roughness) prediction that it falls with the posterior sd of sigma_s is moot
  (that sd is 0.042, 0.043, 0.047 across the rungs, prior-dominated everywhere) and was
  withdrawn with F29's correction; the Occam mechanism (3072 multipliers at every rung)
  predicts exactly this constancy.
- Dust (c): mu_d map 1.030 and bright quartile z +0.34 against the fixed-T reference (nside 8:
  1.029, +0.10; nside 16: 1.029, +0.1): the amplitude-side absorption is independent of nside.
- D025 Gibbs-consistency check, K = 100, dust step 0.5 sd (`conditional_check_gmrf4.txt`, 20:22):
  sync PIT mean 0.478, KS p 0.58, optimum minus draw +0.09 sd; dust 0.516, p 0.23, -0.07 sd:
  passes.
- Analyzer note: the per-quartile lines for the mu maps in `gmrf4_analysis.txt` and
  `gmrf8_analysis.txt` repeat the whole-map numbers (the quartile labels are at the run's nside,
  the mu maps at nside 16); only the whole-map mu-map line is meaningful there.

## Negative control
The comparison is between valid chains of the same target class; the nside-16 rung (E026)
carries the pre-registered reading (a) to (c) as its baseline.

Negative control of the D025 check itself (VALIDATION requirement, run 20:25 to 20:27 on the
gmrf4 chain with K = 100, the same grid as the passing check; `--index-prior-rms` and
`--noise-scale` added to `e026_conditional_check.py` for this): the check evaluated against a
deliberately corrupted target must fail.
- Prior on the field mean narrowed from rms 0.3 to 0.05 (`negctrl_gmrf4_indexprior005.txt`):
  sync PIT mean 0.824, KS p 0.000, optimum minus draw -1.37 conditional sd (predicted from the
  prior pull toward -3.1: about -1.7 sd); dust unchanged (p 0.25, -0.06 sd) because the dust
  field mean 1.58 sits at the prior mean 1.6, so this corruption does not move its conditional.
- Noise rms scaled by 1.5 (`negctrl_gmrf4_noise15.txt`): dust PIT mean 0.249, KS p 0.000,
  optimum +0.83 sd (the Occam and hyperprior terms gain weight against a weaker likelihood);
  sync p 0.033 (the conditional sd widens from 0.018 to 0.026, a mild over-concentration of
  the PIT values).
- cs synchrotron template in place of the exact one (`negctrl_gmrf4_cstemplate.txt`): both
  components fail completely (sync optimum +35 sd from the draw, dust -10.5 sd, KS D 1.0).
A change of the target as small as 1.4 conditional sd along the scanned coordinate is
detected at K = 100; a corruption that leaves the scanned conditional unchanged (the dust case
under the narrowed prior) is not, which is the known limit of a one-coordinate check.

## Conclusion
At the E012 noise depth with exact templates, the index-field resolution between nside 4 and
16 is a second-order choice: it moves the pixel-level index error by 6 to 10 percent, while
the multiplier hyperprior width (F29, +0.12 on the sync level at every rung) and the fixed-T
dust SED (F30, mu_d 1.03 at every rung) set the error budget. The optimum is component
dependent and small: synchrotron is best at nside 8 (0.135 against 0.143 at nside 4 and 0.147
at nside 16, both differences at 6 to 10 sigma by block bootstrap), dust at nside 16, the
finest rung tested (0.060 against 0.063 and 0.067, monotone in every S/N quartile), and the
GLS CMB residual follows the sync field (nside 8 best by 0.07 and 0.22 uK over a 93.2 uK
floor). The pre-registered dust prediction (b), that the faint quartiles would improve from
nside 16 to 8, is refuted: the GMRF prior keeps the finer field competitive even where the
region Fisher sd (0.44 at nside 16) is twelve times the sky sd, because the representation
floor of a coarse region (0.019 at nside 16, 0.028 at nside 4 for dust) is already half the
sky's index sd. The "Fisher sd equals sky sd" crossing used in Q3 is not where the optimum
lies: no rung reaches it for either component (sync 0.65 to 2.7 against 0.063, dust 0.094 to
0.44 against 0.037; extrapolating 1/nside, it sits at nside 1 to 2), yet the finer fields win
or tie. Readings (a) and (c) hold: neither the sync level offset nor the dust multiplier
offset depends on the index resolution, which is what the corrected F29 and F30 mechanisms
predict (both live in the amplitude model, whose 3072 multipliers are the same at every rung).
Cost: 31 gradients per iteration at nside 4 and 8 against 63 at nside 16 (the dense metric's
dimension), wall 1.17, 1.75 h for 1200 iterations against several hours at nside 16 under
load. Recommendation for Paper A: state the resolution rule as "the finest field the budget
allows under the GMRF prior, with the resolution contributing less than 10 percent of the
index error at Planck-like depth", and report the sync level offset and the dust absorption
as amplitude-model properties. All three rungs pass the D025 check (K = 100).

## Follow-ups
