# Paper A: draft tables

Numbers are copied from experiment records; every table names its source. Update only from a
record, never by hand.

## Table 2 (candidate): index-prior comparison on the model-consistent sky (E021 attempt 2)

One sky (seed 2027; 12 bands 23 to 353 GHz, nside 64, white noise at 0.03 x the E012 levels;
truth indices region-constant at nside 16, one MBB at fixed T), 1200 Gibbs iterations per run,
700 post-burn draws. Index errors are RMS over the nside-64 pixels of the posterior-mean index
minus the truth ("pixel"), and RMS over regions of the posterior-mean region index minus the
S/N-weighted truth region mean ("region", the linearized pseudo-true value of a region index;
at nside 16 the two coincide). CMB: RMS of the flat-prior GLS CMB at thinned index draws minus
the input, uK_CMB; the GLS at the true SEDs gives 2.79. Gradients: mean NUTS gradient
evaluations per Gibbs iteration (the dominant cost). Coverage: fraction of regions whose z
against the weighted truth mean is below 1 in absolute value (only meaningful at nside 16,
where the model is exact). Hyperparameters: the realized mean and sd of a/g recovered within 2
sd (pass) or not (fail; the test detects index-model misspecification of any kind).
Source: `experiments/E021_attempt2_analyzer_output.txt`, the ladder table in the E021 record
(third stats review), the hyperparameter lines in the E021 Result.

| nside | prior | pixel RMSE sync / dust | region RMSE sync / dust | CMB rms | grad / it | coverage sync / dust | hyper |
|---|---|---|---|---|---|---|---|
| 16 | GMRF | 0.0333 / 0.0180 | 0.0333 / 0.0180 | 2.89 | 125 | 0.78 / 0.78 | pass (within 1.5 sd) |
| 16 | iid | 0.0505 / 0.0275 | 0.0505 / 0.0275 | 2.93 | 29 | 0.78 / 0.69 | fail (mu_d +6.9 sd) |
| 16 | independent | 0.1424 / 0.1009 | 0.1424 / 0.1009 | 3.16 | 31 | 0.68 / 0.69 | pass |
| 8 | GMRF | 0.0379 / 0.0218 | 0.0318 / 0.0153 | 4.22 | 63 | (0.72 / 0.73) | marginal (mu_d +2.1 sd) |
| 8 | iid | 0.0450 / 0.0274 | 0.0403 / 0.0232 | 4.22 | 15 | (0.71 / 0.63) | fail (mu_d +4.5 sd) |
| 8 | independent | 0.0723 / 0.0489 | 0.0696 / 0.0464 | 4.26 | 15 | (0.67 / 0.62) | pass |
| 4 | GMRF | 0.0474 / 0.0286 | 0.0295 / 0.0138 | 5.56 | 39 | (0.60 / 0.55) | fail (mu_d +5.0 sd) |
| 4 | iid | 0.0497 / 0.0319 | 0.0336 / 0.0211 | 5.56 | 11 | (0.56 / 0.43) | fail (mu_d +7.1 sd) |
| 4 | independent | 0.0548 / 0.0327 | 0.0407 / 0.0208 | 5.57 | 15 | (0.58 / 0.56) | fail (mu_d +4.3 sd) |
| 1 region | single index | 0.0669 / 0.0691 | n/a | 9.03 | 4 | n/a | fail (7 to 32 sd) |

Paired differences at nside 16 (iid minus GMRF, block bootstrap over 192 nside-4 blocks, SEs
understated by block correlation): sync +0.0172 +- 0.0013, dust +0.0095 +- 0.0010, CMB +0.054 +-
0.007 uK. Reading (F24): the spatial prior wins at nside 16 and 8; at nside 4 the gaps are
suggestive only. Against the region reference the GMRF loses 13 / 30 percent going from nside 4
to 16 while independent regions lose factors 3.5 / 4.9; the pixel-level error and the CMB
residual follow what the field can represent. Second noise and template realization (seed 2028,
nside-16 rows): GMRF 0.0323 / 0.0180, iid 0.0501 / 0.0271 (paired +0.0178 +- 0.0013, +0.0091 +-
0.0010), independent 0.1373 / 0.1074; CMB 2.86 / 2.91 / 3.13 uK; iid hyperparameter test fails
(mu_d +5.6 sd), the others pass. The realization scatter of the gaps is about 0.0005; the
table's nside-16 rows can be quoted as the mean of the two realizations.

## Table 3 (candidate): the dust temperature model and the multiplier hyperprior at realistic noise (E026, E029, E028)

One PySM d1s1c1 sky (`/Users/trygvels/data/e023_pysm_n1`; 12 bands 23 to 353 GHz, nside 64,
white noise at the E012 levels), exact amplitude templates unless marked cs (the 30 percent
region-scale template error of E025), GMRF index fields at nside 16,
per-region template multipliers with the hyperprior N(1, sd^2), 1200 Gibbs iterations, 700
post-burn draws; every chain passes the D025 Gibbs-consistency check. "Dust minus beta_d" is
the mean over the regions of a dust S/N quartile (Q1 bright to Q4 faint) of the posterior-mean
region index minus `bd_trueT` (the region fit at the true per-pixel T, E027 `refs_nside16.npz`).
mu_d: mean of the posterior-mean multiplier map. Sync: field mean (truth -3.000). RMSE: posterior
mean against the PySM truth over the nside-64 pixels (floors 0.013 sync, 0.019 dust). CMB: the
flat-prior GLS CMB residual rms in uK_CMB (true-SED floor 93.18; not an index-quality criterion,
D027). Exact-marginal optima: the expected-target shift of each dust quartile under the same
model with the multipliers integrated (`E028/stageA_scan_mapT_vs_fixT_sd05.txt`, `_sd01.txt`).
Sources: `E028/stageA_dust_vs_refs.txt`, `stageAprime_dust_vs_refs.txt`, `stageB_dust_vs_refs.txt`,
`E029/dust_vs_refs.txt`, `E028/stageAprime_e020_analysis.txt`, `E028/fixedT_baseline_vs_trueT.txt`.

| chain | sd | T model | dust minus beta_d Q1 / Q2 / Q3 / Q4 | mu_d | sync mean | RMSE sync / dust | CMB rms |
|---|---|---|---|---|---|---|---|
| E026 `gmrf16_exact_mumap_nc` | 0.5 | fixed 20.93 K | +0.006 / +0.061 / +0.059 / +0.030 | 1.029 | -2.868 | 0.147 / 0.060 | 93.60 |
| E028 A `..._Tmap` | 0.5 | oracle map | +0.047 / +0.095 / +0.098 / +0.066 | 1.025 | -2.880 | 0.133 / 0.085 | 92.45 |
| E029 `..._mu01` | 0.1 | fixed 20.93 K | -0.048 / -0.017 / -0.022 / -0.045 | 1.003 | -2.981 | 0.061 / 0.056 | 94.52 |
| E028 A' `..._Tmap_mu01` | 0.1 | oracle map | +0.010 / +0.033 / +0.031 / +0.003 | 1.002 | -2.967 | 0.066 / 0.038 | 93.08 |
| E028 B `..._Tscalar` | 0.1 | sampled scalar (19.19 +- 0.54 K) | -0.004 / +0.029 / +0.024 / -0.002 | 1.002 | -2.977 | 0.062 / 0.046 | 93.42 |
| E026 `gmrf16_cs_mumap_nc`, cs | 0.5 | fixed 20.93 K | -0.008 / +0.060 / +0.080 / +0.050 | 1.030 | -2.723 | 0.335 / 0.089 | 95.25 (at draws) |
| E030 chain 1 `..._cs_mumap_nc_Tmap_mu03`, cs | 0.3 | oracle map | +0.004 / +0.032 / +0.059 / +0.023 | 1.009 | -2.701 | 0.426 / 0.117 | 96.42 |
| exact marginal | 0.5 | fixed | +0.02 / +0.10 / +0.09 / +0.08 | | | | |
| exact marginal | 0.5 | oracle map | +0.01 / +0.10 / +0.09 / +0.07 | | | | |
| exact marginal | 0.1 | fixed | -0.01 / -0.03 / -0.01 / +0.01 | | | | |
| exact marginal | 0.1 | oracle map | 0.000 / 0.000 / 0.000 / 0.000 | | | | |

Reading (F33 to F36): the sd-0.5 rows carry the Occam push of the marginalized multipliers
(+0.07 to +0.10 in the faint quartiles, F29) with either T model; sd 0.1 removes it and, with a
fixed T, exposes the T misfit in the index (E029); with the oracle T map and sd 0.1 the exact
marginal is unbiased in every quartile and the chain's +0.03 in Q2 and Q3 is the GMRF prior's
interpolation (true quartile means 1.570, 1.523, 1.516, 1.548 against a field mean of 1.558;
shrinkage prediction +0.035, +0.042, +0.010 for Q2 to Q4). The scalar-T row agrees with the
oracle-map row to 0.014 per quartile through a 1.7 K lower T (the beta-T degeneracy, F35).
