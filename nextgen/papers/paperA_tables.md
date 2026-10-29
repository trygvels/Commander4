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
residual follow what the field can represent. A second noise and template realization
(seed 2028) of the nside-16 rows is running; the table is not final until it is in.
