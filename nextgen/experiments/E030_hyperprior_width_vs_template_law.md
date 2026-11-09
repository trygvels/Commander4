# E030: multiplier hyperprior width against the template-error law on component-separation-quality templates

- **Date started / finished:** 2026-11-08 20:31 / scans read 21:55; chain running
- **Branch and commit:** `nextgen`, 8566e0a (plus the uncommitted Gibbs-check and Paper A edits of the same evening)
- **Open question addressed:** Q16 (robustness of the hierarchical amplitude prior to the template law); feeds the real-data default of D017/E029
- **Validation level:** L1 (exact-target scans of the implemented likelihood, no chain)

## Question
When the amplitude template carries a 30 percent region-scale log-normal error (the cs template
of E025 and E026), which hyperprior width on the per-region template multipliers leaves the
collapsed hier index target unbiased, and does the tight N(1, 0.1^2) that removes the Occam
push with exact templates (F33, F36) survive such a template? If sd 0.1 stays within 0.02 of
the reference on cs templates, it is the real-data default candidate; if it breaks, the width
must be matched to the template error and the matched width's Occam push is the price (F29),
which changes what Paper A recommends for real intensity maps.

## Target and priors
The collapsed per-pixel hier likelihood of Commander4's kernel as scanned by
`nextgen/scripts/e023_exact_scan.py` (F31 machinery): data d_p = M_p(beta) a_p + n_p on the 12
bands 23 to 353 GHz at nside 64 with the E012 white-noise levels; CMB amplitude prior N(0, 1e4^2)
uK (flat in practice, F26); for sync and dust a_p ~ N(mu_r g_p, (s_q g_p)^2) with the template
g_p, one multiplier mu_r per nside-16 region with the hyperprior mu_r ~ N(1, sd^2) integrated
analytically, and one scatter s_q per S/N quartile profiled at the reference (held along a scan,
`--marg-fixed-s 1` default). Dust SED: modified blackbody with the oracle nside-16 region-mean
T map (`mapT`, E028 stage A), so the T misfit is absent and the template effect is isolated.
Index fields held at the reference and shifted block-wise (dust quartiles d_Q1 to d_Q4, sync
level s_all); expected-target mode (noise-averaged chi-square, no realization). Widths sd in
{0.1, 0.3, 0.5}. Templates: `cs` (true amplitude times exp(f), f a Gaussian field with rms 0.3
and FWHM 10 degrees, seed 11, as in E025 and E026) and `exact` at sd 0.3 as the Occam-only
control. The cs template's generative law is log-normal with rms 0.3, so sd 0.3 is the matched
(approximately correct) hyperprior and sd 0.1 is misspecified by a factor 3.

## Setup
Sky and data: `/Users/trygvels/data/e023_pysm_n1` (PySM d1s1c1, nside 64, noise scale 1.0, the
E026 to E029 data). Cells `cs_mapT_perregmarg` and `exact_mapT_perregmarg`; `--mu-sd 0.1 0.3 0.5`
in three runs. Mac mini M4, 8 threads; about 5 to 10 minutes per run. Driver and log under
`/Users/trygvels/data/e023_pysm_n1/e030/` (D027: not in the session scratchpad). Outputs in
`E030/scan_sd{01,03,05}.txt` with the `_refs.npz` references.

## Reference
The noise-free per-region fit of the two indices at the mapT SED with the exact amplitudes
(`refs["mapT"]`, saved to `*_refs.npz`): the pseudo-true region values of this model. A scan
reports the argmax shift of each block relative to it with the curvature sd (0.003 in dust Q1,
0.05 to 0.06 in Q4; sync level about 0.01) and the gain in nats. Baselines in hand: exact
template, mapT, sd 0.1 gives 0.000 in every dust quartile and +0.01 for sync; sd 0.5 gives
+0.01, +0.10, +0.09, +0.07 and +0.08 (`E028/stageA_scan_mapT_vs_fixT_sd01.txt`, `_sd05.txt`).

Pre-registered predictions (written before the runs):
1. exact, sd 0.3: dust optima between the sd 0.1 and sd 0.5 values (a few hundredths in the
   faint quartiles), sync +0.01 to +0.08: the Occam push at the matched width.
2. cs, sd 0.3 (matched law): the same as (1) within 0.02 per quartile (F28, F31: cs minus
   exact at most 0.04 at 1 x noise).
3. cs, sd 0.1 (misspecified by 3x): a shift beyond 0.02 in at least the faint dust quartiles or
   the sync level, of unknown sign (the multipliers cannot follow a 30 percent error inside a
   10 percent prior; the misfit goes into the index or the scatters). If it stays within 0.02
   everywhere, the tight hyperprior is robust to a realistic template and becomes the
   real-data default candidate.
4. cs, sd 0.5: the E025/E026 push (+0.07 to +0.10 dust, about +0.1 sync).

## Command
```bash
/Users/trygvels/data/e023_pysm_n1/e030/run_e030_scans.sh   # three e023_exact_scan.py runs, see the script
```

## Result
**Setup error in passes 1 and 2 (found 21:15).** Both passes ran with the scan script's default
per-pixel CMB prior sd of 100 uK (`--cmb-sd` not given; header line "CMB prior sd 100"), not
the 1e4 uK of the chains, of every E025 to E028 scan and of the Target section above. At
realistic noise a finite per-pixel CMB prior leaks the correlated CMB realization into the
collapsed index target (F26), and the exact-template sd 0.1 cell, which gives 0.000 in every
dust quartile with sd 1e4 (E028 stage A), gives +0.04, +0.12, +0.10 in Q2 to Q4 with sd 100
(`scan2_sd01.txt`): that is the size of the F26 leakage with exact templates and the true T on
this realization, not a property of the hyperprior. Passes 1 and 2 are therefore not readings
of the stated target; their files are kept (`scan_sd*.txt`, `scan2_sd*.txt`) as the F26
illustration and nothing below is concluded from them. Pass 3 (`run_e030_scans3.sh`, launched
21:16) repeats pass 2 with `--cmb-sd 10000`; the scan script now prints a warning when the CMB
prior is finite at realistic noise.

**Pass 1 (20:31 to 20:44; CMB prior sd 100, invalid, see above; `scan_sd03.txt`, `scan_sd01.txt`, `scan_sd05.txt`; the script's default
scan set, per-quartile shifts of both components over [-0.1, 0.1]).** Argmax shift of the
expected target against the mapT reference per S/N quartile (Q1 bright to Q4 faint); "e" marks
an argmax at the scan edge (the optimum is beyond 0.1), with the gain over zero in nats:

| cell | sd | d_Q1 | d_Q2 | d_Q3 | d_Q4 | s_Q1 | s_Q2 | s_Q3 | s_Q4 |
|---|---|---|---|---|---|---|---|---|---|
| exact (E028 stage A) | 0.1 | 0.000 | 0.000 | 0.000 | 0.000 | (s_all +0.01) | | | |
| exact | 0.3 | 0.000 | +0.05 (3.6) | +0.10e (4.8) | +0.06 (0.7) | +0.10e (7.2) | +0.10e (1.2) | -0.09 (0.4) | -0.10e (2.7) |
| exact (E028 stage A) | 0.5 | +0.01 | +0.10e (12) | +0.09 | +0.07 | (s_all +0.08) | | | |
| cs | 0.1 | **-0.020 (24.5; 8 curvature sd)** | +0.09 (21) | +0.10e (18) | +0.10e (8.5) | **+0.10e (31)** | +0.10e (6.6) | +0.10e (3.2) | +0.02 (0.05) |
| cs | 0.3 | 0.000 | +0.04 (2.5) | +0.10e (8.2) | +0.10e (4.6) | +0.10e (11) | +0.10e (4.0) | +0.10e (1.4) | -0.10e (0.7) |
| cs | 0.5 | 0.000 | +0.06 (4.2) | +0.10e (5.9) | +0.10e (2.1) | +0.10e (14) | +0.10e (4.6) | +0.10e (1.2) | -0.10e (1.4) |

Profiled scatters s_q at the reference: exact sd 0.3 dust [0.025, 0.057, 0.108, 0.146], cs sd 0.3
dust [0.067, 0.060, 0.108, 0.144], cs sd 0.1 dust [0.071, 0.146, 0.201, 0.177] and sync [0.077,
0.239, 0.283, 0.409] (the misspecified tight prior pushes the template error into the scatters).

Reading of pass 1 against the predictions (superseded: CMB prior sd 100). (3) appeared confirmed in the strongest form: with the cs
template and sd 0.1 the data-dominated dust quartile is biased by -0.020 (8 curvature sd, 24.5
nats; the only cell in E025 to E030 that moves d_Q1 off zero with the true T), the faint dust
quartiles are pushed +0.09 and beyond, and the bright sync quartile is pushed beyond +0.1 (31
nats). (2) holds for the bright dust quartile (0.000 at sd 0.3 with either template) and the
faint quartiles are within 0.02 of the exact cell where both are inside the scan range (d_Q2
+0.04 against +0.05); the faint quartiles hit the edge. (1): the exact template at sd 0.3 carries
an Occam push of +0.05 to +0.10 in the faint dust quartiles, between sd 0.1 (zero) and sd 0.5
(+0.07 to +0.10), and the sync quartiles move by 0.1 or more in both directions, which the
whole-sky s_all of E028 (+0.01 at sd 0.1, +0.08 at sd 0.5) did not show. Several optima lie
beyond the scan range, so a second pass over [-0.3, 0.3] with the stage A scan set and the
per-quartile sync scans at every width (both templates) was launched at the time stamped in
the log (`run_e030_scans2.sh`) before any conclusion is drawn.

**Pass 3 (21:16 to 21:46; `scan3_sd01.txt`, `scan3_sd03.txt`, `scan3_sd05.txt`; CMB prior sd
1e4, both templates at every width, per-quartile shifts over [-0.3, 0.3] plus the whole-sky
s_all and d_all).** The exact-template cells at sd 0.1 and 0.5 reproduce the E028 stage A scans
to the last digit (same profiled scatters, chi-square and optima), so the machinery is the same
and the pass-1 and pass-2 differences were the CMB prior alone. Argmax shift against the mapT
reference per quartile (Q1 bright to Q4 faint), gain over zero in nats in parentheses where it
exceeds 1; "e" marks an optimum beyond +0.3:

| template | sd | d_Q1 | d_Q2 | d_Q3 | d_Q4 | s_Q1 | s_Q2 | s_Q3 | s_Q4 | s_all |
|---|---|---|---|---|---|---|---|---|---|---|
| exact | 0.1 | 0.000 | 0.000 | 0.000 | 0.000 | +0.01 | 0.00 | 0.00 | 0.00 | +0.01 |
| exact | 0.3 | +0.01 (2.9) | +0.04 (1.9) | +0.02 | +0.02 | +0.03 | +0.05 | +0.04 | +0.03 | +0.03 |
| exact | 0.5 | +0.01 (3.9) | +0.13 (12.8) | +0.09 (2.3) | +0.07 | +0.06 (2.4) | +0.20 (2.6) | +0.17 (1.2) | +0.15 | +0.08 (5.2) |
| cs | 0.1 | **-0.020 (24.2)** | +0.03 (1.9) | +0.08 (2.7) | +0.12 (2.8) | **+0.21 (37.6)** | e (13.2) | e (10.3) | e (8.4) | +0.25 (64.7) |
| cs | 0.3 | 0.000 | +0.03 | +0.08 (1.8) | +0.13 (2.8) | +0.07 (3.0) | e (6.6) | e (6.5) | e (5.5) | +0.12 (12.1) |
| cs | 0.5 | +0.01 (2.1) | +0.10 (7.6) | +0.13 (4.4) | +0.19 (4.6) | +0.07 (3.1) | e (7.4) | e (6.7) | e (5.3) | +0.13 (12.4) |

Curvature sds: dust 0.003 (Q1), 0.017 to 0.026 (Q2), 0.035 to 0.046 (Q3), 0.053 to 0.065 (Q4);
sync 0.024 to 0.029 (Q1), 0.08 to 0.13 (Q2 to Q4), 0.022 to 0.026 (s_all). Profiled scatters
with the cs template at sd 0.1: sync [0.088, 0.279, 0.300, 0.345], dust [0.074, 0.198, 0.280,
0.276] against [0.044, 0.159, 0.191, 0.213] and [0.023, 0.098, 0.182, 0.221] for the exact
template: the tight prior pushes the template error it cannot absorb into the scatters.

Reading against the pre-registered predictions. (1) holds: the exact template at sd 0.3 carries
an Occam push of at most +0.04 (dust) and +0.05 (sync), none above 3 nats, between sd 0.1 (zero)
and sd 0.5 (+0.07 to +0.13 dust, +0.15 to +0.20 sync). (2) holds in the bright half and fails in
the faint half: at sd 0.3 the cs and exact optima agree to 0.01 in dust Q1 and Q2 and in the
sync Q1 to 0.04 (1.5 curvature sd), while dust Q3 and Q4 differ by +0.06 and +0.11 and the
faint sync quartiles run beyond +0.3 against exact's +0.03 to +0.05: at realistic noise the
30 percent template error has an effect of its own in the prior-dominated half, at every width.
This qualifies F28's "cs minus exact at most 0.04 in every quartile at 1 x", which was read on a
grid of +-0.16 with both templates at the edge for the faint sync and against the fixed-T model;
it holds for dust Q1 and Q2 and sync Q1, where the data decide. (3) is confirmed: with the cs
template, sd 0.1 shifts the data-dominated dust quartile by -0.020 (6 curvature sd, 24 nats) and
the data-dominated sync quartile by +0.21 (9 sd, 38 nats); the whole-sky sync level moves +0.25
(65 nats). The misspecified tight prior is not admissible with a 30 percent template. (4) holds:
cs at sd 0.5 gives +0.10, +0.13, +0.19 in the faint dust quartiles (E026's cs_mumap chain sits
+0.06 to +0.10 above bd_trueT there with the field prior) and +0.07 in sync Q1.

The decision-relevant reading: the hyperprior width must be at least the template error's
scale. With the width matched to a 30 percent template the bright dust index is unbiased
(0.000), the bright sync index carries +0.07 (2.5 sd, 3 nats; +0.04 of it the template, +0.03
the Occam push), and the faint half of both fields is set by the template error and the GMRF
prior, not by the data. With an exact template the matched width costs under 0.05 anywhere, so
the width is a property of the template, not of the sky; sd 0.1 is admissible only for a
template good to about 10 percent at the region scale.

**Chain (launched 2026-11-08 21:49, `gmrf16_cs_mumap_nc_Tmap_mu03`, driver and log in
`/Users/trygvels/data/e023_pysm_n1/e030/`):** E028 stage A' with the cs templates and
`mu_prior: [1.0, 0.3]` on both components (the matched width, the oracle T map, non-centred,
1200 iterations, 12 ranks; about 6 h). Pre-registered, written before the chain is read:
(a) validity: the D025 check with the cs templates, the T map and the hyperprior width read
from the parameter file passes on sync and dust (KS p above 0.01); (b) bright dust against
`bd_trueT`: z median within 0.5 (scan optimum 0.000); (c) bright sync against the region-fit
reference: between 0 and +0.09 (the scan optimum +0.07 with its 0.025 curvature sd), to be read
as the template-induced shift of a data-dominated index at the matched width, with stage A'
(-0.015) and E029 (-0.028) as the exact-template baselines; (d) the faint quartiles are read
against the independent-region shrinkage prediction (D027) and reported as prior-dominated, with
the expectation that they lie between the field mean and the scan optima; (e) pixel RMSE against
the PySM truth well below E026's cs_mumap chain at sd 0.5 with the fixed T (sync 0.335, dust
0.089) and above stage A' (0.066, 0.038); (f) the multiplier maps track the realized template
error (exp(-f) region means) with a correlation above 0.5 in the bright quartiles, as in E026.
Control: E026 `gmrf16_cs_mumap_nc` (sd 0.5, fixed T; D025 passed) and stage A' (exact, sd 0.1,
oracle T).

**Chain 1 result (read 2026-11-09 04:35; the chain finished 04:29:57 after 6 h 41 min for 1200
iterations, 0 divergences, 31 gradients per iteration post-burn; check set in
`E030/chain_*.txt`; shrinkage prediction in `chain_quartile_shrinkage.txt` from the new
`nextgen/scripts/quartile_shrinkage.py`, which reproduces the stage A' numbers to four
decimals).** Against the pre-registration:

(a) passes. D025 with the cs templates, the T map and sd 0.3 read from the parameter file: sync
KS p 0.93 (CDF mean 0.508), dust p 0.40 (0.527); conditional optimum minus draw -0.01 and -0.10
conditional sds. The chain samples its target; everything below is a property of the model.

(b) passes. Bright dust against `bd_trueT` +0.004 (z median +0.17, rmse 0.095); the scan
optimum was 0.000. Dust multiplier map mean 1.009 +- 0.005 (exp(-f) averages 1.046 over the
sky; the hyperprior shrinks the map toward 1 where the data are weak).

(c) fails as written, and the criterion was wrongly posed. Sync Q1 against `bs_trueT` +0.19
(against `bs_fixT` +0.16, median +0.15; z median +0.48), outside 0 to +0.09. The region level
resolves it: the top 10 regions by Fisher weight (31 percent of the weight, Fisher sd median
0.10) sit at -0.008 (z median 0.00), the top 30 (50 percent) at +0.058, the top 100 (67
percent) at +0.086, the top 300 (82 percent) at +0.112, and the Fisher-weighted level of the
posterior-mean field is -2.912 against the reference -3.009 (+0.097; against the truth -2.988,
+0.076): the data-dominated part of the sky lands where the scans put it (+0.07 with its 0.025
sd; s_all +0.12). The Q1 average is inflated because the median Q1 region has Fisher sd 1.00
against the learned field sigma 0.51 (median shrinkage weight on the data 0.21): sync Q1 is
prior-dominated as a quartile, and the independent-region shrinkage prediction gives +0.206
against the observed +0.194. Stage A' and E029 showed -0.015 and -0.028 there only because their
exact-template field mean sat near the truth, so the same shrinkage was invisible.

(d) passes. The faint quartiles follow the shrinkage prediction: dust Q2 to Q4 +0.032, +0.059,
+0.023 against +0.033, +0.049, +0.021; sync Q2 to Q4 +0.31, +0.37, +0.32 against +0.29, +0.30,
+0.29 (per-region rms of chain minus prediction 0.09 to 0.12 dust, 0.22 to 0.33 sync; the
prediction ignores the spatial correlation).

The mechanism behind (c) and (d) is new (F38): the field-mean hyperparameter is set by the faint
majority of regions. The sync field mean is -2.701 +- 0.050 against the truth -3.000
(`bs_trueT` mean -2.9995), 6 posterior sd off, exactly where the template-error-dominated faint
half sits (scan optima beyond +0.3 at every width); E026's cs chain at sd 0.5 with the fixed T
had -2.723 +- 0.055. Dust: 1.569 against the `bd_trueT` mean 1.539 (+0.030), smaller because
three dust quartiles have data weights 0.89, 0.29, 0.07 and the faint offsets are +0.02 to
+0.05. With the field mean displaced, the GMRF pulls every region whose Fisher sd exceeds sigma
toward it, which for sync is the whole sky except the hundred brightest regions. With exact
templates (stage A') the field means were -2.967 +- 0.029 and 1.558 +- 0.012.

(e) fails, in the opposite direction. Pixel RMSE against the PySM truth: sync 0.426, dust 0.117
(|b| > 20: 0.459, 0.120) against E026 `cs_mumap_nc` 0.335, 0.089 and stage A' 0.066, 0.038;
per S/N quartile (low to high) sync 0.40, 0.48, 0.46, 0.35 and dust 0.096, 0.132, 0.135, 0.097
against E026's 0.34, 0.40, 0.35, 0.23 and 0.068, 0.105, 0.101, 0.076. Paired against stage A'
by block bootstrap: +0.360 +- 0.021 (sync), +0.064 +- 0.005 (dust). The learned field is
rougher: sigma_s 0.512 +- 0.077 against E026's 0.346 +- 0.080, sigma_d 0.174 +- 0.028 against
0.118 +- 0.032, rho 30 and 47 deg against 59 and 57 deg. GLS CMB residual at the posterior-mean
fields 96.42 uK (|b| > 20: 95.72) against E026's 95.25 and stage A' 93.08; with stage A' dust
and this chain's sync 96.25, with this chain's dust and stage A' sync 93.13, so the excess is
the sync field through its +0.3 faint half. The E026 comparison is confounded three ways: the
T map (oracle here, fixed in E026), the width (0.3 against 0.5), and E026's mixing (its
hyperparameter ESS was 4 to 13, so its sigma and RMSE come from a chain that had not explored
the field prior). Chain 2 below breaks the confound.

(f) passes in Q1, as in E026. Multiplier maps against the realized exp(-f) region means: dust
Q1 corr +0.86 (E026 +0.90), Q2 +0.44 (+0.47), Q3 +0.26, Q4 +0.15; sync Q1 +0.66 (+0.70), Q2
+0.37 (+0.41); rms of map minus a/g 0.275 for both components (E026 0.284, 0.297); posterior sd
of a multiplier median 0.10 (dust Q1) to 0.30 (faint half), i.e. the hyperprior width where
the data say nothing.

Mixing: dust region ESS min 38, median 429; sync min 101, median 692; dust hyperparameters ESS
40 to 659 (`index_field/u` 40, `log_sigma` 55), sync 65 to 395; multiplier maps ESS min 49 and
107.

**Chain 2 (launched 2026-11-09 04:36:00, `gmrf16_cs_mumap_nc_Tmap_mu05`; driver, log and
check set `run_e030_chain2*.sh` in the data dir, outputs `E030/chain2_*.txt`):** cs templates,
`mu_prior: [1.0, 0.5]` on both components, oracle T map, otherwise chain 1's parameter file;
identical to chain 1 except the width and to E026 `cs_mumap_nc` except the T map. Pre-registered
before it is read: (i) D025 passes; (ii) bright dust against `bd_trueT` within 0.5 z median
(scan optimum at sd 0.5 +0.01); (iii) the top-100 Fisher-weight sync regions between +0.04 and
+0.12 (scan +0.07 at sd 0.5), the Q1 average again at the shrinkage prediction; (iv) the sync
field mean again at -2.70 +- 0.06, since the faint half is set by the template at every width
(F37); (v) the discriminating prediction: if the width drives the roughness and RMSE difference,
sigma_s falls to about 0.35, the RMSE to E026's 0.33 and 0.09, and the multiplier maps absorb
more of the error (rms against exp(-f) below 0.275); if the T map or E026's unmixed chain is the
cause, chain 2 lands at chain 1's 0.43 and 0.12 with sigma_s about 0.5, and E026's cs_mumap RMSE
is then read as an unconverged chain's; (vi) the faint dust quartiles between the shrinkage
prediction and the sd-0.5 scan optima (+0.10, +0.13, +0.19).

**Mechanism of the faint-half push, from the pass 3 scan files (read 2026-11-09 04:42):** along
every faint sync scan the expected target's logdet and trace terms cancel to within 0.05 nats
at every offset for both templates (for example cs s_Q4 at +0.3: logdet +3.25, trace -3.22;
exact: +2.08, -2.07), so the push is in the chi2 term alone, the expected data misfit of the
collapsed model. With the exact template chi2 is nearly flat (within 4 nats over +-0.3 in
every faint quartile, mildly favouring flatter indices: -2.3 to -3.9 nats at +0.3 at sd 0.3);
with the cs template it rises monotonically toward flatter indices, +5.4, +6.5, +6.6 nats at
+0.3 in Q4, Q3, Q2 at sd 0.3 and +8.2, +10.0, +12.8 at sd 0.1, i.e. the tighter the multiplier
prior, the larger the misfit a flatter index relieves. So the template error's effect is a
misfit effect, not a volume effect, which is why it survives every hyperprior width (F37). A
candidate reading: with a linear-Gaussian multiplier prior a multiplicative error is absorbed
asymmetrically (a too-low template needs mu 1.35 for f = -0.3, a too-high one 0.74 for +0.3,
and the amplitude prior width s g is tighter where the template is low), so the regions with a
too-low template leave the larger residual; with the sync index information dominated by the
bands above the 30 GHz reference, an under-predicted amplitude is relieved by a flatter index.

**Pass 4 (launched 2026-11-09 04:43, niced beside chain 2; driver `run_e030_scans4.sh` in the data
dir; outputs `E030/scan4_*.txt`):** three scans at the matched width with the expected target,
oracle T map, CMB sd 1e4, all ten per-quartile scans. Pre-registered before they are read:
(1) `scan4_sd03_csneg`: the cs error field with its sign flipped (`--cs-scatter -0.3`). If the
asymmetry reading holds, the faint sync optima stay positive (beyond +0.3); if the push is
linear in f, they flip to beyond -0.3, and the F37/F38 effect would average out over template
realizations. (2) `scan4_sd03_seed12`: a second realization of the error field. Prediction:
faint sync optima beyond +0.3 again, sync Q1 within +0.07 +- 0.05, dust Q1 within 0.01 of zero,
faint dust +0.03 to +0.13 as for seed 11 within a factor two. (3) `scan4_sd30`: a near-flat
multiplier prior (sd 3.0) with the exact and cs templates, the test of Q18 candidate (c)
(let the template say nothing where the data cannot check it). Prediction: with the exact
template the faint sync optima move beyond the sd-0.5 values (+0.15 to +0.20; the Occam push
grows with the width, F29), with the cs template they stay beyond +0.3; i.e. a wide multiplier
prior does not rescue the faint half, and Q18 needs the GMRF on the multipliers or a different
law, not a wider width.

## Negative control
Not needed for the scan machinery (validated in F31 and E028 stage A, where the exact-template
sd 0.1 cell reproduces zero in every quartile); the exact-template cells are the control for
the template effect.

## Conclusion
(scans read 2026-11-08 21:55; chain 1 read 2026-11-09 04:35; chain 2 pending) The hyperprior width of the per-region template
multipliers must be at least the template error's scale: with a 30 percent template the tight
sd 0.1 biases the data-dominated indices (sync +0.21, dust -0.02, tens of nats) and is not
admissible; the matched sd 0.3 leaves the bright dust index unbiased and the bright sync index
+0.07 high, and costs under 0.05 of Occam push when the template is exact. The faint half of
both fields is set by the template error and the GMRF prior at every width (qualifies F28).
FINDINGS: F37, F28 qualified. OPEN_QUESTIONS: Q16 status. The real-data default is a choice
for Trygve (D017 territory): the width follows the template quality.

Chain 1 (cs, sd 0.3, oracle T; D025 passes) confirms the scan reading where the data decide:
bright dust at the reference (+0.004), the Fisher-weighted sync level +0.08 to +0.10 and the
hundred brightest regions at -0.01 to +0.09. It adds F38: the index-field mean hyperparameter is
set by the template-dominated faint majority (-2.70 +- 0.05 against -3.00 for sync), and the
GMRF propagates the displacement into every region fainter than the field's sd, so the sync Q1
average (+0.19) and the whole faint half (+0.31 to +0.37) are the shrinkage prediction, not a
likelihood bias. Pixel RMSE (0.43, 0.12) is worse than E026's sd-0.5 chain (0.34, 0.09), against
prediction (e); the comparison is confounded by the T map, the width and E026's mixing, which
chain 2 separates. With a real-data-quality template the reportable quantities are the
data-dominated regions and the Fisher-weighted level, not the field mean or the region map as a
whole (Q18).

## Follow-ups
- Done: chain 1 (cs, sd 0.3, oracle T) read 2026-11-09; chain 2 (cs, sd 0.5, oracle T) running.
- Q18: a field-mean hyperprior informed by the data-dominated regions only, or a spatially
  correlated (GMRF) multiplier prior that lets the bright neighbours inform a faint region's
  template error (the cs error field has a 10 deg correlation length the iid N(1, sd^2) ignores).
