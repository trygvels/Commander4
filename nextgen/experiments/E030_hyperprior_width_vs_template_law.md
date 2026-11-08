# E030: multiplier hyperprior width against the template-error law on component-separation-quality templates

- **Date started / finished:** 2026-11-08 20:31 / (running)
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

## Negative control
Not needed for the scan machinery (validated in F31 and E028 stage A, where the exact-template
sd 0.1 cell reproduces zero in every quartile); the exact-template cells are the control for
the template effect.

## Conclusion
(pending)

## Follow-ups
- If (3) breaks: one chain with cs templates at the matched width (sd 0.3) and the oracle T
  map, read against the scan optima and E026's cs_mumap chain (sd 0.5).
