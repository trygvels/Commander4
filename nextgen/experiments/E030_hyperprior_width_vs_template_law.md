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
(pending)

## Negative control
Not needed for the scan machinery (validated in F31 and E028 stage A, where the exact-template
sd 0.1 cell reproduces zero in every quartile); the exact-template cells are the control for
the template effect.

## Conclusion
(pending)

## Follow-ups
- If (3) breaks: one chain with cs templates at the matched width (sd 0.3) and the oracle T
  map, read against the scan optima and E026's cs_mumap chain (sd 0.5).
