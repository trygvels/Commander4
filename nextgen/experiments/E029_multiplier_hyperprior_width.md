# E029: Does a tight template-multiplier hyperprior remove the faint-sky synchrotron index bias?

- **Date started / finished:** 2026-11-02 / 2026-11-03
- **Branch and commit:** nextgen, e51c4f0 (no code change; parameter file only)
- **Open question addressed:** Q16 (F29 as corrected: the Occam bias of the marginalized per-region multipliers), Q3 at realistic noise
- **Validation level:** L4 (misspecified PySM sky, exact templates)

## Question
The E026 stats review traced the +0.13 of the nside-16 synchrotron field's unweighted mean to
the marginalization volume of the 3072 per-region multipliers under the hyperprior N(1, 0.5^2)
and predicted, by the exact marginal, +0.06 at sd 0.2 and zero within noise at sd 0.1. Does a
chain with `mu_prior [1.0, 0.1]` (both components) show this, with the data-dominated regions
and the dust side unchanged? If yes, the multiplier hyperprior width is the model choice that
sets the faint-sky index at this depth and Paper A states it; if the offset stays, the review's
mechanism is wrong too and the cause is still open.

## Target and priors
E026 gmrf16_exact_mumap_nc (its Target and priors) with one change: the hyperprior on every
template multiplier mu_r is N(1, 0.1^2) for synchrotron and for dust
(`param_gmrf16_exact_mumap_nc_mu01.yml`, output `chains_gmrf16_exact_mumap_nc_mu01/`).
Everything else identical: exact templates, flat CMB prior, GMRF fields at nside 16,
non-centred, dense metric, 1200 iterations, burn 500.

## Setup
Data `/Users/trygvels/data/e023_pysm_n1` (seed 2027, the E026 realization). Queued first in
`run_e028_queue.sh` after E027 (before the E028 stages), 12 ranks on the Mac mini.

## Reference
1. The truth: unweighted sky mean of beta_s -3.000; Fisher-weighted truth per quartile from
   `E026/sync_level_check.txt`; the region-fit reference `bs_fixT`, `bd_trueT`
   (`E027/refs_nside16.npz`, E027 quartile labels).
2. The exact marginal prediction (review `musd.py`): faint-level bias zero within noise at sd
   0.1 (noise realizations -0.08 to +0.09).
3. The E026 chain (sd 0.5) as the control: faint quartiles +0.12 to +0.14 above the reference,
   mu -2.868 +- 0.042; bright quartile at the reference; mu_d map 1.029; bright dust index at
   bd_trueT (z median -0.06).

Pre-registered reading (written before the run):
- (a) Pass: the sync field mean mu within 0.05 of -3.00 and the faint quartiles' posterior-mean
  minus reference within +-0.05 (control +0.12 to +0.14); the top-300 regions by Fisher weight
  within 0.03 of the reference as before. Fail: faint offset above +0.08.
- (b) Dust: the mu_d map mean moves toward 1 (the 1.029 is 0.3 sd of the new prior; expected
  about 1.02 or the T misfit forces it back; either is recorded); the bright dust quartile stays
  at bd_trueT (z median within 0.5). If the dust index moves by more than 0.03 in the bright
  quartile, the multiplier hyperprior also changes what the dust index estimates.
- (c) The D025 check (K = 100, dust step at most 1 sd) passes; sigma_s is reported (the review's
  mechanism predicts no change in the field scale's prior dominance).
- (d) Sync field RMSE against the reference per quartile is reported with the control's.

## Command
```bash
mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_gmrf16_exact_mumap_nc_mu01.yml
python nextgen/scripts/e026_analyze.py --dir /Users/trygvels/data/e023_pysm_n1 --refs nextgen/experiments/E027/refs_nside16.npz --tags gmrf16_exact_mumap_nc_mu01 gmrf16_exact_mumap_nc
python nextgen/scripts/e026_sync_level_check.py --dir /Users/trygvels/data/e023_pysm_n1 --refs nextgen/experiments/E027/refs_nside16.npz --field-run chains_gmrf16_exact_mumap_nc_mu01
```

## Result
Chain 2026-11-02 20:18 to 2026-11-03 01:43 (5.4 h, 1200 iterations, 700 post-burn, 0
tracebacks, 31 gradients per iteration against 63 for the control: the dense metric built at a
smaller multiplier scale). Outputs in `E029/`: `analysis_fixT_ref.txt` (both chains against
the E027 nside-16 references), `sync_level.txt`, `hyper_ess.txt`, `e020_analysis.txt`,
`dust_vs_refs.txt` (`nextgen/scripts/e028_dust_vs_refs.py`, both chains against `bd_trueT` and
`bd_fixT`), `conditional_check.txt` (D025, K = 100, `--mu-sd 0.1` added to the check script).

- **(a) Sync: pass.** Field mean mu -2.981 +- 0.028 (control -2.868 +- 0.042; truth sky mean
  -3.000). Posterior mean minus the region-fit reference per Fisher-sd quartile: -0.028,
  -0.004, +0.014, +0.016 (control +0.12 to +0.14 in the faint three); Fisher-weighted level
  of the posterior-mean field -3.000 (reference -3.009, truth -2.988). The top 10 / 100 /
  300 regions by Fisher weight sit +0.003 / +0.011 / -0.032 from the reference (z medians
  +0.16 / +0.09 / -0.53); the pre-registered 0.03 for the top 300 is missed by 0.002.
  Pixel-level RMSE against the truth 0.061 (control 0.147; paired difference 0.086 +- 0.004
  by block bootstrap); per S/N quartile 0.058, 0.060, 0.059, 0.067 against truth sds 0.051 to
  0.084. The field scale is no longer inflated: sigma_s 0.077 +- 0.040 (control 0.152 +-
  0.047; truth sd 0.063), region posterior sds 0.076 to 0.086 (control 0.13 to 0.16).
  Multiplier map mean 0.9996 +- 0.0018.
- **(b) Dust: the multiplier hyperprior changes what the dust index estimates.** The mu_d map
  mean moves to 1.003 +- 0.002 (control 1.029 +- 0.007): the 3072 multipliers under
  N(1, 0.1^2) no longer absorb the fixed-T misfit at the 545 GHz pivot, so the bright dust
  index moves from beta_d to the fixed-T effective index: bright quartile z median -1.55
  against `bd_trueT` (mean offset -0.048; control -0.06 / +0.006) and +0.11 against `bd_fixT`
  (+0.006; control +1.46 / +0.060). That is a move of 0.054, above the pre-registered 0.03.
  The faint half of the dust sky shrinks to the field mean 1.506 +- 0.011 (control 1.578; sky
  mean of beta_d 1.540, of beta_eff3 1.535): per quartile (bright to faint) the chain sits at
  1.522, 1.505, 1.494, 1.503 against `bd_fixT` 1.515, 1.504, 1.520, 1.585 and `bd_trueT`
  1.570, 1.523, 1.516, 1.548, so the faintest quartile is -0.083 below the fixed-T reference
  (z -1.27) and -0.045 below the true-T one. Pixel-level dust RMSE 0.056 against beta_d and
  0.059 against beta_eff3 (control 0.060 / 0.065); per S/N quartile (low to high) 0.054,
  0.040, 0.045, 0.078 (control 0.032, 0.064, 0.074, 0.061): better in the middle, worse at
  both ends.
- **CMB.** The flat-prior GLS CMB residual at thinned index draws is worse with the tight
  hyperprior: 94.52 uK against 93.60 (paired +0.91 +- 0.07; true-SED floor 93.18), and the
  high-latitude residual 93.95 against 92.91. The excess over the floor rises from 8.9 to
  15.8 uK rms, with the power fraction at l 30 to 100 up from 0.865 to 0.883. The sync field,
  which improved by a factor 2.4, does not carry this; the faint-sky dust index does: an index
  0.08 low over the high-latitude sky mis-extrapolates the 353 GHz dust to 143 GHz by 7
  percent, which is a uK-level residual there.
- **(c) D025 check passes** (K = 100, steps 0.6 and 0.5 conditional sd): sync PIT mean 0.470,
  KS p 0.63, optimum minus draw +0.16 sd; dust 0.485, p 0.78, +0.04 sd. Mixing: sync mu ESS
  201, log sigma 32, u 38; dust mu 68, log sigma 8, u 7 (the dust field-scale hyperparameters
  mix poorly; sigma_d 0.069 +- 0.013 and rho 42 +- 16 deg agree with the control's 0.062 and
  40, so the dust shrinkage is not an artifact of the scale estimate); region ESS medians 215
  (sync) and 139 (dust).
- **(d)** Sync RMSE against the reference per quartile 0.109, 0.075, 0.059, 0.055 (control
  0.102, 0.140, 0.152, 0.146); the bright quartile is unchanged (its RMSE against the reference
  exceeds that against the region-mean truth, 0.066, because the fixed-T region fit leaks the
  dust misfit into the bright sync regions).

## Negative control
The E026 sd-0.5 chain is the control (the bias present).

## Conclusion
The review's mechanism is confirmed in a chain: with the multiplier hyperprior at N(1, 0.1^2)
the faint-sky synchrotron index bias of E026 is gone (field mean -2.981 +- 0.028 against
-3.000; faint quartiles within 0.02 of the reference against +0.12 to +0.14), the field scale
returns to the sky's (0.077 against 0.152, truth 0.063) and the sync pixel error drops from
0.147 to 0.061. So the hyperprior width on the per-region multipliers is the model choice that
sets the faint-sky index at this depth, and the wide default also inflated the field scale.
The price lands on dust: without free multipliers the fixed-T misfit is carried by the dust
index again (bright sky at beta_eff, 0.05 below beta_d) and, through the GMRF shrinkage, by
the faint sky (0.08 below the fixed-T reference), which costs 0.9 uK of GLS CMB residual at
high latitude. Neither hyperprior gives an unbiased dust index with a fixed T: sd 0.5 hides
the T misfit in the multipliers and adds the Occam push, sd 0.1 exposes it. The multiplier
field is therefore not a substitute for a dust temperature model, and the stage A chain
(oracle T map, running) is the direct test of whether the tight hyperprior plus a correct T
gives both indices right. Recommendation (needs Trygve, a default change): multiplier
hyperprior sd 0.1 for the hier amplitude prior, with a T model for dust; Paper A states the
hyperprior width as part of the model and shows the two chains side by side.

## Follow-ups
- Stage A (oracle T, sd 0.1) reading: bright dust at `bd_trueT` and the faint sky at the
  truth would close the dust side; if the faint-sky shrinkage stays, it is the GMRF field
  mean set by the bright sky and needs a region-dependent reference.
- The dust hyperparameter mixing (ESS 7 to 8 for log sigma and u) wants the dense metric
  rebuilt at the chain's scale or more iterations before the dust field scale is quoted.
- The CMB residual's sensitivity to the faint-sky dust index is a Paper A point: the index
  error that matters for the CMB is the one at high latitude, not the bright-sky one.
