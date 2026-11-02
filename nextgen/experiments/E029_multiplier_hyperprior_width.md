# E029: Does a tight template-multiplier hyperprior remove the faint-sky synchrotron index bias?

- **Date started / finished:** 2026-11-02 /
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
(queued)

## Negative control
The E026 sd-0.5 chain is the control (the bias present).

## Conclusion
(pending)

## Follow-ups
