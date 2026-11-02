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
(running)

## Negative control
None beyond E026's: the comparison is between valid chains of the same target class. The
nside-16 rung (E026) carries the pre-registered reading (a) to (c) as its baseline.

## Conclusion
(pending)

## Follow-ups
