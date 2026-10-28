# E021: Index-field prior comparison on a sky drawn from the model class: hierarchical GMRF with inferred smoothness against iid regions, independent regions and a single index (M2.5 first cut, Q3; redesign of E020)

- **Date started / finished:** 2026-10-28 /
- **Branch and commit:** nextgen, e2c8e3b
- **Open question addressed:** Q3 (optimal spatially varying index model); ROADMAP M2.5
- **Validation level:** L3 (sky drawn from the model class with known truth; the samplers are
  validated at L1 to L3 in E018 and the target was verified on this data by direct computation
  in E020)

## Question
The E020 question, asked first on a sky the model can represent exactly: on an intensity sky at
nside 64 with 12 bands and white noise at 0.03 x the E012 levels, whose synchrotron and dust
indices are constant over nside-16 regions, whose dust is one modified blackbody at a known T,
and whose amplitudes scatter about the templates exactly as the hierarchical amplitude prior
assumes, does the hierarchical GMRF index field with inferred mean, scale and correlation length
recover the index maps and the CMB better than (a) exchangeable (iid) hierarchical regions, (b)
independent regions at a fixed nside, (c) a single full-sky index, and is the GMRF answer stable
when nside_beta is refined from 4 through 8 to 16? If the GMRF wins and is stable, it is the
Phase 2 index model and the remaining Q3 variants are compared against it; if the fixed ladder is
as good, the field prior buys nothing at this S/N. E020 could not answer this because the PySM
sky's misspecifications dominated the posterior (F21); this is the clean comparison the priors
need before the robustness question (Q16).

## Target and priors
Identical to E020 (its "Target and priors" section applies word for word: the collapsed hier
likelihood of E018 over the three per-pixel amplitudes, the amplitude priors a_cmb ~ N(0, 100^2),
a_c ~ N(mu_c g_c, (s_c g_c)^2) with hyperpriors mu_c ~ N(1, 0.5^2), log s_c ~ N(log 0.3, 1^2),
guard log s_c >= -8, and the four index-prior options at nside_beta 4, 8, 16 and 0 with mu0 =
(-3.10, 1.60), sd 0.3, independent-region sd 1.0). What changes is the data, so that the model is
well specified:

Sky (writer `nextgen/scripts/pysm_intensity_sim.py --sky consistent --nside-truth 16
--template-scatter 0.3`, seed 2027, same bands, nside and noise as E020):
- beta_s(p) and beta_d(p) are the plain means of the PySM s1 and d1 index maps (degraded to
  nside 64) over each nside-16 region, broadcast back to nside 64: 3072 distinct values each.
  Truth sd 0.0615 (sync) and 0.0323 (dust) at nside 64; ranges -3.267 to -2.761 and 1.430 to 1.727.
  So the nside_beta 16 options are exactly specified, nside_beta 8 and 4 average over 4 and 16
  true regions (their RMSE floors are the truth's own region means, reported by the analyzer),
  and the single index averages over all 3072.
- Dust: one modified blackbody at T = 20.929 K everywhere (the mean of the d1 temperature map,
  the value the model fixes), amplitude a_dust(p) = the PySM d1 emission at 545 GHz in uK_RJ.
  Hence beta_eff = beta_eff3 = beta_d; the fixed-T misspecification of E020 is gone by
  construction.
- Synchrotron: a_sync(p) (nu/30)^beta_s(p) with a_sync the PySM s1 emission at 30 GHz.
- CMB: the PySM c1 map (rms 84.1 uK_CMB) times the thermodynamic-to-RJ factor, as before.
- Noise: white, rms per band 8.4, 7.2, 7.2, 6.6, 7.2, 7.2, 8.4, 4.8, 2.88, 1.92, 2.88, 9.6 uK_RJ per
  nside-64 pixel (0.03 x the E012 nside-16 levels scaled to nside 64), one realization.
- Templates: g_c(p) = a_c(p) / r_p with r_p = max(1 + 0.3 eps_p, 0.2), eps iid N(0, 1) per pixel
  and component (the 0.2 floor is hit by 0.34 and 0.42 percent of pixels; g/true at |b| > 20,
  16/50/84 percentiles: 0.771/1.001/1.430 for both). So a_c = g_c (1 + 0.3 eps) up to the floor:
  the hier prior a ~ N(mu g, (s g)^2) is correctly specified with mu = 1 and s = 0.3, which the
  sampled hyperparameters must recover. The targets are the realized values, not 1 and 0.3: the
  floor and the finite sky give mean(a/g) 0.9990 (sync) / 0.9976 (dust) and sd(a/g) 0.2986 /
  0.2988 for this seed (analytic with the floor: 1.00035 and 0.29895; the full-information
  posterior sd is about 0.0014 for mu and 0.0010 for s, so the dust mu sits 1.8 sd from 1 and a
  check against 1 would fail with everything correct; review item 3). Acceptance: posterior mean
  of mu_c within 2 posterior sd of the realized mean(a/g), posterior mean of s_c within 2 sd of
  the realized sd(a/g), the same in all ten runs. Unlike E020 the templates are built from the
  truth: this is deliberate (L3); the data-built templates return in Q16.

Per-region Fisher sd (profiled over the three per-pixel amplitudes; the writer's
`sim_report.txt`): nside_beta 4 median 0.041 sync / 0.024 dust (22 / 40 percent of regions above
the truth sd); nside_beta 8 0.087 / 0.054 (69 / 64 percent above); nside_beta 16 0.180 / 0.116
(89 / 78 percent above); full sky 0.0005 / 0.00002. The same regime as E020: nside_beta 4
data-dominated, 8 marginal, 16 prior-dominated for most regions, so the comparison is a
comparison of shrinkage.

Verification of the target on this data before any run (offline evaluation of
`RegionIndexLikelihood.from_arrays` on the band files, scratchpad `probe_e021_like.py`): the
log posterior at the truth state is 53,859 above the flat prior-mean start and 519,360 above a
flat beta_s = -4.0; 1-D scans of one region's beta_s around the truth are symmetric with local sd
0.014 to 0.19, between the amplitudes-known and profiled Fisher sd; the analytic gradient matches
central finite differences in every block (ratios 1.0000 to 1.0029, the latter on a 4.4-unit
gradient with a 1e-5 step). The gradient norm at the truth is 2.8e4 against 2.3e6 on the E020
sky, dominated by the full-sky dust mean (21,886 on a coordinate with Fisher sd 2e-5, a shift of
1e-5): the truth is at the mode to within the noise.

## Setup
Writer, parameter files, runs, order, sampler and analyzer as in E020 with these settings:
`parameterization: centred`, `init_inv_mass: fd_hessian`, `adapt_mass: true` (D014: the windows
were cleared of blame in E020 attempt 3; a 40-iteration timing run on this sky,
`chains_timing_gmrf16/`, measures gradients per transition before and after the first window and
decides; if the windows degrade the metric, `adapt_mass: false` is used for all ten runs, a
sampler setting that leaves the target unchanged). Iteration count: set from the timing run so
that the ten runs fit in about two days on the Mac mini (10 cores, the SBC replication running
alongside), with the same warm-up (500) and burn (600) as E020 if 2000 iterations fit, else a
shorter chain with warm-up 500 and burn 600 and the count recorded here before launch. 12 MPI
ranks, Commander4 always runs two chains; the analyzer reads chain 1.

Metrics: E020's (`nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent`),
with three additions from the design review (2026-10-28, below), all in the analyzer before any
production run:
- Cross-nside comparisons use the region-level error sqrt(RMSE^2 - floor^2) with the floor the
  RMSE of the truth's own unweighted region means at that nside_beta (an exact decomposition of
  the pixel RMSE; on this truth the floors are 0.0343 / 0.0195 / 0 for sync and 0.0204 / 0.0146
  / 0 for dust at nside_beta 4 / 8 / 16, so the pixel RMSE favours 16 by construction), and the
  GLS CMB residual is compared against the GLS at the truth's region means per nside_beta
  (`gls_at_region_means_nside{4,8,16}`), not only against the true-index GLS. Rule (c) of E020
  ("8 and 16 tie") is read on the region-level error and on the CMB residual relative to its
  per-nside reference; on this truth (constant per nside-16 region) a pixel-RMSE tie between 8
  and 16 is impossible and is not asked.
- The Monte Carlo error of the posterior means, sqrt(mean over pixels of sd_r^2 / ESS_r), is
  reported per run as an RMSE inflation (it adds in quadrature). Runs are compared only when
  this is below a quarter of the RMSE difference under test; otherwise the chain is extended.
- The hyperparameter check above (realized mean and sd of a/g, 2 posterior sd).

## Design review
stats-reviewer 2026-10-28 (independent recomputation from the files in the data directory and the
code): the sky is from the model class, CONFIRMED on all four points (indices constant per
`region_map(16, 64)` region and equal to the PySM region means; `_dust_sed` textually equal to
`ThermalDust.get_sed`, constants differing by 3.3e-7 relative and T by 2.4e-5 K, at most 0.017
uK_RJ at 353 GHz in the plane, 0.2 percent of the noise sd; CMB conversion equal to 1.4e-6; sync
pivot 30 GHz, a_sync equal to the PySM emission; noise z-scores sd 0.99996; no template floor
active). Must-fix taken: the negative control was stated against the unweighted mean, which a
correct single index does not reach (rewritten above). Should-fix taken: region-level error and
per-nside GLS references for the cross-nside rules (2); realized hyperparameter targets with a
2-sd band (3); Monte Carlo error of the posterior means reported as an RMSE inflation (4).
Confirmed: the analyzer has no division by beta_eff - beta_d and does not use T_d (5); templates,
hyperpriors, Gibbs seed, noise realization and NUTS settings are shared by all ten runs (6).

## Reference
The truth maps in truth.npz (region-constant indices, the c1 CMB, the amplitudes). No exact
posterior exists at this size; sampler correctness is from E018 (L1 to L3) and the direct
verification above. The E020 chain outputs on the misspecified sky are kept as the contrast.

## Command
```bash
python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e021_consistent --nside 64 --noise-scale 0.03 --seed 2027 --sky consistent --nside-truth 16 --template-scatter 0.3
mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent/param_timing_gmrf16.yml    # 40 iterations
for t in gmrf16 iid16 indep16 scalar gmrf8 iid8 indep8 gmrf4 iid4 indep4; do mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent/param_$t.yml; done
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent --burn 600
```

## Result
**Timing run, gmrf16, 40 iterations, centred, fd_hessian start, adapt_mass true (05:31 to
06:28):** Hessian start 14 min per chain (12304 gradients, range 1.4e-8 to 2). The chain
converges towards the truth (sync index RMSE 0.117 at the start, 0.076 at iteration 40 with the
field sd 0.076 against a truth sd 0.061; dust 0.068 to 0.060; sync log sigma -4.9 to -2.5, u
wandering between 8 and 0) but at 255 to 1023 gradients per transition with step sizes 0.004 to
0.017 and acceptance near 1.0: 66 s per Gibbs iteration (both chains), 37 h per run.
Offline diagnosis (scratchpad `diag_e021_nuts.py`, the same likelihood from the band files, NUTS
from the iteration-40 state, fixed metric, 30 transitions each): with a fresh finite-difference
Hessian at that state the sampler takes 136 gradients per transition at step sizes 0.05 to 0.16
(all 6152 coordinates free); with the field hyperparameters log sigma and u additionally held
fixed, 187. So the cost was the stale start-state metric (range 1.4e-8 to 2 there against
1.4e-8 to 0.014 at iteration 40: the hyperparameter and field scales change by orders of
magnitude in the transient), not the hyperparameter funnel. Fix: `nuts.refresh_inv_mass_at`
(commit c707a34; warm-up only: the Hessian is recomputed at the stated transition counts and the
step size re-found, same target; tests on the exact non-Gaussian target of the fixed-metric
test). The first production attempt with the stale metric (07:10 stop, 40 iterations in 40 min,
kept as `chains_gmrf16_stalemetric_attempt1/`) was restarted with the refreshed metric.

**Production settings (2026-10-28, 07:10):** 1200 Gibbs iterations, warm-up 400 with the metric
refreshed at transitions 50 and 200, burn 500 (700 post-burn samples per run; the E020 plan of
2000/500/600 would take 13 h per nside-16 run at 136 gradients per transition and 83 ms per
gradient; 1200 gives about 8 h), `adapt_mass: true`. Two runs concurrently (lane A: gmrf16,
indep16, gmrf8, indep8, gmrf4; lane B: iid16, scalar, iid8, iid4, indep4) on the 10-core Mac
mini alongside the SBC replication; the per-gradient slowdown of the two lanes is measured at
chain-1 iteration 70 and the lanes are serialized if it exceeds 1.5x. Stats review of the refresh
(2026-10-28, target invariance CONFIRMED: the kernel is fixed after warm-up for any sequence of
positive metrics; `_DualAveraging.restart` matches Stan's window restart) with its items taken:
a refresh must sit at least 50 transitions before the end of warm-up (a refresh at n_warmup - 1
left the step size at the heuristic value, acceptance 0 in the reviewer's check), each count
fires once at the first Gibbs iteration that reaches it (numstep need not divide it), and with
`adapt_mass: true` the window ending at 250 would have replaced the refreshed Hessian by a
50-sample variance estimate, so the production runs `adapt_mass: false` (the metric is the
Hessian at transition 200 for the rest of the run, the setting the 136-gradient diagnostic
measured). The first two-lane attempt (07:10 to 07:16, `adapt_mass: true`, stopped inside the Hessian
start, kept as `chains_{gmrf16,iid16}_windows_attempt2/`) was restarted with that setting at 07:17. The MC inflation of the
RMSE (above) decides whether 700 samples suffice for each comparison; a run is extended if not.

## Negative control
The single index converges to the Fisher-weighted sky mean, not the unweighted mean (review item
1): sync -2.9878 (unweighted -3.0000), dust 1.5973 (unweighted 1.5397). Its index RMSE must
therefore equal the RMSE of the truth about the weighted mean, 0.0627 for sync and 0.066 for
dust (the analyzer prints this as `rmse_of_weighted_sky_mean`), up to its posterior sd (2e-5 for
dust), and its CMB residual must be the largest. Note that the dust prior mean mu0 = 1.60
coincides with the weighted mean 1.597 (its offset from the unweighted mean is 1.9 truth sd); the
scalar dust run is therefore not a test of shrinkage away from mu0, only the sync one is.
`indep16` must be the noisiest index map where the GMRF shrinks. New here: the amplitude
hyperparameters must recover the realized mean and sd of a/g (above) in every option (they are
the same in all ten runs; a run where they do not is a sampler or data problem, not a prior
difference).

## Conclusion
(pending)

## Follow-ups
- Q16: the same comparison on the E020 PySM sky at noise scale 1 and 0.03 with robust template
  priors and sampled T_d.
- Second noise realization for the decisive pair (`--seed 2028`).
