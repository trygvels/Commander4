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
  sampled hyperparameters must recover (a check in itself: mu_c within its posterior sd of 1 and
  s_c of 0.3, up to the small bias from the floor). Unlike E020 the templates are built from the
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

Metrics: E020's, unchanged (`nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent`).
Additional check: the recovered amplitude hyperparameters against their true values (mu 1, s 0.3).

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
(timing run started 05:31 on 2026-10-28)

## Negative control
As in E020: the single index must have index RMSE equal to the truth sd (0.0615 / 0.0323) up to
its posterior sd and the largest CMB residual; `indep16` must be the noisiest index map where the
GMRF shrinks. New here: the amplitude hyperparameters must recover mu = 1 and s = 0.3 in every
option (they are the same in all ten runs; a run where they do not is a sampler or data problem,
not a prior difference).

## Conclusion
(pending)

## Follow-ups
- Q16: the same comparison on the E020 PySM sky at noise scale 1 and 0.03 with robust template
  priors and sampled T_d.
- Second noise realization for the decisive pair (`--seed 2028`).
