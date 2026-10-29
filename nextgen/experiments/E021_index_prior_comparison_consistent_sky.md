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
--template-scatter 0.3 --template-story model`, seed 2027, same bands, nside and noise as E020;
attempt 2, the sky the comparison runs on; the attempt-1 chains in
`/Users/trygvels/data/e021_consistent/` are on a misspecified sky and are excluded from the
comparison):
- beta_s(p) and beta_d(p) are the plain means of the PySM s1 and d1 index maps (degraded to
  nside 64) over each nside-16 region, broadcast back to nside 64: 3072 distinct values each.
  Truth sd 0.0615 (sync) and 0.0323 (dust) at nside 64; ranges -3.267 to -2.761 and 1.430 to 1.727.
  So the nside_beta 16 options are exactly specified, nside_beta 8 and 4 average over 4 and 16
  true regions (their RMSE floors are the truth's own region means, reported by the analyzer),
  and the single index averages over all 3072.
- Dust: one modified blackbody at T = 20.929 K everywhere (the mean of the d1 temperature map,
  the value the model fixes), amplitude a_dust(p) = g_dust(p) r_p with g_dust the PySM d1
  emission at 545 GHz in uK_RJ and r_p the per-pixel scatter of the template bullet below.
  Hence beta_eff = beta_eff3 = beta_d; the fixed-T misspecification of E020 is gone by
  construction.
- Synchrotron: a_sync(p) (nu/30)^beta_s(p) with a_sync the PySM s1 emission at 30 GHz.
- CMB: the PySM c1 map (rms 84.1 uK_CMB) times the thermodynamic-to-RJ factor, as before.
- Noise: white, rms per band 8.4, 7.2, 7.2, 6.6, 7.2, 7.2, 8.4, 4.8, 2.88, 1.92, 2.88, 9.6 uK_RJ per
  nside-64 pixel (0.03 x the E012 nside-16 levels scaled to nside 64), one realization.
- Templates and amplitudes (attempt 2, `--template-story model`, data directory
  `/Users/trygvels/data/e021_consistent_v2`): the templates g_c are the PySM amplitude maps at
  the pivots (floored for positivity at 1e-3 of their mean; no pixel hits it), and the true
  amplitudes are a_c(p) = g_c(p) r_p with r_p = max(1 + 0.3 eps_p, 0.2), eps iid N(0, 1) per
  pixel and component (the 0.2 floor is hit by 0.36 and 0.43 percent of pixels; g/true at
  |b| > 20, 16/50/84 percentiles: 0.767/0.998/1.413 sync, 0.770/1.000/1.434 dust). This is the
  hier prior's own generative story, a | g ~ N(mu g, (s g)^2) with mu = 1 and s = 0.3, which the
  sampled hyperparameters must recover. The targets are the realized values: mean(a/g) 1.00309
  (sync) / 0.99986 (dust) and sd(a/g) 0.29888 / 0.29897 for this seed (the full-information
  posterior sd is about 0.0014 for mu and 0.0010 for s). Acceptance: posterior mean of mu_c
  within 2 posterior sd of the realized mean(a/g), posterior mean of s_c within 2 sd of the
  realized sd(a/g), the same in all ten runs. Verified on the data before any run (numpy
  evaluation of the exact 12 x 12 collapsed marginal at the true index field,
  `nextgen/scripts/e021_truth_scan.py <data dir>`): the marginal's peak in mu_s is 1.0020 (realized 1.0031, sd 0.0014),
  in mu_d 1.0004 (0.9999), in log s within 0.0004 of the realized values, and a common shift
  of the index of the 768 lowest-S/N nside-16 regions peaks at +0.003 (dust, sd 0.005) and
  +0.010 (sync, sd 0.007): the score identity holds on this sky. Unlike E020 the templates are
  the smooth truth: this is deliberate (L3); the data-built templates return in Q16.
- Attempt 1 (2026-10-28 07:17 to 17:11, data directory `/Users/trygvels/data/e021_consistent`,
  kept) built the templates the other way round: g_c = a_c / r_p with a_c the PySM amplitude.
  Then a_c = g_c r_p still holds pixel by pixel, but the conditional distribution of a given g
  is not the prior's (at fixed a, a large r means a small g, so r is anti-correlated with g),
  and the exact collapsed marginal is biased on that sky: see Result, "Attempt 1 outcome".

Per-region Fisher sd (profiled over the three per-pixel amplitudes; the writer's
`sim_report.txt`, attempt-2 sky): nside_beta 4 median 0.037 sync / 0.022 dust (15 / 37 percent
of regions above the truth sd); nside_beta 8 0.079 / 0.049 (65 / 62 percent above); nside_beta
16 0.163 / 0.106 (88 / 77 percent above); full sky 0.0005 / 0.00002 (attempt-1 sky: 0.041 /
0.024, 0.087 / 0.054, 0.180 / 0.116). The same regime as E020: nside_beta 4 data-dominated, 8
marginal, 16 prior-dominated for most regions, so the comparison is a comparison of shrinkage.

Verification of the target on this data before any run (offline evaluation of
`RegionIndexLikelihood.from_arrays` on the band files, scratchpad `probe_e021_v2_like.py`;
attempt-1 numbers in parentheses): the log posterior at the truth state is 79,652 (53,859)
above the flat prior-mean start and 649,017 (519,360) above a flat beta_s = -4.0; 1-D scans of
one region's beta_s around the truth are symmetric with local sd 0.013 to 0.17, between the
amplitudes-known and profiled Fisher sd; the analytic gradient matches central finite
differences in every block (ratios 0.9995 to 1.0001). The gradient norm at the truth is 1.0e4
(2.8e4) against 1.3e6 at the flat start: the truth is at the mode to within the noise.

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
That review's "from the model class" verdict was for the attempt-1 sky and was wrong on the
templates: it checked a = g r pixel by pixel, which is not sufficient; the conditional a | g
must be the prior's, which requires g, not a, to be the fixed quantity (lesson recorded).

Second stats review (2026-10-28, after attempt 1; own marginal code, Cholesky form): the score
identity holds for the model story and fails for the data story, CONFIRMED with the mechanism
made precise: with a fixed and g = a / r, p(a | g) is proportional to p(a) p_r(a/g) over the
pixel ensemble, and since the PySM amplitude distribution falls steeply this is
errors-in-variables attenuation with g as the noisy variable, E[a | g] below g at large g and
above it at small g while mean(a/g) stays 1 (faint dust quartile, E[a|g]/g by g bins from the
lowest to the highest decile: old sky 1.29, 1.16, 1.06, 0.97, 0.84, 0.54 and a slope of a on g
of 0.71; new sky 1.00 in every bin, slope 0.98). The scan numbers reproduced on both skies
(dust +0.1143 and +0.0028, sync +0.0385 and +0.0097, hyperparameters as stated above); an
independent replay with fresh r and noise gave +0.125 and +0.127 under the data story and
-0.001 and -0.009 under the model story (sd 0.0055 / 0.0050): the shift is the exact target's,
systematic across seeds. New sky checked for other departures: the CMB per-pixel prior (84.1 or
1e4 instead of 100 uK moves the dust peak by at most 0.002; the per-pixel likelihood ignores
the CMB's spatial correlation by construction), the ratio floor (0.36 / 0.43 percent of pixels;
six model-story realizations give a dust shift of -0.002 +- 0.002, bound 0.0008), the
positivity floor (inactive), the SEDs (sky recomputed to 0 in all bands, T_fixed 2.4e-5 K),
the noise (z sd per band 0.993 to 1.006): nothing above 0.2 sd of the common shift. The
hyperparameter test is valid with the realized targets; its 2-sd band is lenient (the posterior
sd of mu, 0.0014, equals the sampling sd of the realized mean over 49152 pixels, so the
posterior mean sits within about 0.0004 of it) and a pass is not a tight check. Noted: the
`data` option at the new code does not regenerate the attempt-1 sky at the same seed (RNG
order), so the attempt-1 directory is the only copy.

## Reference
The truth maps in truth.npz (region-constant indices, the c1 CMB, the amplitudes). No exact
posterior exists at this size; sampler correctness is from E018 (L1 to L3) and the direct
verification above. The E020 chain outputs on the misspecified sky are kept as the contrast.

## Command
```bash
# attempt 1 (data-story templates, kept as evidence)
python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e021_consistent --nside 64 --noise-scale 0.03 --seed 2027 --sky consistent --nside-truth 16 --template-scatter 0.3 --template-story data
mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent/param_timing_gmrf16.yml    # 40 iterations
mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent/param_gmrf16.yml            # 1200 iterations, complete
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent --burn 500 --tags gmrf16
# attempt 2 (model-story templates), the ten runs one at a time (scratchpad run_e021_v2_queue.sh)
python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e021_consistent_v2 --nside 64 --noise-scale 0.03 --seed 2027 --sky consistent --nside-truth 16 --template-scatter 0.3 --n-iter 1200 --n-warmup 400 --refresh-mass-at 50 200 --no-adapt-mass --template-story model
for t in gmrf16 iid16 indep16 scalar gmrf8 iid8 indep8 gmrf4 iid4 indep4; do mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent_v2/param_$t.yml; done
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent_v2 --burn 500
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

**Refresh and timing at chain-1 iteration 70 (2026-10-28, 10:27):** the refresh at transition
50 does what the offline diagnostic predicted. gmrf16 went from 400 to 640 gradients per
transition before it (step sizes 0.006 to 0.02) to 63 to 272 after it (step sizes 0.04 to 0.12,
acceptance 0.4 to 0.97 in the first blocks while dual averaging settles); iid16 from 300 to 510
to 15 to 34 (step sizes 0.1 to 0.3). Per-gradient cost from the chain file times (gaps over
600 s, the Hessian starts and refreshes, excluded; transitions with at most 7 gradients give the
per-transition overhead, transitions with at least 200 the per-gradient cost): one run alone
(timing run, SBC replication in the background) 2 s overhead and 82 to 95 ms per gradient; two
runs at once 7 s overhead and 150 to 220 ms per gradient (gmrf16 median 152, iid16 222). The
slowdown is 1.8 to 2.3x, so two lanes deliver the throughput of one and the rule above applied:
at 10:31 the jobs were serialized by pausing lane B with SIGSTOP (scratchpad
`e021_serializer.sh`, one tick per minute: the mpirun job that started earliest runs, the others
wait and resume in start order; progress is kept, no run is restarted; D015 amendment). iid16 is
paused at chain-1 iteration 80 and resumes when gmrf16 ends; lane A's next run (indep16) then
waits for it. Projection at single-lane speed: gmrf16 about 2 s + 100 x 0.085 s per transition,
two chains, 1130 iterations left plus the refresh at 200 (2 x 14 min): about 7 h; iid16 about
3 h; the nside-8 and nside-4 runs have 4 to 16 times fewer index coordinates (Hessian 3084 and
774 gradients) and should be cheaper per transition; about 35 to 40 h for the ten.

**Attempt 1 outcome, gmrf16 complete (2026-10-28, 07:17 to 16:49; 1200 iterations, 700
post-burn; two lanes until 10:31, then alone):** after the refresh at transition 200 the
sampler settles at 63 gradients per transition (step size 0.081, acceptance 0.6 to 0.9, no
divergences), 16.6 s per Gibbs iteration. The analyzer (`--burn 500`): sync region-level RMSE
0.0549 (|b| > 20: 0.0609; truth sd 0.0615; single-index floor 0.0627), region z median/MAD
-0.81/0.56, |z| < 1, < 2 for 57 and 93 percent of regions, posterior sd median 0.047, field
mu -3.037 +- 0.003, sigma 0.079, rho 26 deg, ESS min 54 / median 354 over regions, hyper ESS
25 to 36; CMB residual 2.85 uK against the 2.79 floor. Dust: RMSE 0.0764 (|b| > 20: 0.0928),
worse than a constant at the sky mean (0.0323) and than the single-index floor (0.0661), z
median/MAD +1.62/1.14, |z| < 1, < 2 for 30 and 58 percent, sigma 0.068, rho 44 deg, ESS min
3.7 / median 54, hyper ESS 3 to 8. By S/N quartile the dust index is right where the data
decide (highest quartile: error +0.004, RMSE 0.013, sd 0.011) and wrong where the prior does:
lowest quartile truth 1.552, posterior mean 1.679, error +0.128 with posterior sd 0.033 (RMSE
0.134, ESS median 19). The amplitude hyperparameters miss their targets: amp mu_s 0.965 and
log s_s -1.231 against the realized 0.999 / -1.209, mu_d 0.996 / -1.193 against 0.998 /
-1.208, with posterior sds of 0.0014 and 0.003. (The reduced chi-square in the chain files,
about 70, is the alm-represented sky model through the 110-arcmin beams against the maps and
not the per-pixel fit; it is not a diagnostic here.)

**Diagnosis (exact target, numpy 12 x 12 collapsed marginal per pixel;
`nextgen/scripts/e021_faint_scan.py` (scan at the chain state), `e021_faint_cause.py`
(resimulations), `e021_truth_scan.py` (scan at the truth); chain traces by
`e021_field_trace.py <tag> <burn>`):** the chain is at its
target. A common shift of the 768 lowest-S/N dust regions evaluated at the chain's posterior-mean
state peaks at +0.025, one grid step from zero; evaluated at the truth state with the realized
hyperparameters it peaks at +0.114 with sd 0.0055 (21 sd), the log-determinant part of the
marginal gaining 456 nats against 235 lost by the chi-square part at +0.125; the sync shift
peaks at +0.039 (sd 0.007); the marginal in the amplitude hyperparameters at the true index
field peaks at mu_s 0.981 and mu_d 0.973 (sd 0.0014) and log s_d -1.201. The CMB prior width is
irrelevant (84.1, 100 or 1e4 uK: peak +0.113 to +0.117). Resimulating the same amplitudes with
fresh noise reproduces the pull (+0.13, +0.12); resimulating with a = g r drawn from the
attempt-1 templates removes it (peak 0.000, hyperparameter peaks 1.000). Cause: the attempt-1
construction g = a / r makes r anti-correlated with g at fixed a, so a | g is not N(g, (0.3 g)^2)
and the score identity fails; the determinant term, which weights pixels by g^2, then prefers a
steeper dust SED at the observed bands in the regions where the prior variance (s g M)^2 is
comparable to the noise. Not a sampler fault (the E020 checks of gradient and kernel stand) and
not a prior-volume effect of the kind in F22: a mismatch between the simulated template errors
and the model's.

**Attempt 2 (started 2026-10-28 17:15):** sky regenerated with `--template-story model` (writer
change, same seed; the noise realization differs because the RNG order changed), verified as in
"Target and priors"; the ten runs one at a time in the order gmrf16, iid16, indep16, scalar,
gmrf8, iid8, indep8, gmrf4, iid4, indep4 with the attempt-1 sampler settings (D015, D017). The
attempt-1 gmrf16 chain is kept as the evidence for F23 and as a robustness data point for Q16.

**Attempt 2, gmrf16 complete (2026-10-28 17:15 to 2026-10-29 03:59, 10.7 h alone on the machine
with the SBC replication until 03:02; 1200 iterations, 700 post-burn):** Hessian start 15 min
per chain; after the refresh at 200 the sampler runs at 63 to 127 gradients per transition (step
size 0.067, acceptance 0.7 to 0.98, no divergences), 125 gradients and 25 s per Gibbs iteration
(75 to 83 ms per gradient alone). Analyzer (`--burn 500`): sync region-level RMSE 0.0333
(|b| > 20: 0.0342; truth sd 0.0615; single-index floor 0.0630), region z median/MAD -0.08/0.55,
|z| < 1, < 2 for 78 and 98 percent of regions, posterior sd median 0.042, by S/N quartile 0.037,
0.036, 0.033, 0.027 (truth sd 0.051 to 0.081), ESS min 77 / median 352 over regions; dust
RMSE 0.0180 (0.0197; truth sd 0.0323; floor 0.0686), z +0.06/0.54, 78 and 98 percent,
posterior sd median 0.025, by quartile 0.019, 0.019, 0.020, 0.012, ESS min 8 / median lower
than sync (the dust field mean mu_d has ESS 8.5 with lag-1 autocorrelation 0.94, so its MC
error is 0.0005 against a posterior sd of 0.0016). Fields: sync mu -3.0017 +- 0.0027 (truth
mean -3.0000), sigma 0.071 +- 0.005, rho 26 +- 3 deg; dust mu 1.5410 +- 0.0016 (truth 1.5397),
sigma 0.038 +- 0.002, rho 20 +- 2 deg. Amplitude hyperparameters pass the 2-sd test: mu_s
1.00132 +- 0.00155 (target 1.00309, -1.1 sd), s_s 0.29869 +- 0.00107 (0.29888, -0.2 sd), mu_d
1.00248 +- 0.00180 (0.99986, +1.5 sd; ESS 23), s_d 0.29970 +- 0.00120 (0.29897, +0.6 sd). GLS
CMB residual 2.89 uK against the 2.79 floor (per-draw 2.99), power fractions 1e-4 and 9e-4,
regression coefficient 0.997. MC inflation of the RMSE: sync 0.0024, dust 0.0037. The faint
quartile is right where attempt 1 was wrong: dust truth 1.5516, posterior mean 1.5432 (sd
0.025); the exact-marginal common-shift scan at the chain state peaks at +0.025 and at the
truth state at 0.000. The z MAD of 0.55 against 0.67 for a standard normal says the posterior
sds are about 1.2 times too wide at the region level; the coverage figures (78 / 98 percent
against 68 / 95) say the same. Compared with the single index the GMRF at nside 16 halves the
sync error (0.033 against 0.063) and quarters the dust error (0.018 against 0.069) at this
noise.

**Attempt 2, iid16 complete (2026-10-29 03:59 to 06:50, 2.85 h; 1200 iterations, 700
post-burn):** Hessian start 11 min per chain; after the refresh at 200 the sampler runs at 15 to
31 gradients per transition (step size 0.165, acceptance 0.4 to 0.99, no divergences), 29
gradients and 4 s per Gibbs iteration (87 ms per gradient), 4.3 times fewer gradients per
iteration than gmrf16. Analyzer (`--burn 500`): sync region-level RMSE 0.0505 (|b| > 20: 0.0519;
floor 0.0630), region z median/MAD -0.10/0.56, |z| < 1, < 2 for 78 and 98 percent, posterior sd
median 0.066, ESS per region min 247 / median 929; dust RMSE 0.0275 (0.0302; floor 0.0686), z
+0.40/0.57, 69 and 96 percent, sd median 0.032, ESS min 290 / median 916. Fields: sync mu
-3.0054 +- 0.0024 (truth mean -3.0000), sigma 0.079 +- 0.002; dust mu 1.5541 +- 0.0010 (truth
1.5397), sigma 0.038 +- 0.001; mu_d has ESS 9.6 (lag-1 0.98) but its trace sits between 1.552
and 1.556 from iteration 400 on. The amplitude hyperparameters fail the 2-sd test: mu_s 0.99831
+- 0.00151 (target 1.00309, -3.2 sd), s_s 0.29833 +- 0.00112 (-0.5 sd), mu_d 1.01133 +- 0.00165
(0.99986, +6.9 sd; ESS 147), s_d 0.30236 +- 0.00115 (+2.9 sd). GLS CMB residual 2.93 uK
(per-draw 3.08) against 2.89 (2.99) for gmrf16.

**nside-16 pair (iid16 minus gmrf16, paired, block bootstrap over 192 nside-4 blocks):** sync
RMSE +0.0172 +- 0.0013, dust RMSE +0.0095 +- 0.0010, GLS CMB residual +0.054 +- 0.007 uK.
MC-inflation rule: the sync inflation is 0.0024 in both runs, below a quarter of the difference
(0.0043), so the sync comparison stands. For dust iid16 (0.0011) passes but gmrf16 (0.0037)
exceeds the quarter (0.0024), so the rule is formally not met on the gmrf16 side; that inflation
adds to gmrf16's own RMSE, so removing it can only widen the gap: the sign of the dust
comparison stands and its size is uncertain at the 0.004 level until a longer or second-seed
run. The GMRF beats iid shrinkage on both indices and on the CMB. Cost side (principle 2):
iid16 gives a median 920 ESS per region from 20,300 post-burn gradients, gmrf16 352 from
87,500 (sync), so the GMRF spends about 11 times more gradients per effective sample at nside
16 for the 0.017 / 0.010 RMSE gain and for the correct mid-S/N indices below.

**Where iid loses, and why (checked, `nextgen/scripts/e021_iid_check.py`):** by S/N quartile
(profiled Fisher sd at nside 16, high to low) the iid16 dust error means are -0.000, +0.022,
+0.032, +0.003 (RMSE 0.014, 0.033, 0.038, 0.017) against gmrf16 +0.000, +0.006, +0.007, -0.008
(0.012, 0.021, 0.020, 0.019); the truth quartile means are 1.567, 1.522, 1.519, 1.552. The iid
prior shrinks every region toward one global mu, and the posterior of mu is the
precision-weighted mean of the regions with weights 1 / (sigma^2 + s_r^2), which the bright
Galactic regions dominate: the linear-Gaussian prediction from the truth and the Fisher sds
gives mu 1.5561 (chain 1.5541; unweighted truth mean 1.5397), quartile errors -0.000, +0.026,
+0.035, +0.005 and RMSE 0.0268 (chain 0.0275). The two middle quartiles have true indices 0.035
below the weighted mean and are pulled toward it; the faint quartile is spared only because its
true mean (1.552) happens to sit near mu. The GMRF shrinks toward the neighbours instead and
follows the large-scale structure of the PySM index (rho 20 deg). The exact-marginal block-shift
scan at the iid16 chain state (collapsed likelihood plus the iid prior at the chain
hyperparameters) peaks at 0.000 for both middle dust quartiles and for the faint sync quartile,
with the likelihood alone wanting -0.01: the chain sits at the mode of its target, so this is a
property of the iid target, not a sampler failure. For sync the prediction is rougher (predicted
RMSE 0.041, chain 0.051; chain quartile errors -0.004 to -0.007 against predicted +0.001 to
+0.008): the sync offset of -0.005 (2.2 posterior sd of mu_s) is not explained by the shrinkage
alone and is left as an observation (coupling to the dust bias through the joint amplitude
solve is the candidate, not checked). The amplitude hyperparameters are a consequence: with the
dust index 0.02 to 0.03 too high over half the sky, the amplitude ratio a/g the data prefer
shifts, and mu_d lands 6.9 sd above its realized target: the exact-marginal scan of mu_d at the
iid16 chain betas (all other parameters at their chain means) peaks at 1.010 (the grid edge;
chain 1.0113) and at the truth betas at 1.000 (realized 0.99986), so the shift is implied by
the biased index field. The pre-registered reading in
"Negative control" (hyperparameters off target = sampler or data problem, not a prior
difference) was therefore wrong; corrected there.

**Attempt 2, indep16 complete (2026-10-29 06:50 to 07:33, 44 min; 1200 iterations, 700
post-burn):** independent N(mu0, 1.0^2) per region (near flat), no index hyperparameters, 31
gradients and 1.5 s per iteration. Analyzer: sync region-level RMSE 0.1424 (|b| > 20: 0.1604),
z median/MAD -0.06/0.69, |z| < 1, < 2 for 68 and 95 percent, posterior sd median 0.131, ESS min
341; dust RMSE 0.1009 (0.1216), z +0.08/0.65, 69 and 95 percent, sd median 0.069, ESS min 395.
By S/N quartile (low to high) sync 0.189, 0.160, 0.124, 0.067 and dust 0.172, 0.088, 0.056,
0.020 against gmrf16 0.037, 0.036, 0.033, 0.027 and 0.019, 0.019, 0.020, 0.012: the GMRF gains
in every quartile, by a factor 2.5 (sync) and 1.6 (dust) even in the brightest one, and 4.3 and
5.6 over the whole sky. Paired against gmrf16: sync +0.1091 +- 0.0037, dust +0.0829 +- 0.0053,
GLS CMB residual +0.275 +- 0.012 uK (3.16 uK). Amplitude hyperparameters pass: mu_s 1.00003 +-
0.00164 (-1.9 sd), s_s -0.8 sd, mu_d 1.00255 +- 0.00187 (+1.4 sd), s_d +0.5 sd. With no field
prior the region posteriors are calibrated to the precision of the test (z MAD 0.65 to 0.69
against 0.67, coverage 68 / 95 and 69 / 95 against 68 / 95), so the collapsed likelihood and the
sampler are calibrated region by region and the 1.2-times-too-wide sds of gmrf16 belong to the
field-prior part (not traced further here). The posterior sd medians 0.131 / 0.069 are below
the profiled Fisher sds 0.163 / 0.106 (the Fisher linearization is conservative; the z test
uses the chain sd). The precision-weighted mean of the indep16 region means with the iid16
sigma is 1.5519 for dust and -3.0007 for sync (iid16 mu 1.5541 and -3.0054; truth 1.5397 and
-3.0000): the dust number reproduces the iid pull toward the plane from the data alone, the
sync number does not reproduce the -0.005 iid offset (as noted above). The faint dust quartile
has error +0.033 with RMSE 0.175 under the near-flat prior. The pre-registered prediction
(indep16 the noisiest index map where the GMRF shrinks) holds; the trace scripts that read the
field hyperparameters fail on this run by construction (no `index_field/mu` in the state).

**Attempt 2, scalar complete (2026-10-29 07:33 to 07:52, 19 min; 1200 iterations, 700
post-burn):** one index per component, 4 gradients and 0.9 s per iteration. Chain: beta_s
-2.97371 +- 0.00094, beta_d 1.60077 +- 0.00012; amplitude hyperparameters mu_s 0.99295 (-7.2 sd
from the realized target), s_s +3.1 sd, mu_d 1.05020 (+32 sd), s_d +13.6 sd. Index RMSE against
the truth 0.0669 (sync) and 0.0691 (dust); GLS CMB residual 9.03 uK (per-draw 9.03, power
fractions 9e-3 and 1.4e-2, regression coefficient 0.991); paired against gmrf16: sync +0.0336 +-
0.0035, dust +0.0511 +- 0.0014, CMB +6.27 +- 2.06 uK. Negative-control check
(`nextgen/scripts/e021_scalar_check.py`): the exact collapsed marginal over all 49,152 pixels,
scanned in each scalar index at the chain's amplitude hyperparameters, peaks at -2.97374
(curvature sd 0.00084) and 1.60076 (0.00011), on top of the chain (-2.97371 +- 0.00094, 1.60077
+- 0.00012), and at the same place when the realized amplitude targets replace the chain
hyperparameters. The chain sits at the mode of its target. The pre-registered linearized
prediction (the Fisher-weighted sky mean -2.9878 and 1.5973, with RMSE 0.0627 and 0.0661) is
off by 0.014 and 0.0035, so the "RMSE must equal the weighted-mean RMSE" clause is met only
approximately (0.0669 against 0.0627; 0.0691 against 0.0661): the linearization around the
truth does not capture where the nonlinear 12-band marginal with coupled components puts a
single index. The CMB residual is the largest of all runs as required, and the hyperparameters
are far off because one index per component misfits the sky and the amplitude ratios absorb
it (the corrected reading). Note that the data-only dust optimum 1.6008 coincides with the
hyperprior mean mu0 = 1.60, as anticipated in "Negative control".

**Attempt 2, gmrf8 complete (2026-10-29 07:52 to 08:48, 56 min; 1200 iterations, 700
post-burn):** 3088-gradient Hessians, 63 gradients and 2.5 s per iteration (12 to 16 ms per
gradient), no divergences. Analyzer: pixel-level RMSE against the nside-16 truth sync 0.0379
(|b| > 20: 0.0362) and dust 0.0218 (0.0211), of which 0.0195 / 0.0146 is the unbeatable
within-region spread of the truth at nside 8; region-level RMSE with that floor removed 0.0325
(sync) and 0.0163 (dust) against 0.0333 / 0.0180 for gmrf16. Paired against gmrf16 at the pixel
level: sync +0.0046 +- 0.0007, dust +0.0038 +- 0.0007, GLS CMB residual +1.34 +- 0.51 uK (4.22 uK
against 2.89; the reference GLS at the truth's nside-8 region means is 5.37 uK, so the
resolution misfit in the plane, not the sampler, drives the CMB). Region z median/MAD -0.13/0.63
(sync) and +0.09/0.61 (dust), coverage 72 / 96 and 73 / 96 percent, posterior sd medians 0.037 /
0.020, ESS min 133 / 5 over regions (dust hyperparameters again slow: mu_d ESS 9, ESS of the
dust amplitude mean 39). Fields: sync mu -3.0004 +- 0.0028, sigma 0.066 +- 0.005, rho 34 +- 5
deg; dust mu 1.5421 +- 0.0020, sigma 0.036 +- 0.002, rho 32 +- 4 deg (the coarser basis takes a
longer correlation length). Amplitude hyperparameters: mu_s -1.2 sd, s_s +0.2 sd, s_d +1.4 sd,
mu_d 1.00384 +- 0.00190 (+2.1 sd), a marginal exceedance of the 2-sd band (one of 4 tests in the
fifth run; its direction matches a slight dust misfit of the nside-8 model to the nside-16
truth), noted and not treated as a failure. Reading: at this noise the nside-8 field is as good
as nside 16 per region and loses only what it cannot represent (the pixel-level +0.004 and the
CMB residual); the region-model ladder is therefore about representable structure, which is
what Q3 asks.

**Sampler note:** every run prints 13 numpy RuntimeWarnings at its first transition (12 ranks:
divide by zero in the GMRF prior at sigma = 0; one rank: overflow in exp). They come from the
initial step-size heuristic (eps doubling from 1, Hoffman and Gelman Algorithm 4) visiting an
absurd state (log sigma below -745, index far outside range); the density returns -inf there by
the nuts.py contract and the heuristic backs off. No effect on the target. The mass-refresh
log lines are in `chains_<tag>/logs/`, not in `run_<tag>.out`.

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

Correction (2026-10-29, after iid16): the last reading was wrong. The amplitude hyperparameters
are estimated jointly with the index fields, so a biased index field shifts them (iid16: mu_d
+6.9 sd with the mid-S/N dust index 0.02 to 0.03 too high, see the Result). The hyperparameter
test is kept as a detector of index-prior misspecification, not as a sampler check; the sampler
check is the exact block-shift scan at the chain state.
Outcome (scalar run, 2026-10-29): the single index sits at the exact mode of its target (scan in
the Result) but 0.014 (sync) and 0.0035 (dust) away from the linearized Fisher-weighted mean, so
the RMSE clause holds only approximately (0.0669 / 0.0691 against 0.0627 / 0.0661); the CMB
residual is the largest (9.0 uK); indep16 is the noisiest index map (0.142 / 0.101).

## Conclusion
(pending)

## Follow-ups
- Q16: the same comparison on the E020 PySM sky at noise scale 1 and 0.03 with robust template
  priors and sampled T_d; and the attempt-1 sky (g = a / r) as the mildest template-error
  misspecification: a 30 percent scatter whose conditional law differs from the prior's shifts
  the prior-dominated dust index by +0.13 at 4 posterior sd under hier at nside 16.
- Second noise realization for the decisive pair (`--seed 2028`).
