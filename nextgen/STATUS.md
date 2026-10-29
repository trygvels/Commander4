# Status

Newest entry on top. Keep each entry short: learned (with numbers), running, next, needs Trygve.

## 2026-10-29 (Mac mini, 21:15): E023 stage 1: 0.03 sky read; noise-1 sky exposed the per-pixel CMB prior (F26, D021), rerun with a flat prior in progress

**Learned**
- E023, 0.03 sky: the data-built template under one global scatter drives the hier target's
  bias (0.05 to above 0.1 in both indices in every quartile but the brightest; faint-sync
  runaway reproduced as a block shift); fixed T adds at most 0.03 with an exact template; the
  per-quartile scatter removes the runaway and nothing else; per-quartile mean and scatter fix
  the mid-S/N sync (to 0.01) but leave the dust at 0.05 to 0.10, because the 353 GHz-based
  dust template's error varies with T_d, hence with the index: no scalar-hyperparameter law
  fixes it. Neither candidate goes to stage 2 as a fix.
- E023, noise-1 sky: the baseline cell (exact templates, true T) was biased by 0.05 to 0.1 in
  the fainter dust quartiles; traced (three hypotheses, one confirmed) to the per-pixel CMB
  prior N(0, 100^2) of the E018/E021 target, informative at realistic noise (F26). With a flat
  CMB prior the baseline is clean. All noise-1 cells are being rerun with the flat prior
  (`--cmb-sd 1e4`, from 21:04, about 3 h); production runs at realistic noise use it (D021).

**Running**
- E023 noise-1 rerun (flat CMB prior); E021 seed-2028 triple (gmrf16 near its end, then iid16,
  indep16).

**Next**
- Read the noise-1 cells; E023 Conclusion; Q16 and ROADMAP updates; design of a template-error
  field model (the dust template's error must be a field correlated with T_d) or of a
  component-separation-based template as the Q16 candidate for stage 2.
- Seed-2028 paired differences into E021 and F24.

**Needs Trygve**
- New: the per-pixel CMB prior width in the collapsed likelihood for real-data runs (D021 sets
  it flat; the proper correlated CMB prior belongs to the CG amplitude step, Q5). Unchanged
  items as before.

## 2026-10-29 (Mac mini, 16:10): E024 done, dense Hessian metric adopted for production (F25, D020)

**Learned**
- E024 (gmrf8 with the dense metric, same data and target as the diagonal chain): every
  statistic within 1 combined MC SE of the diagonal chain; mu_d ESS 9 to 331 of 700 (27 times
  per gradient), faint dust regions 6 to 135, dust amplitude mean 39 to 376, chain's worst
  quantity 14 times better per gradient; sync regions 1.6 to 2.2 times worse per gradient
  (pre-registered clause failed, recorded; 84 against 63 gradients per iteration). Adopted for
  production runs; default unchanged.

**Running**
- E023 stage 1 scans (E021 control sky, then PySM at 0.03 and at noise 1).
- E021 seed-2028 gmrf16, iid16, indep16 (from 11:13).

**Next**
- E023 stage-1 reading; stage-2 chain on the noise-1 sky with the dense metric (the first
  nside-16 check of the reviewed code).
- Seed-2028 paired differences into E021 and F24; Paper A Table 2 final.
- Efficiency follow-up: the sync-side loss under the dense metric (per-block step size or a late
  refresh); the two persistently non-positive Hessian directions.

**Needs Trygve**
- Unchanged.

## 2026-10-29 (Mac mini, 12:10): dense Hessian metric for the index-field NUTS merged after stats review (E024 run in progress); E023 scans and the seed-2028 triple running

**Learned**
- Diagonal-metric baseline (E021 gmrf8): dust field mean mu_d ESS 9 of 700, dust amplitude
  mean ESS 39, dust region ESS min 6, median 66; sync fine (mu_s 69, regions median 575). In
  centred coordinates the likelihood curvature is diagonal in pixels but the basis is global, so
  a diagonal metric cannot follow the couplings of the field mean and the hyperparameters.
- Implemented `nuts.init_inv_mass: fd_hessian_dense`: the full finite-difference Hessian (same
  2 dim gradients as the diagonal estimate) as a dense inverse mass, eigenvalues 1/|lam| capped
  at 10, refreshed at the configured transitions; same target for any fixed metric.
  Stats-reviewed (consistency of the Hamiltonian confirmed; a silent corruption of the Hessian
  by rejected stencil points found and fixed in both initializers), 952 tests green. On an
  8-dimensional correlated Gaussian the dense metric takes 3 leapfrogs per transition against 7
  for the best diagonal one.
- E023 code control: the numpy exact-target marginal agrees with Commander4's kernel to 4e-12.

**Running**
- E024: gmrf8 with the dense metric on the E021 data (from 11:56, same target as the diagonal
  gmrf8 chain; pre-registered reading in the record: posterior agreement within 4 MC SE, mu_d
  ESS per gradient up by at least 3).
- E023 stage 1: scans on the E021 control sky, then the PySM sky at 0.03 and at noise 1 (one
  core, slow under load; several hours).
- E021 seed-2028 gmrf16, iid16, indep16 (from 11:13; about 14.5 h unloaded, slower now).

**Next**
- E024 result into the record; if adopted, the dense metric becomes the production setting for
  the E023 stage-2 chains and is noted in Paper A Sect. 5 (sampler).
- E023 stage-1 reading per the pre-registration; then the stage-2 chain on the noise-1 sky.
- Seed-2028 paired differences into E021 and F24; Paper A Table 2 final.

**Needs Trygve**
- Unchanged: Paper A presentation of the E022 diagnostics; `COMP_CL_TYPE = none` semantics
  (Q17 item 1); the hier default (D003/D011); Q15 (zero-level model).

## 2026-10-29 (Mac mini, 11:45): E021 Conclusion corrected after stats review (resolution reading withdrawn, ordering stands); E023 designed and under review; seed-2028 triple running

**Learned**
- Third stats review of E021 (reproduced): the nside-16 and 8 orderings and the iid mechanism
  hold; nside 4 is suggestive only; the "per-region error flat with nside" reading was an
  artefact of the unweighted region reference. Against the S/N-weighted reference the GMRF's
  per-region error falls 0.0333 / 0.0318 / 0.0295 (sync) and 0.0180 / 0.0153 / 0.0138 (dust) at
  nside 16 / 8 / 4, while independent regions lose factors 3.5 / 4.9: the spatial prior keeps
  most of the averaging coarse regions buy. The negative control had quoted attempt-1 numbers
  (dust clause actually met, sync 0.012 off); the MC-inflation rule fails for more pairs than
  recorded but its effect on an RMSE is at most 0.0004; coverage at nside 8 and 4 is
  reference-dependent; block-bootstrap SEs are understated. F24, Conclusion, Q3, M2.5, Paper A
  outline rewritten; D019 withdraws the resolution rule from D018 until a field finer than the
  sky's structure is tested.
- E023 (Q16) designed: exact-target block-shift scans of the profiled hier marginal on the PySM
  sky at noise scale 1 and 0.03 (template built from data or exact; T fixed or true; global or
  per-quartile scatter); script written and smoke-tested (3 to 7 s per state). Stats review of
  the target: formula correct; required before running: an expected-log-likelihood scan on the
  noiseless sky (with the tr(C^-1 N) term) so the thresholds are decidable at noise 1, a Laplace
  term for the profiled hyperparameters, a noise-free region-fit reference instead of region
  means, quartiles from template S/N, a mu-per-quartile candidate, and a kernel-agreement check
  instead of the F21 reproduction as the code control. Being implemented.

**Running**
- E021 seed-2028 gmrf16, iid16, indep16 (from 11:13, about 14.5 h).

**Next**
- E023 script v2 with the review's changes, then the control run on the E021 sky and the two
  PySM skies (about 1.5 h each on one core, alongside the queue).
- When the seed-2028 triple ends: paired differences on the new sky, append to E021 and F24.

**Needs Trygve**
- Unchanged: Paper A presentation of the E022 diagnostics; `COMP_CL_TYPE = none` semantics
  (Q17 item 1); the hier default (D003/D011); Q15 (zero-level model).

## 2026-10-29 (Mac mini, 11:30): E021 closed (F24): the GMRF prior wins at every nside; resolution costs nothing per region under it; seed-2028 triple started

**Learned**
- gmrf4, iid4, indep4 (30, 22, 22 min): region-level RMSE 0.033 / 0.020, 0.036 / 0.025, 0.043 /
  0.026; paired against gmrf4 directly iid4 +0.0023 +- 0.0005 and +0.0034 +- 0.0005, indep4
  +0.0073 +- 0.0013 and +0.0041 +- 0.0008. Same ordering, smaller gaps; the GMRF-over-iid gap is
  inside the MC inflation of 700 draws (not established at nside 4 by the pre-registered rule,
  though the inflation sits on the GMRF side and can only understate it). All three under-cover
  (55 to 60 percent at 1 sd) and fail the amplitude-hyperparameter test, including indep4 with
  no index hyperparameters: one index per 256 pixels misfits the nside-16 truth, so the test
  detects any index-model misspecification. iid4 is still the most biased (dust z median +0.89,
  mu chain unresolved at ESS 2). CMB residual 5.56 uK in all three (resolution).
- Ladder (GMRF): region-level error flat, 0.0333 / 0.0180, 0.0325 / 0.0163, 0.0327 / 0.0200 at
  nside 16, 8, 4; pixel-level error and CMB residual (2.89, 4.22, 5.56 uK) track the floor; cost
  125, 63, 39 gradients per iteration. Resolution rule: as fine as the budget allows under a
  spatial prior (D018). E021 Conclusion written; F24; Q3 core answered; M2.5 updated; Paper A
  outline has the ladder point and a Table 2 candidate.
- Re-ranked OPEN_QUESTIONS: Q16 (robustness on misspecified skies) first, then Q15 (zero levels,
  needs Trygve), then the Q4 efficiency items, then Q8/Q7, then Q14 second half and Q5.

**Running**
- E021 seed-2028 realization of gmrf16, iid16, indep16 in `/Users/trygvels/data/e021_consistent_v2_seed2028`
  (queue started 11:13, about 14.5 h; pre-registered expectation in the E021 Follow-ups).

**Next**
- When the seed-2028 triple ends: paired differences on the new sky, append to E021 and F24.
- Meanwhile (does not need the machine): Q16 design (E020 sky at noise scale 1 and 0.03 with
  robust template priors and sampled T_d; attempt-1 sky as the template-law case), the M2.4
  efficiency item for mu_d (non-centred field mean), and the Paper A Table 2 draft from the ten
  E021 runs.

**Needs Trygve**
- Unchanged: Paper A presentation of the E022 diagnostics; `COMP_CL_TYPE = none` semantics
  (Q17 item 1); the hier default (D003/D011). New: Q15 (zero-level model) is now the highest
  ranked item that needs a model decision.

## 2026-10-29 (Mac mini, 08:15): E021 attempt 2 iid16 done; the GMRF beats iid shrinkage at nside 16 on both indices and the CMB; iid biases the mid-S/N dust index

**Learned**
- iid16 (2.85 h, 29 gradients per iteration against 125 for gmrf16): sync region-level RMSE
  0.0505 and dust 0.0275 against gmrf16 0.0333 and 0.0180; paired differences +0.0172 +- 0.0013
  and +0.0095 +- 0.0010, GLS CMB residual +0.054 +- 0.007 uK. The MC-inflation rule holds for
  sync; for dust the gmrf16 inflation (0.0037) exceeds a quarter of the gap, which can only
  understate it. Fair cost: the GMRF spends about 11 times more gradients per effective sample.
- Mechanism verified (E021 Result, `e021_iid_check.py`): iid shrinks every region toward one
  mu, and mu is the precision-weighted mean dominated by the Galactic plane (predicted 1.556,
  chain 1.554, sky mean 1.540); the two middle S/N quartiles of the dust index lie 0.035 below
  it and come out +0.022 and +0.032 too high (gmrf16: +0.006, +0.007). Exact block-shift scans
  at the chain state peak at zero: a property of the iid target, not the sampler. The amplitude
  hyperparameters fail their 2-sd test in iid16 (mu_d +6.9 sd) as a consequence; the
  pre-registered reading of that test (sampler or data problem) was wrong and is corrected in
  the record: it detects index-prior misspecification.
- The 13 numpy RuntimeWarnings at the first transition of every run come from the initial
  step-size heuristic probing an absurd state; the density returns -inf there by contract.

- indep16 (44 min, near-flat per-region prior): sync RMSE 0.142, dust 0.101, CMB residual 3.16 uK;
  region posteriors calibrated to the precision of the test (z MAD 0.65 to 0.69, coverage 68 / 95
  and 69 / 95), amplitude hyperparameters within 1.9 sd. The GMRF gains in every S/N quartile
  (factor 2.5 sync and 1.6 dust even in the brightest), 4.3 and 5.6 over the sky.

- scalar (19 min): beta_s -2.9737, beta_d 1.6008, each on the exact mode of the collapsed marginal
  (scan peak within 0.1 posterior sd), CMB residual 9.0 uK (largest, as required), amplitude
  hyperparameters 7 to 32 sd off (one index misfits the sky). The linearized Fisher-weighted
  prediction missed the sync optimum by 0.014, so that clause of the negative control holds only
  approximately; recorded.

- gmrf8 (56 min): region-level RMSE 0.033 / 0.016, as good as gmrf16 per region; at the pixel level
  +0.005 / +0.004 worse and the CMB residual 4.2 uK against 2.9, both from structure the nside-8
  field cannot represent (the truth is region-constant at nside 16). Calibration z MAD 0.63 / 0.61,
  coverage 72 / 96 and 73 / 96; mu_d amplitude hyperparameter at +2.1 sd (marginal, noted).

- iid8 and indep8 (35 and 26 min): the nside-16 ordering repeats, GMRF 0.033 / 0.016, iid 0.041 /
  0.023, independent 0.070 / 0.047 region-level RMSE; iid again pulls the dust field mean to the
  plane (1.556 against 1.540) with dust z median +0.44 and fails the amplitude-hyperparameter
  test (mu_d +4.5 sd), and mixes badly (hyperparameter ESS 3 to 15 of 700); indep8 is calibrated
  and passes. At nside 8 the CMB residual (4.2 uK in all three) is set by the resolution misfit.

**Running**
- E021 attempt 2 gmrf4, iid4, indep4 (from 09:49, about an hour). A waiter runs the full
  ten-run analyzer when the queue ends.

**Next**
- indep16 as the no-shrinkage reference (its per-region posterior sds are the s_r of the
  weighted-mean argument), scalar as the negative control, then the nside ladder; after all ten:
  E021 Conclusion, F24, Q3, M2.5; second noise seed for the decisive pair.

**Needs Trygve**
- Unchanged: Paper A presentation of the E022 diagnostics; `COMP_CL_TYPE = none` semantics
  (Q17 item 1); the hier default (D003/D011).

## 2026-10-29 (Mac mini, 03:30): hier SBC replication PASS (F16 established); E021 attempt 2 gmrf16 near its end

**Learned**
- The seeded SBC replication of the hier target (300 fresh replicates, non-centred NUTS, 28 h on
  a shared machine) passes the full pre-registered battery (Bonferroni min p 0.142, KS p 0.33,
  all 16 ECDFs inside the 99 percent band, TARP 8 of 8) and the one-sided test on the pixel that
  failed in the first run (rank mean 0.490 against the 0.528 threshold; pooled over 600: 0.520,
  +1.7 se). By the pre-registered rule the first run's band violation is a fluctuation and the
  SBC is PASS on 600 replicates: F16 loses its "provisional", M2.4a is done, L2 passed for the
  hier target at nside_beta 4 (E016 Conclusion written).

- E021 attempt 2, gmrf16 done (04:00, 10.7 h): on the corrected sky the GMRF index field at
  nside 16 recovers sync to 0.033 and dust to 0.018 region-level RMSE (single-index floors 0.063
  and 0.069; truth sds 0.062 and 0.032), z median -0.08 / +0.06 with MAD 0.55 (posterior sds
  about 1.2x too wide), 78 / 98 percent of regions within 1 / 2 sd, amplitude hyperparameters
  within 1.5 sd of their realized targets, CMB residual 2.89 uK against the 2.79 floor, no
  faint-region bias (exact-marginal scan at the chain state peaks at +0.025 of a 0.025 grid).
  Dust hyperparameter mixing is slow (mu_d ESS 8.5 of 700; MC error 0.0005 against sd 0.0016).

**Running**
- E021 attempt 2 iid16 (from 03:59, about 4 to 5 h), then indep16, scalar, nside 8 and 4
  (about 30 h); a waiter runs the analyzer on the nside-16 pair when iid16 ends.

**Next**
- The nside-16 pair comparison (GMRF against iid shrinkage) with the MC-inflation rule; then the
  rest of the queue.

## 2026-10-28 (Mac mini, 17:30): E021 attempt 1 stopped after gmrf16; the sky's template construction biased the target; attempt 2 running on a corrected sky

**Learned**
- gmrf16 on the attempt-1 sky finished (1200 iterations, 9.5 h, 63 gradients per transition
  after the refresh at 200). Sync is fine (region-level RMSE 0.055 against the 0.063
  single-index floor, z median -0.8, 93 percent of regions within 2 sd, ESS 54 to 640). Dust is
  not: RMSE 0.076 against a 0.032 constant-at-the-mean floor, the faint three quarters of the
  sky 0.13 above the truth at a stated sd of 0.033, ESS 4 to 19 there, amplitude
  hyperparameters 2 to 24 sd off their targets.
- The chain is at its target. The exact collapsed marginal (numpy, 12 x 12 per pixel) at the
  true index field prefers the same +0.114 shift of the faint dust regions (21 sd; determinant
  term +456 nats against chi-square -235) and hyperparameter means 2 percent low. Cause: the
  sky built its templates as g = a / r, so a given g does not follow the prior's N(g, (0.3 g)^2)
  law and the score identity fails. Redrawing a = g r removes the shift exactly; the CMB prior
  width is irrelevant (F23, D017; stats review CONFIRMED with the mechanism made precise:
  errors-in-variables attenuation, E[a | g]/g from 1.29 to 0.54 across template deciles on the
  old sky, 1.00 on the new; replays with fresh scatter reproduce +0.125 / -0.001).
- The corrected sky (`e021_consistent_v2`, `--template-story model`) passes: shift peaks
  +0.003 +- 0.005 (dust) and +0.010 +- 0.007 (sync), hyperparameters within 1 sd of the
  realized values, log posterior at the truth 79,652 above the flat start, gradients to 5e-4.

**Running**
- E021 attempt 2: the ten runs one at a time from 17:15 (gmrf16 first, about 8 h; iid16,
  indep16, scalar, then nside 8 and 4; about 40 h in all), scratchpad `run_e021_v2_queue.sh`.
- SBC replication (seed 2027), about 15 h left.

**Next**
- After the nside-16 pair: analyzer with `--burn 500`, MC inflation check, hyperparameter
  targets (1.0031 / 0.2989 sync, 0.9999 / 0.2990 dust), E021 Result, F24, Q3, M2.5.
- Q16 gains the attempt-1 sky as the mildest template-law misspecification.

**Needs Trygve**
- Unchanged: Paper A presentation of E022, `COMP_CL_TYPE = none` semantics (Q17), hier default
  (D003/D011).

## 2026-10-28 (Mac mini, 12:45): E022 closed with the released BeyondPlanck chains; a decision for Trygve on Paper A

**Learned**
- The released BP chain files store the per-region polarized beta_s per sample; two chains
  (1001 samples each) were read over HTTP range requests from the 1.8 TB files (613 MB each).
  From the exact two-block identity var(delta beta)/2 = E[var_cond]: r = 2.4 (Spur) and 1.75
  (Plane), lag-1 0.83 and 0.67 in both chains, matching the nside-16 surrogate with the flat
  CMB prior (2.35 / 2.11, 0.82 / 0.78). IPS tau 17 to 43 (the paper's 30 to 50) against 5 to
  11 from lag-1 (a slow component); 23 to 59 effective samples per chain. Both
  likelihood-sampled regions drift towards flatter indices in both chains (+0.065 to +0.092
  per 1000 samples, 1.5 to 2.0 sigma per chain-region, none in the prior-only regions), the
  direction of the flat-prior volume pull, but the four tests share a start and pair within
  chains, and chain 1's movement is a last-quarter joint excursion: consistent with an
  unsettled chain and with a transient, establishing neither (third review round). The
  reported sds (0.06, 0.07) are the spread of chains with 23 to 59 effective samples that
  cannot be shown to have settled; that the flat-prior target lies at flatter indices is an
  expectation from the surrogate and a toy, not a chain measurement (E022, F22; stats review
  in three rounds: targets, code, inputs and chain numbers CONFIRMED). Chains 3 and 4 share
  the same start, mix better (tau 9 to 11, 89 to 114 effective samples) and move the same way
  more weakly; across the four chains the movement averages +0.04 (Spur) and +0.07 (Plane) per
  1000 samples and stays an indication.
- Method point for Paper A Sect. 4: the coupling and the runaway are carried by the free
  per-pixel CMB amplitude (0.19 uK of signal against 1.5 to 2 uK of noise); a C_l-sized CMB
  prior removes both on the surrogate. Resolution scaling on the real data (LFI and 353,
  nside 16 to 64): with flat priors the Spur's collapsed target goes -2.88, -2.13, boundary,
  the Plane -3.16, -2.88, -2.65 (r 2.6 to 6, tau 18 to 77); the C_l-sized CMB prior
  stabilizes Spur, High Latitude and Galactic Center but not the Plane (flat dust and
  synchrotron priors): data-scaled priors on every component (E016 form) or the profile are
  the hypotheses to test. Sixth review round, Q17 item 2 measured on a released sample: the
  conditional width at BP's resolution is 0.016 (Spur) and 0.0085 (Plane), the synchrotron
  amplitude is bound by BP's gauss C_l prior, the amplitude block is 30 and 4 percent of the
  chain's per-step motion, and the flat-prior volume pull in BP's chain is at most 1e-3 of
  the surrogate's: the reported sds are consistent with posterior widths under BP's full
  model at 23 to 114 effective samples per chain; the mechanism that neutralizes the volume
  term (the CMB's treatment in the CG step, most likely) is not identified. The Plane's
  resolution dependence is the misspecification of one
  index for a region whose index spans 1.5: its nside-4 sub-regions are resolution-stable
  (weighted mean -3.25 at nside 16, 32 and 64) while the single fit moves -3.30 to -2.91,
  the real-data case for index fields at the data's resolution. Fifth review round: BP's effective amplitude element count is undetermined
  (7k to 4M per region), so the released chains cannot distinguish an unsettled flat-prior
  chain from a stationary one driven by the other Gibbs blocks; the record states the
  dichotomy and its discriminator (Commander3 prior semantics, or one conditional-width
  measurement on a released sample).

**Running**
- E021 gmrf16 alone (ends about 17:30), then iid16 and the other eight in start order.
- BP chains 3 and 4 reader; E022 final stats check; E016 hier SBC replication (about 18 h).

**Next**
- E021 analysis once the nside-16 pair is in; SBC verdict; E022 follow-up on the resolution
  scaling of the volume term (nside 32, 64) for Paper A Sect. 4; Q8; Paper A outline update.

**Needs Trygve**
- How to present the E022 BP XV chain diagnostics (23 to 114 effective samples per chain,
  reported sds consistent with posterior widths, amplitude block a minor part of the chain
  motion, flat-prior volume pull absent for an unidentified reason) in Paper A, given BP XV
  is his paper, and whether he can say from the Commander3 code what `COMP_CL_TYPE = none`
  does to the CMB amplitude step (Q17 item 1): as a
  worked example of the prior-volume effect with the numbers above, or only as the method
  recommendation (C_l-sized CMB prior or profile) with the diagnostics in an appendix. The
  record and F22 hold the numbers either way; nothing is published without his decision.
- The hier default question (D003/D011) stands.

## 2026-10-28 (Mac mini, 11:50): E022 (BP XV polarized beta_s Gibbs coupling) done on a surrogate, stats-reviewed; Q7 data point; E021 running serialized

**Learned**
- E022 (F22, draft until the re-check): on an nside-16 surrogate of the BeyondPlanck XV
  polarization data (LFI 30 to 70, WMAP Ka/Q/V nine-year nside-16 maps with their dense
  covariances fetched from LAMBDA, NPIPE 353, the four BP XV regions, prior N(-3.3, 0.1)) the
  collapsed-to-conditional sd ratio of beta_s is 2.1 to 2.35 (Spur, Plane) with the released
  configuration's flat CMB amplitude prior, and a two-block Gibbs sampler with exact
  conditional draws has autocorrelation times of 7 to 11 (Gaussian prediction 8 to 10; lag-1
  equals the exact identity to 0.005). With a Gaussian CMB prior of the size the Planck EE
  spectrum implies (0.19 uK) the coupling vanishes (ratio 1.02 to 1.04). BP XV's reported 30 to
  50 does not understate this coupling; F14's ratio 180 does not transfer to polarization at
  this S/N (synchrotron S/N 2 to 3 per pixel). The flat-prior marginal is shifted from the
  profile by +0.27 (Spur) and runs to the flattest allowed index in the High Latitude region;
  the C_l-sized CMB prior removes the runaway: the prior-volume effect on real data, carried by
  the free per-pixel CMB amplitude. Stats review CONFIRMED WITH ITEMS; all must-fix and
  should-fix items applied (CMB-prior and noise variants, MC errors, lag-1 identity);
  re-check requested.
- Q7 data point: BP's simulated LFI noise rms at nside 16 and 60 arcmin is 1.6 times white-noise
  propagation in all three bands (ratio 1.4 to 1.9 across the sky); a white per-pixel
  likelihood there overstates index precision by 1.6.
- The released BP chain files store the per-region index per sample compactly; a reader over
  HTTP range requests (no download of the 1.8 TB file) is fetching BP XV's actual beta_s chain
  for its autocorrelation.

**Running**
- E021 gmrf16 alone (chain-1 iteration 160 at 10:56, about 20 s per iteration, ends about
  17:30 with the refresh at 200), then iid16 (paused at 80), then the other eight in start
  order (serializer).
- BP chain-1 beta_s reader; E022 re-check by the stats-reviewer; E016 hier SBC replication
  (about 19 h left).

**Next**
- E022: BP chain autocorrelation against 30 to 50 and against 7 to 11; close the record and
  F22 after the re-check. Follow-up on the resolution scaling of the volume term (nside 32, 64)
  for Paper A Sect. 4.
- E021 analysis once the nside-16 pair is in; SBC verdict; E014 BP XV cell is E022 now
  (E014 Follow-ups to be pointed there); M2.6 (Q7 partly measured, Q8); Paper A outline update.

**Needs Trygve**
- Nothing new. The hier default question (D003/D011) stands.

## 2026-10-28 (Mac mini, 10:40): E021 production on the refreshed metric; two lanes give no throughput, jobs serialized

**Learned**
- The warm-up metric refresh (`nuts.refresh_inv_mass_at`, c707a34 plus review items 5b07e7d;
  stats-reviewed, target invariant) does what the offline diagnostic predicted. At transition 50
  gmrf16 dropped from 400 to 640 gradients per transition (step 0.006 to 0.02) to 63 to 272
  (step 0.04 to 0.12); iid16 from 300 to 510 to 15 to 34 (step 0.1 to 0.3). Without it the ten
  runs would have cost about 37 h each.
- Two concurrent runs on the 10-core Mac mini cost 150 to 220 ms per gradient against 82 to 95 ms
  for one run (chain file times, Hessian gaps excluded; other work per transition 2 s alone,
  7 s with two). Two lanes deliver the throughput of one, so the pre-registered rule applied:
  from 10:31 one job runs at a time, the others are paused with SIGSTOP and resume in start
  order (scratchpad `e021_serializer.sh`, D015 amendment); nothing is restarted.
- E021 production settings: 1200 iterations, warm-up 400, refresh at 50 and 200,
  `adapt_mass: false` (a window ending at 250 would have overwritten the refreshed Hessian),
  burn 500. Stats review of the E021 design CONFIRMED the model class; the must-fix negative
  control against the Fisher-weighted sky mean is in the record (sync 0.0627, dust 0.066).

**Running**
- E021 gmrf16 (chain-1 iteration 80 at 10:30; about 7 h at single-lane speed including the
  refresh at 200), then iid16 (paused at iteration 80, about 3 h), then indep16, scalar, gmrf8,
  iid8, indep8, iid4, gmrf4, indep4 in start order; about 35 to 40 h for the ten.
- E016 hier SBC replication (seed 2027; about 20 h left).

**Next**
- Per run: trace, chi-square, hyperparameters against the realized values. Analyzer with burn 500
  once the nside-16 pair is in; MC inflation check; second seed for the decisive pair; Result,
  Conclusion, F22, Q3, ROADMAP M2.5.
- SBC verdict by the pre-registered rule; F16 final or withdrawn.
- E014 BP XV cell; M2.6 (Q7, Q8); Paper A outline update with E017 to E021 numbers.

**Needs Trygve**
- Nothing new. The hier default question (D003/D011) stands.

## 2026-10-28 (Mac mini, 05:45): E020 closed as a design flaw (F21), sampler exonerated by direct computation; E021 on a model-consistent sky designed, under stats review, timing run started

**Learned**
- The 1023-gradient saturation of every E020 attempt was not a mass-matrix problem. Attempt 3
  (fixed Hessian metric, `adapt_mass: false`) reproduced attempt 2 transition for transition until
  the first window at iteration 100, and saturated from iteration 30 like the others. The saved
  chain states, evaluated offline with the same likelihood built from the band files, have a log
  posterior 104,122 above the state at the PySM truth by iteration 182 (log-likelihood part
  +123,001), with beta_s mean -4.05, sd 3.84, range -26 to -0.6 and beta_d shifted by +0.08 at the
  Galactic centre (30 Fisher sd). The sampler was climbing the target. The target is right: the
  analytic gradient matches finite differences to five digits in every block at the truth and at
  the runaway state, and the C++ kernel equals a direct 12 x 12 Gaussian marginal per pixel in
  every region checked, including regions at beta_s = -9.7 (constant offset only).
- Two misspecifications of the E020 design drive it (F21): the fixed-T region-constant dust model
  misfits the PySM sky in the plane (the fit gains 101,841 in chi-square by moving both indices),
  and in the faint half of the sky the data-built template is wrong by factors 0.41 to 1.66 under
  one global scatter s, so killing the synchrotron in all bands but 23 GHz lowers the
  amplitude-prior cost more than it costs in fit. At 0.03 x E012 noise both dominate; every prior
  option would have been scored on the same runaway. E020 production stopped at 05:10 (D014); the
  data and chains are kept as evidence. The two NUTS changes (relative floor, `adapt_mass`) are
  correct and kept but were not the fix; the attempt-2 diagnosis is withdrawn in the record.
- E021: the same ten options, metrics and analyzer on a sky drawn from the model class
  (`pysm_intensity_sim.py --sky consistent`: region-constant indices at nside 16 from the PySM
  region means, one MBB at the fixed T, PySM amplitudes, templates equal to truth over
  1 + 0.3 eps). Verified offline before any run: the truth state is 53,859 above the flat start and
  519,360 above a flat beta_s = -4; 1-D scans are symmetric about the truth; gradient against
  finite differences agrees. Q16 opened for the robustness question (non-stationary or
  heavy-tailed template scatter, sampled T_d) on the kept PySM sky at realistic noise.

**Running**
- E021 timing run, gmrf16, 40 iterations, `adapt_mass: true` (started 05:31; decides the metric
  setting and the iteration count for the ten production runs, target about two days in total).
- stats-reviewer on the E021 design and sky (sky from the model class, template construction,
  fairness of the metrics across nside_beta, analyzer edge cases).
- E016 hier SBC replication (seed 2027, 300 replicates; about 24 h left).

**Next**
- Timing result, then the E021 production loop (ten runs); analysis, Result, Conclusion, F22,
  Q3, ROADMAP M2.5. Apply the reviewer's must-fix items before launch.
- SBC verdict by the pre-registered rule; F16 final or withdrawn.
- E014 BP XV cell; BP XVI; M2.6 (Q7, Q8); Paper A outline update with E017 to E021 numbers.

**Needs Trygve**
- Nothing new. The hier default question (D003/D011) stands.

## 2026-10-27 (Mac mini, 22:55): E020 relaunched on a fixed Hessian metric; hier SBC borderline FAIL, replication running; BeyondPlanck papers read; monopole design decided

**Learned**
- E020 gmrf16 needed two fixes to run at all. Attempt 1: Stan's absolute mass-floor (1e-4)
  inflated the mass of the 1e-8-variance bright-region coordinates after the Hessian start and
  pinned the step size (fix 0e7c282: floor relative to the starting metric). Attempt 2: the warm-up
  windows themselves replaced the Hessian metric by variance estimates from 25 to 200 correlated
  draws over 6152 coordinates (precision range 5e6): 7 gradients per transition before the first
  window, 1023 after, 1.6 min per iteration, 55 h per run. Fix: `nuts.adapt_mass: false` (fadec92,
  b93592e; the Hessian metric is kept, only the step size adapts; stats review CONFIRMED, any
  fixed metric leaves the target unchanged; `inv_mass` is now a copying property). Attempt 3
  launched 22:26; first transitions on the metric take 1 to 7 gradients.
- E016 hier SBC (300 replicates, all clean): Bonferroni, KS and TARP pass, but the ECDF
  simultaneous band fails on one of 16 quantities, `beta_faint1` (+2.9 se rank mean, standardized
  error +0.15 +- 0.06), at the band edge. No correlation with any sampler covariate; the paired
  faint pixel with the same template amplitude is clean. FAIL by the pre-registered rule;
  replication with seed 2027 (300 replicates) running with a pre-registered one-sided test
  (rank mean above 0.528 rejects). F16 stays provisional.
- BeyondPlanck I, III, XIII, XV read in full (research note with section references): BP XV
  sampled beta_s conditionally on amplitudes with enough inner Metropolis steps for lag-n
  correlation below 0.1 and reports a Gibbs correlation length of 30 to 50 over about 500 samples;
  BP I ran 4 x 1000 samples with no Gelman-Rubin; BP XIII disabled the intensity beta_s likelihood
  and sampled nu_p full-sky with band monopoles re-estimated per proposal. F14's scope is now
  explicit: one intensity configuration, not transferable to BP XV without the planned E014 cell.
  Paper A does not claim published results were under-mixed.
- D013: per-band zero levels will be offset columns with Gaussian priors in the collapsed
  per-pixel likelihood (M3.3), not a CG-step monopole component. Q15 item 2 closed.

**Running**
- E020 production attempt 3 (gmrf16 since 22:26, then nine more runs); waiter at chain-1
  iteration 175 gives the per-run budget.
- E016 SBC replication, seed 2027 (since 22:47, about 31 h).

**Next**
- E020 budget; shorten the schedule if a run exceeds about 12 h; analysis, second seed, Result,
  Conclusion, F21, Q3, ROADMAP M2.5.
- SBC replication verdict; E016 Conclusion; F16 final or withdrawn; M2.4a.
- E014 BP XV cell (conditional-versus-collapsed width per region) before any statement about BP
  mixing; BP XVI read; BP I Sects. 9.2 to 11 need the PDF.

**Needs Trygve**
- Nothing new. The hier default (D003/D011) stays open.

## 2026-10-27 (Mac mini, 10:55): E020 (M2.5 first cut) redesigned twice, timed, and launched; non-centred coordinates fail at nside 64 S/N

**Learned**
- E020 design passed two stats reviews (second verdict: correct, every should-fix taken). PySM
  d1+s1+c1 at nside 64, 12 bands, noise 0.03 x E012 levels, fixed T_d, prior means offset from
  the truth, 5-degree smoothed data-built templates floored at the 1st percentile; ten runs
  (GMRF, iid and independent region priors at nside_beta 4, 8, 16, plus a single index); metrics
  pre-registered incl. the GLS CMB at index draws with a true-SED floor and a paired block
  bootstrap; the decision rule split into same-nside and cross-nside parts.
- Profiled Fisher sd of a region index at nside_beta 16 on this sky spans 4e-3 to 0.61 (sync)
  and 5e-4 to 1.16 (dust), a precision range of 5e6 (E020 Result). In the non-centred
  coordinates the strongly informed regions are tied to mu and sigma and NUTS collapses (step
  size 2e-5, state frozen over 40 iterations even with a finite-difference Hessian mass start).
  The centred parameterization moves (step size 1e-2, 0 divergences after the first warm-up
  transitions). Non-centred remains right for the weak-data regime (E016 to E018).
- Two cost bugs fixed (same target): the graph-Laplacian eigenbasis (18 s per component at
  nside_beta 16) was rebuilt every Gibbs iteration on every rank (36 s of each iteration); the
  serving ranks computed the dense field-to-pixel maps for nothing. In-process gradient 47 ms at
  4 threads (66 ms at 1); the finite-difference start costs 12304 gradients per chain.
- New sampler options, all tested: `spectral_index_field: {type: iid | independent}`,
  `nuts.init_inv_mass: fd_hessian`; slow L1 test of the independent prior against exact
  per-region grids (PASS).

**Running**
- E020 production, sequential under nohup (`run_e020_prod.sh` in the session scratchpad; logs in
  `/Users/trygvels/data/e020_d1s1/run_<tag>.log`): gmrf16 started 10:51, then iid16, indep16,
  scalar, gmrf8, iid8, indep8, gmrf4, iid4, indep4; 2000 iterations each, two chains. The cost
  per iteration after warm-up decides whether the 2000-iteration budget stands; a waiter reports
  at chain-1 iteration 150.
- E016 hier SBC (L2 gate for M2.4a): 166 of 300 replicas.

**Next**
- Read the gmrf16 budget; shorten the schedule if a run exceeds about 12 h; analyze with
  `nextgen/scripts/e020_analyze.py`; second noise seed for the decisive pair; Result, Conclusion,
  F21, Q3, ROADMAP M2.5.
- SBC analysis when it finishes; E016 Conclusion; M2.4a done or withdrawn.

**Needs Trygve**
- Nothing new. Q15 (per-band zero levels) and the hier default (D003/D011) still open.

## 2026-10-27 (Mac mini, 05:40): M2.4 stage 1 validated at L3 (E018 PASS); M2.2 full-sky row done (E017); MAMS benchmarks nearly done

**Learned**
- E018 L3 PASS (F19): the 12-rank end-to-end Commander4 run with `likelihood: hier` (sync index
  field at nside_beta 2, K 47, plus field and amplitude hyperparameters; 3000 Gibbs iterations)
  agrees with the SBC-validated E016 prototype chain on the same maps in all 55 scalars: max 1.9
  pooled SE, sd ratios 0.86 to 1.09, 0 divergences, 15 gradients per iteration (same as the
  prototype), rms z against the simulation truth 1.10. Scalar configuration likewise (max 1.8 SE).
  Negative control: the same comparison separates the hier chain from the exact profile posterior
  by 19 MC SE and from the collapsed flat-prior grid by over 1000 sd. M2.4 stage 1 is validated at
  L1 and L3; the L2 gate is the E016 SBC (134 of 300).
- E017 final (F18, scalar index on template data): flat marginal biased +6 sd at full-sky S/N
  680, +18 sd at 230, collapsed to the support edge in all replicates at S/N 68 and 23; fixed
  Gaussian 4 to 10 sd off; profile unbiased everywhere but sd z 1.14 at S/N 68 and 1.90 at S/N 23
  (coverage 0.56/0.93 and 0.33/0.66); hier unbiased (|z| <= 0.3) and calibrated at all four
  scales, also with a 30 percent template distortion, hyperparameters recovered. Eight hier cells
  that had collapsed at the unguarded log s floor (workers predated the guard) were rerun with the
  guard; all healthy. The E014 marginal-minus-ridge gaps are this marginal-mode bias.
- E019 (M2.3, Q4) so far: NUTS non-centred 0.034 to 0.046 ESS per gradient (worst index pixel),
  flat in K; centred 10 to 60 times worse. MAMS non-centred 0.034 / 0.040 (nside_beta 2, noise 1 /
  0.3), 0.009 / 0.048 (nside_beta 4); centred 0.001 to 0.008. So MAMS is within a factor 1 to 1.4
  of NUTS on this target (one centred nside-4 cell left); no case for a second sampler yet.

**Running:** E019 last MAMS cell (centred nside_beta 4, to about 06:00); E016 hier SBC (134 of
300, about 5 per hour now that the machine is less loaded).

**Next**
1. E019 table and Conclusion, Q4 status, registry (when the last MAMS cell finishes).
2. E016 SBC analysis when done: F16 final or withdrawn; then M2.4a declared done.
3. Then M2.5 preparation (Q3: index resolution versus S/N on PySM skies) and the Q15 design note
   for Trygve; Paper A outline gets the E017/E018 numbers.

**Needs Trygve:** unchanged (Q15 per-band zero-level parameter design; hier stays non-default per
D003/D011, with E017/E018 now the evidence for making it the default below full-sky S/N 100 if
Trygve agrees).

## 2026-10-27 (Mac mini, 01:50): pass 3 and pass 4 grids done; M2.4 stage 1 code complete (commits 1 to 8, stats-reviewed)

**Learned**
- E014 pass 3 (20 bands): both chains exact against their grids (ridge -0.09 / +0.03 sd, marginal
  -0.03 / -0.02 sd, widths within 6 percent). Marginal minus ridge +16.1 +- 0.1 sd in beta_s and
  +2.1 in beta_d, against +56 / +5.5 on the pass-2 inputs: the flat-prior volume shift is itself
  data-set dependent.
- E014 pass 4 (one DR1-model zero-level convention for all 20 bands, grids only): ridge
  (-2.87545, 1.39914), marginal (-2.87236, 1.39924); against pass 3 the zero levels alone move
  beta_s by +0.035 (260 sd) and beta_d by -0.016 (460 sd). Marginal minus ridge +22.7 / +2.9 sd.
  Real-data full-sky indices are set by the zero-level convention at the 0.03 to 0.4 level, far
  above the statistical error; Q15 stands (per-band monopole parameter needed).
- M2.4 stage 1 implemented on nextgen, 8 commits, full suite 932 tests green: C++
  `solve_perpix_regions` and `want_extras` (bitwise parity with `solve_perpix`); array-valued
  SEDs and `beta_map` on the components; `AmpPrior` (none / gaussian / hier) as pseudo-bands and
  the prior- and region-aware amplitude draw; `RegionIndexLikelihood` (joint target over index
  fields and hier hyperparameters; stats review 2026-10-27 CONFIRMED against an independent
  covariance-form marginal, FD gradients in every block, 4-rank MPI reduction; its five findings
  applied); `IndexFieldSamplingGroup` (NUTS over MPI, state across Gibbs iterations);
  `likelihood: hier` with a `nuts` block, component keys `spectral_index_nside`,
  `spectral_index_field`, `amp_prior`; chain output `index_field/*`, `amp_prior/*`,
  `mcmc/<group>/{nuts,state}`; restart. Example `params/sims/simparam_compsep_hier_field.yml`.
  Not yet validated end to end (L1 to L3 below); the hier SBC of the prototype (L2 gate) is at
  86 of 300 replicates.

**Running:** pass-4 Haslam-dipole variant grid (01:40 to about 02:10); E017 hier modes (71 of 96
cells); E016 hier SBC (86 of 300, about 2.5 per hour: it will take another day on this machine).

**Next**
1. M2.4 commit 9: `aux/validate_mcmc_likelihood.py make --mode hier_field` (L3, 7 ranks, scalar +
   ridge must reproduce today's numbers; nside_beta 2 + hier against the prototype chain and
   per-region grids), the L1 in-process test, E018 record; stats-reviewer on the harness.
2. pass4hdip grid into E014 and Q15; FINDINGS entry on the zero-level dependence (F18); DECISIONS
   entry on the pass-4 convention and on the M2.4 start-state choice.
3. E017 table and Conclusion when the hier modes finish; E016 SBC result when done.

**Needs Trygve:** unchanged (Q15 per-band zero-level parameter; hier as a non-default likelihood
option is in, ridge stays the default per D003).

## 2026-10-26 (Mac mini, night): pass-3 ridge on 20 bands; zero-level inconsistency found (Q15); M2.4 plan drafted

**Learned**
- E014 pass 3 (all six PR4 HFI bands, 20 bands): ridge beta_s = -2.9103 +- 0.0001, beta_d =
  1.4150 +- 0.00003 (ESS 111 / 316), against -3.3425 / 1.3654 in pass 2 (17 bands): the full-sky
  indices move by 0.43 and 0.05, about 1250 posterior sd, when the HFI set changes. The chain is
  healthy (no bound, acceptance 0.04 because the pass-2 proposals are now too wide).
- Cause candidate, measured on the inputs: two zero-level conventions. WMAP/LFI had the
  |b| > 30 median removed (40.9 uK at 30 GHz, 78 uK at K); HFI were referenced to the DR1 CMB plus
  dust, and the DR1 CMB map has mean +28.6 uK at |b| > 30 (+41.9 at |b| > 60). The pass-3 100 GHz
  map sits at +45 uK at high latitude, the pass-2 one at 0. Commander4 has no monopole sampling, so
  a 29 uK_CMB step between the band sets goes into the foregrounds. E014's method results (chain
  versus grid on the same inputs) are unaffected; the index values are not science. Q15 opened.
- M2.4 implementation plan written and code-checked (`nextgen/plans/M2.4_commander4_index_field.md`):
  Gaussian and hier amplitude priors enter the existing region kernel as pseudo-bands (A' = A +
  T^-1, b' = b + T^-1 m), so one kernel serves likelihood, gradient and the exact amplitude draw;
  NUTS over the non-centred field plus hyperparameters; about 1000 source lines in 8 commits.
  Stats-reviewer on its section 1 running.

**Running:** E014 pass-3 marginal chain (23:16 to about 00:40) then the exact grid; E017 hier
modes (60 percent at 23:30); E016 SBC (66 of 100 replicates); pass-4 consistent zero levels
(subagent: `dr1_model_zero_levels.py`, offset = mean of map minus DR1 sky model at |b| > 30 for
all 20 bands); stats review of the M2.4 target.

**Next**
1. Pass-3 marginal and grid: record the 20-band ridge/marginal gap and the chain-versus-grid
   check in E014; then the exact grid on pass-4 inputs (minutes) to measure the zero-level
   sensitivity of both indices before any further chain.
2. E017 table and Conclusion; E016 SBC result and F16 final.
3. M2.4 stage 1 after the stats review (first commits: C++ `solve_perpix_regions` and the
   pseudo-band augmentation, with the bitwise parity test).

**Needs Trygve:** unchanged, plus (Q15) whether a per-band zero-level (monopole) parameter with a
prior should go into Commander4's per-pixel path as an offset column of the mixing matrix (it is
linear, so it collapses like an amplitude); without it, real-data spectral indices from the fork
depend on the input convention at the 0.4 level.

## 2026-10-26 (Mac mini, evening): E016 residual resolved; tuned legacy chain; scalar-index coverage (E017) running

**Learned**
- E016 residual: the +0.3 sd mean bias at noise 1 was a fluctuation of the 12 shared truths, not
  misspecification (Gaussian-scatter generator: +0.34) and not sampling (non-centred, ESS 394:
  +0.30). The pooled pixel SE (0.04) was wrong because pixels of a replicate share truth and noise;
  the replicate-clustered SE is 0.15. 24 fresh replicates give +0.01 +- 0.07; all 36: +0.10 +- 0.07,
  coverage 0.71 / 0.975. `coverage_analyze.py` prints the clustered SE for every beta row; all
  summaries regenerated (E012 conclusions unchanged: gauss at noise 1 is +2.08 with clustered SE
  0.11). F16 stays provisional only on the SBC.
- E014 tuned legacy chain (7.3 h, 1000 iterations): proposals sized for the ridge width give ESS
  17 to 20 and a width ratio above 11 (lower bound, chain has not traversed its posterior); the
  untuned chain's 16 to 17x stays the measurement. Tuned ridge: chain and exact grid agree to 0.06 sd
  and 5 percent in width.
- E017 (M2.2 scalar index): new `scalar_coverage.py` computes the w0 / w1 / gauss posteriors of one
  full-sky index exactly on adaptive grids and samples the hier target (5 parameters) with NUTS;
  pinned-hier-versus-gauss-grid check passes (ESS 618, -0.02 sd, width ratio 1.07). First cells: at
  full-sky S/N 23 the flat-prior marginal piles up at the boundary of the support (z +1700), at
  S/N 220 it is +17 sd off, gauss +10 sd, w0 within 0.3 sd.

**Running:** E017 grid modes (24 replicates x 4 noise scales, 2 workers, half done at 19:40) and
hier modes (launched 19:38, 2 workers). E016 SBC non-centred (about 25 of 300 at 19:40; speeds up
when the E014 chains end). E014 tuned `conditional` (since 17:44), then `marginal`.

**Next**
1. E017 table and Conclusion (expected: w0 and hier calibrated at full-sky S/N, w1 and gauss not;
   quantify the w1 shift in w0 sd against the 54 sd of E014 on real data); F17 if it holds.
2. E016 SBC result, Conclusion, F16 final; E014 tuned conditional / marginal, Conclusion.
3. Commander4 implementation plan for the hierarchical amplitude prior (D-entry: where mu_c, s_c
   live in the parameter file, which template on real data), stats-reviewer on the target first.
4. Paper A outline.

**Needs Trygve:** unchanged (Cosmoglobe contact on F1, F10, F11, F13 and the Jeffreys default; Q13;
cluster; DR2 parameter file; ridge default for spatially varying indices with F16 as the candidate
replacement).

## 2026-10-26 (Mac mini, afternoon): hierarchical amplitude prior calibrated at every S/N on template data

**Learned**
- E016 coverage study (108 cells, same replicates and S/N bins as E012 study_v1): with the
  amplitude scale and relative scatter per component learned from the data, the index-field
  posteriors are calibrated at every S/N tested: mean z +0.10 / +0.35 / +0.24 sd over all pixels at
  noise 0.3 / 1 / 3 with coverage 0.66 / 0.95, 0.64 / 0.95, 0.67 / 0.93, against w0 z +1 to +9 and
  gauss z +14 / +2 / +0.6 on the same replicates. In the bins below S/N 1 (where w0 is biased by 8
  to 13 sd): z -0.02 and +0.35, coverage 0.56 / 0.96 and 0.67 / 0.91. Shape templates wrong by 30
  percent (large scale) or 23 percent (per pixel) change nothing: the learned scatter absorbs it
  (log s recovered to its posterior sd). F16, provisional until the SBC passes. Residuals: +0.3 sd
  mean bias at noise 1 in all variants; one 58-pixel bin with coverage 0.48 / 0.86 (same as w0).
- F15 (joint flat-prior posterior has the w1 collapse as its exact beta marginal) confirmed by the
  stats-reviewer and recorded; wording restricted to the flat-prior joint, not Commander3/4
  production targets.
- E016 L1 checks all pass (direct form 7e-13, gradients 4e-9, w1 limit scales as 1/s^2, gauss limit
  exact); `HierNonCentredPosterior` added (same target, Jacobian identity to 5e-14) because the
  centred SBC at K 192 could not complete one replicate in 4.5 CPU-hours per worker (stopped,
  recorded in E012; not needed: the centred parameterization is used only at K 48).

**Running:** E014 tuned round: `resolve_chisq` (since 10:27, 667 iterations of 1000 at 15:30),
then `conditional` and `marginal`. E016 SBC of the hier target (non-centred, 300 replicates, 4
workers) starts after its smoke test; about 10 h at the current load.

**Next**
1. E016: SBC result, Conclusion, F16 final or withdrawn; understand the +0.3 sd residual at noise 1
   (log-normal versus Gaussian scatter is the first suspect: rerun one noise-1 cell with Gaussian
   generative scatter).
2. E014: tuned legacy / conditional / marginal results, Conclusion, registry row.
3. M2.2 scalar-index coverage table with the hier mode included; then the Commander4
   implementation plan for the hierarchical amplitude prior (D-entry: where mu_c, s_c live in the
   parameter file, which template on real data) with stats-reviewer on the target before code.
4. Paper A outline: F1/F2 (sampler), F3/F12/F16 (priors), F10/F11 (Commander3), F13 (DR1 dust).

**Needs Trygve:** unchanged (Cosmoglobe contact on F1, F10, F11, F13 and the Jeffreys default; Q13;
cluster; DR2 parameter file; ridge default for spatially varying indices, now with F16 as the
candidate replacement: learned amplitude prior rather than profile below S/N 15).

## 2026-10-26 (Mac mini, midday): both exact modes verified on real data; coverage study: no safe prior-free option for a varying index below S/N 15

**Learned**
- E014: the `marginal` chain lands on the exact flat-prior marginal grid target within 0.05 sd
  (width ratios 1.08 / 1.23 at ESS 76 / 45), and the tuned `ridge` chain (ESS 559 / 422) matches
  its target to 0.06 sd with width ratios 0.96 / 0.94. Two exact chains measure the prior effect
  directly on the data: 54 ridge sd in beta_s, 5 sd in beta_d, same widths, chi-square differs by
  0.07 percent. Q1 is answered at nside 64 and the ranking changed (D010: Q2 first).
- E012 coverage study (108 cells, 12 replicates x noise 0.3 / 1 / 3 x gauss / w0 / w1): the
  flat-prior marginal (w1) is unusable at every S/N (z 100 to 300, zero coverage); the profile (w0)
  is calibrated at index-pixel S/N above about 15 (coverage 0.70 / 0.96) and biased with heavy tails
  below S/N 3 (z +5 to +13, field sd biased +5 to +28 sd, ESS 4 to 6); the Gaussian amplitude
  prior is calibrated only when it is the generative one (diagnostic: z +0.08 +- 0.06, coverage
  0.69 / 0.95), biased by 2 to 14 sd on template data with the default width, and 4x worse with a
  10x wider prior, because the Gaussian marginal tends to the w1 volume term as T grows. F12 final.
  There is no uninformative safe amplitude prior for a spatially varying index: learn it from the
  data (hierarchical amplitude model) or use the profile at an index resolution set by S/N.
- A corollary (joint posterior with flat amplitude priors has the w1 collapse as its exact beta
  marginal; Gibbs chains avoid it only by not mixing) is with the stats-reviewer before it is
  claimed anywhere.

**Running:** E014 tuned round: `resolve_chisq` (since 10:27, about 5 h), then `conditional` and
`marginal` (about 1 h each). E012 centred SBC (300 replicates, since 10:07, slow under load 28 on
10 cores). Stats-reviewer on the Gibbs corollary.

**Next**
1. E014: tuned legacy / conditional / marginal results, Conclusion, registry row; then IQU and the
   nside 128 scaling check (lower priority).
2. E012: centred SBC analysis; record the reviewer's verdict; then the hierarchical amplitude-prior
   experiment (Q2 next step) with its target written and reviewed first; M2.2 scalar-index coverage
   table.
3. Paper A outline: F1/F2 (sampler), F3/F12 (priors), F10/F11 (Commander3), F13 (DR1 dust).

**Needs Trygve:** unchanged from the morning entry (Cosmoglobe contact on F1, F10, F11, F13 and the
Jeffreys default; Q13; cluster; DR2 parameter file). New for the ridge-default decision (D003): on
a spatially varying index the profile is safe only above index-pixel S/N 15; the default for
Phase 2 will need either an S/N-driven index resolution or a learned amplitude prior. No default
changed yet.

## 2026-10-26 (Mac mini, morning): Q1 answered at nside 64; index-field SBC passes; DR1 dust header error

**Learned**
- Q1 (E014, pass 2 on NPIPE 353/857 + PR3 100 + DR1 WMAP/LFI + Haslam, nside 64, dx12 mask): the
  legacy `resolve_chisq` index posterior is 16.5x (beta_d) and 17.1x (beta_s) too wide against an
  exact grid evaluation of the ridge target on the same data, and beta_d is shifted by +3.0 +- 0.5
  ridge sd. The ridge chain itself matches the grid within 0.2 sd (width ratios 1.01 / 1.18). The
  inflation follows the F1 scaling (k^(1/4) predicts 18.5x); nside 512 to 1024 would give above 50x.
- The flat-prior marginal target on the same data sits 54 ridge sd from the ridge target in beta_s
  (+0.019): the amplitude prior (Q2) matters more than any sampler detail on real low-resolution maps.
- E012 L2: SBC of the hierarchical index-field sampler (gauss, non-centred, 300 replicates, 14 h)
  PASSES on all 12 quantities (min p 0.055, Bonferroni, ECDF band, TARP); M2.1 acceptance met for
  the default configuration.
- Commander3 branch survey (E009, E015): all four Commander3 defects are present at the master tip
  (2026-03-17) and on every active branch (akari, spider, hfi, hfi_Tbol, litebird). The marginal
  sign is corrected only on an unmerged test branch; `hfi_Tbol`'s 2026-10-23 MH rewrite still
  accepts on chi2 after CG. Our earlier readings were on dr1_repro-era code; conclusions unchanged.
- F14 (E014): the Gibbs `conditional` chain cannot mix on real data (conditional sd of beta_s 2e-6
  against marginal 3.6e-4); the collapsed modes do.
- F13: the Cosmoglobe DR1 dust map is referenced at 70 GHz, not the 545 GHz in its header
  (parameter file `COMP_NU_REF_T03 = 70 857`, Commander3 keeps the first token); the BeyondPlanck
  dust file has the same header error. The truncated local DR1 synchrotron file was re-fetched.

**Running:** E014 `conditional` (since 05:56) then `marginal` (untuned, 2000 iterations), then the
tuned-proposal round for all four modes (1000 iterations each). E012 coverage study (gauss / w0 /
w1 x noise 0.3 / 1 / 3, 12 replicates, 4 jobs) starting as the SBC released its workers; the centred
SBC is queued behind it.

**Next**
1. E014: conditional and marginal against the exact marginal target; tuned round; write the Q1
   conclusion and the E014 record; then IQU and a nside 128 scaling check (lower priority).
2. E012: coverage analysis (promote or drop F12), centred SBC; then M2.2 scalar-index coverage table
   (Q2) with the same code.
3. Re-rank OPEN_QUESTIONS (Q1 answered; Q2 now has a real-data number) and log in DECISIONS.

**Needs Trygve:** the Q1 answer sharpens the Cosmoglobe question (F1 is a 16x effect on real data at
nside 64; the `perpix_betasamp` example config runs at nside 512). Unchanged otherwise (F1 and the
sign error, F11, the disabled Jeffreys prior, Q13, cluster access, DR2 parameter file, ridge
default revisit once F12 is quantified). New: F13 is a documentation defect in a public data
release; whether to report it to Cosmoglobe.

## 2026-10-25 (Mac mini, night): ridge field pathology found; Q1 pass 2 switched to NPIPE

**Learned**
- Ridge (w0) chains in the index-field coverage harness mix at ESS 4 to 5 per 1500 s in all four
  parameterization and initialization settings, and infer field sigma about 1.2 against a truth of
  0.065 (z +22). A profile scan shows the cause is the target: per-region profile likelihoods are
  rough and multimodal at S/N below a few (ML scatter rms 2.2 about truth). The flat-prior marginal
  collapses every region to beta about -0.1 (sd 0.01). The proper Gaussian amplitude prior behaves.
  Recorded as F12 (preliminary, one replicate) in E012.
- Q1 pass 2 inputs: the Planck archive drops PR3 transfers above about 60 MB for 143 and 217 GHz
  (100 GHz verified). The second pass uses the NPIPE A-split 353 and 857 GHz maps hosted at UiO
  (BeyondPlanck inputs) plus PR3 100 GHz, dx12 chi-square mask, proposals 0.002 / 0.004, numstep 20.

- E014 pass 2 ridge (NPIPE 353/857 + PR3 100, dx12 mask, 2000 iterations, 1.6 h): beta_d 1.36539
  +- 0.00005, beta_s -3.3425 +- 0.0004, but acceptance 0.4 percent (proposals 80x / 6x the posterior
  sd), ESS 44 / 97. Means are solid, widths to 15 percent. A tuned-proposal round (1.2e-4 / 9e-4,
  1000 iterations, tag `npipetuned`) is queued for all four modes.
- DR1 dust amplitude anomaly resolved (F13): the DR1 dust map is referenced at 70 GHz, not the
  545 GHz in its FITS header (parameter file `COMP_NU_REF_T03 = 70 857`, Commander3 keeps the first
  token). Predicted 857 / 353 GHz slopes 11.3 / 7.3 against 12.2 / 7.2 measured. Zero levels
  unaffected; prep script fixed; the BeyondPlanck dust header has the same error.

**Running:** E014 pass 2 untuned driver (resolve_chisq since 20:08, then conditional, marginal;
about 2 h each), then the tuned round (four modes, about 1 h each). SBC default, 96 of 300
replicates at 20:12 (load 27, rate dropping toward 20 per hour). Coverage study queued behind the
SBC run.

**Next**
1. E014 pass 2 analysis (analyze_q1_chains.py --suffix _npipemask, then _npipetuned) as modes
   finish; resolve_chisq is the first comparison against ridge.
2. SBC analysis (sbc_analyze) when the 300 replicates finish; then the centred SBC.
3. Coverage analysis; promote or drop F12; M2.2 scalar-index coverage table.

**Needs Trygve:** unchanged (Cosmoglobe contact about F1, the sign error, F11, the disabled
Jeffreys prior; Q13; cluster access; the DR2 parameter file). New: whether the ridge default of the
fix branch (F3, F4) should be revisited for spatially varying indices once F12 is quantified.

## 2026-10-25 (Mac mini, evening): index-field prototype passes L1; SBC running

**Learned**
- M2.1 prototype (`nextgen/prototypes/indexfield/`, E012 v2): collapsed likelihood with gauss / w0 /
  w1 modes, GMRF and harmonic Matern fields (constant mode removed, logit rho), NUTS; an
  independently written reference path (plain-loop likelihood, own bases, centred Gibbs, closed
  forms). L1 complete: 1a gradients (1.7e-8), C++ kernel agreement (7e-12); 1b exact SED grids,
  12 cases, pooled ESS 3450 to 10700, |z| <= 2.1; 1c hyperparameters vs Gibbs reference, |z| <= 1.7
  at ESS > 800; 1d closed form vs NUTS with independent bases, 96000 draws, |z| <= 3.5; 1e all
  three wrong targets fail as expected (dropped logit Jacobian makes the target improper).
- Centred parameterization mixes 10x better than non-centred at strong data (E-BFMI 1.0 vs 0.8).
- Flat-prior marginal (w1) collapses the index field to beta ~ -0.2 even at S/N 30 (field version
  of F3); ridge and the proper Gaussian prior behave. To be quantified in the coverage study.
- Independent review of the implementation: correct with caveats; caveats fixed (no replicate ever
  dropped from the SBC, singular-A handling, dispersed starts).

**Running:** SBC, default configuration (gauss, GMRF nside_beta 4, non-centred), 300 replicates,
4 workers (`sbc_out/`). HFI PR3 maps from the Planck Legacy Archive (1 of 6 verified; the archive
does not resume, so each file must complete in one transfer). Q1 pass 2 armed behind it.

**Next**
1. SBC analysis; then the centred-parameterization SBC and the coverage study (w0 chains need a
   budget or geometry decision: ESS 3 to 4 at 300 + 400 transitions).
2. E014 pass 2 (HFI, Galactic chi-square mask, narrow proposals).
3. FINDINGS F12 (w1 collapse) once quantified; M2.2 scalar-index coverage table with the same code.

**Needs Trygve:** unchanged (Cosmoglobe contact about F1, the sign error, F11, the disabled
Jeffreys prior; Q13; cluster access; the DR2 parameter file).

## 2026-10-25 (Mac mini, afternoon): Commander3 reviews done, Q14 answered for BeyondPlanck, real maps load

**Learned**
- Ten independent stats-reviewer passes (two per item, code paths and stated targets only) confirm
  E009/E010: BP-era index samplers valid conditional steps, `marginal` sign flip (quadrature by
  both), 2024 MH module has F1 (toy inflation 1.3x to 4.6x). Corrections: BP10's Jeffreys prior was
  hard-disabled in the code; BP10 sampled AME nu_p as one full-sky value with monopoles profiled.
- New defect F11 (E015, verified line by line): BP's absolute bandpass chi^2 is never MPI-reduced,
  so the 30 GHz absolute shift was accepted on rank 0's scans only. Present at HEAD and in the
  WMAP/DIRBE/SPIDER/LiteBIRD/QUIET TOD modules.
- Q14 (E013): the sign error moved BP's nu_p by 0.008 GHz (0.8 of a 0.01 GHz conditional sd);
  negligible because the masked AME pixels have S/N above 30. DR1 never used the `marginal` path.
- E012 target revised (v2) after review: constant mode dropped, radians, logit rho, SBC with a
  proper Gaussian amplitude prior (D008, D009); F4 wording corrected (conditional Jeffreys on a).
- M1.3: DR1 WMAP+LFI, Haslam and (partially) PR4 maps are local (~15 GB); Commander4 maps-only runs
  on them at nside 64 after dipole and monopole removal (reduced chi^2 ~150 with spinning dust).

**Running:** E014 (Q1) first pass without HFI: ridge, resolve_chisq, conditional, marginal, 2000
Gibbs iterations each, ~50 min per mode, 14 ranks. PR4 (NERSC) and PR3 (IRSA) HFI downloads, slow
(~0.3 MB/s).

**Next**
1. M2.1: implement the index-field prototype per E012 v2 (L1 checks 1a to 1e, then the SBC).
2. E014 analysis when the four chains finish; then the HFI variant and IQU.
3. M2.2 coverage study with the same prototype at N_beta = 1.

**Needs Trygve**
- Whether to tell Cosmoglobe about F1, the `marginal` sign error, and now F11 (bandpass reduce)
  and the disabled Jeffreys prior. All are code-reading findings with written derivations.
- Unchanged: Q13 upstream pull request, cluster access for T3/T4.

## 2026-10-25 (Mac mini): M1.1 and M1.2 done

**Learned**
- The nextgen branch builds and passes on the Mac mini (Apple M4, 10 cores, 24 GB, macOS 27.0,
  Apple clang 17, CPython 3.13.11, pip MPICH 5.0.2): 859 passed, 4 skipped in 25 s (E011).
- The macOS setup script needed no Homebrew for the build: every dependency has an arm64 wheel,
  cmake and ninja install as pip wheels, gfortran is unnecessary. The script was fixed accordingly
  and the machine record and fixes are in COMPUTE.md. tmux and Julia (juliaup) were installed with
  Trygve's approval; the specidx-bench prototype project is instantiated.
- E005 reproduces to the quoted precision (resolve_chisq 3.22x / 4.35x too wide at 150 / 600
  pixels, mean shifts +0.50 / +0.83 sd; marginal, ridge and conditional within 0.1 sd and 10% in
  width). E006 reproduces: ridge within 0.03 sd of the exact ridge target, width ratios 0.94 / 0.97;
  exact flat-prior target on the bounds; marginal stuck in the beta_s = -3.5 mode. F1, F2, F3 and F8
  therefore hold on an independent toolchain. One E006 mode takes about 80 s on 7 ranks here (15
  min in the cloud); resolve_chisq takes 5 min.
- Disk is the binding constraint on this machine (~100 GB free): low-resolution maps only; T2 CG
  runs and anything at nside 1024 go to a cluster (COMPUTE.md tiers adjusted).

**Running:** nothing.

**Next**
1. Independent `stats-reviewer` pass on E009 and E010 (Commander3 code paths and stated targets,
   not our conclusions). Clone Commander3 next to this repository first (KICKOFF.md).
2. Q14: rerun the BeyondPlanck / DR1 AME nu_p step with the sign corrected and with `ridge`.
3. M1.3 / M1.4 (Q1 on real low-resolution maps, within the disk budget) and M2.1 in parallel.

**Needs Trygve**
- Unchanged: Q13 (upstream pull request), whether to tell Cosmoglobe about F1 and the `marginal`
  sign error, cluster access for T3/T4.

## 2026-10-25 (later): Commander3 review, handover to the Mac mini

**Learned**
- Commander3's BeyondPlanck-era polarized-index sampler (alm / pixel regions) and bandpass sampler
  are valid conditional MCMC steps in the BP10 configuration (E010). Main approximations: prior
  draws for frozen regions, diagonal spurious-map chi^2 for relative bandpass shifts.
- Commander3's per-pixel `marginal` likelihood has the log-det sign flipped (+1/2 ln A instead of
  -1/2 ln A). BeyondPlanck final used it for the AME peak frequency (E009, F10).
- Commander3's 2024 MH module (`sample_specind_mh` and siblings) has the same acceptance bug as
  Commander4 when CG groups are updated inside the MH loop (E009).

**Running:** nothing.

**Next**
1. M1.1 / M1.2 on the Mac (build, tests, reproduce E005 and E006).
2. Independent `stats-reviewer` pass on E009 and E010 (give it the Commander3 code paths and the
   stated targets, not our conclusions). Required before anything about Commander3 leaves this repo.
3. Q14: rerun the BeyondPlanck / DR1 AME nu_p step with the sign corrected and with `ridge`.
4. Then the plan in the entry below (Q1, M2.1, M2.2).

**Needs Trygve**
- Whether and when to tell Cosmoglobe about F1 (Commander4 and the 2024 Commander3 MH module) and
  the `marginal` sign error.

## 2026-10-25: handover from the exploratory cloud session

**Learned**
- Commander4's spectral-index MCMC move is not a valid MCMC step; posteriors 1.2x to 5.7x too wide
  in tests, growing with resolution (F1). Fixed and validated end to end (F2).
- The flat-amplitude-prior marginal can be dominated by its volume term; ridge (profile, joint
  Jeffreys) is the default now (F3, F4). This prior question is open and central (Q2).
- Collapsing amplitudes beats Gibbs + MH by 9x to 22x in amplitude solves per effective sample
  (F5). Language is not a bottleneck (F6).
- Two per-pixel model issues: alm storage of amplitudes (F7) and smoothed noise treated as white
  (F8).

**Pushed to `trygvels/Commander4`**
- `fix/mcmc-marginal-acceptance` (PR-ready, not opened upstream).
- `feat/spatially-varying-beta-foundation` (region kernel with gradients, NUTS).
- `nextgen` (this scaffolding, prototypes, research reports).

**Running:** nothing.

**Next** (see ROADMAP Phase 0 and 1)
1. M1.1 and M1.2: set up the Mac (`nextgen/scripts/setup_macos.sh`), fill in the machine record in
   COMPUTE.md, run the test suite, reproduce E005 and E006.
2. M1.3 and M1.4: fetch low-resolution DR1 / PR4 / WMAP maps and answer Q1 (real-data impact of F1).
3. M2.1 in parallel: prototype the hierarchical GMRF index field outside Commander4.
4. M2.2: SBC coverage study of ridge vs flat-prior marginal vs Jeffreys (Q2).

**Needs Trygve**
- Q13: whether and when to open the fix as an upstream pull request.
- Access: ITA / Sigma2 cluster account and allocation for T3/T4 work.
- Which Mac mini spec (memory, cores, disk) is available.
