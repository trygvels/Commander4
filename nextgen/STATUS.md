# Status

Newest entry on top. Keep each entry short: learned (with numbers), running, next, needs Trygve.

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
