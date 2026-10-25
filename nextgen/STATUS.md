# Status

Newest entry on top. Keep each entry short: learned (with numbers), running, next, needs Trygve.

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
