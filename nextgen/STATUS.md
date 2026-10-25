# Status

Newest entry on top. Keep each entry short: learned (with numbers), running, next, needs Trygve.

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

**Pushed to `trygvels/commander4`**
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
