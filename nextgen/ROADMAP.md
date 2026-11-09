# Roadmap

Phases overlap. Each milestone has acceptance criteria; a milestone is done only when its criteria
are met and logged in `experiments/registry.md`. Revisit this file when findings change the plan.

## Proposed revision after the 2026-11-09 review (D029, pending Trygve)

The phases below are the plan of 2026-11-02 and are kept as written. The review
(`reviews/2026-11-09_program_review.md`) proposes this order instead:

| Step | Weeks | Acceptance |
|---|---|---|
| R1. E032 read; no follow-ups | 0 to 1 | Pre-registered items of the E032 record, four chains |
| R2. E036 completed (noise ladder, off-centre truths, free CMB and dust, a CMB prior, the global-scale arm) | 0 to 1 | Predictions in the E036 record |
| R3. Gate 0 for polarization: WMAP K-band res 4 map and inverse covariance fetched and checked; an LFI 30 GHz low-resolution covariance that passes a half-split null test; one window for WMAP and LFI | 1 | Challenge report Sect. F, Gate 0 |
| R4. K band against 30 GHz null test: region indices from {K, Ka}, {K, 30}, {30, 44, Ka}, profile with sandwich error, dense covariances; WMAP9 with PR3, PR4 and DR1 30 GHz | 1 to 3 | The three pairs agree within 2 combined sd, or a model with at most five named modes brings the residual to a PTE above 0.01, or the disagreement is reported as the result |
| R5. Paper A1 draft: legacy move, incidental parameters, limits, polarization section | 1 to 6 | Every claim A or B grade; exact grids; the 20-seed F8 test; SBC of the Commander4 ridge and marginal modes |
| R6. Commander4 fixes after E032 (reader order and floor behind a switch, `get_sky` with an index field, maps instead of alms for per-pixel amplitudes) and the upstream reports | when cores are free; reports need Trygve | Tests green; E014 grid rerun if the real-intensity section is kept |
| R7. Paper B, planned on the outcome of R4: an index field at nside 8 data and nside 2 field on Gates 1 to 6, or a systematics paper | after R4 | Challenge report Sect. F |

**Paper A1 (methods):** spectral indices and incidental amplitudes in pixel-based Bayesian
component separation. **Paper B (science):** the polarized synchrotron index at the resolution
the public data determine, positioned against Dunkley et al. 2009. **tau:** a paired difference
between index models on the same WMAP and LFI data, as a by-product of Paper B.

## Phase 0: groundwork (done, 2026-10-25)

- Correctness fix for spectral-index MCMC with exact likelihood options and a ridge default.
- Per-region likelihood kernel with analytic gradients; NUTS sampler.
- Literature survey (four reports in `research/`), findings, open questions, validation protocol.

## Phase 1: environment and real-data check (weeks 1 to 2)

| Milestone | Acceptance |
|---|---|
| M1.1 Build and full test suite pass on the Mac mini | All tests green; setup issues recorded in COMPUTE.md |
| M1.2 Reproduce E005/E006 on the Mac | Same numbers within MC error |
| M1.3 Fetch DR1 / PR4 / WMAP maps at low resolution (DATA.md) | Maps load in a Commander4 maps-only config |
| M1.4 Answer Q1 at nside 64 to 128 | Report of resolve_chisq vs ridge vs conditional on real maps |
| M1.5 Decide on the upstream pull request (Trygve) | Decision logged in DECISIONS.md |

## Phase 2: optimal spatially varying spectral parameters (weeks 2 to 10)

| Milestone | Acceptance |
|---|---|
| M2.1 Prototype hierarchical GMRF index field (outside Commander4, NumPy/JAX or Julia) on per-pixel sims | Exact-target checks pass; SBC ranks uniform for beta field and hyperparameters |
| M2.2 Answer Q2 (prior calibration) | Coverage table for marginal / ridge / conditional / Jeffreys, full-sky and spatially varying. Done 2026-10-27 (E012, E016, E017; F12, F16, F18): hier prior calibrated, profile only at high S/N, flat marginal never |
| M2.3 Sampler choice (Q4) | ESS per gradient table for NUTS vs MAMS with 2 to 3 preconditioners. Done 2026-10-27 for the diagonal preconditioner (E019, F20): parity, NUTS stays; other preconditioners deferred to Q5 |
| M2.4 Implement in Commander4 per-pixel path: index maps on components, per-region mixing in `solve_perpix`, chain output, config, tests | Full suite green; end-to-end runs match exact targets (ladder L3). Stage 1 done 2026-10-27 (E018, F19): L1 and L3 pass; L2 gate passed 2026-10-29 (E016 SBC, 600 replicates, F16) |
| M2.5 Compare against alm GRF, clustering, fixed ladder (Q3) on PySM d1s1 / d10s5 / d12s7 | Table of index error, coverage, CMB residuals, cost. 2026-10-28: first cut E020 on the PySM sky closed (F21, misspecification dominates at 0.03 x noise); the comparison runs as E021 on a model-consistent sky first, PySM skies return under Q16. 2026-10-29: E021 done (F24): GMRF beats iid and independent regions at nside 16, 8, 4 on the model-consistent sky; the nside-16 GMRF keeps most of the per-region averaging of coarse regions (stats-reviewed, D019); GMRF is the prior of record (D018). Remaining for this milestone: alm GRF and edge-preserving variants, clustering baselines, PySM skies (Q16) |
| M2.6 Moments (Q8) and smoothed-noise model (Q7) | Quantified effect; decision logged. 2026-10-30: sampled T_d per region scheduled here as the structural fix for the dust template error (E023, F27, D022). 2026-11-02: code for a sampled log T (scalar or field) and a fixed T map in place (E028); staged as oracle T map, scalar T, then a field with 545 and 857 GHz (D026) |

**Paper A (methods):** exact spectral-parameter sampling in Bayesian component separation; the
prior-volume problem and the ridge/Jeffreys resolution; hierarchical GMRF index fields with inferred
smoothness; collapsed sampling efficiency. Target journal A&A or JCAP.

## Phase 3: beams, real data, degeneracy breaking (months 2 to 6)

| Milestone | Acceptance |
|---|---|
| M3.1 Extend collapse to CG amplitudes (Q5): low-res likelihood or stochastic log-det | Bias and efficiency quantified on sims |
| M3.2 Ingest low-frequency data (Q9): C-BASS, QUIJOTE MFI, S-PASS (with Faraday handling), CLASS 40 GHz where public | Each dataset validated alone, then jointly |
| M3.3 Per-band calibration/monopole nuisances where missing in Commander4 | Implemented and tested |
| M3.4 Science run: polarized synchrotron spectral-index field from WMAP + LFI (+ external) with full uncertainties | Posterior maps; comparison with BeyondPlanck XV, Weiland 2022, QUIJOTE, CLASS |

**Paper B (science):** the polarized synchrotron spectral-index field at the resolution the data
support, with calibrated uncertainties, and its consequences for foreground cleaning of large-scale
polarization (and AME polarization limits).

## Phase 4: performance and cosmology (months 4 onwards)

- CG cost (Q12): recycling, warm starts, multilevel preconditioners; target 10x fewer SHTs per
  independent sample.
- MAMS/MCLMC for index fields and possibly joint amplitude moves.
- C(l) sampling with rescaling and overrelaxation (Q11).
- Toward tau: Gibbs correlation lengths (Q10), HFI systematics templates (ADC nonlinearity,
  polarization efficiency, bandpass mismatch, transfer functions), coordinated with Cosmoglobe's
  OpenHFI (Q13).

**Paper C (long-term, high impact):** end-to-end tau with full error propagation, addressing the
low-l (~0.055) versus DESI+CMB (~0.090) tension.
