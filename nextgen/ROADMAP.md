# Roadmap

Phases overlap. Each milestone has acceptance criteria; a milestone is done only when its criteria
are met and logged in `experiments/registry.md`. Revisit this file when findings change the plan.

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
| M2.2 Answer Q2 (prior calibration) | Coverage table for marginal / ridge / conditional / Jeffreys, full-sky and spatially varying |
| M2.3 Sampler choice (Q4) | ESS per gradient table for NUTS vs MAMS with 2 to 3 preconditioners |
| M2.4 Implement in Commander4 per-pixel path: index maps on components, per-region mixing in `solve_perpix`, chain output, config, tests | Full suite green; end-to-end runs match exact targets (ladder L3) |
| M2.5 Compare against alm GRF, clustering, fixed ladder (Q3) on PySM d1s1 / d10s5 / d12s7 | Table of index error, coverage, CMB residuals, cost |
| M2.6 Moments (Q8) and smoothed-noise model (Q7) | Quantified effect; decision logged |

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
