# Validation ladder

Every new sampler, likelihood, prior or model component climbs this ladder in order. A level is
passed only with a logged experiment. Efficiency claims are allowed from level 1 upwards, but only
in units of effective samples per dominant cost.

## L0: unit tests
Kernels against independent references (NumPy, finite differences); bookkeeping (state restore on
reject, when amplitudes are redrawn); config validation. Part of the pytest suite, fast.

## L1: exact-target statistical tests
Small problems where the target is known exactly: analytic Gaussians, 1D/2D grids of the exact
log posterior computed independently of the code under test, or an independent valid sampler.
Pass: means within 4 MC standard errors (using ESS), widths within ~15 to 20%. Include at least one
negative control (a known-wrong method must fail the same test, as `resolve_chisq` does).
Chain-length L1 tests live in the pytest suite behind `C4_SLOW_TESTS=1` (they take minutes; e.g.
`tests/test_index_field_l1.py`, E018) and are run explicitly when the experiment record says so.
A reference sampler is only a reference where it mixes: check its ESS and, when a parameter has a
funnel (a scale with a long tail towards zero), compare where the data identify it or use a
sampler that reaches the tail (E018 L1 note).

## L2: simulation-based calibration
Draw parameters from the prior, simulate data from the model, sample, record the rank of the truth
among posterior draws. Ranks must be uniform (Talts et al. 2018). Use the joint log likelihood as an
extra test quantity to catch posteriors that equal the prior (Modrak et al. 2022). For
high-dimensional fields use TARP coverage (Lemos et al. 2023). Low resolution (nside 8 to 32), a few
hundred replicates. Required for every prior choice in Q2/Q3.

## L3: end-to-end Commander4 runs on model-consistent simulations
Full Commander4 (MPI, real I/O and chain writing) on simulations drawn from the model being sampled,
compared to the exact or independently computed target. Harness: `aux/validate_mcmc_likelihood.py`
(MH likelihood modes; `--mode hier_field` for the index field with a hierarchical amplitude prior,
E018). Mind F7/F8: per-pixel runs need common-resolution smoothing, and smoothed noise is not white.

## L4: misspecified realistic skies
PySM 3 skies of increasing complexity (d1s1, d9s4, d10s5, d12s7, and the PEGS skies of 2502.20452),
realistic noise and beams. Metrics: index-map error and coverage, CMB residual power (TT/EE/BB at
l < 100), r bias where relevant, chi^2 and posterior-predictive maps, cost. Keep the BeyondPlanck IV
distinction: Bayesian (constrained) simulations test self-consistency; frequentist (random-input)
simulations test bias and coverage.

## L5: real-data regression
Reproduce, or explain differences from, published products: Cosmoglobe DR1 maps and chains,
BeyondPlanck synchrotron regions, Planck PR4 component maps. Differences are findings, but each must
be traced to a cause (method, prior, model, data) before it is claimed.

## Reporting checklist for any result
- Target distribution and priors stated explicitly.
- Exact or independent reference used, and how it was computed.
- Number of independent samples (ESS) behind every mean and width.
- Negative control, where one exists.
- Experiment ID, commit hash, command line.

## Standard checks on production chains (D025)
- Gibbs-consistency (conditional CDF) check: for K >= 20 post-burn draws, scan the full target
  along one coordinate (at least the field mean of every index field) with every other
  coordinate at the draw's value; the conditional CDF at the draw must be consistent with
  uniform (KS test) and the conditional optimum with the draw. Exact for the implemented
  target, independent of the sampler; 0.1 s per evaluation at nside 64
  (`nextgen/scripts/e026_conditional_check.py`). A chain that fails is not read. Stats-review
  additions (E026): K >= 100 draws (K = 20 detects only KS D > 0.3), a scan step of at most
  one conditional sd, and a negative control (a deliberately corrupted target or chain must
  fail the check) before the check is quoted as evidence. Negative control done (E027, 2026-11-02,
  `E027/negctrl_gmrf4_*.txt`): a narrowed prior on the field mean (optimum 1.4 conditional sd
  off), 1.5 times the noise rms, and a wrong template all fail at K = 100 (`--index-prior-rms`,
  `--noise-scale` in the script); a corruption that leaves the scanned conditional unchanged is
  not detected, so scan every field mean and, where a claim rests on a hyperparameter, that
  coordinate too.
- Marginal, not profile, for nuisance multipliers in exact-target scans (F31): per-region
  amplitude multipliers are integrated analytically (closed form in mu_r); profiled scans are
  biased (Neyman-Scott) and global-multiplier scans are a different model. The roughness
  check (`e026_roughness_check.py`) is kept as a diagnostic only: its block-shift optimum is
  linear in the coefficient scale and measures the pull of the data-pinned regions, not a
  roughness effect (E026 review).

## Rules added after the 2026-11-09 review (D028)

Each rule answers a failure that the review found in the record. They apply to every experiment
from E032 on.

- **Regime gate (F39).** Before an experiment runs, its record states which real data it stands
  for: the noise per pixel relative to the real maps of that band, and the per-pixel S/N in the
  leading band. Words such as "realistic" or "Planck-like" are not used without those numbers.
  A simulation with no real counterpart is labelled as a stress test.
- **Trivial baselines (F40, referee objection 7).** Every index-field result is reported next to
  a constant index at the truth mean and next to the truth's own region means, for the index
  error and for the CMB. CMB comparisons use the bias term of the E033 split at high latitude
  and its low-multipole band power, never the total residual.
- **Four chains (F42).** A production configuration is run with both in-flight chains written
  (`output.chains.write: [1, 2]`) and at least two runs with different seeds and dispersed
  starts. Split R-hat over the eight half-chains must be below 1.05 for every quantity a claim
  rests on, and the effective sample size is reported per chain. The D025 check stays necessary
  and is not sufficient. A single-chain result is graded B at best.
- **Replication.** A number enters a paper as a property of a method only after it has been
  reproduced on at least three noise seeds, five for a coverage statement. Otherwise it is
  labelled as one realization.
- **Calibration at the configuration reported.** Simulation-based calibration is run for the
  Commander4 implementation at the configuration a paper reports, at reduced resolution if
  needed, not only for the prototype.
- **Reading exact-target scans (E030 chain 3).** A set of scans is summarized by the optimum
  pooled over the regions a hyperparameter pools, with its error. "Flat" is not a reading for
  scans whose gains add. A pre-registration built on a scan states the pooled number.
- **Pre-registration.** Predictions state the number that would falsify them and the trivial
  baseline they must beat. The time is the record file's modification time, and it must precede
  the first output file.
- **Real maps (F41).** No real-map number is quoted in units of posterior sd unless the
  high-latitude reduced chi-square of the fit, in the noise model the chain used, is below 2.
  The rms that reached the likelihood is printed by the reader, not reconstructed.
