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

## L2: simulation-based calibration
Draw parameters from the prior, simulate data from the model, sample, record the rank of the truth
among posterior draws. Ranks must be uniform (Talts et al. 2018). Use the joint log likelihood as an
extra test quantity to catch posteriors that equal the prior (Modrak et al. 2022). For
high-dimensional fields use TARP coverage (Lemos et al. 2023). Low resolution (nside 8 to 32), a few
hundred replicates. Required for every prior choice in Q2/Q3.

## L3: end-to-end Commander4 runs on model-consistent simulations
Full Commander4 (MPI, real I/O and chain writing) on simulations drawn from the model being sampled,
compared to the exact or independently computed target. Harness: `aux/validate_mcmc_likelihood.py`.
Mind F7/F8: per-pixel runs need common-resolution smoothing, and smoothed noise is not white.

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
