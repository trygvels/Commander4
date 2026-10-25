# E010: Review of Commander3's polarized spectral-index and bandpass samplers (BeyondPlanck)

- **Date:** 2026-10-25
- **Code reviewed:** Cosmoglobe/Commander tag BP10 (3aa493e) and HEAD (f760333):
  `comm_nonlin_mod.f90` (`sample_nonlin_params`, `sample_specind_alm`), `commander.f90`
  (`process_TOD`), `comm_tod_bandpass_mod.f90` (`sample_bp`), `comm_tod_lfi_smod.f90`,
  `comm_tod_mapmaking_mod.f90` (`finalize_binned_map`), `comm_tod_driver_mod.f90`
  (`compute_chisq_abs_bp`); `param_BP10_final_v1.txt`, `LFI_tod.defaults`, DR1v2 defaults.
- **Validation level:** code review (L0). Independent `stats-reviewer` pass still to do.

## Polarized spectral indices (BeyondPlanck: synchrotron beta_s in 4 regions, dust beta_d full sky)
Settings: `COMP_BETA_POL_LMAX` 100 (synch) / 0 (dust), `ALMSAMP_PIXREG = .true.`,
`ALMSAMP_PRIORSAMP_FROZEN_REGIONS = .true.`, synch regions 1 and 3 frozen with prior N(-3.3, 0.1),
`ALMSAMP_OPTIMIZE_ALM = .false.`, `ALMSAMP_BURN_IN = 2`, 100 MH steps per Gibbs iteration, Jeffreys
prior on for dust.

Target of the MH chain: p(theta_regions | a, d, other params) with every component amplitude held
fixed, full-resolution chi^2 over the index mask, Gaussian prior per region. Region values map
deterministically to the index field (10 deg smoothing, clipping to the hard bounds, band limit
lmax 100). Proposals: correlated Gaussian with a fixed Cholesky factor from an earlier run; step
length adapted only in the first 2 Gibbs iterations. **Valid conditional Gibbs step. No F1 bug.**

Caveats, all approximations rather than errors:
- Frozen regions are drawn from the prior each iteration and the synchrotron amplitude is then
  redrawn by CG. Those regions report prior draws, not posterior samples.
- The index likelihood uses only pixels inside `COMP_INDMASK` (a cut likelihood; the amplitudes see
  all pixels).
- Proposals outside the hard bounds are redrawn without a normalization correction; negligible with
  bounds [-4.5, -1.5] and prior width 0.1.
- The dust Jeffreys prior uses the Fisher information at fixed amplitudes, which depends on the
  amplitudes; the amplitude step ignores that dependence. Negligible.

## Bandpass shifts (BeyondPlanck LFI)
Proposals are drawn in `process_TOD` around the current state: relative shifts as
x + P L eta with P projecting onto zero mean across detectors (symmetric random walk on the
constraint surface), absolute shifts (30 GHz only, odd iterations) as a common N(0, s) offset. The sky
signal for each proposal uses the current component amplitudes; n_corr and gains are drawn once with
the current bandpass and shared by current and proposed chi^2.
- **Absolute shift:** TOD-level white-noise chi^2 given n_corr, gains, sidelobes and orbital dipole.
  Exact conditional likelihood, flat prior. Valid.
- **Relative shifts:** chi^2 of the per-pixel ML spurious maps, sum_j S_j^2 / Var(S_j), with I, Q, U
  solved freely per pixel. A deterministic function of the proposal, so MH is valid for it. It is an
  approximation of the I/Q/U-profiled likelihood: it uses only the diagonal of the S covariance,
  ignoring correlations between the detectors' S maps, which can misstate the width of the
  relative-shift posterior. Worth quantifying.
- **`sample_bp` structure:** all proposals are centred on the starting value, but each is compared
  with the currently accepted one. With `TOD_NUM_BP_PROPOSALS_PER_ITER = 1` (BP10, all DR1v2
  defaults) this is a single symmetric random-walk step and valid. With more than one proposal it
  would be an independence-type sampler without the proposal-density correction and would not
  target the posterior. The parameter comment already says "1 for sampling". Unchanged at HEAD.

## Conclusion
Both methods are valid MCMC steps for their stated conditional targets in the BeyondPlanck
configuration. Open items: an independent re-derivation (stats-reviewer), and a check of how much
the diagonal spurious-map chi^2 misstates the relative bandpass uncertainty.
