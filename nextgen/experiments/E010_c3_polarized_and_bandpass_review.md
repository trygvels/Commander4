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

## Independent review outcome (2026-10-25, two stats-reviewer agents, code paths and stated targets only)
Full reports: `E009_c3_nonlin_review/independent_reviews_2026-10-25.md`.
Polarized index sampler: confirmed valid conditional MH for its low-resolution surrogate
likelihood; see the E009 addendum for the new facts (Jeffreys prior hard-disabled, Jensen-biased
low-resolution weights, frozen regions as prior draws make the joint chain a cut-type sampler).
Bandpass sampler, confirmed: single symmetric proposal per iteration is a valid MH step for a
FLAT-prior conditional (no Gaussian prior enters the acceptance, `comm_tod_bandpass_mod.f90:46-58`;
the "PROP MEAN" entry in `bp_init_*.dat` is a proposal width, not a prior); multi-proposal mode is
not invariant (exact 5-state enumeration and Gaussian MC: variance off by 10 to 40 percent); the
relative-shift statistic uses only the diagonal of the spurious-map covariance (toy: 8.6 vs exact
17.1), and is a per-pixel profile over I, Q, U with flat priors, not the conditional given the sky
model.
NEW, verified by reading the code myself after the review (E015): in the ABSOLUTE-shift branch
(30 GHz, odd iterations) `compute_chisq_abs_bp` accumulates `chisq_S` per MPI rank from that
rank's scans and `finalize_binned_map` is called without `chisq_S`, so the only reduce
(`comm_tod_mapmaking_mod.f90:565-570`) never runs; `sample_bp` then sums `chisq_S` on rank 0 only
(`comm_tod_lfi_smod.f90:796, 846-855` at BP10; `:837, 900-919` at HEAD f760333, same structure
in the WMAP, DIRBE, SPIDER, LiteBIRD and QUIET TOD modules). The 30 GHz absolute bandpass shift was
therefore sampled from the conditional given only rank 0's scans. Also noted: on acceptance the
frequency map handed to component separation still carries the old bandpass correction while the
mixing matrices switch to the new shifts (one TOD iteration out of sync).
