# E009: Does Commander3's nonlinear sampler share the Commander4 bug?

- **Date:** 2026-10-25
- **Code reviewed:** Cosmoglobe/Commander tags BP10 (3aa493e, 2021-11-19) and HEAD (f760333, 2026-03-17):
  `comm_nonlin_mod.f90` / `comm_nonlin_smod.f90`, `comm_mh_specind_mod.f90`, `comm_chisq_mod.f90`,
  `commander.f90`; BeyondPlanck final parameter file `param_BP10_final_v1.txt` and its defaults.
- **Validation level:** code review plus L1 toy (`E009_c3_nonlin_review/ame_nup_sign.py`).

## Question
Does any Commander3 spectral-parameter sampler accept on chi^2 after redrawing amplitudes (F1), and are
the per-pixel likelihood options exact for their stated targets?

## Results

1. **Alm / pixel-region sampler (`sample_specind_alm`)**, used in BeyondPlanck for polarized beta_s
   (4 regions, regions 1 and 3 frozen) and polarized beta_d (one region, Jeffreys prior): amplitudes are
   held fixed for the whole MH chain and the acceptance uses full-resolution chi^2 plus the region
   prior. This is a valid conditional Gibbs step. **No F1 bug.**
2. **Local sampler, `chisq` lnL**: residual evaluated with the smoothed amplitude held fixed. Valid
   conditional step on the low-resolution approximate likelihood. **No F1 bug.**
3. **Local sampler, `ridge` and `marginal` lnL**: followed by a CG redraw of the component's amplitude
   (`samp_cg` in BP10 `sample_nonlin_params`), which is the right partially collapsed structure.
   **But the `marginal` log-determinant has the wrong sign**: both `comp_lnL_marginal_diagonal` and
   `comp_lnL_max_chisq_diagonal` add `-0.5*log(invMNM)` = +1/2 ln A, while the flat-prior marginal is
   -chi2_min/2 - 1/2 ln A (checked by quadrature). The code therefore samples ridge x A^(+1/2) instead
   of ridge x A^(-1/2): the prior-volume pull has the right size and the opposite direction. Still
   present at HEAD. `ridge` and `chisq` are unaffected.
   Used in: BeyondPlanck final and DR1v2 defaults for the AME peak frequency nu_p (intensity,
   `ame_LFI_localsampler.defaults`); BP v7/v8/8.1 parameter files also used it for intensity beta_s
   and beta_d (BP10 final switched those to `prior`).
4. **Newer MH module (`sample_specind_mh`, `sample_template_mh`, `sample_mbbtab_mh`,
   `sample_gain_firas`, added 2024-02-08)**: propose, run `sample_all_amps_by_CG` (which adds
   fluctuations when OPERATION = sample), accept on chi^2. **Same structure as F1** whenever
   `MCMC_SAMPLING_GROUP_UPDATE_CG_GROUPS` is not `none`. `sample_specind_mh` also adds no prior term
   and no bounds check. Commander4's `mcmc.py` is evidently a port of this. No parameter file in the
   repository uses these groups, so which production runs (likely DR2-era) used them with CG updates
   is unknown.

## Toy (AME-like nu_p, 1344 nside-16 pixels, 8 bands 23 to 94 GHz, prior N(21, 3), truth 22 GHz)

| Noise per band | ridge bias | correct marginal bias | C3 code 'marginal' bias |
|---|---|---|---|
| 2 uK | -0.21 sd | -0.43 sd | +0.01 sd |
| 20 uK | -0.05 sd | -2.25 sd | +2.17 sd |
| 100 uK | +0.14 sd | -18.1 sd | +27.1 sd (at the 30 GHz grid edge) |

Means over 20 noise realizations (MC error on each bias about 0.2 sd). At high S/N the sign barely
matters; once many low-S/N pixels enter the sum, the volume term dominates and the sign decides which
bound the chain drifts to. The impact on the published BeyondPlanck / DR1 nu_p depends on the S/N of
the pixels in their sampling mask and is not measured here.

## Smaller issues found (approximations, not F1-type errors)
- Frozen regions and `prior` lnL draw indices from the prior each iteration (BP: polarized beta_s in
  regions 1 and 3, intensity beta_s and beta_d). Deliberate, but those are prior draws, not posterior
  samples, and the amplitudes are conditioned on them.
- Region proposals are redrawn until inside the hard bounds without the proposal-normalization
  correction; negligible unless the posterior sits near a bound.
- Low-resolution likelihood treats smoothed noise as white (as F8).
- Combined monopole sampling sets monopoles to their conditional mean given the old amplitude, then
  marginalizes or profiles the amplitude: a deterministic heuristic, not an exact conditional.
- Alm mode without pixel regions applies the Gaussian prior to a_00 with the index-space rms, which is
  sqrt(4 pi) = 3.5x tighter on the mean index than stated. BeyondPlanck used pixel-region mode, so it
  is not affected.
- Jeffreys prior for beta_d uses the Fisher information at fixed amplitudes (depends on a); the
  amplitude step ignores that dependence. Negligible.

## Conclusion
The BeyondPlanck-era index samplers do not have the Commander4 acceptance bug. The `marginal`
likelihood option has a sign error in its volume term, which affected the AME nu_p step in BeyondPlanck
(and likely DR1). The 2024 MH module has the Commander4 bug.

## Independent review outcome (2026-10-25, two stats-reviewer agents per item, given code paths and stated targets only)
Full reports: `E009_c3_nonlin_review/independent_reviews_2026-10-25.md` (also covers E010 and E012).
Confirmed: (1) alm / pixel-region sampler is a valid MH step for its target with amplitudes fixed,
no F1 structure (both reviewers; 2e6-step toy check of the bound handling). (2) `marginal` lnL
sign flip at `comm_nonlin_mod.f90:3513` (and the unused `comp_lnL_marginal_diagonal`), verified by
quadrature by both reviewers independently. (3) 2024 MH module has the F1 structure; toy
inflation 1.3x (k = 4) to 4.6x (k = 1024) amplitude degrees of freedom, acceptance pinned near 0.5
independent of the step size; the no-refresh variant converges to the profile likelihood.
New facts the review added (not in the text above):
- BP10 sampled nu_p as ONE full-sky pixel region (sum of the lnL over all unmasked nside-16
  pixels), not pixel by pixel. The released BP nu_p map is constant (25.60 GHz). E013 uses this.
- With `COMBINED_MONOPOLE_SAMPLING` the band monopoles are refit deterministically (prior-shrunk,
  noise-unweighted least squares) inside every nu_p proposal: monopoles are profiled, not sampled.
- The low-resolution chi-square used by the alm / region sampler weights each nside-16 pixel by
  N_sub * mean(1/sigma^2) instead of N_sub / mean(sigma^2) (Jensen; 1.56x for a 2x variance
  contrast, 2.7x for lognormal(0.5) scatter), and masked sub-pixels enter as zeros with full weight.
- The Jeffreys prior was hard-disabled in BP10: `comm_param_mod.f90:933` sets
  `cs_apply_jeffreys = .false.` after reading the parameter. E009's statement that BP10 used a
  Jeffreys prior for polarized beta_d is therefore WRONG; the dust beta_pol target was likelihood
  times N(1.56, 0.1^2). Even when enabled, `compute_jeffreys_prior` returns the Fisher information
  of one global index shift, not sqrt(det F) over the region parameters.
- The 2024 MH module also: ignores every prior and bound on theta (flat improper prior); perturbs
  only pixel region 1; `sample_mbbtab_mh` recomputes chi2_old from a fresh CG draw every step;
  `sample_template_mh` uses a multiplicative proposal without the Hastings correction (biased toward
  zero, chain drifts); the chi-square mask and band list can differ from the CG ones, so even the
  `UPDATE_CG_GROUPS = none` Gibbs chain has no common joint in that case.
- Hard-bound handling by redraw-until-inside makes the stationary density pi(theta) Z(theta);
  negligible at BP10 settings (Z = 1 to double precision).
