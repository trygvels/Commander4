# Findings so far

Established on 2026-10-25 in a single exploratory session. Experiment IDs refer to
`experiments/registry.md`. Numbers are from simulations unless stated.

## 1. Sampler correctness

**F1. Commander4's original spectral-index move is not a valid MCMC step** (E002, E005, E006).
It proposes beta', re-solves the amplitudes as a constrained realization, and accepts on chi^2 at
the drawn amplitudes. For "propose beta', then draw a' ~ p(a | beta', d)" the acceptance that
leaves the joint posterior invariant is p(beta' | d) / p(beta | d), which does not depend on the
drawn amplitudes. Accepting on chi^2 at the draws adds a random term (X' - X), X ~ chi^2 over the
amplitude degrees of freedom, and drops the -1/2 ln det A term.

| Setup | Posterior width vs exact | Mean shift |
|---|---|---|
| Per-pixel indices, 300 pixels (E002) | 1.24x | small |
| Full-sky scalar index, 256 pixels (E002) | 3.9x | -0.43 sd (+-0.09) |
| Full-sky scalar index, 1024 pixels (E002) | 5.7x | -0.99 sd (+-0.20) |
| Commander4 unit test, 150 / 600 pixels (E005) | 3.2x / 4.4x | +0.50 / +0.83 sd |

The inflation grows roughly as k^(1/4), with k the number of amplitude degrees of freedom, so at
full resolution (k ~ 10^7) it would be far larger. The shifted means match the profile ("ridge")
posterior. Commander4's `perpix_betasamp` example config uses exactly this setup at nside 512. An
independent subagent reading the upstream code reached the same conclusion. Real-data impact is
not yet measured (see OPEN_QUESTIONS Q1).

**F2. Fixed** on branch `fix/mcmc-marginal-acceptance` (E005, E006). Options `marginal`, `ridge`,
`conditional` are exact for their stated targets; `resolve_chisq` keeps the old behaviour with a
warning. Unit tests reproduce the exact posterior; a full Commander4 run on 7 MPI ranks matched the
exact ridge target within 0.03 sd with width ratios 0.94 and 0.97.

## 2. Priors matter as much as the sampler

**F3. The flat-amplitude-prior marginal can be scientifically useless for full-sky indices** (E006).
Its volume term, -1/2 ln det A summed over every pixel, acts as an implicit prior on beta. With dust
normalized at 545 GHz (above all bands), the log density rose by ~1800 towards the beta_d bound and
the exact posterior sat on the hard limits (beta_s = -2.008 at the bound, input -3.0), with a second
local mode near beta_s = -3.5 separated by a barrier of ~22 in log density. Chi^2 at the best-fit
amplitudes was slightly better at the true beta.

**F4. Ridge is the principled default.** Ridge is the profile ("spectral") likelihood of parametric
component separation (Stompor et al. 2009). As an MH target it equals the flat-prior marginal times
det(A)^(1/2), i.e. the marginal under the Jeffreys prior of the amplitude block conditional on beta,
p(a | beta) proportional to det A(beta)^(1/2), with a flat prior on beta. (Corrected 2026-10-25 after
the E012 review: it is NOT the full joint Jeffreys prior sqrt(det I(a, beta)), which carries an
a-dependent factor and gives a different beta posterior; quadrature check in
`experiments/E009_c3_nonlin_review/check_a_e.py`.) Eriksen et al. 2008 and Planck 2015 X used
Jeffreys priors to remove prior-volume bias of more than 4 sigma. `auto` now resolves to ridge where available. Conditional Gibbs (Commander3
'chisq' default) targets the flat-prior posterior and therefore carries the same volume pull, only
hidden by slow mixing. This needs a proper coverage study (Q2).

## 3. Efficiency

**F5. The Gibbs alternation between amplitudes and indices is the mixing bottleneck** (E001).
In a per-pixel toy (300 pixels, CMB + synchrotron + dust, 8 bands):

| Sampler | Amplitude solves per effective sample (noise x0.3 / x1 / x3) |
|---|---|
| Commander-style Gibbs + MH (chisq) | 42 / 14 / 3.2 |
| Same, 10x more MH steps | 41 / 14.5 / 3.2 |
| Commander3 'marginal' lnL | 20 / 7.8 / 3.4 |
| Same Gibbs loop with NUTS instead of MH | 41 / 13 / 3.1 |
| Collapsed (amplitudes integrated out) + NUTS | 1.9 / 1.5 / 1.3 |

More MH steps, or gradients inside the same Gibbs loop, do not help. Collapsing does, by 9x at
nominal noise and 22x at high S/N (20 to 28x per second in bright plane pixels). The gain grows
with S/N. The toy has independent pixels; with beams the collapse needs either low-resolution
likelihoods or log-determinants of the CG system (Q5).

**F6. Language is not the bottleneck** (E004). Full-sky per-pixel marginal likelihood at nside 512,
8 bands, 2 cores: Julia 0.070 s, Numba 0.060 s, C++ comparable, plain NumPy 1.8 s. Commander4's
Python plus C++ (nanobind) design is adequate; algorithmic gains dominate.

## 4. Commander4 implementation issues found

**F7. Per-pixel amplitudes are stored as alms with lmax = 3 nside - 1** (E007). For unsmoothed maps
the round trip alone raises reduced chi^2 from 0.57 to ~794 (matching Commander4's own logged value);
raising lmax makes it unstable. With common-resolution smoothing the mean is fine, but the white
per-pixel fluctuation draws still lose variance above lmax. Not changed yet.

**F8. The per-pixel model treats smoothed noise as white** (E006). After smoothing to a common beam
the noise is correlated; the per-pixel likelihood ignores that, and in the validation simulation
the exact ridge target sat 1.9 sd (beta_s) and 3.7 sd (beta_d) from the input. This is a model
issue, not a sampler issue (Q7).

**F9. The two-chain pipelining in `cli.py` is valid**: each chain's TOD step uses its own chain's
latest sky.

**F10. Commander3: the BeyondPlanck-era index samplers do not have the F1 bug, but `marginal` has a
sign error, and the 2024 MH module has F1** (E009). The alm / pixel-region sampler and the local
`chisq` sampler hold amplitudes fixed during MH: valid conditional Gibbs. The local `marginal` lnL adds
+1/2 ln A instead of -1/2 ln A (both `comp_lnL_*_diagonal` functions, still at HEAD), so it samples
ridge x A^(1/2) rather than the flat-prior marginal; BeyondPlanck final used it for AME nu_p (confirmed in `BP_param_c0001_v2.txt`). Cosmoglobe DR1 did not: its parameter file models AME with the exponential parameterization and `prior` index steps (E013 scoping). In a
toy the sign is irrelevant at high S/N and flips a 2 to 27 sd volume bias at low S/N. The 2024
`sample_specind_mh` family redraws amplitudes by CG and accepts on chi^2, like Commander4.
The BeyondPlanck polarized-index (alm / pixel-region) and bandpass samplers are valid conditional
steps in the BP10 configuration (E010); their main approximations are prior-drawn frozen regions
and a diagonal spurious-map chi^2 for relative bandpass shifts.
Independent review (two stats-reviewer agents per item, 2026-10-25, reports in
`experiments/E009_c3_nonlin_review/`) confirmed all of the above and added: BP10 sampled nu_p as
ONE full-sky value with the band monopoles profiled inside each proposal; the BP10 Jeffreys prior
was hard-disabled in the code (`comm_param_mod.f90:933`), so the dust beta_pol target was
likelihood times N(1.56, 0.1^2) (E009's earlier statement is corrected); the alm / region sampler's
low-resolution chi^2 weights nside-16 pixels by N_sub mean(1/sigma^2) instead of N_sub/mean(sigma^2)
(1.56x for a 2x variance contrast); the 2024 MH module also ignores all priors and bounds, perturbs
only pixel region 1, and `sample_template_mh` uses a multiplicative proposal without the Hastings
term. **Impact of the sign error on BeyondPlanck's nu_p is negligible (E013):** the AME pixels in
the sampling mask have median S/N above 30 and nu_p was a single parameter; the flipped volume term
moved the full-sky conditional posterior by 0.8 of its 0.01 GHz sd (0.008 GHz), the correct
marginal would move it equally the other way; per pixel the median shift is 0.04 ridge sd, 0.14 sd
at S/N 3 to 10. Cosmoglobe DR1 did not use the `marginal` path at all.

**F11. BeyondPlanck's absolute bandpass shift was accepted on rank 0's chi^2 only** (E015, code
reading confirmed line by line). In the absolute-shift branch `compute_chisq_abs_bp` accumulates the
TOD chi^2 per MPI rank and `finalize_binned_map` is called without `chisq_S`, so the only reduce never
runs; `sample_bp` sums on rank 0. The 30 GHz absolute shift posterior is the conditional given rank
0's scans (wider by about sqrt(N_ranks), centred on those scans). Present at BP10 and HEAD, and in
the WMAP, DIRBE, SPIDER, LiteBIRD and QUIET TOD modules. Also: no Gaussian prior enters the
bandpass acceptance (flat prior; the `bp_init` "PROP MEAN" is a proposal width), and after an
accepted shift the frequency map handed to component separation still carries the old correction
while the mixing matrices use the new one.

## 5. Commander4 versus Commander3 coverage

Commander4 covers Planck LFI/HFI TOD processing (gain, correlated noise, PSD, orbital dipole,
sidelobes, jumps, glitches, HFI demodulation), CG and per-pixel component separation with
CMB/synchrotron/dust/free-free/spinning dust/radio sources/templates, SO and LiteBIRD simulations,
chain restarts and a large test suite. Missing relative to Commander3: WMAP, DIRBE, SPIDER, QUIJOTE
TOD; zodiacal light (`dev-zodi` branch unmerged); per-band monopoles/dipoles; bandpass-shift
sampling; ADC corrections; crosstalk (`dev-cross-talk-sim`); C(l) sampling; spatially varying
spectral parameters; several components (CO lines, physical dust, tabulated MBB, curved/broken
power laws, lognormal AME, PAH, dust extinction). Commander3 is still needed for full Cosmoglobe-style
runs until these exist and a published Commander3 result is reproduced.

## 6. From the literature (2026-10-25 survey; verify before citing)

- No published method infers the beta power spectrum or correlation length as a free hyperparameter
  inside a CMB Gibbs sampler, and no SPDE/GMRF prior has been used on HEALPix for foreground
  parameters. Grumitt et al. 2020 (NUTS, hierarchical per-region beta) is the closest; it hit
  funnel divergences. This is an open gap.
- Commander3's alm-based beta sampler clipped values "in case of ringing" and fixed its prior
  hyperparameters; ringing was a real problem.
- Field-level CMB inference at nside 512 (Flinch, 2510.26691) found microcanonical samplers
  (MCLMC/MAMS) 80 to 100x more efficient per gradient than NUTS.
- Full Gibbs samples cost 170 CPU-h (BeyondPlanck), 812 CPU-h (Cosmoglobe DR1 WMAP+LFI), ~3000
  CPU-h (LiteBIRD); Cosmoglobe's own cosmological-parameter sampler needs ~2000 CPU-h per
  independent sample and states a 10x reduction is needed.
- The tau debate (low-l polarization ~0.054 to 0.058 vs DESI+CMB 0.090 +- 0.012) is the clearest
  high-impact opening; it needs HFI end to end (Cosmoglobe is building "OpenHFI").
- Several groups find beta_s steepening and spatial variation at scales finer than BeyondPlanck's
  four regions (Weiland 2022, QUIJOTE, CLASS rejects uniform beta_s at PTE < 0.001).

## 7. Process

- Commander4 builds on Linux with pip-installed MPICH; ~860 tests pass on all three branches.
- The Commander4 lead developer openly uses Claude and Copilot for helper code; AI-assisted
  contributions are normal in this project.

## F12 (preliminary): at low per-pixel S/N the ridge likelihood inflates a hierarchical index field and the flat-prior marginal collapses it

E012 probe and profile scan (one replicate, nside_d 16, 48 index regions, S/N 0.4 to 16): with
amplitudes profiled (ridge, w0) the per-region likelihood in beta is rough and multimodal at S/N
below a few, the per-region maxima scatter with rms 2.2 about the truth, the inferred field sigma is
about 1.2 against a truth of 0.065 (z = +22 to +25) and NUTS mixes at ESS 4 to 5 per 1500 s in every
parameterization. With the flat amplitude prior (w1) every region collapses to beta about -0.1 with
conditional sd 0.01 (field version of F3). With a proper Gaussian amplitude prior the conditionals are
unimodal and prior-dominated at low S/N and correct at S/N 16. Preliminary until the E012 coverage
study (queued) gives rates over replicates. Bearing on Q2: for spatially varying indices the amplitude
prior is not a nuisance choice; neither the profile nor the flat-prior marginal is safe at low S/N.
