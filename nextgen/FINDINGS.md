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
det(A)^(1/2), which is the contribution of a joint Jeffreys prior on (a, beta) (literature agent's
derivation; Eriksen et al. 2008 and Planck 2015 X used Jeffreys priors to remove prior-volume bias of
more than 4 sigma). `auto` now resolves to ridge where available. Conditional Gibbs (Commander3
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
ridge x A^(1/2) rather than the flat-prior marginal; BeyondPlanck final used it for AME nu_p. In a
toy the sign is irrelevant at high S/N and flips a 2 to 27 sd volume bias at low S/N. The 2024
`sample_specind_mh` family redraws amplitudes by CG and accepts on chi^2, like Commander4.

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
