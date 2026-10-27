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
| Real maps, nside 64, 38165 masked pixels, 17 bands, 5 components (E014) | 16.5x (beta_d) / 17.1x (beta_s) | +3.0 +- 0.5 sd (beta_d), -2.7 +- 2.0 sd (beta_s) |

The inflation grows roughly as k^(1/4), with k the number of amplitude degrees of freedom, so at
full resolution (k ~ 10^7) it would be far larger. The shifted means match the profile ("ridge")
posterior. Commander4's `perpix_betasamp` example config uses exactly this setup at nside 512. An
independent subagent reading the upstream code reached the same conclusion. Real-data impact
(E014, Q1): on public maps at nside 64 the legacy posterior is 16 to 17x too wide against the exact
ridge target computed on the same data, consistent with the k^(1/4) scaling (predicted 18.5x), and
the dust index is shifted by 3 ridge sd. At nside 512 to 1024 the scaling predicts factors above 50.

**F2. Fixed** on branch `fix/mcmc-marginal-acceptance` (E005, E006). Options `marginal`, `ridge`,
`conditional` are exact for their stated targets; `resolve_chisq` keeps the old behaviour with a
warning. Unit tests reproduce the exact posterior; a full Commander4 run on 7 MPI ranks matched the
exact ridge target within 0.03 sd with width ratios 0.94 and 0.97. On real maps at nside 64 (E014)
the ridge chain matches an independent grid evaluation of its target within 0.2 sd (width ratios
1.01 and 1.18 at ESS 97 and 44; with tuned proposals and ESS 559 and 422: shifts -0.06 +- 0.04
and +0.06 +- 0.05 sd, width ratios 0.96 and 0.94), and the marginal chain matches the exact flat-prior marginal
target within 0.05 sd (width ratios 1.08 and 1.23).

## 2. Priors matter as much as the sampler

**F3. The flat-amplitude-prior marginal can be scientifically useless for full-sky indices** (E006).
Its volume term, -1/2 ln det A summed over every pixel, acts as an implicit prior on beta. With dust
normalized at 545 GHz (above all bands), the log density rose by ~1800 towards the beta_d bound and
the exact posterior sat on the hard limits (beta_s = -2.008 at the bound, input -3.0), with a second
local mode near beta_s = -3.5 separated by a barrier of ~22 in log density. Chi^2 at the best-fit
amplitudes was slightly better at the true beta.
On real nside-64 maps (E014, dust referenced at 70 GHz, 17 bands, 38165 pixels) the volume term
shifts the exact beta_s target by +0.019 (54 ridge sd) and beta_d by +0.0002 (5 sd) with the same
widths; two exact chains (`ridge`, `marginal`) reproduce both targets within 0.2 sd, and the
reduced chi-square differs by 0.07 percent, so the data do not pick the prior.

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
at S/N 3 to 10. Cosmoglobe DR1 did not use the `marginal` path at all. Branch status (2026-10-26 survey, E009): present at the master tip
(2026-03-17) and on all active branches; corrected only on an unmerged test branch
(`fable-test-SOLAT`); the 2024 MH module's acceptance is still wrong on `hfi_Tbol` after its
2026-10-23 "Fixed MH sampler" rewrite.

**F11. BeyondPlanck's absolute bandpass shift was accepted on rank 0's chi^2 only** (E015, code
reading confirmed line by line). In the absolute-shift branch `compute_chisq_abs_bp` accumulates the
TOD chi^2 per MPI rank and `finalize_binned_map` is called without `chisq_S`, so the only reduce never
runs; `sample_bp` sums on rank 0. The 30 GHz absolute shift posterior is the conditional given rank
0's scans (wider by about sqrt(N_ranks), centred on those scans). Present at BP10 and HEAD, and in
the WMAP, DIRBE, SPIDER, LiteBIRD and QUIET TOD modules. Also: no Gaussian prior enters the
bandpass acceptance (flat prior; the `bp_init` "PROP MEAN" is a proposal width), and after an
accepted shift the frequency map handed to component separation still carries the old correction
while the mixing matrices use the new one. Present on all 37 Commander3 branches with source (2026-10-26
survey, E015).

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

**F12. For a spatially varying index the flat-prior marginal is unusable at every S/N, the profile
likelihood (ridge) is calibrated only above a per-index-pixel S/N of about 15, and a proper Gaussian
amplitude prior gives calibrated indices only if it is well specified** (E012 coverage study
`study_v1`: 12 replicates x noise 0.3 / 1 / 3 x three modes, nside_d 16, 48 index pixels, template
amplitudes with 0.1 dex scatter, 108 NUTS cells). Flat-prior marginal (w1): every index pixel pulled
to beta about -0.2 with sd 0.01, standardized bias z of +100 to +300, zero 68 / 95 percent coverage,
field mean off by 400 to 1400 sd (field version of F3, now quantified). Ridge (w0): at S/N >= 15.6, z
+0.09 +- 0.06 and +0.17 +- 0.11 with coverage 0.70 / 0.96 and 0.67 / 0.91; at S/N below 3, z +5 to
+13 with sd(z) 20 to 35, coverage 0.3 to 0.6, field sd biased +5 to +28 sd, ESS 4 to 6 (the profile
maxima of low-S/N regions scatter by about 2 in beta and the hierarchy reads that as field variance).
Gaussian amplitude prior N(m, T) with T^(1/2) = max(0.5 |m|, 20) on template data with 0.23 |m|
scatter: z +0.6 (coverage 0.59 / 0.93) at noise 3, +2.1 at noise 1, +14 at noise 0.3 (median error
1.5 in beta). With the generative prior in the likelihood the same cells are calibrated (z +0.08 +-
0.06 and -0.22 +- 0.06, coverage 0.69 / 0.95 and 0.68 / 0.92), so the bias is prior misspecification,
not a sampler error; widening the prior to 10 abs(m) makes it 4x worse (z +60 and +26), because the
Gaussian marginal tends to the w1 volume term as T grows. There is no uninformative safe Gaussian
amplitude prior: it must be learned from the data or replaced by the profile at adequate S/N. Bearing on
Q2/Q3: the index resolution must be chosen from the S/N map so index pixels stay above the ridge
threshold, or a proper amplitude prior must be used and its specification treated as a modelled
systematic.

**F13. The Cosmoglobe DR1 dust map is referenced at 70 GHz, not the 545 GHz its FITS header states**
(E014; DATA.md section "Cosmoglobe DR1 thermal dust reference frequency"): the parameter file
overrides `COMP_NU_REF_T03` with `70 857` and Commander3 keeps the first token. Verified by
regressing the 857 and 353 GHz maps on the DR1 amplitude at |b| > 30 deg (slopes 12.2 and 7.2
against 11.3 and 7.3 predicted for 70 GHz, correlation 0.99). Consequence for us: DR1 dust can be
used as a template or simulation truth only with the 70 GHz reference; the BeyondPlanck dust file
has the same header error. The residual 7 to 14 percent is the delta-bandpass versus
bandpass-integrated SED and is not yet quantified.

**F14. On real maps the Gibbs conditional for a full-sky spectral index cannot mix** (E014, nside
64, 38165 masked pixels, 17 bands): the conditional sd of beta_s given the per-pixel amplitudes is
2.0e-6 against a marginal sd of 3.6e-4 (ratio 180; for beta_d 2.8e-5 against 4.3e-5, ratio 1.5),
so Gibbs with an MH step on chi2 at fixed amplitudes needs of order 3e4 iterations per independent
beta_s sample. The Commander4 `conditional` chain with proposals of 0.002 stayed at its start value
for 2000 iterations. The collapsed (`ridge`, `marginal`) chains mix in tens of iterations on the
same data. This is F5 measured on real data and is the quantitative case for collapsing amplitudes
in any index sampler at Planck resolution.

**F15. The joint posterior with flat amplitude priors has the collapsed flat-prior marginal (w1, with
-1/2 ln det A) as its exact spectral-parameter marginal; a Gibbs sampler that draws amplitudes with
their fluctuation and moves beta on chi2 at fixed amplitudes targets exactly that** (E012 corollary,
independently derived and checked numerically by the stats-reviewer, 2026-10-26). Consequences: the
F3 and F12 pathologies are properties of the flat-prior joint, not of the collapsed implementation;
the profile ("ridge") and the fixed-amplitude conditional are not marginals of that joint; and a
Gibbs chain with that target either shows the collapse or is not mixing. This does not transfer to
Commander3/4 production chains as they stand, whose targets include a C_l prior on the CMB
amplitudes (which regularises the CMB/synchrotron direction but not synchrotron/dust), index priors
and bounds, and index steps on smoothed data; what it does say is that "sample amplitudes with flat
priors, then sample the index on chi2" is not a route around the prior question of Q2.

**F16 (provisional until the E016 SBC passes). A hierarchical amplitude prior, with the per-component
scale and relative scatter learned from the data, gives calibrated index-field posteriors on
template data at every index-pixel S/N tested, including below S/N 1 where the profile is biased by
5 to 13 sd and the fixed Gaussian prior by 0.6 to 14 sd** (E016 coverage study, same 108 replicates
and S/N bins as E012 study_v1). Over all index pixels, mean z is +0.10 +- 0.04, +0.35 +- 0.04 and
+0.24 +- 0.05 at noise 0.3 / 1 / 3 with coverage 0.66 / 0.95, 0.64 / 0.95 and 0.67 / 0.93; in the
lowest populated S/N bins (S/N 1.6 to 3.3 at noise 0.3, below 0.54 at noise 1 and 3) z is +0.03,
-0.02 and +0.35 with coverage 0.70 / 0.96, 0.56 / 0.96 and 0.67 / 0.91. The result survives a wrong
shape template: a 30 percent large-scale distortion or 23 percent independent per-pixel scatter of
the template gives z within 0.3 sd and coverage 0.64 to 0.75 / 0.92 to 0.98 everywhere, because the
learned scatter absorbs the shape error (log s recovered to its posterior width in both cases, e.g.
-1.10 +- 0.04 against the expected -1.11 for the scattered template at noise 0.3). Standard errors
above are pooled over pixels; with the replicate-clustered SE (pixels of a replicate share truth
and noise) the all-pixel biases are +0.10 +- 0.07, +0.35 +- 0.16 and +0.24 +- 0.25, so the apparent
+0.3 sd bias at noise 1, common to all variants because they share the 12 truths, is a 2 sigma
fluctuation: 24 fresh replicates at noise 1 give +0.01 +- 0.07, and all 36 give +0.10 +- 0.07
with coverage 0.71 / 0.975. One 58-pixel bin (S/N 3.3 to 15.6 at noise 0.3) has
coverage 0.48 / 0.86, the same as the profile there. This is the
first mode in this programme usable below index-pixel S/N 15 without a hand-set amplitude prior;
it answers the Q2 design question in favour of learning the amplitude prior rather than fixing it or
removing it (F12).

**F17. On real low-resolution maps the full-sky spectral indices are set by the zero-level
(monopole) convention of the input maps, at the 0.03 to 0.4 level, hundreds of posterior sd**
(E014 passes 2 to 4, exact 2D grids of the ridge target, nside 64, |b| > 30 chisq mask, 17 or 20
bands). With the WMAP/LFI maps median-subtracted at high latitude and the HFI maps referenced to
the DR1 CMB plus dust (two conventions in one run) the ridge beta_s is -3.3424 (17 bands) or
-2.9103 (20 bands); with one convention for all 20 bands (offset = mean over |b| > 30 of map minus
the DR1 sky model, `nextgen/scripts/dr1_model_zero_levels.py`) it is -2.8755 and beta_d moves from
1.4150 to 1.3991; removing the 3.4 K non-solar dipole of the Haslam map moves beta_s by a further
+0.011. The posterior sd is 1e-4 (beta_s) and 3e-5 (beta_d), so the statistical error of a full-sky
index on these maps is irrelevant next to the convention. Commander4 has no monopole sampling; until a
per-band zero level with a prior is part of the model (Q15), spectral indices from maps-only runs on
real data are convention-dependent and are not science results. The method results of E014 (chain
against exact grid on identical inputs) are unaffected.
