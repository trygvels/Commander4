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
in any index sampler at Planck resolution. Scope (added 2026-10-27): measured for one
configuration only (full-sky intensity indices, 17 bands, nside 64, high S/N); the ratio already
varies from 180 (beta_s) to 1.5 (beta_d) within it, so it does not transfer to other published
configurations without measurement. In particular nothing is established about the mixing of the
BeyondPlanck XV polarized beta_s step (four regions, LFI 30 to 70 plus WMAP Ka/Q/V and HFI 353 in
Q/U, prior N(-3.3, 0.1)); that cell is a planned E014 follow-up. BP XV itself ran many Metropolis
steps per Gibbs iteration at fixed amplitudes (lag-n correlation below 0.1) and reports a beta_s
Gibbs correlation length of 30 to 50 over about 500 samples (2011.08503 Sect. 5.3.2), so its
within-conditional mixing was handled and the amplitude-index Gibbs coupling was measured, not
ignored; the open question is only whether 30 to 50 is the full coupling length.

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

**F16 (SBC PASS on 600 replicates: the first run of 300 failed the pre-registered ECDF band on one of 16 quantities at the band edge, +2.9 se on one faint pixel's rank mean; the seeded replication of 300 passed the full battery and the one-sided rank-mean test, -0.6 se, pooled +1.7 se, so the deviation is a fluctuation by the pre-registered rule; E016). A hierarchical amplitude prior, with the per-component
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

**F18. For one full-sky spectral index on template data the flat-prior marginal is biased or
collapsed at every S/N, the profile is unbiased but up to two times too narrow below full-sky S/N
100, and the hierarchical amplitude prior is unbiased and calibrated at every S/N tested** (E017,
nside 16, 12 bands, three components, 96 profile and 24 other replicates per noise scale, exact 1D
grids for the three non-hier modes). Flat marginal: z +6.0 sd at S/N 680, +17.7 at 230, and at the
upper edge of the support in all 48 replicates at S/N 68 and 23. Fixed Gaussian with the E012
misspecified template: +5.2 / +10.1 / +4.4 / +0.4 sd. Profile: z within 0.1 +- 0.1 at all scales,
sd z 0.93 / 0.89 / 1.14 / 1.90, coverage 0.72/0.98, 0.72/0.98, 0.56/0.93, 0.33/0.66. Hier: z
+0.10 / +0.31 / -0.12 / -0.07 (+- 0.2), sd z 0.8 to 1.1, coverage 0.50 to 0.83 / 0.88 to 1.00, the
same with a 30 percent large-scale template distortion, hyperparameters recovered where detected;
hier minus profile mean within 0.3 profile sd. Eight hier chains that had collapsed at the
unguarded log s floor were rerun with the guard. The marginal-mode gaps of E014 (54 sd, 16 to 23
sd) are this bias; with hier or ridge the real-data index would agree to a fraction of a sd.

**F19. Commander4's `likelihood: hier` (index field at nside_beta plus hierarchical amplitude
priors, NUTS over MPI) samples its stated target at L1 and L3** (E018). L1: exact-target tests in
process (scalar chain against a grid; nside-1 field against an independent centred Gibbs
reference). L3: a 12-rank end-to-end run (nside 16, 12 bands, synchrotron index field at nside_beta
2 with K 47 non-centred coefficients, mu, log sigma, log rho, and hier hyperparameters for sync and
dust; 3000 Gibbs iterations) agrees with the SBC-validated E016 prototype chain on the same maps
in all 55 scalars: max 1.9 pooled SE (rms 0.75 over the 48 index pixels), sd ratios 0.86 to 1.09,
0 divergences, 15 gradients per Gibbs iteration (prototype 15 per sample), rms z against the
simulation truth 1.10; the scalar configuration likewise (max 1.8 SE). Negative control: the same
comparison puts the Commander4 hier chain 19 MC SE (0.6 sd) from the exact profile posterior on
the same data and over 1000 sd from the collapsed flat-prior grid. The L2 gate (E016 SBC of the
prototype target) is still running; hier stays non-default (D003, D011).

**F20. On the non-centred hierarchical index-field target (K 11 to 191 plus hyperparameters,
diagonal preconditioning) MAMS and NUTS have the same effective samples per gradient** (E019:
ratios 0.9 to 1.2 at the median index pixel in all five cells, 0.3 to 1.1 at the worst pixel,
MAMS warm-up 2 to 20 times more gradients; both samplers 5 to 240 times worse on the centred
funnel). The literature's 3 to 100 times gains for MAMS are reported at d of order 1e4 and above;
at Phase 2 dimensions there is no case for a second sampler, and NUTS stays the production sampler
(M2.3 closed; revisit for the joint amplitude-plus-index target of Q5).

**F21. At high S/N the collapsed hierarchical index posterior runs away under SED or template
misspecification; the sampler is not at fault** (E020, PySM d1s1c1 at nside 64, 12 intensity bands,
white noise at 0.03 x the E012 levels, `likelihood: hier`, GMRF index fields at nside_beta 16 for
beta_s and beta_d with fixed T_d = 20.93 K). The NUTS chain climbed from the flat prior-mean start
(log posterior 217,720 below the state at the PySM truth) to states 104,122 above it by iteration
182, with beta_s mean -4.05, sd 3.84 and range -26 to -0.6 (truth -3.00, sd 0.063) and beta_d
shifted by up to +0.08 at the Galactic centre (30 Fisher sd). The target and its gradient were
verified offline (analytic gradient against central finite differences to five digits in every
block; the C++ kernel against a direct 12 x 12 Gaussian marginal per pixel including regions at
beta_s = -9.7). Two mechanisms, both properties of the stated target on this sky: (i) in the plane
the region-constant fixed-T dust model misfits the sky (noise-free reduced chi-square 0.25 at the
region-mean truth, 1.8 per dof at 353 GHz) and the fit buys 101,841 in chi-square by moving both
indices; (ii) in the faint half of the sky the data-built template is wrong by factors 0.41 to 1.66
(16 to 84 percent at |b| > 20) under a prior with one global scatter s (inferred 0.36), and an index
that removes the synchrotron from every band but 23 GHz, or flattens it into the CMB, lowers the
amplitude-prior cost by more than the data cost (gain 23,818 against a GLS loss of 24,877 in the
two faintest quartiles). Consequences: a prior comparison needs a sky from the model class first
(E021); real-data runs at high S/N need a template-error model that is not a single global
scatter, and T_d sampled or marginalized (Q16); the 1023-gradient saturation seen in every E020
attempt was the funnel climb of this runaway, not a mass-matrix defect (the relative floor and the
`adapt_mass` switch are correct, tested and kept, but were not the fix). Scope: one sky, one noise
level, intensity only.

**F22. The amplitude-index coupling in the BeyondPlanck XV polarized beta_s configuration is
carried by the free CMB amplitude; the released chains are consistent with it and cannot be
shown to have settled** (E022). (a) Surrogate (LFI 30 to 70, WMAP Ka/Q/V nine-year nside-16
maps with dense covariances, NPIPE 353, four BP XV regions, prior N(-3.3, 0.1), flat per-pixel
CMB prior as in the released configuration): collapsed-to-conditional sd ratio 2.35 (Spur) and
2.11 (Plane), two-block Gibbs autocorrelation time 7 to 11 (lag-1 equals 1 - E[var_cond]/var
to 0.005 with E[var_cond] from the grid), for either LFI noise model; with a Gaussian CMB
prior sized from the Planck EE spectrum (0.19 uK at this resolution) the ratio is 1.02 to
1.04. The polarized synchrotron S/N per pixel is 2 to 3 (hundreds in the F14 intensity case);
the volume term and the stationary distribution grow with the number of amplitude elements,
so the surrogate bounds BP's native-resolution coupling from below. (b) The released chains (four, 1001 samples each, one shared start, read over HTTP from the
1.8 TB files): lag-1 0.79 to 0.83 (Spur) and 0.56 to 0.67 (Plane), IPS tau 9 to 43
consistent with the paper's 30 to 50, 23 to 114 effective samples per chain, chain means
agreeing to 0.02, and a movement towards flatter indices of +0.04 (Spur) and +0.07 (Plane)
per 1000 samples (0.3 to 2.8 sigma per chain-region, shared start). The conditional width at
BP's own resolution, from a released sample's synchrotron amplitudes through the LFI beams,
is 0.016 (Spur) and 0.0085 (Plane): the gauss C_l prior binds the synchrotron amplitude (700
to 2500 effective Stokes-elements, comparable to the surrogate's), the amplitude block is 30
(Spur) and 4 (Plane) percent of the chain's per-step motion, and a flat-prior volume term of
the surrogate's strength would move the chain by 0.03 to 0.06 per sample against the observed
5e-5, so the flat-prior shift is absent in BP's chain (pull at most 1e-3 of the surrogate's);
the mechanism is not identified, the leading candidate being that the CMB amplitude is not a
free per-pixel parameter in BP's CG step. The reported sds of 0.06 and 0.07 are consistent
with posterior widths under BP's full model at 23 to 114 effective samples per chain; the
slow movement remains unexplained. On the surrogate the flat-prior marginal is shifted from the profile by +0.27
(Spur) and +0.13 (Plane) and runs to the flattest allowed index in the High Latitude region;
the C_l-sized CMB prior removes the runaway. (c) Resolution scaling on the real data (LFI and
353, nside 16, 32, 64, white LFI noise): with flat amplitude priors the Spur's collapsed
target moves -2.88, -2.13, uniform-prior boundary (7434 pixels), the Plane -3.16, -2.88,
-2.65, with r 2.6 to 6.0 and tau 18 to 77 before pinning; the Galactic Center (high S/N) and
the profile are stable. With the C_l-sized CMB prior Spur, High Latitude and Galactic Center
move by 0.03 to 0.16 with tau 1.2 to 1.6 (High Latitude's 0.16 with a stable profile is the
synchrotron volume term), the Plane moves -3.29 to -2.88 with its profile moving the same
way: not a volume term but the misspecification of one index for a region whose index spans
1.5 (its nside-4 sub-regions are resolution-stable with a weighted mean of -3.25 at every
resolution), the real-data case for index fields at the resolution the data support. By element count
(inference, labelled; a direction, not a value) the BP XV chain's flat-prior stationary
target for the Spur lies at the boundary and the released chains are far from it.
Stats-reviewed in four rounds.

**F23. The collapsed hierarchical amplitude prior is exact only for its own template-error law;
template errors of the same size but a different conditional law bias the index where the prior
dominates** (E021 attempt 1). On a sky with region-constant indices, the model's SEDs and
noise, and templates g = a / r with r = max(1 + 0.3 eps, 0.2) iid per pixel (so that a = g r
holds pixel by pixel), the exact collapsed marginal at the true index field prefers a common
shift of the 768 lowest-S/N nside-16 dust regions of +0.114 (sd 0.0055; the determinant term
+456 nats against the chi-square term -235) and amplitude-hyperparameter means 2 to 2.5
percent below the realized mean of a/g (mu_s 0.981, mu_d 0.973, sd 0.0014); the gmrf16 chain
reproduces it (faint-quartile dust index +0.128 at posterior sd 0.033, amp mu_s 0.965) and sits
at the target's mode. The CMB prior width does not matter; redrawing a = g r from the same
templates removes the shift (peak 0.000, hyperparameters 1.000), as does the corrected sky.
Cause: with a fixed and g = a / r, p(a | g) is proportional to p(a) p_r(a/g), and since the
amplitude distribution falls steeply this is errors-in-variables attenuation with g as the
noisy variable (E[a | g]/g from 1.29 in the faintest template decile to 0.54 in the brightest,
slope of a on g 0.71, while mean(a/g) stays 1), so a | g is not N(g, (s g)^2) and the score
identity fails for the per-pixel marginal. Consequence for real data (Q16): the hier posterior
in prior-dominated regions depends on how the template errors are distributed, not only on
their size; a 30 percent scatter of the wrong law moves the faint-sky dust index by 0.11 at
four stated sigma. Exact-target computation (numpy 12 x 12 marginal); stats-reviewed
(independent marginal code, replays with fresh scatter and noise: +0.125 / +0.127 under the
data story, -0.001 / -0.009 under the model story; the corrected sky checked for other
departures from the model class, none above 0.2 sd).

**F24. On a model-consistent sky a spatial (GMRF) prior on the index field beats exchangeable
shrinkage and independent regions at nside 16 and 8, and keeps most of the per-region averaging
that coarse regions buy while retaining the structure they cannot represent** (E021 attempt 2:
ten runs on one sky, 12 bands 23 to 353 GHz at 0.03 x E012 noise, truth region-constant and
spatially smooth at nside 16, 700 post-burn draws each, every chain at the mode of its exact
target where checked; stats-reviewed, corrections taken). Index RMSE at nside 16 (sync / dust):
GMRF 0.0333 / 0.0180, iid 0.0505 / 0.0275 (paired +0.0172 +- 0.0013 and +0.0095 +- 0.0010 over
192 sky blocks, SEs understated by block correlation), independent N(mu0, 1) regions 0.1424 /
0.1009, single index 0.0669 / 0.0691; GLS CMB residual 2.89, 2.93, 3.16 and 9.03 uK against the
true-SED floor 2.79. The GMRF gains in every S/N quartile (factors 2.5 and 1.6 over independent
regions in the brightest). The same ordering at nside 8 (iid +0.0072 +- 0.0012 / +0.0052 +-
0.0008 region-level, independent +0.036 / +0.029) and, suggestively, at nside 4 (gaps 0.002 to
0.007, inside the pre-registered MC-inflation rule's margin). Mechanism of the iid loss: it
shrinks every region toward one mean whose posterior is the precision-weighted mean of the
regions, dominated by the Galactic plane (dust 1.554 against the sky mean 1.540; linear-Gaussian
prediction 1.556), so the mid-S/N regions whose truth differs are pulled 0.02 to 0.03 off, and
the biased field drags the amplitude hyperparameters off their realized targets (mu_d +6.9 sd);
that hyperparameter test detects index-model misspecification of any kind (at nside 4 every
prior fails it from the resolution misfit) and is not a sampler check. Calibration: the GMRF
over-covers at nside 16 (78 / 98 percent at 1 / 2 sd, sds 1.2 times wide), independent regions
cover 68 / 95; iid dust is biased (z median +0.40, +0.44, +0.89 at nside 16, 8, 4); coverage at
nside 8 and 4 is reference-dependent and not a calibration test. Resolution: against the
S/N-weighted region reference the GMRF's per-region error falls with coarser regions (0.0333 /
0.0318 / 0.0295 sync, 0.0180 / 0.0153 / 0.0138 dust at nside 16 / 8 / 4), so the nside-16 field
is within 13 and 30 percent of the nside-4 model that has 16 times the data per region, while
independent regions lose factors 3.5 and 4.9; the pixel-level error and the CMB residual (2.89,
4.22, 5.56 uK) follow what the field can represent; cost 125, 63, 39 gradients per iteration.
No field finer than the sky's structure was run, so the resolution rule for real data is still
open (Q16). Cost of the GMRF against iid at nside 16: 11 times the gradients per effective
sample, with the dust field mean mixing slowly (hyperparameter ESS 3 to 14 of 700) in the
centred parameterization. One sky, one noise and template realization (the seed-2028 nside-16
triple is running), a spatially smooth truth (on an uncorrelated truth the iid prior is the
correct one), fixed T_d, no template-law misspecification (F23).

**F25. A dense finite-difference Hessian metric fixes the slow mixing of the GMRF field mean in
the centred index-field NUTS at high S/N without changing the posterior** (E024, the E021 gmrf8
target and data, 1200 iterations each; stats-reviewed kernel change). The full Hessian costs the
same 2 dim gradients as the diagonal estimate (plus an O(dim^3) eigendecomposition); used as the
inverse mass (eigenvalues 1/|lam| capped at 10, refreshed at transitions 50 and 200) it gives
every index, field and amplitude statistic within 1 combined MC SE of the diagonal-metric chain
and the analyzer metrics within the MC inflation, with the dust field mean mu_d at ESS 331 of
700 against 9 (27 times more per gradient), the faint dust regions 135 against 6 (17 times), the
dust amplitude mean 376 against 39, and the chain's worst quantity 14 times better per gradient;
the sync regions, already well sampled by the diagonal metric, lose a factor 1.6 to 2.2 per
gradient (84 against 63 gradients per iteration, tree depth 7 against 6). Mechanism: in centred
coordinates the likelihood curvature is diagonal in pixel space but the basis is global, so the
couplings of the field mean and the hyperparameters with the coefficients that carry the S/N
pattern are invisible to a diagonal metric. Scope: nside 8, one sky, high S/N; nside 16 and the
low-S/N funnel untested.

**F26. A finite per-pixel CMB amplitude prior in the collapsed hier likelihood leaks the fixed,
spatially correlated CMB realization into the index likelihood at realistic noise; the induced
index shift has a realization-dependent sign and is several block sd because the per-pixel
likelihood counts correlated leakage as independent** (E023 stage 1, PySM d1s1c1 sky at the
E012 noise levels, exact-target expected-log-likelihood scans with exact templates and the
true T_d; stats-reviewed with independent code, sign-flip, shuffle and roll tests reproduced:
`experiments/E023/f26_check.py`, `f26_shuffle.py`, `f26_scaling.py`). With the E018/E021 prior
N(0, 100^2 uK^2) the target prefers the dust index of the second S/N quartile 0.05 from the
noise-free region fit (7.7 nats, 3.8 of the model's own block sd), the two fainter quartiles
at or beyond 0.10 and the faintest sync quartile at 0.10, all in the term of the expected
chi-square that is linear in the fixed CMB realization: flipping the sign of the CMB flips
every shift, rolling the CMB over the sky gives dust Q2 shifts between -0.02 and +0.02 and
dust Q3 shifts between -0.10 and +0.10, and shuffling the CMB pixels (same marginal, no
spatial correlation) gives 0.002 and 0.010 (0.2 to 0.4 sd): an iid per-pixel prior is
calibrated for an iid field and not for the real CMB (neighbour correlation 0.56 at nside 64).
The leak fraction sigma_c^2 / (tau^2 + sigma_c^2) is 0.25 at tau = 100 uK with sigma_c = 58 uK
(the per-pixel CMB noise at noise scale 1) and 1.7 uK at 0.03 x noise; the shift scales as
noise^2 at fixed tau (dust Q2: 0.046, 0.008, 0.001 at noise 1, 0.3, 0.1; about 1e-4 at 0.03),
so E018, E021 and E024 are unaffected. The true per-pixel variance (tau = 84) makes it larger
(dust Q2 0.058, 11 nats); the flat CMB prior (`amp_prior: none`, or sd 1e4, leak 3e-5) projects
the CMB out exactly per pixel, is independent of the CMB realization and its correlation, and
costs the CMB column's one degree of freedom (block sd 20 to 25 percent wider in the faint dust
quartiles); with it every noise-1 baseline scan peaks within 0.03 of the reference at gains
below 0.05 nats (dust Q2 +0.01 at 0.6 sd by shift). The correlated CMB prior of the CG step
(Q5) is the only choice that is both calibrated and uses the CMB's smoothness. Scope: one sky,
one CMB realization, template-Fisher quartiles, the noise-free region-fit reference, scan-edge
peaks in the fainter quartiles; the per-pixel kernel on real maps at this depth must not carry
a finite CMB prior (D021).
