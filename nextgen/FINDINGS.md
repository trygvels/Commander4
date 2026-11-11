# Findings so far

F1 to F11 date from the exploratory session of 2026-10-25. F12 to F45 were added as experiments
closed, through 2026-11-09. Experiment IDs refer to `experiments/registry.md`. Numbers are from
simulations unless stated.

## 0. Reading guide after the 2026-11-09 review

The review is `reviews/2026-11-09_program_review.md`, with the referee report and the pipeline
audit beside it. Five things apply to everything below.

- **Noise regime.** "Realistic noise", "Planck-like depth" and "noise scale 1" in F26 to F38 mean
  scale 1 of the PySM ladder, which is 24 to 90 times real WMAP/LFI intensity white noise (F39).
  Scale 0.03 (F21, F23, F24) is 0.7 to 2.7 times the real depth. Only F22 concerns polarization.
- **Baselines.** On the PySM d1s1 sky a constant index at the truth mean has pixel RMSE 0.0629
  (sync) and 0.0374 (dust); the truth's own nside-16 region means have 0.0130 and 0.0187. An
  index-field RMSE means something only against these. Every scale 1 synchrotron field with the
  sd 0.5 or cs-template configurations is worse than the constant.
- **The CMB residual.** The total flat-prior GLS CMB residual quoted in F32 to F38 is the noise
  term at scale 1 and ranks nothing (F40).
- **One realization, one chain.** F26 to F38 rest on one sky, one noise realization and single
  chains. Each run computed a second chain that was never written (F42).
- **Real maps.** Every E014 number expressed in posterior sd is the width of a fit whose
  high-latitude reduced chi-square was of order 13 to 60 in the chain's own noise model (F41).

Grades: A, established against an exact or independent reference with replication. B, supported
on one sky, seed or configuration. C, suggestive. D, possibly an artefact. W, withdrawn. They
follow the referee report except where a later experiment changed them (F1, F38).

| Finding | Grade | Status after the review |
|---|---|---|
| F1 | A for the mechanism and the law (derived twice independently), B for the real-map size | Restated (E035): Monte Carlo within Metropolis; temperature 1/2 + 0.443 times the sd of the chi-square at a fresh draw, sqrt(pi k / 8) + 1/2 for exact stored draws; the real-map width and mean shift explained by the alm storage (F7) |
| F2 | A | Stands |
| F3 | B | Mechanism stands and is in the literature (0709.1058, 0811.4280, 1608.00551); "54 ridge sd" dropped (F41) |
| F4 | C | "Principled default" is a convention, not a result; literature sentence corrected |
| F5, F6 | C | One seed per cell; the qualitative result is known |
| F7 | B | Amended: the stored mean maps are not band-limited either; 1.9 k of chi-square inside the E014 mask; stored amplitude samples under-dispersed; maps instead of alms (Q6) is now a correctness item |
| F8 | B | Amended by E037: over 40 seeds the offset is a zero-mean scatter of 4.4 and 4.8 posterior sd; the white-noise control alone gives 2.6 and 7.8, from the per-pixel S/N of 0.04 to 0.06; no bug |
| F9 | C | Code reading, confirmed by a second reading (audit G1); no experiment |
| F10, F11 | B for the code facts, C for their impact | Commander3 was never run |
| F12, F15, F18 | A on simulations | Stand |
| F13 | B | Stands |
| F14 | B for the ratio, C for iterations per sample | One configuration |
| F16 | B | Calibrated on data drawn from its own hyperprior; says nothing about misspecified skies |
| F17 | B | A statement about means; the like-for-like shift is +0.035 (beta_s) and -0.016 (beta_d); "0.4" compares two band sets; "hundreds of sd" dropped |
| F19 | B | One replicate; no SBC of the Commander4 implementation |
| F20, F25 | B minus, C | Efficiency on one seed |
| F21 | B plus | The most important negative result: a runaway at real depth with a data-built template |
| F22 | B for (a) and (b), D for (c) | (c) is an inference, not a measurement |
| F23 | A for the mechanism, B for the size | Stands |
| F24 | B plus | The one place an index field beats a scalar; needs the same comparison on PySM at scale 0.03 with longer chains |
| F26 | A | Algebra |
| F27 | B to C | Scale 1 |
| F28 | W at scale 1 | Contradicted by F37, F38 and chain 3. At scale 0.03 the level is unbiased (F45) but the field does not beat a constant index |
| F29 | B to C | Wording of the slope corrected; sizes from 6 to 8 realizations in a review note |
| F30 | W | Withdrawn in 2026-11-03 |
| F31 | A | Stands |
| F32 | B for the numbers, W for the rule | All three rungs lose to a constant index for synchrotron at scale 1 |
| F33 | C | One chain with exact templates; not a rule (F37) |
| F34 | C | One chain |
| F35 | C | "Validated" withdrawn: the reference was refined after pre-registration, 22 divergences |
| F36 | C | Holds for one sky with exact templates, the oracle T map and a width chosen knowing the template is exact |
| F37 | B minus for the scans, C for the rule | Scale 1 |
| F38 | B for the description, W for the cause and the remedy | Chain 3 amendment below; absent at real depth (F45) |
| F39, F40 | B | One sky; Fisher at the truth |
| F41 | A for the floor, B for the order | Real maps and one simulated CMB |
| F42, F43 | A | Code facts |
| F44 | A for the flat-prior shift and the calibration of the profile at real depth, B for the rest | Exact grids, 100 noise realizations per cell, one sky; the zero-mean prior is qualified by F46 (e) |
| F45 | B | Four chains, one sky; closes the intensity template-prior sequence |
| F46 | A for (a) to (d) in one dimension, B for (e) | Exact one-dimensional targets, 400 realizations per cell, pre-registered; the coupled dense case is E039 |
| F47 | B | Deterministic input checks; the LFI noise model is empirical |
| F48 | A for the shift, B for the values | Pre-registered; real maps; see F49 for the values |
| F49 | B for the size of the scale dependence, W for "beam aliasing" as its cause | Amended 2026-11-11 after the referee read: cause not established, estimator not yet validated on simulations, stability criterion failed in two regions |
| F50 | B, C for the leakage reading | Band pairs pre-registered; the PR4 against PR3 part is exploratory |
| F51 | B for the numbers, D for the reading of published chains | Not their maps, pixelization or chains |
| F52 | A | Measured with a guard allocator and a malloc canary; reruns compared digit by digit |
| F53 | B | One data set; the white prior was the pre-registered baseline and the pre-registered criterion failed |
| F54 | B | Estimator checked on simulations with and without a true spread (60 realizations each, nside 8); one data set; 5 to 9 effective patches per region |
| F55 | A for the differences between estimators, B for the index values | Pre-registered (H1 to H3); regions by geometry with a mask; one data set; values depend on the choices of F58 |
| F56 | B | 30 simulated realizations at nside 8 and one map fit at nside 16; the size on the maps is not predicted |
| F57 | A for the Tau A test, B for the reading of the block errors | Pre-registered (H4, H6); null of 120 realizations at nside 8, the map fit at nside 16 |
| F58 | A for the calculation and the scan, B for its truth model | Pre-registered (H5, H11); first order in the mismatch; the LCDM covariance of the paper taken as the truth |
| F59 | B | Pre-registered (H7, H8); one data set; 24 correlated comparisons, largest 2.4 expected sd |
| F60 | B for the deficit, C for its reading | Not pre-registered (found on reading R1 again); a simulated null; the literature has no statement either way |

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

**Restated 2026-11-09 (E035).** Taken alone the move is a valid pseudo-marginal
Metropolis-Hastings kernel, but for an extended target whose index marginal is the profile (ridge)
posterior and whose amplitudes have half their posterior covariance. The amplitude step of the
Gibbs sweep redraws the amplitudes with their full covariance, which is not that target's
conditional, so the composite chain leaves neither target invariant. With one index step per
amplitude draw it converges to the ridge posterior tempered by T = sqrt(2 pi k) / 4, with k the
number of amplitude degrees of freedom re-solved per proposal: the ridge mean, and a width of
0.7916 k^(1/4) times the ridge width. The law is derived in the E035 record and holds in toy
chains to 0.7 percent for k from 1000 to 190825; without the amplitude refresh, or with the
refresh drawn from the extended target's conditional, the same toy returns the ridge posterior
exactly. More index steps per amplitude draw lower the coefficient (0.64 at five, 0.54 at
twenty). The simulated rows of the table agree with the toy at their own settings to 5 percent
(E002: 4.0 to 4.1 and 5.7; E005: 3.0 and 4.3). The real-map row sits 22 to 31 percent above the
toy at E014's settings (13.1 to 13.6 against 16.4 to 17.0 in two chains, beta_d known to 2
percent), and the independent stats review traced that to the alm storage of the per-pixel
amplitudes (F7): the chi-square the move sees at a fresh draw has 2.44 +- 0.06 times the exact
variance (measured from the chain's own log; 2.69 +- 0.18 computed directly on the inputs), and
the toy with that factor gives 15.95 and 16.57 with the chain's acceptance. In general the
temperature is 1/2 + 0.443 times the sd of the chi-square at a fresh draw, which equals
sqrt(pi k / 8) + 1/2 when the stored amplitudes are the drawn ones. The exact small-step law at
every k is T = 1/2 + 1/(4 c_k), c_k = 2^-(k+1) C(k, k/2). The temperature law for a noisy
acceptance is in print for the sign rule (Ball, Fink and Bowler 2003; unverified in the PDF by
the main session), and the scheme is the "Monte Carlo within Metropolis" of Andrieu and Roberts
2009; new here are the identification of the move, the chi-square form and the k^(1/4)
consequence. Commander4's
`perpix_betasamp` example config uses this move at nside 512, and an independent subagent
reading the upstream code reached the same conclusion about it. On the real maps of E014 the
dust index is also shifted by 3 ridge sd, which is not a property of the move: the gradient of
the storage residual predicts +3.1 (F7). Scope: flat amplitude priors and the per-pixel
amplitude path. In the CG path the re-solve is a truncated solve (audit A1) and nothing is
established. Superseded wording, kept for the record: "not a valid MCMC step", an empirical
k^(1/4) fit with a predicted 18.5 for E014, and "factors above 50" at nside 512 to 1024.

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
**Caveat (2026-11-09, F41):** the "54 ridge sd" and "5 sd" are widths of a fit that does not hold
at high latitude in the chain's own noise model, and are dropped. The absolute shift of +0.019
belongs to that error model. The mechanism is published (Eriksen et al. 2008, Dunkley et al. 2009,
Alonso et al. 2017); what is not published is a coverage curve against S/N (novelty check).

**F4. Ridge is the principled default.** Ridge is the profile ("spectral") likelihood of parametric
component separation (Stompor et al. 2009). As an MH target it equals the flat-prior marginal times
det(A)^(1/2), i.e. the marginal under the Jeffreys prior of the amplitude block conditional on beta,
p(a | beta) proportional to det A(beta)^(1/2), with a flat prior on beta. (Corrected 2026-10-25 after
the E012 review: it is NOT the full joint Jeffreys prior sqrt(det I(a, beta)), which carries an
a-dependent factor and gives a different beta posterior; quadrature check in
`experiments/E009_c3_nonlin_review/check_a_e.py`.) Eriksen et al. 2008 used a Jeffreys prior to
remove a bias of more than 4 sigma in the global offset driven by a uniform index prior (their
Sect. IV.2), and Alonso et al. 2017 (1608.00551, Fig. 11) report index biases of more than 4 sigma
from the volume factor and apply a volume prior that makes their target the ridge likelihood
(corrected 2026-11-09 from the novelty check; verify both in the PDFs before citing). `auto` now resolves to ridge where available. Conditional Gibbs (Commander3
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

**F7. Per-pixel amplitudes are stored as alms with lmax = 3 nside - 1** (E007; amended
2026-11-09, see the end of this entry). For unsmoothed maps
the round trip alone raises reduced chi^2 from 0.57 to ~794 (matching Commander4's own logged value);
raising lmax makes it unstable. With common-resolution smoothing the mean is fine, but the white
per-pixel fluctuation draws still lose variance above lmax. Not changed yet.
**Amended (2026-11-09, stats review of the legacy move, `reviews/2026-11-09_stats_review_legacy_move.md`
Sect. 7d):** "with common-resolution smoothing the mean is fine" does not hold at the chi-square
level. Computed on the E014 inputs (nside 64, 120 arcmin, lmax 191, 38165 masked pixels, five
components), the stored amplitudes are P(a_hat + delta) with P the alm round trip, and the
chi-square at stored amplitudes exceeds the minimum by a deterministic 1.90 k (9.5 per kept pixel,
k = 190825), from the round-trip residual of the mean amplitude maps, which are not band-limited
at lmax = 3 nside - 1 (residual 0.5 to 1.6 marginal sd per component; not an iteration artefact
and not confined to the mask edge). The fluctuation keeps 0.74 of its chi-square variance, so
every stored amplitude sample is under-dispersed. Consequences: any acceptance evaluated through
`global_loglike` sees chi2_min plus this offset. For the legacy move that multiplies the
acceptance noise by 2.44 and shifts beta_d by +3 ridge sd through the offset's gradient (F1).
The `conditional` mode accepts through the same function and is probably affected; unchecked.
`ridge`, `marginal` and `hier` accept on the pixel-space kernel and are not affected, so no
index result of F12 to F38 changes. The amplitudes written to the chain files are affected in
every per-pixel run. This re-implements the draw and the round trip with Commander4's library
calls and has not been confirmed by a Commander4 run with map storage (Q6, Q23).

**F8. The per-pixel model treats smoothed noise as white** (E006; amended 2026-11-10 by E037, see
the end of this entry). After smoothing to a common beam
the noise is correlated; the per-pixel likelihood ignores that, and in the validation simulation
the exact ridge target sat 1.9 sd (beta_s) and 3.7 sd (beta_d) from the input. This is a model
issue, not a sampler issue (Q7).
**Amended (2026-11-10, E037: 40 seeds of the harness simulation, exact grids, a white-noise
control):** the 1.9 and 3.7 sd of the single realization are typical draws of a zero-mean scatter
with sd 4.4 (beta_s) and 4.8 (beta_d); there is no significant bias (-0.18 +- 0.69 and -0.51 +-
0.76). The smoothing is not the only cause: with white noise and no smoothing the same ridge
target is overconfident by 2.6 and 7.8, because the harness has a per-pixel amplitude S/N of 0.04
to 0.06, where the curvature of the profile likelihood understates the error. With the noise
divided by 10 and by 30 the white control is calibrated (sd of z 0.82 and 1.07, then 0.80 and
0.97), so there is no bug in the harness or the target. The harness tests a sampler against its
target (F2 is untouched); its target is not a calibrated statement about the input indices. The
split between the smoothing and the per-pixel S/N is not measured (Q7).

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
**Caveat (2026-11-09, F41):** both widths belong to an error model with the HFI bands at 235 to
760 times their intended weight. The ratio of the two is the finding; neither width is an
uncertainty of the sky's index.

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
**Corrected (2026-11-09, referee Sect. 2.5 and F41):** the like-for-like effect of the zero-level
convention, with the same 20 bands, is +0.035 in beta_s (-2.9103 to -2.8755) and -0.016 in beta_d,
plus +0.011 from the Haslam dipole. The "0.4" of the headline compares two band sets (17 against
20 bands) and is not a zero-level effect. "Hundreds of posterior sd" is dropped: the sd is that of
a fit with a high-latitude reduced chi-square of tens. The ingest errors of F41 act on the same
maps and were not separated from the zero levels.

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
centred parameterization. Replicated on a second noise and template realization (seed 2028):
gmrf16 0.0323 / 0.0180, iid16 0.0501 / 0.0271 (paired +0.0178 +- 0.0013, +0.0091 +- 0.0010),
indep16 0.1373 / 0.1074, iid dust mean 1.553, iid hyperparameter test failing (mu_d +5.6 sd)
and the GMRF passing: the realization scatter of the gaps is about 0.0005. One sky, a
spatially smooth truth (on an uncorrelated truth the iid prior is the correct one), fixed T_d,
no template-law misspecification (F23), the per-pixel CMB prior inert at this noise (F26).

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

**F27. Under the collapsed hier amplitude prior, a template built from one scaled, smoothed
frequency map biases the spatially varying index field at 0.03 x and at 1 x the E012 noise,
and no scalar or region-scale hyperparameter law repairs the dust template** (E023 stage 1,
PySM d1s1c1 sky, exact-target expected-log-likelihood scans of Commander4's kernel, model
control on the E021 sky PASS, kernel agreement 4e-12; stats-reviewed design). With the E020
templates (23 GHz scaled to 30 GHz for sync, 353 GHz scaled at a fixed T to 545 GHz for dust,
5 deg smoothing, floored) and global hyperparameters, the target prefers index shifts of 0.05
to beyond 0.1 in every S/N quartile but the brightest at 0.03 x noise (faintest sync quartile
toward beta_s -3.6) and, at realistic noise with a flat CMB prior, -0.09 in the brightest sync
quartile (3 sd), beyond +0.1 in sync Q2 and dust Q2 and Q3, +0.01 in the brightest dust
quartile (3 sd), with the faintest sync quartile running to beta_s near -5; with exact
templates the same target is unbiased at realistic noise and off by at most 0.03 (fixed T) at
0.03 x. A scatter per quartile removes the 0.03 x runaway only; a mean and scatter per quartile
fixes the mid-S/N sync and leaves the dust at 0.05 to 0.10 (0.03 x) or biases the brightest
dust quartile by +0.03 at 8 sd (noise 1); a free multiplier per nside-16 region (3072 per
component) fixes the sync and leaves the dust untouched: the dust template's error is
pixel-scale and tied to T_d, hence to the index itself. Fixed T_d is at most a 0.03 effect at
0.03 x and invisible at noise 1. At realistic noise the faint half of the sky is
prior-dominated at nside 16 (block sds 0.08 to 0.12 sync, 0.03 to 0.05 dust), so a wrong
template is the posterior there. Consequence (D022): hier for synchrotron with a region-scale
template-error field; for dust only with a component-separation template or sampled T_d, to
be tested in E025. Scope: one sky, one realization per cell, the data term without the index
prior, fixed against known T only.

**F28. With a template-error field at the region scale the hier amplitude prior is close to
unbiased for a component-separation-quality template at 0.03 x and 1 x noise, and the dust
temperature misfit is below everything else** (E025, exact-target scans of Commander4's kernel
on the PySM d1s1c1 sky, flat CMB prior, template g = a exp(f) with f a 10-deg, 30 percent
Gaussian field uncoupled to the sky). Under the global law (one mu and s per component) the
target shifts the two fainter dust quartiles by 0.07 to 0.09 at 0.03 x noise and beyond 0.1 at
realistic noise (the F23 law mismatch on a realistic error); with one free template multiplier
per nside-16 region (3072 per component, N(1, 0.5^2) hyperprior each) plus a scatter per S/N
quartile every shift is within 0.02 at 0.03 x (dust Q4 -0.02 to -0.03 at 2 sd, the only
excess over the pre-registered 0.01) and no scan at realistic noise has a gain above 1.6 nats
(14 of 18 inside the block sd, four at 1.0 to 1.7 sd). The profiled residual scatter drops to
0.06, the within-region part of the error. The dust SED at the true T_d against the fixed T
changes no shift by more than 0.01 in any of eight cells. Consequence (D023): Commander4's
hier prior needs mu as a region map (a template-error field), after which a template whose
within-region error is at the few-percent level is admissible; sampled T_d is deprioritized.
Scope: one sky, one error realization, one error scale (10 deg), the data term without the
index prior.
**Recheck under the analytic multiplier marginal (F31; E025 second pass, 2026-11-02):** the
template effect (cs minus exact) is at most 0.02 at 0.03 x and 0.04 at 1 x in every quartile,
so this finding stands; the absolute shifts of the hier target against the true-T region fit
are model properties (the fixed-T misfit at 0.03 x, the F29 Occam push at 1 x: +0.10 to +0.16
sync, +0.02 to +0.10 dust in the expected target), which the profiled scans had hidden.
**Qualified (2026-11-08, E030, F37):** "at most 0.04 in every quartile at 1 x" was read on a
+-0.16 grid with both templates at the edge for the faint sync quartiles. Scanned over +-0.3
with the oracle T map, the template effect at 1 x is at most 0.01 in dust Q1 and Q2 and 0.04 in
sync Q1 (where the data decide), but +0.06 to +0.11 in dust Q3 and Q4 and beyond +0.25 in the
faint sync quartiles at every hyperprior width 0.1 to 0.5: in the prior-dominated half of the
sky a 30 percent template error is not removed by the multiplier field, and the GMRF index
prior sets those regions in a chain.
**Withdrawn at noise scale 1 (2026-11-09 review):** four chains with the multiplier field and a 30
percent template put the synchrotron field mean at -2.72, -2.70, -2.72 and -3.16 against -3.00
(F38). The scale 0.03 statement rests on exact-target scans (E032 scan: every quartile within
0.02) and awaits the E032 chains.

**F29. At realistic noise the nside-16 GMRF synchrotron index field is prior-dominated over
three quarters of the sky, and there the hierarchical amplitude prior's per-region multipliers
bias the index flatter through their marginalization volume** (E026 gmrf16_exact_mumap_nc:
PySM d1s1c1 at the E012 noise, exact templates, one multiplier per nside-16 region with
hyperprior N(1, 0.5^2), flat CMB prior, non-centred NUTS, chain validated by the D025 check;
independent stats review in `E026/stats_review_sync_level.md`). The bright template-S/N quartile
holds 90 percent of the sync Fisher weight and agrees with the noise-free region-fit reference
(the 10 to 300 regions with the largest weight: -0.02 to +0.03 from the reference, z medians
-0.6 to +0.2, 77 to 90 percent inside |z| < 1); the Fisher-weighted level of the field is
-2.978 against -2.988 for the truth. The unweighted field mean is -2.868 +- 0.042, +0.13 above
the sky mean, because the 2300 regions with region Fisher sd 2.5 to 3.7 sit at +0.12 to +0.14.
Cause (exact marginal with the multipliers integrated analytically): the volume term -1/2
sum_r log(-Q_r) of the 3072 marginalized multipliers rises with a flatter index (slope +130
nats per unit beta over the faint regions; 66 percent of the sync amplitude information is at
23 GHz, so with the reference at 30 GHz a steeper index gives each multiplier more leverage,
which the volume term penalizes; wording corrected 2026-11-09, the text said "rises with a
steeper index" against its own slope sign), which puts the faint level
+0.15 +- 0.05 above the truth over simulated noise realizations (+0.19 on this one); the chain's
+0.13 is that bias held back by the field prior. The bias scales with the hyperprior variance:
sd 0.2 gives +0.06, sd 0.1 zero within noise (E029 tests this in a chain). It is not a property
of the field's roughness (a first version of this finding said so; refuted by the review: the
block-shift optimum is linear in the coefficient scale, and the faint-region data level is the
same for constant, reference, shuffled and iid patterns). A single sync index on the same data
is -2.944 +- 0.031 against a model pseudo-true value of -2.956 +- 0.003 (the same Occam push
diluted by the full-sky curvature) and an exact 1-d marginal of -2.951 +- 0.025: consistent.
Consequences: the unweighted mean of a prior-dominated index field is not a reportable
quantity under this hyperprior; the multiplier hyperprior width is a scientific choice that
sets the faint-sky index; the dust state is a sync confounder at this depth (with the dust at
the truth and mu_d at 1 the scalar pseudo-true value moves to -2.83). The same push acts on
the dust temperature: with the multipliers integrated, the expected target over a scalar T
peaks at 27 K (sd 1.7) under the sd-0.5 hyperprior and at 20.2 K (sd 0.5; sky 20.9) under sd
0.1 (E028 T grids). Scope: one sky, one noise realization; six to eight simulated noise
realizations for the pseudo-true values.

**F30. With the per-region template multiplier field at realistic noise, a
component-separation-quality template is usable for dust and not for synchrotron, the scalar
law fails, and the fixed-T dust misfit becomes a bright-sky dust-index bias through the
multipliers** (E026 non-centred chains, same setting as F29). Against exact templates the cs
template (true amplitude times a 10-degree, 30 percent log-normal error field) costs +0.19
(sync) and +0.02 (dust) pixel-level RMSE; dust's bright quartile is as good as with exact
templates (0.069 against 0.072), the fainter quartiles 1.3 to 1.7 times worse with no extra
bias, while sync is 2 to 3 times worse everywhere with the field scale at 0.35; the
multiplier maps track the realized error with correlation 0.7 to 0.9 in the bright quartiles
and 0.2 to 0.4 elsewhere (posterior sds 0.25 to 0.48). The scalar law (one multiplier and
scatter per component) gives +0.48 and +0.24 with a sync field scale of 0.87: the predicted
failure (F27, F28), now in a valid chain. With exact templates the bright dust quartile sits
+0.06 above the reference with the dust mu map at 1.029 +- 0.007 (should be 1): the fixed-T
misfit (beta_eff3 minus beta_d rms 0.04) is absorbed jointly by the multipliers and the index,
1191 + 375 nats better than the reference, a path the global law could not take (E025: +0.5
percent, no index shift). So the region-scale multiplier field, which removes the template
bias (F28), converts an SED misfit into an index bias where the data are strong: region-wise
or sampled T_d comes before the multiplier field is used for dust on real data (M2.6).
**Correction (2026-11-02, E028 baseline, `E028/fixedT_baseline_vs_trueT.txt`):** the +0.06 is
relative to the fixed-T region fit (`bd_fixT`, whose pseudo-true value is the effective index
beta_eff3 of a fixed-T fit). Against the region fit at the true per-pixel T (`bd_trueT`, the
true beta_d) the same chain's bright quartile has z median -0.06 (mean offset +0.006, RMSE
0.055): with the per-region multipliers taking the T misfit at the 545 GHz pivot (the 1.029),
the chain's bright-sky dust index is the true beta_d, not beta_eff. The fainter quartiles sit
+0.06 above bd_trueT (z +1.0 to +1.2; the field mean 1.578 against a sky mean of 1.540), a
prior-side effect to be read with E027. So the multiplier field does not bias the bright-sky
index; it changes which index the chain estimates (beta_d instead of beta_eff) and moves the
T misfit into the multipliers. The recommendation stands for a different reason: on real data
the multiplier absorbs any pivot-frequency mismatch, so its 3 percent is not interpretable as
a template error, and a T model is needed to separate the two.
**Withdrawn (2026-11-03, E028 stage A, F34):** the mechanism statement "the multipliers absorb
the fixed-T misfit" and the correction's "the chain's bright-sky dust index is the true
beta_d" do not survive the oracle-T chain: with the true region T map the mu_d map stays at
1.025 and the bright quartile sits +0.047 above beta_d. The fixed-T chain's agreement with
beta_d was a cancellation of two offsets of about 0.05 with opposite signs. The recommendation
(a T model before the multiplier field is used for dust on real data) stands on the E029 and
E028 evidence that neither hyperprior width gives an unbiased dust index with a fixed T.

**F31. Profiling per-region nuisance multipliers in the exact-target scans is not a valid
stand-in for the chain's marginal; the multipliers must be integrated analytically** (E026
stats review, `E026/review_scripts/marg.py`, `occam.py`). The per-pixel hier marginal is
Gaussian in the region multiplier mu_r, so the integral over mu_r under its hyperprior is
closed-form (F(mu_hat) - 1/2 log(-Q_r) per region; checked against brute force to 1e-8). The
E025 and E026 scan cells with 3072 profiled multipliers carry a Neyman-Scott bias toward
steeper indices (-0.2 to -0.4 on the E026 chain state against -0.003 +- 0.025 for the
marginal), and the global-multiplier cells are a different model (no per-region volume term;
-0.109, 7.3 nats, where the marginal gives 0). Conclusions drawn from profiled per-region
cells (F28's "within 0.02 at 0.03 x noise, insignificant at 1 x" for the cs template; the E026
whole-sky and chain-state scans) are to be redone with the marginal before they are cited;
the E026 chain results (F30) do not depend on the scans. The scan script's `perregmarg` law
now does the integration (checked against a 2-d quadrature to 5e-12); at the E026 chain's
posterior-mean state it gives a whole-sky optimum of +0.000 (sd 0.026) where the profiled law
gave -0.12 to -0.18, reproducing the review.

**F32. At Planck-like depth with exact templates the index-field resolution between nside 4
and 16 changes the pixel-level index error by 6 to 10 percent, with the optimum at nside 8
for synchrotron and at the finest rung (nside 16) for dust; the sync level offset (F29) and
the dust multiplier offset (F30) do not depend on the resolution** (E027, three valid chains
by the D025 check, paired block-bootstrap differences). Sync RMSE 0.143, 0.135, 0.147 and
dust 0.067, 0.063, 0.060 at nside 4, 8, 16 (differences 0.008 +- 0.001, 0.012 +- 0.001;
0.006 +- 0.001, 0.002 +- 0.001); the GLS CMB residual follows the sync field (nside 8 best by
0.07 and 0.22 uK over the 93.2 uK floor). The region Fisher sd exceeds the sky's index sd at
every rung for both components (sync 0.65 to 2.7 against 0.063, dust 0.094 to 0.44 against
0.037), so "the resolution where the per-region information crosses the sky scatter" is not
the rule; under the GMRF prior the finer field stays competitive and the coarse field pays the
representation floor (within-region truth sd, half the sky sd at nside 4). The sync field mean
sits +0.13, +0.12, +0.13 above the sky mean and the mu_d map at 1.030, 1.029, 1.029 at the
three rungs: both offsets are amplitude-model properties (the 3072 multipliers are the same at
every rung), as F29 and F30 predict. Pre-registered dust reading (b) (faint quartiles better
at nside 8 than 16) refuted: the finer rung is better in every dust quartile, most in the
faintest (0.046, 0.039, 0.032). Costs 31, 31, 63 gradients per iteration.

**F33. A tight multiplier hyperprior N(1, 0.1^2) removes the faint-sky synchrotron index bias
of the hier model at Planck-like depth and restores the field scale, but moves the fixed-T
dust misfit from the multipliers onto the dust index and the high-latitude CMB** (E029 against
E026, same data, both chains pass D025). Sync: field mean -2.981 +- 0.028 (sd 0.5: -2.868;
truth -3.000), faint quartiles within 0.02 of the region-fit reference (sd 0.5: +0.12 to
+0.14), sigma_s 0.077 +- 0.040 (sd 0.5: 0.152; truth 0.063), pixel RMSE 0.061 against 0.147.
Dust: the mu_d map falls from 1.029 to 1.003, the bright quartile moves from beta_d to the
fixed-T effective index (z -1.55 against bd_trueT, +0.11 against bd_fixT; a move of 0.054),
the faint half shrinks to the field mean 1.506 (0.08 below the fixed-T reference), and the
flat-prior GLS CMB residual rises by 0.91 +- 0.07 uK (high latitude 92.91 to 93.95 uK; the
excess over the true-SED floor from 8.9 to 15.8 uK rms; by hybrid fields the better sync field
is worth 1.0 uK and the dust field costs 2.1 uK, 0.75 of it the faint half; the chain's own
CMB draws degrade by 0.5 uK). The hyperprior width is a model
choice that sets the faint-sky index (F29 confirmed in a chain); it does not replace a dust
temperature model (F30): with a fixed T, sd 0.5 hides the misfit in the multipliers and sd 0.1
exposes it in the index. The stage A chain (E028, oracle T) tests the combination.

**F34. The fixed dust temperature is not the cause of the dust multiplier offset (1.029) and
the bright-sky dust-index offset of the sd-0.5 hier chain: with the true region T map the
multipliers stay at 1.025 and the dust index rises everywhere, so the sd-0.5 model's
pseudo-true dust index sits about 0.05 above its own reference under either T model** (E028
stage A, `gmrf16_exact_mumap_nc_Tmap`: E026's chain with the oracle nside-16 region-mean T map
in the dust SED, everything else equal; the chain passes the D025 check against its T-map
target, KS p 0.12 and 0.39, and fails it against the fixed-T target, dust KS D 0.98, 1269 nats
lower). Both pre-registered criteria fail: mu_d map 1.0254 +- 0.0073 (required 1.00 +- 0.01)
and bright-quartile z median +1.31 against `bd_trueT` (required within 0.5). The dust index
moves up by 0.04 in every quartile (offsets against beta_d +0.047, +0.095, +0.098, +0.066
against +0.006, +0.061, +0.059, +0.030 with the fixed T; field mean 1.616 against 1.578, sky
mean 1.540; pixel RMSE 0.085 against 0.060) while the sync field is unchanged (field mean
-2.880 against -2.868, quartile offsets within 0.02). The multiplier offset lives in the bright
half (dust mu map per quartile 1.038, 1.041, 1.017, 1.005; fixed T 1.062, 1.038, 1.012, 1.002;
sd 0.1 with fixed T 1.012, 1.000, 1.000, 1.000). Against each model's own reference the
bright quartile is +0.047 (oracle T) and +0.060 (fixed T, where beta_eff3 minus beta_d is
-0.054): F30's reading that the fixed-T chain estimates beta_d was a cancellation, now
withdrawn. The flat-prior GLS CMB residual rms is not an index-quality criterion: the oracle
chain's biased-high dust field gives 92.3 uK, 0.7 uK below the true-SED GLS floor of 93.2
(a wrong mixing matrix trades bias for variance), and putting `bd_trueT` in its place raises it
to 94.4. Consequence: with the N(1, 0.5^2) multiplier hyperprior the dust index and multipliers
carry an offset that a correct T model does not remove, and the sd-0.1 hyperprior with the T
model (stage A') is the configuration that can be unbiased. Mechanism (exact-marginal scans
in the E028 record, multipliers integrated, read against beta_d): under sd 0.5 the F29 Occam
push puts the three fainter dust quartiles at +0.07 to +0.10 with either T model (the stage A
chain shows +0.07 to +0.10, the fixed-T chain +0.03 to +0.06 held back by the field prior),
the bright quartile's pseudo-true offset is +0.01 to +0.02, and the chain's +0.047 there is
the GMRF prior transmitting the faint-sky level (field mean 1.616 against a bright-quartile
level of 1.57). Under sd 0.1 with the true T map the optimum is 0.000 in every dust quartile
and +0.01 for the sync level: the stage A' prediction. Scope: one sky, one noise realization,
exact templates.

**F35. A sampled scalar dust temperature in the collapsed hier index-field likelihood is
exact for its stated target and lands at the exact-marginal optimum of its own law: T =
19.19 +- 0.54 K against a pseudo-true 19.90 +- 0.47 K (1.5 sd), on a sky whose unweighted
mean T is 20.9 K and dust-weighted mean 22.6 K; the dust index field moves in step along the
beta-T degeneracy and ends within 0.03 of beta_d in every quartile by cancellation** (E028
stage B, `gmrf16_exact_mumap_nc_Tscalar`: E029's chain with `sample_T: true`, `T_nside: 0`,
prior N(log 20, 0.2^2), sd-0.1 multiplier hyperprior; 7.0 h, 700 post-burn draws). The D025
conditional-CDF check passes along sync, dust and log T (KS p 0.59, 0.24, 0.41; conditional
sd of log T 0.0024, so the 0.54 K posterior width is the beta-T-multiplier degeneracy, not
noise). The pre-registered reference (20.2 +- 0.5 K, a coarse read of a three-point grid) was
refined by the same construction to 19.90 +- 0.47 K before the comparison; the chain's 0.72 K
offset below it is not Monte Carlo noise (0.07 K) and lies in the index fields: with the
chain's own posterior-mean fields held fixed the same exact marginal peaks at 19.14 +- 0.12 K,
0.4 sd from the chain (the GMRF-shrunk faint half is 0.03 to 0.05 shallower than the free
refit at 19.9 K). The implementation was stats-reviewed after the run (target and gradient
sound; two hazards in unused configurations, a fixed T_map beside the scalar and the scalar
restart, fixed with tests). Against the region fit at the posterior T the dust field
is at its reference in the bright half and 0.03 to 0.05 low in the faint half (the E029 and
sd-0.1 scan pattern); against beta_d it is -0.004, +0.029, +0.024, -0.002 (pixel RMSE 0.046
against E029's 0.056) because the region fit at 19.2 K lies above beta_d by -0.016 to +0.063
per quartile and the faint-half shrinkage cancels it. Sync is unchanged from E029 (field mean
-2.977 +- 0.028, RMSE 0.062), the multiplier maps are at 1.002 and 0.999, the dust region ESS
improves (min 76, median 602 against 35, 136) at the same gradient cost, with 22 divergences.
Consequence: the sampled-T code path is validated at realistic noise (the L3 rung for the
stage C T field), and a scalar T is a degeneracy-weighted compromise, not the sky's dust
temperature; the spatial T model with 545 and 857 GHz remains the structural fix. Scope: one
sky, one noise realization, exact templates.

**F36. With the oracle dust temperature map and the tight multiplier hyperprior N(1, 0.1^2),
the hier collapsed index-field chain has no likelihood bias in either component at realistic
noise: the bright dust quartile is at beta_d (+0.010, z +0.42), the multiplier maps are at
1.002 and 0.999, sync is at E029's level plus the predicted +0.013; the fainter dust quartiles
sit up to +0.03 above beta_d where the GMRF prior interpolates, and an independent-region
shrinkage prediction reproduces the pattern (+0.035, +0.042, +0.010 predicted against +0.033,
+0.031, +0.003 observed in Q2 to Q4)** (E028 stage A', `gmrf16_exact_mumap_nc_Tmap_mu01`:
E029's chain with the nside-16 truth T map in the dust SED; 5.5 h, 700 post-burn draws, 31
gradients per iteration, 0 divergences; D025 passes, KS p 0.68 and 0.69). Dust pixel RMSE
against the PySM truth 0.038 against E029's 0.056 (floor 0.019); the GLS CMB residual 93.08 uK
(E029 94.52, floor 93.18; not a criterion). Pre-registered criteria: sync within 0.03 of E029
pass, bright dust z within 0.5 pass, mu_d map 1.00 +- 0.01 pass, "every dust quartile within
0.03 of beta_d" missed by 0.003 and 0.001 in Q2 and Q3. The miss is the field prior's: the true
quartile means are 1.570, 1.523, 1.516, 1.548 while the field mean is 1.558 +- 0.012 and the
region Fisher sd in Q2 to Q4 (0.27 to 2.1) is far above the field sigma (0.049), so those
regions are the GMRF's stationary interpolation; the exact marginal of the same model has its
optimum at zero in every quartile (stage A scan, sd 0.1). The scalar-T chain (F35) agrees
with this one to 0.014 per quartile. Consequence: on this sky the T model plus the tight
hyperprior is the unbiased configuration (the T model is what moves the dust RMSE from 0.056
to 0.038 against E029), what remains is prior-dominated interpolation at the field mean, and
that is a Paper A reporting item (D025), not a method defect; a per-quartile tolerance of 0.03
is not a usable criterion where the quartile's true mean sits 0.04 from the field mean. Scope:
one sky, one noise realization, exact templates, oracle T map.

**F37. The hyperprior width of the per-region template multipliers must be at least the
template error's scale: with a 30 percent region-scale template error the tight N(1, 0.1^2)
that is unbiased with exact templates (F33, F36) biases the data-dominated sync quartile by
+0.21 (9 curvature sd, 38 nats) and the data-dominated dust quartile by -0.020 (6 sd, 24 nats)
in the exact expected target, while the matched width 0.3 leaves the bright dust quartile at
zero and the bright sync quartile at +0.07 (2.5 sd, 3 nats; +0.04 of it the template, +0.03
the Occam push) and costs at most +0.04 (dust) and +0.05 (sync) when the template is exact**
(E030 pass 3, exact scans of the collapsed hier target of Commander4's kernel with the oracle
nside-16 T map, multipliers integrated analytically, scatters profiled, per-quartile shifts
over +-0.3, CMB prior sd 1e4; `E030/scan3_sd{01,03,05}.txt`; the exact-template cells at sd 0.1
and 0.5 reproduce the E028 stage A scans exactly). Optima per dust quartile Q1 to Q4 and sync
Q1, s_all: exact sd 0.3: +0.01, +0.04, +0.02, +0.02; +0.03, +0.03. cs sd 0.1: -0.020, +0.03,
+0.08, +0.12; +0.21, +0.25. cs sd 0.3: 0.000, +0.03, +0.08, +0.13; +0.07, +0.12. cs sd 0.5:
+0.01, +0.10, +0.13, +0.19; +0.07, +0.13. The faint sync quartiles run beyond +0.3 with the cs
template at every width (exact: +0.03 to +0.05 at sd 0.3, +0.15 to +0.20 at sd 0.5): the
template error's own effect dominates the Occam push in the prior-dominated half (qualifies
F28), where the GMRF index prior sets the field in a chain (E026's cs_mumap at sd 0.5 with the
fixed T: sync RMSE 0.335, field scale 0.35). The misspecified tight prior pushes the template
error into the profiled scatters (sync 0.09 to 0.35 against 0.04 to 0.21). Consequence: the
hyperprior width is a property of the template, not of the sky; on real data it follows the
template quality (sd 0.1 only for a template good to about 10 percent at the region scale), and
at the matched width a data-dominated sync index carries a template-induced shift of order
+0.07 for a 30 percent template (the chain `gmrf16_cs_mumap_nc_Tmap_mu03` confirms it with the
field prior: the hundred brightest regions at -0.01 to +0.09, the Fisher-weighted level +0.08
to +0.10; the region map as a whole follows the displaced field mean, F38). Pass 4 scans
(2026-11-09; chain 2 at sd 0.5 read 12:20: the bright dust quartile is +0.036 there, of which
+0.018 is shrinkage toward the field mean through the smoother field, sigma_d 0.08, so "leaves
bright dust at zero" holds at the matched width and not above it; the hundred brightest sync
regions sit at +0.06 at sd 0.5 against +0.09 at sd 0.3, both at the scan's +0.07 within the
region scatter): over three realizations of the 30 percent error (seed 11, its negative, seed 12)
the data-dominated sync shift at the matched width is +0.07, -0.06, +0.01, a random shift with
scale 0.06 and mean near zero, so the "+0.07" is one draw and not a bias; the faint push stays
positive in all three (sync +0.2 to +0.3, dust +0.03 to +0.14; even in the error field's sign);
and a near-flat multiplier prior (sd 3) pushes every prior-dominated quartile beyond +0.3 with
40 to 136 nats for exact and cs templates alike (the F29 Occam push at full strength), so the
width must be matched and not wider. Also
recorded: two earlier passes ran with the scan script's default per-pixel CMB prior sd 100 and
showed +0.04 to +0.12 in the faint dust quartiles with exact templates and the true T, the F26
leakage at realistic noise; invalid as readings of the target, kept as the F26 illustration.
Scope: one sky, one template-error realization (seed 11, 10 degrees, log-normal rms 0.3),
expected target without the index prior.

**F38. With a 30 percent region-scale template error at realistic noise the index-field mean
hyperparameter is set by the prior-dominated faint majority of regions, and the GMRF prior
propagates its displacement into every region fainter than the field's own sd; only the
data-dominated regions and the Fisher-weighted level stay where the exact target puts them**
(E030 chain 1, `gmrf16_cs_mumap_nc_Tmap_mu03`: cs templates, the matched hyperprior sd 0.3,
oracle nside-16 T map, 1200 iterations, 0 divergences; D025 passes with KS p 0.93 sync and 0.40
dust, so the chain samples its target; `E030/chain_*.txt`). Sync field mean -2.701 +- 0.050
against the truth -3.000, 6 posterior sd off, where the template-error-dominated faint half
sits (scan optima beyond +0.3 at every width, F37); E026's cs chain at sd 0.5 with the fixed T
had -2.723 +- 0.055; with exact templates the means are -2.967 (stage A') and -2.981 (E029).
Dust field mean 1.569 against the `bd_trueT` mean 1.539, smaller because three dust quartiles
carry data weights 0.89, 0.29, 0.07. The faint sync quartiles sit +0.31, +0.37, +0.32 above the
region-fit reference, at the independent-region shrinkage prediction (+0.29, +0.30, +0.29;
`quartile_shrinkage.py`), and so does the bright quartile's average +0.19 (prediction +0.21: the
median Q1 region has Fisher sd 1.00 against the learned sigma_s 0.51, data weight 0.21). The
ten, thirty and hundred brightest regions by Fisher weight (31, 50, 67 percent of the weight)
sit at -0.01, +0.06, +0.09 and the Fisher-weighted level of the field at +0.08 (truth) to +0.10
(reference), where the scans put a data-dominated index at this width (+0.07, F37). Bright dust
is at the reference (+0.004, z median +0.17); dust Q2 to Q4 +0.032, +0.059, +0.023 against the
prediction +0.033, +0.049, +0.021. The multiplier maps track the realized template error in the
bright quartile (corr +0.86 dust, +0.66 sync) and not beyond (Q2 +0.44, +0.37), their mean
1.009 and 0.982 against exp(-f)'s 1.046. Pixel RMSE 0.426 sync, 0.117 dust against stage A' 's
0.066, 0.038 (exact templates and sd 0.1, so template and width differ at once; label corrected
2026-11-09, the text said stage A; paired +0.360 +- 0.021, +0.064 +- 0.005) and against E026's
cs chain 0.335, 0.089; the field is rougher (sigma_s 0.51 against 0.35, sigma_d 0.17 against
0.12), and the E026 comparison is confounded by the T map, the width and E026's hyperparameter
ESS of 4 to 13 (chain 2, cs at sd 0.5 with the oracle T, separates them). GLS CMB residual at
the posterior-mean fields 96.4 uK against stage A' 's 93.1 (the difference is the noise term, F40).
Consequence: with a real-data-quality template the reportable index quantities are the
data-dominated regions and the Fisher-weighted level, not the field mean or the region map as a
whole; the field-mean hyperprior, or the multiplier prior's law, must be changed so that the
template-dominated regions do not set the field mean (Q18). Scope: one sky, one template-error
realization, the oracle T; the sd 0.5 chain is pending.

**Amended (2026-11-09 12:20, E030 chain 2, `gmrf16_cs_mumap_nc_Tmap_mu05`: cs, sd 0.5, oracle T;
D025 passes, KS p 0.78 and 0.72; 20 divergences):** the width drives the roughness. sigma_s 0.32
against 0.51 at sd 0.3, sigma_d 0.081 against 0.174; pixel RMSE 0.327 sync, 0.091 dust against
0.426, 0.117 at sd 0.3 and E026's 0.335, 0.089, so the T map changes the RMSE by under 0.01 and
E026's cs numbers are a converged chain's. GLS CMB residual 94.45 uK at the draws against 97.66
(floor 93.18). The sync field mean is -2.722 +- 0.051 at sd 0.5 as at sd 0.3, the faint sync
half +0.26 to +0.34 at both, the hundred brightest sync regions +0.06 against +0.09 (scan
+0.07), so the F38 displacement does not depend on the width. The width trades coverage for
RMSE: at sd 0.5 the region map is one posterior sd high almost everywhere (z median against
the reference +1.03 sync, +1.00 dust; |z| < 1 for 47 and 49 percent of the regions against 69
and 65 at sd 0.3), the dust field mean +0.073 above the `bd_trueT` mean against +0.030, and the
bright dust quartile +0.036 (z +0.73) against +0.004, of which +0.018 is shrinkage toward the
field mean (data weight 0.64 at sigma_d 0.08 against 0.89 at 0.17). The multiplier maps absorb
the same share of the error at either width (rms against exp(-f) 0.28 to 0.30, the faint
half's own error rms; the width sets only their posterior sd there, 0.46 to 0.49 against 0.29).
Neither width gives a calibrated region map under a 30 percent template; what holds at both is
the data-dominated regions and the Fisher-weighted level (+0.07 to +0.09 against truth), D027.
Scope unchanged: one sky, one error realization, the oracle T.
**Amended (2026-11-09 evening, E030 chain 3, `gmrf16_cs_mumap_nc_Tmap_mu03_ref23`: chain 1 with
the amplitude and template reference at the leading band, 23 and 353 GHz; 1200 iterations, 0
divergences; D025 passes, KS p 0.72 and 0.43):** the displacement changes sign and the remedy
fails. Sync field mean -3.164 +- 0.049 against the truth -3.000 (3.4 posterior sd low;
pre-registered within 0.10), sigma_s 0.121. The field is a near-constant at the displaced mean:
every quartile sits 0.15 to 0.18 below the true-T region fit (all sky -0.164 against the
independent-region shrinkage prediction -0.162, rms 0.024), with |z| < 1 for 27 to 33 percent of
the regions. The hundred brightest regions by Fisher weight sit at -0.126, where chains 1 and 2
had +0.09 and +0.06, so "the data-dominated regions and the Fisher-weighted level stay where the
exact target puts them" does not hold across reference frequencies and is withdrawn as a general
statement. The sync pixel RMSE of 0.177 equals that of a constant at -3.164 (0.176) and carries
no information. Dust: field mean 1.591, sigma_d 0.110, RMSE 0.084. The pre-registration rested on
a misreading of the pass 7 scan, whose faint-quartile optima were -0.18 to -0.30 and not flat;
the four quartile parabolas sum to -0.06 for the level and -0.23 +- 0.09 for the faint three.
What stands of F38 is the description: the field mean follows the prior-dominated majority, the
GMRF carries it into every region, and the regions follow the shrinkage prediction to rms 0.024.
What does not stand is the causal reading through the lever arm and the remedy of referencing the
template at the leading band. Four chains with a 30 percent template put the sync field mean at
-2.72, -2.70, -2.72 and -3.16 against -3.00. At noise scale 1 the hier prior with a
component-separation-quality template has no configuration with a calibrated field mean. Whether
that matters at real intensity depth is E032.

## 8. Findings of the 2026-11-09 review

**F39. Noise scale 1 of the PySM ladder is 24 to 90 times real WMAP/LFI intensity white noise (88
to 1090 times HFI, 17 to 64 times WMAP/LFI Q/U), and scale 0.03 is 0.7 to 2.7 times the real
WMAP/LFI depth. Real intensity synchrotron is limited by the component model and the 408 MHz
anchor, not by noise; real polarized synchrotron is where the weak-data machinery applies**
(E031: per-pixel white noise of the ladder against the real maps at nside 64; Fisher sd of region
indices at the truth with the amplitudes profiled). Median synchrotron index sd per nside-16
region: ladder scale 1 5.99, scale 0.03 0.180, real three-component intensity 0.078, real
five-component intensity 1.36, with 408 MHz 0.10; real dust 0.012. Polarized synchrotron per
nside-16 region and Stokes parameter: 2.18 with the BP XV band set, 0.77 with K band added (0.37,
0.18, 0.087, 0.035 at nside 8, 4, 2, 1), 0.17 with an 11 GHz channel at an assumed depth. The
truth's index sd on this sky is 0.063 (sync) and 0.037 (dust). The median synchrotron S/N per
nside-64 pixel at 23 GHz is 0.31 at scale 1 and 16 at real depth. Consequences: the label
"realistic noise" for scale 1 is withdrawn in every record; F26 to F38 describe a weak-data
regime whose real counterpart is polarization, where none of them was obtained; at real
intensity depth the binding terms are the component model (free-free, spinning dust), the zero
levels (F17) and the ingest path (F41). Scope: Fisher numbers at the truth with white noise; real
low-multipole polarization noise is correlated and worse.

**F40. The total flat-prior GLS CMB residual at noise scale 1 is the noise term and ranks no
model. Its bias part is a Galactic-plane number, and at high latitude no chain with a
component-separation-quality template beats a constant index** (E033: the residual of the GLS CMB
estimate at assumed SEDs split analytically into bias, w^T sky minus cmb, and noise; the weights
reproduce `gls_cmb` to 4e-12 uK). The total is 90.9 to 96.4 uK for every model against a noise
term of 90.5 to 96.0. Bias over the sky and at abs(b) > 20, in uK: the truth's nside-16 region
means 5.95 and 0.42; a constant at the truth means 16.3 and 0.95; the exact-template chain (stage
A') 5.60 and 0.84; the cs-template chains 6.7 to 9.8 and 1.56 to 3.34. Every chain beats the
constant over the full sky by 30 to 60 percent, through dust and the plane. At abs(b) > 20 the
constant is matched by the exact-template chain and beaten by no cs-template chain, and at l = 2
to 10 chain 1 has six times the constant's bias amplitude (1.92 against 0.32 uK). Withdrawn with
this: every ranking of models by the total residual in E026 to E030 and Paper A Table 3, "chain 3
improves the CMB residual", and the E029 attribution. Consequence: in intensity the stake of the
index model for the CMB is 0.3 to 0.5 uK per low-multipole band-power bin against a signal near
30 uK, so index fields in intensity are a foreground-science gain first. The same bias in
polarization would be compared with a signal of 0.1 to 1 uK and has not been measured. Scope: one
sky, pseudo-C_l over f_sky, band powers exploratory.

**F41. Two steps of the map ingest make the real-map error model wrong: the model-error floor
reaches the likelihood 27.6 times below its nominal value, and maps are degraded before they are
smoothed, which leaves band-dependent CMB aliasing of 0.6 to 1.7 uK per nside-64 pixel** (E034;
raised by the pipeline audit as B14 and B1, expectations written before the runs).
`add_signal_fraction_to_rms` adds a white term per native pixel, which the degradation (8 times)
and the smoothing (3.45 times) average down, so the E014 chains had 100, 353 and 857 GHz at 0.37,
0.49 and 0.83 uK_RJ per pixel, the tightest bands, where the E014 table lists 5.1, 13.4 and 23.
That table's rms column is the degraded rms with the nominal floor and without the smoothing
reduction; it reproduces to 2 percent. In the chain's own noise model the high-latitude reduced
chi-squares are of order 13 to 60 in every band, where the table shows 0.0 to 6.4. For a CMB-only
simulated sky the code path (`ud_grade`, then smoothing to 120 arcmin) leaves differences between
bands of 0.79 uK_CMB (13 against 32 arcmin native beams), 1.5 (13 against 60) and 0.64 (32
against 53), and the other order leaves 0.0000. Against the per-pixel noise that reaches the
likelihood these are 0.65 to 3 times the noise for the CMB-dominated bands. What stands of E014:
the chains agree with the exact targets of the model as implemented (F2), the ratio of F14, and
the mean shifts of F17. What does not: any number in posterior sd, the reading of which bands
carry the fit, and the dust index as a measurement. Upstream is affected by both items. The
reader fix is queued behind the E032 chains. Scope: the aliasing test is CMB only and is a lower
bound; the effect on the indices needs the chains rerun.

**F42. Commander4 computes two chains in every run and writes one. Every chain of E020 to E030
had an unwritten twin, and no between-chain diagnostic has been computed for any result** (code
reading 2026-11-09, audit G9 and G7: `gibbs_schedule` in `cli.py` walks `chain in (1, 2)` inside
every iteration, the random stream is keyed on seed, chain, iteration and rank, the NUTS state is
kept per chain and group in `compsep/index_field_sampler.py`, and `file_io/chain_writer.py`
writes only the chains in `output.chains.write`, `[1,]` in every parameter file used so far). The
two chains share their start and differ in random stream. Half of the compute of every run was
discarded, and wall times per written iteration cover both chains. From E032 on both are written.
Consequence: split R-hat over at least four chains from two starts becomes a standard check
(VALIDATION.md), at no extra compute.

**F43. `get_sky` evaluates a component that carries an index field at the field's mean index**
(code reading 2026-11-09, audit B3: `sky/diffuse_components.py`, `get_sky` multiplies the
amplitude map by the scalar `get_sed(nu)` at the mean that `set_beta_map` stores). The
index-field sampler and both amplitude solvers use the per-pixel SED, so sampled posteriors are
not affected. Affected: Commander4's own chi-square in every index-field run (the E020 and E021
records attribute its value of 50 to 100 to the alm representation alone, which is incomplete),
the residual passed to a group when an index-field component is held fixed, the `conditional`
likelihood, and the sky handed to the TOD side. No recorded result uses any of the four. The fix
is queued behind the E032 chains.

**F44. For signed Q/U amplitudes at real WMAP/LFI depth the profile likelihood and a zero-mean
amplitude prior with one free scale per index region give calibrated region indices without any
template, while the flat-prior marginal is shifted by ln(nu_ref / nu_pivot) times the index
variance of a single amplitude, 0.06 to 0.9 with a sign and size set by the reference frequency
and by the CMB prior** (E036: exact grid posteriors for one index per region, no sampler; truth
amplitudes from the DR1 polarized synchrotron map at 5 deg and nside 16, one index over the sky,
real Q/U white noise of WMAP K, Ka, Q and LFI 30, 44; 100 noise realizations per cell; predictions
written before each stage). At real depth with synchrotron only (median S/N per element 4.45) the
profile and the zero-mean prior have mean z within 0.11, sd of z 0.97 to 1.10 and 68 percent
coverage 0.64 to 0.72, for 48 regions of 128 amplitudes and for 192 regions of 32, down to region
index errors of 0.4. The flat-prior marginal sits at -0.061 (reference 23 GHz), +0.062 (30 GHz)
and +0.34 (60 GHz) for an index error of 0.06, and at -0.40 and +0.34 for 0.16; the formula gives
0.118 for the difference between 23 and 30 GHz against the measured 0.122, and the shift does not
average down with the number of amplitudes that share the index. Noise ladder: at twice and four
times the noise (S/N 2.2 and 1.1) the profile stays calibrated (sd of z 1.00 to 1.09) except in
the faintest class of the 30 degree regions (1.27 and 1.40, coverage 0.60 and 0.58), and the
free-scale prior shows mean z -0.30 in one class at four times the noise, which a scale map
removes even when it carries a 30 percent error. With a free CMB and a free dust amplitude per
element under flat priors the index error grows by 1.63, the profile and the zero-mean prior stay
calibrated (mean z within 0.08, sd of z 0.91 to 1.08), and the pivot of the flat-prior shift
moves to the lowest band: +0.006 at 23 GHz and +0.17 at 30 GHz for an index error of 0.064, +0.86
for 0.13. With a Gaussian CMB prior of 0.2 uK the index error returns to the synchrotron-only
value (0.059 against 0.060) and the pivot to 26.4 GHz. With the truth moved off the grid centre
(-2.8, -3.4) nothing changes below a region index error of 0.25; above it the posterior mean
reverts 13 percent toward the centre. Three of my predictions failed, all in the same direction:
I expected the profile to under-cover at real depth, at twice and at four times the noise, and it
did so only in one class. Consequences: the flat-prior marginal, which is what a Gibbs sampler
with flat amplitude priors targets (F15), is not a usable target for a polarized index at any
reference except by accident of the pivot; E022's +0.27 in the Spur at a 30 GHz reference with a
free CMB amplitude is this effect on real data, and its recorded numbers obey the formula: for
two noise models the Spur shifts of +0.271 and +0.172 give the same pivot, 26.04 GHz, to 0.2
percent, while with a 0.19 uK CMB prior the pivot is 29.5 to 30.2 GHz, BeyondPlanck's reference,
which is a candidate explanation of why its chains show no flat-prior pull (Q17, the open point
of F22); the amplitude prior of the polarized track is
the profile with a sandwich width or the per-region zero-mean prior, with a CMB prior the size of
the signal (D029); no template is needed, so the template-error findings F21 to F38 do not arise
there. Scope: white noise, one index over the sky, a truth amplitude field that is a posterior
mean of real data, no dense covariance. The profile's interval is too narrow by the sandwich
factor sqrt(1 + sum V / sum u V), with u the amplitude S/N squared and V the weighted variance of
ln(nu) per element; computed from the truth alone it matches the observed sd of z within 0.05 in
15 of 16 classes below a region index error of 0.25 and gives 1.16 and 1.32 against the observed
1.27 and 1.40 in the faint 30 degree regions at two and four times the noise.
A zero-mean prior with one scale shared by all regions on a scale map that is 30 percent high at
high latitude moves the high-latitude index by -0.026 at a 23 GHz reference, 0.000 at the pivot
and +0.029 at 30 GHz (30 degree regions; half that for 15 degree regions): the sign structure the
independent challenge derived, small at this depth for regions the data inform.

**F45. At WMAP and LFI white-noise depth the collapsed hier prior with a 30 percent template error
shows no field-mean displacement, but its index fields do not beat a constant index, and its
field hyperparameters fail a four-chain convergence criterion that every chain's own check
passes** (E032: the E030 chain 1 model at noise scale 0.03, two runs with both in-flight chains
written, seeds 7 and 8 and dispersed starts, 700 post-burn draws per chain; eleven pre-registered
items plus four for the chain pairs). Sync field mean -3.016 to -3.018 against -3.000 (+0.30 and
-0.16 at scale 1, F38), dust +0.005 to +0.009 from the truth mean, every quartile within 0.03 of
the region-fit reference, multiplier maps at correlation 0.91 to 1.00 with the realized template
error (0.25 to 0.72 at scale 1), CMB bias 4.55 uK over the sky and 0.76 to 0.79 at abs(b) > 20
against 16.3 and 0.95 for a constant index. Against that: pixel RMSE 0.065 to 0.067 (sync) and
0.059 to 0.062 (dust) against 0.0629 and 0.0374 for a constant at the truth mean; inferred field
amplitudes 0.093 to 0.096 and 0.067 to 0.079 against truth scatters of 0.062 and 0.032; dust
regions inside one posterior sd in 34 to 49 percent of the faint half. Convergence: all four
chains pass the D025 conditional check (KS p 0.09 to 0.81), yet split R-hat between the first
chains of the two runs is 1.11 (sync field mean), 1.34 (dust field mean) and 1.11 (dust log
sigma), and over the eight half-chains 1.05, 1.18, 1.42 and 1.41 (dust range). Chains with the
same start differ as much as chains with different starts, so this is slow mixing of the
hyperparameters and their within-chain ESS (10 to 130) is overestimated. Region indices,
amplitude scatters, multiplier maps and the CMB bias agree across the four chains. Consequences:
F29 to F38 are weak-data effects of a regime with no real intensity counterpart; the model has no
advantage over a constant index at real depth on this sky, which answers the referee's objection
7; every single-chain statement about a field mean, amplitude or range (F29 to F38) carries a
Monte Carlo error of about one posterior sd; D025 is necessary and not sufficient. Scope: one
sky, one noise and one template-error realization, oracle T map; for the HFI bands this noise is
2.6 to 33 times the real one, so the dust numbers do not describe real dust data.

**F46. For a profile-likelihood index the curvature width is too narrow by sqrt(1 + T / J), the
sandwich evaluated at the fitted amplitudes is narrower still, and the debiased form J(a_hat) - T
gives nominal coverage from one data set; the zero-mean prior with a free scale is the same
estimator under a common noise pattern and a biased one otherwise** (E041: one index for 32 to
8192 heavy-tailed amplitudes, five bands, diagonal noise with a per-amplitude pattern, per-amplitude
S/N from 9 to 0.07, 400 realizations per cell, six predictions written before the run; unit test
in `nextgen/lowres/tests`). J = sum u V is the profiled Fisher information and T = sum V the
noise-noise term of the score variance (general forms in `lowres/likelihood.py`). (a) With the
curvature width the sd of z is sqrt(1 + T / J) within 0.054 in all 30 measured cells and follows
it to 2.6 at T / J = 6. (b) J evaluated at the fitted amplitudes has expectation J + T, so the
sandwich built from it is too narrow by up to a factor 5 in sd of z: worse than no correction.
(c) With J_c = J(a_hat) - T and variance J(a_hat) / J_c^2 the sd of z is 0.946 to 1.052 and the
68 percent coverage 0.650 to 0.718 in 29 of 30 measured cells (0.635 in one), for T / J up to
1.6. (d) The flat-prior shift equals sum_p ln(nu_ref / nu_p) / J within 4 percent at a 23 GHz
reference and 15 percent at 30 GHz up to shifts of 0.3, and is below 1e-3 at the pivot reference.
(e) With a noise pattern common to the amplitudes the zero-mean prior with a maximized scale has
the profile optimum exactly and the debiased sandwich as its curvature variance (algebra,
confirmed to 0.008 in sd of z). With a pattern that varies between amplitudes it is biased by
-0.6 to +0.3 posterior sd depending on the amplitude draw, because a Gaussian prior is
misspecified for heavy-tailed amplitudes and the bright ones carry whatever pattern they happen
to have. This qualifies F44: its "calibrated within 0.11 except -0.30 in one cell" was this
effect. The estimator of record is the profile with the debiased sandwich; the zero-mean prior is
reported beside it. Scope: one region, one component, diagonal noise; the coupled dense case is
E039.

**F47. Two input facts for the polarized track: the white-noise columns of the PR4 LFI map files
understate the noise of an nside-16 pixel average at every scale, and the 353 GHz band cannot be
treated as a delta band** (E038 b and b2; E040 amendment). (a) The ring-half difference of the
PR4 30, 44 and 70 GHz Q/U maps, averaged to nside 16 with inverse-covariance weights, has a
chi-square per element of 1.77, 1.40 and 1.21 against the covariance columns, with a ratio of
1.67, 1.36 and 1.15 at multipoles 30 to 47 and 3 to 7 at multipoles 2 to 6 (two of three
expectations failed: I expected the pixel scale to be white). A three-parameter model per band,
white level times r plus a power law in multipole in the whitened frame, fitted to the half
difference, passes 22 sub-tests per band with none beyond 3 sd; the 70 GHz model is the weakest.
It is isotropic, fitted to one realization, and blind to errors common to the halves. (b) A
modified blackbody with index 1.53 and 19.6 K seen through the HFI 353 GHz band is 1.137 (all
bolometers) or 1.121 (PSBs) times its value for a delta band at 353.0 GHz; the bandpass reader
reproduces the published unit conversions to 0.3 percent and the LFI colour-correction ratios to
3e-4. A delta band at 353 GHz therefore overstates the dust amplitude at every lower frequency
by 12 to 14 percent, more than a dust-index change of 0.05 does at V band. (c) WMAP K band and
PR4 30 GHz agree in sign and units (correlation 0.985 in Q, 0.970 in U on elements above 10 uK),
and their unweighted ratio is 0.453 and 0.470 against 0.502 for an index of -3.1 with
band-averaged spectra. That is not a measurement: an unweighted fit gives the poorly measured
large-scale modes of the K-band map full weight.

## 9. Real polarization data (2026-11-10 and 11)

**F48. On the WMAP nine-year K to V polarization maps the flat-prior reference shift is there at
the predicted size, and the profile likelihood gives indices near -3.1 with errors of 0.01 to
0.02** (E040, four BP XV regions, nside 16, released dense covariances, 353 GHz for dust through
its measured bandpass, CMB prior 0.2 uK, bands brought to a common window; expectations written
before the first fit; E039 stage 1 passed first). Profile optimum with debiased sandwich error,
all scales: High Latitude -3.094 +- 0.012, Spur -3.097 +- 0.015, Galactic Center -3.087 +- 0.009,
Plane -2.995 +- 0.012; chi-square per element and band 0.58, 0.59, 1.82, 0.78 against an expected
0.6. The same likelihood under flat amplitude priors: -3.102 (reference 23 GHz) and -2.944
(reference 30 GHz) at high latitude, a move of +0.158 against a first-order prediction of +0.172
(Spur +0.075 against +0.079, Galactic Center +0.025 against +0.027, Plane +0.037 against +0.038);
12 sandwich sd at high latitude. At nside 8 the move was +0.062. The zero-mean prior is within
0.4 sd of the profile in the three regions that fit. Grade A for the shift (an identity,
confirmed on simulations with the real covariances in E039 and on the maps), B for the index
values, which carry F49.

**F49. The pixel scale of the nside-16 maps carries a beam-aliasing bias toward a flatter index,
0.02 to 0.07 for WMAP K to V** (E042 negative control and low-pass arm; E038 c2). Removing
multipoles 2 to 20 from one band moved the K, Ka index by +0.08, +0.16, +0.11 and +0.37 where the
pre-registered control required no movement. Every amplitude ratio to K band rises by 6 to 10
percent at multipoles 35 to 47 (Ka over K: 0.32 to 0.33 up to multipole 34, 0.35 above). The maps
are pixel averages of skies that are not band-limited and K band has the widest beam; the window
operator of E038 works to multipole 47 only, and its validating simulation had no compact
structure. With the lowest band and 353 GHz at all scales and the other bands informing the index
through multipoles 2 to 20 only, K, Ka, Q, V gives -3.113 +- 0.014, -3.112 +- 0.017, -3.140 +-
0.011 and -3.061 +- 0.014: steeper by 0.019, 0.015, 0.053 and 0.066, errors 10 to 20 percent
larger. A separate failure of the same family: degrading to nside 8 with the dense weights biases
the Galactic Center by +0.22 in simulation and +0.23 on the maps (E040 R8). Grade B: one cut
value read so far; stability across L = 12, 20, 30 pending.

**F50. Planck 30 GHz polarization cannot be combined with WMAP K to V in one power law at this
precision, and the answer depends on the release** (E040 stage 2, E042). Index from multipoles 2
to 20: K with PR4 30 GHz -3.590 +- 0.023, -3.244 +- 0.028, -3.553 +- 0.017, -3.217 +- 0.025
against K with Ka -3.092, -3.085, -3.109, -3.042; 30 GHz with Ka -2.31 to -2.79. PR4 30 GHz is
fainter than the extrapolation of K band by 4 to 11 percent. With the PR3 map (all scales) the
K to 30 GHz index is -2.921 at high latitude where PR4 gives -3.428, with a reduced chi-square of
6. Observed minus predicted for PR4 30 GHz from a WMAP-only fit: 16 to 43 times the noise model
at multipoles 2 to 6 in B as in E, against 2 to 4 for the Ka control. Projecting out the 20
K-band modes with more than ten times the median variance moves the high-latitude index by 0.03.
Exploratory (found after looking at the residual map): for abs(b) < 5 in the inner Galaxy PR4 is
15 percent fainter than PR3 (slope 0.854), 4 percent for 5 to 20 degrees, equal above 20, and the
plane residual follows the intensity at 0.05 to 0.13 percent, the size of a difference between
two bandpass-leakage corrections. Grade B for the band-pair indices, C for the leakage reading.

**F51. For the BeyondPlanck band set and reference (30 to 70 GHz, Ka, Q, V, 353; no K band; flat
priors on every amplitude; amplitude defined at 30 GHz) the flat-prior optimum is at the upper
bound of the allowed range at high latitude** (E040 stage 3, PR4 and WMAP maps, nside 16). Profile
-2.695 +- 0.040, -2.921, -2.819, -2.780; flat prior at 30 GHz -1.600 (bound), -2.314, -2.752,
-2.479; flat prior at 23 GHz -4.036, -3.432, -2.885, -2.997. Svalheim et al. (2023) and Watts et
al. (2023) report "unreasonably flat" high-latitude indices for a conditional index step at a 30
GHz reference and constrain them with N(-3.3, 0.1^2). The volume term has that sign and more than
that size in this configuration. Grade B for the numbers, D for any statement about those
chains, which were not examined. The profile index of this band set moves by 0.3 with the CMB
prior, which is F50 from another side.


**Amendments of 2026-11-11 to F48 to F51, after the independent referee read of Paper A1 draft 1
(`reviews/2026-11-11_referee_report_paperA1_draft1.md`).**

- F49: the heading says more than was shown. What is established is the size of the effect: with the
  index informed by multipoles up to 30, 20 and 12 only, the four-band index moves by -0.014 to -0.048
  (High Latitude), -0.008 to -0.015 (Spur), -0.043 to -0.093 (Galactic Center) and -0.043 to -0.086
  (Plane), and the pre-registered stability criterion (Z4) holds for the Spur, barely for the High
  Latitude region, and fails for the other two. "Beam aliasing" is one candidate cause; compact
  structure with a flatter spectrum is another, and Fuskeland et al. (2014) match our all-scale value in
  the plane, which favours the sky. The low-pass estimator had not been run on simulations when these
  numbers were written (E039 addendum B does that). The low-pass indices are sizes of an effect, not
  results, and Paper A1 quotes none of them as an index.
- F50: stands as a description of these maps under this model. The statement that the two releases
  differ "by 0.5" rests on PR3 with white-noise columns and a reduced chi-square of 6 and is not a
  measurement of a difference between releases.
- F51: the reading of the published analyses is reduced to what the numbers support: in this
  configuration the flat-prior optimum depends on the reference by more than a unit of index at high
  latitude. Whether that is what happened in any published chain is not claimed.
- F48: see F53 for the dependence of the values on the CMB prior.

**F52. On macOS 27.0 (build 26A5388g) Accelerate's Cholesky factorization with the upper-triangle
option writes up to 768 kB past an internal buffer for matrix orders of about 1500 to 6000, and
scipy's `cho_factor(lower=True)` and `cholesky(lower=True)` on C-ordered arrays take that path**
(E043). Measured with Apple's guard allocator (a byte write one past a heap block inside
`dpotrf$NEWLAPACK`, called with no scipy wrapper in between) and with a malloc canary
(`nextgen/tools/malloc_canary.c`): 786432 bytes past a 524288-byte block at order 1536, 352000 bytes
at 4608, none at 6144 to 18432, none ever with the lower-triangle option. The factor returned is
correct; the damage is to neighbouring heap memory, and it showed as crashes in unrelated code. Every
nside-8 dense fit and noise draw of E038 to E040 was exposed; the nside-16 fits were not (the baseline
rerun on the old path under the canary shows no overflow and reproduces every printed digit). Fix: the
other triangle, at the same cost. Reruns on the sound path under the canary reproduce the recorded
outputs as printed: plug-in addendum (40 of 40 realizations), nside-8 window checks, nside-8 real-data
arm; E039 stages 1 and 2 in progress at the time of writing (see the E043 record for the final count).
Grade A.

**F53. The profile index of the WMAP K to V fit depends on the shape of the CMB prior at the level
of its statistical error** (E040 stage 4). With the LCDM pixel covariance as the prior in place of a
white 0.2 uK prior of the same total variance, the four indices are -3.082 +- 0.012, -3.080 +- 0.016,
-3.081 +- 0.009 and -2.964 +- 0.013: flatter by 0.012, 0.017, 0.006 and 0.031 (0.9, 1.0, 0.7 and 2.3
sd). The pre-registered criterion (within 1 sd) failed, so the LCDM version is the central table of
the paper. The flat-prior reference shift holds in both versions (+0.167, +0.095, +0.029, +0.049 from
23 to 30 GHz with the LCDM prior, 5 to 9 percent below the first-order prediction). Whether the move
is what a LCDM sky does to a white-prior fit, or large-scale achromatic structure in the maps that the
LCDM prior absorbs, was tested on simulations (C4, 60 skies with a LCDM CMB at nside 8): with the matched
prior the estimate is unbiased and calibrated (sd of z 0.93 to 1.14); **with the white prior it is unbiased
on average and over-confident**, its quoted error too small by 32, 24 and 30 percent at High Latitude, in the
Spur and in the Plane (coverage of one sd 0.55, 0.48, 0.65), and it differs from the matched fit by a random
0.007, 0.009, 0.003 and 0.007 per sky; with a flat prior it is unbiased and calibrated at 2.2 to 2.5 times
the error. Against that scatter the moves on the maps are 1.6, 1.8 and 2.2 sd at High Latitude, in the Spur
and in the Galactic Center (jointly unremarkable: 22 percent of simulated skies further out) and 4.4 sd in
the Plane, which a LCDM sky does not do. With a flat CMB prior the map indices are -3.028 +- 0.029, -3.055 +-
0.036, -3.008 +- 0.020 and -2.988 +- 0.030, and the flat-amplitude-prior optimum of the same four WMAP bands
is at -2.51 (reference 23 GHz) and -1.71 (30 GHz) at High Latitude; the first-order prediction of the
difference holds within 10 percent for every region and prior, for shifts from 0.03 to 0.8. Grade B.

**F54. The four regions are not one index each: block errors are 2.3 to 7.5 times the model errors, and
half a dozen sky patches carry the information of each region** (E040 stage 5). The scatter of the profile
score over the 48 nside-2 pixels, with a leverage correction (first-order delete-one-patch jackknife), gives
errors of 0.044 (High Latitude), 0.038 (Spur) and 0.100 (Plane) with the LCDM prior, where the sandwich
errors are 0.012, 0.016 and 0.013; the effective number of patches is 5.7, 8.8, 2.5 (Galactic Center, no
block error quoted) and 5.1, and one patch carries 37 percent of the High Latitude information. On
simulations the estimator follows the true scatter with the model true (rms of z 0.82, 1.20, 1.02) and with
an index spread of 0.1 between patches (0.59, 0.99, 0.91), where the uncorrected form is too small by up
to 39 percent; the choice of form was made by a rule written before the tests were read. A spread of 0.1
reproduces the High Latitude and Spur ratios of the maps; the Plane needs about twice that. With the block
errors the four region indices (-3.08, -3.08, -3.08, -2.96) do not differ. Any misfit scatters the block
scores, so "a spread of the index" is one reading, not a measurement. Grade B.

**Amendment of 2026-11-11 08:50 to F48 to F54, after the second referee read of Paper A1
(`reviews/2026-11-11_referee_report_paperA1_draft3.md`).** Every region name in F48 to F54 refers to the E022
surrogate regions with no mask, and those are not what they are named:

- "High Latitude" contains the inner Galactic plane at 25 to 70 degrees from the centre on both sides (140 of
  2194 pixels within 10 degrees of the plane, holding 55 percent of the region's squared fitted amplitude;
  pixels above 30 degrees hold 27 percent). Its index (-3.08) is mostly an inner-plane index.
- Tau A is in "Plane", unmasked: fitted polarized amplitude 175 uK against a region median of 38; by the
  referee's first-order calculation its patch carries 58 percent of the Plane's block variance and moves the
  Plane index by 0.07 (from -2.96 to about -3.04).
- **F54 is withdrawn for High Latitude.** 65 percent of that region's block variance involves the Galactic
  Center's scores, through a 0.20 correlation of the two indices in J_c^-1; from the region's own scores the
  block error is 0.026 (2.1 sandwich errors), which 4 of 60 null simulations exceed. For the Plane the scatter
  is real and mostly one source. What stands of F54: the estimator, its tests, and the statement that half a
  dozen patches carry each region's information.
- F53: "a quarter to a third" is from a simulation at nside 8 with a prior variance 1.9 times the simulated
  one; relative to the matched prior on the same skies the understatement is 14 to 27 percent in three regions,
  with a Monte Carlo error of 10 points, and the independent review found the opposite sign in another setup.
  What stands: a white prior of the wrong shape miscalibrates the error by tens of percent, with a sign that
  depends on the configuration. E044 replaces the simulation by an exact first-order calculation.
- F48, F53: the estimator comparison is differential and stands as measured (flat prior at 23 GHz -3.092, at
  30 GHz -2.925, profile -3.082 on identical maps and regions), with two repairs: the change of convention is
  13.4 model sd (12.6 is the distance of the 30 GHz optimum from the profile optimum), and "predicted to 9
  percent" is a first-order step with the curvature of the fitted marginal; the prediction that needs no
  flat-prior fit (J_c) is off by 11 percent, and 18 percent with a free CMB.
- E044 redoes the real-data results on regions defined by geometry with a processing mask (D035).

**F55. On regions defined by geometry with a processing mask the simple model fits the WMAP K to V
polarization maps, and under flat amplitude priors the index still moves by 0.07 to 0.48 with the
reference frequency** (E044 R1; pre-registered). Regions: High Latitude (|b| > 30), Spur, Intermediate
(10 to 30), Plane (|b| < 10); mask: brightest 5 percent of the DR1 polarized synchrotron model and a 7.5
degree disc around Tau A; LCDM covariance as the CMB prior; nside 16, dense covariances, window matching.
Chi-square per element and band 0.577, 0.591, 0.572, 0.614 against 0.6 (the old regions without a mask
gave 0.58, 0.59, 1.78, 0.77). Profile index -3.147 +- 0.029, -3.139 +- 0.028, -3.091 +- 0.019,
-3.039 +- 0.015. Flat prior at 23 GHz: within one sandwich error of the profile optimum in every region
(pivot 23.3 GHz). Flat prior at 30 GHz minus at 23 GHz: +0.478, +0.154, +0.204, +0.072, which is 16.4,
5.6, 10.5 and 4.8 sandwich errors with no change of data or model. First-order predictions are always
high: measured over predicted 0.75, 0.84, 0.86, 0.93 with the curvature of the fitted marginal and 0.69,
0.79, 0.82, 0.90 with J_c. The index values are answers of one model on one set of regions and are not
quoted as a measurement (F58 will list what moves them). Grade A for the differences, B for the values.

**F56. The zero-mean amplitude prior with a free scale per region is biased steep at the configuration
of the map fit** (E044 R1 and `sim_n8_zm`; not expected). In 30 simulated realizations at nside 8 with
these regions, mask, covariances and prior, in which the profile optimum is unbiased (pooled mean z
-0.05, -0.20, +0.02, +0.03 over 120 realizations, se 0.09), the zero-mean optimum is steeper than the
profile optimum of the same realization by 1.05, 0.54, 0.44 and 0.22 sandwich errors, with a scatter of
0.02 or less. On the maps at nside 16 the difference is -0.104, -0.043, -0.030, -0.004 (3.6, 1.6, 1.5,
0.2 sandwich errors): same sign and order. F44 and F46 (e) had it as the same estimator for a common
noise pattern and up to 0.6 sd off otherwise; with three components, dense noise and a faint region it
is larger. Consequence: the zero-mean prior is not a substitute for the profile likelihood on these
maps. Grade B.

**F57. The misfit and the large block errors of the earlier real-data fits were the brightest emission
and one source; with the mask the block scores give no evidence that any region is more than one index**
(E044 R1, R2; pre-registered H4, H6). With only a 7.5 degree disc around Tau A masked on the old regions,
the old Plane index moves from -2.964 to -3.048 +- 0.014 (predicted -0.05 to -0.10; the referee's
first-order estimate -0.074) and its chi-square per element and band falls from 0.770 to 0.596. On the
new regions the own-score block error (CR3) over the sandwich error is 1.63, 1.00, 1.44, 1.28; under the
model (120 simulations at nside 8 with the regions, mask, covariances and prior of the fit) that ratio
has median 1.33, 1.30, 1.11, 1.12 and 97.5th percentile 2.02, 3.11, 1.69, 1.72; 15, 85, 12 and 31 of the
120 lie above the map values. H4 (ratios of 1.5 to 4 for Intermediate and Plane) failed on the low side.
This supersedes what was left of F54 for the regions. Grade A for the Tau A test, B for the reading.

**F58. On the masked regions a white CMB prior makes the index error 9 to 22 percent too small while quoting
a smaller error than the LCDM prior, and with a free CMB amplitude the flat-prior marginal has no optimum
at high latitude** (E044 R3, R5, R7; pre-registered H5, H11). White prior of 0.2 uK judged under the LCDM
covariance by the first-order formula of `lowres.mismatch` (tested against 1500 simulated fits): true over
quoted error 1.22, 1.21, 1.09, 1.10; bias below 0.1 of the quoted error; quoted error 0.91 to 0.97 of the LCDM
fit's, true error 1.05 to 1.10 of it. Flat CMB prior: profile index -3.041 +- 0.062, -3.122 +- 0.060,
-3.017 +- 0.044, -3.038 +- 0.035 (errors 2.1 to 2.4 times larger); with flat priors on all three amplitudes the
High Latitude optimum is at the upper bound (-1.6) for references of 23 and 30 GHz alike, and a scan shows the
marginal rising monotonically over the whole range -4.6 to -1.6 (gain of 463 and 1601 in ln L between the
profile optimum and the edge, where the profile likelihood loses 445). The index variance carried by one
amplitude is 11 there. Unexplained and reported: the index at High and Intermediate latitude steepens the more
tightly the CMB amplitude is held (white, LCDM, flat: -3.173, -3.147, -3.041 and -3.104, -3.091, -3.017), each
step 1.7 to 1.95 of its expected standard deviation (Q29). This replaces F53. Grade A for the calculation and
the scan, B for the truth model.

**F59. With the mask no test rejects the one-index model, and the dependence on angular scale is gone** (E044
R6; pre-registered H7, H8). Index from multipoles up to 30, 20, 12 only: at most 1.3, 1.4, 1.3 of the standard
deviation that the model predicts for the change (without the mask, on the old regions: 2.8 to 10.4 at l <= 20).
K and Ka alone: +0.6, -0.2, +2.4, +2.0. Model changes: dust index 0.05 away, 0.004 to 0.021 (most in the
Plane); delta bandpasses, -0.046 to -0.051 everywhere; drift correction 0.001; 353 GHz bandpass choice 0.003 or
less; floor of the CMB prior 0.002 or less. By-product: chi2_min at dust indices 1.48, 1.53, 1.58 puts the dust
index of these maps at 1.57 +- 0.02 (Planck: 1.53 +- 0.02). This supersedes F49 for these regions and answers
Q25 at the level of this experiment. Failed in the pre-registration: the sign of the Plane's change at large
scales (H7) and where the dust index matters most (H8). Grade B.

**F60. The chi-square of the masked fit against the released WMAP covariances is low, and not uniformly**
(E044 R4, R5; not pre-registered). chi2_min over degrees of freedom: 0.973 with the LCDM prior (2.5 sd), 0.956
with a flat CMB prior (3.4 sd; two residual degrees of freedom per element, exact count), 0.962 for K and Ka
alone (2.1 sd); a simulated null of the LCDM fit gives 1.002 +- 0.003 (40 realizations at nside 8). The deficit
is in the two faint regions. The one degree of freedom per element along the CMB spectrum (the difference of the
two nested fits) carries 1.008 +- 0.019. Read as noise it is 1.3 to 2.2 percent in rms, and every error of E044
would be too large by that. Literature (`research/2026-11-11_wmap_lowres_covariance_normalization.md`): no
statement of the accuracy of the r4 normalization, no claim that the matrices overstate the noise, two
uncommented published V-band values of 0.96 to 0.98 with the same matrices. Not rescaled. Q28. Grade B for the
deficit, C for any reading of it.

**Amendment of 2026-11-11 12:25: what E044 does to F48 to F54.** F48 (the reference shift on real maps) is
replaced in its numbers by F55 and stands in kind. F49 (scale dependence) described the unmasked regions and is
superseded by F59 for the masked ones: the effect was the masked emission. F53 (white CMB prior) is replaced by
F58. F54 (block errors) is superseded by F57. F50 and F51 (Planck 30 GHz, the band set without K band) were not
redone on the new regions and are Paper B material with the labels of the old regions.

