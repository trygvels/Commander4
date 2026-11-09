# Way-forward challenge: independent review of the proposed plan

- **Written:** 2026-11-09, reading from 23:03 CET, written from 23:22 CET, finished 23:33 CET.
- **Role:** independent reviewer (Bayesian component separation, Gibbs sampling, low-multipole
  polarization likelihoods), asked to find what is wrong or unsound in the plan, not to confirm it.
- **Repository state:** branch `nextgen`, read from 1eaac33 on. Three commits landed while I read
  (60ca7f0, 1c839d0, 63f811c: revised Paper A thesis, D028, new `VALIDATION.md` rules, Q19 to Q23)
  and E036 wrote its first output at 23:17. I read those too and say where the plan differs.
- **What I ran:** four numpy checks on one thread, about 30 s in total, piped to the interpreter
  so that no file was written (Appendix). No MPI, tests or builds. I ran `git log` and `git status`
  (read-only) three times before registering that the brief excludes git; nothing was changed.
- **Marks:** [abs] verified in the arXiv abstract today; [search] a web search summary today, not
  the paper; [rep] from a repository report whose numbers went through a summarizer; [mem] from
  memory, unverified. "My algebra" marks a derivation that no experiment record holds.

## Summary

1. Verdict: no-go on the plan as worded, go on a reduced one. As worded it yields two half-papers.
2. Paper A holds as (a) plus (b) with F21, F23 and a Q/U section. (c) and (d) lack an A or B base.
3. "Validated at real WMAP/LFI depth" is white noise on a three-component sky. E031 says so.
4. (a) is the known inexact MCWM scheme. The law is new and misses the real-map width by a quarter.
5. A Q/U prior with one global scale on a template is biased: +-0.12 in the faint sky, 3.5 sd.
6. One free scale per index region is sound and equals the profile with a sandwich error.
7. E036: profile and zero-mean prior agree at real depth. The prior that matters is the CMB's.
8. Dense NUTS is feasible at nside 8 (1 s per factorization), not at nside 16 (about a minute).
9. K band against 30 GHz is degenerate with the index. Test it first, on grids, before any field.
10. Instead: submit (a)+(b) in six weeks, run the K/30 null test, decide Paper B on the gates in F.

## A. The five elements: assumptions, support, cheapest falsifier

### A1. Paper A, narrowed

**(a) The legacy move and the width law 0.79 k^(1/4).**

1. *The derivation is right.* Supported: E035 (toy to 0.7 percent for k >= 1000, exact controls;
   E035:84-102) and the parallel stats review (`reviews/2026-11-09_stats_review_legacy_move.md`
   Sect. 1 and 2). I re-derived the kernel, its target and the law, and agree.
2. *It is new.* Not as framed. A chain that refreshes the estimate of the current state between
   accept steps is "Monte Carlo within Metropolis", named and shown inexact in Andrieu and Roberts
   2009 (arXiv 0903.5480) and said to give conservative, too-wide approximations [search]. That is
   F1. The stats review says so (its line 54); E035 and the thesis do not. New: the closed form.
3. *The law describes Commander4 on real maps.* Contradicted today: toy 13.1 to 13.6 against a
   measured 16.5 and 17.1 (E035:125-131, 155-161; Q23). Cheapest test, one thread, no chain: 20
   amplitude draws at fixed beta on the E014 inputs, chi-square with and without the alm round
   trip of F7. Its variance over 2k must be 2.4 +- 0.5 for that candidate to live.
4. *Small k is covered.* The first correction is T = sqrt(2 pi k)/4 + 1/2 (my algebra, and the
   stats review's alpha(0) = 1/2 + c_k). It explains the toy's excess at k = 10 and 100 (check 2).
5. *The move matters to a reader.* Unsupported. No published result is shown affected (referee
   objection 8), the Commander3 2024 module was never run (F10), there is no Commander4 paper
   ([rep]), and the full-resolution path is CG, where nothing is established (F1 scope, Q22).

**(b) Flat, profile and hierarchical amplitude priors compared by coverage.**

1. *The simulation results.* Supported: A for F12, F15, F18; B for F16.
2. *The thresholds are general.* Untested. "S/N about 15" belongs to one configuration
   (`FINDINGS.md:260-279`). The flat-prior shift is ln(nu_ref / nu_pivot) times the
   single-amplitude index variance (revised outline; E036: -0.061 at 23 GHz, +0.062 at 30 GHz), so
   every flat-prior number in F3, F12, F18 belongs to one reference. Test: E017 grids at three.
3. *The profile must be replaced below the threshold.* Weaker than claimed. For known noise and
   linear amplitudes the profile score is unbiased at any S/N: E[b^T A^-1 b] is the signal term
   plus the number of components, and the noise part does not depend on beta. What fails is the
   width, by the sandwich factor sqrt(1 + 1/u), u the mean amplitude S/N squared (my algebra).
   Check 4: sd z 1.09 against 1.14 at u = 3.35 and 1.04 against 1.02 at u = 30, mean z -0.06; at
   u = 0.38, where the index sd is 0.6, the asymptotics fail (2.44 against 1.91). F18's sd z of
   1.90 and 1.14 at S/N 23 and 68 fit one u scaled with the noise. A referee will ask why the
   paper gives an empirical threshold instead of this factor. Test: overlay it on E012 and E017.
4. *"Hierarchical" is the calibrated arm.* Only inside its model class (`FINDINGS.md:44`); F21, F23
   and F38 are what happens outside. The plan's list for (b) omits F21 and F23, and the revised
   thesis moves them to a "limits" item. They are the other half of (b).
5. *The SBC.* The first 300 replicates failed the pre-registered rule and the replication passed;
   one index was sampled with the dust SED fixed; E-BFMI was 0.03 to 0.17 (referee Sect. 2.3).

**(c) GMRF index fields "validated at real WMAP/LFI depth".**

1. *The GMRF beats exchangeable and independent priors.* B plus, on a truth that is region-constant
   and smooth at the model's own resolution, with exact templates, a fixed T, hyperparameter ESS of
   3 to 14 in 700 and intervals that over-cover at 78 / 98 (`FINDINGS.md:486-521`).
2. *Scale 0.03 is real depth.* Contradicted by E031. It matches WMAP and LFI white noise and is 2.6
   to 33 times above HFI. With free-free and spinning dust at known SEDs the real-noise region sd
   of beta_s goes from 0.078 to 1.36 (E031 Table 2), and the record concludes that the
   three-component sky "tests neither the noise nor the degeneracy structure of a real intensity
   analysis" (E031:250-254). Honest label: WMAP and LFI white noise on a three-component sky.
3. *E032 validates it.* One sky, one noise and one template-error realization, an oracle T map and
   a template that is the truth times a smooth error (E032:31-32, 41-42, 146-148): the benign law
   of F23, not the data-built law of F21. A pass earns B. A fail removes the only at-depth
   evidence, and the plan does not say what Paper A is then.
4. *The Commander4 implementation is validated at the reported configuration.* No SBC exists (D028
   now asks for one), and the L1 field test at template brightness is unresolved (log sigma -1.83
   against -1.31 to -1.52; referee Sect. 2.4).
5. Tests: the E032 follow-ups (PySM index maps with exact templates, a scalar baseline), and an
   exact-target scan with free-free and AME columns added. E031 predicts the regime fails it.

**(d) Real-map demonstrations that survive F41.** "The sampler agrees with its exact target" is
supported (F2), but on a model with a high-latitude reduced chi-square of 13 to 60 it is an L1
test with real maps as arbitrary inputs. The zero-level shifts (+0.035, -0.016) were measured
with the HFI bands at 235 to 760 times their intended weight (`FINDINGS.md:307-309, 362-363`);
means depend on weights, so they are untested after the fix. Test: the exact 2-d ridge grid on
the corrected reader output (Q19). The effect is published at 0.03 to 0.07 ([rep]: 2606.21334).

### A2. Freezing template laws, the dust temperature field and sampler efficiency

1. *None of them gates a paper.* Mostly supported (F39, F40).
2. *The freeze is consistent with A1(c).* It is not: E032 and its two follow-ups are intensity
   template-law experiments. A stop rule is missing. Mine: E032 is the last, pass or fail.
3. *It is consistent with element 5.* It is not. Hyperparameter ESS of 3 to 14 in 700 draws (F24)
   is an efficiency problem, and split R-hat below 1.05, let alone the 1.01 a referee expects
   ([mem]: Vehtari et al. 2021), needs longer chains or better mixing (E032:107-111 says so).

### A3. The polarized track

1. *It builds on E022.* Less than it sounds. E022's collapsed targets used the diagonal of the
   WMAP covariance; the dense N^-1 entered a conditional-only check (E022:80-82, 481-482). K band
   is not on disk (`~/data/wmap9_lowres/` holds Ka, Q, V). The WMAP maps were unsmoothed against 60
   arcmin LFI maps (E022:412-413). No dense LFI covariance is on disk; the PR4 one has a cosine
   window, a transfer function and an unread scaling step ([rep], science report lines 207, 287).
2. *A zero-mean prior with a global scale tied to a template is sound.* Unsound by my analysis (E).
   E036 tests a different model, one free scale per index region, and its first output
   (`E036/coverage_real_depth.txt`, read 23:20) says something the plan does not expect: at nside
   16 and real white-noise depth the profile and the zero-mean prior agree to 0.02 or better in
   shift, with sd z of 0.97 to 1.08 and 0.97 to 1.10, and the template adds nothing. The amplitude
   prior is not the problem at this depth. The reference frequency of a flat prior is (0.06 to 0.6).
3. *tau is a by-product.* The prior the index is sensitive to is the CMB's. In E022, with a flat
   per-pixel CMB prior the Spur marginal sits at -2.87 against a profile of -3.14 and High Latitude
   runs to the boundary; with a C(l)-sized prior of 0.19 uK the marginal is -3.197 against -3.193
   (E022:168-177). In a dense analysis the proper CMB prior is the signal covariance
   S(C_l^EE(tau)): tau is a hyperparameter of the one prior that matters, and comes for free.
4. *The data support an index field with inferred smoothness.* Weakly. With K band the profiled
   Fisher sd per region is 0.77, 0.37, 0.18, 0.087 at nside 16, 8, 4, 2 against a PySM scatter of
   0.06, for white noise, which is optimistic (E031 Table 3, E031:222-229). At nside 2 that is 48
   numbers with an S/N near 1 each: a field amplitude may be detectable in the bright third, a
   correlation length hardly. "Inferred smoothness" would be the hyperprior, the regime of F29 to
   F38. Test, minutes: the exact posterior of (sigma, rho) for a field seen with those errors.
5. *An explicit K against 30 GHz model can be identified.* Unspecified, and degenerate with the
   measurand. K and 30 GHz carry 91 percent of the index information (weights 0.52 and 0.39, check
   1). A relative factor g between them shifts the index by 4.6 ln g, so the 0.37 difference in
   the Fan region ([rep]: Weiland et al. 2022 Table 1; its abstract confirms "inconsistencies
   present in these public data" [abs]) is an 8 percent discrepancy, far above calibration. Only
   Ka (5 percent of the information) tells a K/30 systematic from the index. Test: add a relative
   gain and a few large-scale additive modes in one band to the E031 Fisher script with the real
   covariances. If the region sd grows by more than 3, the model eats the signal.
6. *The analysis is exact and cheap at nside 8 to 16* (Q21). Half supported. One Cholesky of order
   4608 (nside 8, Q and U, three components) took 0.98 s on one thread beside the E032 chain; order
   18432 (nside 16) scales to a minute. A gradient needs the inverse, three to four factorizations.
   At 15 (the K 47 field of E018) to 125 (E032:138) gradients per iteration a 1200-iteration chain
   costs 20 to 167 hours on one thread at nside 8, and SBC with 300 replicates of 300 iterations
   1500 to 3000 core-hours. The feasible rung is nside 8 data with an nside-2 field. Nside 16 needs
   a reduced algebra (marginalize synchrotron only) or a cluster, which is Trygve's call.
7. *QUIJOTE can be added in the north.* Its declination filter acts below l of about 15, most of
   an nside-8 analysis; no pixel covariance is released, 11 and 13 GHz share 1/f noise, and the 5
   percent calibration is 0.068 in beta ([rep], science report lines 58-61, 213). "Dense and
   exact" does not extend to it. Test, seconds: a declination-ring filter on a PySM s1 map.
8. *It is new.* The architecture is WMAP5's. Dunkley et al. 2009 (0811.4280) ran a pixel-based
   Metropolis-within-Gibbs estimate of polarized maps with spatially varying synchrotron indices
   and reported tau = 0.090 +- 0.019 and beta = -3.03 +- 0.04 with an informative prior [abs]; at
   nside 8 with the full N^-1, 48 index pixels and a zero-mean dust Q/U amplitude prior scaled by
   the dust intensity ([rep], science report line 75, novelty report lines 102-111, 179-180). The
   plan does not mention it. Paper B has to be positioned against it.
9. *WMAP plus LFI lack the power, so tau is secondary.* True for settling 0.055 against 0.090. Not
   true for the paired shift of tau between index models on the same data, where the noise
   cancels. Nobody has published that for WMAP plus LFI ([rep], science report Sect. 3.3).

### A4. Real intensity as engineering for upstream

1. *Paper A needs no real intensity.* Contradicted by A1(d), which needs the Q19 rerun.
2. *Offsets and calibration factors are intensity engineering.* The same linear nuisance columns
   with priors (D013, decided 2026-10-27, not implemented; audit Sect. 10.1 item 1) are what a K
   against 30 GHz model is made of. Deprioritizing them delays the polarized track.
3. *Fixing B1 and B14 makes the maps usable.* Untested: D2 (up to 3.3 in width by hand estimate),
   D3, D7 and D10 remain. Test: Q19's own acceptance, by the exact grid.
4. *Upstream will take care of it.* Nobody has told upstream, and both ingest items affect its
   users today. Reporting them costs an hour, needs Trygve's yes, and is the most useful item.

### A5. Process

1. *The regime gate is enough.* It asks for the noise per pixel and the S/N in the leading band
   (D028). E032 passes that and is still not a real regime (A1c item 2). Add the components the
   simulation lacks, the noise correlation, the window, and how the template was built.
2. *Four chains.* Two starts, each shared by a pair (F42). State an ESS floor (400, [mem]) next to
   the R-hat bound and exclude warm-up draws by a recorded count (audit G7).
3. *Five seeds support a coverage statement.* They do not: the binomial error of a 68 percent
   coverage on five independent units is 0.21. Require a replicate-clustered error below 0.03.
4. Missing rules: a time box (two follow-ups on a failed pre-registration, then a design review;
   E023 to E030); literature numbers checked in the PDF; no truth at the centre of a grid (E).

## B. The most likely failure within three months, and what I would do instead

In three months there is a Paper A draft whose (a) and (b) are finished and whose (c) and (d)
still wait on experiments in a regime chosen for convenience, and a polarized prototype that
passes SBC inside its own model and has not met real data, because its first real number depends
on whether K band or 30 GHz is believed. Two half-papers. The cause is the one the audit named:
the validated region and the real-data region do not overlap (audit lines 51-69). The plan freezes
intensity where it would touch data (A4) and starts polarization from new code, new priors and
data that are not on disk (A3), while (c) keeps Paper A open.

What I would do:

1. Split Paper A now. Paper A1 is (a), (b), F21 and F23 as the limits, and a Q/U section made of
   E036 and a dense-N upgrade of E022 with K band: flat against profile against proper priors on
   the four BP XV regions and on nside-2 regions. Everything left is exact grids on one thread,
   nothing in it is graded below B, and its real-data section sits where the prior decides the
   answer. Six weeks.
2. Move (c) to Paper B, where a field is validated in the configuration that uses it. Cut (d).
3. Before any NUTS code for polarization, run the K against 30 GHz null test with dense
   covariances (C2): two weeks, grids only. It decides whether Paper B is about a field or about
   systematics.
4. Report the ingest bugs upstream (Trygve's decision).

## C. Alternatives with more value per unit of effort

- **C1. Paper A1** as in B.
- **C2. K band against 30 GHz as a result of its own.** Region indices from {K, Ka}, {K, 30} and
  {30, 44, Ka} with dense covariances and the profile likelihood with a sandwich error; the
  discrepancy projected on a few named modes (relative gain, WMAP loss imbalance, LFI gain or
  bandpass leakage); WMAP9 with PR3, PR4 and the DR1 maps. The discrepancy is known, unmodelled,
  and moves tau at 70 GHz from 0.053 to 0.094 with the tracer ([rep]: Weiland et al. 2018).
- **C3. tau against the index model as a paired difference.** The dense likelihood with S(tau) on
  a tau grid, with the index constant, in four regions, in nside-2 regions, and under two
  high-latitude priors (-3.0 and -3.3). No NUTS, no field prior. The 2009 version exists (a prior
  of -2.7 +- 0.3 moved tau from 0.090 to 0.092, [rep]); BP XV asserts an effect without a number.
- **C4. Reproduce Dunkley et al. 2009 on WMAP9** as the L5 gate of the polarized code, and add its
  phase-space correction prior p_s / L_0 ([rep]) as an arm of the (b) comparison.
- **C5. The profile likelihood with a sandwich error as the default index estimator.** No prior,
  no reference frequency, no template, valid for dense N, and to first order as informative as a
  correctly specified zero-mean prior (E). It changes the width of a mode that is grade A (F2).
- **C6. A filtered northern analysis.** One declination filter applied to WMAP K and Ka, LFI 30
  and QUIJOTE 11 and 13 GHz; the index of the structure above l = 15 per region, with diagonal
  noise from the weights and half-mission differences. It removes the large-scale K/30 modes by
  construction, has 4.5 times the index information of set B (E031 Table 3, assumed depth), and
  addresses the gap between map-domain values of -3.0 to -3.1 ([rep]) and -3.26 +- 0.05 at |b| >
  20 deg from spectra, with no evidence for spatial variation there ([abs]: 2610.31469).

My order by value per effort: C1, C2, C5, C4, C3, C6.

## D. The strongest abstract from A-grade and B-grade evidence, and what a referee still rejects

> Spectral indices in pixel-based Bayesian component separation share the data with per-pixel
> amplitudes, which are incidental parameters. We analyse a natural joint move, "propose the
> index, re-solve the amplitudes, accept on chi-square": alone it is a pseudo-marginal kernel for
> the profile posterior, and combined with the Gibbs amplitude draw it is a Monte Carlo within
> Metropolis scheme that converges to the profile posterior tempered by sqrt(2 pi k)/4, with
> widths 0.79 k^(1/4) times too large for k re-solved amplitudes (toy chains to 1 percent for k
> above 1000, five simulated configurations to 5 percent); exact collapsed alternatives reproduce
> independently computed targets to 0.1 to 0.2 sd, also on real nside-64 maps. With flat amplitude
> priors the exact index marginal carries a volume term that depends on the reference frequency:
> it biases a full-sky index by 6 to 18 sd, collapses a spatially varying index at every S/N
> tested, and on real WMAP and LFI polarization data shifts the index of the North Polar Spur by
> +0.27 and sends the high-latitude index to the prior boundary. The profile likelihood is
> calibrated above an index-pixel S/N of about 15; below it a full-sky index stays unbiased with
> intervals up to two times too narrow, and an index field is biased by 5 to 13 sd. A hierarchical
> Gaussian amplitude prior with learned scale passes simulation-based calibration on 600
> replicates, but is exact only for its own template-error law: an errors-in-variables template
> error of 30 percent shifts the faint-sky index by 0.11, and a data-built template makes the
> posterior run away at the white-noise depth of WMAP and LFI. On a model-consistent sky at that
> depth a Gaussian Markov random field prior with inferred amplitude and correlation length lowers
> the index-field error by a third against exchangeable shrinkage and by at least a factor of two
> against a single index.

Evidence used: F1 (A), F2 (A), F12, F15, F18 (A), F22 a and b (B), F16 (B), F23 (A, B), F21 (B
plus), F24 (B plus). What a referee would still reject:

1. The move belongs to unreleased code, no published number is shown to be affected, and the one
   real-map width is 22 to 31 percent above the law, unexplained.
2. The law holds for one index step per redraw, small proposals, flat priors, the per-pixel solve
   and k above 1000. The scheme is MCWM, which the text must say and cite.
3. The volume term and its remedies are in print since 2008 to 2017 ([rep]: 0709.1058, 0811.4280,
   1608.00551). New is the coverage, in one configuration, with a threshold where a formula exists.
4. SBC tests the code, not the prior: a correct sampler passes for any prior. The first run failed
   its rule, one index was sampled, and the Commander4 implementation has no SBC.
5. The paper recommends a prior that it shows running away on the first realistic sky. It has no
   positive recommendation for real data.
6. The GMRF result has the truth inside the model class, exact templates, a fixed T, single
   chains, hyperparameter ESS of 3 to 14 and over-coverage, on three components with white noise.
7. The polarization numbers come from a surrogate with diagonal noise, mismatched windows and
   delta bandpasses; its BeyondPlanck comparison ends in "mechanism not identified" (E022:453-454).
8. No real-data fit with an acceptable chi-square, no beams, no CG path, and no second code: every
   reference is the authors' own prototype of the same derivation.

## E. Is a zero-mean hierarchical amplitude prior for Q/U sound when the sky is not drawn from it?

**Setting.** One Stokes element with data d = a f(beta) + n, n ~ N(0, N), f_b = (nu_b /
nu_ref)^beta. Write A = f^T N^-1 f, a_hat the GLS amplitude, u = a^2 A (true amplitude S/N
squared), x = sigma_p^2 A (prior variance over noise variance), and L and V the mean and variance
of ln(nu_b / nu_ref) under weights f_b^2 / N_b. L is the log distance from the reference to the
information pivot. For band set B at beta = -3.1 with the E031 white noise the pivot is 25.9 GHz,
L = -0.146 at a 30 GHz reference and +0.128 at 22.8 GHz, and V = 0.025.

**Three identities** (my algebra; check 1 confirms the second numerically, check 4 the third):

1. The zero-mean marginal is the flat-prior marginal times N(a_hat(beta); 0, sigma_p^2 + 1/A(beta)).
2. For a fixed true amplitude, not a draw from the prior, the expected score at the true index is
   L (u - x) x / (1 + x)^2 per element, and that of a common log scale is (1/2) (u - x) x /
   (1 + x)^2. The model's information on beta is x^2 [2 L^2 / (1 + x)^2 + V / (1 + x)].
3. The profile score is unbiased for any amplitudes, with curvature u V and variance (u + 1) V.

**What follows.**

- The zero-mean prior estimates the index from second moments: the data covariance along f must
  equal sigma_p^2. It is the errors-in-variables problem of F23 moved up one moment. What must be
  right is E[a_p^2 given the scale map]. Heavy tails and spatial correlation of the amplitudes do
  not bias the index ([mem]: robustness of normal-theory factor loadings, as E036 P1 cites).
- **If the index and the scale share a region, it is sound.** The two scores are proportional
  (factor 2 L), so the scale absorbs any misfit of the level and the index is consistent for any
  sky. With a flat prior on the log scale the target does not depend on the reference at all: a
  change of reference multiplies every amplitude of the region by one factor, which the scale
  takes (check 1: shift 0.00000 at both references). The information is then x^2 V / (1 + x) per
  element, which is the sandwich-corrected profile, (u V)^2 / ((u + 1) V), at u = x. They are the
  Bayesian and the frequentist form of one estimator, and E036's first output shows them agreeing.
- **If the scale is shared across index regions, it is not.** That is the plan as worded: one
  global scale times a template. Region r keeps the expected score L sum_p (u_p - x_p) x_p /
  (1 + x_p)^2. A region whose mean-square amplitude differs from the scale map's by a factor rho
  moves its index by about L (rho - 1) / (x V) at high x and L (rho - 1) / (2 L^2 + V) at low x.
  Check 1, two zones of 24 nside-2 regions, the bright zone at x = 300, the faint zone at x = 10
  with a scale map 30 percent high, truth -3.1 everywhere: the faint-zone index moves +0.125 at a
  30 GHz reference, -0.120 at 22.8 GHz and 0.000 at 25.9 GHz, against a pooled model sd of 0.035.
  At x = 3 it is +0.31 and -0.34; at x = 30, +0.046 and -0.042 (2.2 sd pooled). A map 30 percent
  low gives -0.22 and +0.17 at x = 10. The bright zone moves by under 0.007.
- This is F38 and E030 chain 3 again (a displaced faint sky, a sign that flips with the
  reference), in polarization and of the size of the signal: the published high-latitude
  steepening is 0.25 to 0.4 ([abs]: Weiland et al. 2022, 2203.11445).
- **A data-built scale map is wrong in exactly this way.** sqrt(Q^2 + U^2) of a noisy map is biased
  high where the S/N is low ([mem]: the Rician floor), so rho is below 1 in the faint sky: a
  spurious flattening at a 30 GHz reference, a spurious steepening at a K-band reference. A map
  built from K band is also the data used twice (F21).
- **K band creates the lever.** For the BP XV band set the pivot is 29.8 GHz and L at a 30 GHz
  reference is -0.006 (check 1): a prior stated at 30 GHz was nearly neutral there, and is not now.
- **What remains exposed in the per-region form.** (i) A GMRF field finer than the scale partition
  restores the lever inside each scale unit. (ii) The pivot varies inside a region with the WMAP
  and LFI depth patterns, and by mode once the covariance is dense, so the two scores are
  proportional only on average. (iii) Synchrotron and dust are correlated at 0.1 to 0.25 ([abs]:
  2610.31469) and independent priors ignore it: a second-moment misfit of the same kind, size not
  computed. (iv) At low region S/N the scale hits its lower bound and the index reverts to its
  prior.

**Against the alternatives.**

- *Flat priors with conditional Gibbs.* The target is the flat-prior marginal (F15), with a shift
  of -L / (u V) per region whatever the pixel count: -0.06 and +0.06 at 23 and 30 GHz in E036's
  middle class, +0.27 for the Spur in E022. Unsound here. Keep it as the negative control.
- *A C(l)-type Gaussian field prior.* Right for the CMB (A3 item 3). For synchrotron it is the
  zero-mean prior with a global scale and no scale map, on a sky whose amplitude varies by an
  order of magnitude from plane to poles: the lever at full strength unless the amplitude is
  stated at the pivot. A modulated version needs the same per-region scale freedom.
- *Profile with a sandwich error.* Unbiased score, no reference, no template, valid for dense N,
  widths corrected by 1 + n / sum u. It loses to a correct hierarchical model only through the
  weights when the S/N varies inside a region, and it fails where the index posterior is not
  locally Gaussian (sd of order 0.5: sd z 2.4 in check 4).

My order: the profile with a sandwich error or the per-region zero-mean prior for synchrotron
(tied to a template only softly, with a log-normal width no tighter than the template error), a
C(l) prior for the CMB, and never a global scale on a template.

**The minimal simulation that would show the problem.** E036 with four changes, all grids on one
thread. (1) Add the estimator "one global scale times a scale map" at references of 23, 26 and 30
GHz. (2) Scale maps: the oracle; the oracle 30 percent high above |b| = 30 deg; sqrt(Q^2 + U^2) of
the simulated K map itself smoothed to 10 deg; the same from an independent noise draw. (3) Truths
of -2.8 and -3.4 and a latitude step of 0.3. E036's truth of -3.1 is the centre of its support
[-4.6, -1.6] (E036:29, 47-49), so reversion to the prior reads as zero bias; that matters in its
class above 0.25, where 14 to 17 percent of the posteriors sit at the edge. (4) A scale partition
coarser than the index partition (scales at nside 1, indices at nside 2 and 4). Identity 2
predicts a sign flip between 23 and 30 GHz, zero at 26 GHz, and the size L (rho - 1) / (x V).

## F. Go or no-go for the first polarized milestone

The milestone: a simulation at real depth with the real WMAP K, Ka, Q, V covariances and LFI,
nside 8 data, an nside-2 index field, calibrated at the configuration reported. Thresholds are
mine. The Fisher expectation behind Gate 3 is a latitude-step sd of about 0.03 (E031 Table 3: 24
regions of sd 0.087 per zone, LFI variance times 2.65).

- **Gate 0, inputs (one week).** (a) The K-band res-4 map and N^-1 are fetched; N^-1 is symmetric
  to 1e-10 and positive semi-definite, with the processing mask as its null space; E022's
  conditional with the dense Ka, Q, V matrices is reproduced within 5 percent. (b) An LFI 30 GHz
  dense covariance passes a half-split null test with a PTE between 0.01 and 0.99, and its
  diagonal is within 20 percent of BP's simulated rms (1.65 uK, E022:129). If not, LFI enters with
  a diagonal inflated by 2.65 in variance and "exact" is dropped for LFI. (c) A beta = -3.1 sky
  pushed through the WMAP res-4 path and the LFI low-resolution path returns region indices
  within 0.01 with the noise off. If not, one window for both sets comes first.
- **Gate 1, exactness.** Dense collapsed log likelihood against brute force at nside 4 to 1e-8;
  gradient against central differences to 1e-5; a scalar chain against its grid within 4 Monte
  Carlo standard errors and 15 percent in width; the target at references of 23 and 30 GHz equal
  to 1e-6 nats up to a constant.
- **Gate 2, SBC at the reported configuration.** At least 300 replicates, ESS above 200 per tested
  quantity, ranks uniform by a band fixed in advance for the field mean, log sigma, the range, ten
  pre-chosen regions and the amplitude hyperparameters. After a failure the cause is found first.
- **Gate 3, calibration on a sky that is not drawn from the model** (this gate decides). Truth: a
  data-derived Q/U amplitude field as in E036, plus the CMB and dust with a 0.2 synchrotron-dust
  correlation; indices constant at -2.9, constant at -3.3, and a step of 0.3 at |b| = 30 deg; 100
  noise draws from the real N. Required: bias of the |b| > 30 deg mean index below 0.03 and below
  0.5 of its posterior sd; the step recovered with bias below 0.05 and posterior sd below 0.08,
  and the null truths giving a step inside 2 sd in at least 95 of 100 draws; coverage of region
  indices of 0.60 to 0.76 at 68 percent and 0.91 to 0.98 at 95 percent, clustered by draw; the
  |b| > 30 deg mean on the same data at references of 23 and 30 GHz within 0.02, and within 0.03
  across scale maps (none, oracle, data-built, independent 30 percent error); region-map RMSE
  below the constant at the truth mean and below independent nside-2 regions with an N(-3.1,
  0.3^2) prior. If the field does not beat the regions, it is dropped for regions.
- **Gate 4, chains.** Split R-hat below 1.01 and bulk ESS above 400 for every hyperparameter and
  for the |b| > 30 deg mean, over four chains from two dispersed starts; no divergences after
  warm-up; E-BFMI above 0.3; warm-up excluded by a recorded count.
- **Gate 5, cost.** One gradient in at most 5 s on one thread at the reported configuration.
  Otherwise the configuration is reduced before the SBC, not after.
- **Gate 6, before any real index is read.** The |b| > 30 deg mean index from {K, Ka}, {K, 30} and
  {30, 44, Ka} agrees within 2 combined sd, or a model with at most five named nuisance modes
  brings the K against 30 GHz residual chi-square to a PTE above 0.01. If neither holds, the
  result is reported as limited by systematics, with the three values: a different paper.
- **No-go.** A Gate 3 bias above 0.05 that no parameterization removes; a posterior sd of the
  latitude step above 0.10 in simulation, since the data then cannot answer the question; a Gate 6
  disagreement above 3 sd with no identifiable mode set.

## G. What the plan ignores and should not

1. **Dunkley et al. 2009** did the proposed architecture on WMAP5 (A3 item 8). It is in the
   science report's table and absent from the plan.
2. **F8 and Q7.** On the validation simulation the exact ridge target sits 1.9 and 3.7 sd from
   the input (`FINDINGS.md:178-181`). By principle 8 that is a bug until shown otherwise, and it
   sits under the harness that carries F2. The referee's test (20 seeds, grids only) is unplanned.
3. **E018 and F7.** The L1 field test disagrees with its reference at template brightness and was
   moved to a brighter sky (referee Sect. 2.4); (c) leans on that implementation. Per-pixel
   fluctuation draws stored as alms lose variance above lmax (F7): a candidate for Q23, and every
   delivered amplitude sample is under-dispersed.
4. **Q17.** E022 could not say why the flat-prior pull is absent from the BeyondPlanck chains, and
   Paper B will be compared with BP XV. Reading the Commander3 source for `COMP_CL_TYPE = none`
   is cheap and has not been done.
5. **The polarized track runs outside Commander4.** Commander4 has no Q/U covariance (audit D1),
   one index for I and Q/U, a hierarchical prior for intensity only (audit S2;
   `compsep/processing.py:373-376`) and scalar smoothing of Q and U (B2, unchecked). Paper B as
   planned exercises none of the fork. That may be right, and it should be a recorded decision.
6. **Unverified inputs.** The science report's own risk list for WMAP plus LFI (NPIPE transfer
   function at 44 and 70 GHz, an LFI covariance that misses part of the error budget, simulations
   behind a NERSC account; its lines 44-47). The 11 GHz depth and the 143 and 217 GHz
   sensitivities, from memory (E031:46-49, 54-58). The 0.37 in the Fan region and tau by tracer,
   summarizer extracts. PySM's index scatter of 0.06, which sets every "prior-dominated" statement
   of E031 while published map dispersions are 0.09 to 0.21 ([rep]; E031:257-258). Check the five
   or six load-bearing numbers in the PDFs before they gate anything.
7. **C-BASS.** A public 5 GHz release would multiply the lever arm by about six and move the
   pivot. I found none today ([search]). Paper B's design should be able to take it.
8. **Records behind the review.** `STATUS.md` (last entry 12:22) and `ROADMAP.md` (2 November)
   predate the review. At 23:32 the program review that Q21 cites for go or no-go numbers had
   placeholders and E036 had no registry row. The E034 and E035 records end in duplicated
   "Conclusion: Pending" stubs (E034:191-200, E035:176-185).
9. **E032's exact-template faint dust quartile** already sits beyond its scan edge (+0.039, 8
   nats; E032:93-97). The pre-registration reads the cs cells only, and the exact-template chain
   listed as a follow-up would be read against a scan that has already missed the 0.03 tolerance.

## Appendix: the four checks

One thread, numpy and scipy, script on stdin, nothing saved. Band set B: K, LFI 30, Ka, Q, LFI 44,
V, LFI 70 with the Q/U white noise of E031 Table 1 divided by 4 for nside 16, beta = -3.1.

1. **Zero-mean prior, wrong scale map.** Expected -2 log marginal per element, tr(C^-1 Sigma_true)
   + ln det C with C = N + sigma^2 f f^T, minimized over (beta_bright, beta_faint, log scale); two
   zones of 24 regions of 128 elements; sds are model curvatures for independent elements. Pivot
   25.92 GHz; L = -0.1462 and +0.1282 at references of 30 and 22.8 GHz; V = 0.0251; band weights
   0.522, 0.391, 0.052, then below 0.02. Without K band: pivot 29.81 GHz, L = -0.0064. With one
   scale per zone the shift is 0.00000 at both references.
2. **First correction to the E035 law.** T = 1/(4G) + 1/2, G = E[exp(-Z/2); Z > 0], Z the
   difference of two chi-square_k variates. Widths at k = 10, 100, 1000, 1e4: closed form
   sqrt(sqrt(2 pi k)/4 + 1/2) 1.575, 2.601, 4.507, 7.948; Monte Carlo G 1.591, 2.606, 4.505,
   7.960; the E035 toy 1.577, 2.584, 4.463, 7.870; the leading law 1.408, 2.503, 4.451, 7.916.
3. **Dense factorization.** Cholesky of order 4608: 0.98 s beside the E032 chain; 18432: 63 s.
4. **Profile sandwich.** One index, 256 amplitudes of equal S/N, 400 noise draws: mean z -0.09,
   -0.06, -0.07 and sd z 2.44, 1.09, 1.04 at u = 0.38, 3.35, 30; sqrt(1 + 1/u) = 1.91, 1.14, 1.02.

Limits. The checks are toys: white diagonal noise, one component, independent elements. With CMB
and dust columns L and V become Schur-complement quantities that I have not rederived. I did not
verify the E031 Fisher numbers or the E036 script. Statements marked [rep] or [mem] are unverified.
