# Paper A (methods): exact spectral-parameter sampling in Bayesian component separation

Working outline, 2026-10-27. Every claim below cites an experiment ID (`experiments/registry.md`)
or a finding (`FINDINGS.md`); items marked (pending) have no result yet and name the experiment
that will supply it. Target journal A&A (methods section of the Cosmoglobe series) or JCAP.
Writing preference: no em-dashes.

## Working title
Exact sampling of spectral parameters in pixel-based Bayesian component separation: the
prior-volume problem, collapsed likelihoods, and hierarchical index fields with inferred smoothness.

## Thesis (one paragraph)
In Gibbs-sampled component separation the spectral parameters (synchrotron and dust indices) are
updated conditionally on, or jointly with, foreground amplitudes that carry flat priors defined at
an SED reference frequency. We show that (1) the acceptance step used in production code does not
leave the joint posterior invariant when the amplitude redraw is a constrained realization (F1,
E009/E010 for Commander3, E014 for Commander4), (2) the flat-prior joint posterior itself has a
volume term that moves a full-sky index by tens of posterior standard deviations and collapses a
spatially varying index at low signal to noise (F3, F12, F15), (3) a collapsed (amplitude
marginalized) likelihood with a hierarchical amplitude prior whose scale and scatter are learned
from the data gives calibrated index posteriors down to index-pixel S/N below 1 (F16, E016; E017
for the full-sky case, pending), and (4) this target, with a Gaussian Markov random field index
field whose amplitude and correlation length are inferred, is implemented in Commander4's per-pixel
path and reproduces exact targets end to end (E018, pending L3). As a by-product we show that on
real low-resolution maps the full-sky index is set by the zero-level convention of the inputs at
the 0.03 to 0.4 level, hundreds of posterior sd (F17), which turns the monopole treatment into a
model component rather than a preprocessing choice.

## 1. Introduction
- Pixel-based Bayesian component separation (Commander, BeyondPlanck, Cosmoglobe DR1/DR2): Gibbs
  over amplitudes and spectral parameters; the role of spectral indices for B-mode foreground
  cleaning. Cite the standard references (LITERATURE.md).
- The two questions of this paper: is the index step exact, and under which prior is the answer
  meaningful? Statement of contributions (thesis items 1 to 4).

## 2. The per-pixel model and its collapsed likelihoods
- Data model d_p = M_p(beta) a_p + n_p at a common resolution, diagonal noise; the mixing matrix,
  reference frequencies, units (uK_RJ).
- Collapsed likelihoods: flat-prior marginal (w = 1, includes -1/2 log det A), profile or ridge
  (w = 0), fixed Gaussian prior (A' = A + T^-1, b' = b + T^-1 m) and the hierarchical prior
  a_p,c ~ N(mu_c g_p,c, (s_c g_p,c)^2) with the trailing terms -1/2 m^T T^-1 m - 1/2 log det T
  that carry the dependence on (mu_c, s_c). Exact expressions and the constants dropped
  (plans/M2.4 section 1; stats-reviewed twice).
- The amplitude Gibbs draw under each prior: N(A'^-1 b', A'^-1); the pseudo-band implementation
  (one kernel for likelihood, gradient and draw; E018 L0 tests).
- Analytic gradients with respect to the indices (per region) and the hyperparameters.

## 3. The acceptance step in production samplers is not exact
- Commander4 legacy `resolve_chisq`: accept on chi-square at redrawn amplitudes; posterior 16 to 17
  times too wide and beta_d shifted 3 sd on real maps (F1, E009, E014 pass 2). Commander3: the
  BeyondPlanck-era samplers do not have the bug, but `marginal` has the volume problem (F10, E010).
- The exact alternatives: ridge and marginal chains reproduce independent grid targets to 0.1 sd
  on 17 and 20 bands (E014). The Gibbs conditional step cannot mix for a full-sky index (F14).
- Figure 1: chain against exact grid for the four modes (E014 pass 2 and 3 tables).

## 4. The prior-volume problem
- Full-sky index: marginal minus ridge +56 sd (17 bands), +16 sd (20 bands), +23 sd (pass 4):
  the shift depends on the data set through the noise levels and the mixing matrix (E014).
- Spatially varying index: the flat-prior marginal collapses below index-pixel S/N 15; the
  profile is calibrated above S/N 15 and biased by 5 to 13 sd below (F12, E012); the joint
  flat-prior posterior has exactly this marginal, so any Gibbs chain with that target either shows
  the collapse or is not mixing (F15).
- Full-sky coverage table for w0, w1, gauss and hier at S/N 680 to 23 (E017, pending): profile
  calibrated at high S/N, 2 times too narrow at S/N 23; w1 unusable at every S/N (6 to 20 sd,
  edge pile-up below S/N 70); hier calibrated (provisional numbers in E017).
- Figure 2: z and coverage versus S/N, spatially varying (E012/E016) and full sky (E017).
- Discussion: Jeffreys-type factor of the ridge; why the fixed Gaussian prior is not a fix (4 to 10
  sd bias where the index is well measured; E012, E017).

## 5. Hierarchical amplitude prior and the GMRF index field
- Model: index field beta = mu 1 + sigma S(rho) z on the graph-Laplacian basis of the HEALPix
  index grid, Matern alpha 2, constant mode dropped, hyperpriors on (mu, sigma, rho) (D009);
  amplitude hyperparameters (mu_c, log s_c) with their hyperpriors; the collapsed joint target;
  non-centred and centred parameterizations and their Jacobian.
- Validation ladder: exact-target checks (E012 L1), SBC at K 192 (E012; hier SBC E016 L2,
  pending), coverage with a misspecified template (E016: a 30 percent large-scale distortion or 23
  percent scatter of the template leaves z within 0.3 sd and coverage 0.64 to 0.75 / 0.92 to 0.98).
- Sampler: NUTS with dual averaging and windowed diagonal mass adaptation; efficiency per gradient;
  the funnel at low S/N (E018 L1 diagnostic: the centred reference chain does not reach the
  small-sigma tail where NUTS in non-centred coordinates does). MAMS comparison (M2.3, pending).
- Figure 3: recovered index field and hyperparameters against truth on template skies (E016).

## 6. Implementation in Commander4 and end-to-end validation
- Per-region mixing in the per-pixel solver, priors as pseudo-bands, `RegionIndexLikelihood`
  over MPI, `likelihood: hier` groups; parity statements (bitwise with the flat path when no
  field and no prior; E018 L0).
- End-to-end: 12-band simulation at nside 16, field at nside 2, 3000 Gibbs iterations against the
  prototype chain and per-region references (E018 L3, pending); cost per iteration.
- Table: E018 comparison (means within pooled SE, sd ratios, ESS per gradient, divergences).

## 7. Real data: the zero-level convention dominates the full-sky index
- E014 passes 2 to 4 on nside-64 DR1 + PR4 + Haslam intensity maps: beta_s from -3.34 to -2.86
  depending only on the input zero levels (F17); a consistent convention from the DR1 sky model
  moves beta_s by 260 sd and beta_d by 460 sd; the Haslam dipole alone 80 sd.
- Consequence: per-band monopole (and dipole) parameters with priors belong in the model (Q15; the
  offset is linear and collapses like an amplitude). Implementation and its effect: Phase 3
  (M3.3), out of scope here beyond the measurement.
- Figure 4: the zero-level ladder (table) and the pass-2 to pass-4 posteriors.

## 8. Conclusions
- Exact index sampling is available at the cost of one collapsed-likelihood evaluation per
  proposal; the ridge is the pragmatic default for full-sky indices (D003), the hierarchical prior
  the calibrated choice for spatially varying indices; the monopole is a model component.
- What is not shown: polarization (QU streams with the hier prior are excluded in stage 1),
  CG-solved amplitudes (Q5, M3.1), beams beyond a common resolution.

## Status of the evidence (2026-10-27)
| Section | Evidence in hand | Pending |
|---|---|---|
| 3 | E009, E010, E014 (passes 2, 3) | none |
| 4 | E012, E014, E015, E016 | E017 Conclusion (hier modes finishing) |
| 5 | E012 (L1, SBC), E016 (coverage) | E016 hier SBC (86 of 300), M2.3 MAMS table |
| 6 | E018 L0, L1 | E018 L3 (running), stage-2 efficiency numbers |
| 7 | E014 pass 4, F17 | none for the measurement; Q15 decision for the fix |

## Open decisions that change the paper
- Q15: with a monopole component in the model, section 7 becomes a result (indices with the
  zero-level uncertainty propagated) instead of a caveat.
- Whether to include the Commander3 code-path audit (E009/E010) as an appendix or a separate note
  to the Cosmoglobe team (Trygve's call; contacting the team needs his go-ahead).
