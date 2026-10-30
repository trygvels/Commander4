# Open questions

Ranked by how much the answer changes what we build next. Each entry: why it matters, how to answer
it, and what the answer decides. Mark an entry RESOLVED with the experiment ID when answered.

**Ranking (2026-10-29):** Q16 first (robustness of the hierarchical index posterior to template
law, template scatter and the dust SED on the PySM skies at realistic noise: it decides whether
Paper A's method runs on real intensity maps), then Q15 (zero levels; blocks any real-data index
result and needs Trygve's model decision), then the Q4 efficiency items (slow mixing of the
field mean in the centred parameterization; MAMS and preconditioners), then Q8/Q7 (M2.6), then
the second half of Q14 and Q5 (CG amplitudes). Answered: Q1 (E014), Q2 (E012, E016, E017), the
core of Q3 (E021, F24: GMRF prior of record; the resolution rule remains open, D019). The entries below keep their
original numbering.

## Q1. Does the legacy acceptance matter on real-data configurations? (ANSWERED at nside 64, E014)
- **Answer (2026-10-26):** yes. On DR1 WMAP + LFI, NPIPE 353 / 857, PR3 100 and Haslam at nside 64
  with the dx12 mask, the legacy `resolve_chisq` index posterior is 16 to 17x too wide and the dust
  index is shifted by 3 ridge sd (E014, exact grid reference). The scaling of F1 predicts factors
  above 50 at nside 512 to 1024. The `marginal` chain lands on its own exact target within 0.05
  sd and sits 54 ridge sd from the ridge target in beta_s (prior effect, see Q2); the `conditional`
  chain cannot mix on real data (F14). Remaining: the tuned round (running), IQU, and a direct run
  at nside 128 to confirm the scaling (open, lower priority). Ranking after this answer: Q2 first,
  then Q3/Q4 (M2.x), Q14 second half; D010.
- **Why:** F1 shows the original move inflates and shifts beta posteriors in simulations. We need
  to know how much it changed real Commander4 results before raising it upstream.
- **How:** Run `params/Planck_mapsonly/param_LFI+HFI+WMAP_perpix_betasamp.yml` (and the CG variant)
  on public Cosmoglobe DR1 / Planck PR4 / WMAP maps, first downgraded to nside 64 to 128, with
  `likelihood: resolve_chisq`, `ridge` and `conditional`. Compare beta posteriors, CMB residual
  maps, chi^2 maps, and acceptance rates.
- **Decides:** urgency and framing of the upstream pull request; whether published Commander4
  preliminary results need rechecking.

## Q2. Which prior gives calibrated posteriors for spectral parameters?
- **Why:** F3/F4: flat amplitude priors can make the index posterior useless; ridge is equivalent
  to a joint Jeffreys prior. Coverage has not been tested for either, and the answer differs between
  full-sky and spatially varying indices.
- **How:** Simulation-based calibration at low resolution (nside 8 to 32): draw beta (and amplitudes)
  from the prior, simulate data, sample with marginal, ridge, conditional and an explicit Jeffreys
  prior; check rank uniformity (Talts et al. 2018) with the joint log likelihood as an extra test
  quantity (Modrak et al. 2022). Then misspecified tests on PySM skies.
- **Decides:** the default prior for every nonlinear parameter in the new sampler; a central result
  for the methods paper.
- **Design note (E012 v2, D008):** SBC can only test the sampler when the generative amplitude
  prior is the one in the sampled target. Sampler correctness is therefore tested with a proper
  Gaussian amplitude prior and its exact marginal; ridge, flat-prior marginal and Jeffreys are
  compared as a frequentist coverage study (bias and coverage versus S/N), as principle 4 requires.
- **Real-data size of the effect (E014, 2026-10-26):** at nside 64 with 17 bands the flat-prior
  marginal and ridge targets for beta_s differ by 0.019, which is 54 posterior sd, with identical
  widths and chi-square. The prior is the dominant systematic for full-sky indices on real data;
  Q2 ranks above sampler efficiency until the coverage study settles it.
- **Coverage study result (E012 study_v1 and diagnostics, 2026-10-26, F12):** for a spatially
  varying index the flat-prior marginal is unusable at every S/N, the profile likelihood is
  calibrated only above index-pixel S/N of about 15, and a Gaussian amplitude prior is calibrated
  only when it is the generative one; widening it moves the answer toward the flat-prior collapse.
- **Hierarchical amplitude prior (E016, 2026-10-26, F16 provisional):** learning the per-component
  amplitude scale and relative scatter from the data gives calibrated index posteriors at every S/N
  tested (z within 0.35 sd, coverage 0.64 to 0.73 / 0.92 to 0.97 over all pixels), also with a shape
  template wrong by 30 percent on large scales or 23 percent per pixel. Pending: the SBC of the
  hier target (sampler correctness; running), the scalar-index version (M2.2, E017 running), and
  the real-data step (which template, and the Commander4 implementation).
- **Scalar index (E017, 2026-10-27, F18):** flat marginal biased 6 to 18 sd and collapsed to the
  support edge below full-sky S/N 70; fixed Gaussian 4 to 10 sd off; profile unbiased but 15
  percent too narrow at S/N 68 and two times at S/N 23; hier unbiased and calibrated at all four
  S/N, also with a distorted template. Answer for Q2: learn the amplitude prior (hier); ridge is
  acceptable only where the index S/N is high (full-sky S/N above a few hundred, index-pixel S/N
  above 15). In Commander4 since E018 (F19). Remaining for Q2: the hier SBC (L2 gate, running), the
  template choice on real data, and the index-resolution rule from the S/N map (Q3).

## Q3. What is the optimal spatially varying index model? (CORE ANSWERED, E021)
- **Answer (2026-10-29, E021, F24, stats-reviewed):** on a model-consistent, spatially smooth
  sky the hierarchical GMRF prior beats exchangeable (iid) shrinkage and independent regions on
  both indices and on the CMB at nside 16 and 8 (nside 16: 0.033 / 0.018 against 0.051 / 0.028
  and 0.142 / 0.101; suggestive at 4); iid is biased where the sky differs from its
  plane-weighted global mean and is dropped as a candidate (D018). Against the S/N-weighted
  region reference the nside-16 GMRF is within 13 / 30 percent of the per-region error of the
  nside-4 model with 16 times the data, while independent regions lose factors 3.5 / 4.9, and the
  finer field keeps the structure. The resolution rule (as fine as the budget allows under a
  spatial prior) is untested in the direction that matters, a field finer than the sky's
  structure (D019; E023 and Q16). Still open inside Q3: the alm GRF and heavy-tailed
  (edge-preserving) variants and the clustering baselines were not run; joint fields for beta_d
  and T_d.
- **Status 2026-10-28:** E020 (PySM sky at 0.03 x E012 noise) closed as a design flaw (F21: the
  misspecified sky dominates); the ten-option comparison runs as E021 on a model-consistent sky.
- **Why:** Trygve's direction: indices as fine-grained as S/N allows, statistically optimal rather
  than hand-drawn regions; hard region boundaries are unphysical and alms ring.
- **Candidate (leading):** hierarchical Gaussian Markov random field on HEALPix pixels at nside_beta
  (SPDE/Matern form, sparse precision on pixel neighbours), with mean, amplitude and correlation
  length sampled as hyperparameters in non-centred coordinates. The posterior mean is then a
  Bayesian Wiener filter of the index field: data-driven where S/N is high, shrunk where it is low.
- **Variants to compare:** alm Gaussian random field with sampled C_l^beta (Commander3 style, but
  with hyperparameters inferred); heavy-tailed increments (Student-t / horseshoe on neighbour
  differences) to allow edges; non-stationary amplitude; K-means clustering (FURAX) and the fixed
  nside [64, 8, 4] ladder (fgbuster) as baselines.
- **How:** validation ladder levels 2 to 4: exact-target checks, then sims drawn from each prior,
  then PySM d1s1, d10s5, d12s7 and the PEGS skies. Metrics: index-map error and coverage, CMB
  residual power at l < 100 (and r bias where relevant), chi^2 maps, ESS per CPU-hour.
- **Sub-questions:** nside_beta per parameter; whether the funnel (Grumitt et al. 2020) appears and
  whether non-centring plus MAMS fixes it; joint vs separate fields for beta_s, beta_d, T_d.
- **Decides:** the core of the methods paper and of the science application.

## Q4. Which gradient sampler for the index field?
- **Why:** Literature reports MCLMC/MAMS 3 to 100x more ESS per gradient than NUTS at field level.
- **How:** Same targets as Q3 at nside_beta 16 to 128; compare NUTS (ours), MAMS (BlackJAX or a
  NumPy port), with diagonal, Fisher and low-rank-plus-diagonal preconditioning; carry adaptation
  across Gibbs iterations.
- **Decides:** the production sampler.
- **Answer for Phase 2 (E019, 2026-10-27, F20):** on the non-centred hier target at K 11 to 191
  MAMS equals NUTS in ESS per gradient (0.9 to 1.2 at the median pixel) with a 2 to 20 times
  dearer warm-up; NUTS stays the production sampler. Reopen only for d above a few thousand (the
  joint target of Q5, or nside_beta 16 and above), with the reference MAMS implementation and a
  non-diagonal preconditioner.

## Q5. How to keep the collapse benefit when beams couple pixels?
- **Why:** F5's 9 to 22x gain is for independent pixels. With CG amplitudes (beams, C_l prior), the
  exact marginal needs ln det of the CG operator.
- **Options:** (a) low-resolution collapsed likelihood for the index field plus full-resolution CG
  amplitudes, quantifying the approximation; (b) stochastic log-determinants (Hutch++, XTrace, SLQ)
  using the existing CG preconditioner, with noise-aware MH; (c) joint HMC/MAMS over amplitudes and
  indices with hand-written SHT adjoints (as Flinch does).
- **Decides:** whether the method extends to full-resolution polarization analyses.

## Q6. Should per-pixel amplitudes be stored as maps instead of alms?
- **Why:** F7. The alm round trip at lmax = 3 nside - 1 loses pixel-scale content and inflates chi^2.
- **How:** measure variance loss and chi^2 bias at realistic settings (nside 512, 60 arcmin); prototype
  map-space storage for per-pixel groups.
- **Decides:** a small upstream fix, and whether per-pixel chains are trustworthy at small scales.

## Q7. How wrong is the white-noise assumption on smoothed maps?
- **Why:** F8. Smoothing correlates noise; the per-pixel likelihood ignores it.
- **How:** compare with a low-resolution likelihood using dense or banded noise covariance; quantify
  bias and coverage versus smoothing scale.
- **Decides:** whether the per-pixel path is acceptable for science, or needs a noise-covariance
  model.
- **Status (2026-10-28, E022 side measurement):** on BeyondPlanck's LFI maps at nside 16 and
  60 arcmin, the noise rms from BP's noise simulations is 1.6 times the white-noise-propagated
  rms in all three bands (ratio 1.4 to 1.9 across the sky): correlated noise raises the
  large-scale noise variance by 2.5 to 2.8. A white-noise per-pixel likelihood at that
  resolution overstates the index precision by 1.6 unless it uses simulated or measured
  low-resolution rms maps (as BP did). The bias as a function of smoothing scale and the effect
  on coverage are still to be measured.

## Q8. Do moment components absorb sub-resolution variation?
- **Why:** Line-of-sight and in-beam averaging distort SEDs (Chluba et al. 2017); moment templates
  are linear and fit naturally in the CG amplitude step.
- **How:** add first-order moment components (d f / d beta times an amplitude map) with Gaussian
  priors; test on d10s5 and d12s7.

## Q9. Which external data are actually available, and in what form?
- **Why:** low-frequency data (C-BASS 5 GHz, QUIJOTE MFI 11 to 19 GHz, S-PASS 2.3 GHz, CLASS 40 GHz)
  are what breaks synchrotron degeneracies; high-frequency (HFI 545/857, DIRBE) constrains T_d.
- **How:** verify download locations, formats, beams, calibration uncertainties, masks and transfer
  functions; record in DATA.md.
- **Status (2026-10-25, lit-scout with fetched URLs, DATA.md "verified" section):** public: QUIJOTE
  MFI DR1 (LAMBDA, nside 512 IQU with covariance; FDEC declination filtering must be applied to
  external data), S-PASS DR1 (LAMBDA/INAF, IAU convention, dec < -1 deg), CLASS 40 GHz DR1 (LAMBDA,
  nside 128 with a C_l-domain transfer matrix only), Haslam 2014 (LAMBDA). NOT public: C-BASS North
  maps (the research/02 claim was wrong). Open: how to implement the QUIJOTE FDEC and CLASS filters
  in a joint likelihood (transfer-function problem); PR4 mirrors (NERSC slow, CMB-S4 Globus links
  dead on 2026-10-25, Hugging Face Parquet copy exists).

## Q10. How long are Commander4's Gibbs correlation lengths for gain against sky?
- **Why:** gain/sky and noise/sky degeneracies drive BeyondPlanck's long chains.
- **How:** run Commander4's simulation configs with TOD steps; measure autocorrelation of gains,
  sky modes, indices.

## Q11. C(l) sampling (missing in Commander4)
- **How:** implement with Racine et al. 2016 rescaling plus Ducrocq et al. 2021 auxiliary-variable /
  overrelaxed moves; validate against exact low-l posteriors.

## Q12. Can CG solves be made ~10x cheaper?
- **How:** warm starts and Krylov recycling across Gibbs iterations, multilevel or pseudo-inverse
  preconditioners, learned preconditioners inside CG; measure iterations and wall time on
  `simparam_compsep_only.yml`.

## Q13. Strategic coordination
- Cosmoglobe is building OpenHFI (end-to-end HFI) and plans AKARI, FIRAS, SO, LiteBIRD work. Which
  of our results should be upstreamed, which published independently, and when to talk to Jonas
  Lunde / Hans Kristian Eriksen? This is Trygve's call.

## Q14. How much did the Commander3 `marginal` sign error and 2024 MH module change published results? (BeyondPlanck part RESOLVED, E013)
- **Why:** E009/F10. BeyondPlanck final sampled AME nu_p with the sign-flipped `marginal` lnL
  (DR1 did not: it used the exponential AME model with `prior` index steps, see E013); DR2-era runs may have used `sample_specind_mh` with CG updates (F1 structure).
- **How:** rerun the AME nu_p step on BeyondPlanck / DR1 inputs (or on DR1 chain samples as fixed
  inputs) with the sign corrected and with `ridge`; compare nu_p and the AME amplitude map. For DR2,
  find the production parameter files and check `MCMC_SAMPLING_GROUP_UPDATE_CG_GROUPS`.
- **Decides:** whether to report this to Cosmoglobe (Trygve decides), and how to frame it in Paper A.
- **Status (2026-10-25):** BeyondPlanck: negligible (E013: 0.008 GHz, 0.8 conditional sd, because
  nu_p was one full-sky parameter at high S/N). DR1: not affected (E013 scoping). Open: which
  DR2-era runs used `sample_specind_mh` with CG updates (no parameter file in the repository does;
  the DR2 product page's link to `CG_param_c0001_DR2.txt` on tsih3.uio.no returns 404 as of
  2026-10-25, so this needs the file from Cosmoglobe or from Trygve's own access),
  and the size of the E015 bandpass effect (needs Commander3 on simulated TOD, cluster).

## Q15. Zero levels (monopoles) for maps-only Commander4 runs on real data
- **Why:** Commander4 has no per-band monopole sampling. E014 pass 2 and 3 mixed two zero-level
  conventions (WMAP/LFI: high-latitude median removed; HFI: referenced to the DR1 CMB, which has a
  +28.6 uK mean at |b| > 30, plus the DR1 dust template); the full-sky synchrotron index moved by
  0.43 (1250 posterior sd) between the two HFI sets and the dust index by 0.05. Any real-data
  index (Paper B) is meaningless until the zero levels are consistent, and the question of which
  reference (DR1 model, HI correlation, Commander3 `md` samples) is right is a systematic that must
  be propagated, not fixed by hand.
- **How:** (1) Pass 4 of E014: one convention for all bands, offset_b = mean over |b| > 30 of
  (map_b minus DR1 sky model at band b) with `nextgen/scripts/dr1_model_zero_levels.py`; compare
  with the DR1 chain `md` values read remotely; exact 2D grid (minutes) on pass-2, pass-3 and
  pass-4 inputs to measure the index sensitivity to the convention. (2) Decide whether Commander4
  needs a monopole component (per band, with the Commander3 priors: HI cross-correlation for HFI,
  monopole-dipole from a mask for WMAP/LFI) or whether offsets with a prior can be sampled as one
  extra parameter per band in the collapsed likelihood (the offset is linear, so it can be
  profiled or marginalized exactly like an amplitude: an "offset template" column in the mixing
  matrix). (3) Propagate the zero-level uncertainty into the index posterior (prior width per band).
- **Decides:** whether real-data spectral-index results from this fork can be reported at all;
  the design of the mixing-matrix offset column (ties into M2.4).
- **Status (2026-10-27):** sensitivity measured (E014 pass 4, exact grids on the same 20 maps):
  one DR1-model convention for all bands moves the ridge beta_s from -2.9103 (mixed conventions)
  to -2.8755 (260 posterior sd) and beta_d from 1.4150 to 1.3991 (460 sd); removing the 3.4 K
  non-solar Haslam dipole moves beta_s by a further +0.011 (80 sd). The convention is the
  dominant systematic of a full-sky index (F17). Pass 4 (`dr1_model_zero_levels.py`) is the
  working convention for E014-style inputs (D011); item (2) decided 2026-10-27 (D013): per-band
  offset columns with Gaussian priors in the collapsed per-pixel likelihood, implemented at M3.3;
  item (3) follows the implementation.

## Q16. How must the hierarchical amplitude prior and the dust SED be made robust to misspecification at high S/N?
- **Status 2026-10-30 (E023 stage 1, F26, F27, D022):** the driver is the template, not the SED:
  a template from one scaled frequency map biases the index field at 0.03 x and at 1 x noise,
  and no scalar or region-scale hyperparameter law repairs it for dust (pixel-scale error tied
  to T_d); fixed T_d is a 0.03 effect at most; the per-pixel CMB prior must be flat at realistic
  noise (leakage of the correlated CMB). Next (E025): component-separation-quality templates
  and the true-T dust SED by the same scans, then the first stage-2 chain; sampled T_d as the
  structural fix. At realistic noise the faint half of the sky is prior-dominated at nside 16.
- **Status 2026-10-30, later (E025, F28, D023):** a free template multiplier per nside-16 region
  makes the hier target close to unbiased for a smooth 30 percent template error at both
  depths and both components (within 0.02 at 0.03 x, insignificant at 1 x); the global law is
  biased for dust (0.07 to beyond 0.1). Fixed against true T_d changes nothing: sampled T_d is
  deprioritized. Next: the multiplier field in `AmpPrior`, then the stage-2 chain.
- **Why:** F21 (E020): with a data-built template wrong by factors 0.4 to 1.7 at high latitude
  under one global scatter s, and with T_d fixed, the collapsed hier posterior for beta_s runs
  away (to -26) and beta_d shifts by 30 Fisher sd in the plane. Real Planck plus WMAP intensity
  data at nside 64 are in this regime for dust and near it for synchrotron. F23 (E021 attempt
  1) adds the mildest case: a 30 percent template scatter whose conditional law differs from
  the prior's (g = a / r instead of a = g r) shifts the prior-dominated dust index by +0.13 at
  four stated sigma and the amplitude hyperparameters by 2 percent, with the data-dominated
  regions untouched; the determinant term of the collapsed marginal carries it.
- **Candidates:** a non-stationary scatter (s per latitude band or per region, or s as a smooth
  field); a heavy-tailed scatter (Student-t pseudo-band, which the collapsed Gaussian kernel does
  not give in closed form, so a scale-mixture with one latent scale per pixel sampled by Gibbs);
  an explicit template-error model (template = truth times a smooth field with its own prior);
  sampling or marginalizing T_d per region (a second spectral parameter in the region kernel).
- **How:** the E020 sky (kept) at realistic noise (noise scale 1) and at 0.03, each candidate
  against the same metrics as E021; a well-specified control from E021. Pre-register the
  acceptance: index RMSE within two bootstrap SE of the E021 result on the consistent sky, and
  no region with |beta_s - truth| above 0.5 on the PySM sky.
- **Depends on:** E021 (the clean baseline); M3 (real-data model) before any real-data claim.
- **Decides:** whether Paper A can run hier on real intensity maps or must restrict the template
  prior to polarization and low-S/N intensity.

## Q17. Is the BeyondPlanck XV polarized beta_s chain stationary, and what is its effective number of free amplitude elements?
- **Why:** E022 / F22. The released chains (four, 1001 samples, one shared start) have 23 to
  114 effective samples per chain and move towards flatter indices; on an nside-16 surrogate
  the flat-prior collapsed target is shifted by +0.27 (Spur) and runs to the prior boundary
  at nside 64, but BP's effective element count (7k to 4M per region) decides whether its
  own target is shifted or neutralized by the run's component priors. Paper A's statement
  about BP XV depends on it; so does the method recommendation (data-scaled priors on every
  amplitude).
- **How:** (1) the Commander3 semantics of `COMP_CL_TYPE = gauss` (D_l = 200, synchrotron
  and dust) and `none` (CMB): does the amplitude step apply S^-1 with those spectra, and how
  many effective elements does that leave per region (Trygve knows the code; otherwise read
  comm_Cl_mod.f90 and comm_diffuse_comp_mod.f90); (2) one conditional-width measurement: the
  index chi2 curvature at a released sample's amplitudes (BP_c0001_v2.h5 holds the
  amplitude alms per sample; a 1-D scan of beta_s for the Spur at fixed amplitudes gives
  sd_cond directly, 0.01 against 5e-4 is the discriminating range); (3) if sd_cond is below
  1e-3, the chains never approached their flat-prior target and the reported sds are not
  posterior widths; if it is near 0.01 the target is unshifted and the reported values are
  those of a stationary chain forced by the other Gibbs blocks.
- **Decides:** whether Paper A reports BP XV's polarized index uncertainties as
  underestimated or as consistent; Needs Trygve for the presentation either way (STATUS).
- **Status (2026-10-28, 13:50):** item 2 measured on sample 000500 of chain 1
  (`bp_sample_cond_width.py`): sd_cond 0.016 (Spur), 0.0085 (Plane), the wide end; the
  synchrotron amplitude is bound by the gauss C_l prior, the chains sit at their target, the
  amplitude block is 30 and 4 percent of the per-step motion, r 4 to 9, two-block tau 36
  (Spur) and 150 (Plane). Item 1 (the semantics of `none` for the CMB) remains; the
  measurement implies the CMB column is not a free per-pixel amplitude in BP's chain.

