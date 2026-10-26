# Open questions

Ranked by how much the answer changes what we build next. Each entry: why it matters, how to answer
it, and what the answer decides. Mark an entry RESOLVED with the experiment ID when answered.

## Q1. Does the legacy acceptance matter on real-data configurations? (ANSWERED at nside 64, E014)
- **Answer (2026-10-26):** yes. On DR1 WMAP + LFI, NPIPE 353 / 857, PR3 100 and Haslam at nside 64
  with the dx12 mask, the legacy `resolve_chisq` index posterior is 16 to 17x too wide and the dust
  index is shifted by 3 ridge sd (E014, exact grid reference). The scaling of F1 predicts factors
  above 50 at nside 512 to 1024. Remaining: the `conditional` / `marginal` chains (prior effect,
  running), IQU, and a direct run at nside 128 to confirm the scaling (open, lower priority).
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

## Q3. What is the optimal spatially varying index model?
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

