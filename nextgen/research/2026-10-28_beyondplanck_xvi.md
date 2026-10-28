# BeyondPlanck XVI (Herman et al.), polarized AME limits: facts for spectral-index sampling

Date: 2026-10-28. Author: lit-scout agent for Trygve Leithe Svalheim.

Paper: arXiv:2201.03530v2, "BeyondPlanck XVI. Limits on Large-Scale Polarized Anomalous Microwave Emission from Planck LFI and WMAP", D. Herman, B. Hensley, K. J. Andersen et al., A&A 675, A15 (2023). Title verified on the arXiv abstract page. Read level: full text via arXiv HTML and ar5iv (all sections 1 to 6, figure captions, Eqs. 1 to 18), extracted through an HTML-to-text summarizer instructed to quote verbatim. No tables exist in the paper. Internal inconsistencies noted below are in the rendered text, not extraction gaps. The code is a separate Gibbs sampler, footnote in Sect. 2 pointing to github.com/hermda02/dang, not Commander (inference: the text never names Commander).

## 1. Spectral parameters, priors, spatial model

Polarization-only model, Eq. 3: synchrotron plus AME plus thermal dust. CMB is excluded from the model (Sect. 3.2, "known to be small"), tested in Sect. 5.4. Free-free: not stated.

- Synchrotron beta_s: Eq. 4, power law with nu_0 = 30 GHz, index beta_{s,p} "per pixel, but common for Q and U" (Sect. 3.2); amplitudes a_s per pixel and per Stokes parameter. Pixelization Nside 64 at 1 degree FWHM (Sect. 2). Nominal prior Gaussian beta_s = -3.1 +- 0.1 "following the results of Svalheim et al. (2022b)" (Sect. 4.3). Sect. 5.2 then treats the prior mean "as a free hyperparameter" and repeats the analysis for a grid of prior means from -3.3 to -3.0 (Sect. 5.3, Fig. 7); the grid spacing and the prior width at non-nominal means are not stated. No regions, no proposal smoothing, no curvature (Sect. 6: curvature "has not been addressed").
- Dust beta_d: one full-sky value, drawn from the Gaussian prior beta_d = 1.62 +- 0.04 each Gibbs iteration, no likelihood term (Sect. 3.2, Eq. 15, Sect. 4.1). T_d fixed at 19.6 K, scalar (Sect. 3.2). Dust amplitude fixed to the Planck DR4 353 GHz Q/U map, no CMB subtraction (Sects. 2, 3.2, Eq. 5).
- AME: no SED parameter. AME is assumed perfectly correlated with 353 GHz polarized dust; one full-sky template amplitude a_{AME,nu} per band, Eq. 6, jointly in Q and U; fixed to zero at WMAP V and LFI 70 GHz to avoid "a perfect degeneracy between the spatial synchrotron template and the AME template amplitudes" (at most N_band - 1 free, Sect. 3.2). Conversion to polarization fraction: Eqs. 8 to 11, using the Planck 2015 AME intensity map, p_353^max = 22% and <I_353/A_AME> = 2.5 (+0.2, -0.7) (Sect. 3.3).
- Intensity: no intensity fit in this paper.

## 2. Sampling algorithm

Gibbs chain, Eqs. 13 to 15: (a_AME, a_s) from P(a | d, beta_s, beta_d) by CG on the linear system after Eq. 17 with a fluctuation term (Sect. 4.2); then beta_s from P(beta_s | d, a_AME, a_s, beta_d) (Eq. 14), i.e. conditional on amplitudes, not marginalized; then beta_d from its prior. Amplitude priors: none stated (Sect. 4.2). beta_s step: Metropolis with a Gaussian proposal of width sigma "determined by the standard deviation of the Gaussian prior distribution", acceptance ratio Eq. 18 (likelihood ratio times prior ratio), "This procedure is repeated until convergence" (Sect. 4.3). Number of Metropolis steps per Gibbs iteration: not stated. Whether the Metropolis chain runs per pixel independently or on a global chi-square: not stated (Sect. 3.2 writes beta_{s,p} per pixel, Sect. 4.3 writes a single beta_j without pixel index). Noise covariance N_nu: "Gaussian distributed with a covariance matrix" (Sect. 4.2); diagonal or dense not stated. Between Gibbs samples the LFI maps are replaced by a different one of 1000 BeyondPlanck realizations to marginalize systematics (Sects. 2, 5). TOD parameters are not sampled here.

Mixing, acceptance rate, correlation length, trace plots: not stated anywhere. Convergence criterion: not stated.

## 3. Chain lengths

"For each chain, 500 full Gibbs samples are generated" (Sect. 5 intro). Number of chains, burn-in, thinning: not stated. Gelman-Rubin, ESS: not stated. Whether each beta_s prior mean on the Fig. 7 grid is a separate run: implied by "we repeat the above calculations for a grid" (Sect. 5.3) but not explicit.

## 4. Posterior widths and degeneracies

- beta_s posterior: Fig. 3 mean map only; no numerical mean or sd, no sd map (Sect. 5.2). "beta_s is very close to the prior mean over large portions of the sky", deviations along the Northern Galactic Spur, Fan Region, Galactic Plane "in agreement with Svalheim et al. 2022b" (direction not stated). "even the combination of LFI and WMAP has very little constraining power compared to the prior" (Sect. 5.2).
- AME versus beta_s: nominal prior -3.1 +- 0.1 gives posterior peak below 0% with p_AME^max < 1.0% (95%); for beta_s <~ -3.1 "typical upper limit" p_AME^max <~ 1.0%; for beta_s = -3.0 a "nominal detection" p_AME = 2.5 +- 1.0% (95%) (Sect. 5.3, Fig. 7). Mechanism: "synchrotron emission is oversubtracted by the flat spectral index" and "countered by adding a spurious polarized AME component" (Sect. 5.3). Inconsistency: the abstract and Sect. 6 quote 0.6% for beta_s <~ -3.1, and Sect. 6 quotes 3.5% for beta_s = -3.0; Sect. 5.3 gives 1.0% and 2.5 +- 1.0%. The paper does not reconcile these.
- CMB: adding synfast LCDM fluctuations "increases the upper limits on p_AME^max by 10 to 15%", so the main limits are "conservative" (Sect. 5.4, Fig. 8). Number of realizations: not stated.
- Synchrotron amplitude agrees with BeyondPlanck and Planck DR4 "to around the 3 uK_RJ level" (Sect. 5.2, Fig. 5); dust-synchrotron degeneracy along the plane (Fig. 5 caption).
- Monopoles, zero levels, dipoles: no occurrence of these words in the paper (polarization only).

## 5. Data

Polarization Q/U at 1 degree FWHM, Nside 64, uK_RJ, full bandpass integration (Sect. 2): BeyondPlanck LFI 30, 44, 70 GHz (1000-sample ensemble); WMAP K 22.8, Ka 33.0, Q 40.6, V 60.8 GHz (W excluded for "a high level of systematic residuals"); Planck DR4 353 GHz as dust template; Planck 2015 AME intensity map for the fraction conversion. Mask: Svalheim et al. central Galactic mask plus compact high chi-square regions including Tau A, f_sky = 91% (Sect. 2). Figure 1 and 5 maps smoothed to 3 degrees for display.

## 6. Limitations and recommendations

Sect. 6: data "have limited constraining power with respect to the spectral index of polarized synchrotron emission"; "not able to robustly and simultaneously constrain polarized AME and synchrotron emission"; LFI-only analyses "susceptible to even stronger degeneracies"; "It is conceivable that the steep spectral index they observe is partially caused by polarized AME"; curvature not addressed; "The only way to actually break these degeneracies is through additional high signal-to-noise observations" (C-BASS, S-PASS, QUIJOTE).

## Relation to repo notes

LITERATURE.md line 12 (0.6% vs 2.5 +- 1.0%) matches the abstract; add that Sect. 5.3 gives 1.0%. research/01_spatially_varying_parameters.md line 71 ("steeper than -3.1") matches Sect. 6. New: BP XVI samples beta_s per pixel at Nside 64 (not the four BP XV regions), conditional on amplitudes, prior-dominated, with no mixing diagnostics.

Not found despite searching the full text: Metropolis step count, chain count, burn-in, acceptance rate, correlation length, beta_s sd values, noise covariance form, grid spacing of prior means.
