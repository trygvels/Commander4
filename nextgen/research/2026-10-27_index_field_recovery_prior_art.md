# Prior art for spectral-index map recovery on PySM skies (input to E020 / M2.5 / Q3)

Lit-scout report, 2026-10-27. Read levels: "full" = arXiv HTML/ar5iv full text or PDF read;
"abs" = abstract page only. Section/table references are to the paper named. Items marked
(unverified) were not confirmed in a primary source. (inference) marks my own synthesis.

Scope of what was searched: arXiv HTML and ar5iv full texts of the papers below; the PySM 3
docs `models/` page and `presets.cfg` on GitHub; two arXiv API queries and two web searches for
GMRF / GP / random-field priors on index maps (listed in Section 6).

## 1. fgbuster multi-patch, multi-resolution and adaptive configurations

**arXiv:2510.08534, Rizzieri, Leloup, Errard, Poletti 2025, "Cleaning Galactic foregrounds with
spatially varying spectral dependence from CMB observations with fgbuster"** (full)
- Setup (Sect. 4.1): PySM 3 d1s1 and d10s5 (plus d0s0 reference), CMB c1 with r = 0, simulated
  directly at nside 64 (not downgraded, to keep per-pixel MBB/power-law exact). LiteBIRD-like
  (15 bands 40 to 400 GHz, 2.16 uK-arcmin combined) and PICO-like (21 bands 21 to 799 GHz,
  0.86 uK-arcmin). Delta bandpasses, no beams, f_sky = 60 %, white noise, 100 noise realizations.
  Fit bounds (Sect. 4.3): beta_d in [0, 5], T_d in [10, 40] K, beta_s in [-6, 0]. d12s7 not used.
- Configurations (Sect. 4.2): multi-patch "nsN" (one HEALPix nside for all three parameters,
  ns2/16/32/64; ns64 gives numerically zero systematics by construction); multi-resolution
  "mres" nside = [64, 8, 4] for [beta_d, T_d, beta_s]; "adaptive" = latitude-dependent mres from
  the LiteBIRD 2023 forecast; "ideal" subsets = true template binned into n_bins equal-width bins
  (200 bins gives [95, 147, 126] subsets for [beta_d, T_d, beta_s], Fig. 16 caption);
  "realistic" subsets = bins of a first-pass recovered map.
- Metrics: systematic residual = CMB weights at the fitted parameters applied to the noiseless
  foreground maps, delta s_syst = W^CMB(beta_bar) d (Sect. 3.1, 4.4), power spectrum via anafast;
  statistical residual = first-order propagation of delta beta from the Fisher matrix, Sigma^-1
  (Eq. 37), drawn by sparse Cholesky (Eqs. 41, 42), plus post-separation noise
  N_p = [(A^T N^-1 A)^-1]_CMB (Eq. 43); total C_l = C_l^syst + C_l^stat (Eq. 48); r likelihood
  Eq. 45. Index-map error is decomposed per pixel as delta beta_tot = delta beta_stat +
  delta beta_syst (Sect. 3.1, Fig. 7). No single RMS summary of the index map is reported.
- Headline numbers (Table 1, Sect. 4.6; r_bias +- sigma_r in units of 1e-3, simulation-based):
  LB d1s1 adaptive 0.8 +- 1, 200 ideal 0.0004 +- 0.2; LB d10s5 adaptive 3 +- 1, 200 ideal
  0.003 +- 0.2; PICO d1s1 adaptive 23 +- 0.7; PICO d10s5 adaptive 20 +- 0.7. Pearson correlation
  of beta_d and T_d statistical residual maps: -0.995 (ns16) vs -0.805 (mres) (Table 2).
- Relevance to E020: this is the reference definition of "fixed ladder [64, 8, 4]" and of the
  systematic vs statistical residual split. Reuse the delta s_syst definition and report
  C_l^syst and C_l^stat of the CMB map separately. Note they never report index-map RMS, so
  our index-error metric must be defined by us (inference).

**arXiv:1811.00479, Errard & Stompor 2019, PRD 99, 043529** (full)
- Setup (Sect. II.3): d1s1 templates smoothed/degraded to nside_beta 2 to 64; maps at nside 64;
  15 bands 40 to 402 GHz, ~1 deg beams, 4 to 38 uK-arcmin white noise, 60 % HFI dust mask,
  r_in = 0.001, 10 CMB+noise realizations. Sim#2 uses a2d6s3 with 2 % polarized AME.
- Metrics: residual Delta_CMB (Eq. 9), its spectrum C_l^Sigma (Eq. 14), Fisher Sigma (Eq. 7),
  r likelihood (Eq. 23). Bias = offset of the r peak from 0.001.
- Numbers: ignoring spatial variation biases r by ~0.01 (abstract, Sect. III.2 at
  nbar_beta = 4); with the residual term included the fit is unbiased with sigma(r) from
  ~4e-4 (nbar = 2) to ~7e-4 (nbar = 16) (Sect. III.1, Fig. 7). Patch counts 48/192/768/3072
  for nside 2/4/8/16 are my arithmetic.
- Relevance: origin of the HEALPix-patch "nside ladder"; shows statistical residual grows with
  patch count while bias falls. Use the same trade-off plot (bias and sigma vs nside_beta).

## 2. Clustering-based patches

**arXiv:2109.11562, Puglisi et al. 2022, MNRAS (doi 10.1093/mnras/stac069)** (full)
- Setup (Sect. 4.1): PySM d1s1 at nside 64, clustering at nside 32. LiteBIRD 15 channels
  40 to 402 GHz, f_sky 60 %; SO-SAT 6 channels, f_sky 27 %. White noise, 20 signal+noise
  realizations, noise-free runs for bias.
- Clustering (Sect. 3.1, 3.2): spectral clustering on beta_s, beta_d, T_d maps with per-pixel
  sigma; affinity exp[-(1-Theta)^2 / 2 sigma_pix^2] with sigma_pix ~ 46 arcmin, deformed by a
  KS-test similarity with weight alpha (0 < alpha <= 0.5); agglomerative cut at threshold
  delta; objective W = sqrt(V_w^2 + V_b^2) (Eq. 6). Median pixels per cluster ~44 (beta_d),
  ~88 (T_d), ~16 (beta_s), angular sizes ~3.5, 7, 14 deg (Sect. 4.1). ~4000 CPU-h at nside 32.
- Index metric (Fig. 10, Sect. 4.1): median relative error 0.002 (beta_s), 0.006 (beta_d),
  0.017 (T_d); worst errors <= 20 % in T_d. This is a median relative error, not RMS.
- r (Table 2) with template mis-modelling fraction f: LiteBIRD (0.22 +- 0.66)e-3 at f = 0 to
  (2.18 +- 0.68)e-3 at f = 1; SO-SAT (0.08 +- 1.21)e-3 to (0.20 +- 1.19)e-3. HI cloud-count
  patches (Sect. 4.2, K = 1481 at nside 32): r = (0.09 +- 0.26)e-3 (LiteBIRD).
- Relevance: a clustering baseline for E020 and the only paper here that reports an index-map
  error statistic (median relative error). Their "f" test is a cheap way to probe prior
  misspecification.

**arXiv:2604.08463, Kabalan et al. 2026, FURAX K-means patches** (full)
- Setup (Sect. 4, App. B.1): c1d1s1 at nside 64, LiteBIRD-like white noise (levels cited to
  LiteBIRD 2023, not printed), no beams, 40 noise realizations, 3 latitude regions of 20 % sky.
- Patches: spherical K-means on pixel positions only (equal-area; the text does not cluster on
  beta maps). K for each parameter from a 24-value grid (50 to 10000), 24^3 configurations.
  Selection minimizes r_hat + sigma(r) (Sect. 2.4) where r_hat uses the true CMB and noise-free
  foreground leakage (Eqs. 14, 15), so selection needs simulation truth.
- Metrics: sigma^2_CMB map variance (Eq. 10); r from a Wishart-type likelihood l = 2 to 150
  (Eq. 12); systematic, statistical and total BB residuals (Eqs. 14 to 16). No index-map error
  statistic; index maps shown only visually (Fig. 9).
- Numbers: r_hat = 0.9e-4, 68 % interval [-2.7, +6.5]e-4; 68 % upper limit 7.4e-4 vs 1.0e-3 for
  the mres baseline (Sect. 4.3); r = 3e-3 recovered as 2.9 (+0.37/-0.34)e-3 (Sect. 4.4).
  Optimal K (Table 2): beta_d 7000 to 9685, T_d 300 to 2500, beta_s 150 to 300. Baseline
  ladder (Table 3, from LiteBIRD 2023): beta_d nside 64 everywhere; T_d nside 8/4/0 for
  low/mid/high latitude; beta_s nside 4/2/2. ~640 GPU-hours for the grid search.
- Relevance: the "optimal patch count" ordering (beta_d finest, beta_s coarsest) is the
  empirical target a hierarchical field should reproduce without a grid search. Their Table 3
  ladder is the latitude-adaptive version of [64, 8, 4].

## 3. Commander3 / BeyondPlanck / Cosmoglobe practice

**arXiv:1502.01588, Planck 2015 X** (full, first 200k chars)
- Intensity: baseline fit at 1 deg FWHM, nside 256; beta_d and T_d both free per pixel there
  (Table 5). High-resolution 7.5 arcmin, nside 2048 fit keeps T_d fixed to the 1 deg solution
  and refits beta_d per pixel (Sect. 5.4). Priors (Table 4): beta_d ~ N(1.55, 0.1),
  T_d ~ N(23, 3) K, T_e ~ N(7000, 500) K with 2 deg smoothness, nu_p ~ N(19, 3) GHz. "We adopt
  a Jeffreys prior for the spectral parameters in order to suppress prior volume effects."
- Output statistics (Table 5, 1 deg): beta_d max-post 1.53 +- 0.05, mean 1.51 +- 0.06, rms
  0.05 +- 0.03; T_d 21 +- 2 K, rms 1.1 +- 0.7 K. These are the statistics of the map that PySM
  d1 uses (see Section 5).
- Polarization: spectral parameters fixed to the intensity result; 40 arcmin / nside 256 and
  10 arcmin / nside 1024 products.
- Relevance: gives the d1 truth map's own posterior scatter (rms 0.05 in beta_d), i.e. a
  floor on how much structure in d1 beta_d is real (inference).

**arXiv:1807.06208, Planck 2018 IV** (full, first 200k chars; Sect. 5.3 not reached)
- Sect. 2.1: polarized beta_s is a single spatially constant value, beta_s = -3.1 +- 0.1
  (abstract); polarized beta_d fitted per 3 deg FWHM region, beta_d = 1.55 +- 0.05 (abstract);
  T_d fixed from intensity. A web snippet claims beta_s was estimated at 5 deg smoothing
  (unverified; contradicts the Sect. 2.1 reading, which says spatially constant).

**arXiv:2011.08503, BeyondPlanck XV (Svalheim et al.)** (full)
- beta_s in 4 disjoint regions (High Latitude, Spur, Centre, Plane; merged from the 24-region
  Fuskeland 2014 partition). Proposed region map smoothed with a 10 deg FWHM Gaussian before
  evaluating the mixing matrix (Sect. 4.3). chi^2 evaluated at nside 16 at each channel's native
  resolution (Eq. 26). Prior N(-3.3, 0.1) (Sect. 4.4); width 0.2 let beta_s reach -2.7 and
  produced foreground excess. Metropolis-within-Gibbs with tuned Gaussian proposal, acceptance
  0.2 to 0.8, correlation length 30 to 50 samples (Sect. 5.3.2).
- Results: Spur -3.17 +- 0.06, Plane -3.03 +- 0.07; other two regions prior dominated. Prior
  mean shift of 0.2 moves the Spur by <= 0.03 (Fig. 11). beta_d single full-sky 1.62 +- 0.04
  with prior 1.56 +- 0.10 (Sect. 5.3.3). No simulation recovery test of beta_s is in the paper.
- Relevance: the "smooth the region map with a wide beam" trick is a fixed-prior substitute for
  a correlated field; ringing/clipping in Commander3 (research/01) is the known failure mode.

**arXiv:2310.13740, Cosmoglobe DR1 III** (abs): single full-sky beta_s = -3.07 +- 0.07 from
T-T plots. **arXiv:2507.05324, LiteBIRD in Commander3** (full, PDF pp 1 to 14): sky is PySM
d1s1a1f1co1 (Sect. 4.1); spectral indices "kept at their true input values for simplicity"
(Sect. 3.1), so no index recovery test exists there.

Published use of Commander3's alm index sampler (ALMSAMP_PRIOR_FWHM): not found. Searches:
web "Commander spectral index alm sampling smoothing prior FWHM Cosmoglobe DR2 / Planck 2018 /
BeyondPlanck" (results: 0903.4311, 1103.2554, 2310.13740, 1807.06208, Gjerlow FIR2025 slides,
none describe an alm sampler). A FIR 2025 talk (Gjerlow) says DR2 fits per-pixel beta_d and
T_d in a one-MBB model (slides, unverified).

## 4. GMRF / GP / hierarchical priors on index maps

**arXiv:1910.14170, Grumitt, Jew, Dickinson 2020, MNRAS (doi 10.1093/mnras/staa1857)** (full)
- Setup (Sect. 3, Table 1): PySM s1 + d1 (Q, U only), CAMB lensed CMB with r = 5e-3, all maps
  smoothed to 70 arcmin, nside 64. C-BASS 5 GHz (4320 uK-arcmin) + 15 LiteBIRD channels
  (40 GHz 27.9, 402 GHz 19.6 uK-arcmin); white noise, Gaussian likelihood (Eq. 10). Input
  index maps: beta_s mean -3.0, sd 0.06; beta_d 1.54, sd 0.04; T_d 20.9 K, sd 2.2 K.
- Regions: mean-shift on pixel coordinates plus two naive two-frequency indices (Eqs. 4, 5),
  zeta = 0.5, bandwidth 0.3, min 10 pixels: 171 regions of 138 to 645 pixels (mean 290).
- Hierarchical model (Sect. 5.4, Eqs. 19 to 21, Table 3): per pixel beta_s^p ~ N(mu, sigma)
  with hyperpriors mu_bs ~ N(-3, 0.3), sigma_bs ~ HalfNormal(0.35); mu_bd ~ N(1.6, 0.3),
  sigma_bd ~ HalfNormal(0.35); mu_Td ~ N(21, 2), sigma_Td ~ HalfNormal(2.5); non-centred
  (beta = mu + Gamma sigma, Gamma ~ N(0, 1)); flat amplitude priors. The pooling model uses
  N(-3, 0.3), N(1.6, 0.3), N(21, 2) times Jeffreys. No spatial correlation inside a region.
- Metrics (Eq. 24): normalized deviation eta = (input - output)/sigma per pixel, compared to
  N(0, 1) via median and MAD (Table 4, Fig. 3); R_hat <= 1.1; n_eff ~ 1000; NaMaster pseudo-C_l
  with 60 % mask, cross-split; r from a Gaussian likelihood 30 <= l < 180 (Eq. 26).
- Numbers: r = (12.9 +- 1.4)e-3 CP(L), (9.0 +- 1.1)e-3 CP(LC), (5.2 +- 1.0)e-3 H(LC). Table 4
  (Q, median / MAD of eta): beta_s CP(LC) -1.03 / 16.8 vs H(LC) 0.09 / 1.46; beta_d -0.35 /
  4.96 vs -0.35 / 1.42; T_d 1.23 / 5.11 vs 0.63 / 1.96. Dust beta residuals ~50 % lower and
  T_d ~25 % lower for H(LC) (Fig. 8). NUTS ~1e3 x MH in ESS per time for pooling (Sect. 7).
  H(L) (LiteBIRD alone) was not run: divergences, convergence "suspect".
- Relevance: the closest published analogue of our hierarchical model, minus spatial
  correlation. Their eta median/MAD table is a ready-made coverage metric; adopt it alongside
  SBC/TARP. Note MAD of eta is not an RMS error in beta.

**arXiv:2205.13968, de Belsunce, Gratton, Efstathiou 2022** (full, PDF)
- Prior (Eqs. 4, 5): C_beta,ij = sigma(beta_i) sigma(beta_j) exp(-theta_ij^2 / 2 theta_corr^2),
  theta_corr = 5 deg fixed by choice (Sect. 2.1.2), sigma(beta_s) = 2.0, sigma(beta_d) = 1.0,
  means -3.0 and 1.5. Newton-Raphson MAP, not sampling. Real Planck 2018 LFI+HFI (SRoll2) Q, U
  at nside 16 (3072 pixels), f_sky 0.52; no simulation recovery test.
- Numbers (Table 1, corrNCM+corrPrior, SRoll2): beta_s = -2.833 +- 0.620, beta_d = 1.429 +-
  0.236 (mean and sd across pixels). With tight priors (sigma 0.1) they recover
  -3.000 +- 0.008 and 1.494 +- 0.028 and state "previous results are strongly prior dominated".
- Relevance: the only Gaussian-kernel correlated prior on a beta map in CMB component
  separation that I found; correlation length fixed, not inferred; no truth comparison.

No paper inferring a correlation length, power spectrum or GMRF precision for a spectral index
map inside a CMB component-separation sampler was found (see Section 6 for the queries run).
This agrees with research/01 line 37.

## 5. PySM 3 model definitions for the truth maps

Source: `presets.cfg` (galsci/pysm main) and the docs `models/` page (both read), Thorne et
al. 2017 arXiv:1608.02841 (full, ar5iv) and PEGS arXiv:2502.20452 (full, HTML, first 200k).

- **d1** (`ModifiedBlackBody`): I at 545 GHz, Q/U at 353 GHz, `pysm_2/dust_beta.fits` and
  `pysm_2/dust_temp.fits`, uK_RJ. Templates from Planck 2015 Commander; intensity nside 2048
  degraded to 512; polarization smoothed to 2.6 deg (docs) or 2 deg (Thorne Sect. 2.2.1; the
  two values are not reconciled) with Gaussian small scales added. beta_d and T_d are the
  Planck 2015 X Commander maps (1 deg, nside 256, Table 5 stats above: beta_d 1.51 +- 0.06
  rms 0.05, T_d 21 +- 2 K rms 1.1 K). Same index in I and P. Grumitt Table 1 quotes beta_d
  1.54 sd 0.04, T_d 20.9 K sd 2.2 K for the d1 maps at nside 64. Exact stored nside of
  dust_beta.fits not stated in the docs (unverified; 512 by analogy with the template).
- **s1** (`PowerLaw`): I = Haslam 408 MHz (Remazeilles 2015 reprocessing, 56 arcmin),
  freq_ref_I 408 MHz; Q/U = WMAP 9-yr 23 GHz smoothed to 3 deg (Thorne) or 5 deg (docs);
  `pysm_2/synch_beta.fits` = Miville-Deschenes et al. 2008 "Model 4" from Haslam + WMAP 23
  GHz with a magnetic-field model; mean -3.00, sd 0.06 (Thorne Sect. 2.1.3); resolution not
  stated in Thorne, PEGS says the Miville-Deschenes map is about 7 deg (Sect. III.2.3).
- **d10** (`ModifiedBlackBody`): GNILC PR2 353 GHz template at 21.8 arcmin uniform
  resolution, CIB monopole 0.13 MJy/sr removed, freq_ref 353 GHz, available nside 2048/4096/
  8192, lower nside via ud_grade from 2048. beta_d and T_d = GNILC PR2 maps (Planck Int. XLVIII,
  arXiv:1605.09387, abs: full-sky means T = 19.4 +- 1.3 K, beta = 1.6 +- 0.1) plus injected
  small-scale fluctuations: power-law fits alpha = 0.04 (beta_d, 200 < l < 400) and
  alpha = -0.47 (T_d, 100 < l < 400), l_1 = 400, modulated by the intensity modulation map,
  "below 10 %" relative and a few percent at high latitude (PEGS Sect. III.2.2, Table 2).
  Caveat: ~3 % of sky inside GAL097 replaced by the input GNILC map (docs). Native FWHM of the
  GNILC beta_d/T_d maps not found (unverified; the GNILC products are variable resolution).
- **s5** (`PowerLaw`): I and P templates `synch_template_nside{nside}`, freq_ref 23 GHz;
  `synch_beta` = Miville-Deschenes map rescaled to the S-PASS dispersion: subtract mean -3.00,
  multiply by 1.572, add mean back (PEGS Sect. III.2.3); small scales from a power-law fit
  over 10 < l < 36 with alpha = -0.61, modulated by the synchrotron intensity. Brightest 7 %
  of sky replaced by input maps (docs). No printed mean/sd after rescaling; sd ~0.06 x 1.572
  = 0.094 at large scales is my arithmetic (inference).
- **s7** (`CurvedPowerLaw`): same template and beta as s5 plus `synch_curvature`, freq_curve
  23 GHz; large-scale c_s = Haslam smoothed to 5 deg linearly rescaled to span -0.076 to
  -0.041 (ARCADE patch mean -0.0517, sd 0.0054), small scales alpha = -0.61 (PEGS
  Sect. III.2.3). Note freq_ref_I = 23 GHz, unlike s1 (408 MHz).
- **d12** (`ModifiedBlackBodyLayers`): 6 MKD layers (arXiv:1706.04162), each with amplitude,
  beta_d^(k), T_d^(k) at nside 2048, MJy/sr, color_correction 0.911, max_nside 2048. Layer
  spectral maps are random per layer and "do not match the observed beta_d and T_d maps"
  (PEGS Sect. IV.1). There is no single truth beta_d map: any one-MBB "truth" must be defined
  by us, e.g. the best-fit effective index at our bands (inference).
- PEGS composite skies (Table 7): Low = d9 s4 f1 a1 co1; Medium = d10 s5 f1 a1 co3;
  High = d12 s7 f1 a2 co3. This confirms the mapping flagged "(verify)" in nextgen/DATA.md.

## 6. Gaps, contradictions, and suggested follow-ups

- No published index-map RMS or coverage benchmark on PySM skies exists for nside 16 to 64 in
  intensity with Planck bands. All the papers above are polarization-only (Q, U) at nside 64
  with LiteBIRD/PICO/SO bands. Our intensity setup has no direct comparator; we should still
  report Rizzieri's C_l^syst / C_l^stat split and Grumitt's eta median/MAD so numbers can be
  placed next to theirs (inference).
- Queries run that returned nothing relevant: arXiv API `abs:"spectral index" AND
  abs:"component separation" AND (abs:"Gaussian process" OR abs:"Markov random field" OR
  abs:"Gaussian random field" OR abs:hierarchical)` (0 results); arXiv API `("spectral index"
  OR "spectral indices") AND foreground AND ("Gaussian process" OR "random field" OR
  "smoothness prior" OR "correlation length") AND cat:astro-ph.CO` (2 results: 2205.13968 and a
  21 cm paper 2211.10448); web search for GMRF/GP index priors (only 2205.13968 and 21 cm GP
  work 2407.11296). 
- LITERATURE.md check: entries for 2510.08534, 2604.08463, 2109.11562, 2205.13968, 1910.14170,
  1811.00479 are consistent with the full texts. One web snippet described de Belsunce's
  correlation length as "inferred"; the paper fixes it at 5 deg (Sect. 2.1.2), so LITERATURE.md
  is right.
- Unverified items to close: GNILC beta_d/T_d native resolution (read 1605.09387 full text);
  stored nside of `pysm_2/dust_beta.fits` and `synch_beta.fits` (inspect the FITS headers);
  Planck 2018 IV Sect. 5.3 details for beta_s in polarization; Kabalan first half for any
  index-map statistic (not found in the second half).
- Suggested follow-ups: Rizzieri et al. 2025 beams paper arXiv:2409.12944 for mixed-beam
  treatment at the beta resolution; LiteBIRD 2023 forecast (the source of the latitude-adaptive
  ladder, not read here); Martinez-Solaeche et al. 2018 arXiv:1706.04162 for the d12 layer index
  distributions.
