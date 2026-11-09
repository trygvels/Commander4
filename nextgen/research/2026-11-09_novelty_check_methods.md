# Novelty check for the methods paper: six claims against the literature to 2026-11-09

Date: 2026-11-09. Author: lit-scout agent for Trygve Leithe Svalheim.

Question: for each of six method claims, what is the closest published prior art, and is the claim
new, partly known, or known? Special attention to 2025 and 2026.

Read levels used below:
- FT: full text read through an HTML rendering (ar5iv or arXiv HTML) via a text summarizer.
- FT-part: as FT, but only the stated character range was read (pages over 100 000 characters).
- ABS: abstract only (arXiv API). BIB: title and authors verified only, content not read.
- NOTE: taken from an existing note in `nextgen/research/`, which records its own read level.

Caveat that applies to every quote and number: they were extracted by a summarizer instructed to
quote in pieces under 120 characters. Greek letters, minus signs and ranges are transliterated to
ASCII. Verify against the PDF before citing. Statements marked (inference) are mine, not a paper's.

## 0. Short answer

| # | Claim | Verdict | Closest prior art |
|---|---|---|---|
| 1 | Re-solve amplitudes per proposal, accept on chi-squared, is not a valid sampler | Partly known | 0804.2645, 2201.08188, 1608.00551 |
| 2 | Flat-amplitude-prior marginal has a volume term that biases low-S/N indices; coverage threshold | Partly known (mechanism and remedies known; coverage threshold not found) | 0811.4280, 1608.00551, 0709.1058 |
| 3 | Hierarchical, template-tied amplitude prior as the remedy, validated by SBC | Partly known (new in CMB component separation as far as searched) | 2201.08188, 1203.0025, 1907.11426, 2211.10448 |
| 4 | GMRF or SPDE prior on the index map, hyperparameters inferred, NUTS inside Gibbs or collapsed | New for CMB component separation as far as searched; partly known across fields | 2606.04094, 1401.4711, 2011.11575 |
| 5 | Conditional Gibbs mixes slowly at high S/N; collapsed sampling removes it | Partly known (qualitatively known, cure published; S/N-resolved ESS not found) | 2201.08188, 1608.00551, 2509.11894 |
| 6 | Zero-level convention sets large-region indices at many sd | Partly known (degeneracy and treatments known; "many sd" statement not found) | 0709.1058, 1411.7616, 2606.21334 |

No 2025 or 2026 paper was found that overlaps substantially with the proposed methods paper. The
whole arXiv listing for `abs:"component separation"` from 2025-01-01 to 2026-11-09 (165 entries)
was scanned by title and summarizer gist. Partial overlaps are in Sect. 3.

## 1. Claim by claim

### Claim 1. Validity of the index step (re-solved amplitudes, chi-squared acceptance)

Verdict: partly known. The profile likelihood, its difference from the flat-prior marginal, and the
valid collapsed step are all published. Not found anywhere: an analysis of the validity of a
Metropolis step that re-solves or redraws amplitudes for every proposal and accepts on chi-squared,
a statement of what such a step targets, or a quantification of the width errors (project F1).

Closest papers:
- 0804.2645, Stompor, Leach, Stivoli, Baccigalupi 2009 (FT, both halves). Sect. 3.1.1: the
  spectral likelihood, Eq. (9), follows from substituting the amplitude solution Eq. (8) into
  Eq. (5); it "is not a proper likelihood in a common sense, and is thus often referred to as a
  'ridge' or 'profile' likelihood". Eq. (10) is the flat-prior marginal with the extra term
  ln|(A^t N^-1 A)^-1|; "the marginal likelihood peaks at their values, which do not coincide with
  their maximum likelihood estimates". The authors prefer the profile form and say it gives
  "useful, if not strictly precise, estimates of their uncertainty". Sect. 5: maximized "direct or
  via MCMC sampling". No discussion of Gibbs or Metropolis validity.
- 2201.08188, BeyondPlanck XIII (FT-part, first 100 000 of about 167 000 characters this session;
  the 2026-10-27 repo note records a full read). Sect. 4.1, Eqs. 19 to 24: P(a, beta) = P(beta)
  P(a|beta); "a short (typically a few hundred steps) standard Metropolis sampler" on the marginal
  Eq. (23); "Immediately following the last Metropolis step, we draw one sample from P(a|m_nu,
  beta) using Eq. (24)"; "it is critically important that no other parameters are updated between
  beta and a". The ridge option "sets a equal to its maximum likelihood value for a given value of
  beta" and is "identical to the above expression, with the exception of the last determinant term
  being excluded"; "We implemented support for both options in our codes". The paper does not say
  that the ridge option changes the target, does not say which option produced the released
  results, and never discusses a chi-squared step with refitted amplitudes.
- 1608.00551, Alonso, Dunkley, Naess, Thorne 2017, BFoRe (FT, three passes). Sect. II.3.1 and
  App. A.2, Eq. 29: indices are sampled from the amplitude-marginal. App. B: a "volume prior"
  |N_T|^(-1/2) is applied to cancel the determinant. (inference) Marginal times volume prior is
  algebraically the Stompor profile likelihood times the index prior, so BFoRe's production
  configuration samples the ridge likelihood. The paper does not say so and never uses the words
  profile likelihood.

Other relevant: 0709.1058 (FT-part), indices sampled conditional on amplitudes, a valid Gibbs
conditional (Sect. III.2, Eq. 17). 2405.18365, Morshed et al. 2024 (FT-part, first 100 000
characters): profile likelihood Eq. (7); "This likelihood may be significantly biased and requires
a special correction procedure"; foreground amplitudes maximized, not sampled; (inference) the bias
corrected there comes from the CMB-covariance term of the minimally informed method, not from the
flat-prior volume term. 2511.04546, Sathyanathan, Errard, Basak 2025 (FT-part): a bias-corrected
"ridge likelihood" (Eqs. 10, 12 to 17) sampled with emcee as if it were a posterior, no marginal
versus profile discussion. 0805.0093, Gratton 2008 (FT): foregrounds marginalized analytically
under a Gaussian prior taken to the flat limit (Sect. III, Eqs. 13 to 15), with a determinant in
Eq. 28 and a saddle-point prefactor (Eq. 51); footnote 6: "this term would be missed in any
approach that attempts to use a standard CMB likelihood formula"; no numerical application. This
is the "Gratton" candidate of the task. 1409.0858, Vansyngel et al. 2016 (FT): mixing matrix
sampled by Metropolis-Hastings from the component-marginal (Sect. 2.4). 2301.11895 and 2408.02063
(both ABS): profile versus marginal in cosmological parameter estimation.

(inference) Maximizing per proposal and accepting on the minimum chi-squared samples the profile
posterior, a well-defined target that is not the marginal of the stated model. Redrawing a from
P(a|beta, d) per proposal and accepting on the chi-squared at the drawn a omits the proposal
density from the Hastings ratio and has no stated target. Neither statement was found in a paper.

Referee must-cite: 0804.2645, 2201.08188, 1608.00551, 0709.1058, 1004.4756 (Stivoli et al. 2010,
BIB), 2405.18365 with 2308.03412 (BIB), 2511.04546, 0805.0093. Statistics (BIB only, DOIs
verified on Crossref, relevance is my inference): Liu 1994 (10.1080/01621459.1994.10476829), van
Dyk and Park 2008 (10.1198/016214508000000409), Andrieu and Roberts 2009 (10.1214/07-AOS574), Cox
and Reid 1987 (10.1111/j.2517-6161.1987.tb01422.x), Berger, Liseo, Wolpert 1999
(10.1214/ss/1009211804), Golub and Pereyra 1973 (10.1137/0710036).

### Claim 2. Prior-volume term at low S/N and coverage versus S/N

Verdict: partly known. The mechanism, its growth with the number of amplitude pixels per index, and
several remedies are in print since 2008 to 2017. Not found: any coverage-versus-S/N curve, any
S/N threshold for calibrated coverage, any rank-histogram or SBC test of index posteriors.

Closest papers:
- 0811.4280, Dunkley et al. 2009, WMAP5 polarization (FT, two passes). Not in LITERATURE.md.
  Gibbs with Metropolis for the index conditional on amplitudes (Sect. II.2.2, Eqs. 19 to 20);
  48 index pixels against 768 amplitude pixels. Sect. II.2.2: "In the limit of low
  signal-to-noise, there is a larger volume of models with a flatter (i.e. beta tending to 0)
  synchrotron spectrum." Sect. III.1: a noise-only simulation gives marginal index distributions
  "centered on beta=-2.6" under a -3.0 +- 0.3 prior; the Jeffreys prior is "a common approach" but
  "this would not overcome the phase space effect"; remedy p(beta) = p_s(beta)/L_0(beta) with L_0
  measured on a noise-only simulation; "the phase-space volume ... increases for more pixels
  sharing a common index parameter"; caveat "the prior is informative in the presence of a signal
  and could therefore lead to parameter bias". Pixels with sigma(beta) > 0.25 masked (Sect. IV.2).
- 1608.00551, Alonso et al. 2017 (FT). Not in LITERATURE.md. App. B, Eq. 34: the volume factor
  sqrt|N_T| "becomes arbitrarily large for beta_s -> -infinity"; the volume prior is "equivalent
  to the widely used Jeffreys prior for the spectral parameters" (Sect. II.3.2) "except for the
  subdominant logarithmic term" (App. B). Fig. 11: at high S/N the factor "biases the estimate of
  the spectral indices by more than 4 sigma" (text; the caption says more than 3 sigma), and at
  low S/N the means sit "hundreds of" sigma from truth even with a broad Gaussian prior. No S/N
  threshold, no coverage test. (inference) This anticipates project finding F4 (ridge equals the
  marginal under a conditional Jeffreys-type prior of the amplitude block).
- 0709.1058, Eriksen et al. 2008 (FT-part, first 100 000 of about 156 000 characters). Jeffreys
  prior, Sect. III.1 and Eq. (36). Sect. IV.2, Figs. 3 and 4: with per-pixel amplitude and index
  plus one global offset m, a uniform index prior gives m = "-1.8 +- 0.4, more than 4 sigma away
  from zero"; "there is a larger allowed parameter volume with steep indices than with shallow";
  "a proper ignorance prior leads to unbiased estimates". No S/N scaling.

Other relevant: 1502.01588, Planck 2015 X (FT-part, 0 to 200 000 of about 249 000), Sect. 2.1:
"we adopt a Jeffreys prior for the spectral parameters in order to suppress prior volume effects".
0903.4311, Dickinson et al. 2009 (FT-part, first 100 000 of 111 149), Sect. III.1: priors
"multiplied by the appropriate Jeffreys prior to account for the different volume of likelihood
space"; Sect. IV.4.1: "the prior is important in regions of weak emission". 1203.0152,
Armitage-Caplan et al. 2012 (FT): flat amplitude and Gaussian index priors (Sect. 4); bias on r up
to about 3 sigma (abstract); "These biases arise from over-parameterizing the model in low
signal-to-noise regions" (Fig. 3 caption); "the result is being driven by the prior" (Sect. 5.1);
no volume language. 1103.2554 (FT), Sect. 4.1.1: "the prior is stronger than the likelihood" away
from the plane. 1907.11426, Jew and Grumitt 2019 (FT), Sect. 4: in the faintest pixels "the means
of the posterior distributions are biased towards shallow spectral indices". 1907.11642, Jew et
al. 2019 (FT), Sect. 3.1: "Flat priors are not always uninformative". 2509.11894 (FT-part),
Sect. 2.5: "For non-linear parameters such as spectral indices, uniform priors are generally
informative". 2606.21334 (FT-part), Sect. 5.4: "in low S/N regions posterior distributions are
partially prior-driven". 2301.05117, QUIJOTE VIII (FT-part), Sect. 5.1.2: a cut where "beta_s is
estimated with a signal-to-noise ratio larger than 15"; this is a cut on beta_s/sigma(beta_s), a
different quantity from the project's index-pixel S/N of about 15. Dust SED analogues: 0902.0636
and 0902.3477 (ABS, noise produces a spurious T-beta anticorrelation in chi-squared fits),
1203.0025 (FT-part, hierarchical remedy), 1305.2130 (ABS, the summarizer's gist: all methods
retain some bias, worst for hierarchical models at high noise).

(inference) The papers disagree on the direction of the push (shallow in 0811.4280 and 1907.11426,
steep in 1608.00551 and 0709.1058). The determinant is that of the whole mixing matrix, so the
direction depends on the reference frequency and on which other components are present. This is
consistent with the project's E030 pass 7 and should be said explicitly in the paper.

Referee must-cite: 0811.4280, 1608.00551, 0709.1058, 1502.01588, 0903.4311, 1203.0152, 1907.11426,
1907.11642, 1203.0025, and Neyman and Scott 1948 (10.2307/1914288, BIB, relevance is inference).

### Claim 3. Hierarchical, template-tied amplitude prior validated by SBC

Verdict: partly known. Not found in CMB component separation: a Gaussian prior on foreground
amplitudes with mean and scale learned from the data, per-region template multipliers, presented
as the cure for the index volume bias, or any simulation-based calibration of a CMB component
separation sampler. Known: every ingredient separately.

Closest papers:
- 2201.08188, BP XIII (FT-part). Sect. 4.4.2, Eqs. 25, 28 to 29: informative Gaussian AME
  amplitude prior with mean equal to the Planck DR4 857 GHz map scaled by alpha (2 to 4 x 10^-5,
  chosen by a scan, Fig. 7) and a power-law width chosen by a chi-squared scan (Fig. 9).
  Sect. 4.4.1: harmonic smoothing priors with widths set by hand. No hyperparameter is sampled.
- 1203.0025, Kelly et al. 2012 (FT-part, first 100 000 of about 133 000). Per-pixel (log N, log T,
  beta) follow a multivariate Student-t whose mean and covariance are inferred (Sect. II.3,
  Eqs. 10 to 11), so the learned prior includes the amplitude; calibration errors sampled jointly;
  "convergence of our MCMC sampler proves to be extremely slow due to the strong dependencies".
- 1907.11426, Jew and Grumitt 2019 (FT). Amplitude prior = reference prior proportional to A1 A2
  (Eq. 6) times a maximum-entropy exponential prior whose mean is the Haslam map extrapolated with
  beta = -3.11 and a 20 percent polarization fraction (Sect. 2.2, Eqs. 7 to 8). A
  template-informed amplitude prior on real Planck data, hyperparameters fixed.
- 2211.10448, Pagano et al. 2024 (ABS). Global 21 cm: a monopole offset and per-region scale
  factors multiply the base map, number of regions chosen by evidence. With 2010.09644 (ABS,
  per-region index scaling a base map) this is the adjacent-field form of template multipliers.

Other relevant: 0811.4280 (FT), Sect. II.2.1: dust Q/U amplitude prior 0 +- 0.2 times the dust
intensity. 0902.4803, Efstathiou, Gratton, Paci 2009 (FT), Sect. 4.2, Eqs. 33 to 35: foregrounds
integrated out under a Gaussian prior with mean F_0 and covariance P, then P is neglected.
0804.2645 (FT), Sect. 3.1.1: Gaussian amplitude priors "can be easily incorporated here".
2204.09360, Platz et al. (FT-part): a template "imprints ... into the prior mean of the component"
(Eq. 9) with a learned modification field, gamma rays, NIFTy. 2501.01417, B-GSM 2025 (FT-part):
Gaussian component prior with fixed hyperparameters (Eqs. 28 to 30, 39 to 40) and stated prior
sensitivity (Sect. 9). 2509.11894 (FT-part), Sect. 2.3: "Sampling the prior covariance S is left
for future work". 1910.14170, Grumitt et al. 2020 (FT-part): amplitudes Unif(-inf, inf) in both
models (Tables 2 and 3), hierarchy on spectral parameters only, validation by normalized
deviations (Eq. 24, Fig. 3), no SBC. 1910.08077, HGMCA (FT-part): blind, MAP, hierarchy on mixing
matrices. 2408.10952, Cosmoglobe DR2 I (FT-part), Sect. 2.1.3: "we adopt existing spatial
templates as delta function priors on the spatial morphology". 1801.06660, HerBIE (ABS).

(inference) The CMB power-spectrum Gibbs sampler (astro-ph/0310080, astro-ph/0407028, both ABS)
is a learned-variance Gaussian prior on the CMB amplitudes. The new element is giving foreground
amplitudes a learned prior and showing that it removes the index volume bias.

Referee must-cite: 2201.08188, 0811.4280, 1907.11426, 1203.0025, 1801.06660, 1910.14170,
2211.10448, 2501.01417, 2204.09360, 1804.06788 (SBC, in LITERATURE.md, not re-read), and
1305.2130 as the caution that hierarchical models can themselves be biased at high noise.

### Claim 4. GMRF or SPDE prior on the index map with inferred hyperparameters, NUTS

Verdict: new for CMB component separation as far as searched. Partly known across fields.

Closest papers:
- 2606.04094, Chavez et al. 2026, CHIBI (FT-part, first 100 000 of about 168 000). VLBI
  multifrequency imaging: zero-mean GMRF prior on the spectral index map (Sect. II.2, Eq. 19),
  correlation length under an inverse-gamma hyperprior (Eq. 18), "we fit a unique value of rho for
  each map"; Adam then NUTS (Sect. II.3); real VLBA data and synthetic EHT data; flat images up to
  256 x 256. Not CMB, not on the sphere, no Gibbs or collapsed structure.
- 1401.4711, Junklewitz, Bell, Ensslin 2015 (FT). Gaussian prior on the index map (Eq. 13); "The
  power spectrum of the spectral index is reconstructed together with the field itself"
  (Sect. 3.3, App. A.3); MAP with alternating conditional solves; VLA simulations only.
- 2011.11575, Azzoni et al. 2021 (FT-part, first 100 000 of about 134 000). Index fluctuations as
  a Gaussian field with C_l^beta = B (l/l_0)^gamma, B and gamma free (Eq. 2.21, Table 1), inferred
  at power-spectrum level with emcee; no index map reconstructed; SO-like simulations and BK15.

Other relevant: 2205.13968 (NOTE): fixed 5 degree Gaussian-kernel prior, MAP, real Planck Q/U.
1910.14170 (FT-part): exchangeable hierarchy, NUTS, simulations. 1907.11426 (FT): global mean and
intrinsic scatter (whole sky -2.99 +- 0.03, scatter 1.12 +- 0.02, Sects. 4 and 5.1), pixels
independent, index-map power spectrum only as a diagnostic (Sect. 5.2). 0902.4803 (FT), Sect. 5:
"nobody has yet developed a model fitting technique that incorporates prior information on the
angular variation of the spectral indices" (a 2009 statement). 2504.07090 (ABS): field-level HMC
with foreground parameters for 21 cm. 2510.08534 and 2604.08463 (NOTE): patch-based indices.

Not found: a GMRF, SPDE, Matern or Gaussian-process prior with inferred hyperparameters on a
spectral-index map inside a CMB component-separation sampler (searches in Sect. 5).

Referee must-cite: 2606.04094, 1401.4711, 2011.11575, 2205.13968, 1910.14170, Lindgren, Rue,
Lindstrom 2011, and (BIB only) 1403.4630, 1503.00256, 0708.3797, 1312.0906.

### Claim 5. Conditional Gibbs mixing at high S/N, collapsed sampling as the cure

Verdict: partly known. The slow mixing and the collapsed cure are both published. Not found: an
ESS or correlation-length comparison of conditional against collapsed index sampling as a function
of S/N in CMB component separation.

Closest papers:
- 2201.08188, BP XIII (FT-part), Sect. 4.1 and Fig. 1: Gibbs has an "inability to probe
  degenerate distributions", with "poor Markov chain mixing and a very long correlation length";
  toy with correlation 0.99: "less than half of the distribution is actually explored", "many
  thousands of samples will be required"; the marginal-then-conditional joint step is the cure.
- 1608.00551 (FT), App. A.2: "By skipping the intermediate sampling of the amplitudes ..., we
  reduce the correlation length of the the MCMC chains for b" (Fig. 10). No numbers.
- 2509.11894, Nasirudin and Bull 2025 (FT-part, first 100 000 characters), Sect. 4: with the index
  sampled conditional on the sky, integrated autocorrelation times of about 112 (sky), 22 (gains)
  and 72 (index), ESS 9, 44 and 14 from 1000 samples, "which is rather low"; analytic
  marginalization named only as future work (Sect. 4.3). Radio, six index regions.

Other relevant: 2011.08503, BP XV (NOTE): correlation length 30 to 50 (Sect. 5.3.2). 0709.1058
(FT-part), Sect. III.2.2: "hundreds of MCMC steps may be required to generate two uncorrelated
samples". astro-ph/0508268 (FT-part), Sect. V.4.2: about 10^7 samples per pixel. 0811.4280 (FT),
Sect. II.4: "about 10,000 iterations are required" with the index step. 1103.2554 (FT),
Sect. 2.2.2: three index iterations per loop "to reduce correlations between consecutive samples".
1502.01588 (FT-part), Sect. 2.2: "this algorithm is notorious for its slow convergence for nearly
degenerate parameters". 1409.0858 (FT), Sect. 2.3.1: full-posterior schemes "inefficient because
of high correlation length in the chains". 2303.08095, Cosmoglobe DR1 (FT-part, characters
100 000 to 300 000 of 314 116), Sect. 6.1: strong SED priors motivated by foreground and gain
degeneracies "which otherwise can lead to very long Markov chain correlation lengths"; this is
about gains, not amplitude and index. 1203.0025 (FT-part): interweaving needed.

Referee must-cite: 2201.08188, 2011.08503, 1608.00551, 2509.11894, 0709.1058, 1409.0858, and
(BIB only) Roberts and Sahu 1997 (10.1111/1467-9868.00070), Liu 1994, van Dyk and Park 2008.

### Claim 6. Zero-level convention and large-region indices

Verdict: partly known. The degeneracy, its treatment as a nuisance, and index errors of 0.03 to
0.07 from monopoles are published. Not found: a statement that the zero-level convention moves
full-sky or large-region indices by many posterior standard deviations (project F17).

Closest papers:
- 0709.1058 (FT-part), Sect. IV.1: "This degeneracy is a crucial issue for CMB component
  separation"; Sect. IV.4, Eqs. 41 to 43: "It is thus impossible to determine individually the
  absolute zero-level of foreground component and offsets"; "Spatial variations in the spectral
  index partially resolve this degeneracy". The more-than-4-sigma example of Sect. IV.2 is a bias
  of the offset m driven by the index prior.
- 1411.7616, Wehus et al. 2017 (FT). T-T regression with one index per region (Sect. 2.2);
  "While spurious monopoles do not change the net slope of a T-T plot, spurious dipoles do"
  (Sect. 2.3); dipole-induced index bias typically +-0.2 for 408 against 1420 MHz (Fig. 9,
  Sect. 4); "The offsets are nearly degenerate with the effective zero-level of each foreground
  component" (Sects. 1 and 5); Haslam 8.9 +- 1.3 K (Table 1).
- 2606.21334, Hoerning et al. 2026 (FT-part, 0 to 200 000 of about 205 000). Monopoles set by
  anchoring and T-T fits (Sect. 3.2, App. B), not shown to be sampled jointly; Sect. 4.1.2:
  sigma_beta from a 1 mK monopole uncertainty at 4.76 GHz is about 0.03 in bright regions and up
  to about 0.06 at high latitude, combined floor 0.05 to 0.07; Sect. 5.4: "the derived synchrotron
  spectral index is sensitive to the adopted relative offsets between frequency channels".

Other relevant: 2201.08188 (FT-part), Sect. 4.3: "for each proposed spectral parameter value, we
estimate a new band monopole value"; "by marginalizing over the band monopoles the log-likelihood
function widens by a factor of 1.5-2". 0804.2645 (FT), Sect. 4: "arbitrary global offsets" as
marginalized nuisance parameters (Eqs. 14, 19 to 22). 1502.01588 (FT-part), Sect. 4.1, Table 6:
several monopoles fixed; the analysis "relies critically on external priors for these values".
0903.4311 (FT-part), Sect. III.4: without monopole and dipole fitting "there is a large-scale
gradient running across the map, which is unlikely to be real"; monopoles shift by up to about
10 uK across the index prior range (Sect. IV.4.1). 2501.01417 (FT-part), Sect. 8.2: the index
"shows a strong positive correlation with the scale correction ... and a strong negative
correlation with the zero correction". 2303.08095 (FT-part): "the monopoles are in general
determined by the foreground model" (Sect. 5.4 or 6.2, the two passes disagree on the section).
2408.10952 (FT-part, first 100 000 of 148 511): one monopole per band sampled (Eq. 1); about 10
percent correlations with zodiacal parameters (Sect. 5). 1907.11426 (FT), Sect. 5.3.5: removing
8.9 K from the prior map leaves mean indices "hardly affected".

Referee must-cite: 0709.1058, 1411.7616, 0804.2645, 2201.08188, 1502.01588, 0903.4311,
2606.21334, 2501.01417, 2303.08095, 2408.10952, 2012.11019 (BIB).

## 2. Table

| arXiv ID | Read | Finding | Relevance |
|---|---|---|---|
| 0804.2645 | FT | Profile likelihood Eq. 9 "not a proper likelihood"; marginal Eq. 10; offsets Sect. 4 | 1, 2, 6 |
| 2201.08188 | FT-part | Marginal Metropolis then one amplitude draw; ridge option; AME template prior; monopoles | 1, 3, 5, 6 |
| 1608.00551 | FT | Marginal index sampling; volume prior; biases of 4 to hundreds of sigma | 1, 2, 5 |
| 0811.4280 | FT | Phase-space effect, noise-only centre -2.6, correction prior p_s/L_0 | 2, 3, 5 |
| 0709.1058 | FT-part | Jeffreys prior; offset biased by more than 4 sigma; offset degeneracy | 2, 5, 6 |
| 1502.01588 | FT-part | Jeffreys prior "to suppress prior volume effects"; fixed monopoles | 2, 5, 6 |
| 0903.4311 | FT-part | Jeffreys factor; prior-driven high latitude; monopole and dipole fit | 2, 6 |
| 1203.0152 | FT | r biased up to 3 sigma by over-parameterized low-S/N index pixels | 2 |
| 1103.2554 | FT | Flat amplitude priors, prior dominates off-plane, 10 000 samples | 2, 5 |
| 1907.11426 | FT | Reference times Haslam-based amplitude prior; shallow bias in faint pixels | 2, 3, 4, 6 |
| 1907.11642 | FT | Jeffreys-rule priors; "Flat priors are not always uninformative" | 2 |
| 2301.05117 | FT-part | S/N > 15 cut defined on beta_s/sigma(beta_s) | 2 (different quantity) |
| 0805.0093 | FT | Analytic foreground marginalization, determinant prefactor, no application | 1, 3 |
| 0902.4803 | FT | Gaussian foreground prior integrated out then neglected; 2009 gap statement | 3, 4 |
| 1203.0025 | FT-part | Learned Student-t prior over (log N, log T, beta); slow mixing | 2, 3, 5 |
| 1910.14170 | FT-part | Flat amplitudes, hierarchy on indices only, NUTS, no SBC | 3, 4 |
| 1910.08077 | FT-part | Blind MAP, hierarchy on mixing matrices | 3 |
| 2501.01417 | FT-part | Collapsed Gaussian-prior likelihood, zero level and scale sampled | 3, 5, 6 |
| 2509.11894 | FT-part | Conditional index step, autocorrelation 72 to 112, prior covariance fixed | 2, 3, 5 |
| 2204.09360 | FT-part | Template in the prior mean with learned modification field | 3, 4 |
| 2211.10448, 2010.09644 | ABS | Per-region scale factors and indices on a base map, evidence selection | 3, 4 |
| 2606.04094 | FT-part | GMRF on index map, correlation length inferred, NUTS, VLBI | 4 |
| 1401.4711 | FT | Index-field power spectrum inferred with the field, MAP | 4 |
| 2011.11575 | FT-part | Amplitude and tilt of index-field spectrum inferred in C_l space | 4 |
| 1409.0858 | FT | Collapsed sampling of the mixing matrix | 1, 5 |
| 1411.7616 | FT | Offsets from T-T regression; dipoles bias indices by about 0.2 | 6 |
| 2606.21334 | FT-part | Commander low-frequency model; index floor 0.05 to 0.07 from offsets | 2, 5, 6 |
| 2303.08095 | FT-part | SED parameters prior-driven; monopoles set by the foreground model | 5, 6 |
| 2408.10952 | FT-part | Monopole per band sampled; global dust indices conditional on amplitudes | 5, 6 |
| 2405.18365 | FT-part | Profile likelihood with a separate bias correction, amplitudes maximized | 1 |
| 2511.04546 | FT-part | Bias-corrected ridge likelihood sampled with emcee | 1 |
| 2609.04510 | FT-part | Commander2 with correlated noise, amplitudes only | none (checked) |
| 2507.05324 | FT | LiteBIRD with Commander3, indices held at input values (Sect. 3.1) | none (checked) |
| 2604.16901 | FT | Cosmoglobe overview; Commander4 named, no methods paper cited | none (checked) |
| 0811.3915, astro-ph/0508268 | FT-part | Typical index priors; per-pixel MCMC with 10^7 samples | 2, 5 |
| 0902.0636, 0902.3477, 1305.2130, 1801.06660 | ABS | Noise-induced T-beta bias and hierarchical remedies | 2, 3 |
| 2301.11895, 2408.02063 | ABS | Profile versus marginal and volume effects in cosmology | 1, 2 |

## 3. 2025 and 2026 papers that overlap with the proposed methods paper

None overlaps substantially. Partial overlaps that must be cited:
- 2606.21334 (2026): the same code family on real data; quantifies claim 6 at the 0.03 to 0.07
  level; no sampler analysis, no convergence diagnostics, Commander4 not mentioned.
- 2509.11894 (2025): measured slow mixing of a conditional index step (claim 5) and a remark on
  uniform index priors (claim 2), radio.
- 2501.01417 (2025): collapsed likelihood with zero-level and scale sampling, radio, one global
  spectrum per component.
- 2606.04094 (2026): claim 4 in VLBI imaging.
- 2511.04546 and 2502.15511 (2025, the second ABS): corrected ridge and minimally informed
  likelihoods; a different bias from the one in claim 2.
- 2510.08534 and 2604.08463 (already in LITERATURE.md): patch-based indices.
- 2604.26890 and 2504.07090 (ABS): Gibbs and HMC foreground models for 21 cm.

Checked and not overlapping: 2609.04510, 2507.05324, 2604.16901, 2512.07235, 2512.08522,
2607.00299, 2603.02177, 2502.04300. `all:Commander4` returns 0 arXiv results.

## 4. Contradictions and nuances against the repo notes

1. LITERATURE.md line 35 ("Jeffreys prior removes >4 sigma prior-volume bias"): in 0709.1058
   Sect. IV.2 the more-than-4-sigma bias is in the global offset m, not in the index.
2. LITERATURE.md lacks 0811.4280 and 1608.00551, the closest prior art for claim 2. The second
   anticipates F4 (inference). It also lacks 1907.11426, 1203.0152, 1203.0025, 0805.0093,
   2606.04094, 1401.4711, 2011.11575, 2501.01417, 2509.11894, 1411.7616, 1502.01588, 0903.4311.
3. `2026-10-27_index_field_recovery_prior_art.md` says no paper infers a correlation length or
   spectrum for an index map inside a CMB component-separation sampler. That still stands for
   pixel-domain CMB samplers, but 2606.04094, 1401.4711 and 2011.11575 should be added.
4. `2026-10-27_beyondplanck_chain_lengths_and_index_sampling.md` says no paper reports a
   numerical amplitude-index correlation length. 2509.11894 now does, outside CMB.
5. `04_adjacent_field_methods.md` lists 2509.11894 without its mixing numbers.
6. LITERATURE.md line 22: 2507.05324 keeps the indices "at their true input values" (Sect. 3.1)
   although its Table 2 lists a spectral-index sampling cost. Not a contradiction.
7. `01_spatially_varying_parameters.md`: no contradiction found.

## 5. Searches run and what was not found

arXiv API (`export.arxiv.org/api/query`), result counts in brackets:
- `abs:"component separation"`, 2025-01-01 to 2026-11-09 [165, all pages scanned].
- Same with `AND abs:Bayesian` [7], `AND abs:parametric` [21], `AND abs:"spectral index"` [6].
- `abs:foreground AND abs:"spectral index" AND abs:CMB`, 2025 to 2026 [6].
- `abs:foreground AND abs:Gibbs`, 2025 to 2026 [2]. `abs:"spectral index" AND abs:Gibbs` [1].
- `abs:foreground AND abs:hierarchical AND abs:CMB`, since 2024-06 [3, none relevant].
- `abs:foreground AND abs:"Hamiltonian Monte Carlo"`, since 2024-06 [4, one adjacent].
- `abs:"spectral index" AND abs:Bayesian AND abs:synchrotron`, 2025 to 2026 [6].
- `abs:"Markov random field" AND abs:"spectral index"` [0]; `abs:"spectral index" AND
  abs:"spatial prior"` [0]; `abs:dust AND abs:"spectral index" AND abs:"Gaussian process"` [1,
  irrelevant]; Gaussian process, Matern or correlated field with spectral index [3, irrelevant].
- `abs:"spectral index" AND abs:hierarchical AND abs:foreground` [0, although 1910.14170 exists];
  `abs:"spectral index" AND abs:"hierarchical Bayesian" AND abs:synchrotron` [1, HerBIE use].
- Prior volume with foreground [4, none relevant]; profile likelihood with component separation
  and CMB [2, none relevant]; zero level, monopole or offset with spectral index [9];
  `abs:Cosmoglobe` [16]; `all:Commander4` [0]; `au:Gratton_S AND abs:foreground` [8].
- `id_list` verification of every ID cited as BIB or ABS.

INSPIRE citation lists, scanned by title: 0804.2645, 0709.1058, 1910.14170, 2201.08188,
2205.13968, 1411.7616 (counts not retained), 1608.00551 since 2022 [11], 0811.4280 since 2021
[24]. WebSearch (6 queries): hierarchical Bayesian component separation 2026; Gaussian process or
Markov random field prior on a spectral index map; prior volume and spectral index bias;
"Commander4"; simulation-based calibration with CMB component separation; Commander3 parameter
documentation for the ridge, marginal and chisq options (nothing found). Crossref for the eight
statistics DOIs.

Not found in any of the above: (a) an analysis of the validity or target of a chi-squared index
step with re-solved or redrawn amplitudes; (b) coverage of index posteriors as a function of S/N,
or an S/N threshold; (c) a learned hierarchical prior on foreground amplitudes in CMB component
separation, or SBC of such a sampler; (d) a GMRF, SPDE or Gaussian-process index-map prior with
inferred hyperparameters in CMB component separation; (e) conditional against collapsed ESS as a
function of S/N; (f) index shifts from the zero-level convention expressed in posterior sd; (g) a
Commander4 methods paper.

## 6. Gaps and suggested follow-up

Gaps:
- Every FT-part entry has an unread remainder (ranges given above). BP XIII beyond 100 000
  characters was read only in the 2026-10-27 note.
- Statistics references were verified bibliographically only. Their relevance is my inference.
- Not read beyond title or abstract: 1004.4756, 2308.03412, 2502.15511, 1801.06660, 2308.09176,
  2211.10448, 2010.09644, 2504.07090, 2604.26890, 1305.2130, the maximum-entropy separation
  papers astro-ph/9806387 and astro-ph/0405494, 1404.5323, 1909.05923, 1802.01145, 1509.04714,
  1709.07897, 2012.11019.
- Abstract-term searches miss papers that use other words (shown by the Grumitt null result).
  INSPIRE lists were scanned by title only. Journals outside arXiv were not searched.
- Summarizer inconsistencies left unresolved: 1608.00551 Fig. 11 text (4 sigma, 1 sigma) against
  caption (3 sigma, 2 sigma); 0811.4280 "92" against "96" index pixels; 2606.21334 Table 5 C-BASS
  monopole; 2303.08095 section number of the monopole sentence; 1203.0152 Test G bias (0.6 against
  1 sigma).

Suggested follow-up:
1. Check in the PDFs the eight quotes the paper will lean on most: 0804.2645 Sect. 3.1.1,
   2201.08188 Sect. 4.1 and 4.3, 0811.4280 Sect. III.1, 1608.00551 App. B and Fig. 11,
   0709.1058 Sect. IV.2 and IV.4, 2606.21334 Sect. 4.1.2.
2. INSPIRE citers of 0811.4280 and 1608.00551 before 2021, for any coverage study.
3. Commander3 source and documentation for the `chisq`, `ridge` and `marginal` options, to state
   exactly what the released BeyondPlanck and Cosmoglobe runs used.
4. Fuskeland et al. (1404.5323, 1909.05923) and Krachmalnicoff et al. (1802.01145) full text, for
   low-S/N bias of T-T index estimators.
5. Second half of 2606.04094 and the Comrade.jl documentation for the GMRF hyperprior details.
6. Read 1305.2130 in full: its bias of hierarchical models at high noise is the obvious referee
   objection to claim 3.
