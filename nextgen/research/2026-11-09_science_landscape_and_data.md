# Science landscape and public data: polarized synchrotron index and tau

Date: 2026-11-09. Author: lit-scout agent for Trygve Leithe Svalheim.

Question: where could a statistically rigorous Bayesian analysis of spatially varying foreground
spectral parameters, with calibrated uncertainties, make the largest scientific difference using
public data, in (a) the polarized synchrotron index beta_s and (b) the optical depth tau?

How to read this file:
- Read levels: FT = full text; FT-part = stated part of the full text; ABS = abstract only;
  PAGE = landing or documentation page; LOCAL = file on this machine.
- Every number was extracted through a small summarizer model that returns short quotes. Section,
  table and equation references are given so each can be checked. Check the PDF before citing.
  Items where two extractions disagree are marked UNRESOLVED.
- "(inference)" marks my own synthesis or arithmetic, never a result of the cited paper.
- No data files were downloaded. Data facts come from landing pages and files already on disk.

## 1. Short answer

State of the field. (a) Away from the bright structures every published map of polarized beta_s is
prior dominated or noise dominated, and the authors say so (2301.05117, 2201.03530, 2310.13740,
2511.14572). High-latitude values range from about -3.0 to -3.3 with quoted errors of 0.01 to 0.10
that mostly omit calibration. A new QUIJOTE power-spectrum result (2610.31469, 2026-10-22) gives
-3.26 +- 0.05 at |b| > 20 deg, 0.15 to 0.25 steeper than the same collaboration's map-domain
values. No analysis infers the smoothness of the index field and none validates coverage of an
index map. (b) Low-multipole tau falls in two groups: HFI based, 0.051 to 0.060 with errors near
0.006, and WMAP plus LFI, 0.062 to 0.069 with errors 0.011 to 0.013. DESI DR2 BAO plus CMB without
low-l polarization gives 0.090 +- 0.012 (2504.16932). All of them treat synchrotron with one
template coefficient per band or with fixed or prior-drawn indices. The statement that the beta_s
prior moves tau exists as one qualitative sentence in BeyondPlanck XV Sect. 4.4, with no number.
Cosmoglobe DR1 has no tau.

Ranked candidate analyses (all inference; ranking is value times workstation feasibility with
maps-only inputs):

1. **tau from WMAP plus LFI with the synchrotron index field marginalized in the pixel domain.**
   Joint nside 16 Q/U model of WMAP K to V (res 4 maps and N^-1, including loss-imbalance modes),
   LFI 30/44/70 and 353 GHz (PR4 nside 16 maps and dense covariances), with beta_s as constant,
   regional, and smooth field with inferred smoothness, and tau sampled with it. Deliverable: tau
   as a function of the index model and prior, with simulation-based coverage. Why first: it
   quantifies what BP XV only asserts, sits exactly between 0.058 (HFI), 0.066 to 0.069 (WMAP plus
   LFI) and 0.09 (DESI implied), and no published analysis does it (Sect. 3.4). Feasible:
   matrices are 144 MB to 648 MiB each; marginalizing three components per pixel is one dense
   Cholesky of order 18 432, about 2e12 flops (operation count, not a benchmark). Risks: NPIPE
   low-l transfer function at 44 and 70 GHz; LFI covariance "is not capturing the full error
   budget" (1907.12875 Sect. 2.3.5); the 600 PR4 simulations need a NERSC account; band-specific
   smoothing differs between WMAP res 4 and PR4 low-resolution maps.
2. **K band versus 30 GHz large-scale consistency as an explicit low-dimensional systematic
   model.** Prerequisite for 1 and 3 and the cheapest. Evidence that the two synchrotron tracers
   disagree: Fan region index -2.90 +- 0.020 from [23, 33] GHz against -3.27 +- 0.022 from
   [23, 30] GHz (2203.11445 Table 1); LFI 70 GHz tau of 0.094 +- 0.016 cleaned with K against
   0.053 +- 0.023 or 0.076 +- 0.020 cleaned with LFI 30 (1801.01226); template coefficient ratio
   30 GHz over K of 1.54 +- 0.13 in Planck 2018 V against 1.96 to 1.98 expected for beta = -3.1
   (Sect. 3.3). Compare WMAP9 plus PR3/PR4 with the Cosmoglobe DR1 reprocessed maps.
3. **Northern beta_s map from QUIJOTE 11 and 13 GHz, WMAP K and Ka, and Planck 30 GHz** with
   calibration, polarization angle and the FDEC filter inside the model, a spatial prior with
   inferred smoothness, and simulation-based coverage. Resolves map domain (-3.0 to -3.1) against
   power spectrum (-3.26 +- 0.05) and the documented prior dominance. A 5 percent calibration
   error is ln(1.05)/ln(22.8/11.1) = 0.068 in beta_s, several times the 0.01 to 0.02 regional
   errors quoted in 2510.17761. All maps are public; noise simulations are on request only; no
   pixel covariance is released; 11 and 13 GHz share 1/f noise.
4. **Southern sky with S-PASS 2.3 GHz, Faraday rotation and depolarization modelled jointly, plus
   CLASS 40 GHz.** Gives the long lever arm and a curvature test, but risk is high: ground
   subtraction in declination rings can remove sky signal (1903.09420 Sect. 5), rotation measures
   are underestimated where depolarized, and CLASS has only a C_l-domain transfer matrix.
5. **Curvature or running in bright regions.** Lowest value: the papers themselves cannot
   distinguish models with present data (2301.05117 Table 2; 2610.31469 App. C).

## 2. Part A: polarized synchrotron index

### 2.1 Evidence table

| arXiv, read level | Data, method, resolution | Values and where from | Uncertainty derivation and caveats |
|---|---|---|---|
| 0811.4280v2 Dunkley 2009, FT | WMAP5 K to V, nside 8, full N^-1; Gibbs; index in 48 nside-2 pixels, prior -3.0 +- 0.3 | mean -3.03 +- 0.04; low latitude -3.00 +- 0.04, high -3.08 +- 0.06 | posterior; see Sect. 3.3 for tau |
| 1807.06208 Planck 2018 IV, ABS + sections | Planck 30 to 353; Commander, SMICA | polarization fit has "only a spatially-constant spectral index for synchrotron", fixed -3.1; SMICA -3.10 +- 0.06 | "the Planck measurements by themselves have little sensitivity to the synchrotron spectral index" |
| 1801.04945 Planck 2018 XI, FT-part (first 100k of 227k chars) | Planck 30, 100 to 353 and WMAP 23, 33 cross-spectra; six regions LR24 to LR71; l 4 to 159 | beta_s = -3.13 +- 0.13 (Sect. 4.2), weighted mean and sd over regions and bins, then used as a Gaussian prior; dust-synchrotron correlation rho fitted, detected only at l <~ 40 | spread over fits; per-region values in App. C not read |
| 1802.01145v3 Krachmalnicoff 2018, FT-part (first 100k of 121k) | S-PASS 2.3 GHz, WMAP, Planck; harmonic SED fit l 20 to 140, masks abs(b) > 20 to 50 deg; pixel map at 2 deg | -3.22 +- 0.08 (Sect. 4.2.1, Table 2); map mean -3.25, dispersion 0.15, "the spread of βs we find is a real feature of the sky emission and not only due to noise" (Sect. 5.2); index-map spectrum slope -2.6 +- 0.2 (Sect. 5.3); s_run "compatible with zero at 1.5σ at most" (Table 3) | 100 simulations plus 5 percent calibration |
| 1909.05923v2 Fuskeland 2021, FT | Faraday-corrected S-PASS vs WMAP K; T-T in 95 regions of about 15 deg; 1 deg, nside 64 | mean -3.24 +- 0.01 corrected (std 0.09, 65 regions), -3.32 +- 0.02 uncorrected (Table 1); plane about -2.8, high latitude about -3.3; flattening 0.1 +- 0.2 from 2.3 to 23 against 23 to 33 GHz (Sect. 4.4) | fit, bootstrap (10 000), rotation-angle half range (Eq. 9); no calibration term; no simulations; "The flat spectral index at low Galactic latitude are likely partly due to depolarization effects" |
| 2110.12803v2 Martire 2022, FT-part (first 100k) | WMAP K x PR3 30; NaMaster, l 30 to 300 | abstract: beta_EE -3.00 +- 0.10, beta_BB -3.05 +- 0.36, alpha_EE -2.95 +- 0.04, B/E 0.22 +- 0.02; Table 4: -3.22 +- 0.09 (94 percent of sky) to -2.98 +- 0.10 (brightest 30 percent) | derivation not recorded |
| 2203.11445v2 Weiland 2022, FT-part (first 100k of 121k) | WMAP K, Ka, LFI 30; zero-intercept Q/U correlation in nside 8 superpixels from 2 deg, nside 32 maps (Sect. III); per-pixel fit with free RM at nside 16 when adding 1.4 or 2.3 GHz (Sect. IV.2) | steepening off the plane 0.4 (south), 0.25 (north) (abstract); BICEP2/Keck field -3.25 +- 0.04 (stat) +- 0.02 (sys) (Sect. VII); Fan region [23,33] -2.90 +- 0.020 vs [23,30] -3.27 +- 0.022; NPS -3.11 +- 0.036 (Table 1) | "The discrepancy traces to large-scale systematic differences between the WMAP 23 GHz and LFI 30 GHz maps"; loss-imbalance bias up to 0.05 not in the errors (Sect. IV.3); "Current public data are insufficient to characterize βpol for the future" |
| 2011.08503v2 BeyondPlanck XV, FT | LFI 30/44/70, WMAP Ka/Q/V, 353 (no K); Gibbs, 4 regions, prior N(-3.3, 0.1) | Spur -3.17 +- 0.06, Plane -3.03 +- 0.07, difference "only statistically significant at the 2σ level" (Sect. 5.3.2); high latitude and Galactic centre are prior draws; nine-region test with prior -2.8 +- 0.1 gives Galactic centre -4.15 +- 0.05 (Sect. 5.3.1) | posterior sd over about 500 samples, 2 chains; high latitude "particularly sensitive to relative gain uncertainties and degeneracies with the CMB quadrupole" |
| 2201.03530v2 BeyondPlanck XVI, FT (see `2026-10-28_beyondplanck_xvi.md`) | per-pixel beta_s, nside 64, prior -3.1 +- 0.1 | no numerical value; prior dominated | not stated |
| 2301.05113 QUIJOTE IV, ABS + Sects. 1 to 2.6, 5 to 11 | MFI 11 GHz vs WMAP K; per-pixel likelihood, prior -3.1 +- 0.3 (Sect. 8.1.2) | "peaks at β=−3.09 with a standard deviation of 0.14"; SED beta_EE -2.99 +- 0.14, beta_BB -2.79 +- 0.35; B/E 0.36 +- 0.04 (abs(b) > 5), 0.26 +- 0.07 (abs(b) > 10) (Table 22) | 1000 noise simulations; calibration 5 percent (11, 13 GHz), 6 percent (17, 19); "the FDEC correction is affecting low multipoles (ℓ≲15)" (Sect. 2.5) |
| 2301.05117 QUIJOTE VIII, FT | MFI 11/13, WMAP K/Ka, Planck; per-pixel B-SeCRET + emcee; 2 deg, nside 64, f_sky 51 percent; priors beta_s N(-3.1, 0.3), c_s N(0, 0.1) | mean -3.08, dispersion 0.13 (Sect. 5.1.3); NPS -3.028 +- 0.002 (Table 1); c_s = -0.0797 +- 0.0012 in RC1 (Table 2), models not distinguishable | "the prior governs the β_s estimation" outside bright regions; "we have considered statistical uncertainties but not calibration errors"; prior N(-3.0, 0.3) moves faint regions by about 0.1 (Sect. 5.5.2); Q-only and U-only maps "distinctly different" (App. A) |
| 2309.00675 Eimer 2024, FT | CLASS 40 GHz with WMAP, Planck; PolSpice l 5 to 125 | s9 mask: -2.95 +- 0.09 (EE), -3.28 +- 0.24 (BB); B/E 0.33 to 0.42; map median -3.09 (+0.23/-0.11) | transfer function keeps 85, 67, 34, 6 percent of power at l = 30, 20, 10, 5 |
| 2404.17567v2 Shi 2024, FT-part (first 100k of 115k) | WMAP K, Ka and CLASS + WMAP Q combined; total least squares at nside 8 | both pairs mean -3.08; uniform-index PTE 0.01 (K, Ka) and < 0.001 (K, combined 40 GHz); Pearson r = 0.04 between the two index maps (Fig. 10) | Fisher errors inflated by simulations (App. C); 5 percent calibration "equivalent to shifting the β_s measurement in all N_side=8 pixels by a common ±0.1" |
| 2310.13740 Cosmoglobe DR1 III, FT | DR1 K, Ka, 30 GHz maps; T-T in 24 regions, 1 deg, nside 64 (Sect. 4.1.1) | K/30: all sky -3.07 +- 0.07, high latitude -3.31 +- 0.07; K/Ka: -2.95 +- 0.07, -3.20 +- 0.10 (Sect. 4.1.2); B/E 0.40 +- 0.02 | rotation-angle systematic and sample sd; the Gibbs chain "adopts an informative βs∼𝒩(−3.15,0.05) prior sampling"; with a N(-3.1, 0.1) prior "the posterior uncertainty is limited by the prior" (Sect. 4.2) |
| 2510.17761v2 QUIJOTE XIX, FT | MFI 11, 13, WMAP K/Ka, PR3 30; moment expansion + cPILC + T-T; 2 deg, beta at nside 32 | weighted mean -3.11, scatter 0.21; plane -3.05 +- 0.01; high latitude (NPS dominated) -3.13 +- 0.02; Table 1: plane -2.96 +- 0.01, high-latitude regions 1 to 14 -3.02 +- 0.02, range -2.67 to -3.40 | effective variance with 100 noise realisations plus rotation-angle term; no calibration term; PySM check: input sd 0.06, recovered sd 0.16, 81 percent correlation; no coverage test |
| 2610.31469v1 QUIJOTE XII, ABS + FT (three windows) | MFI 11.1, 12.9, WMAP K/Ka, PR4; PolSpice l 30 to 129; abs(b) > 10, 15, 20 deg; Stan | "β_s^{EE} = -3.26 ± 0.05", "β_s^{BB} = -3.25 ± 0.10" at abs(b) > 20 deg (abstract); WMAP + Planck alone -3.00 +- 0.09, -3.05 +- 0.11, -3.12 +- 0.12 (Table 4); alpha_EE -3.24 +- 0.12, alpha_BB -3.81 +- 0.34; "abs(s_run)<0.21 at the 3σ level" (App. C); "We find no evidence for spatial variations of the spectral index βs in any of the three Galactic cuts" | Sect. 4.2: errors may be underestimated (single noise realizations); calibration treatment not found. UNRESOLVED: Table 4 extraction gave +- 0.03 for beta_EE |
| 2510.12860 Casas, ABS + first 100k | QUIJOTE MFI only; neural network; angular slope, not beta_s | alpha_EE -3.04 +- 0.18, alpha_BB -3.00 +- 0.26; B/E 0.28 +- 0.06 | claims spatial variation while plane and high latitude agree within errors |
| 2205.13968v1 de Belsunce 2022, FT | Planck 30, 70, 100 to 353 (no K); per-pixel at nside 16, f_sky 0.52; priors -3.0 +- 2.0 with a fixed 5 deg Gaussian kernel; MAP + Laplace | mean beta_s -2.833 +- 0.620; "previous results are strongly prior dominated" | no simulation validation |
| 2511.14572 MFI2 forecast, ABS + first 100k | simulated WMAP + Planck + MFI + MFI2 | plane sigma(beta_s) 0.23, 0.10, 0.06; cosmological fields 0.29, 0.21, 0.12; "improvement factors reach ∼10" | "Most of the sky remains prior-dominated"; normalized bias -0.67 to -0.77 with WMAP + Planck (Table 2); no coverage |
| 2605.23489v1 Shaw 2026, FT | C-BASS north 4.76 GHz ("uses pre-publication data"), S-PASS scaled with -3.1 | no beta_s measured; T-T slope C-BASS vs scaled S-PASS 1.13/1.14 at b > 0, 0.31/0.29 at b <= 0 (App. B) | depolarization <~ 1 to 2 percent at 4.76 GHz against 5 to 16 percent at 2.3 GHz |
| 2403.18123 Delabrouille 2024, ABS | WMAP K/Ka/Q + LFI 30/44 synchrotron Q/U at 30 GHz, 1 deg | assumes index -3.1 | not a measurement |

### 2.2 Where they disagree, and to what it is attributed

- **High-latitude steepening.** Steep: -3.31 +- 0.07 (2310.13740, K/30), -3.26 +- 0.05
  (2610.31469), -3.25 +- 0.04 (2203.11445), about -3.3 (1909.05923), -3.22 +- 0.08 (1802.01145).
  Shallow: -3.02 +- 0.02 (2510.17761 Table 1), -3.08 +- 0.06 (0811.4280), means -3.08 to -3.09
  (2301.05117, 2301.05113). The shallow group is map domain with priors centred near -3.0 to -3.1
  and stated prior dominance. (inference) The split follows estimator and prior more than data.
  Note 2610.31469 itself: WMAP plus Planck alone give -3.12 +- 0.12 and adding QUIJOTE steepens.
- **North Polar Spur.** -3.028 +- 0.002 (2301.05117), -3.11 +- 0.036 (2203.11445), -3.17 +- 0.06
  (2011.08503), -3.13 +- 0.02 (2510.17761). (inference) The spread, 0.14, is up to 70 times the
  smallest quoted error; region definitions were not checked to coincide.
- **K versus 30 GHz.** 2203.11445 attributes the Fan-region 0.37 difference to large-scale
  systematics between WMAP 23 and LFI 30 GHz. 2310.13740 finds K/30 steeper than K/Ka by about
  0.11 on the reprocessed maps.
- **Plane versus high latitude.** Flat plane at 2.3 GHz attributed partly to depolarization
  (1909.05923); 2310.13740: "S-PASS is much more strongly affected by Faraday rotation".
- **Curvature.** Not established: 2301.05117 cannot distinguish models; 2610.31469 and
  1802.01145 find running consistent with zero; 1909.05923 finds 0.1 +- 0.2.
- **Q versus U.** Maps differ (2301.05117 App. A); T-T analyses fold the difference into a
  rotation-angle systematic (1909.05923 Eq. 9, 2310.13740, 2510.17761).
- **E versus B.** B/E ranges from about 0.2 (2610.31469, 2110.12803) to 0.40 +- 0.02
  (2310.13740).
- **Calibration** is absent from the errors of 2301.05117, 2510.17761 and 1909.05923.
- **Internal Cosmoglobe inconsistency.** DR1 I (2303.08095, FT, Sect. 2.4) says the BeyondPlanck
  prior N(-3.30, 0.1^2) is applied; DR1 III and the local parameter file
  `/Users/trygvels/data/cosmoglobe_DR1/param_CG_v09.txt` give N(-3.15, 0.05) with
  `COMP_BETA_POL_LNLTYPE02 = prior` as the last assignment. (inference) DR1 polarized beta_s is a
  prior draw with no data likelihood.

### 2.3 Spatial priors and coverage validation

No analysis uses a spatial prior with inferred smoothness. Closest: fixed 5 deg kernel
(2205.13968), 48 coarse pixels (0811.4280), four regions (2011.08503). 1802.01145 measures the
index-map power spectrum but does not use it as a prior. No analysis validates coverage of an
index map: 2510.17761 reports a correlation with a PySM input, 2511.14572 a normalized bias,
1103.2554 (FT) four simulations with "A thorough test for biases should combine results from
thousands of simulations." This agrees with
`/Users/trygvels/code/Commander4/nextgen/research/2026-11-09_novelty_check_methods.md` claim 4.

## 3. Part B: optical depth

### 3.1 Evidence table

| arXiv, read level | Data and likelihood | tau and where from | Foreground treatment |
|---|---|---|---|
| 1907.12875 Planck 2018 V, FT (windows 100k to 300k) | HFI 100 x 143 SimAll; LFI 70 pixel | lowE + lowl "τ=0.0506±0.0086"; LFI 70 "τ=0.063±0.020"; WMAP with 353 dust template "τ=0.062±0.012" (Sect. 2.3.6) | scalar templates; 100 GHz alpha = 0.0183 +- 0.0012 (30 GHz tracer, Table 2), 0.0119 +- 0.0006 (K tracer, Sect. 2.2.5): "the posterior on τ is not substantially affected by these variations" |
| 1801.01226 Weiland 2018, FT | pixel Q/U likelihood, nside 8 | WMAP KaQV cleaned with K: 0.045 +- 0.019 (46 percent sky), 0.058 +- 0.017 (75 percent); LFI 70 cleaned with LFI 30: 0.053 +- 0.023 or 0.076 +- 0.020; cleaned with K: 0.094 +- 0.016; WMAP cleaned with LFI 30: 0.37 +- 0.02 | attributed to a common-mode LFI systematic, "not consistent with spectral index variations over the sky" |
| 1908.09856 Pagano 2020, FT | SRoll2 100 x 143 | 0.0566 (+0.0053/-0.0062) (Eq. 3); headline 0.059 +- 0.006 | tracer variants "extremely consistent"; "we start to be limited by the quality of the dust template" |
| 2005.05600v1 Natale 2020, FT-part (first 100k) | LFI 70 + WMAP Ka, Q, V, nside 16, pixel Gaussian, l 2 to 29 | 0.069 (+0.011/-0.012), 54 percent sky (Eq. 15); WMAP only 0.064 (+0.017/-0.015) | scalar templates, coefficient variance in the covariance (Eq. 5); "a clear trend towards high values of τ for large sky fractions" |
| 2007.04997 NPIPE, ABS + FT-part | PR4 | "τ= 0.051 ± 0.006" | polarization prior "requires a measurable and correctable transfer function affecting large-scale CMB polarization" (Sect. 2.1) |
| 2010.01139, 2309.10034v2 Tristram 2021, 2024, FT and FT-part | LoLLiPoP on Commander-cleaned PR4 | 0.0577 +- 0.0056; 0.058 +- 0.006 (Eq. 14) | foregrounds "kept fixed in the simulations"; "we cannot include uncertainties arising from an imperfect foreground model" |
| 2103.14378 de Belsunce 2021, FT | SRoll2 momento | 0.0581 +- 0.0055; SRoll1 0.0507 +- 0.0063 | no synchrotron removal: 100 x 143 0.0582 to 0.0585, autos 0.064 to 0.072 (Table 4) |
| 2205.13968v1 de Belsunce 2022, FT | Planck only, two stage | 0.0598 +- 0.0059 (momento), 0.0599 +- 0.0066 (pixel) | tight index priors push tau "higher than the TC-cleaned ones by several σ" (Sect. 4.3) |
| 2205.10104 BeyondPlanck XII, FT | LFI + WMAP Ka to V, 225-mode KL Gaussian, f_sky 0.68 | 0.066 +- 0.013 (Table 3) | "the uncertainty of the fully marginalized case is 46 % larger than for white noise" (Fig. 6), 32 percent over astrophysical-only; no synchrotron prior test; about 2000 Gibbs samples needed |
| 2501.11904v2 Li 2025 (CLASS), FT | CLASS 90 x Planck, xQML nside 16 | 0.053 (+0.018/-0.019) (Eq. 12); 0.055 +- 0.013 (Eq. 13); CLASS x WMAP < 0.070 at 95 percent (Eq. 14) | scalar templates (Table 1, x 1e-2): Planck 100 with K 1.14 +- 0.07, with PR4 30 2.26 +- 0.14; CLASS 90 with K 1.24 +- 0.08, with 30 3.43 +- 0.17; transfer function "essentially zero for ℓ<5" |
| 2504.13254v2 Ilic 2025, FT | PR4 | 0.0576 +- 0.0060 +- 0.0006 (Eq. 6) | Planck l 6 to 25 "may be affected by systematic residuals (arising from foreground or instrumental contamination)" |
| 2603.22454v1 elica 2026, FT-part (first 100k) | SRoll2 100, 143 + Natale WMAP/LFI, HL cross-spectra | 0.0575 (+0.0048/-0.0058) (Eq. 10) | scalar templates; "mild tension (2.3σ)" with DESI DR2; "A deviation is observed at ℓ = 5" |
| 2504.16932 Sailer 2025, FT; PRL 136, 081002 (2026) | DESI DR2 BAO + PR4 high-l + lensing, no low-l | 0.090 +- 0.012 | "≃3−5σ tension with Planck low-ℓ polarization data when taken at face value"; "at present we have no evidence to doubt the CMB polarization measurements of τ" |
| 2602.06115 Hergt 2026, FT-part | tension quantification | CMB against BAO mostly below 2 sigma | "do not interpret this as evidence for a standalone τ-driven discrepancy" |
| 1801.01481 Watts 2018, ABS | CLASS forecast | within a factor of two of cosmic variance | "We use a power spectrum-based likelihood"; not pixel based |

Abstract-level only, for orientation (verify before citing): 2508.21069 (4 sigma), 2601.09644
(2.4 sigma), 2607.26373 (0.067 +- 0.011), 2504.21813 (0.091 +- 0.011 without low-l EE),
2606.05166 (compressed SRoll2 likelihood), 2505.24829 (likelihood modelling), 2606.26576
(second CLASS 90 GHz receiver, mapping speed up 41 percent).

### 3.2 Status of the tension

Published assessments range from "3 to 5 sigma at face value" (2504.16932) through 2.3 sigma
(2603.22454) to below 2 sigma (2602.06115). (inference) In naive quadrature against
0.090 +- 0.012, the WMAP plus LFI values are only 1.3 sigma (2005.05600) and 1.4 sigma
(2205.10104) away, so the tension is carried by HFI. Cosmoglobe DR1 offers no tau: 2303.08095
Sect. 6 says "estimating the optical depth of reionization will require about nine months of
continuous runtime" and Sect. 6.2 leaves parameters for future work.

### 3.3 Attribution of tau uncertainty or bias to foreground modelling

- **BeyondPlanck XV, exactly.** Sect. 4.4 of 2011.08503v2: a prior of width 0.2, which lets the
  high-latitude index reach about -2.7, "leads to large foreground excesses", "an obvious
  foreground bias in the large-scale CMB extraction", and "a corresponding excess in estimates of
  the optical depth of reionization", citing BP XII. There is no tau number, table or figure for
  this test. BP XII (2205.10104) does not contain the words synchrotron or beta_s in any such
  test. The claim is therefore qualitative and unquantified.
- **Quantified, on older or different data.** 0811.4280: tau 0.090 +- 0.019 against 0.086 +- 0.016
  template cleaned; prior -2.7 +- 0.3 gives 0.092 +- 0.019; a curved-synchrotron simulation fitted
  with a power law returns 0.111 +- 0.019 for input 0.1. 2205.13968: several sigma for tight
  priors, Planck only.
- **Tracer dependence.** 1801.01226 (table above) is the largest documented effect and is
  attributed to an LFI systematic, not to the index.
- (inference) Template coefficient ratios, 30 GHz over K, expected (28.4/22.8)^3.1 = 1.98
  (1.96 in thermodynamic units): Planck 2018 V at 100 GHz 1.54 +- 0.13; Li 2025 with SRoll2
  100 GHz 1.98 +- 0.17; Li 2025 CLASS 90 GHz 2.77 +- 0.23 (filtered templates). Errors are naive
  ratios of the quoted errors. Template noise attenuation may contribute. No conclusion is drawn
  beyond showing that a joint model with an explicit index would be informative here.
- (inference) BP XII Fig. 6 numbers imply astrophysical marginalization alone adds about
  11 percent to sigma(tau) (1.46/1.32) under tight index priors.

### 3.4 A pixel-based likelihood with spatially varying indices and tau together

Not found as a joint posterior. Two-stage versions exist: 0811.4280 (WMAP5, 48 index pixels,
marginal CMB map and covariance passed to a tau likelihood), 1103.2554 (simulated Planck,
"exact pixel likelihood given the estimated CMB map and CMB covariance"), 2205.13968 (Planck,
fixed kernel, MAP plus Laplace), and 2205.10104 (Gibbs CMB samples compressed to a Gaussian). None
uses WMAP K with LFI on current maps, none infers smoothness, none reports coverage.

## 4. Part C: public data (landing pages only)

| Dataset, URL | Format and resolution | Noise products | Large-scale caveats |
|---|---|---|---|
| WMAP9, `https://lambda.gsfc.nasa.gov/product/wmap/dr5/m_products.html` | res 4 (nside 16) nine-year band Q,U `wmap_band_qumap_r4_9yr_{K,Ka,Q,V,W}_v5.fits`; nested, Galactic, mK | `wmap_band_quninv_r4_9yr_{band}_v5.fits`, 6144 x 6144, 144 MB each, sigma0^2 divided out, singular in processing-mask pixels (`ninv_info.html`, as extracted); loss-imbalance templates (20 files, 56.3 MB each) | N^-1 includes loss-imbalance modes; no band-to-band correlation stated. LOCAL: Ka, Q, V in `~/data/wmap9_lowres/`, K band not held |
| Planck PR4 low resolution, `https://data.cmb-s4.org/planck_pr4-lowres.html` | 63 files, 155.4 GiB; `npipe6v20_{freq}_map_ns0016.fits` plus two subset maps | dense `npipe6v20_ncm_ns0016_smoothed_{freq}_bin.dat` and inverse, 648.0 MiB each; nside 32 unsmoothed, 10.1 GiB each | cosine smoothing; transfer function at 44 to 217 GHz; 600 simulations need a NERSC account; NCM scaling step in NPIPE Sect. 4.2 not fully read |
| Planck 2018, PLA likelihood page | `COM_Likelihood_Data-baseline_R3.00.tar.gz`, bflike extension, SimAll lowE | projected quantities only | no nside 16 HFI maps or simulations found; last 64k chars of the `Frequency_maps` page unread |
| SRoll2, `https://wiki.cosmos.esa.int/planck-legacy-archive/index.php/Sroll2` | 28 nside 2048 maps `HFI_SkyMap_<freq>-<detset>_2048_RSR.00_<ringset>.fits` | 500 simulations each for 100, 143, 353 psb and 300 for 217, at nside 32, three tar files of about 1.4 GB | project site `sroll20.ias.u-psud.fr` refused connections on 2026-11-09 |
| Natale WMAP + LFI dataset, `http://www.fe.infn.it/u/pagano/low_ell_datasets/wmap_lfi_legacy` | from the arXiv comment of 2005.05600 | not verified | all fetch attempts failed on 2026-11-09 |
| BeyondPlanck v2, `https://beyondplanck.science/products/files_v2/` | chains 4 x 1.6 TB; `BP_c000x_Presamp_v2.h5` (2.2 to 2.8 GB); `BP_CMB_likelihood_v2.tar.gz` (18 MB) | CMB-level only | no band covariances or per-sample band maps; tarball contents not inspected |
| Cosmoglobe DR1, `https://www.cosmoglobe.uio.no/products/cosmoglobe-dr1.html` | IQU maps nside 512 (70 GHz at 1024), half-mission sets, chains 2 x 841 GB, `synch_beta/` | rms only | no dense low-resolution covariance, no per-sample maps, no tau product; 500 samples total (2303.08095 Sect. 3.3) |
| QUIJOTE MFI DR1, `https://lambda.gsfc.nasa.gov/product/quijote/quijote_mfi_data_get.html` | nside 512, mK_CMB, 11/13/17/19 GHz, beams 39 to 56 arcmin | weights and QU covariance columns; half-mission maps; no simulations | FDEC filter (l <~ 15); 5 to 6 percent calibration; dec -12 to 70 deg |
| S-PASS DR1, `https://lambda.gsfc.nasa.gov/product/foreground/fg_s_pass_get.html` | nside 1024, 8.9 arcmin, 2303 MHz, dec < -1 deg | sensitivity map; RM and RM error maps | ground residual rms 0.40 mK (Q), 0.51 mK (U) (1903.09420 Table 2); "Flux density calibration accuracy 5%" |
| CLASS, `https://lambda.gsfc.nasa.gov/product/class/class_prod_table.html` | 40 GHz nside 128; combined CLASS + WMAP Q map with 200 noise simulations; 90 GHz spectra only | noise simulation tarballs | filtered; only a C_l-domain transfer matrix; no 90 GHz maps |
| C-BASS | not public | none | no data page at `https://cbass.site.ox.ac.uk/`; `https://lambda.gsfc.nasa.gov/product/cbass/` is 404; release cited as "Taylor et al., in preparation" in 2605.23489 |
| Synchrotron revisited, `https://lambda.gsfc.nasa.gov/product/planck/p1/planck_revisited_pol_synch_map_get.html` | `synchrotron_map_QU_60arcmin.fits` (100.7 MB) | noise and scaling alm archives (2.2 GB, 505 MB) | assumes index -3.1; simulation counts not stated |

## 5. Flags against existing repo notes

`/Users/trygvels/code/Commander4/nextgen/LITERATURE.md`:
- Line 10, "~30% larger errors": BP XII Fig. 6 gives 46 percent over white noise and 32 percent
  over astrophysical-only; say which.
- Line 11, "widening the prior biases tau": one qualitative sentence, no number (Sect. 3.3).
- Line 16: add the DR1 I against DR1 III prior inconsistency (Sect. 2.2).
- Line 17: -3.07 +- 0.07 is the K/30 all-sky weighted mean only; the chain draws from the prior.
- Line 46: the method is a Q/U correlation fit in nside 8 superpixels, with a documented
  K against 30 GHz systematic.
- Line 47: the -3.13 is NPS dominated; Table 1 gives -3.02 +- 0.02 for regions 1 to 14.
- Line 48: this is Shi et al. 2024; PTE < 0.001 holds for one pair only (0.01 for the other);
  the two index maps correlate at r = 0.04. Eimer 2024 is 2309.00675.
- Line 49: "~10x" is the maximum; typical per-pixel gain is 2 to 4.
- Line 50: pre-publication C-BASS data, no beta_s measured.
- Line 90: now PRL 136, 081002 (2026). Line 94: add the NPIPE low-l transfer-function caveat.
- Missing: 2610.31469, 2510.12860, 2403.18123, 2301.05117, 2301.05113, 1909.05923, 1802.01145,
  1903.09420, 2309.00675, 2110.12803, 0811.4280, 1103.2554, 1801.01226, 1801.04945, 2005.05600,
  1908.09856, 2103.14378, 2010.01139, 1907.12875, 2205.10104, 2602.06115, 1801.01481.

`/Users/trygvels/code/Commander4/nextgen/DATA.md`:
- Lines 31 to 37 still say "verify" although lines 102 onward verify them.
- Line 78 (CMB-S4 Globus mirror returns 404) conflicts with lines 329 to 337; the direct links
  are still listed on the CMB-S4 pages on 2026-11-09; not tested here.
- Line 14 says about 1.6 TB per BeyondPlanck chain, line 59 says 1.8 TB; the page says 1.6 TB.
- No entries for the PR4 low-resolution set, SRoll2 maps and simulations, the Natale dataset,
  WMAP K-band res 4 N^-1 and loss-imbalance templates, or the BeyondPlanck likelihood tarball.

`/Users/trygvels/code/Commander4/nextgen/research/05_external_data_verified.md`: QUIJOTE, S-PASS,
CLASS and C-BASS entries are still correct on 2026-11-09. Its "(A)" read levels are superseded by
the full-text reads above.

## 6. Searches run and what was not found

arXiv API (`export.arxiv.org/api/query`); counts given where retained:
- `all:"C-BASS"` (22); `all:"C-BASS" AND (polarization OR polarisation OR survey)`.
- `(Cosmoglobe OR BeyondPlanck OR Commander) AND abs:"optical depth"` (5);
  `all:Cosmoglobe AND (reionization OR "optical depth" OR likelihood)` (1).
- `abs:synchrotron AND abs:"spectral index" AND polarized terms AND ("Gaussian process" OR
  hierarchical OR "spatial prior" OR smoothness OR "correlation length" OR "random field")` (7).
- Same stem with `("simulation-based calibration" OR "coverage probability" OR "posterior
  coverage" OR "credible interval" OR "rank statistics")`.
- Same stem restricted to 2026-01-01 to 2026-11-09 (12; found 2610.31469), and from 2025-01-01.
- `abs:"optical depth" AND abs:reionization` with: large-scale or low-multipole terms from
  2025-07-01 (13); `(foreground OR synchrotron) AND polarization`; `("spectral index" OR
  "component separation")`; `("spatially varying" OR "pixel-based" OR Gibbs) AND (synchrotron OR
  foreground)`; `synchrotron AND "spectral index"`; `(pixel OR "map-based")`; and
  `foreground AND (forecast OR simulated) AND (synchrotron OR "spectral index")` (2: 1507.08853,
  1103.2554).
- INSPIRE citing lists of 2504.16932 (117), 2603.22454 (4), 2205.13968 (15), 2205.10104 (22),
  2005.05600 (25).
- Web search for the SRoll2 data site and for a C-BASS release.

Not found:
- any public C-BASS polarization release, or a C-BASS polarized beta_s;
- any Cosmoglobe tau;
- any tau number for a widened beta_s prior in BP XV or BP XII;
- any spatial prior with inferred smoothness on an index map;
- any simulation-based coverage validation of a beta_s map;
- any joint sampler of tau with spatially varying indices;
- any public dense low-resolution covariance or per-sample maps for Cosmoglobe DR1;
- CLASS 90 GHz maps; QUIJOTE noise simulations; the Weiland 2022 composite index map on LAMBDA.

## 7. Gaps and follow-up

- Resolve QUIJOTE XII Table 4 (+- 0.03 against the abstract's +- 0.05) and its calibration
  treatment from the PDF. Check its INSPIRE citing list in a few weeks.
- Planck 2018 XI App. C per-region indices; remaining 127k characters unread.
- NPIPE Sect. 4.2: the scaling step required for the low-resolution covariances.
- Whether the public `CG_param_c0001_v1.txt` equals the local `param_CG_v09.txt`.
- Contents of `BP_CMB_likelihood_v2.tar.gz` and of the Natale dataset (HTTP-only host).
- Whether the CMB-S4 Globus HTTPS links work (needs a header request from a shell).
- ACT DR6 and SPT-3G parameter papers were not read for their tau statements.
- Result counts were not retained for several searches listed without a number.
- Several abstract-level 2025 to 2026 tension papers are listed by ID only; read before citing.
