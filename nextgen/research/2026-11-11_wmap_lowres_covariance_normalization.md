# WMAP nine-year low-resolution Q/U covariance: normalization and published chi-square values

Date: 2026-11-11. Author: lit-scout agent. Conventions follow
`nextgen/research/2026-11-10_citation_verification_paperA1_round2.md` (Sect. 2a of that report already
covers how the matrices are constructed; it is extended here, not repeated).

Question: a three-component fit to the nine-year r4 (N_side 16) Q/U maps of K, Ka, Q, V with the
released dense N^-1 gives a minimum chi-square of 0.956 of expectation (3.4 sd low; 0.973 and 2.5 sd
with a LambdaCDM prior on the CMB), located on the faint sky. Read as a noise miscalibration this is
about 2 percent in rms. Is that known or expected?

Read levels:
- PDF: I read the stated pages of the arXiv (or LAMBDA) PDF as page images and transcribed the
  quotes and table values myself.
- HTML-S: full-text HTML (ar5iv or a LAMBDA page) read through the fetch tool's small model, which
  returns quotes of at most about 125 characters, so long sentences are given in pieces. "x2" means
  the passage was queried twice with different prompts and gave the same numbers.
- SEC: secondary only (a search-engine summary); not verified, do not cite.
- BIB: bibliographic record only (arXiv API).

Conventions: mathematics in ASCII (the approximately-equal and similar-to signs are both written
"~"); typographic dashes inside quotes written as hyphens; "(sic)" for misprints in the source;
"(inference)" marks my own arithmetic or synthesis, never a result of the cited source. For a
chi-square with n degrees of freedom I use sd = sqrt(2/n) for the reduced value (inference
throughout).

Access limits: arXiv PDFs above 10 MB could not be fetched (1507.02704, 1907.12875, 1303.5075,
2303.08095), so those four are HTML-S. The ar5iv rendering of 1303.5075 ends inside App. A, so its
App. E was not read. A&A returned HTTP 403. Larson et al. 2011 was read only to character 100 000
of 153 175 (its Sect. 5 not read).

## 0. Short answer and summary table

Not known, and not expected from the documentation. No WMAP paper or the Explanatory Supplement
states an accuracy for the normalization of the r4 polarization N^-1, and none reports a reduced
chi-square below one for pixel-space tests. The team's own pixel-space numbers are at or above one:
1.016 to 1.022 (three-year, cleaned Q and V, r4), 1.000 to 1.018 (three-year pixel model, Ka, Q, V),
1.009 +- 0.040 (nine-year null maps, r3). The only direct calibration of the r4 noise against the
high-resolution sigma_0 (three-year) found the r4 noise 1.5 to 5.7 percent HIGHER, the opposite sign
to the present finding. Two things do point the same way as the present finding: (a) the five-year
paper says of the Fisher errors built on these matrices "If anything, there is a slight tendency to
overestimate the uncertainties at higher l"; (b) two published re-analyses that use exactly the
nine-year raw r4 maps and matrices print V-band reduced chi-squares of 0.96 to 0.98 (and a WMAP
minus LFI 70 GHz half-difference of 0.973), each called consistent with unity and none discussed.
BeyondPlanck XV saw a "uniform and slightly negative bias at high Galactic latitudes" in its
chi-square map and attributed it to the count of degrees of freedom.

| Source | Read | Finding | Relevance |
|---|---|---|---|
| Jarosik et al. 2007, astro-ph/0603452v2 | PDF pp. 32-54 | r4 matrix built from a per-year stationary noise filter normalized to one at zero lag, with grouped observations; maps degraded from r9, "noise properties not precisely described by the noise matrix". r4 noise from year differences is 1.015 to 1.057 times the r9 sigma_0(I) (Table 5) | Q1, Q3: the only r4 calibration; sign opposite to a chi-square deficit |
| Page et al. 2007, astro-ph/0603450v2 | PDF pp. 22-32, 69-73 | Table 4: cleaned maps outside P06, chi2/nu = 1.096, 1.018, 1.016, 1.050 (Ka, Q, V, W), nu = 4534 | Q1: three-year values at or above one |
| Kogut et al. 2007, 0704.3991v1 | PDF pp. 3-9 | Table 2: pixel model outside P06, chi2/DOF = 1.018, 1.000, 1.017 (Ka, Q, V), DOF 4524 | Q1: closest team analogue of a per-pixel fit |
| Hinshaw et al. 2009, 0803.0732v2 | PDF pp. 5-10, 16-22, 47-54 | Sect. 6.1: Fisher errors describe the scatter at l >= 6, "slight tendency to overestimate the uncertainties at higher l"; Table 5 | Q1: the one team statement with the sign of the present finding |
| Jarosik et al. 2011, 1001.4744v1 | PDF pp. 3-16, 25-32 | Sigma^-1 scaled by sigma_0(I), sigma_0(Q,U) of Table 1; W-band measured over Fisher variance agrees at l >= 8, exceeds it at l <= 7; Sigma^-1 omits baseline and calibration noise | Q1, Q3 |
| Bennett et al. 2013, 1212.5225v3 | PDF pp. 16-31, 97-111; HTML-S rest | sigma_0 from r9 maps by a diagonal method, white-noise model; "the r9 values are adopted for all WMAP nine-year analysis". Sect. 7.1: r3 likelihood chi2 1321 for 1170 dof; null maps 1180 +- 47 | Q1, Q3 |
| Nine-Year Explanatory Supplement (Greason et al., v5.0, 2012-12-19) | PDF front matter, pp. 71-78 | Sect. 2.2.3: imbalance "projected out" of the matrices; processing-mask elements zero; Table 2.13 gives band sigma_0. Contradicts the LAMBDA page on the imbalance term | Q1 |
| LAMBDA DR5 pages (ninv_info, skymap_info, loss_imbal_template_r4_info, product list) | HTML-S | "a true 1-sigma uncertainty" for the raw-map matrices; sigma_0^2 divided out; r4 maps "created from the res 9 maps" | Q1 |
| Lattanzi et al. 2017, 1611.01123v1 | PDF pp. 1-6 | Raw nine-year r4 maps and matrices, K and 353 GHz templates: chi2/N_pix = 1.045, 1.004, 0.977 (Ka, Q, V; noise only) and 1.026, 0.983, 0.963 (signal plus noise) | Q2: V band 2 to 4 percent low, uncommented |
| Planck 2018 V, 1907.12875v2 | HTML-S x2 | Table 12: reduced chi2 1.004, 1.048, 0.970 (Ka, Q, V), f_sky 73 percent | Q2 |
| Planck 2015 XI, 1507.02704v3 | HTML-S x2 | Sect. 2.6: WMAP and LFI 70 half-difference reduced chi2 0.973 (PTE 0.84, 2786 dof); half-sum 1.116 | Q2 |
| Weiland et al. 2018, 1801.01226v3 | PDF pp. 3-9 | Table 2: KaQV with K template, chi2/dof = 754/734, 752/734, 1190/1162 at N_side 8 | Q2: no deficit at N_side 8 with signal covariance |
| Svalheim et al. 2023 (BeyondPlanck XV), 2011.08503v2 | HTML-S x2 | Sect. 5.1 footnote 5: slightly negative chi-square offset at high latitude, attributed to n_dof | Q2: same sign and place |
| Watts et al. 2023, 2303.08095v1 | HTML-S | No chi-square against the r4 matrices; Sect. 5.3 says the template treatment is not "a complete description of the full large-scale uncertainties" | Q2: published claim of the opposite sign |
| Natale et al. 2020, 2005.05600v1; Paradiso et al. (BP XII), 2205.10104; Eimer et al. 2024, 2309.00675; Fuskeland et al. 2014, 1404.5323v2; Vidal et al. 2015, 1410.4438v2 | HTML-S | Use the matrices or WMAP noise; none reports an absolute chi-square against them or a rescaling | Q2: nothing found |

## 1. Question 1: what the WMAP team says

### 1.1 What the nine-year r4 N^-1 files are documented to contain

Nine-Year Explanatory Supplement, Version 5.0, "December 19, 2012", "editor M. R. Greason, M. Limon,
E. Wollack, et al." (PDF, title page and Preface). Sect. 2.2.3 "Reduced Resolution Sky Maps" (PDF,
document p. 76):
- "These low resolution maps are supplied at HEALPix Res 4 (Nside=16; 3072 pixels) and are used
  primarily for foreground analysis on the polarization component of the WMAP data."
- "However, these noise covariance matrices do not fully characterize the pixel-to-pixel noise in
  these lower resolution maps. Therefore a accurate (sic) and much larger inverse covariance matrix
  is supplied for each map as a separate product."
- "The first, more common, form consists of a 6144x6144 matrix containing four 3072x3072 blocks
  representing Q-U covariance".
- "The elements corresponding to the Galactic plane pixels as defined by the standard processing
  mask are set to zero. herefore (sic), the matrices as supplied are singular. To make them
  non-singular, the rows and the columns corresponding to each zero element of the diagonal should
  be deleted."
- "A small effect due to transmission imbalance between the A side and the B side of the instrument
  has been projected out of the inverse covariance matrices as described in Section 3.5.1 of
  (Jarosik et al., 2007)."
- Product table (p. 77): "Nine Year Coadded Q,U Maps per Freq. Band" 5 files; "Inverse Covariance
  Matrices for Nine Year Coadded Q,U Maps per Freq. Band" 5 files; single-year per-DA Q,U maps and
  their inverse covariance matrices, 90 files each.
- Sect. 2.2.4 (p. 77): the imbalance "has been measured for each radiometer using WMAP flight data,
  with an uncertainty of 20% (Jarosik et al., 2011)"; "Section 3.5.1 of Jarosik et al. (2007)
  describes how these templates can be projected out of the inverse covariance matrices so that
  transmission imbalance uncertainty will not affect subsequent analyses." Page 78: "where Sigma^-1
  is the 6144x6144 QU inverse covariance matrix for a given DA in mK^-2 and v is the 6144 element
  QU transmission imbalance template for the DA in mK."
- Table 2.13 "Pixel Noise Coefficients: cleaned indicates foreground cleaned" (p. 72), sigma_0 in
  mK. Per DA, Q,U uncleaned: K1 1.435, Ka1 1.472, Q1 2.254, Q2 2.140, V1 3.324, V2 2.958, W1 5.912,
  W2 6.577, W3 6.958, W4 6.795. Per frequency band, I: K 1.429, Ka 1.466, Q 2.188, V 3.131,
  W 6.544; Q,U uncleaned: K 1.435, Ka 1.472, Q 2.197, V 3.141, W 6.560; Q,U cleaned: Ka 2.166,
  Q 2.641, V 3.339, W 6.832.

LAMBDA DR5 pages (HTML-S):
- `ninv_info.html`: "All the inverse covariance matrices have had sigma0^2 divided out of them."
  "These factors may be found in" (link to `skymap_info.html`). "The matrices include a
  contribution due to uncertainty in the loss imbalance coefficients" "In previous releases, this
  uncertainty was fully marginalized--equivalent to assuming no information" "In the latest release
  of the matrices for raw Stokes Q and U maps (those without foreground cleaning)," "a true 1-sigma
  uncertainty was incorporated into this term."
- `skymap_info.html`: "The low resolution (res 4, NSide=16) maps were created from the res 9 maps."
  "The Nobs vectors were replaced with the appropriate" N^-1 diagonal elements, "scaled by sigma0^2."
  Same sigma_0 tables as Table 2.13.
- `loss_imbal_template_r4_info.html`: same text as Supplement Sect. 2.2.4; no statement on whether
  the released nine-year matrices already contain the treatment.
- `maps_band_r4_quninv_9yr_get.html`: file list only (five files of 144.00 MB, 2012-07-29); no
  description of normalization.

Not settled by these documents:
1. The Supplement says the imbalance modes are projected out; the LAMBDA page says the raw-map
   matrices carry a 1-sigma term instead. The two official texts disagree. (inference) The LAMBDA
   wording is the more specific and matches later practice (Weiland et al. 2018 fit the templates
   explicitly); the Supplement sentence reads as carried over from earlier releases.
2. Which sigma_0 goes with a band matrix. The pages say only that sigma_0^2 is divided out. Table
   2.13 gives band values that are the arithmetic means of the DA values (inference: Q (2.254 +
   2.140)/2 = 2.197; V (3.324 + 2.958)/2 = 3.141). Jarosik et al. 2011 Sect. 2.11 says the DA
   matrices are combined "scaled to reflect the noise level of each DA". How the band file is then
   brought to unit sigma_0 is not documented. (inference) If a band file were the plain sum of
   unit-sigma_0 DA matrices and were scaled by the band sigma_0, the noise variance would be
   overstated by 1.0 percent in V and 0.2 percent in Q and not at all in K and Ka: the right sign,
   but too small and absent in Ka. This can be checked on the files (Sect. 6).
3. Whether the QU matrix is the inverse of the QU block of the IQUS covariance or the QU block of
   its inverse (already open in the 2026-11-10 report). Not found.

### 1.2 How sigma_0 is set, and the stated accuracies

- Three-year, high resolution. Jarosik et al. 2007 Sect. 4.3.1 (PDF pp. 51-52): sigma_0 from
  year-to-year differences by two methods, results "agree within 0.3% with the values obtained
  using the first technique"; "the sensitivity of the radiometers has decreased at most by 0.1%
  over the 1 year period". Sect. 4.3 (pp. 47, 51): the diagonal approximation "uses different
  values of the noise per observation parameter for maps of the different Stokes components,
  parameterizing the effects of the ignored, off diagonal components of the noise correlation
  matrices."
- Three-year, r4. Jarosik et al. 2007 Sect. 4.3.2 "Noise Evaluation the r4 Full Sky Maps" (sic;
  PDF pp. 52-54): "The noise of the full sky reduced resolution (r4) maps may be determined using
  the pixel-pixel noise matrices. This analysis is performed only on the high Galactic latitude
  regions, those pixels originating from degradation of the high resolution spm maps." Eqs. (46)
  and (47): (sigma_0^{x,y})^2 = (t_x - t_y)^T Sigma_{x,y}^-1 (t_x - t_y), Sigma_{x,y}^-1 =
  (Sigma_x + Sigma_y)^-1, "where t_x and t_y are reduced resolution sky map sets comprising I, Q, U
  and S components". "The ratio of these average values to the noise values of the r9 I maps are
  presented in Table 5. The differences between the two noise estimates arise from ignoring the
  pixel-pixel noise correlations in both the r9 - r4 map degradation, and the evaluation of the
  sigma_0 of the r9 maps."
- Table 5 "Noise Levels for r9 Maps" (p. 53), columns sigma_0^I (methods 1, 2), sigma_0^P (method
  1), sigma_0^Q, sigma_0^U (method 2), sigma_0^r4/sigma_0^I:
  K1 1.439 1.439 1.455 1.453 1.453 1.057; Ka1 1.464 1.462 1.483 1.484 1.480 1.019;
  Q1 2.245 2.244 2.269 2.273 2.272 1.027; Q2 2.135 2.140 2.156 2.157 2.156 1.016;
  V1 3.304 3.311 3.330 3.333 3.333 1.023; V2 2.946 2.951 2.970 2.970 2.974 1.015;
  W1 5.883 5.896 5.918 5.935 5.930 1.023; W2 6.532 6.550 6.571 6.587 6.578 1.024;
  W3 6.885 6.885 6.925 6.918 6.931 1.020; W4 6.744 6.765 6.780 6.793 6.800 1.033.
  Note to the table: "The values in the last column are the ratios of the low resolution (r4) maps
  noise estimates, obtained using the noise correlation matrices (4.3.2), to the high resolution
  noise of the I maps (4.3.1). Differences in the noise estimates arise from ignoring pixel-pixel
  noise correlations both in the map degradation and the noise estimation of the high resolution
  maps."
  (inference) Relative to the r9 polarization value the ratio is 1.045 (K1), 1.006 (Ka1), 1.016
  (Q1), 1.006 (Q2), 1.015 (V1), 1.007 (V2). So in the three-year data the r4 noise measured with the
  matrices was HIGHER than the r9 sigma_0, by 0.6 to 1.6 percent in Ka to V. A matrix scaled with
  the r9 sigma_0 would then give chi-square above one. The estimate uses I, Q, U and S together, so
  it is not a pure Q/U number.
- Which sigma_0 the released three-year matrices carry: "The covariance matrices are the
  correlation matrices scaled by a multiplicative factor, 1/sigma_0^2, where sigma_0 is the noise
  per observation corresponding to each DA." (Sect. 3.5.1, p. 40). It does not say whether this is
  the r9 or the r4 value. Not found.
- Five-year. Hinshaw et al. 2009 Table 1 (PDF p. 49), sigma_0(Q,U) in mK: K1 1.453, Ka1 1.488,
  Q1 2.278, Q2 2.163, V1 3.341, V2 2.975; note d: "Noise per observation for resolution 9 and 10
  I, Q, & U maps, to ~0.1% uncertainty."
- Seven-year. Jarosik et al. 2011 Table 1 (PDF p. 5), sigma_0(Q,U) uncleaned: K1 1.456, Ka1 1.490,
  Q1 2.280, Q2 2.164, V1 3.348, V2 2.979 (sigma_0(I): 1.437, 1.470, 2.254, 2.140, 3.319, 2.955).
  Sect. 2.7 (p. 12): "In practice the terms M_am^T N^-1 M and M_am^T N^-1 M_am are evaluated at r4
  with the N^-1 normalized to unity at zero lag", and the result is "converted to noise units by
  dividing by the appropriate combinations of sigma_0(I) and sigma_0(Q, U) given in Table 1." So the
  seven-year r4 matrices are normalized with the high-resolution sigma_0. No r4 calibration like the
  three-year Sect. 4.3.2 appears in pp. 3-16 or 25-32.
- Nine-year. Bennett et al. 2013 Sect. 4.1.3 "Noise Characterization of the High Resolution Maps"
  (PDF pp. 24-26): "The noise in the r9 and r10 maps is described assuming the radiometer noise
  distribution is stationary, has a white spectrum and is normally distributed." "Ideally the value
  of sigma_0 is obtained by evaluating sigma_0^2 N_pix = <t_n^T Sigma^-1 t_n>," "but such a
  calculation is intractable with high resolution maps. In practice only the diagonal elements of
  Equation (15) are evaluated." "The omission of the off-diagonal terms in the evaluation of
  Equation (15) ignores the contribution to the noise from the reference beam pixels in the
  evaluation of sigma_0." "The sigma_0 values from the r10 maps have values from 0.3% (W-band) to
  1.5% (K-band) higher than those obtained form (sic) the corresponding r9 sky maps." "A more
  accurate determination of sigma_0 can be made by equating the diagonal elements of Equation
  (14)". "The values of sigma_0 obtained from r9 and r10 Stokes I maps, evaluated using diagonal
  elements of Equation (14), agree to ~0.05% with the worst discrepancy being ~0.1%." "The N_obs
  fields of the nine-year r9 and r10 sky maps now contain the reciprocals of the diagonal element of
  the Sigma matrix". "Because the sigma_0 values computed from r10 maps differ by less than 0.1%
  from those computed from r9 maps, the r9 values are adopted for all WMAP nine-year analysis."
  "These methods have been extended and applied to the Stokes Q and U maps and the spurious response
  map S. This change improved the agreement between the sigma_0 values for the temperatures and
  polarization maps to ~0.5% from ~1.1% in earlier data releases."
  No section of Bennett et al. 2013 treats the noise of the r4 maps (PDF pp. 16-31; HTML-S for the
  rest; same conclusion as the 2026-11-10 report).
- (inference) Nine-year over seven-year sigma_0(Q,U): 0.986 (K1), 0.988 (Ka1), 0.989 (Q1), 0.989
  (Q2), 0.993 (V1), 0.993 (V2). If the r4 matrices are scaled with these, as in the seven-year
  recipe, the nine-year normalization assigns 1.4 to 2.9 percent LESS noise variance than the
  seven-year one. The nine-year change therefore goes against, not towards, a conservative matrix.

Answer to "how accurately normalized": the stated accuracies (0.3 percent, about 0.1 percent,
about 0.05 percent) all refer to sigma_0 of the r9 and r10 maps. No statement gives the accuracy of
the r4 N^-1 normalization for the five-, seven- or nine-year release. The one r4 measurement
(three-year) differs from the r9 value by 1.5 to 5.7 percent.

### 1.3 Chi-square values the team reports at low resolution

Three-year, r4, template-cleaned maps. Page et al. 2007 Sect. 4.3 (PDF pp. 23-25): "In both cases
the full pixel covariance matrix was used to compute chi^2 for Stokes Q and U simultaneously. For
the full sky the number of degrees of freedom, nu, is 6144 (twice the number of pixels in an r4
map) and outside the P06 mask nu = 4534." "Since chi^2/nu for Q and V bands is so close to unity
for the cleaned maps, it is no longer an effective measure of cleaning." Table 4 "Comparison of
chi^2 Between Pre-cleaned and cleaned Maps" (p. 72), columns chi2/nu pre-cleaned, chi2/nu cleaned,
nu, Delta chi2:
- full sky: Ka 10.65, 1.20, 6144, 58061; Q 3.91, 1.09, 6144, 17326; V 1.36, 1.19, 6144, 1045;
  W 1.38, 1.58, 6144, -1229;
- outside P06: Ka 2.142, 1.096, 4534, 4743; Q 1.289, 1.018, 4534, 1229; V 1.048, 1.016, 4534, 145;
  W 1.061, 1.050, 4534, 50.
The text adds (pp. 24-25) "the reduced chi^2 in the Q and V band maps is 1.022 and 1.016 as
compared to 1.014 and 1.018 for the starlight-based directions" (the second pair does not match the
table's 1.018 and 1.016; as printed). The chi-square is "relative to a null signal", the cleaned
covariance is scaled "by a factor of 1/(1 - a_0)^2" and the K-band pattern is projected out.
(inference) sd = 0.021, so Q and V are 0.8 sd above one; CMB and residual foregrounds add to it.

Three-year, r4, pixel model. Kogut et al. 2007 (arXiv:0704.3991v1; PDF p. 6), Eq. (9): chi^2 =
Sum R N^-1 R, "where N^-1 is the inverse noise matrix and R_i,alpha is the (data-model) residual in
pixel i and Stokes parameter alpha. The synchrotron amplitude P_s and polarization angle gamma_s in
the pixel model are derived almost entirely from the K-band data, leaving nearly zero residuals at K
band. Since the synchrotron amplitude and polarization angle may be treated as an external (K-band)
template, we assess the goodness of fit using only the 4 higher-frequency bands, counting only the
parameters fitted using data at these bands." Table 2 "chi^2 of Residual Maps" (chi2, DOF,
chi2/DOF):
- Pixel Model, Outside P06: Ka 4604, 4524, 1.018; Q 4521, 4524, 1.000; V 4600, 4524, 1.017;
  W 4774, 4524, 1.055; All 18505, 18098, 1.023.
- Pixel Model, Inside P06: Ka 1626, 1601, 1.016; Q 1541, 1601, 0.963; V 1670, 1601, 1.043;
  W 1773, 1601, 1.107; All 6611, 6404, 1.032.
- Template Model, Outside P06: Ka 4968, 4532, 1.096; Q 4616, 4532, 1.019; V 4608, 4532, 1.017;
  W 4760, 4532, 1.050; All 18952, 18128, 1.045.
- Template Model, Inside P06: Ka 2403, 1608, 1.494; Q 2081, 1608, 1.294; V 2704, 1608, 1.682;
  W 4946, 1608, 3.076; All 12034, 6432, 1.871.
"Outside the P06 mask, the variance of the cleaned maps is close to the value expected for
instrument noise alone." The paper does not comment on the level of chi2/DOF.

Five-year, l-space internal scatter. Hinshaw et al. 2009 Sect. 6.1 (PDF pp. 19-20): "We have
evaluated these individual spectra from the 5 year data and have assigned noise uncertainties to
each estimate using the Fisher formalism described in Page et al. (2007). We subject the ensemble
to an internal consistency test by computing the reduced chi^2 of the data at each multipole l
within each band or band pair, under the hypothesis that the data at each multipole and band
measures the same number from DA to DA and year to year." "For l >= 6, the most significant
deviation from 1 in reduced chi^2, in any spectrum or band, is 1.594 in the l = 7 BB spectrum for
KaQ." "Thus we conclude that the Fisher-based errors provide a good description of the DA-to-DA and
year-to-year scatter in the l >= 6 polarization data. If anything, there is a slight tendency to
overestimate the uncertainties at higher l." "For l <= 5, we find 37 out of 120 points where the
reduced chi^2 deviates from 1 at more than 4 sigma significance, indicating excessive internal
scatter in the data relative to the Fisher errors. However, all but 5 of these occur in cross-power
spectra in which one or both of the bands contain W band data." "For Ka-V bands, we believe that
the Fisher errors provide an adequate description of the scatter in this l <= 5 polarization data".
Table 5 "Polarization chi^2 Consistency Tests" (pp. 52-53), from "the template-cleaned maps";
columns KaKa (10), KaQ (50), KaV (50), KaW (100), QQ (45), QV (100), QW (200), VV (45), VW (200),
WW (190), the bracket being "the number of degrees of freedom in the reduced chi^2". Rows l = 6 to
10 for the six pairs without W (KaKa, KaQ, KaV, QQ, QV, VV):
- EE: l6 1.003 0.847 0.688 0.740 0.856 1.384; l7 0.600 0.671 0.689 0.936 0.780 0.900;
  l8 1.578 1.262 1.337 1.080 0.763 1.025; l9 0.760 0.710 0.891 0.582 0.726 0.791;
  l10 0.494 0.821 0.996 0.656 0.763 0.676.
- EB: l6 0.475 1.183 0.651 0.829 1.023 1.201; l7 1.014 1.007 0.829 0.817 0.759 0.616;
  l8 0.849 0.897 1.279 0.681 0.689 1.021; l9 0.743 0.734 1.007 0.820 0.798 0.882;
  l10 0.413 1.003 1.316 0.722 0.900 1.124.
- BB: l6 0.379 0.805 0.483 1.009 0.861 0.800; l7 1.925 1.594 0.967 1.074 0.817 0.772;
  l8 0.804 1.005 0.999 1.069 0.782 0.997; l9 0.320 0.489 0.502 0.884 0.491 0.748;
  l10 1.181 1.162 1.028 1.218 1.165 1.079.
(inference) The dof-weighted mean of these 90 entries is 0.879 (EE 0.849, EB 0.893, BB 0.896).
Caveats: the bracketed dof equal the number of spectra, so if the mean was fitted the expectation
is about 0.98, not 1; spectra within a band pair share maps and are not independent; these are
five-year cleaned maps and five-year matrices. For noise-dominated cross-spectra the variance
scales as the product of two noise variances, so 0.88 corresponds roughly to 6 percent in map-noise
variance, 3 percent in rms. This is the same order as the present finding. It is my arithmetic on
their table, not a number the paper states.

Seven-year, W band only. Jarosik et al. 2011 Sect. 4.1.2 (PDF pp. 25-31): cross spectra were
calculated "utilizing the full Sigma^-1 covariance matrix to weight the input maps. The
polarization analysis mask was applied by marginalizing the Sigma^-1 over pixels excluded by the
mask". Fig. 7 caption: "The bottom panel is the ratio of the measured variance between the
individual power spectra estimates (DA year x DA year) and the variance predicted by the Fisher
matrix calculation using the full Sigma^-1 inverse noise covariance matrix. Note the good agreement
for l >= 8." Text (p. 29): "For l >= 8 there is good agreement between the measured variance and
analytical estimate, but the measured variance slowly grows larger with decreasing l for l <= 7 for
both the cleaned and uncleaned sky maps." No equivalent figure for Ka, Q or V. (My reading of the
plotted points for l = 8 to 20: scatter between about 0.7 and 1.4, no number in the text.)

Nine-year, r3, cleaned KaQV. Bennett et al. 2013 Sect. 7.1 (PDF pp. 106-107): "The low-l polarized
pixel-based likelihood contains 585 unmasked res 3 pixels each with a Q and U Stokes parameter, for
1170 degrees of freedom. The chi^2 value for this part of the likelihood is 1321. The probability
to exceed this chi^2 value is 0.13%, which is unusually low." "we computed chi^2 from six
combinations of difference maps involving Ka-, Q-, and V-bands: Ka-Q, Ka-V, Ka-QV, Q-V, Q-KaV, and
V-KaQ". "We find that none of these combinations show an anomalous chi^2. The average and standard
deviation of chi^2 is 1180 +- 47 for 1170 degrees of freedom. The largest value of chi^2 is 1236
from Ka-QV, and the probability to exceed (PTE) is 8.8%." "We find sigma_r3 = 0.27 uK brings the
reduced chi^2 to unity. The instrumental noise per r3 pixel of the co-added KaQV map ranges from
0.43 to 1.57 uK, with the average and standard deviation of 0.86 +- 0.17 uK." "the reduced chi^2"
"drops from 1.13 to 1.12, 1.10, and 1.09" when the mask is enlarged. "Therefore, we conclude that
the excess chi^2 likely to be (sic) at least partially due to residual foregrounds, which we do not
include in the noise covariance matrix."
(inference) The null maps give 1.009 with sd 0.041 per combination. A 4 percent deficit is one sd
of a single combination, so this test neither shows nor excludes it. It is for cleaned maps at r3
with the cleaned-map matrices, not the raw r4 matrices.

Other team papers:
- Gold et al. 2011 (1001.4555v3, HTML-S, characters 0 to 120 000 of 144 616), Sect. 4.5,
  polarization-only fit at N_side 16: "Residuals after the fit are small compared to the noise, and
  over all bands the mean reduced chi2 per pixel is 1.1." Sect. 2.4 (polarization template
  cleaning): no chi-square value and no noise weighting described.
- Gold et al. 2009 (0803.0715v2, HTML-S): no polarization-only chi-square; the pixel fits use
  N_obs noise and "an overall multiplicative rescaling factor was determined for each frequency
  band" for smoothed maps.
- Dunkley et al. 2009, Gibbs paper (0811.4280v2, HTML-S): chi-square values are given for
  simulations only (for example "chi2=1497 (chi2/pixel=0.98)" for the CMB map without the dust
  prior); none for the flight data.
- Dunkley et al. 2009, likelihood paper (0803.0586v2, HTML-S, all 188 285 characters in two
  portions): no chi-square for the low-l polarization maps.
- Larson et al. 2011 (1001.4635v2, HTML-S, first 100 000 characters): no pixel-space chi-square for
  low-l polarization. Sect. 5 not read.
- Bennett et al. 2013 Sect. 5.3.2.2 (HTML-S): polarization cleaning, "none of them gave significant
  improvement in chi2"; Table 11 lists coefficients only.

### 1.4 Statements that the matrices carry extra variance

- Transmission imbalance. Three-, five- and seven-year matrices have one or two modes per DA
  projected out (Jarosik et al. 2007 Sect. 3.5.1, Eqs. 29 and 30, PDF pp. 39-40; Hinshaw et al.
  2009 Sect. 6.2, PDF p. 22: "we marginalize over residual uncertainties in the imbalance
  coefficients when we form the pixel-pixel inverse covariance matrices"; Jarosik et al. 2011 Sect.
  2.9). Nine-year: see the contradiction in Sect. 1.1. (inference) Either way this touches at most
  two modes per DA and year, a few tens of modes in a band matrix of 6144. It cannot produce a 4
  percent deficit over several thousand degrees of freedom.
- Poorly constrained modes. Jarosik et al. 2011 Sect. 4.1.1 (PDF p. 25): "they are properly
  de-weighted when these maps are analyzed using their corresponding Sigma^-1 matrices". Sect. 2.11
  (p. 14): for K and Ka the smallest eigenvalue of the summed QU matrix is set to zero; "No modes
  are removed from the Q, V and W-band maps."
- Spurious map S. Jarosik et al. 2007 Sect. 3.6 (PDF p. 45): "Including the spurious term in the
  map degradation and likelihood evaluation was found to have a negligible effect on the likelihood
  results (Page et al. 2006)." No statement that the S degrees of freedom make the QU matrices
  conservative. Not found.
- Partially filled pixels. From the five-year release the processing mask is defined on r4 pixels:
  "The new mask is defined in terms of low-resolution (r4) HEALPix sky pixels (Gorski et al. 2005)
  to facilitate a cleaner definition of the pixel-pixel inverse covariance matrices, N^-1."
  (Hinshaw et al. 2009 Sect. 2, PDF pp. 7-8). No statement of extra variance from partial pixels.
- Not found anywhere: a sentence calling the low-resolution matrices conservative, or a rescaling
  of the polarization noise for the raw maps.

## 2. Question 2: independent analyses that used these matrices

Lattanzi et al. 2017 (arXiv:1611.01123v1, JCAP DOI 10.1088/1475-7516/2017/02/041; PDF pp. 1-6).
Sect. 2.1: "we consider the WMAP9 raw low resolution frequency maps for the Ka, Q, and V channels
and the related noise covariance matrices that are publicly available." Each map "is masked by the
WMAP9 processing mask", and two scalings are fitted by minimizing chi^2_S+N = m~^T C^-1 m~ with
templates "chosen as the WMAP9 K-band and the Planck 353 GHz maps"; "The covariance matrix C = S +
N is taken as the sum of the signal, S, and noise, N, matrices". "We nonetheless also report" "the
chi^2_N figures we find assuming C = N only." "The reduced chi^2 values we find are consistent with
unity for all bands, in all masks considered." Table 1 (p. 5), columns alpha, beta, N_pix,
chi2_S+N/N_pix, chi2_N/N_pix, and the same two outside P06 ("leaving N_pix = 4510 good pixels"):
- Ka: 0.315 +- 0.003, 0.0031 +- 0.0004, 5720, 1.026, 1.045, 1.012, 1.029;
- Q: 0.163 +- 0.003, 0.0039 +- 0.0004, 5742, 0.9830, 1.004, 0.9778, 0.995;
- V: 0.047 +- 0.003, 0.0076 +- 0.0004, 5894, 0.9627, 0.977, 0.9606, 0.974.
(inference) N_pix equals twice (3072 minus the masked pixels of Bennett et al. 2013 Table 4: 212,
201, 125), so it counts Q and U. In sd: V -2.0 (S+N) and -1.3 (N); V outside P06 -1.9 and -1.2;
Q -0.9 and +0.2; Ka +1.4 and +2.4. The noise-only form should exceed one by the CMB term, so the V
value of 0.977 is low by at least 2.3 percent. Summed over the three bands the noise-only total is
1.008 and the S+N total 0.990. The model has one index per band over the sky, so foreground
mismatch adds chi-square in Ka and Q.

Planck 2018 V (arXiv:1907.12875v2, A&A 641, A5; HTML-S x2). Sect. 2.2.5: "As covariance matrices we
adopt the official WMAP covariance matrices and smooth them in harmonic space". Sect. 2.3.6: "The
Ka, Q, and V channel maps have first been smoothed with the same cosine window used for LFI" "then
foreground-cleaned individually on the WMAP analysis mask P06." "The estimated foreground
coefficients and associated chi2 values are given in Table 12." Table 12 (band, f_sky, alpha, beta,
reduced chi2): Ka 73 %, 0.318 +- 0.003, 0.0036 +- 0.0005, 1.004; Q 73 %, 0.161 +- 0.003, 0.0040 +-
0.0005, 1.048; V 73 %, 0.048 +- 0.003, 0.0079 +- 0.0005, 0.970. No degrees of freedom printed and no
comment on the values. (inference) With about 4500 elements, sd is about 0.021: V -1.4, Q +2.3.

Planck 2015 XI (arXiv:1507.02704v3, A&A 594, A11; HTML-S x2). Sect. 2.6 "Comparison with WMAP-9
polarization cleaned with Planck 353 GHz": Ka, Q, V cleaned with K and 353 GHz, "we use Nside=16";
"To exactly mimic the procedure followed by the WMAP team, we exclude the signal correlation matrix
from the noise component of the chi2 form." For the WMAP and Planck 70 GHz combination on "the
intersection of the WMAP P06 and Planck R1.50 masks (fsky=45.3%)": "Assuming that chi_sd^2 is chi2
distributed with 2786 degrees of freedom we find a PTE(chi2>chi_sd^2)=1.3x10^-5" "(reduced
chi2=1.116) for the half-sum, and PTE=0.84 (reduced chi2=0.973) for the half-difference." "This
strongly suggests that the latter is consistent with the assumed noise". (inference) 0.973 is 1.0
sd low; the statistic mixes the WMAP and LFI covariances.

Weiland et al. 2018 (arXiv:1801.01226v3, DOI 10.3847/1538-4357/aad18b; PDF pp. 3-9). "we obtain the
9-year coadded WMAP low-resolution (Nside = 16) Stokes Q and U maps and inverse covariance matrices
from LAMBDA"; "We further downgrade the Nside = 16 products to Nside = 8 using inverse noise
weighting". Likelihood with M = S + N. Table 2 "Derived values for tau, EE only", WMAP KaQV with
the K template, chi2_maxlike/dof: 754/734 (row 4, coefficients free, 46 %), 752/734 (row 5, fixed),
1190/1162 (row 6, 75 %). Also LFI 70 with LFI 30: 728/734. On the poorly measured modes: "within
the context of a low-multipole likelihood, where full covariance matrices are included, poorly
measured modes are correctly weighted in the analysis (Page et al. 2007)." (inference) 1.027,
1.025, 1.024 with sd 0.052 and 0.041: no deficit at N_side 8.

BeyondPlanck XV, Svalheim et al. 2023 (arXiv:2011.08503v2, A&A 675, A14; HTML-S x2). Data: WMAP Ka,
Q, V at N_side 16 with "the precomputed dense low-resolution noise covariance matrices provided by
the WMAP team", plus LFI and 353 GHz. Sect. 5.1, on Fig. 2 ("Normalized and reduced chi2 per pixel,
summed over Stokes Q and U and averaged over all Gibbs samples"), n_dof = 20 798: "Overall, we see
that our model appears to perform well at high Galactic latitudes," "while the Galactic plane shows
clear excess chi2." Footnote 5: "We note that it is very difficult to compute n_dof rigorously due
to the presence of active priors," "As such, the important feature of Fig. 2 is its spatial
structure, not the absolute zero-level." "The uniform and slightly negative bias at high Galactic
latitudes is thus simply an indication that" "n_dof is very slightly over-estimated;" "if it were
due to over-estimating the white noise rms level," "which is the only other possible explanation
for a chi2 deficit," "the Planck scanning strategy would have been visually apparent." No number
is given for the offset. This is the closest published precedent: same sign, same part of the sky,
explained by the count of degrees of freedom and argued only for the Planck noise, not for WMAP.

Cosmoglobe DR1, Watts et al. 2023 (arXiv:2303.08095v1; HTML-S, characters 60 000 to 310 693 in
four portions). No chi-square is computed against the official r4 matrices. Sect. 5.3 (quoted in
the 2026-11-10 report): "the template-based correction of the low-resolution covariance matrices"
"does not provide a complete description" "of the full large-scale uncertainties", that is,
understated errors on a few large-scale modes, the opposite sign. Table 5 gives template amplitudes
and the fractional decrease in standard deviation per DA (0.30, 0.25, 0.11, 0.20, 0.06, 0.16 for K1,
Ka1, Q1, Q2, V1, V2), not a chi-square. Sect. 4.4, time-domain white noise against WMAP: "about
half of the Cosmoglobe values lie between the two WMAP results, while the other half are higher."
"we have not yet been able to identify the cause of the major difference in reported white noise
levels" "it is more likely just a matter of different conventions". So DR1 does not confirm the
WMAP white-noise level at the percent level in either direction.

Others:
- Natale et al. 2020 (2005.05600v1, HTML-S): builds covariances "starting from the polarization
  pixel-pixel inverse covariance matrices at Nside=16 (Res 4)", smooths them, adds (20 nK)^2
  regularization. Only an excess chi-square in units of its dispersion is plotted; no absolute
  value in the text, no statement of a negative excess (figures not readable by me).
- BeyondPlanck XII, Paradiso et al. (2205.10104, HTML-S): chi-square PTE 0.32 for 225 modes of its
  own LFI plus WMAP likelihood; nothing on the WMAP matrices alone.
- Eimer et al. 2024 (2309.00675, HTML-S, first half): WMAP noise simulations "by sampling according
  to the full noise covariance matrix" (footnote names `wmap_band_quninv_r4_9yr_Ka_v5.fits`); no
  validation number.
- Fuskeland et al. 2014 (1404.5323v2, HTML-S): works at N_side 64 with white noise from the
  full-resolution maps ("Correlated pixel noise covariance matrices are only in Nside=16."); does
  not use the r4 matrices; no comment on their level.
- Vidal et al. 2015: the filament paper is arXiv:1410.4438v2 (DOI 10.1093/mnras/stv1328), not
  1410.4436 (which is "A new polarisation amplitude bias reduction method"). HTML-S: errors from 500
  noise realizations per band, and "We note that this is an underestimate of the uncertainty due to
  the presence of correlated noise in large scales." Does not use the r4 matrices.
- Planck 2013 XV App. E (1303.5075v2): SEC only. A search-engine summary gives reduced chi-squares
  of about 1.13, 1.13, 0.99 for Ka, Q, V with the 353 GHz template and says only diagonal noise
  elements were used. Not verified; do not cite.

Summary of Q2 (inference): with the nine-year r4 products, published per-band reduced chi-squares
span 0.96 to 1.05. V band is low in both tables that print it (0.961 to 0.977 and 0.970), by 1.2 to
2.0 sd each; these are the same data, not independent confirmations. No paper rescales the
covariance, and none calls it conservative or optimistic in the white-noise sense.

## 3. Question 3: are the matrices computed like the maps?

They are known approximations, documented in Jarosik et al. 2007 (PDF) and not revised in print
since.

- Time-domain noise model. Sect. 3.4.2 (p. 35), Eq. (23): the filter N_tt^-1(Delta t) is built
  from "the parametrized noise correlation function as described in 2.5", set to zero beyond
  Delta t_max; "a constant K is added to the values of N_tt^-1 for lags less than Delta t_max to
  force the mean of this portion of the filter to zero. This ensures that very low frequencies,
  those with periods longer than the filter extent, are filtered out. Finally, the filter is
  normalized by an overall multiplicative constant, C, such that the value at Delta t = 0 is unity."
  "The radiometer noise properties are treated as stationary for each year of data, making the
  inverse noise matrix, N^-1, circulant." Hinshaw et al. 2009 Sect. 2 (p. 7): "The functional form
  of the filter is unchanged (Jarosik et al. 2007), but the fits have been updated to cover years 4
  and 5 of the flight data."
- Evaluation at r4 with grouped observations. Sect. 3.5 (pp. 38-39): "These matrices were evaluated
  at r4 on a year-by year basis for all 10 DAs, resulting in 30 matrices. One approximation was
  used in the evaluation of these matrices to speed processing." "The approximation consists of
  grouping time contiguous observations for which both the A-side beam and B-side beam remain in
  the same r4 pixels. Each such group comprises ~30 rows of the mapping matrix (about 15
  observations)". "Using the group average values for the data in the columns of M rather than the
  exact row by row values in calculation of the inverse noise matrix is a good approximation for
  the vast majority of the sky, since the values being averaged are smooth and slowly varying.
  However these assumptions fail for pixels very close to the Galactic poles since the
  polarizations (sic) angles gamma_A and gamma_B vary significantly within a r4 pixel. The effect
  of these approximations has been investigated and is found to have negligible affect (sic) on the
  calculated power spectra."
- Maps are not made at r4. Sect. 3.6 (p. 40): "Since the matrices are produced at r4 due to
  computational constraints, corresponding low resolution sky maps must be produced that preserve
  the signal and noise properties as well as possible." "Even without application of a filter, the
  omission of the inter-pixel noise correlations results in the degraded maps having noise
  properties not precisely described by the noise matrix as calculated in 3.5." Page 42: maps made
  directly at r4 "will have noise properties accurately described by the pixel-pixel noise matrix,
  but will have a higher degree of signal aliasing"; "Low resolution maps of Stokes Q and U were
  produced both directly at r4 and by degrading r9 maps without application of an anti-alias
  filter. Analysis of both sets of maps produced very similar low-l polarization results for both
  simulated and flight data. Based on these results it was determined that the noise properties of
  the degraded form of low resolution polarization maps was adequately described by the pixel noise
  matrices, so the degraded form is used for all three Stokes parameters." "Although the
  inter-pixel noise correlations are ignored in the map degradation, the noise correlations between
  the Stokes Q and U signals within each high resolution pixel are included."
- Seven-year form. Jarosik et al. 2011 Sect. 2.10 (p. 13): "The low resolution sky maps were
  generated by performing an inverse noise weighted degradation of the high resolution maps that
  takes into account the intra-pixel noise correlations of the high resolution input maps." Sect.
  2.7, Eq. (14), evaluated "at r4 with the N^-1 normalized to unity at zero lag". Sect. 4.1.2
  (p. 29): the matrix "incorporates information regarding the sky scanning pattern, through the
  mapping matrices M and M_am, and correlations in the radiometer noise, through the use of inverse
  radiometer noise covariance matrix N^-1. This method does not, however, include the effect of any
  noise correlations that might be introduced through the baseline fitting and calibration
  procedure." Sect. 4.1.3 (pp. 29-31): with gain and baseline marginalized, "For E-mode
  polarization the maximum increase in uncertainty in a_lm is about 3.5% occurring for l = 3 and is
  less than 1% for other multipoles. For the B-mode a_lm uncertainties for l = 2, 3 increase by 40%
  and 90% respectively".
- Nine-year. Bennett et al. 2013 Sect. 2.2.1: "The standard WMAP map-making procedure is unchanged
  from the previous release" (quoted from the 2026-11-10 report, PDF p. 8). The new elements are
  the per-band processing masks (Table 4, p. 19: masked pixels K 312 or 270, Ka 212, Q 201, V 125,
  W 98) and the sigma_0 method of Sect. 1.2. Nothing on how the nine-year r4 matrices were
  recomputed.
- Flagged-data fraction and weights. Not found: no sentence says the r4 matrices use the same data
  flags as the maps. (inference) Since Sigma^-1 is accumulated from the mapping matrix of each
  year's observations, cut data would drop out of both, but the documents do not say so. The
  weights do differ by construction: the matrix is a direct r4 evaluation with the full filter, the
  map is an N_obs-weighted average of r9 pixels.

(inference) Sign of the documented approximations. Each one that can be signed makes the matrix
understate the noise of the released map: the reference-pixel and inter-pixel terms ignored in the
degradation (three-year measurement: r4 noise above the r9 sigma_0), the baseline and calibration
noise left out of Sigma^-1, and the residual-foreground and imbalance arguments of Bennett et al.
2013 and Watts et al. 2023. None predicts a chi-square deficit. What the documents leave open is
the absolute scale: the filter is normalized to one at zero lag by fiat, so the scale of Sigma^-1
rests entirely on sigma_0, and since the seven-year release that sigma_0 is taken from the r9 maps
under a white-noise model rather than measured at r4. A mismatch of 2 percent in rms in that
identification is neither documented nor excluded.

## 4. Flags against existing project notes

- `research/2026-11-10_citation_verification_paperA1_round2.md`, Sect. 2a, open point 3 ("The
  nine-year change from projection to a 1-sigma term is documented only on the LAMBDA page"): the
  Nine-Year Explanatory Supplement has now been read (Sect. 2.2.3) and it says the modes are
  "projected out", so it contradicts the LAMBDA page instead of supporting it. Paper A1 should not
  state the nine-year treatment as a fact without checking the files (for example whether
  N^-1 v = 0 for the released templates).
- Same report, open point 1 (mask): the Supplement repeats that processing-mask elements are zero
  and the matrices singular, against the E038 finding of positive definite band matrices. Still a
  documentation conflict; the file check stands.
- Same report, Sect. 6: its reading of Watts et al. 2023 Sect. 5.3 as "a published claim" "that
  those matrices understate the large-scale errors" stands, and is the opposite sign to the present
  deficit. The two are not in conflict: one concerns a few large-scale modes, the other the overall
  scale.
- `nextgen/DATA.md` lines 309-311 list sigma_0 per DA only. For the band-averaged r4 files the band
  values of Table 2.13 apply (Q,U uncleaned: K 1.435, Ka 1.472, Q 2.197, V 3.141).
- Task statement: "Vidal et al. 2015 (arXiv:1410.4436)" is the wrong identifier for the filament
  paper; it is 1410.4438.
- Nothing in `nextgen/LITERATURE.md` is contradicted.

## 5. Searches run, and what was not found

Searches and fetches (generic queries only):
- arXiv API (`export.arxiv.org/api/query?id_list=...`) for 22 identifiers (records of all papers
  above).
- INSPIRE API: `fulltext:"quninv"` (1 hit: 2309.00675); `refersto:arxiv:1212.5225 and
  fulltext:"inverse covariance matrices" and fulltext:"Nside=16" and fulltext:"WMAP" and
  fulltext:"reduced chi"` (response truncated to the first record, unusable).
- Web search: "WMAP nine-year low resolution polarization inverse noise covariance matrix reduced
  chi-square below one noise overestimated Ka Q V Nside 16"; "WMAP res 4 Q U "inverse covariance"
  matrices "sigma0" normalization chi-square template fitting 353 GHz "reduced chi" Ka Q V raw
  maps". Neither returned a source that says the matrices overstate the noise.
- Full-text reads listed in the table and sections above. The PDFs were read from the fetch
  tool's session cache, which is not kept.

Looked for and not found:
- A stated accuracy of the r4 N^-1 normalization in the five-, seven- or nine-year papers or the
  Nine-Year Supplement.
- A year-difference or DA-difference chi-square at r4 for raw Ka, Q or V maps in the seven- or
  nine-year papers (the three-year paper has it, for I, Q, U, S jointly).
- Any rescaling of sigma_0 for the raw r4 polarization maps, or any sentence calling the matrices
  conservative.
- Any independent paper that reports a chi-square deficit against these matrices as a finding, or
  rescales them.
- A description of how the nine-year r4 matrices were computed (same recipe as seven-year is
  assumed, not stated).

Not read: Larson et al. 2011 Sect. 5; Planck 2013 XV App. E; the three-, five- and seven-year
Explanatory Supplements; Gruppuso et al. and Paci et al. (BolPol analyses of WMAP low-resolution
polarization); the FITS headers of the matrix files.

## 6. Suggested follow-up

Tests on public products that need no sky model:
1. Year differences per DA at r4 with the single-year Q,U maps and matrices
   (`maps_da_r4_qu_1yr`, `maps_da_r4_quninv_1yr`): the nine-year repeat of Jarosik et al. 2007
   Eq. (46), Q and U only. A 2 percent rms miscalibration would show as about 0.96 in every pair.
2. DA differences within a band and year (Q1 minus Q2, V1 minus V2).
3. Band-file convention: compare sigma_0,band^-2 times the band matrix with the sum over years and
   DAs of sigma_0,DA^-2 times the single-year matrices.
4. Apply the released imbalance templates v to the band matrices and report v^T N^-1 v against
   v^T v times the mean diagonal, to settle projected versus 1-sigma.
5. When cutting a mask, marginalize (invert, cut, invert back), as in Jarosik et al. 2011 Sect.
   4.1.2, and state it.

Literature still to read: Larson et al. 2011 Sect. 5; Planck 2013 XV App. E from the journal PDF;
the seven-year Explanatory Supplement (wording on the imbalance term, to date the Supplement
sentence); Gruppuso et al. 2011 and later BolPol papers for per-band noise checks; the FITS header
comments of `wmap_band_quninv_r4_9yr_*_v5.fits`.

## 7. What is established, what is not

Established (read levels as given above):
- No WMAP document gives the accuracy of the r4 polarization N^-1 normalization. Since the
  seven-year release the matrices are a unit-normalized r4 filter product divided by a sigma_0
  measured on r9 maps (Jarosik et al. 2011 Sect. 2.7; Bennett et al. 2013 Sect. 4.1.3).
- The matrices are stated approximations of the map noise: per-year stationary noise fit, grouped
  observations, maps degraded from r9 (Jarosik et al. 2007 Sects. 3.4.2 to 3.6).
- The only r4 calibration (three-year) found r4 noise 1.5 to 5.7 percent above the r9 sigma_0(I).
- Team pixel-space chi-squares are at or above one (1.000 to 1.022; nulls 1.009 +- 0.040).
- Hinshaw et al. 2009: "slight tendency to overestimate the uncertainties at higher l".
- Published V-band values with the nine-year raw matrices: 0.961 to 0.977 (Lattanzi et al. 2017),
  0.970 (Planck 2018 V); WMAP minus LFI 70 half-difference 0.973 (Planck 2015 XI). None discussed.
- BeyondPlanck XV saw a slightly negative high-latitude chi-square offset and blamed n_dof.

Not established:
- That the nine-year matrices overstate the noise by 2 percent in rms. No source says so; the
  published hints are 1 to 2 sd each and share the same data.
- The nine-year imbalance treatment (Supplement and LAMBDA page disagree) and the band sigma_0
  convention.
