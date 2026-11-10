# Citation verification for Paper A1, round 2: 7 items checked against the sources themselves

Date: 2026-11-10. Author: lit-scout agent. Round 1 is
`nextgen/research/2026-11-10_citation_verification_paperA1.md`; format and conventions are the same.

Question: for each of the 7 items below, is the statement supported by the cited source, where
exactly, in what words, and what is the correct bibliographic record?

Read levels:
- PDF: I read the stated pages of the PDF as page images and transcribed the quotes myself.
- HTML-S: full-text HTML (ar5iv, or an official documentation page) read through the fetch tool's
  small model, which returns quotes of at most about 125 characters. Long sentences are therefore
  given as consecutive pieces, each in its own quotation marks. Where a passage was queried more
  than once with different prompts, this is stated.
- ABS: abstract only. BIB: bibliographic record only (Crossref, arXiv API, JSTOR cover page).

Conventions: Greek letters and mathematics are written in ASCII. Typographic dashes inside quotes
are written as hyphens. "(sic)" marks a misprint that is in the source. "(inference)" marks my own
arithmetic or synthesis, never a result of the cited source. Section numbers are those of the
version read.

Access limits that shaped the read levels: arXiv PDFs above 10 MB could not be fetched (2007.04997
and 2303.08095). The A&A site, the Cardiff ORCA copy and the Padova copy of the NPIPE paper returned
HTTP 403. The ar5iv rendering of 2007.04997 ends inside Sect. 7.4, so Sects. 8 and 9 and all
appendices of the NPIPE paper were not read at any level. Project Euclid pages for Huber (1967)
returned empty content; the paper was read from the UC Berkeley library scan instead.

## 0. Summary table

| # | Source | Read | Verdict | Correction or caveat |
|---|---|---|---|---|
| 1a | Huber 1967 | PDF, pp. 221, 226-233 of 221-233 | CONFIRMED | Corollary to Theorem 3 (p. 231): covariance Lambda^-1 C (Lambda')^-1. i.i.d. data, fixed dimension. |
| 1b | White 1982 | PDF, cover, pp. 1, 4-7 of 1-25 | CONFIRMED | Theorem 3.2 (p. 6); C = A^-1 B A^-1 defined, unnumbered, on p. 5. Pages 1-25 (JSTOR). A corrigendum exists (1983), not read. |
| 2a | Bennett et al. 2013, 1212.5225 | PDF, pp. 1-26, 133-138 of 177; HTML-S keyword search | CONFIRMED WITH CORRECTION | The r4 Q/U N^-1 products are not described in this paper. Construction is in Jarosik et al. 2011 (Sects. 2.6-2.11) and 2007 (Sect. 3.5); nine-year specifics only on the LAMBDA page. |
| 2b | same | PDF, App. A, pp. 136-138 | CONFIRMED | App. A: reduce pre-flight centre frequencies by 0.13, 0.12, 0.11, 0.06 percent for K, Ka, Q, V; no variation detected in W. |
| 2c | same; Page et al. 2007 | HTML-S (Bennett); PDF pp. 6-8 (Page) | CONFIRMED WITH CORRECTION | Not stated in Bennett 2013. Page et al. 2007 Sect. 2: HEALPix convention, not the IAU position angle. |
| 3 | Jarosik et al. 2003 | PDF, pp. 13-24 of 24 | CONFIRMED WITH CORRECTION | The record in the task mixes two papers. Table 11 is in "Design, Implementation, and Testing of the MAP Radiometers", ApJS 145, 413. Equations are (48) and (49); there is no Eq. (50). |
| 4a | NPIPE, 2007.04997 | HTML-S, to Sect. 7.4 | CONFIRMED WITH CORRECTION | The paper's term is "half-ring" (HR1, HR2). 30 and 44 GHz A/B maps are split by year, not by ring half. White-noise covariance is documented on distribution pages, not found in the paper text read. |
| 4b | same | HTML-S | CONFIRMED WITH CORRECTION | Only 30 and 353 GHz are calibrated without the polarization prior. 217 GHz is calibrated with 30 and 353 GHz templates. No EE suppression at 30 or 353 GHz. |
| 4c | same | HTML-S | NOT FOUND (explicit) | No sentence says which bolometers enter the 353 GHz maps. Tables 4 and 5 imply the four SWBs are included (inference). PR3 was PSB only. |
| 4d | same; Planck 2018 II | HTML-S; PDF | CONFIRMED WITH CORRECTION | NPIPE: 167 ms baselines everywhere, no knee frequencies in the text read. Knee frequencies 113.9, 53.0, 19.6 mHz are from Planck 2018 II Table 1. |
| 5 | Planck 2018 XI, 1801.04945 | PDF, pp. 1-3, 8-15 | CONFIRMED | T_d fixed at 19.6 K, beta_d fitted. High-latitude mean from cross-spectra, PR3, PSB-only 353 GHz. MBB holds "to below 70 GHz". |
| 6 | Watts et al. 2023, 2303.08095 | HTML-S, arXiv v1, 0 to 300 000 of 310 693 chars | CONFIRMED | Prior N(-3.30, 0.1^2) for the high-latitude index. K to 30 scaling 0.495 for index -3.1. No numerical nu_0,s in the text; 30 GHz from a figure caption. |
| 7 | Planck 2018 II, 1807.06206 | PDF, pp. 1-4, 11, 19-20, 22-23, 26-27, 30-31; HTML-S | CONFIRMED WITH CORRECTION | Public N_side 16 NCVMs exist (`LFI_NoiseCovMat_0??_0016_R3.??.tgz`) but are documented in the PLA Explanatory Supplement, not in the paper. Large-scale residuals at 30 and 44 GHz: Sect. 8.6 and App. A. |

## 1. Sandwich covariance: Huber 1967 and White 1982

### 1a. Huber 1967

Identifier: Peter J. Huber (Swiss Federal Institute of Technology), "The behavior of maximum
likelihood estimates under nonstandard conditions", Proceedings of the Fifth Berkeley Symposium on
Mathematical Statistics and Probability, Vol. 1, pp. 221-233. Title, author and page numbers are
printed on the pages read (first page 221, references end on 233; the next article, by Kempthorne,
starts on 235). Year 1967, editors L. M. Le Cam and J. Neyman, University of California Press: this
is the customary record and is NOT confirmed from a record I could load (Project Euclid page
`euclid.bsmsp/1200512988` returned empty content; a Crossref bibliographic query returned no
record, so no DOI is given).

Read: PDF (scan, `digitalassets.lib.berkeley.edu/math/ucb/text/math_s5_v1_article-13.pdf`),
pp. 221 and 226 to 233. Pages 222 to 225 (Sects. 2 and 3, consistency) not read.

Where: Sect. 4 "Asymptotic normality", pp. 226 to 232. Theorem 3 on p. 230, its Corollary on p. 231.

Quotes:
- p. 221: "This paper proves consistency and asymptotic normality of maximum likelihood (ML)
  estimators under weaker conditions than usual. In particular, (i) it is not assumed that the true
  distribution underlying the observations belongs to the parametric family defining the ML
  estimator, and (ii) the regularity conditions do not involve the second and higher derivatives of
  the likelihood function."
- p. 226: "Assume that x_1, x_2, ... are independent identically distributed random variables with
  values in X and common distribution P. We want to give sufficient conditions ensuring that every
  sequence T_n = T_n(x_1, ..., x_n) satisfying (1/sqrt(n)) sum_{i=1}^n psi(x_i, T_n) -> 0 in
  probability (27) is asymptotically normal."
- p. 227, Eq. (28): lambda(theta) = E psi(x, theta). "Expectations are always taken with respect to
  the true underlying distribution P."
- p. 230: "THEOREM 3. Assume that (N-1) to (N-4) hold and that T_n satisfies (27). If
  P(|T_n - theta_0| <= d_0) -> 1, then (1/sqrt(n)) sum_{i=1}^n psi(x_i, theta_0) +
  sqrt(n) lambda(T_n) -> 0 (56) in probability."
- p. 231, the decisive sentence: "COROLLARY. Under the conditions of theorem 3, assume that lambda
  has a nonsingular derivative Lambda at theta_0 (that is, |lambda(theta) - lambda(theta_0) -
  Lambda (theta - theta_0)| = o(|theta - theta_0|)). Then sqrt(n)(T_n - theta_0) is asymptotically
  normal with mean 0 and covariance matrix Lambda^-1 C (Lambda')^-1, (where C stands for the
  covariance matrix of psi(x, theta_0), and Lambda' is the transpose of Lambda)."
- p. 231, the correctly specified case: "We want to show that lambda(theta_0) = 0 and Lambda = -C =
  -I(theta_0), where I(theta_0) is the information matrix. This implies, by the corollary to
  theorem 3, that the asymptotic variance of the ML estimator is I(theta_0)^-1."

Verdict: CONFIRMED. The form is in the Corollary to Theorem 3, p. 231, in an unnumbered display.
In the draft's notation A is Lambda (the derivative of the expected estimating function) and B is C
(the covariance of the estimating function). The estimator is any root of an estimating equation
(27), not only a maximum-likelihood estimator, which is what the draft says.

### 1b. White 1982

Identifier: Halbert White, "Maximum Likelihood Estimation of Misspecified Models", Econometrica,
Vol. 50, No. 1 (Jan., 1982), pp. 1-25 (JSTOR cover page of the copy read; journal header "Volume
50, January, 1982, Number 1" on p. 1). DOI 10.2307/1912526 (Crossref gives the first page only).
The Econometric Society web page lists pages 1-26 (HTML-S); the JSTOR cover and RePEc say 1-25.

Read: PDF (scanned JSTOR copy on a university course page,
`faculty.utrgv.edu/diego.escobari/teaching/Econ8370/Papers/White(1982)Econometrica-Corey.pdf`),
cover page and journal pp. 1 and 4 to 7. Pages 2 to 3 (Assumptions A1 to A3) and 8 to 25 not read.

Where: Sect. 3 "Asymptotic normality", pp. 5 to 7.

Quotes:
- p. 1 (abstract): "The quasi-maximum likelihood estimator (QMLE) converges to a well defined limit,
  and may or may not be consistent for particular parameters of interest. Standard tests (Wald,
  Lagrange Multiplier, or Likelihood Ratio) are invalid in the presence of misspecification, but
  more general statistics are given which allow inferences to be drawn robustly."
- p. 1: "The consistency question was apparently first considered independently by Berk [7, 8] and
  Huber [20]."
- p. 5, definitions (unnumbered displays): A_n(theta) = {n^-1 sum_t d^2 log f(U_t, theta) /
  d theta_i d theta_j}; B_n(theta) = {n^-1 sum_t d log f(U_t, theta)/d theta_i . d log f(U_t,
  theta)/d theta_j}; A(theta) = {E(d^2 log f(U_t, theta)/d theta_i d theta_j)}; B(theta) =
  {E(d log f(U_t, theta)/d theta_i . d log f(U_t, theta)/d theta_j)}. "When the appropriate inverses
  exist, define C_n(theta) = A_n(theta)^-1 B_n(theta) A_n(theta)^-1, C(theta) = A(theta)^-1
  B(theta) A(theta)^-1."
- p. 6, the decisive statement: "THEOREM 3.2 (Asymptotic Normality): Given Assumptions A1-A6
  sqrt(n)(theta^_n - theta_*) ~A N(0, C(theta_*)). Moreover, C_n(theta^_n) -> C(theta_*) a.s.,
  element by element."
- p. 7: "THEOREM 3.3 (Information Matrix Equivalence): Given Assumptions A1-A7, if g(u) = f(u,
  theta_0) for theta_0 in Theta, then theta_* = theta_0 and A(theta_0) = -B(theta_0), so that
  C(theta_0) = -A(theta_0)^-1 = B(theta_0)^-1, where -A(theta_0) is Fisher's information matrix."
- p. 4, footnote 3: "Strictly speaking, the regularity conditions given here allow our results to be
  applied only to regression equations with i.i.d. regressors (e.g., as in White [33, 36])."

Verdict: CONFIRMED. Theorem 3.2, p. 6, with C = A^-1 B A^-1 defined on p. 5. The only numbered
equation nearby is (3.1), the information-matrix equality.

Notes for the draft (both papers):
- Neither paper uses the word "sandwich" on the pages read. Writing "the covariance now called a
  sandwich" is accurate; attributing the name to them is not.
- A corrigendum to White (1982) is listed by the Econometric Society (Econometrica 51, No. 2,
  March 1983, reported as p. 513 by a search result). Seen only as a title and a search snippet;
  page number and content not verified.
- (inference) Both theorems are stated for independent identically distributed observations and a
  parameter vector of fixed dimension. The profile likelihood of Paper A1 has one amplitude per
  element, so the number of nuisance parameters grows with the data and the elements are not
  identically distributed. The two citations support the form A^-1 B A^-1; they do not by
  themselves establish that the sandwich is the right covariance of the profile estimator. That
  rests on the draft's own derivation (App. on score and curvature) and on E041, which is how
  `sec_profile.tex` already presents it.

Bibliographic corrections for `refs.tex`: Huber, add the full page range 221-233 (and editors and
publisher once confirmed from a library record). White, pages 1-25.

## 2. Bennett et al. 2013, nine-year WMAP final maps and results

Identifier: arXiv 1212.5225v3 (stamped 4 Jun 2013). ApJS 208, 20 (2013), DOI
10.1088/0067-0049/208/2/20 (Crossref: volume 208, issue 2, article 20, online 2013-09-20, issue
date 2013-10-01). 21 authors, first five C. L. Bennett, D. Larson, J. L. Weiland, N. Jarosik,
G. Hinshaw (p. 1 of the PDF and Crossref). The arXiv API record has no journal reference. The
bibitem in `refs.tex` is correct and lacks only the arXiv number.

Read: PDF of arXiv v3, pp. 1 to 26 and 133 to 138 of 177 (pp. 139 to 165 are the point-source
tables and App. D). HTML-S keyword search of the ar5iv rendering, characters 0 to 300 000 of
314 674, for the polarization convention and for any description of the r4 N^-1 matrices.

Earlier WMAP papers read for item (a):
- Jarosik et al. 2011, arXiv 1001.4744v1, "Seven-Year Wilkinson Microwave Anisotropy Probe (WMAP)
  Observations: Sky Maps, Systematic Errors, and Basic Results", DOI 10.1088/0067-0049/192/2/14
  (ApJS 192, 14, read off the DOI; Crossref not queried). PDF, pp. 3 to 15 of 42.
- Jarosik et al. 2007, arXiv astro-ph/0603452v2, "Three-Year WMAP Observations: Beam Profiles, Data
  Processing, Radiometer Characterization and Systematic Error Limits", ApJS 170, 263 (2007), DOI
  10.1086/513697 (arXiv API). PDF, pp. 22 to 41 of 58.
- Page et al. 2007, arXiv astro-ph/0603450v2, "Three Year WMAP Observations: Polarization
  Analysis", ApJS 170, 335 (2007), DOI 10.1086/513699 (arXiv API). PDF, pp. 6 to 8 of 105; HTML-S
  for characters 0 to 200 000 of 232 537.
- LAMBDA product pages `product/wmap/dr5/ninv_info.html` and `m_products.html`, HTML-S, the first
  queried twice.

### 2a. Low-resolution Q/U maps and dense N^-1 matrices

What Bennett et al. 2013 itself says (PDF):
- Sect. 2.2.1 (p. 8): "The standard WMAP map-making procedure is unchanged from the previous release
  and the resulting maps are used for the core cosmological analyses. Progress has been made on the
  algorithm for estimating the noise properties of the maps."
- Sect. 2.2.1 (p. 8): "Another difference from previous analyses is that this procedure now
  determines the noise in the polarized maps from the Stokes Q and U year-to-year differences while
  including a spurious ("S") map term, and a mean monopole is subtracted from each S map, as is done
  separately for Stokes I in the temperature map analysis."
- Footnote 2 (p. 8): "The map resolution levels refer to the HEALPix pixelization scheme (Gorski et
  al. 2005) where r4, r5, r9, and r10 refer to N_side values of 16, 32, 512, and 1024,
  respectively."
- Sect. 2.2.1 (p. 8): "In the current processing we have changed to masking based on the brightness
  in each individual band." Table 4 (p. 19), masked pixels of 3072: K (year not 2) 312, K (year 2)
  270, Ka 212, Q 201, V 125, W 98.
- Captions of Figs. 4 and 5 (pp. 21 to 22): "The smooth large angular scale features visible in
  W-band, and to a lesser extent in V-band, are the result of a pair of modes that are poorly
  constrained in map-making, yet properly de-weighted in the analysis."
- Sect. 4.1.2 (pp. 16, 19), Eqs. (8) to (10): t~_0 = M_am^T N^-1 d, "N^-1 is the inverse of the
  radiometer noise covariance matrix", t~ = (M_am^T N^-1 M)^-1 t~_0.
- Sect. 4.1.3 is titled "Noise Characterization of the High Resolution Maps" and treats r9 and r10
  only. No section on the pages read describes the r4 Q/U maps or their matrices. The HTML-S search
  found no sentence that ties "loss imbalance" or "transmission imbalance" to a covariance matrix
  in characters 100 000 to 300 000; Sect. 7.1 concerns the res 3 polarization likelihood.

What the matrices contain, from the earlier papers (PDF):
- Time-correlated radiometer noise. Jarosik et al. 2007, Sect. 3.4.2, Eq. (23), defines the inverse
  noise correlation function N_tt^-1(Delta t) from "the parametrized noise correlation function",
  and: "The radiometer noise properties are treated as stationary for each year of data, making the
  inverse noise matrix, N^-1, circulant." Sect. 3.5, Eq. (28): Sigma^-1 = M^T N^-1 M =
  sum_{t1,t2} M(t1, p1) N_tt^-1(t1 - t2) M(t2, p2), where "p1 and p2 are pixel indices spanning the
  I, Q, U and S maps". "These matrices were evaluated at r4 on a year-by year basis for all 10 DAs".
  "The result is a 4N_p x 4N_p matrix, Sigma^-1, which may be inverted by standard numerical
  techniques to form the noise matrix Sigma."
- The spurious map S. Jarosik et al. 2011, Sect. 2.6 (p. 10): "polarization degrees of freedom N_p
  is four times the number of map pixels, corresponding to the Stokes I, Q and U components and a
  spurious component, S, used to absorb effects arising from bandpass differences between the two
  radiometers comprising each differencing assembly (DA)."
- Seven-year form with asymmetric masking. Jarosik et al. 2011, Sect. 2.7, Eq. (14): Sigma^-1 =
  (M^T N^-1 M_am) (M_am^T N^-1 M_am)^-1 (M_am^T N^-1 M). "In practice the terms M_am^T N^-1 M and
  M_am^T N^-1 M_am are evaluated at r4 with the N^-1 normalized to unity at zero lag". "Regions of
  the Sigma^-1 matrix corresponding to combinations of Stokes I, Q and U are then converted to noise
  units by dividing by the appropriate combinations of sigma_0(I) and sigma_0(Q, U) given in
  Table 1."
- Transmission (loss) imbalance modes, three-year. Jarosik et al. 2007, Sect. 3.5.1, Eqs. (29) and
  (30): Sigma~^-1 = Sigma^-1 - (Sigma^-1 v (x) Sigma^-1 v)/(v^T Sigma^-1 v), "so that Sigma~^-1 v =
  0", where "v is the map mode to be removed, obtained from differences between the simulated maps
  made with the correct transmission imbalance factors and the incorrect factors." And: "Inverse
  noise covariance matrices based on these correlation matrices are contained in the data release.
  The covariance matrices are the correlation matrices scaled by a multiplicative factor,
  1/sigma_0^2, where sigma_0 is the noise per observation corresponding to each DA."
- Transmission imbalance modes, seven-year. Jarosik et al. 2011, Sect. 2.9 (pp. 12 to 13): "To
  prevent biasing the cosmological analyses, sky map modes that can be excited by these measurement
  errors are identified and projected from the Sigma^-1 matrices." "(1) The simulated sky maps are
  formed at r4 by direct evaluation of equations (5, 6) at r4. This change was adopted to eliminate
  contamination of the transmission imbalance templates by the poorly measured sky map modes
  associated with monopoles in the Stokes I and spurious mode S sky maps." "(2) Transmission
  imbalance modes were evaluated both for the case when the x_im for both radiometers comprising
  each DA were increased 20% above their measured values, and for the case when the x_im value for
  radiometer 1 was increased by 10% and that of radiometer 2 decreased by 10%." "For six of the DAs
  very similar sky map modes were generated by both sets of x_im, while for the V2, W1, W2, and W4
  somewhat different modes were generated for the two different sets. As a result, only one mode was
  projected from the K1, Ka1, Q1, Q2, and W3 matrices while two modes were projected out of the the
  (sic) V2, W1, W2 and W4 matrices." (The sentence says six DAs and names five; V1 is not named.)
- Band maps and a de-weighted mode in K and Ka. Jarosik et al. 2011, Sect. 2.10: "The low
  resolution sky maps were generated by performing an inverse noise weighted degradation of the high
  resolution maps that takes into account the intra-pixel noise correlations of the high resolution
  input maps." Sect. 2.11, Eqs. (18) and (19): (Q, U) = Sigma~^tot sum_{yr,DA} Sigma~^-1_{yr,DA}
  (Q, U)_{yr,DA}, Sigma~^tot = (sum_{yr,DA} Sigma~^-1_{yr,DA})^-. "For K and Ka-bands the
  pseudo-inverse is calculated by performing a singular value decomposition of the sum of the QU x
  QU Sigma~^-1 matrices and inverting all its eigenvalues except for the smallest, which is set to
  zero. This procedure de-weights the one nearly singular mode associated with the monopoles in the
  I and S maps. No modes are removed from the Q, V and W-band maps."
- Why the I and S monopoles matter. Jarosik et al. 2007, Sect. 3.4.3 (p. 37): "For an ideal
  differential radiometer the mean value of the I and S maps would be undetermined, corresponding to
  singular modes in Sigma. However, inclusion of the non-zero valued transmission imbalance factors,
  x_im, converts these previously singular modes into modes with very small eigenvalues, indicating
  that they are very poorly constrained by the measurement."
- Mask. Jarosik et al. 2007, Sect. 3.5 (three-year): "The r4 noise matrices contained in the data
  release are calculated from mapping matrices corresponding to cut sky coverage and are intended
  for use in regions outside of the region excluded by the processing mask. Values corresponding to
  r4 pixels entirely contained within the processing mask are set to zero." Jarosik et al. 2011,
  Sect. 2.6 (p. 11, seven-year): "Utilizing asymmetric masking eliminates both of the aforementioned
  shortcomings of the simple conjugate gradient method - each DA only requires a single map
  solution, and a pixel-pixel inverse noise covariance matrix can be generated which describes the
  noise correlation between all the pixels in the map."

Nine-year specifics, LAMBDA `ninv_info.html` for data release 5 (HTML-S, pieces):
- "The inverse covariance (N-1) matrices supplied with the higher resolution maps do not fullly
  (sic) characterize" "the pixel-to-pixel noise in the res 4 maps." "Therefore a more accurate and
  much larger inverse covariance matrix is supplied for each map as a separate product."
- "The matrices include a contribution due to uncertainty in the loss imbalance coefficients in
  each differencing assembly" "(Jarosik et al. 2010)." "In previous releases, this uncertainty was
  fully marginalized--equivalent to assuming no information" "about these coefficients." "In the
  latest release of the matrices for raw Stokes Q and U maps (those without foreground cleaning),"
  "a true 1-sigma uncertainty was incorporated into this term."
- "The elements corresponding to the Galatic (sic) plane pixels as defined by the processing mask
  were set to zero." "Therefore, the matrices as supplied are singular."
- "All the inverse covariance matrices have had sigma0^2 divided out of them."
- Ordering and coordinates: "HEALPix, nested, Res 4/NSide=16", "Galactic".
- `m_products.html`: "Polarization (Q,U) coadded maps for each frequency at HEALPixe (sic) Nside=16
  (Res 4) degraded from Nside=512 (Res 9)." (5 files) and "Polarization (Q,U) pixel-to-pixel inverse
  covariance matrices at HEALPixe (sic) Nside=16 (Res 4)." (5 files, 144 MB each). Loss imbalance
  templates, 20 files: "Polarization templates (Q,U) that trace the effect of a small transmission
  imbalance between the A & B sides", "Two per differencing assembly."

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: the nine-year release (Bennett et al. 2013; LAMBDA data release 5) provides
N_side 16 Stokes Q and U maps per band with dense 6144 x 6144 inverse noise covariance matrices.
Bennett et al. 2013 does not describe these products. Their construction is given by Jarosik et al.
2007 (Sect. 3.5, Eqs. 23 and 28 to 30) and Jarosik et al. 2011 (Sects. 2.6 to 2.11). The matrices
contain the time-correlated radiometer noise propagated through the differential scan, are built
in the joint space of I, Q, U and the spurious map S, and carry the uncertainty of the transmission
imbalance: projected out completely in the three-, five- and seven-year matrices (one or two modes
per DA in the seven-year set), and, per the LAMBDA page, entered with "a true 1-sigma uncertainty"
in the nine-year matrices of the raw Q and U maps. sigma_0^2 is divided out and must be restored.

Three points are not settled by what I read:
1. Mask. The LAMBDA page says processing-mask pixels are zeroed and the matrices singular. That
   wording matches the three-year procedure (Jarosik et al. 2007). Jarosik et al. 2011 says the
   seven-year matrices describe all pixels. The draft (`sec_data.tex`, E038) reports the released
   nine-year band matrices as positive definite with no null space. (inference) The LAMBDA sentence
   is probably carried over from earlier releases for the band Q/U matrices; the E038 check of the
   files is the stronger evidence. State in the paper that this was checked on the files.
2. Whether the released QU matrix is the QU block of the joint inverse covariance (conditional on I
   and S) or the inverse of the QU block of the covariance (marginal over I and S). Jarosik et al.
   2011 Sect. 2.11 speaks of "the QU x QU Sigma~^-1 matrices" without saying which. Not found in
   Page et al. 2007 App. D either (HTML-S).
3. The nine-year change from projection to a 1-sigma term is documented only on the LAMBDA page; I
   found no description in Bennett et al. 2013 (pages and characters read as stated) and did not
   read the nine-year Explanatory Supplement (Greason et al. 2012).

The draft sentence "The matrices carry the correlated noise of the differential scan and the poorly
constrained large-scale modes" is supported, but by Jarosik et al. 2007 and 2011 and the LAMBDA
page, not by Bennett et al. 2013 alone. Cite them together.

### 2b. Drift of the effective band centre frequencies

Where: App. A "Band Center Frequencies", pp. 136 to 137; Fig. 44, p. 138; Sect. 3, p. 12; Table 3,
p. 18.

Quotes (PDF):
- "Figure 44 shows small year-to-year variations of Galactic plane brightness measured from yearly
  maps in K-, Ka-, Q-, and V-bands. Each yearly map was correlated against the nine-year map for
  pixels at |b| < 10 deg. A linear slope and offset was fit to each correlation, and the slope
  values are shown in Figure 44. Results for W-band are not shown because the scatter in the yearly
  slopes is large and no significant variation was detected."
- "The K-Q band brightness variations were previously presented for the seven-year data in Jarosik
  et al. (2011), where they were described as variations in the WMAP calibration. Further analysis
  has shown that the CMB signal in yearly maps does not show such variation." "We have also found
  that the Galactic plane brightness variations depend on spectral index, with greater variation for
  regions of steeper spectral index, so we conclude that they are caused by variations in the
  effective center frequencies of the WMAP bandpasses over the mission."
- Eqs. (A1) to (A3): T_i(p) = T_0(p) (nu_i/nu_0)^beta_p, "with thermodynamic temperature spectral
  index beta_p"; Delta nu_i/<nu_i> = (T_i(p)/<T_i(p)> - 1)/beta_p.
- The decisive paragraph: "No correction for bandpass drift is applied in our map-making. Since the
  WMAP observations are made simultaneously in the different bands, the map-making always forms band
  maps that have a common epoch, and each band map can be treated as having a single effective band
  center frequency valid for that epoch. Our previously published band center frequencies (see Table
  4 of Jarosik et al. (2011) for point sources and Table 11 of Jarosik et al. (2003a) for diffuse
  emission) are based on pre-flight measurements, so presumably are valid for year 1 of the flight
  data. For nine-year data, a correction based on Figure 44 should be applied. The correction is a
  reduction of the pre-flight center frequency by 0.13, 0.12, 0.11, and 0.06% for K-, Ka-, Q-, and
  V-band, respectively. This correction is included in the center frequencies for point sources
  listed in Table 3."
- Sect. 3 (p. 12): the point-source quantities of Table 3 "were determined as described in Jarosik
  et al. (2011), except a small correction for bandpass drift was included in the calculation of
  effective frequency for K-, Ka-, Q-, and V-bands as described in Appendix A."
- Table 3, note e: "The effective center frequency for a point source with flux spectral index
  alpha = -0.1. The estimated uncertainty, due to uncertainties in the pre-flight passband response
  measurements, is 0.1% for all DAs."

Verdict: CONFIRMED. The numbers 0.0013, 0.0012, 0.0011, 0.0006 and no correction for W agree with
App. A.

(inference) The paper says "a reduction ... by 0.13 ... %", that is a factor 1 - 0.0013; the IDL
routine divides by 1 + 0.0013. The two differ by 1.7e-6 in relative terms. In the lower panel of
Fig. 44 the year-1 points sit at about +0.0013 (K), +0.0012 (Ka), +0.0011 (Q) and +0.0006 (V)
relative to the mission mean, so the correction moves the pre-flight (year 1) value to the
nine-year mean. The correction is the same size as the 0.1 percent uncertainty the paper assigns to
the pre-flight passbands (Table 3, note e); the draft should say so when it quotes effective
frequencies to four digits. `effective_freq.pro` itself was not fetched.

### 2c. Polarization convention

Bennett et al. 2013: no statement found. The HTML-S keyword search (characters 0 to 300 000) found
no occurrence of "IAU" and no use of "convention" for the Stokes parameters.

Page et al. 2007, Sect. 2, p. 6 (PDF): "gamma is the angle between the polarization direction of
the electric field and the Galactic meridian (Kogut et al. 2003)." "From these, we form maps of
polarization intensity, P = sqrt(Q^2 + U^2), and direction, gamma = (1/2) tan^-1(U/Q). This
convention has gamma positive for North through West and follows the convention in Zaldarriaga &
Seljak (1997) and HEALPix (Gorski et al. 1998). However, it differs from the standard astronomical
position angle (PA) which has gamma_PA = (1/2) tan^-1(-U/Q) with gamma_PA positive for North
through East. The choice of convention does not affect the plots."

LAMBDA `product/about/pol_convention.html` (HTML-S): "In the polarization convention used
historically in the CMB community (COSMO (HEALPix))," "the local x-axis points South, the local
y-axis points East," ... "As a result, there is a difference in the sign of the U Stokes
parameter." "Most HEALPix files will follow the COSMO convention."

Verdict: CONFIRMED WITH CORRECTION. WMAP maps follow the HEALPix (COSMO) convention, not the IAU
one; the source is Page et al. 2007, Sect. 2. Planck states the same for its own maps (item 7).
Not checked: the `POLCCONV` keyword in the nine-year r4 FITS headers; look at the local files.

## 3. Jarosik et al. 2003: Table 11 and the effective frequency

Bibliographic correction: the record in the task joins the title of one paper to the journal
reference of another.
- ApJS 145, 413 is N. Jarosik, C. L. Bennett, M. Halpern, G. Hinshaw, A. Kogut, M. Limon,
  S. S. Meyer, L. Page, M. Pospieszalski, D. N. Spergel, G. S. Tucker, D. T. Wilkinson, E. Wollack,
  E. L. Wright, Z. Zhang, "Design, Implementation, and Testing of the Microwave Anisotropy Probe
  Radiometers", ApJS 145, 413-436 (April 2003), DOI 10.1086/346080 (Crossref), arXiv
  astro-ph/0301164 (arXiv API; arXiv title "Design, Implementation and Testing of the MAP
  Radiometers").
- "First-Year Wilkinson Microwave Anisotropy Probe (WMAP) Observations: On-Orbit Radiometer
  Characterization" is Jarosik, Barnes, Bennett, et al., ApJS 148, 29-37 (September 2003), DOI
  10.1086/377221 (Crossref). Not read; its arXiv number was not checked.
- Table 11 is in the first of these. Bennett et al. 2013 App. A cites "Table 11 of Jarosik et al.
  (2003a) for diffuse emission", and its introduction (p. 3) identifies 2003a as the paper on "the
  radiometer design and characterization" and 2003b as the one on "on-orbit radiometer
  characteristics".

Read: PDF of arXiv astro-ph/0301164v2, pp. 13 to 24 of 24. Pages 1 to 12 not read.

Where: Sect. 7.4.1 "Bandpass Measurement", pp. 20 to 21, Eqs. (36) to (49); Table 11, p. 23.

Quotes:
- "Measurements of the radiometers' bandpasses are necessary for the analysis of in-flight data.
  Unlike the radiometric gains and offset, which are calibrated in flight, bandpass measurements
  must be made before launch."
- "The bandpass was measured by stepping the synthesized source through 201 frequencies, nu_i, and
  recording the output of the radiometer at each frequency with the RF output of the synthesizer
  enabled and disabled."
- "The procedure resulted in a bandpass measurement relating the response from each detector to the
  input microwave power, referenced to the input port of the cold hybrid tee, so it does not include
  possible frequency dependences of the feed horn, OMT and waveguides connecting the OMT to the cold
  hybrid tees. The transmission coefficients of the OMT and waveguide sections were measured with
  the VNA and small (< 0.3 dB) corrections applied to the measured responses, resulting in corrected
  responses, r_l(nu_i) and r_r(nu_i), for each detector."
- Combination of the two detectors: "The calibrated data of the two detectors are combined with a
  simple average to yield the calibrated radiometer output," Eq. (44): T^_th,avg = (1/2)(r^_l y_l +
  r^_r y_r). "The frequency response of this combined signal is an average of the two detector
  frequency responses", Eq. (45): r_avg(nu_i) = (1/2) [r_l(nu_i)/sum_i r_l(nu_i) w'(nu_i) +
  r_r(nu_i)/sum_i r_r(nu_i) w'(nu_i)], with w' from Eq. (40), w'(nu) = x^2 e^x/(e^x - 1)^2 at
  T_CMB = 2.725 K.
- Eq. (46): d T^_th,avg/d T_RJ = sum_i r_avg(nu_i) / sum_i r_avg(nu_i) w'(nu_i). "The thermal center
  frequency for each radiometer, nu_th, is the frequency at which the Rayleigh-Jeans to
  thermodynamic conversion equals that calculated for the radiometer," Eq. (47): d T^_th,avg/d T_RJ
  = 1/w'(nu_th).
- The decisive passage: "The factors need (sic) to convert an observed Rayleigh-Jeans temperature
  difference to a Rayleigh-Jeans temperature difference at another frequency, nu_0, for objects with
  antenna temperature scaling as T_A proportional to nu^beta were calculated from the relation"
  Eq. (48): [dT_RJ/dT_RJ(nu_0)](beta) = sum_i r_avg(nu_i) (nu_i/nu_0)^beta / sum_i r_avg(nu_i). "In
  these units a source with a constant Rayleigh-Jeans temperature has beta = 0. The resulting
  coefficients were fit to a function of the form" Eq. (49): [dT_RJ/dT_RJ(nu_0)](beta) =
  (nu_avg/nu_0)^beta (1 + c_1 beta + c_2 beta^2) "for -3 < beta < 1 and the resultant parameters,
  nu_avg, c_1 and c_2 are given in Table 11."
- "These quantities apply to diffuse sources that completely fill the beams of the instrument. For
  resolved sources additional terms resulting from the frequency dependent beam sizes must be
  included (Page 2003). Tabulated bandpass measurements will be made available as part of the MAP
  data release."
- Table 11, "MAP Radiometer Bandpass Parameters", one row per radiometer (20 rows). Columns:
  nu_avg (GHz), c_1, c_2, a fourth column printed as "nu_eff (GHz)", nu_th (GHz), dT_th/dT_RJ. Rows
  K to V as I read them from the page image (check against the journal table before use):
  K11 22.29, -2.43E-03, 2.78E-03, 5.26, 22.36, 1.013; K12 23.12, -2.35E-03, 2.83E-03, 4.09, 23.18,
  1.014; Ka11 32.78, -1.58E-03, 1.75E-03, 6.75, 32.84, 1.028; Ka12 33.12, -1.85E-03, 2.06E-03, 7.04,
  33.19, 1.029; Q11 40.71, -1.85E-03, 2.10E-03, 8.62, 40.79, 1.044; Q12 40.81, -1.48E-03, 1.70E-03,
  7.66, 40.88, 1.044; Q21 40.16, -1.58E-03, 1.78E-03, 7.55, 40.23, 1.043; Q22 40.92, -1.53E-03,
  1.78E-03, 7.68, 40.99, 1.044; V11 59.20, -2.13E-03, 2.20E-03, 8.37, 59.32, 1.095; V12 61.10,
  -2.00E-03, 2.23E-03, 12.53, 61.22, 1.101; V21 61.61, -1.95E-03, 2.22E-03, 12.81, 61.72, 1.103; V22
  60.63, -2.12E-03, 2.35E-03, 11.66, 60.75, 1.099.

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: Table 11 of Jarosik et al. 2003 (ApJS 145, 413, the radiometer design paper)
gives per radiometer the parameters nu_avg, c_1 and c_2 of Eq. (49), a fit to the bandpass-averaged
factor of Eq. (48) that scales a Rayleigh-Jeans temperature with spectrum nu^beta from nu_0 to the
band. The fit is stated for -3 < beta < 1. The paper does not define an "effective frequency for a
power law" as such, and the arXiv version has no Eq. (50) (the main-text equations end at 49). The
two detectors of a radiometer are combined by a simple average (Eqs. 44 and 45). Nothing on pp. 13
to 24 says how the two radiometers of a differencing assembly, or the assemblies of a band, are
combined.

Flags:
- The fourth column of Table 11 is printed "nu_eff (GHz)" in the arXiv PDF but holds 4 to 21 GHz.
  (inference) It is the effective bandwidth, with a Delta lost in the rendering; the column is not
  an effective frequency. Check the journal table.
- (inference) From Eq. (49), the frequency at which a delta band reproduces the band response to
  nu^beta is nu_eff(beta) = nu_avg (1 + c_1 beta + c_2 beta^2)^(1/beta). With beta = -3.1 and an
  unweighted mean over the radiometers of a band, the Table 11 values give K 22.46, Ka 32.70, Q
  40.36, V 60.10 GHz, and 22.43, 32.66, 40.31, 60.06 GHz after the App. A drift correction. The
  draft quotes 22.44, 32.70, 40.34 and 60.03 GHz from the bandpass profiles. The agreement is within
  0.11 percent, which is the size of the drift and of the stated 0.1 percent pre-flight
  uncertainty. Two caveats: beta = -3.1 is just outside the fitted range of Eq. (49), and the
  unweighted mean is my choice, not the paper's.

## 4. Planck Collaboration Int. LVII 2020 (NPIPE)

Identifier: arXiv 2007.04997v1 (2020-07-09), "Planck intermediate results. LVII. Joint Planck LFI
and HFI data processing". A&A 643, A42 (2020), DOI 10.1051/0004-6361/202038073, online 2020-11-11
(arXiv API and Crossref). The bibitem in `refs.tex` is correct and lacks only the arXiv number.

Read: HTML-S only. ar5iv, characters 0 to 292 899 in three portions, each queried three or four
times with different prompts. The rendering ends inside Sect. 7.4 ("Global parameters"); Sect. 8
(optical depth), Sect. 9 and the appendices are absent and were not read. The published PDF could
not be obtained (see access limits). Supporting official pages, HTML-S: PLA Explanatory Supplement
`NPIPE_Introduction` and `Frequency_maps` (section "2020 - NPIPE", first 100 000 characters), the
CMB-S4 mirror listing `data.cmb-s4.org/planck_pr4-half_ring.html`.

### 4a. Split maps and what a half-difference captures

Quotes (HTML-S, pieces):
- Footnote 5, Sect. 2.1: "Half-ring data sets are built from the first and second halves of the
  pointing periods."
- Sect. 4.1.2, the whole paragraph as returned: "The repetitive Planck scanning strategy (each
  scanning circle is repeated 30-75 times) allows us to build subset maps" "that have effectively
  identical systematics but independent instrument noise." "We split each pointing period into two
  half-rings, and assign the halves to separate subsets, which we call "HR1" and "HR2."" "The
  half-ring difference (i.e., HR1 - HR2) is a useful measure of the instrumental noise," "but should
  not be assumed to be unbiased for large-scale cross-spectral analysis." "The half-ring maps share
  all the calibration, bandpass mismatch, and other template residuals." "The noise in the HFI
  half-ring maps contains a small amount of correlated error between the two ring halves" "(cf.
  Sect. 2.3.4 and Fig. 18)." "This error comes from a noisy signal estimate used in the removal of
  the glitches."
- Sect. 2.3.4, on earlier HFI processing: "This bias was evident in the noise estimates derived
  from half-ring difference maps," "which were about a percent short of the full-data noise due to
  correlated noise cancellation".
- Sect. 4.1.5, the A/B split: "We have split the horns in the focal plane into two independent sets,
  A and B, and performed reprocessing independently on each set." "The detector-set splitting was
  not possible at either 30 or 44 GHz because of the lack of redundant polarized horns." "Instead,
  for these two frequencies, the split is done time-wise: set A comprises operational years 1 and
  3;" "while set B has years 2, 4, and the month of integration time from the final, fifth observing
  year." "The time-wise split is not expected to have independent systematics," "since the two
  subsets will share initial beam and bandpass mismatch;" "NPIPE provides only one data split where
  the systematics between the splits are expected to be uncorrelated." "A half-mission split is not
  feasible because the second half mission does not have full sky coverage."
- Table 5, "Details of the A/B subset split." (as returned): 30 and 44 GHz, years 1 and 3 against
  years 2, 4 and start of 5; 70 GHz, horns 18, 20, 23 against 19, 21, 22; 100 GHz, horns 1, 4
  against 2, 3; 143 GHz, horns 1, 3, 5, 7 against 2, 4, 6; 217 and 353 GHz, horns 1, 3, 5, 7 against
  2, 4, 6, 8.
- Sect. 6.2: "the detector-set difference does not cancel systematics that vary detector by
  detector."
- PLA Explanatory Supplement, "2020 - NPIPE": "Half-ring frequency maps can be used to gauge
  instrumental noise." "The large-scale systematics residuals in this split are fully correlated and
  cancel in the difference." "The half-ring split data were destriped independently with Madam to
  avoid correlated noise residuals." File name patterns: `?FI_SkyMap_fff{-tag}_????_R4.??_full.fits`
  and `?FI_SkyMap_fff{-tag}_????_R4.??_full-ringhalf-{1|2}.fits`. FITS columns listed: I_STOKES,
  Q_STOKES, U_STOKES, HITS, II_COV, IQ_COV, IU_COV, QQ_COV, QU_COV, UU_COV; the NPIPE section does
  not define the covariance columns.
- CMB-S4 mirror, half-ring page: "Half-ring maps, hitmaps and white noise covariance matrices."
  Files per frequency: `npipe6v20_030_map_sub{1,2}of2.fits`, `..._hmap_...`, `..._wcov_...`,
  `..._wcov_inv_...` and `..._wcov_sub{1,2}of2_hrscaled.fits` (same pattern for 044, 070, 353).
- Low-resolution products, Sect. 4.2: "The maps are downgraded to HEALPix resolution Nside = 16 and
  convolved with a cosine apodizing kernel." "The low-resolution maps are accompanied by pixel-pixel
  noise-covariance matrices that reflect the baseline uncertainties" "and the high-frequency
  instrumental noise approximated as uniform white noise." "Since the high-frequency noise in NPIPE
  is not uniform for any of the Planck detectors, the approximation is inexact and," "as a result,
  requires an extra scaling step performed on the diagonal of the matrix to match the half-ring
  difference maps." "The off-diagonal elements correspond to low frequency noise and need not be
  scaled."

Verdict: CONFIRMED WITH CORRECTION.

Corrections:
1. The paper's term is "half-ring" (subsets HR1 and HR2: the first and second half of each
   pointing period). "ring-half" occurs only in the PLA file names (`full-ringhalf-1`). The draft
   should say "half-ring" and may give the file-name form once.
2. The "A/B" maps at 30 and 44 GHz are a split by year, not by ring half and not by horn, and the
   paper says that split "is not expected to have independent systematics". If Paper A1 ever uses
   A/B differences at 30 or 44 GHz, this applies.
3. What the half-difference captures, in the paper's words: it "is a useful measure of the
   instrumental noise" but "should not be assumed to be unbiased for large-scale cross-spectral
   analysis", and the two halves "share all the calibration, bandpass mismatch, and other template
   residuals". So the draft's sentence that the difference "is blind to errors common to the two
   halves" is supported.
4. Not found in the text read: any statement on whether the half-ring difference over- or
   underestimates the LFI noise of the full map, or on destriping residuals in the half-difference.
   The only quantified statement concerns earlier HFI processing ("about a percent short"). Planck
   2018 II (item 7) says of the PR3 LFI maps that "The actual noise in the LFI data is given
   directly by the half-ring difference maps."
5. "Per-pixel white-noise covariance" is supported by the PLA column list and the CMB-S4 listing,
   not by a sentence of the paper that I could read. The mirror also carries `_hrscaled` variants
   of the half-ring white-noise covariance, and the round-1 landscape note lists
   `..._wcov_mcscaled.fits` for the full maps. (inference) These are covariances rescaled to
   half-ring differences or to Monte Carlo scatter. The draft should name the file it used.

(inference) Two effects push the half-ring difference in opposite directions as an estimate of the
noise of the full map: noise correlated over more than half a pointing period partly cancels in the
difference (too low), and independent destriping of half the data leaves larger baseline errors
than the full solution (too high). Neither is quantified for LFI in what I read.

### 4b. Calibration without the polarization prior; transfer function; 30 GHz

Quotes (HTML-S, pieces):
- Sect. 2.1: "We use a hybrid calibration scheme to work around a degeneracy in the Planck
  calibration process." "The 30- and 353-GHz data are first calibrated along with the polarized sky
  much like PR3 HFI is calibrated." "We then mimic the PR3 LFI calibration and calibrate the
  intermediate CMB frequencies (44-217 GHz)," "using an approximation that the sky polarization can
  be captured with templates derived from the foreground frequencies" "(Sect. 2.4.13)." "The
  calibration for the CMB frequencies is thus based on temperature-only sky and polarized foreground
  templates" "The use of this polarization prior significantly reduces the uncertainty in gain" "and
  the large-scale polarization systematics," "but also requires a measurable and correctable
  transfer function affecting large-scale CMB polarization" "(described in Sect. 4.3)."
- Sect. 2.4.13: "To this end, the 217-GHz reprocessing fits for 30- and 353-GHz polarization
  templates smoothed to 1 deg," "while 44-143 GHz reprocessing fits for 30, 217, and 353 GHz." "For
  each template, we only allow for a single template amplitude across all detectors."
- Sect. 4.3: "The suppression of EE power is limited to frequencies that are calibrated with the
  polarization prior." "There is no suppression at 30 GHz or 353 GHz". "We measure the transfer
  function from our simulations by comparing the input CMB map to the output" "foreground-cleaned
  map." "We find that the EE quadrupole and octupole are significantly suppressed by our
  calibration, with about 64 %" "of the quadrupole and 73 % of the octupole missing in the NPIPE
  power spectra." "Multipoles 4 and 5 are missing about 30 % of the power, and l = 6 is missing
  about 20 %." "In these files, the transfer function is unity everywhere except for the l < 42" "EE
  modes." At 44 GHz the function is closer to unity because "the 30-GHz template shields about 22 %
  of the CMB polarization."
- Sect. 4.3.1: "measurement of these signals is driven by the 30- and 353-GHz channels," "neither of
  which is subject to the polarization prior in the calibration."
- Fig. 21 caption, on transfer functions computed for 30 and 353 GHz: "These two transfer functions
  are merely of diagnostic value" and "are not applied in any analysis."
- Fig. 41 caption: "353 GHz is not shown because it is always calibrated without the polarization
  prior."
- Sect. 6.3.1: "polarized Galactic emission is not affected by the same transfer function at low
  multipoles as the CMB."
- Sect. 6.3.2, comparison with WMAP K band at 5 degrees: "The residuals for 30 GHz are shown in Fig.
  50, and for 44 GHz in Fig. 51." "They also show that the NPIPE large-scale polarization is more
  compatible with WMAP than is that of PR3."

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: in NPIPE the 30 and 353 GHz channels are calibrated with the polarized sky
solved alongside, without the polarization prior. The 217 GHz channel is calibrated WITH the prior,
using 30 and 353 GHz templates, and 44 to 143 GHz with 30, 217 and 353 GHz templates. The EE
suppression is limited to the channels calibrated with the prior; the paper says there is none at
30 or 353 GHz.

(inference) Two consequences for Paper A1, neither stated by the paper in these words. First, the
30 GHz map does not benefit from the reduction of gain uncertainty that the prior brings to 44 to
217 GHz, so its large-scale polarization errors are the larger ones. Second, the 44 and 70 GHz maps
were calibrated assuming that their polarized foreground follows the 30, 217 and 353 GHz templates
with one amplitude per template, which is a constant-index assumption at the calibration stage. A
synchrotron index fitted between PR4 30 GHz and PR4 44 or 70 GHz on large scales is therefore not
independent of that assumption. The paper's remark that Galactic emission "is not affected by the
same transfer function" as the CMB does not quantify this. Worth one sentence in `sec_data.tex`.

### 4c. Which bolometers enter the 353 GHz polarization maps

Not found as an explicit statement, in the three portions and with four differently worded
queries. What the text read does contain:
- Sect. 2.2: "For each polarized HFI horn there are two polarization-sensitive bolometers (PSBs
  labelled "a" and "b")" "at orthogonal angles of polarization sensitivity." "Some of the HFI horns
  instead house nearly unpolarized "spiderweb" bolometers (SWBs)." "With two of the SWBs
  irrecoverably compromised by random telegraphic signal," "we have 52 bolometer timelines that are
  preprocessed into 50 optical timelines."
- Table 4 (Sect. 2.4.15) lists a fitted polarization efficiency and angle for all twelve 353 GHz
  bolometers, including 353-1 (0.23 +- 0.13 percent), 353-2 (7.49 +- 0.05), 353-7 (5.31 +- 0.12) and
  353-8 (9.88 +- 0.08), against 85 to 93 percent for 353-3a to 353-6b. Caption: "For NPIPE,
  statistical uncertainties in the PSB polarization angles are small," "while those for the SWBs are
  up to 2.3 deg."
- Table 5: the 353 GHz A and B sets are horns 1, 3, 5, 7 and 2, 4, 6, 8.
- Sect. 6.3: "In temperature, maps from 100 to 217 GHz have diverged in the Galactic plane, owing
  to the exclusion of all SWBs from SRoll2 processing."
- Secondary source, Rosenberg, Gratton, Efstathiou 2022 (arXiv 2205.10869, Sect. 2.3, HTML-S): "The
  NPIPE detector sets are similar to the detector-set splits explored extensively in EG21 but also
  group" "one or two spiderweb bolometers with each set of polarization sensitive bolometers."
- For PR3: Planck 2018 III (arXiv 1807.06207, Sect. 1, HTML-S): "However, for polarization studies,
  baseline maps at 353 GHz are based on polarization-sensitive bolometer (PSB)" "observations only
  (for reasons explained later in this paper), although maps with spider-web bolometers (SWBs)" "are
  also made available for intensity studies." Planck 2018 XI, Sect. 2 (PDF, item 5) says the same.

Verdict: NOT FOUND as an explicit statement. (inference) The NPIPE 353 GHz maps include the four
SWBs (horns 1, 2, 7, 8) together with the eight PSBs: the A/B sets contain those horns and the
pipeline fits polarization parameters for them. With polarization efficiencies of 0.2 to 10
percent, the SWBs carry little weight in Q and U, so the Q/U bandpass is close to the PSB one, but
it is not identical.

Consequence for the draft: `sec_data.tex` says the dust model is "seen through the band of the
polarization-sensitive bolometers". For PR4 that is an approximation, not what the map is. Either
use the NPIPE instrument model weights or say that the SWB contribution to Q and U is neglected.

### 4d. LFI noise

- NPIPE, Sect. 2.1: "We destripe all data with Madam, using extremely short 167-ms baselines
  (Sect. 3)." Sect. 3: "fitting each detector with 167-ms baseline offsets that correspond to 1 deg
  steps on the sky." Sect. 6.2: "Finally, we fit the 1/f noise fluctuations with 167-ms baseline
  steps, compared to the 0.25-s (30 GHz)," "1-s (44 and 70 GHz), and full-pointing-period (35-75 min
  at HFI frequencies) baseline steps used in PR3."
- NPIPE, Sect. 5.1: "The 1/f instrument noise is simulated using the inverse fast Fourier transform
  (FFT) technique." "all polarized horns show significant internal noise correlations at low
  frequency," "but only the HFI horns have significant noise correlations at high frequency." Fig.
  27 caption: "The PSDs are used for simulating the 1/f noise fluctuations, and as inputs to the
  Madam noise filter for destriping."
- No knee frequency is given in the NPIPE text read (Sects. 1 to 7.4).
- Planck 2018 II (PDF, item 7), Table 1, "LFI performance parameters", for 30, 44, 70 GHz: centre
  frequency 28.4, 44.1, 70.4 GHz; white-noise level in timelines 147.9, 174.0, 151.9 uK_CMB s^1/2;
  f_knee 113.9, 53.0, 19.6 mHz; 1/f slope -0.92, -0.88, -1.20. Note c: "Typical values derived from
  fitting noise spectra (see Sect. 5)." Sect. 5, Eq. (6): P(f) = P_0^2 [1 + (f/f_knee)^beta].
  Sect. 1 (p. 3): "As in the 2015 release, we adopt baselines of 1 s at 44 and 70 GHz, and 0.25 s at
  30 GHz. The shorter baseline at 30 GHz is appropriate for the higher 1/f noise of the radiometers
  at this frequency (see Table 1), which introduces a larger correlated component in the noise."

Verdict: CONFIRMED WITH CORRECTION. NPIPE uses 167 ms baselines at all frequencies; the 0.25 s and
1 s baselines are PR3. Knee frequencies must be cited from Planck 2018 II, Table 1, and are PR3
values.

(inference) The ordering of the knee frequencies, 30 above 44 above 70 GHz, matches the ordering of
the large-scale excess the draft measures in the half-ring differences (1.77, 1.40, 1.21). That is
a consistency remark, not a derivation.

## 5. Planck 2018 results XI, polarized dust foregrounds

Identifier: arXiv 1801.04945v3 (12 Nov 2018). A&A 641, A11 (2020), DOI
10.1051/0004-6361/201832618 (arXiv API). "Received 11 January 2018 / Accepted 14 September 2018"
(p. 1). The bibitem in `refs.tex` is correct and lacks only the arXiv number.

Read: PDF of arXiv v3, pp. 1 to 3 and 8 to 15. Sects. 3 and 6 to 7 and the appendices not read.

Where: abstract; Sect. 2 (p. 3); Sect. 4.1 and 4.2 (pp. 8 to 11), Eq. (2), Table 2; Sect. 5.1 to
5.4 (pp. 12 to 15), Table 3, Figs. 14 to 16.

Quotes:
- Abstract: "The SED of polarized dust emission is fit well by a single-temperature modified
  black-body emission law from 353 GHz to below 70 GHz. For a dust temperature of 19.6 K, the mean
  dust spectral index for dust polarization is beta_d^P = 1.53 +- 0.02. The difference between
  indices for polarization and total intensity is beta_d^P - beta_d^I = 0.05 +- 0.03."
- Abstract, sky: "for a set of six nested high-Galactic-latitude sky regions covering from 24 to
  71 % of the sky."
- Sect. 2: "The 100-, 143-, and 217-GHz HFI maps are made using data from all bolometers, while the
  353-GHz maps are constructed using only data from the polarization-sensitive bolometers (PSBs), as
  recommended in Planck Collaboration III (2018)." "Thus, we consider that there is a 1.5 %
  photometric uncertainty on the 353 GHz polarized emission." "In addition, in Sect. 4 we use
  polarization maps from LFI at 30 GHz, and the K and Ka WMAP channels (Bennett et al. 2013) to
  separate dust and synchrotron polarized emission and quantify the correlation between the two
  sources of emission."
- Sect. 4.1 (p. 8): "All 28 spectra are computed for each of the six sky regions described in Sect.
  3.1, within nine multipole bins in the range 4 <= l <= 159."
- Sect. 4.2 (p. 8), what is fixed: "The Planck function B_nu(T_d) is computed for a fixed dust
  temperature T_d = 19.6 K, derived from the fit of the SED of dust total intensity at high Galactic
  latitude in Planck Collaboration Int. XXII (2015). We use a fixed temperature because, over
  microwave frequencies, the dust SED depends mainly on the dust spectral index of the modified
  black-body (or MBB) emission law and the temperature cannot be determined independently of the
  spectral index."
- Sect. 4.2 (p. 9), what is fitted: "The spectral model has five parameters: the two amplitudes A_s
  and A_d and the two spectral indices beta_s and beta_d, characterizing the synchrotron and dust
  SEDs, respectively; and the correlation factor rho quantifying the spatial correlation between
  synchrotron and dust polarized emission." "We fit our spectral model to the EE and BB spectra
  separately, for each sky region and for each multipole bin independently." "We then compute the
  weighted mean and standard deviation of the MPFIT results for all sky regions and multipole bins,
  finding beta_s = -3.13 +- 0.13." "We use it as a Gaussian prior for a second fit of the same data
  with the same model."
- Sect. 4.2 (p. 9), band integration: "These are very close to the factors in Table 3 of Planck
  Collaboration Int. XXII (2015), though here at 353 GHz the evaluation is for the PSBs only. The
  values of C are evaluated for the following SED. For the LFI and WMAP channels used, the
  synchrotron component dominates, for which we assume beta_s = -3, while for the Planck HFI
  channels the polarized dust MBB spectrum dominates, for which we assume beta_d = 1.5 and T_d =
  19.6 K." Table 2 (p. 10), "Unit conversion factors and colour corrections", at 353 GHz: U = 0.075,
  C = 1.098. "Accordingly, the spectrum is converted into units of the data by multiplication by
  C/U".
- Sect. 5.1 (p. 13), the decisive passage: "The difference between beta_d^EE and beta_d^BB is small
  and so we averaged them. Specifically, the uncertainty-weighted average of the fit results for all
  multipole bins and sky regions is beta_d^P = 0.5 (beta_d^EE + beta_d^BB) = 1.53 +- 0.02, where the
  error bar includes the uncertainty from the polarization efficiencies of HFI (Sect. 2) and the
  uncertainty from the CMB subtraction, which affects the determination of beta_d^EE. This is the
  uncertainty of the mean; the weighted dispersions of individual measurements are 0.046 and 0.034
  for E and B modes, respectively. This value of beta_d^P is lower than the mean polarization index
  1.59 +- 0.02 derived from the analysis of earlier (PR2) Planck data (Planck Collaboration Int.
  XXII 2015)."
- Sect. 5.2 (p. 13), SMICA on the LR71 region: "These SEDs, coming from blind component separation,
  are remarkably close to a single-temperature MBB over the full range of Planck polarization
  observations, despite the fact that an MBB spectral shape was not a prior assumption. Performing
  MBB fits after the fact to the SMICA dust spectral data in Fig. 14 (again with T_d = 19.6 K and
  using colour corrections as described in Sect. 4.2), we find a mean spectral index of beta_d^P =
  1.53 +- 0.02, taking into account the 1.5 % uncertainty on the polarization efficiency at 353
  GHz."
- Table 3 (p. 14), SMICA dust SEDs normalized to 1 at 353 GHz, in thermodynamic units. At 44.1
  GHz: EE 6.0 +- 1.1 x 10^-3, BB 2.7 +- 1.4 x 10^-3, MBB 4.6 +- 0.2 x 10^-3. At 70.4 GHz: EE 9.46 +-
  0.75 x 10^-3, BB 9.97 +- 0.96 x 10^-3, MBB 9.76 +- 0.31 x 10^-3. Note b: "Corresponding intensities
  for an MBB SED with T_d = 19.6 K and beta_d = 1.53 +- 0.02."
- Sect. 5.3 (p. 14): "The corresponding value for total intensity is beta_d^I = beta_d^TT = 1.48,
  with much smaller uncertainty. The spectral indices for polarization and total intensity differ by
  0.05 +- 0.03."
- Sect. 5.4 (p. 14): "Therefore, to the sensitivity of the Planck data, a single temperature MBB
  emission law is a satisfactory model of the polarized dust emission." "There is no evidence for a
  flattening or steepening of the dust SED".

Verdict: CONFIRMED. beta_d = 1.53 +- 0.02 with T_d fixed at 19.6 K. The temperature is not fitted;
it comes from the intensity SED of Planck Int. XXII (2015).

Scope to state when citing:
- The value is a weighted mean over six nested high-latitude regions and nine multipole bins (l = 4
  to 159) of fits to EE and BB cross-spectra of PR3 maps at 23, 28.4, 33, 100, 143, 217 and 353 GHz
  (28 spectra from 7 maps (inference)), with the synchrotron index under a Gaussian prior. It is not
  a map-domain or Galactic-plane measurement. Two of the four regions of Paper A1 (Plane, Galactic
  Center) lie outside the sky it was measured on.
- Below 100 GHz: yes, something is said. The abstract claims validity "to below 70 GHz", backed by
  the SMICA SEDs of Table 3 at 70.4 and 44.1 GHz on the LR71 region. At 44.1 GHz the EE and BB
  points have errors of 18 and 52 percent and differ from each other by a factor 2.2 (inference from
  Table 3), so the test below 70 GHz is weak.
- The 353 GHz maps are PSB only and carry a 1.5 percent photometric uncertainty in polarization.
- Dispersion of individual measurements: 0.046 (E) and 0.034 (B).

(inference) Cross-check of the draft's factor 1.12 for the 353 GHz band. Table 2 gives the
conversion from a model in brightness temperature at 353 GHz to the data in K_CMB as C/U =
1.098/0.075 = 14.64. For a delta band at 353.0 GHz the same conversion is (e^x - 1)^2/(x^2 e^x) =
12.91 with x = h nu/(k 2.7255 K) = 6.216. The ratio is 1.134, or 1.127 to 1.142 given the rounding
of U to two digits. The draft's 1.12 is consistent with this to about 1 percent. The Planck factors
are for beta_d = 1.5, PSBs only and the PR3 bandpass, so exact agreement is not expected.

## 6. Watts et al. 2023, Cosmoglobe DR1 results I

Identifier: arXiv 2303.08095v1 (2023-03-14), "Cosmoglobe DR1 results. I. Improved Wilkinson
Microwave Anisotropy Probe maps through Bayesian end-to-end analysis". A&A 679, A143 (2023), DOI
10.1051/0004-6361/202346414, online 2023-11-29 (arXiv API and Crossref; Crossref prints the title
as "COSMOGLOBE DR1 results: I. Improved Wilkinson Microwave Anisotropy Probe maps through Bayesian
end-to-end analysis"). 33 authors. First five on arXiv v1: D. J. Watts, A. Basyrov, J. R. Eskilt,
M. Galloway, L. T. Hergt; in the Crossref record of the published paper: Watts, Basyrov, Eskilt,
Galloway, E. Gjerlow. The bibitem in `refs.tex` (three names and et al.) is correct for both.

Read: HTML-S of arXiv v1 on ar5iv, characters 0 to 300 000 of 310 693 in three portions; the first
portion queried twice. The PDF exceeds the fetch size limit. The published version was not read, so
section and figure numbers are those of v1.

### 6b. K band against 30 GHz in polarization

Quotes (HTML-S, pieces):
- Abstract: "Similarly, the long standing discrepancy between the WMAP K-band and LFI 30 GHz maps is
  finally resolved," "and the difference between the two maps appears consistent with instrumental
  noise at high Galactic latitudes." Also: "our WMAP polarization maps differ more notably from the
  WMAP9 results," "and in general exhibit significantly lower large-scale residuals." "We attribute
  this to a better constrained gain and transmission imbalance model."
- Sect. 2.4: "the polarization modes that are poorly measured by WMAP alone are well measured by
  Planck" "and vice-versa, resulting in an overall better constrained fit." "One important example
  of this regarding the WMAP data is a strong degeneracy between the" "transmission imbalance
  factors and the polarized sky signal;" "it is exceedingly difficult to break this degeneracy using
  data from only one DA alone,"
- Sect. 5.1: "These differences do not match known Galactic component morphologies, but are more
  reminiscent of the poorly measured transmission imbalance modes" (citing Jarosik et al. 2011).
- Sect. 5.2: "we can compare them by scaling K-band assuming a polarized synchrotron power law SED
  of beta_s=-3.1". "the K-30 difference has received significant attention ever since the Planck
  2015 data release" "showing clear signatures of instrumental systematics." "Still, even after all
  these developments, large-scale residuals remained" "When comparing the official WMAP maps,
  particularly for K-band, we see structures sweeping across large angular scales across the sky,"
  "likely due to the poorly measured modes in K-band."
- Sect. 5.3: "The WMAP9 likelihood deweights estimates of the transmission imbalance modes" "in the
  low resolution polarization likelihood covariance matrix". "the two linear templates are able to
  account for between 6 and 52 % of the full large-scale difference" "between Cosmoglobe and WMAP".
  "Overall, it seems likely that the template-based correction of the low-resolution covariance
  matrices" "used in the cosmological WMAP likelihood does not provide a complete description" "of
  the full large-scale uncertainties."
- Sect. 6.4: "for polarized synchrotron emission, K-band is stronger than LFI 30 GHz by about 25 %."
- Sect. 7.3: "the WMAP9 solution does not rely on any knowledge of the sky," "but as a result induces
  poorly measured modes with a 2.5 uK amplitude."
- Sect. 8: "This work provides a natural resolution to the long-standing discrepancy" "between the
  polarized WMAP K-band and Planck LFI 30 GHz observations". "These two channels are sufficiently
  close in frequency that any uncertainty in the synchrotron SED" "should be subdominant to
  instrumental noise." "However, when differencing these two maps (after scaling one by the
  synchrotron SED)," "large-scale residuals correlated with the scanning strategy appeared." "It is
  only now, when both LFI and WMAP are processed from scratch," "that the two datasets agree to a
  level compatible with instrumental noise." "This progress is illustrated in Fig. 46, which shows
  the K-30 difference maps" "for Planck 2018, Planck PR4, BeyondPlanck, and Cosmoglobe." "In all but
  the last case, the Planck maps are differenced with the WMAP9 K-band map." Also: "In
  polarization, we see differences of up to ~10 uK." "we have also shown that the bilinear template
  marginalization method used in WMAP9" "is unable to account for the full observed differences."
- Fig. 46 caption: "the K-band map has been scaled by 0.495 to account for different central
  frequencies," "assuming a synchrotron spectral index beta_s=-3.1."

Summary of (b): before reprocessing, the K minus 30 GHz difference (K scaled by 0.495 for an
assumed index of -3.1) shows large-scale residuals "correlated with the scanning strategy" for
Planck 2018, PR4 and BeyondPlanck against WMAP9 K band. After the joint reprocessing the paper says
the two "agree to a level compatible with instrumental noise", in the abstract qualified to "at
high Galactic latitudes". The index -3.1 is assumed for the scaling; the paper does not fit an
index between the two bands. It attributes the WMAP9 side of the discrepancy to "poorly measured
modes in K-band" and to transmission imbalance, and states that the nine-year covariance treatment
"does not provide a complete description of the full large-scale uncertainties".

That last statement bears directly on Paper A1, which uses the official nine-year maps and N^-1
matrices: it is a published claim, by a different analysis, that those matrices understate the
large-scale errors. The draft should cite it as such and not only as background.

(inference) A scale factor of 0.495 for index -3.1 corresponds to a frequency ratio of 1.2546 for
a pure power law in brightness temperature, for example 28.4 over 22.64 GHz. The text read gives
neither the two frequencies nor the temperature unit of the scaling.

### 6c. Reference frequency and prior of the polarized synchrotron index

Quotes (HTML-S, pieces):
- Sect. 2.4: "when fitting both the polarized synchrotron amplitude and calibration freely without
  priors, the synchrotron" "spectral index at high Galactic latitudes tend (sic) to be biased toward
  unreasonably flat values, beta_s <~ -2.5," "which was likely due to a low level of unmodeled
  systematics, for instance temperature-to-polarization leakage," "rather than true polarized
  synchrotron emission." "To break this degeneracy, Svalheim et al. 2023b chose to marginalize the
  high-latitude synchrotron spectral index" "over a Gaussian prior of N(-3.30, 0.1^2), informed by
  Planck Collaboration V 2020," "rather than estimate it from the data themselves." "We observe the
  same degeneracy, and the introduction of the K-band data is not sufficient to break it on its
  own." The paper then applies "the same informative prior" and says: "These strong priors imply
  that the foreground spectral parameters sampled within the chain do not carry independent
  significance in the traditional posterior sense."
- Sect. 6.1: "this does imply that the SED parameters that are sampled as part of the Gibbs chain
  are non-informative." "Rather, spectral parameters must be estimated through external analyses
  from the frequency maps," "and this is for instance done for polarized synchrotron emission in a
  companion paper by Fuskeland et al. 2023."
- Reference frequency: the sky model has nu_0,s with "nu_0,c is a reference frequency for component
  c" and no number in the text read. Fig. 31 caption: "Polarized synchrotron maps and their standard
  deviations evaluated at 30 GHz." Fig. 29 caption (intensity): "Synchrotron emission evaluated at
  408 MHz."

Verdict: CONFIRMED on the record and title (a); (b) as summarized; (c) prior N(-3.30, 0.1^2) on the
high-latitude polarized synchrotron index, taken over from the BeyondPlanck analysis, with the
index sampled inside the chain but prior dominated. A reference frequency of 30 GHz for the
polarized synchrotron amplitude is supported only by the caption of Fig. 31. This agrees with
`nextgen/LITERATURE.md` line 112 (DR1 I states N(-3.30, 0.1^2); the DR1 III value was not checked
here).

## 7. Planck 2018 results II (LFI processing): low-resolution NCVM and 30 GHz systematics

Identifier: arXiv 1807.06206v2 (11 Sep 2018), "Planck 2018 results. II. Low Frequency Instrument
data processing". A&A 641, A2 (2020), DOI 10.1051/0004-6361/201833293 (arXiv API).

Read: PDF of arXiv v2, pp. 1 to 4, 11, 19 to 20, 22 to 23, 26 to 27 and 30 to 31 (the table of
contents puts App. B.3 on p. 33). HTML-S keyword search of the ar5iv rendering, characters 0 to
200 000 of 205 994, for "noise covariance". Official documentation, HTML-S: PLA Explanatory
Supplement page `Map-making_LFI` (2018 version) and the IRSA page "Planck Public Data Release 3
Mission Ancillary Data".

Is there a public low-resolution NCVM for the 2018 LFI maps? Yes, but the paper does not describe
it.
- The paper mentions NCVMs once, Sect. 8.3 (p. 19): "The actual noise in the LFI data is given
  directly by the half-ring difference maps. Detailed noise characterization is of paramount
  importance for the creation of adequate noise-covariance matrices (NCVMs), as well as for the
  noise MC realizations that are required in subsequent steps of the data analysis. Any noise model
  has to be validated against such half-ring difference maps."
- Sect. 10, "LFI data products available through the Planck Legacy Archive" (pp. 26 to 27), lists
  pointing timelines, time timelines, housekeeping timelines, timelines in volts, cleaned and
  calibrated timelines, scanning beam, effective beam, full-sky maps, bandpass correction maps, the
  gain correction template and the RIMO. Low-resolution maps and NCVMs are not in the list.
- PLA Explanatory Supplement, `Map-making_LFI` (HTML-S): "This data set has been packed into three
  different files, one per frequency", named "LFI_NoiseCovMat_0??_0016_R3.00.tgz" in the body text
  and "Contents of LFI_NoiseCovMat_0??_0016_R3.20.tgz" in the table heading, "that can be downloaded
  from the Cosmology section of the Planck Legacy Archive." "The low-resolution data set can
  currently be utilized efficiently only at resolution" N_side = 16 or lower. The matrices "are C
  binary format files;" "are organized in block form," (II, IQ, IU, QI, QQ, QU, UI, UQ, UU) "are in
  the HEALPix nested pixelisation scheme" "have KCMB units." (as written) "are in Galactic
  coordinates;". "The output of the NCVM module of MADAM mapmaker are inverse NCVMs." "The noise
  model was determined by three noise parameters: the white noise level sigma; slope; and knee
  frequency" "baseline lengths were 0.25s for 30GHzand (sic) 1.0s for 44GHz, and 70GHz". For the
  maps: "subsequently the temperature part was smoothed with a symmetric Gaussian beam with FWHM =
  440arcmin." Files with and without regularization noise are listed.
- IRSA, PR3 ancillary data (HTML-S): `LFI_NoiseCovMat_030_0016_R3.20.tgz`,
  `LFI_NoiseCovMat_044_0016_R3.20.tgz`, `LFI_NoiseCovMat_070_0016_R3.20.tgz`; "NSIDE=16 tarfiles, 1
  per frequency. File sizes are 2.5GB."

The 30 GHz map and large-scale systematic modes (PDF):
- Sect. 8.6, "Validation Summary" (pp. 22 to 23): "At the end of this detailed validation process,
  the improved quality of the 2018 data release, especially at 30 and 44 GHz, is clear. This
  improvement is mainly due to the new calibration procedure implemented. It is also evident that
  residuals at very low multipoles (large angular scales) are still present at 30 and at 44 GHz. One
  of the reasons could be the fact that the initial foreground model used for the iterative
  calibration, mapmaking and component-separation process was entirely based on LFI data. We took
  this approach to avoid dependence on WMAP and to avoid ingesting any HFI systematic effects.
  However, we realize that this results in a foreground model with less power than expected. This
  choice should in principle be investigated further, but the actual level of residuals is small and
  does not prevent the use of the 30-GHz data as a synchrotron template in cosmological analyses
  involving the 70-GHz LFI channel (see Planck Collaboration V 2018, for a detailed analysis of the
  low-l likelihood)".
- App. A, "Comparison of LFI 30 GHz with WMAP K and Ka bands" (pp. 30 to 31): "In order to provide a
  further check on systematic effects, we present here a comparison in polarization at large scales
  between the LFI 30-GHz channel and the K and Ka bands from WMAP (Bennett et al. 2013)." "We then
  filter Planck and WMAP maps with a 10 deg-FWHM Gaussian beam in order to suppress high-frequency
  noise and highlight the large-scale structures. We rescale the K and Ka bands to the LFI 30-GHz
  effective frequency, assuming synchrotron emission with a spectral index of -3." "The amplitude
  of the large-scale structure in the difference maps decreases by roughly a factor of two between
  the 2015 and the 2018 maps, reflecting the improvements in the new calibration scheme." "However,
  as pointed out in Sect. 3.2, the calibration could not be run to convergence for practical
  reasons, and we expect residuals in the 2018 maps with a pattern similar to that of the 2015
  maps." "There is also some agreement between features in the difference maps and the 70-GHz gain
  correction template (Fig. 6), which is based in part on the 30-GHz maps. Despite the substantial
  improvement in the 2018 30-GHz polarization maps, residual systematics are still present, and
  included in our estimated error budget. Additional improvements are planned in a future paper."
- Sect. 8.3 (p. 20), half-ring differences against the white-noise matrices at 1150 <= l <= 1800:
  "there is still an excess of 1/f noise in this high-l regime, meaning that both the real data and
  the noise MCs predict slightly larger noise than the WNCVM." "Residuals are <~ 1.4 % at 30 GHz,
  <~ 1 % at 44 GHz, and <~ 0.6 % at 70 GHz, for both temperature and polarization."
- Sect. 10 (pp. 26 to 27), convention: "We note that the Planck adopted polarization convention is
  not the one proposed by "IAU" but the one used more generally in CMB full-sky maps and referred to
  as "COSMO" (see Planck Collaboration ES 2018)."
- Table 1: see item 4d.

Verdict: CONFIRMED WITH CORRECTION.

Corrected statements:
- Public N_side 16 noise covariance matrices exist for the 2018 LFI maps at 30, 44 and 70 GHz,
  distributed through the PLA (Cosmology section) and IRSA as `LFI_NoiseCovMat_0??_0016_R3.??.tgz`
  (version string R3.00 in the PLA text, R3.20 in the PLA table heading and at IRSA), about 2.5 GB
  each, IQU block form, nested, Galactic. They are documented in the PLA Explanatory Supplement.
  Planck 2018 II neither lists nor describes them, so cite the Explanatory Supplement (and Planck
  2015 II or 2018 V for the method) for the product.
- On the 30 GHz map, the paper's own words are that residuals at very low multipoles "are still
  present at 30 and at 44 GHz" (Sect. 8.6) and that in the 2018 30-GHz polarization maps "residual
  systematics are still present, and included in our estimated error budget" (App. A). Weiland et
  al. 2022 write that "the large-scale systematic modes are primarily attributable to the 30 GHz
  map". (inference) That is their reading of App. A, where the K minus 30 structure halves when only
  the LFI calibration changes; the paper does not say "primarily" and does not assign the residual
  between WMAP and LFI.

Relevance: these matrices belong to the PR3 maps (PR3 destriping baselines, PR3 calibration), not
to the PR4 maps the draft uses. PR4 has its own N_side 16 matrices for cosine-smoothed maps (NPIPE
Sect. 4.2, item 4a).

## 8. Flags against the repo notes and the draft

1. Task statement, item 3, and any future bibitem: "Jarosik et al. 2003, First-Year WMAP
   Observations: On-Orbit Radiometer Characterization, ApJS 145, 413" is a mixed record. Table 11
   is in "Design, Implementation, and Testing of the Microwave Anisotropy Probe Radiometers", ApJS
   145, 413. `refs.tex` has no Jarosik entry yet.
2. Task statement, item 4b: "30 GHz (and 217, 353 GHz) maps are calibrated without the
   polarized-sky prior" is wrong for 217 GHz. Only 30 and 353 GHz are prior-free.
3. `sec_data.tex`, WMAP paragraph: cites Bennett et al. 2013 alone for the dense matrices and what
   they carry. The content comes from Jarosik et al. 2007 and 2011 and the LAMBDA page (item 2a).
4. `sec_data.tex`, Planck paragraph: "ring-half" should be "half-ring" (item 4a).
5. `sec_data.tex`, spectra paragraph: "the band of the polarization-sensitive bolometers" at 353 GHz
   is the PR3 description. For PR4 the SWBs appear to be included (item 4c, inference).
6. `sec_data.tex`, spectra paragraph: the modified blackbody of Planck 2018 XI is a high-latitude
   power-spectrum mean with T_d fixed, not fitted (item 5). The sentence is accurate as written;
   the scope is missing.
7. `2026-11-09_science_landscape_and_data.md`, Sect. 4, WMAP9 row: "singular in processing-mask
   pixels (ninv_info.html, as extracted)". The sentence is on the LAMBDA page, but Jarosik et al.
   2011 and E038 indicate that the band matrices from the seven-year release on cover all pixels
   (item 2a). The same row says "N^-1 includes loss-imbalance modes"; more exactly, the nine-year
   raw Q/U matrices carry them with "a true 1-sigma uncertainty", where earlier releases projected
   them out.
8. `03_instrument_and_external_data.md`, line 13 ("foreground polarization priors (30/217/353 GHz
   templates)"): consistent with NPIPE Sect. 2.4.13. Add that 217 GHz is itself calibrated with the
   prior.
9. `nextgen/LITERATURE.md`, line 120 ("low-multipole transfer-function caveat at 44 and 70 GHz"):
   the transfer function applies to all channels calibrated with the prior, 44 to 217 GHz, and not
   to 30 or 353 GHz (item 4b).
10. `refs.tex`: Bennett, NPIPE and Planck 2018 XI entries lack arXiv numbers (1212.5225, 2007.04997,
    1801.04945); Huber lacks the end page (233); White lacks the end page (25).

No statement in the repo notes that I re-read is contradicted outright by the sources.

## 9. What was looked for and not found

- A description of the r4 Q/U maps and their N^-1 matrices in Bennett et al. 2013 (PDF pp. 1 to 26
  and 133 to 138; HTML-S keyword search of characters 0 to 300 000).
- A statement of the Stokes convention in Bennett et al. 2013 (same search; no "IAU", no
  "convention" in that sense).
- The nine-year change from projecting loss-imbalance modes to a 1-sigma term, in any paper. Found
  only on the LAMBDA page.
- Whether the released QU matrix is conditional on or marginal over I and S (Jarosik et al. 2007,
  2011; Page et al. 2007 App. D by HTML-S).
- An Eq. (50), a definition of "effective frequency for a power law", or a rule for combining
  radiometers into a band in Jarosik et al. 2003 (pp. 13 to 24).
- An explicit sentence on which bolometers enter the NPIPE 353 GHz maps (four queries over the
  three portions of the ar5iv text).
- LFI knee frequencies in the NPIPE text (Sects. 1 to 7.4).
- A statement in the NPIPE text that the half-ring difference over- or underestimates LFI noise, or
  a description of per-pixel white-noise covariance columns.
- Low-resolution NCVMs in the product list of Planck 2018 II (Sect. 10, PDF).
- A numerical nu_0,s for synchrotron in the text of Watts et al. 2023 (characters 0 to 300 000).
- A fitted spectral index between K band and 30 GHz in Watts et al. 2023; -3.1 is assumed.
- The word "sandwich" in Huber (1967) or White (1982), on the pages read.
- A DOI for Huber (1967).

Not attempted or not reached: the IDL routine `effective_freq.pro`; the WMAP nine-year Explanatory
Supplement; NPIPE Sects. 8 and 9 and appendices; the published versions of 2303.08095 and
2007.04997; White (1982) pp. 2 to 3 and 8 to 25 and its corrigendum; Huber (1967) pp. 222 to 225;
the last 10 693 characters of 2303.08095 and the last 14 674 of 1212.5225 on ar5iv.

Sources tried without success: `arxiv.org/pdf/2007.04997` and `arxiv.org/pdf/2303.08095` (above
10 MB); `www.aanda.org/articles/aa/full_html/2020/11/aa38073-20/aa38073-20.html` (403);
`orca.cardiff.ac.uk/id/eprint/136666/1/aa38073-20.pdf` (403);
`www.research.unipd.it/retrieve/.../aa38073-20_part_1.pdf` (403);
`projecteuclid.org/euclid.bsmsp/1200512988` and the Project Euclid ebook page (empty content);
`users.ssc.wisc.edu/~bhansen/718/White1982.pdf` (404).

Searches run. arXiv API (`export.arxiv.org/api/query`), `id_list` for 1212.5225, 2007.04997,
1801.04945, 2303.08095, 1807.06206, astro-ph/0301164, 1807.06205, astro-ph/0603450,
astro-ph/0603452, 1001.4744. Crossref (`api.crossref.org/works/<DOI>`) for 10.1086/346080,
10.1086/377221, 10.1088/0067-0049/208/2/20, 10.2307/1912526, 10.1051/0004-6361/202346414,
10.1051/0004-6361/202038073, and one bibliographic query for Huber (1967). Web searches: (1) NPIPE
353 GHz SWB inclusion, twice with different wording; (2) PLA wiki pages for NPIPE frequency maps
and half-ring files; (3) Planck 2018 LFI low-resolution NCVM product names; (4) full text of White
(1982), twice. INSPIRE was not needed.

## 10. Gaps and suggested follow-up

1. NPIPE from the published PDF (a local or library copy): confirm every HTML-S quote of item 4,
   read Sect. 4.1.1 and the appendix on data products for the definition of the covariance columns
   and of the `wcov`, `hrscaled` and `mcscaled` files, and find the sentence, if any, on SWBs in the
   353 GHz maps. The NPIPE RIMO or the reduced instrument model weights would settle 4c without the
   paper.
2. WMAP nine-year Explanatory Supplement (Greason et al. 2012) for the r4 N^-1 products: the
   1-sigma loss-imbalance term, whether the QU matrix is conditional on or marginal over I and S,
   and the mask. Check `POLCCONV` and the diagonal of the K-band matrix in the local files.
3. Fetch `effective_freq.pro` from the LAMBDA IDL library and record its header and constants
   against App. A of Bennett et al. 2013.
4. Check the journal version of Table 11 of Jarosik et al. 2003 (column heading of the fourth
   column, and the K to V values transcribed in item 3), and read pp. 1 to 12 for the definition of
   the effective bandwidth.
5. Decide what the draft says about the calibration-stage constant-index assumption for PR4 44 and
   70 GHz (item 4b, inference). A test that would answer it: repeat the region index with PR3 or
   BeyondPlanck 44 and 70 GHz maps and compare.
6. Compare the draft's fitted three-parameter LFI noise model with the released PR4 N_side 16
   noise covariance matrices (NPIPE Sect. 4.2), allowing for their cosine smoothing; say in the
   paper why the released matrices are or are not used.
7. Cite Watts et al. 2023 Sect. 5.3 where the draft relies on the nine-year N^-1 matrices at the
   largest scales: it is a published claim that those matrices do not fully describe the
   large-scale uncertainty.
8. Library reads: White (1982) Assumptions A1 to A6 and the 1983 corrigendum; the publication
   record of the Fifth Berkeley Symposium, Vol. 1, for editors, publisher and year.
9. Read the published version of Watts et al. 2023 for the section and figure numbers used in
   item 6, and its sky-model table for nu_0,s.
