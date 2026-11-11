# Citation verification: cluster-robust (block) standard errors, 5 references and 3 claims

Date: 2026-11-11. Author: lit-scout agent. For Paper A1, `sec_data.tex`, Eq. `eq:block` and its
`\pend{references}`. Nothing on this topic was in `nextgen/LITERATURE.md` or `nextgen/research/`.
Read levels: PDF = pages read as images, quotes transcribed by me. HTML-S = official page read through the
fetch tool's small model. ABS = abstract only. BIB = record only. "(inference)" marks my own synthesis.

## 1. References

| # | Verified record | Source URL | Read | Status |
|---|---|---|---|---|
| 1 | K.-Y. Liang, S. L. Zeger (1986), "Longitudinal data analysis using generalized linear models", Biometrika 73(1), 13-22 (April 1986), DOI 10.1093/biomet/73.1.13 | https://academic.oup.com/biomet/article-abstract/73/1/13/246001 ; https://api.crossref.org/works/10.1093/biomet/73.1.13 | ABS, BIB | Record verified. Full text NOT read (OUP paywall; a hosted copy returned 403) |
| 2 | R. M. Bell, D. F. McCaffrey (2002), "Bias reduction in standard errors for linear regression with multi-stage samples", Survey Methodology 28(2), 169-181 (December 2002), Statistics Canada Catalogue No. 12-001-XPB, online article 12-001-X20020029058. No DOI found | https://www150.statcan.gc.ca/n1/en/pub/12-001-x/2002002/article/9058-eng.pdf ; https://www150.statcan.gc.ca/n1/en/catalogue/12-001-X20020029058 | PDF, all 13 pages | Verified from the article header. Ref. 4's own reference list prints 169-179; the article says 169-181, use that |
| 3 | J. G. MacKinnon, H. White (1985), "Some heteroskedasticity-consistent covariance matrix estimators with improved finite sample properties", J. Econometrics 29(3), 305-325 (September 1985), DOI 10.1016/0304-4076(85)90158-7 | https://api.crossref.org/works/10.1016/0304-4076(85)90158-7 (DOI resolves to Elsevier PII 0304407685901587) | BIB | Record verified (Crossref; same in the reference lists of refs. 2 and 4). Paper NOT read: ScienceDirect 403, Queen's WP 537 TLS error |
| 4 | A. C. Cameron, D. L. Miller (2015), "A Practitioner's Guide to Cluster-Robust Inference", J. Human Resources 50(2), 317-372 (March 2015), DOI 10.3368/jhr.50.2.317 | https://jhr.uwpress.org/content/50/2/317 ; text read: https://cameron.econ.ucdavis.edu/research/Cameron_Miller_JHR_2015_February.pdf | BIB + PDF of the authors' final manuscript (ms pp. 1-12, 23-34, 45-52) | Record verified. Page numbers below are manuscript pages; journal pages of the passages could not be confirmed |
| 5 | A. V. Carter, K. T. Schnepel, D. G. Steigerwald (2017), "Asymptotic Behavior of a t-Test Robust to Cluster Heterogeneity", Rev. Econ. Stat. 99(4), 698-709 (October 2017), DOI 10.1162/REST_a_00639 | https://api.crossref.org/works/10.1162/REST_a_00639 ; text read: authors' draft of 2016-02-01 on faculty.econ.ucsb.edu/~doug/researchpapers/ | BIB + PDF of the draft (pp. 1-13, 17-21) | Record verified (Crossref, EconPapers; MIT Press page 403). Equation numbers of the printed version not confirmed |
| 6 | A. C. Cameron, J. B. Gelbach, D. L. Miller, "Bootstrap-Based Improvements for Inference with Clustered Errors", NBER Technical Working Paper 344 (September 2007); published Rev. Econ. Stat. 90(3), 414-427 (2008) per the reference lists of ref. 4 and the sandwich docs | https://www.nber.org/system/files/working_papers/t0344/t0344.pdf | PDF of the NBER version, pp. 1-10 | Supplementary. Journal record not checked at the publisher |
| 7 | J. G. MacKinnon, M. O. Nielsen, M. D. Webb, "Cluster-Robust Inference: A Guide to Empirical Practice", arXiv:2205.03285v1 (2022-05-06), arXiv comment "Accepted at Journal of Econometrics" | https://arxiv.org/abs/2205.03285 | PDF pp. 5-11 (Sects. 2.1 to 3.3) | Supplementary. Journal volume and pages not checked |
| 8 | Software docs: clubSandwich `vcovCR`; sandwich `vcovCL` | https://jepusto.github.io/clubSandwich/reference/vcovCR.html ; https://sandwich.r-forge.r-project.org/reference/vcovCL.html | HTML-S | Supplementary |

## 2. Claims

### A. The cluster sandwich originates with, or is standardly attributed to, Liang and Zeger (1986)

Verdict: PARTLY CONFIRMED. "Standardly attributed to" holds; "originates with" does not; the formula inside
ref. 1 was not read (abstract only: estimating equations with consistent estimates of the variance).
- Ref. 4, Sect. I (ms p. 2): cluster-robust standard errors were "proposed by White (1984, p.134-142) for OLS
  with a multivariate dependent variable (directly applicable to balanced clusters); by Liang and Zeger (1986)
  for linear and nonlinear models; and by Arellano (1987) for the fixed effects estimator in linear panel
  models." Sect. II.C (ms p. 8), Eqs. (10)-(11): (X'X)^-1 [sum_g X_g' u_g u_g' X_g] (X'X)^-1, "Initial
  derivations of this estimator, by White (1984, p.134-142) for balanced clusters, and by Liang and Zeger
  (1986) for unbalanced". Its reference list gives Biometrika 73(1): 13-22. Heteroskedastic case: White
  (1980), Econometrica 48(4), 817-838 (Sect. II.A; record as printed there, not checked).
- Ref. 6, Sect. 1 (p. 3): "In the statistics literature these are called sandwich standard errors, proposed by
  Liang and Zeger (1986) for generalized estimating equations". Ref. 8 (clubSandwich): CR0 "is the original
  form of the sandwich estimator (Liang & Zeger, 1986)".
- Against "originates": ref. 5 (draft p. 5), "The cluster-robust variance estimator, which Shah, Holt and
  Folsom (1977) are among the first to use"; ref. 2 (p. 169) treats it as the survey "linearization" estimator
  (Fuller 1975) and cites for GEE the other 1986 paper, Zeger and Liang, Biometrics 42, 121-130.

### B. Downward bias of CR0; CR2 = (I - H_gg)^(-1/2); CR3 = (I - H_gg)^(-1) and the jackknife; names

Verdict: CONFIRMED, with the naming and the jackknife relation made precise. H_gg = X_g (X'X)^-1 X_g'.
- Bias: ref. 2 abstract (p. 169): "when the number of PSUs is small or a coefficient depends primarily on data
  from a small number of PSUs, linearization estimators can have large negative bias." Theorem 1 (p. 170):
  for V = sigma^2 I and c = n/(n-1), E(v_L) <= Var, with equality only if l'(X'X)^-1 X_i'X_i is constant
  across clusters. Table 1 (p. 175), 20 clusters, iid errors: bias of the linearization variance -9.6, -13.2,
  -32.5, -13.3, -1.8 percent for the five coefficients; jackknife +11.7, +17.2, +51.2, +17.6, +2.1; BRL 0.0.
  Unequal leverage (p. 172): "Even greater bias potential appears to occur when certain PSUs account for most
  of the variability in the covariates".
- Forms: ref. 2, Sect. 4 (p. 173): v_L, v_BRL, v_JK are all c l'(X'X)^-1 [sum_i X_i' A_i r_i r_i' A_i X_i]
  (X'X)^-1 l with "c = n/(n-1), 1, or (n-1)/n, respectively, and A_i = I_i, (I_i - H_ii)^-1/2, or (I_i -
  H_ii)^-1"; "v_BRL can be viewed as a compromise between v_L and v_JK". BRL: Eq. (7), A_i'A_i = (I_i -
  H_ii)^-1, symmetric root; "a generalization of the method of MacKinnon and White (1985)" (p. 172).
- Jackknife: ref. 2, Eqs. (5)-(6) (p. 171). In linear regression the identity is EXACT, not approximate, when
  the delete-one estimates are centred on the full-sample estimate, and it carries the factor (n-1)/n.
  Centring on their mean gives "very similar estimates". Theorem 2: E(v_JK) >= Var under iid errors. The form
  WITHOUT the factor is Mancl and DeRouen (2001): it "equals n/(n-1) v_JK" (p. 172).
- Names: ref. 2 uses none of CR0, CR2, CR3, HC2, HC3 (whole paper read); it says linearization, BRL, jackknife.
  "CR3": earliest use I could confirm is ref. 6, Sect. 2.2 (p. 10): "This correction generalizes the HC3
  measure (jackknife) of MacKinnon and White (1985), and so we refer to this correction as the CR3 variance
  estimator." "CR2VE", "CR3VE": ref. 4, Sect. VI.B, Eqs. (24)-(25) (ms p. 25), "so we refer to it as the
  CR2VE"; its Tables 2-3 print "CR2 bias correction", "CR3 bias correction". "CR0, CR1, CR1p, CR1S, CR2,
  CR3": ref. 8 (clubSandwich); sandwich `vcovCL` keeps HC0 to HC3. Ref. 7, Sect. 2.2, Eqs. (8)-(10) uses CV1
  to CV3: "CV3 is actually a jackknife estimator which generalizes the familiar HC3 estimator".
- The factor differs by source: ref. 4 Eq. (25) and ref. 7 Eq. (10) put (G-1)/G in CR3/CV3 and call it equal
  to the jackknife; clubSandwich says CR3 "approximates the leave-one-cluster-out jackknife" (formula not on
  the page read). Preference: ref. 4 (ms p. 25), "The literature has gravitated to using the CR2VE adjustment
  rather than the CR3VE adjustment"; ref. 7 (p. 7), "CV3 typically yields the most reliable tests".
- Not confirmed at source: the HC1, HC2, HC3 definitions of ref. 3 (known only through refs. 2, 4, 6, 7).

### C. Use t critical values with few degrees of freedom (G - 1 or Satterthwaite) when clusters are few

Verdict: CONFIRMED.
- G - 1: ref. 2, Sect. 4 (p. 173) traces it to Shah, Holt and Folsom (1977), "now the default in Stata ...,
  SUDAAN ... and SAS", and says such tests "would tend to have Type I error rates that exceed the nominal
  value, even if v_L were unbiased". Ref. 4, Sect. VI (ms p. 23): "One should at a minimum use T(G-1)
  critical values"; Sect. VI.D.1 (ms p. 29) repeats it against the standard normal.
- Satterthwaite: ref. 2, Sect. 4 (p. 173): f = (sum_i lambda_i)^2 / sum_i lambda_i^2, lambda_i the eigenvalues
  of G = {g_i' V g_j}, computed with V = I; "Tests based on reference t-distributions with f degrees of freedom
  would be expected to provide better Type I error rates than tests based on n - 1". Table 2 (p. 176): f = 2.6
  to 16.7 for 20 clusters. Caveat (pp. 174, 178): it overcompensates (conservative) when f is "less than 4 or
  5". Ref. 4, Sect. VI.D.2, Eq. (26) restates it; further developed by Imbens and Kolesar (2012, NBER WP 18478).
- t with G* degrees of freedom is NOT recommended by ref. 5: it is "the logical implication" but "could lead to
  difficulty in controlling the size of the test" (draft pp. 18-19). Ref. 4, Sect. VI.D.3: "Carter et al.
  (2013) do not actually propose using critical values based on the T(G*) distribution"; T(G*) is ref. 4's own.

## 3. Effective number of clusters

Yes: the standard term and reference is ref. 5. Definition (draft Lemma 1, Sect. 3, Eq. (5)):
G* = G / (1 + Gamma), Gamma = (1/G) sum_g (gamma_g - gamma_bar)^2 / gamma_bar^2,
gamma_g = a' (X'X)^-1 X_g' Omega_g X_g (X'X)^-1 a, a the selection vector of the tested coefficient.
"We refer to G* as the effective number of clusters" (p. 10). Feasible version G*A: Omega_g replaced by a
matrix of ones. Ref. 4 Eq. (27) reproduces it (delta for Gamma) and reports size distortion for G* < 20 with
normal critical values. Ref. 5 advises reporting G* and using Gaussian critical values if it is large.

## 4. Flags for Paper A1

1. (inference) G / (1 + Gamma) = (sum_g gamma_g)^2 / sum_g gamma_g^2, the same form as the paper's (sum_B
   D_B)^2 / sum_B D_B^2 and as the Satterthwaite f of ref. 2. But gamma_g is the cluster's share of the
   VARIANCE of the estimate, not of the information. The paper's number equals G* only if Var(G_B) is
   proportional to D_B; the text gives the score covariance as J + T, so check before citing ref. 5 for it.
2. Eq. `eq:block` divides by (1 - h_B)^2 with no (G-1)/G factor. By ref. 2 (p. 172) that is the Mancl and
   DeRouen form, G/(G-1) times the jackknife variance in the linear case: 1.33 in variance, 1.15 in error for
   4 patches (inference: the arithmetic, and carrying a linear result to a nonlinear estimator).
3. Every result above is for linear regression. Ref. 2 (p. 179): "Theorem 3 does not extend to GLMs because
   the weights are estimated from the data". The sandwich docs (HTML-S) say HC3 and the jackknife coincide in
   linear models and not in nonlinear ones. "Jackknife to first order" is the safe wording.

## 5. Not found, and follow-up

Searches: Crossref for four DOIs; StatCan catalogue; OUP page and PDF; ScienceDirect and doi.org for ref. 3;
MIT Press for ref. 5; EconPapers for refs. 3 and 5; arXiv API (one HTTP 429) and abs page; web searches for
full text of refs. 1, 3, 6 and for the G* definition. Not found or not read: Theorem 2 of ref. 1; the text of
ref. 3; a DOI for ref. 2; any use of "CR2" before ref. 4 (the 2006 working paper of ref. 6 and Pustejovsky and
Tipton 2018, arXiv:1601.01981, not read); the clubSandwich CR3 formula; Sects. 3.4 onward of ref. 7. Follow-up:
refs. 1 and 3 from a library copy; Mancl and DeRouen (2001), Biometrics 57, 126-134 (as printed in ref. 2);
Imbens and Kolesar (2016), Rev. Econ. Stat. 98(4), 701-712 (as printed in ref. 8).
