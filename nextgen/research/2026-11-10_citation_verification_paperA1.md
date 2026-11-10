# Citation verification for Paper A1: 14 statements checked against the papers themselves

Date: 2026-11-10. Author: lit-scout agent.

Question: is each of the 14 statements below supported by the cited paper, where exactly, and in
what words?

Read levels:
- PDF: I read the stated pages of the arXiv PDF as page images and transcribed the quotes myself.
- HTML-S: full-text HTML (ar5iv or arxiv.org/html) read through the fetch tool's small model, which
  returns quotes of at most about 125 characters. Long sentences are therefore given as consecutive
  pieces. Where two renderings were queried independently, this is stated.
- ABS: abstract only. BIB: bibliographic record only (Crossref, arXiv API, INSPIRE, RePEc).

Conventions: Greek letters and mathematics are written in ASCII. Typographic dashes inside quotes
are written as hyphens. "(sic)" marks a misprint that is in the source. "(inference)" marks my own
arithmetic or synthesis, never a result of the cited paper. Section numbers are those of the arXiv
version read, which uses Arabic numerals where the task statement used Roman ones.

Access limits that shaped the read levels: arXiv PDFs above 10 MB could not be fetched
(2011.08503 v1 and v2, 1404.5323v2, 1907.12875, 2201.08188, 2011.05609). The A&A site and the
University of Milan repository returned HTTP 403. JSTOR and Project Euclid pages did not load.

## 0. Summary table

| # | Paper | Read | Verdict | Correction or caveat |
|---|---|---|---|---|
| 1 | Ball, Fink, Bowler 2003, cond-mat/0301179 | PDF, all | CONFIRMED | Eq. (7) holds for the sign rule (accept if and only if the estimate is negative). arXiv v1 read, PRL text not. |
| 2 | Andrieu, Roberts 2009, 0903.5480 | PDF, pp. 1-5, 18-20, 28-29 | CONFIRMED | None. |
| 3 | Bardenet, Doucet, Holmes, 1505.02827 | PDF, pp. 1-2, 20-23 | CONFIRMED WITH CORRECTION | Not a Gaussian-noise statement: Berry-Esseen bound, sign rule, Barker target pi^beta with beta = lambda sqrt(t)/(n sigma), lambda = 1.702, C = 0.0095. |
| 4 | Stompor et al. 2009, 0804.2645 | PDF, pp. 1-6, 15-17 | CONFIRMED | Eq. (9) spectral, Eq. (10) marginal. The printed sign of the determinant term in Eq. (10) looks reversed (inference). |
| 5 | Eriksen et al. 2008, 0709.1058 | PDF, pp. 1-12 | CONFIRMED | Section is 4.2. Jeffreys prior is Eq. (32) in general and Eq. (36) for the index. |
| 6 | Dunkley et al. 2009, 0811.4280 | PDF, all | CONFIRMED | None. Maps at Nside 8 from Nside 16 inputs; remedy is p_s(beta)/L_0(beta). |
| 7 | Alonso et al. 2017, 1608.00551 | PDF, pp. 4-7, 15-19 | CONFIRMED WITH CORRECTION | Text says more than 4 sigma, Fig. 11 caption says more than 3 sigma. Both are in the paper. |
| 8 | Weiland et al. 2022, 2203.11445 | PDF, pp. 1-19 | CONFIRMED WITH CORRECTION | Superpixel Q/U correlation is used only for plane and spurs. High latitude is a per-pixel fit at Nside 16 with 1.4 or 2.3 GHz. |
| 9 | Weiland et al. 2018, 1801.01226 | PDF, pp. 1-10 | CONFIRMED | 2015 Planck data, EE only. The three rows differ in dust template and in fixed or free coefficients. |
| 10 | Svalheim et al., 2011.08503 | HTML-S, two renderings | CONFIRMED WITH CORRECTION | Published title is "BeyondPlanck XIV" (A&A 675, A14). Values are those of v2; v1 had different numbers. No numerical nu_0 in the text; 30 GHz from a figure caption and BP XIII Table 2. |
| 11 | Fuskeland et al. 2014, 2021 | PDF v1 pp. 1-9; PDF v2 pp. 1-3, 7-13 | CONFIRMED WITH CORRECTION | 2014 values are the mean of T-T and ML. 2021 uses S-PASS and WMAP K only, no QUIJOTE. |
| 12 | QUIJOTE XII 2610.31469; de la Hoz et al. 2023, 2301.05117 | PDF pp. 1-3, 7-14; PDF pp. 1, 4-11, 16-19 | CONFIRMED WITH CORRECTION | XII: abstract and Sect. 6 give +-0.05, Table 4 gives +-0.03. 2301.05117 is paper VIII, not IV. |
| 13a | Neyman, Scott 1948 | BIB | CONFIRMED (bibliographic) | Page range not in the Crossref record. |
| 13b | Errors-in-variables statement | PDF of a 1985 review; ABS, BIB for the primaries | CONFIRMED WITH CORRECTION | Fuller 1987 and Kendall and Stuart: NOT FOUND (text not accessible). Statement verified in Engel 1985, which cites Patefield 1977, 1978. |
| 14 | Planck 2018 V, 1907.12875 | HTML-S, two renderings | CONFIRMED WITH CORRECTION | The ratio 1.54 +- 0.13 and the expected 1.96 to 1.98 are not in the paper. Only the two coefficients are. |

## 1. Ball, Fink, Bowler 2003, "Stochastic Annealing"

Identifier: arXiv cond-mat/0301179v1 (13 Jan 2003). Phys. Rev. Lett. 91, 030201, published
2003-07-14, DOI 10.1103/PhysRevLett.91.030201 (Crossref, BIB). The arXiv PDF carries no journal
reference; its arXiv comment reads "4 pages, 2 figures, submitting to PRL".

Read: PDF, all 5 pages of arXiv v1. The PRL text was not read, so the equation number in the
journal version is not verified.

Where: pp. 1 to 2, Eqs. (3) to (7) and Fig. 1; p. 4 for the definition of sigma and Eq. (20).

Quotes:
- p. 1: "We now consider the case when the true energy change is not known exactly, and we must
  accept moves with probabilty (sic) A(x) where x is only an estimate of the energy change. We
  will assume that these estimates have statistically independent errors."
- Eq. (6): K(Delta E) = int_{-inf}^{+inf} f(x | Delta E) A(x) dx.
- p. 1: "We can gain insight by considering the crude choice A(x) = 1 for x < 0 and A(x) = 0
  otherwise. This simple strategy can give a good aproximation (sic) to eq. (3) of great value in
  optimisation. The resulting acceptance function K when f is given by a Gaussian distribution is
  simply related to the standard Error Function and is graphed in Fig. 1."
- pp. 1 to 2, the decisive sentences: "The resemblance between this and the Glauber acceptance
  function (whose symmetry it shares) is striking, showing how the random energy errors make the
  selection look thermal - although in this case the match is not exact, so detailed balance is
  not strictly achieved. The standard deviation sigma of the Gaussian distribution controls an
  approximate effective temperature using this rule, as inverting eq. (3) gives

  T ~= sqrt(pi/8) sigma (1 - 0.018 (Delta E/T)^2 + order (Delta E/T)^4). (7)"
- p. 2: "For many optimisation purposes the departure from detailed balance, due to the energy
  dependent terms at large |Delta E/T|, is not a problem."
- Fig. 1 caption: "A simple approximate stochastic annealing is obtained by accepting moves on the
  basis of the sign of the estimated energy change x."
- p. 4: "sigma is the standard deviation of x".
- Eq. (3), which is inverted: K(mu -> nu)/K(nu -> mu) = e^{-beta Delta E}. Eq. (20), the exact rule
  for Gaussian errors: A(x) = min(1, e^{-beta (x + beta sigma^2/2)}).

Verdict: CONFIRMED.

Scope to state when citing: Eq. (7) is for the sign rule with Gaussian errors of standard deviation
sigma on the estimated energy difference, not for a Metropolis rule applied to a noisy difference
at a finite temperature. The paper says the match to a thermal rule is approximate.

(inference) Check of the constants: with K(D) = Phi(-D/sigma), ln[K(-D)/K(D)] =
2 sqrt(2/pi) (D/sigma) [1 + (2/(3 pi) - 1/6)(D/sigma)^2 + ...]. This gives T = sqrt(pi/8) sigma at
leading order and a correction coefficient (2/(3 pi) - 1/6)(pi/8) = 0.0179, which rounds to the
printed 0.018.

## 2. Andrieu and Roberts 2009, "The pseudo-marginal approach for efficient Monte Carlo computations"

Identifier: arXiv 0903.5480v1. Ann. Statist. 2009, Vol. 37, No. 2, 697-725, DOI 10.1214/07-AOS574
(printed on p. 1 of the PDF, which is the journal's electronic reprint: "This reprint differs from
the original in pagination and typographic detail").

Read: PDF, pp. 1 to 5, 18 to 20, 28 to 29 of 29. Pages 6 to 17 and 21 to 27 not read.

Where: Sect. 1 (pp. 2 to 5), Table 1 (p. 3), Sect. 6 (pp. 18 to 20), reference list (p. 29).

Quotes:
- Table 1 (p. 3), column MCWM: sample "theta* ~ q(theta, .)", "Z ~ q^N_theta(.)", "Z* ~
  q^N_theta*(.)"; compute both estimates pi~^N(theta) and pi~^N(theta*); accept with probability
  1 ^ r, r = pi~^N(theta*) q(theta*, theta) / [pi~^N(theta) q(theta, theta*)]. Column GIMH: only
  Z* is sampled, and Z is carried with the state.
- Definition (p. 3): "The first approach considered here, which corresponds to the middle column,
  is to attempt to approximate P, independently at each iteration using the importance ratio
  averages given in equation (1.1). More precisely both Z and Z* are "refreshed" at each iteration
  independently of previously sampled auxiliary variables given theta and theta*. This algorithm
  is referred to as the Monte Carlo within Metropolis (MCWM) in [1] following the terminology of
  [7]."
- Invariant distribution (p. 3): "However MCWM and the marginal algorithm P are not equivalent. In
  particular, pi(theta) is typically not the invariant distribution density of P~^MCWM_N and
  therefore will not produce samples from pi(theta) even in steady state. However, intuitively,
  provided that P~^MCWM_N is ergodic, the samples generated by this procedure will asymptotically
  be distributed according to an approximation of pi(theta), which should be all the more precise
  that N is large."
- GIMH (p. 4): "Rather, GIMH can be interpreted as a form of MCWM where Z is "recycled" from the
  previous iteration".
- Sect. 6, "Epilogue: MCWM" (p. 18): "However the existence of an invariant distribution for
  P~^noisy_N is not obvious in general (it is not a MH update). This is in contrast with
  P~^exact_N, for which the invariant distribution as well as its marginal distribution are known
  to be pi~^N(d theta, dZ) and pi(d theta), respectively."
- Theorem 9 (p. 19): "Assume (A1), (A4), (A5) and that for any N >= 1 there exists a probability
  distribution pi°^N on (Theta, B(Theta)) such that pi°^N P~_N = pi°^N, with P~_N = P~^noisy_N
  defined in equation (2.7)." Conclusions: ||P~^k_N(theta_0, .) - pi°^N(.)|| <= C~ rho~^k (6.1)
  and ||pi - pi°^N|| <= C epsilon / (1 - rho) (6.2), for N large enough.
- Reference list (p. 29): [1] is Beaumont, M. A. (2003), Genetics 164, 1139-1160. [7] is O'Neill,
  Balding, Becker, Eerola, Mollison (2000), Appl. Statist. 49, 517-542.

Verdict: CONFIRMED. MCWM is the scheme in which both estimates are refreshed at every iteration.
The paper says the target is typically not its invariant density, that the existence of an
invariant distribution is not obvious, and bounds the total-variation distance when one exists.
It gives no explicit form of the invariant law and uses no tempering language (on the pages read).

Note: p. 5 says that Theorem 9 shows that "generalizations of GIMH can inherit this property with
arbitrary precision". (inference) "GIMH" there is a slip for MCWM, since Theorem 9 is about
P~^noisy_N.

## 3. Bardenet, Doucet, Holmes, "On Markov chain Monte Carlo methods for tall data"

Identifier: arXiv 1505.02827v1 (11 May 2015). J. Mach. Learn. Res. 18(47):1-43, 2017 (JMLR page,
BIB). The arXiv record has no journal reference.

Read: PDF of arXiv v1, pp. 1 to 2 (title, contents) and 20 to 23 (end of Sect. 6.2.2 and all of
Sect. 6.3). Section numbers of the JMLR version were not checked.

Where: Sect. 6.3, "Exchanging acceptance noise for subsampling noise", pp. 22 to 23, Eqs. (22) to
(24). Related sentence in Sect. 6.2.2, p. 20.

Quotes:
- p. 22: "This section is an original contribution, which illustrates a way to obtain subsampling
  algorithms with guarantees under weaker assumptions than Gaussianity. This approach is
  impractical, but it is of methodological and illustrative interest."
- p. 22: "Why, then, not rely on the subsampling noise to guarantee exploration, and accept a move
  if and only if Lambda*_t(theta, theta') + (1/n) log[p(theta') q(theta|theta') / (p(theta)
  q(theta'|theta))] > 0 (22) instead of (21)? This idea has been first used by Branke et al.
  (2008) to develop heuristics for simulated annealing in the presence of noise."
- p. 23: "We do not assume that the Lambda*_t(theta, theta')'s are Gaussianly distributed, but we
  make the parametric assumption that the second and third absolute moments sigma^2 and rho of
  -log p(x*_i|theta') + log p(x*_i|theta) are known and independent of theta, theta'."
- p. 23: "Now let C, lambda > 0 be such that for any x in R, |Phi(x) - 1/(1 + e^{-lambda x})| <=
  C. Bowling et al. (2009) for instance, empirically found C = 0.0095 and lambda = 1.702."
- p. 23, the decisive sentence: "Hence, the acceptance probability of an algorithm that would
  accept the move from theta to theta' by checking whether Lambda*(theta, theta') > 0 is close to
  the acceptance probability of an MCMC algorithm with a Baker (sic) acceptance criterion (Robert
  & Casella, 2004, Section 7.8.1) that targets pi^beta with temperature beta = lambda sqrt(t) /
  (n sigma)."
- p. 23: "This reasoning shows again a close relation between subsampling and tempering, as in
  Section 6.1, with a clear link between the variance of the subsampled log likelihood ratios and
  the temperature."
- Sect. 6.2.2, p. 20, on Austerity MH in a Gaussian example: "In the Gaussian case, the posterior
  is rightly centered, but is slightly too wide. This is a tempering effect due to too small
  subsamples".

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: Bardenet et al. (Sect. 6.3) show that accepting on the sign of a subsampled
log-likelihood ratio has an acceptance probability close to that of a Barker-type chain targeting
pi^beta with beta = lambda sqrt(t)/(n sigma), using a Berry-Esseen bound and the logistic fit to
the normal distribution function with lambda = 1.702 and C = 0.0095 (Bowling et al. 2009). They
explicitly do not assume Gaussian noise. They call beta the "temperature"; it is the exponent of
the target. sigma is defined through the second absolute moment of the single-datum log-likelihood
ratio, and Lambda is the average over data.

(inference) The estimated total log ratio n Lambda*_t has standard deviation s = n sigma/sqrt(t),
so beta = lambda/s and the effective temperature in the sense of item 1 is s/1.702 = 0.588 s,
against sqrt(pi/8) s = 0.627 s from Ball et al.

## 4. Stompor, Leach, Stivoli, Baccigalupi 2009

Identifier: arXiv 0804.2645v2. MNRAS 392, 216 (2009), DOI 10.1111/j.1365-2966.2008.14023.x (arXiv
API record; the PDF header carries a placeholder).

Read: PDF of arXiv v2, pp. 1 to 6 and 15 to 17 of 18. Pages 7 to 14 (Sects. 3.2 and 4) not read.

Where: Sect. 3, Eqs. (5) and (8); Sect. 3.1.1, Eqs. (9) to (12), Fig. 2; Sect. 3.1.2; Sect. 5.

Equations as printed:
- Eq. (5): -2 ln L_data(s, beta) = CONST + (d - A s)^t N^-1 (d - A s).
- Eq. (8): s = (A^t N^-1 A)^-1 A^t N^-1 d.
- Eq. (9): -2 ln L_spec(beta) = CONST - (A^t N^-1 d)^t (A^t N^-1 A)^-1 (A^t N^-1 d).
- Eq. (10): -2 ln L_marg(beta) = -2 ln int ds exp[-(1/2)(d - A s)^t N^-1 (d - A s)]
  = CONST - (A^t N^-1 d)^t (A^t N^-1 A)^-1 (A^t N^-1 d) + ln|(A^t N^-1 A)^-1|.
- Eq. (12): -2 d ln L_marg / d beta_i at beta = beta_ML equals -2 tr[N^ A^t N^-1 A_{,beta_i}], with
  N^ = (A^t N^-1 A)^-1.

Quotes (Sect. 3.1.1, pp. 4 to 6):
- "We can use that equation to eliminate the sky signals from the full data likelihood expression,
  equation (5), obtaining," Eq. (9).
- "The above expression permits a calculation of the 'ridge' values of the full data likelihood,
  i.e., the maximum value of the likelihood as attained for any chosen value of the spectral
  parameters and, therefore, is maximized by the true (global) maximum likelihood value of the
  spectral parameters. Hereafter, we will call this function a spectral likelihood and use it to
  determine the maximum likelihood estimates of the spectral parameters and obtain an insight into
  the resulting errors. We point out that the spectral likelihood is not a proper likelihood in a
  common sense, and is thus often referred to as a 'ridge' or 'profile' likelihood to emphasize
  this fact (e.g., Barndoff-Nielsen & Cox 1994). Nevertheless, it permits to retrieve the maximum
  likelihood values of the spectral parameters, and provides useful, if not strictly precise,
  estimates of their uncertainty."
- "In the alternate approach as implemented by Eriksen et al. (2006) the spectral parameters are
  estimated with help of the marginal distribution derived for the full data likelihood, equation
  (5). The marginalisation over the sky amplitudes is performed numerically, via MCMC, and assumes
  flat priors for the sky signals. We will show here that that approach is not only more
  computationally involved but also fails to reproduce the maximum likelihood estimates of the
  spectral parameters derived by maximizing equation (5)."
- "We first observe that the marginalisation procedure can be done analytically and that the
  integration over the sky amplitudes leads to," Eq. (10), "where |...| denotes the matrix
  determinant."
- "(We remark here in passing that though this equation is derived under the assumption of the
  flat priors for the component amplitudes. However, any kind of a spatial template-like prior
  constraint with known, and Gaussian, uncertainty can be easily incorporated here, and treated as
  yet another data set in equation (5).)"
- "The marginalised and spectral likelihoods, equations (10) & (9), clearly differ. As the
  difference depends on the spectral parameters, we can conclude that the marginal likelihood
  peaks at their values, which do not coincide with their maximum likelihood estimates."
- Sect. 3.1.2 (p. 6): the comparison shows "the macroscopic gain in adopting the expression in
  equation (9) throughout the analysis."
- Sect. 5 (pp. 15 to 16): "on the first step, we perform a maximization of the spectral likelihood
  as given by formula equation (33); on the second step, we use the derived parameters to estimate
  the component maps via equation (34). We have shown that this procedure indeed leads to the
  maximum likelihood estimates derived directly from the full data likelihood (equation 5)." And:
  "We note that it is possible to treat the method described here just as an algorithm for an
  efficient maximization of the full data likelihood, equation (5), avoiding therefore a need of
  interpreting the expression in equation (33) as a common likelihood."
- App. A (p. 17): "We will interpret the inverse of the curvature matrix, computed at the peak of
  the likelihood, as a good approximation to the covariance matrix for the parameters estimated via
  maximization of the likelihood."

Verdict: CONFIRMED. They define the spectral (ridge, profile) likelihood in Eq. (9), give the
flat-prior amplitude-marginal with a determinant term in Eq. (10), and choose Eq. (9), on the
ground that the marginal does not peak at the maximum-likelihood values. Their criterion is the
location of the peak; uncertainties come from the curvature of the profile likelihood, which they
themselves call "useful, if not strictly precise".

Flag (inference, my derivation; check before quoting Eq. (10)): the Gaussian integral gives
int ds exp[-(1/2) s^t M s + b^t s] = (2 pi)^{n/2} |M|^{-1/2} exp[(1/2) b^t M^-1 b] with
M = A^t N^-1 A, hence -2 ln L_marg = const - b^t M^-1 b - ln|M^-1|. The arXiv v2 prints the
determinant term with a plus sign, and Eqs. (11) and (12) carry the same sign consistently. Alonso
et al. Eq. (A7) (item 7), p(b|d) proportional to sqrt(det N^_T) exp[...], has the sign of my
derivation. The MNRAS version was not checked. Paper A1 should write the determinant term from its
own derivation and not copy the sign from the arXiv version.

Not found on the pages read: the phrase "direct or via MCMC sampling" that
`2026-11-09_novelty_check_methods.md` attributes to Sect. 5. It may be on pp. 7 to 14.

## 5. Eriksen et al. 2008, "Joint Bayesian component separation and CMB power spectrum estimation"

Identifier: arXiv 0709.1058v2. ApJ 676, 10-32 (2008), DOI 10.1086/525277 (arXiv API record).

Read: PDF of arXiv v2, pp. 1 to 12 of 23.

Where: Sect. 3.1 (p. 4) for the prior adopted; Sect. 4.1 (p. 8), Eq. (30); Sect. 4.2 "Uniform vs.
Jeffreys' prior" (pp. 9 to 10), Eqs. (31) to (36), Figs. 3 and 4; Sect. 4.5 (p. 12).

Quotes:
- Sect. 3.1: "For this paper, we choose the prior most widely accepted in the statistical
  community, namely Jeffreys' ignorance prior (e.g., Box & Tiao 1992). This prior is given by the
  square root of the Fisher information measure, P_theta proportional to sqrt(|F|_theta theta) =
  sqrt(|-d^2 ln L / d^2 theta|). Its effect is essentially to "normalize" the parameter volume
  relative to the likelihood, and make the likelihood so-called "data translated"."
- Sect. 3.1: "For amplitude-type degrees of freedom, the ignorance prior works out be (sic) the
  usual flat prior, but for non-linear parameters, e.g., spectral indices, it is non-uniform. In
  particular, for the power law spectrum parametrized by a spectral index beta described above, it
  reads P(beta) proportional to sqrt(sum_nu ((G(nu)/sigma_nu)(nu/nu_0)^beta ln(nu/nu_0))^2)."
- Eq. (32): P_J(theta) ~ sqrt(F_theta theta) = sqrt(-<d^2 ln L / d^2 theta>).
- Eq. (36): P_J(beta) ~ sqrt(sum_nu ((1/sigma_nu)(nu/nu_1)^beta ln(nu/nu_1))^2), with Eqs. (34)
  and (35): P_J(A) ~ 1, P_J(m) ~ 1.
- Setup (Sect. 4.2): "Consider a full-sky data set pixelized at HEALPix resolution N_pix = 16
  (sic) (3072 independent pixels). Reduce the noise to 1 muK RMS per pixel. Use the same signal
  model as before, but with an offset common to all pixels," Eq. (31): T_nu(p) = m delta_{nu,nu_1}
  + A(p) (nu/nu_1)^{beta(p)}. The true offset is m = 0 (Sect. 4.1).
- The example: "Second, however, the distribution has a mean and standard deviation of -1.8 +-
  0.4, more than 4 sigma away from zero! Repeated experiments with different noise seeds gave
  similar results."
- "This behaviour is a result of the choice of prior. We initially adopted a uniform prior on the
  offset, the amplitudes and the spectral indices, with little thought to why we should do so.
  This was a poor choice."
- Prior volume: "Figure 4 shows Jeffreys' prior for the spectral index beta, limited to -4 <= beta
  <= -2. beta = -2 is given about two and half times more weight than beta = -4. Intuitively,
  this is necessary because there is an asymmetry between a steep and a shallow spectrum: A steep
  spectrum means that the signal dies off quickly with frequency, while a shallow spectrum implies
  that it maintains its strength longer. Thus, there is a larger allowed parameter volume with
  steep indices than with shallow, leading to an imbalance in terms of marginal probabilities.
  This parametrization effect is countered by the Jeffreys' prior."
- "The conclusion is clear: a proper ignorance prior leads to unbiased estimates, while a naive
  uniform prior leads to biased estimates."
- "In the rest of the paper, when we say that a Gaussian prior is adopted for the spectral
  indices, we mean a product of a Gaussian and Jeffreys' prior."
- Sect. 4.5: "The uniform prior is a proper ignorance prior for Gaussian variables only. In
  general, Jeffreys' rule should be used in the absence of informative priors."

Verdict: CONFIRMED. The quantity biased by more than 4 sigma is the global offset m (-1.8 +- 0.4
against a true value of 0), under uniform priors on offset, amplitudes and indices. The unit of
-1.8 +- 0.4 is not stated in the sentence (the noise is 1 muK per pixel and Fig. 3 has an unlabeled
unit). The paper shows the Jeffreys result only as a curve in Fig. 3, without a number in the text.

## 6. Dunkley et al. 2009, WMAP five-year Bayesian polarization maps

Identifier: arXiv 0811.4280v2. ApJ 701, 1804-1813 (2009), DOI 10.1088/0004-637X/701/2/1804 (arXiv
API record).

Read: PDF of arXiv v2, all 11 pages.

Where and quotes:
- Resolution, Sect. 2.1: "In this analysis we will use HEALPix Nside=8 with Np = 768." Sect. 4:
  "We apply the sampling method to the five-year WMAP data, using low resolution coadded Nside=16
  maps and inverse noise matrices as the inputs, and then degrading to Nside=8 as for the
  simulations."
- Index pixels, Sect. 2.1: "Third, we define Ni < Np pixels within which beta_2 takes a common
  value, rather than allow it to take a unique value at each Nside=8 pixel." and "In this case we
  use Nside=2 HEALPix pixels (Ni = 48)." Also: "First, that the spectral indices in Q and U are
  the same in a given pixel".
- Bands, Sect. 3: "We perform the Gibbs sampling using the four WMAP frequency bands K, Ka, Q, V
  (hereafter KKaQV). We do not include W band as standard". Synchrotron pivot 23 GHz, dust pivot
  94 GHz, dust index fixed at 1.7 (Sect. 2.1, Eqs. 8 and 9).
- Sampler, Sect. 2.2.2: index drawn from p(beta|A, d) by Metropolis-Hastings, Eqs. (19) and (20):
  -2 ln p(beta|A, d) = chi^2 - 2 ln p(beta).
- Prior, Sect. 2.2.2: "We choose a Gaussian prior of -3.0 +- 0.3, motivated by understanding of the
  synchrotron emission (Rybicki & Lightman 1979) and allowing for variations of the size observed
  in the synchrotron intensity (e.g., Bennett et al. (2003)). This is combined with a uniform
  prior on the CMB and synchrotron amplitudes, and a Gaussian prior 0 +- 0.2 I_d on the dust Q and
  U amplitudes."
- Volume effect, Sect. 2.2.2: "In the limit of low signal-to-noise, there is a larger volume of
  models with a flatter (i.e. beta tending to 0) synchrotron spectrum, allowing large CMB and
  synchrotron amplitudes of opposite sign. One approach is to modify the prior on the spectral
  indices to include an additional 'phase-space' factor."
- Sect. 3.1: "The resulting marginalized distributions for the Ni = 48 spectral indices are close
  to Gaussian, but centered on beta = -2.6." and "This is not a significant effect for individual
  pixels, but the phase-space volume relative to that of the prior central value increases for
  more pixels sharing a common index parameter."
- Jeffreys, Sect. 3.1: "Imposing the Jeffreys' prior on the spectral index parameters in this model
  would give less weight to steep indices than shallow ones, accounting for their smaller effect
  on the likelihood. However, this would not overcome the phase space effect that arises from the
  increase in volume of models with shallow indices."
- What they did, Sect. 3.1: "While there may be alternative parameterizations that avoid this
  problem, we choose instead to modify the spectral index prior. We use p(beta) = p_s(beta) /
  L_0(beta), where p_s is the Gaussian prior -3.0 +- 0.3, and 1/L_0 is an additional 'phase-space'
  factor. We estimate this factor by evaluating the marginalized posterior, p(beta|d), for a
  noise-only simulation with prior p_s(beta), and defining L_0(beta) = p(beta|d)/p_s(beta)." The
  factors 1/L_0 "are modeled as Gaussian distributions", "with beta = -3.7 +- 0.5 for the Ni = 48
  adopted in this analysis."
- Their own caveat, Sect. 3.1: "This is not the optimal solution, as the prior is informative in
  the presence of a signal and could therefore lead to parameter bias, in particular in the weak
  signal regime."
- Mean index, Sect. 4.2: "In regions of low synchrotron signal-to-noise the index is driven by the
  prior, so we mask pixels with sigma(beta) > 0.25. There is a preference in the North Polar Spur
  and Fan region for an index of ~ -3: the index averaged over regions with sigma(beta) < 0.25 is
  -3.03 +- 0.04."
- Latitude, Sect. 4.2: "Cutting the sky into high and low latitude at b = 20 degrees, the low
  latitude weighted mean, -3.00 +- 0.04, is consistent with the high latitude -3.08 +- 0.06. This
  contrasts with the more significant steepening of the index with latitude observed in Kogut et
  al. (2007) in the analysis of the three-year WMAP polarization data, although here we exclude
  part of the plane."
- Prior test, Sect. 4.2: with p_s = -2.7 +- 0.3, "The sigma(beta) < 0.25 pixel average in this
  case is -2.94 +- 0.04. The estimated CMB power is little affected by this change in prior, with
  tau = 0.092 +- 0.019."
- tau, Sect. 4.1: "Computing the likelihood as a function of tau we find tau = 0.090 +- 0.019 for
  the Gibbs-sampled maps outside the P06 mask, which is consistent with the results obtained
  through template cleaning, which give tau = 0.086 +- 0.016 for the KaQV data combination." With
  W band: "tau = 0.085 +- 0.017."
- Simulations, Sect. 3.2: input tau = 0.1 gives "tau = 0.091 +- 0.019 outside the P06 Galactic
  mask"; and "A larger increase of ~1 sigma in tau, to 0.111 +- 0.019, is seen when the spectral
  index of the simulation shallows from -3.0 to -2.5 between K and W band, but is modeled by a
  pure power law."
- Abstract: "This estimate does not vary significantly with Galactic latitude, although includes
  an informative prior."

Verdict: CONFIRMED on every point: Nside 8 maps (768 pixels, from Nside 16 inputs), 48 index
pixels at Nside 2, Gaussian index prior -3.0 +- 0.3 multiplied by a phase-space factor, mean index
-3.03 +- 0.04 over pixels with sigma(beta) < 0.25, -3.00 +- 0.04 at low and -3.08 +- 0.06 at high
latitude (cut at b = 20 degrees), tau = 0.090 +- 0.019.

Resolved from the repo notes: Sect. 3.1 prints both "increased to 96, and 192" and "Ni=48, 92
(sic) and 192"; the legend of Fig. 1 has N=96. The "92" is a misprint in the paper.

## 7. Alonso, Dunkley, Naess, Thorne 2017 (BFoRe)

Identifier: arXiv 1608.00551v1. Phys. Rev. D 95, 043504 (2017), DOI 10.1103/PhysRevD.95.043504
(arXiv API record).

Read: PDF of arXiv v1, pp. 4 to 7 and 15 to 19 of 19.

Where: Sect. II C 2 "Priors" (p. 7); App. A 2, Eqs. (A6) and (A7) (p. 16); App. B "Volume prior
for spectral parameters", Eqs. (B1) to (B4), Fig. 11 (pp. 17 to 19).

Quotes:
- Sect. II C 2: "In this work we imposed loose Gaussian priors on the spectral parameters, with
  beta_s = -3 +- 0.5, beta_d = 1.54 +- 0.5 and Theta_d = 20 +- 0.5 K." and "Besides these, we
  included a "volume prior" designed to take into account the volume of likelihood space for
  non-linear parameters. This is described in detail in Appendix B, and is equivalent to the
  widely used Jeffreys prior for the spectral parameters [46]."
- Eq. (A7): p(b|d) proportional to sqrt(det(N^_T)) exp[(1/2) T-bar^T N^_T^-1 T-bar] p_p(b), with
  N^_T^-1 = F^T N^-1 F (Eq. A4).
- App. B: "Thus, in this regime, the posterior for the spectral parameters b is dominated by the
  prefactor sqrt(|N^_T|). In the simplified case where the only component is synchrotron, and in
  the absence of polarized channels, this volume factor is given by: sqrt(|N^_T|) proportional to
  [sum_nu (nu/nu_0^s)^{2 beta_s} (1/sigma_nu^2)]^{-1/2}, (B2) which becomes arbitrarily large for
  beta_s -> -infinity."
- App. B text (p. 18): "The presence of the volume factor biases the estimate of the spectral
  indices by more than 4 sigma, and their true value is recovered within 1 sigma after accounting
  for it. The situation in a low-S/N region is ilustrated (sic) in the lower panels of Fig. 11:
  even using a broad Gaussian prior with sigma_beta ~ 1. centered on the true values of the
  spectral indices, the mean estimated indices are hundreds of sigmas away from their true
  values."
- Fig. 11 caption: "The posterior mean is more than 3-sigmas away from the true input value.
  Upper right panel: same plot after cancelling the volume prior. The true value lies within 2
  sigma of the posterior mean."
- App. B: "These two undesirable features (i.e. the non-flat posterior for the spectral indices in
  the signal-less case, and the bias it entails) can be avoided by including a factor
  |N^_T|^{-1/2} in our prior that exactly cancels the volume factor. This problem can also be
  solved by applying a Jeffreys prior, given by the square root of the Fisher information matrix
  [46]", Eq. (B3), and for a power law Eq. (B4): p_Jef(beta_s) proportional to sqrt(sum_nu
  [(nu/nu_0^s)^{beta_s} (1/sigma_nu) ln(nu/nu_0^s)]^2), "which is equivalent to the volume prior
  p_Vol = |N^_T|^{-1/2} defined by Eq. B2 except for the subdominant logarithmic term."

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: Alonso et al. apply a volume prior |N^_T|^{-1/2} that cancels the determinant
factor of the amplitude-marginal posterior (App. B). Without it, the high-S/N example is biased by
"more than 4 sigma" according to the text and by "more than 3-sigmas" according to the caption of
Fig. 11, and the low-S/N example by "hundreds of" sigma. Cite as "more than 3 sigma" or quote both.
The examples are one simulated region each, shown for the dust and synchrotron indices jointly.

## 8. Weiland et al. 2022, "Polarized Synchrotron Foreground Assessment for CMB Experiments"

Identifier: arXiv 2203.11445v2. ApJ 936, 24 (2022), DOI 10.3847/1538-4357/ac83ab (Crossref, BIB).

Read: PDF of arXiv v2, pp. 1 to 19 of 21 (the remainder is the reference list).

Where: abstract; Sect. 3 (pp. 3 to 6), Figs. 2 to 4, Table 1; Sect. 4.2.1 (p. 8); Sect. 4.3
(pp. 12 to 13); Sect. 5 and Table 2 (pp. 14 to 15); Sect. 7 (p. 16); Sect. 8 (pp. 17 to 19).

Quotes:
- Method for plane and spurs, Sect. 3: "We use a linear correlative analysis method similar to
  that described by Fuskeland et al. (2014) for each data frequency pair. Stokes Q and U maps at
  each frequency are first smoothed to 2 degrees FWHM, and degraded to HEALPix Nside = 32 pixels,
  or about 1.8 degrees on a side. We then fit a line of zero intercept through the combined Q and
  U data points within a ~ 7 degree superpixel (HEALPix Nside = 8)." The index follows from the
  slope m as beta^pol = log(m)/log(nu_y/nu_x). The pairs are WMAP 23 and 33 GHz, and WMAP 23 with
  Planck PR3 LFI 30 GHz.
- Method off the plane, Sect. 4.2.1: "We work at HEALPix Nside = 16, or roughly a pixel resolution
  of 3.7 degrees." A three-frequency model, Eqs. (2) to (4), is fit per pixel to [1.4, 23, 33] or
  [2.3, 23, 33] GHz with free RM, Q_0, U_0 and beta^pol.
- Table 1, "beta^pol for Plane and Spur Regions", columns [23, 33] and [23, 30]:
  region 1: -2.90 +- 0.020 and -3.27 +- 0.022; region 2: -3.05 +- 0.010 and -3.01 +- 0.011;
  region 3: -3.11 +- 0.036 and -3.09 +- 0.037; region 4: -3.04 +- 0.062 and -3.15 +- 0.059.
  Regions (Sect. 3): "(1) an outer plane region dominated by the Fan Region, (2) an inner plane
  dominated region, and northern (3) and southern (4) extensions off the plane."
- Inconsistency, Fig. 2 caption: "those in the Fan Region (on the left at l ~ 140 degrees, where
  the colors shift from mostly orange to mostly dark blue) are statistically significant, but are
  attributable to large-scale systematic differences between WMAP and Planck LFI maps. As
  discussed by Planck Collaboration II (2020), the large-scale systematic modes are primarily
  attributable to the 30 GHz map."
- Sect. 3: "There is no likely physical mechanism that suggests that the region 1 results for
  [23,30] and [23,33] are both correct." and "A non-physical origin for the Region 1 discrepancy
  is more compelling."
- Sect. 3: "Given that the polarized spectral index results depend significantly on processing
  differences of the same 30 GHz data set, we favor the use of the more stable and reliable 33 GHz
  data. We adopt the beta^pol determined from the [23,33] frequency pair as our best estimate for
  these high SNR regions."
- Sect. 8: "The discrepancy traces to large-scale systematic differences between the WMAP 23 GHz
  and LFI 30 GHz maps, which have been shown (Planck Collaboration Int. LVII 2020) to change as a
  function of the mapmaking algorithm used to process the 30 GHz data. For this reason, we
  sacrifice the SNR advantages that would be gained from a combination of WMAP and Planck LFI
  data, and our results are based on the WMAP [23, 33] GHz maps only."
- High latitude, abstract: "We find a clear spectral index steepening with increasing latitude
  south of the Galactic plane with Delta beta^pol = 0.4, and a smaller steepening of 0.25 in the
  north. Near the south Galactic pole the polarized synchrotron spectral index is beta^pol ~
  -3.4." and "Within the BICEP2/Keck survey footprint, we find consistency with a constant value,
  beta^pol = -3.25 +- 0.04 (statistical) +- 0.02 (systematic)." Also: "Data quality masking leaves
  44% sky coverage (73% for |b| > 45 degrees)."
- Table 2, latitudinal profile (bin centre, beta^pol): 80: -3.238 +- 0.011; 60: -3.217 +- 0.007;
  40: -3.217 +- 0.012; 20: -3.088 +- 0.016; 5: -2.964 +- 0.014; -5: -2.960 +- 0.013; -20: -3.166
  +- 0.026; -40: -3.198 +- 0.010; -60: -3.245 +- 0.007; -80: -3.378 +- 0.010. Sect. 5 warns that
  the profile "is likely biased" at mid-latitudes, toward larger values.
- Loss imbalance, Sect. 4.3: modes "could possibly produce biases in the recovered beta^pol of up
  to Delta beta ~ 0.05"; "we do not have sufficient evidence of a systematic bias in beta^pol
  resulting from loss imbalance uncertainty, and do not include it in the estimated
  uncertainties."
- Sect. 8, conclusion 2: "Current public data are insufficient to characterize beta^pol for the
  future, for three reasons: a lack of full sky coverage; conflicting results between experiments;
  and a lack of sensitivity to support future B-mode experiments."

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: the zero-intercept Q/U correlation in Nside 8 superpixels is the method for
the Galactic plane and spur regions only. The high-latitude results come from per-pixel fits at
Nside 16 that add 1.4 or 2.3 GHz and a rotation measure, and the composite map uses WMAP 23 and
33 GHz with no Planck 30 GHz. The Fan-region numbers are -2.90 +- 0.020 from [23, 33] GHz and
-3.27 +- 0.022 from [23, 30] GHz (Table 1, region 1, statistical errors, PR3 30 GHz). The paper
also reports that the PR4 and BeyondPlanck 30 GHz maps change the plane results (Fig. 4).

(inference) The Fan difference is 0.37, about 12 times the quadrature sum of the two quoted
errors; the two estimates share the 23 GHz map.

## 9. Weiland et al. 2018, tau and template uncertainties

Identifier: arXiv 1801.01226v3. ApJ 863, 161 (2018), DOI 10.3847/1538-4357/aad18b (Crossref, BIB).

Read: PDF of arXiv v3, pp. 1 to 10 of 13.

Where: Table 2 "Derived values for tau, EE only" (p. 5); Sect. 4 "An unaddressed systematic"
(pp. 4 to 6); Sect. 4.3 (p. 9); Sect. 5 (p. 10); abstract.

Table 2 (CMB data; synchrotron template; dust template; coefficients; f_sky; tau; chi^2/dof):
- row 1: LFI 70 GHz; LFI 30 GHz; plc; fixed; 46%; 0.052 +- 0.023; no chi^2; the Planck value.
- row 2: LFI 70 GHz; LFI 30 GHz; plc; fixed; 46%; 0.053 +- 0.023; 728/734.
- row 3: LFI 70 GHz; LFI 30 GHz; default; free; 46%; 0.076 +- 0.020; 823/734.
- row 4: WMAP KaQV; WMAP K; plc; free; 46%; 0.045 +- 0.019; 754/734.
- row 5: WMAP KaQV; WMAP K; plc; fixed; 46%; 0.040 +- 0.017; 752/734.
- row 6: WMAP KaQV; WMAP K; default; fixed; 75%; 0.058 +- 0.017; 1190/1162.
- row 7: LFI 70 GHz; WMAP K; plc; free; 46%; 0.094 +- 0.016; 846/734.
- row 8: WMAP KaQV; LFI 30 GHz; plc; free; 46%; 0.37 +- 0.02; 933/734.

Quotes:
- Sect. 4: "We next perform foreground cleaning of WMAP KaQV data using the Planck 30 GHz
  synchrotron template, and conversely, clean LFI 70 GHz using the WMAP K-band synchrotron
  template, holding all other analysis assumptions constant. The last two rows of Table 2 show the
  resultant values for tau when using these alternate choices of the synchrotron template. In both
  cases, chi^2 values are significantly increased, as are the derived tau values. The 70 GHz
  result returns tau = 0.094 +- 0.016, in tension with the HFI results listed in Table 1. WMAP
  KaQV cleaned using the LFI 30 GHz synchrotron template returned tau near 0.37".
- Sect. 4: "that the magnitude of the residual is not consistent with spectral index variations
  over the sky".
- Sect. 4.3, with a mask that keeps 11 percent of the sky: "returned tau for WMAP cleaned with LFI
  30 GHz of 0.24 +- 0.03, and 0.09 +- 0.03 for LFI 70 GHz cleaned with WMAP K-band."
- Sect. 5, issue 3: "An unexpectedly high value for tau is derived when WMAP K-band Stokes Q and U
  maps are used to clean synchrotron emission from polarized Planck LFI 70 GHz data, and a
  severely compromised tau value results when LFI 30 GHz Q and U maps are used to clean
  synchrotron from the combined WMAP Ka, Q and V bands."
- Sect. 5, attribution: "We reject both of these assertions and show evidence that the signature
  is a systematic internal to LFI and roughly common-mode to all LFI bands. The spectral behavior
  of the systematic scales like the intensity signal and is thus likely to arise from
  temperature-to-polarization leakage in the LFI channels."
- Abstract: "We demonstrate significant bias in the optical depth derived when using the LFI 30
  GHz map as a template to clean synchrotron from WMAP data, and recommend against the use of the
  2015 LFI 30 GHz polarization data as a foreground template for non-LFI datasets." and "all
  values are consistent with the range of tau = 0.07 +- 0.02."

Verdict: CONFIRMED. LFI 70 GHz cleaned with WMAP K gives 0.094 +- 0.016; cleaned with LFI 30 GHz it
gives 0.053 +- 0.023 (plc dust template, fixed coefficients) or 0.076 +- 0.020 (default dust
template, free coefficients). Say when citing: these are 2015 Planck release data, the fits are EE
only at Nside 8, and no 30 GHz row has exactly the settings of row 7 (plc dust, free
coefficients). The authors attribute the effect to an LFI systematic, not to the spectral index.

## 10. Svalheim et al., polarized foregrounds between 30 and 70 GHz (2011.08503)

Identifier and title: arXiv 2011.08503. The title depends on the version (arXiv API, Crossref,
INSPIRE; BIB):
- arXiv v1 (2020-11-17): "BeyondPlanck XIV. Polarized foreground emission between 30 and 70GHz".
- arXiv v2 (2022-01-18): "BeyondPlanck XV. Polarized foreground emission between 30 and 70 GHz".
- Published (Crossref and the EDP Sciences entry in INSPIRE): "BEYONDPLANCK XIV. Polarized
  foreground emission between 30 and 70 GHz", A&A 675, A14 (2023), DOI
  10.1051/0004-6361/202243160, published online 2023-06-28.

Read: HTML-S. The arXiv PDF (both versions) exceeds the fetch size limit and the A&A pages return
403, so no PDF page was read. Two renderings of arXiv v2 were queried independently with different
prompts: ar5iv (123 189 characters, both halves) and arxiv.org/html/2011.08503v2 (127 923
characters, both halves). The abstract of the published version was read from Crossref. Every
quote below is a short piece returned by the summarizer; where both renderings returned the same
passage they agree.

Where and quotes:
- Bands, Sect. 3 and Table 1: LFI 30, 44 and 70 GHz (time-ordered data), WMAP Ka, Q and V (Nside
  16, full covariance), HFI 353 GHz in polarization. "we include only the WMAP Ka-, Q-, and V-band
  channels in the following, centered on 33, 41 and 61 GHz, respectively." and "the 353 GHz
  channel in polarization, to constrain the amplitude of thermal dust emission." "The W-band
  channel is omitted because it is known to be contaminated by systematic residuals".
- K band, footnote 2 (placed in Sect. 2.1.1 by one rendering and in Sect. 3 by the other): "Note
  that the BeyondPlanck data combination does not include the WMAP K-band sky, as this channel
  would otherwise compete with Planck LFI 30 GHz" (the summarizer's piece ends there). Sect. 3:
  "its statistical sensitivity is so high that it would rival that of the Planck LFI 30 GHz
  channel."
- Index likelihood, Sect. 4.3, Eq. (26): P(beta | d, a) proportional to P(d | a, beta) P(beta),
  where "d={d_nu} is the set of all available frequency maps"; "we instead adopt a standard
  Metropolis-within-Gibbs sampler for beta." The index is sampled conditional on the amplitudes.
  "we first propose beta_s,i independently for each region i, and then smooth the resulting map"
  "with a 10 degree FWHM Gaussian beam to suppress edge effects."
- Regions, Sect. 2.1.1: "We therefore partition the sky into four large disjoint regions, as shown
  in the top panel of Fig. 1"; "We started from the 24-region partitioning defined by Fuskeland et
  al. (2014)"; "If this was found to be strongly prior-dominated, we merged neighboring regions.";
  "This process resulted in four main regions, which we will refer to as 1) High Latitudes;" "2)
  the Galactic Spur or simply Spur; 3) the Galactic Center; and 4) the Galactic plane."
- Prior, Sect. 4.4: "we adopt for the main analysis Gaussian priors of beta_s=-3.3+-0.1 and
  beta_d=1.56+-0.10, respectively." "While the Spur and Galactic plane regions are fitted using
  both the likelihood and prior as described above," "we only include the prior for the
  high-Galactic and Galactic Center regions." "This is because of very strong degeneracies with
  respect to instrumental systematic effects in these regions." "while the High Latitude region is
  particularly sensitive to relative gain uncertainties" "and degeneracies with the CMB
  quadrupole".
- Values, Sect. 5.3.2: "For these, we find posterior mean and standard deviations of
  beta_s^Spur=-3.17+-0.06" "and beta_s^Plane=-3.03+-0.07." "while we find a variation of Delta
  beta_s=0.14 between the Galactic plane and the Spur" "The steepening is however only
  statistically significant at the 2 sigma level" "as determined in the current analysis." The
  same two values are in the abstract of arXiv v2, in Sect. 6 item 1, and in the abstract of the
  published version (Crossref).
- Mixing, Sect. 5.3.2: "Overall, we see that the correlation length is modest, as in 30-50
  samples." Two Gibbs chains, "both after discarding 10 samples for burn-in."
- Wider prior, Sect. 4.4: "Consequently, a broader prior width of sigma_beta=0.2, which permits
  synchrotron indices as shallow as beta_s=-2.7" "leads to large foreground excesses that only can
  be accommodated by the CMB component within the current model." "The result is an obvious
  foreground bias in the large-scale CMB extraction" "and a corresponding excess in estimates of
  the optical depth of reionization". No number is given.
- Reference frequency: the data model (Sect. 2.1, Eq. 6) has a_s (nu/nu_0,s)^beta_s with "nu_0,c
  is the reference frequency for component c", and no numerical value for nu_0,s was found in the
  text. The caption of Fig. 6 reads: "Same as Fig. 5, but for polarized synchrotron emission
  evaluated at 30 GHz."
- Companion papers. 2201.08188 (BeyondPlanck XIII, HTML-S, two renderings, first 100 000
  characters of each): Table 2, "Summary of main parametric signal models for the temperature
  analysis.", lists for synchrotron nu_0,s = 30 GHz, prior beta_s ~ N(-3.3, 0.1^2), and C = 0. The
  paper is about intensity and says polarization is treated in the Svalheim et al. paper.
  2011.05609 (BeyondPlanck I, HTML-S, ar5iv characters 0 to 300 000 of 363 959): Sect. 3.5 Eq.
  (28) has a_s (nu/nu_0,s)^beta_s and "nu_0,i is the reference frequency for component i"; no
  numerical synchrotron reference frequency was found in the part read. The last 64 000
  characters were not read.

Verdict: CONFIRMED WITH CORRECTION.

Corrections:
1. The published paper is "BeyondPlanck XIV" (A&A 675, A14, 2023). "XV" is the title of arXiv v2
   only. The reference list of Paper A1 should use the published number and say which version the
   numbers come from.
2. Spur -3.17 +- 0.06 and Plane -3.03 +- 0.07 are the values of arXiv v2 and of the published
   abstract. The abstract of arXiv v1 gives -3.15 +- 0.07 and -3.12 +- 0.06 (arXiv API, ABS).
3. The reference frequency of 30 GHz for the polarized synchrotron amplitude is supported by the
   Fig. 6 caption and by Table 2 of BeyondPlanck XIII, not by an explicit statement of nu_0,s in
   the text of this paper.
4. Read level: the paper was not read as a PDF. The first author should check the quotes against
   his own copy before they go into print.

Side observation (HTML-S, outside the item): in arXiv v2 the abstract gives beta_d = 1.64 +- 0.03
and Sect. 6 item 2 gives "beta_d=1.62+-0.04". The published abstract has 1.64 +- 0.03.

## 11. Fuskeland et al. 2014 and 2021

### 11a. Fuskeland, Wehus, Eriksen, Naess 2014

Identifier: arXiv 1404.5323. ApJ 790, 104 (2014), DOI 10.1088/0004-637X/790/2/104 (arXiv API).

Read: PDF of arXiv v1, pp. 1 to 9 of 11 (the v2 PDF exceeds the fetch size limit). The v2 abstract
was read from the arXiv abstract page (HTML-S).

Where: abstract; Sect. 2.1 to 2.3 (methods); Sect. 3 (data); Sect. 5, Table 1, Figs. 6 to 8.

Quotes (v1 unless stated):
- Data, Sect. 3: "We therefore focus on the two lowest frequencies, the K- and Ka-bands with
  effective frequencies of 22.45 and 32.64 GHz, respectively, for a synchrotron spectral index of
  beta = -3". Maps smoothed to 1 degree FWHM at Nside = 64; 24 regions.
- Methods: T-T plot with the effective variance method (Sect. 2.1) and a maximum-likelihood fit
  with offsets (Sect. 2.2), both averaged over polarization angles 0 to 85 degrees in steps of 5
  (Sect. 2.3, Eq. 14).
- Combination, Sect. 5: "we adopt the straight mean of the indices derived with the two methods as
  our final point estimate of the spectral index", with the difference between the methods added
  as an "algorithmic uncertainty".
- Result, Sect. 5: "Adopting a weighted average over Galactic and high-latitude regions, we find
  mean spectral indices of beta_plane = -2.98 +- 0.01 and beta_high-lat = -3.12 +- 0.04,
  respectively; the full-sky weighted mean is beta^all-sky = -2.99 +- 0.01, being strongly
  dominated by the Galactic plane regions."
- Longitude, Sect. 5: best-fit "(c_2, a) = (-2.85, 0.17)" for beta_2(l) = c_2 + a sin(2l - 90
  degrees).
- Abstract of v1: "we find that the synchrotron spectral index steepens by ~ 0.15 from low to high
  latitudes". Abstract of v2 (HTML-S): "Averaging over polarization angles, we derive a mean
  spectral index of beta_all-sky=-2.99+-0.01" "in the frequency range of 23-33 GHz." "We find that
  the synchrotron spectral index steepens by 0.14 from low to high Galactic latitudes," "with mean
  spectral indices of beta_plane=-2.98+-0.01 and" "beta_high-lat=-3.12+-0.04."

### 11b. Fuskeland et al. 2021

Identifier: arXiv 1909.05923v2. A&A 646, A69 (2021), DOI 10.1051/0004-6361/202037629 (arXiv API).

Read: PDF of arXiv v2, pp. 1 to 3 and 7 to 13 of 16.

Where: abstract; Sect. 2 (data); Sect. 4.1 and 4.2, Figs. 10 to 12; Sect. 4.3.

Quotes:
- Abstract: "We constrain the spectral index of polarized synchrotron emission, beta_s, by
  correlating the recently released 2.3 GHz S-Band Polarization All Sky Survey (S-PASS) data with
  the 23 GHz 9-year Wilkinson Microwave Anisotropy Probe (WMAP) sky maps. We subdivide the S-PASS
  field, which covers the southern ecliptic (sic) hemisphere, into 95 15 x 15 degree regions and
  estimate the spectral index of polarized synchrotron emission within each region using a simple
  but robust T-T plot technique."
- Abstract: "Adopting this correction model, we find that the mean spectral index of polarized
  synchrotron emission gradually steepens from beta_s ~ -2.8 at low Galactic latitudes to beta_s ~
  -3.3 at high Galactic latitudes, in good agreement with previously published results. The flat
  spectral index at the low Galactic latitudes is likely partly due to depolarization effects."
- Abstract: "A comparison with a similar analysis performed in the 23-33 GHz range suggests a
  flattening of about Delta beta_s ~ 0.1 +- 0.2 from low to higher frequencies, but with no
  statistical significance due to high uncertainties."
- Sect. 4.2: "Inverse-variance weighting of the estimates for all regions yields a full-sky
  average of beta_s = -3.24 +- 0.01 with Faraday corrections and beta_s = -3.32 +- 0.02 without
  Faraday correction. The corresponding standard deviations are 0.09 and 0.19, respectively."
  With correction 65 regions pass the correlation cut. Without correction the text reads 29
  regions and the caption of Fig. 10 reads 19 in the page images; one of the two may be a misprint
  or my misreading at this image resolution.
- Sect. 4.2: "Based on these measurements, we find that the spectral index typically range between
  beta_s = -2.7 and -2.9 along the Galactic plane, and between beta_s = -3.1 and -3.3 at high
  Galactic latitudes."
- Sect. 4.3: "the mean spectral indices are beta_BICEP2 = -3.22 +- 0.06 and beta_SPIDER = -3.21 +-
  0.03."
- Data, Sect. 2: maps at 1 degree FWHM and Nside = 64; effective frequencies 2.303 and 22.45 GHz.

Verdict: CONFIRMED WITH CORRECTION.

Corrected statements:
- 2014: the quoted indices (all sky -2.99 +- 0.01, plane -2.98 +- 0.01, high latitude -3.12 +-
  0.04, between WMAP K and Ka) are the straight mean of the T-T plot and maximum-likelihood
  estimates, not T-T plot values alone. The steepening is 0.14 in the published abstract and
  "~ 0.15" in arXiv v1.
- 2021: the data are S-PASS 2.3 GHz and WMAP 23 GHz only. QUIJOTE is not among the data sets on
  the pages read. The variation is from about -2.8 at low latitude to about -3.3 at high latitude
  (abstract), and the authors attribute the flat low-latitude value partly to depolarization.

## 12. QUIJOTE XII (2610.31469) and QUIJOTE VIII (de la Hoz et al. 2023, 2301.05117)

### 12a. Vansyngel et al., QUIJOTE XII

Identifier: arXiv 2610.31469v1 (2026-10-22), "QUIJOTE scientific results - XII. Characterization of
the polarized synchrotron emission at the power spectrum level in the QUIJOTE-MFI wide survey",
F. Vansyngel, J. A. Rubino-Martin, et al. The PDF header says "Accepted 2026 October 21" for
MNRAS. No journal reference or DOI in the arXiv record.

Read: PDF, pp. 1 to 3 and 7 to 14 of 18. Table 4 was cross-checked in arxiv.org/html (HTML-S).
Appendices not read.

Where: abstract; Sect. 2.1.3 (data sets); Sect. 3 (priors, p. 7); Sect. 4.2 (p. 8); Table 4
(p. 11); Sect. 5.2.3 (p. 12); Sect. 6 (p. 14).

Quotes:
- Abstract: "The synchrotron frequency spectral indices are found to be beta_s^EE = -3.26 +- 0.05
  and beta_s^BB = -3.25 +- 0.10 for the same cut, in agreement with earlier real-space analyses.
  We find no evidence for spatial variations of beta_s across different Galactic latitude cuts,
  within current uncertainties." The cut is |b| > 20 degrees and the multipole range is 30 <= l <
  130.
- Sect. 6: "The best-fit synchrotron spectral index is found to be beta_s^EE = -3.26 +- 0.05 and
  beta_s^BB = -3.25 +- 0.10 for the |b| > 20 degree. These values are consistent with previous
  analyses based on real-space statistics using QUIJOTE-MFI data in combination with Planck and
  WMAP (de la Hoz et al. 2023; Adak et al. 2026)."
- Table 4, column beta_s (EE, BB), read from the PDF and from the HTML with the same result:
  mask +-10: WM+Pl -3.00 +- 0.09, -2.83 +- 0.20; WM+Pl+QJ -3.22 +- 0.03, -3.33 +- 0.07.
  mask +-15: WM+Pl -3.05 +- 0.11, -3.01 +- 0.25; WM+Pl+QJ -3.24 +- 0.03, -3.34 +- 0.08.
  mask +-20: WM+Pl -3.12 +- 0.12, -3.01 +- 0.30; WM+Pl+QJ -3.26 +- 0.03, -3.25 +- 0.10.
- Sect. 5.2.3: "The addition of the QUIJOTE-MFI frequencies leads to the reduction of the error
  bars on beta_s by a factor 2 to 3."
- Sect. 5.2.3: "Finally, we note that the spectral indices recovered from the WM+Pl+QJ and WM+Pl
  datasets show a systematic difference, with the tension increasing with f_sky and reaching 2.3
  sigma for the 10-degree galactic cut where the signal-to-noise is larger." Sect. 6 says "at the
  ~ 2.4 sigma level for the |b| > 10 degree case".
- Sect. 5.2.3: "For example, de la Hoz et al. (2023) report a mean spectral index of -3.08 with a
  dispersion of 0.13 using a parametric component separation method, while a semi-blind approach
  yields an inverse-variance-weighted mean synchrotron index of -3.11 with a standard deviation of
  0.21 (Adak et al. 2026)."
- Prior, Sect. 3: "beta_s ~ N(-3, 1) and beta_d ~ N(1.6, 1)".
- Sect. 4.2 caveat: "In addition, the associated uncertainty might be slightly underestimated in
  the likelihood for two reasons."

### 12b. de la Hoz et al. 2023, QUIJOTE VIII

Identifier: arXiv 2301.05117v1, "QUIJOTE scientific results - VIII. Diffuse polarized foregrounds
from component separation with QUIJOTE-MFI". MNRAS 519, 3504-3525 (2023), DOI
10.1093/mnras/stac3020 (Crossref, BIB). This is paper VIII. Paper IV of the series is arXiv
2301.05113, Rubino-Martin et al., "QUIJOTE scientific results - IV. A northern sky survey in
intensity and polarization at 10-20GHz with the Multi-Frequency Instrument", DOI
10.1093/mnras/stac3439 (arXiv API, BIB; content not read).

Read: PDF, pp. 1, 4 to 11 and 16 to 19 of 25.

Where: abstract; Sect. 3.2 (priors); Sect. 4.1 (data); Sect. 5.1.1 to 5.1.3, Table 1, Fig. 7;
Sect. 5.5.2.

Quotes:
- Abstract: "We obtain an average value of -3.08 and a dispersion of 0.13, considering only pixels
  with reliable goodness-of-fit."
- Abstract: "Finally, we find that the beta_s estimation in the high signal-to-noise synchrotron
  emission areas is prior-independent while, outside these regions, the prior governs the beta_s
  estimation."
- Sect. 3.2: "and use a broad Gaussian distribution N(-3.1, 0.3) as a prior on beta_s." and
  "Finally, flat priors are used in the characterization of the amplitude parameters." Footnote 3:
  "N(x, sigma) represents a normal distribution with mean x and variance sigma^2."
- Sect. 5.1.3: "In the QS95 region, the mean and the standard deviation from the "Model 4" of
  Miville-Deschenes et al. (2008) template are -3.00 +- 0.05 while those from our estimate are
  -3.08 +- 0.13." For pixels with signal-to-noise on beta_s above 15: "we find that the mean value
  and dispersion of the distribution of beta_s are -3.12 +- 0.15 for our analysis".
- Sect. 5.1.3: "One may wonder if this result can be affected by the considered prior, since the
  estimated spectral indices for low signal-to-noise pixels are significantly constrained by it
  (see Section 5.5.2)." and "At high latitudes we cannot make a reasonable comparison since our
  beta_s estimates are driven by the prior."
- Sect. 5.5.2, with the prior N(-3.0, 0.3): "Outside the regions where the synchrotron emission is
  the largest, the difference is close to -0.1 which is the difference between the expected value
  of the priors. In other words, when there is not enough information from the data the recovered
  beta_s is close to the expected value of the prior."
- Sect. 5.1.1: "However, we note that in this analysis we have considered statistical
  uncertainties but not calibration errors, which in QUIJOTE are of 5%."
- Table 1 (uniform index per region): R1 (North Polar Spur), f_sky 4.84 percent, -3.028 +- 0.002;
  R3 -2.945 +- 0.008; R4 -3.319 +- 0.011; R5 -3.228 +- 0.019.
- Data, Sect. 4.1: MFI 11 and 13 GHz, WMAP K and Ka, all Planck polarized channels (PR3 and PR4),
  smoothed to 2 degrees FWHM at Nside = 64; f_sky = 51 percent (Fig. 1).

Verdict: CONFIRMED WITH CORRECTION.

Corrections:
1. QUIJOTE XII is internally inconsistent on the error of beta_s^EE at |b| > 20 degrees: the
   abstract and Sect. 6 give -3.26 +- 0.05, Table 4 gives -3.26 +- 0.03. Cite the abstract value
   and note the table, or ask the authors. BB is -3.25 +- 0.10 everywhere.
2. arXiv 2301.05117 is QUIJOTE VIII. The series paper IV is 2301.05113 and was not read here.
3. The mean -3.08 with dispersion 0.13 is a mean and scatter over pixels within the 95 percent
   goodness-of-fit region, not a value with an error on the mean.

(inference) The 2.3 sigma quoted in Sect. 5.2.3 for the 10 degree cut equals (3.22 - 3.00) /
sqrt(0.09^2 + 0.03^2) = 2.3, which matches the Table 4 error of 0.03 and not 0.05 (which would
give 2.1). This suggests that the table is the source of the text and that the abstract value is
either rounded up or includes something not described on the pages read.

## 13. Incidental parameters and the errors-in-variables model

### 13a. Neyman and Scott 1948

Identifier: J. Neyman and Elizabeth L. Scott, "Consistent Estimates Based on Partially Consistent
Observations", Econometrica, Vol. 16, No. 1, January 1948, first page 1, DOI 10.2307/1914288
(Crossref, BIB).

Read: BIB only. The Crossref record gives the first page only; the page range (usually cited as
1 to 32) is not in the record and the JSTOR page did not load, so it is not verified. The content
was not read.

Verdict: CONFIRMED (bibliographic).

### 13b. Linear functional model with known error variances

No single source could be read that is both standard and open. What was verified:

1. Engel, B. (1985), "The linear functional relationship model for two variables", Kwantitatieve
   Methoden 19, pp. 5-29 (journal of the Dutch statistical society, scanned PDF at
   https://www.vvsor.nl/wp-content/uploads/2020/06/KM1985019001.pdf). Read: PDF, journal pages 5
   to 20 and 27 to 29 (of 5 to 29). This is a review, not a primary source. Model: Y_t1 = x_t1 +
   e_t1, Y_t2 = x_t2 + e_t2, x_t2 = alpha + beta x_t1, Var(e_t) = Omega = sigma^2 [[1, phi rho],
   [phi rho, phi^2]], with the x_t1 unknown constants, t = 1..T (Eqs. 3 to 5).
   - Without knowledge of the error variances (Sect. 0.4, p. 9): "However, beta^_ML is not a
     consistent estimator for beta (T -> infinity) and does not correspond to a maximum but to a
     saddle point of the likelihood surface [Copas (1972), Kendall and Stuart (1968), 29.14,
     p.399-400]." and "with each pair of observations (Y_t1, Y_t2) an unknown and independent
     parameter x_t1 is entered into the model. For infinitely many observations (T -> infinity) we
     end up with infinitely many unknown parameters."
   - Theorem 1 (Sect. 1.1, p. 12): with Omega known, "possibly upto (sic) a constant", maximum
     likelihood under normality is equivalent to Sprent's generalized least squares.
   - Theorem 2 (Sect. 2.1, p. 15): "For the functional relationship model described by (3) and (4)
     with Omega known, the following conditions are sufficient for alpha^ and beta^ to be
     consistent estimators of alpha and beta: lim x-bar_.1 exists and takes a finite value, say
     x_0. lim (1/T) sum (x_t1 - x-bar_.1)^2 exists and is equal to a finite positive value, say
     m_0."
   - Variances (Sect. 2.3, p. 17), the decisive sentences: "Results derived by Patefield (1978)
     are given below. These expressions are not the same as those derived from Fisher's
     information matrix, see Patefield (1977). Therefore the results obtained by Barnett (1973)
     (sic; the reference list has Barnett 1970) for the case rho = 0 are not correct although as
     we shall see in section 3 in practice it will not really matter." Eq. (21): Var_as(beta^) =
     (1/T) (sigma^2/m_0^2) {psi m_0 + (1 - rho^2) phi^2 sigma^2}, with psi = phi^2 - 2 beta phi
     rho + beta^2 (Eq. 13).
   - Size of the error (Sect. 3, p. 18): "A similar observation can be made with respect to the
     results derived from Fisher's information matrix in Barnett (1970). Although these results
     are incorrect the difference with (21), (22) and (23) is O({(sigma^2/K)/m_0}^2) and not
     really important." Eq. (23) is the leading term, Var_as(beta^) = (1/T) psi (sigma^2/K)/m_0,
     for K replicates.
   - Sect. 0.5 (p. 10) states the contrast: with T fixed and replicates K -> infinity, "ordinary
     maximum-likelihood theory applies" and variances "may be derived from Fisher's information
     matrix"; with T -> infinity and K fixed, "we meet the problems discussed in the previous
     section."
2. Patefield, W. M. (1977), "On the Information Matrix in the Linear Functional Relationship
   Problem", J. R. Stat. Soc. C (Applied Statistics) 26(1), 69-70, DOI 10.2307/2346870. Read: ABS
   (RePEc record, through the summarizer, which paraphrased all but one sentence): Barnett (1970)
   derived an asymptotic covariance matrix for the maximum-likelihood intercept and slope by
   partitioning the information matrix; "It is pointed out that his procedure is incorrect."; a
   correct matrix and an estimator of it are given. Full text not accessible.
3. Patefield, W. M. (1978), "The unreplicated ultrastructural relation: large sample properties",
   Biometrika 65, 535-540. BIB only, as listed in Engel's reference list. Not read.
4. Gleser, L. J. (1981), "Estimation in a Multivariate "Errors in Variables" Regression Model:
   Large Sample Results", Ann. Statist. 9(1), 24-44, DOI 10.1214/aos/1176345330. Read: ABS
   (Crossref). For error covariance sigma^2 I: "Under reasonable assumptions concerning the
   unknown sequence {u_1i, i = 1,2,...}, alpha^, B^ and r^{-1}(r + p) sigma^^2 are shown to be
   strongly (with probability one) consistent estimators of alpha, B and sigma^2, respectively, as
   n -> infinity, regardless of the common distribution of the errors e_i." (inference) For p = r
   = 1 this says that the maximum-likelihood estimate of sigma^2 converges to sigma^2/2 while the
   slope estimate is consistent.
5. Fuller, W. A. (1987), Measurement Error Models, Wiley, ISBN 9780471861874, DOI
   10.1002/9780470316665 (Crossref, BIB). Text not accessible: NOT FOUND as a source of a quotable
   statement. Kendall and Stuart, The Advanced Theory of Statistics, vol. 2: not accessible, cited
   only through Engel (Sects. 29.14 and 29.15, pp. 399 to 401 of the 1968 Griffin edition, per
   Engel).

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement, with its source: in the linear functional relationship with a known error
covariance matrix (known up to a constant is enough), the maximum-likelihood estimators of slope
and intercept coincide with generalized least squares and are consistent as the number of points
grows (Engel 1985, Theorems 1 and 2; Gleser 1981, abstract, for equal error variances). Their
asymptotic variances "are not the same as those derived from Fisher's information matrix" (Engel
1985, Sect. 2.3, citing Patefield 1977 and 1978). The discrepancy is of second order in the ratio
of the error variance to the spread of the true values, sigma^2/m_0, so it is negligible at high
signal-to-noise and matters at low signal-to-noise. "Inconsistent standard errors" overstates it:
the information-matrix variance is wrong by a factor that tends to 1 as sigma^2/m_0 tends to 0.

(inference) From Engel's Eqs. (21) and (23) with K = 1, the ratio of the correct asymptotic
variance to the leading term is 1 + (1 - rho^2) phi^2 sigma^2 / (psi m_0). For uncorrelated errors
of equal variance this is 1 + sigma^2 / ((1 + beta^2) m_0). Engel does not print Barnett's
formula, so the sign and size of the error of the information-matrix result beyond "second order"
are not established by what was read. Paper A1 should cite Patefield (1977, 1978) and Gleser
(1981) only after reading them in a library, and may cite Fuller (1987) as the textbook without
attributing a specific sentence to it.

## 14. Planck 2018 results V: template coefficients with 30 GHz and with K band

Identifier: arXiv 1907.12875v2. A&A 641, A5 (2020), DOI 10.1051/0004-6361/201936386 (arXiv API).

Read: HTML-S. ar5iv characters 0 to 200 000 of 297 187; arxiv.org/html/1907.12875v2 characters
100 000 to 200 000 of 301 492 for an independent check of Sect. 2.2.5. The PDF exceeds the fetch
size limit.

Where: Sect. 2.2.1, Eq. (5) and Table 2; Sect. 2.2.5, paragraph "Synchrotron-tracer stability",
Fig. 16.

Quotes:
- Eq. (5): m^_{100,143} = (m_{100,143} - alpha m_30 - beta m_353) / (1 - alpha - beta), "using the
  30-GHz full-channel map as a tracer of polarized synchrotron". The sentence after Eq. (5) says
  "where alpha and beta are the amplitudes of dust and synchrotron templates, respectively", which
  is the reverse of the equation (flagged by the summarizer; not checked in the PDF).
- Table 2, "Template scalings measured on data." "The uncertainties are obtained from chi^2
  minimization. The errors in square brackets are computed from the 300 FFP10 simulations." Row
  100 GHz: alpha x 10^2 = 1.83 +- 0.12 [0.18], beta x 10^2 = 1.950 +- 0.014 [0.015]. Row 143 GHz:
  no alpha, beta x 10^2 = 4.078 +- 0.011 [0.013]. "we fit both synchrotron and dust at 100 GHz,
  but only dust at 143 GHz".
- Sect. 2.2.5, both renderings: "First we process with the usual smoothing procedure the map and
  covariance of the WMAP K-band data." "We substitute WMAP K band for Planck 30 GHz in the
  cleaning of the 100-GHz map," "finding for synchrotron a scaling of alpha=0.0119+-0.0006."
  "Secondly we clean both 100 and 143 GHz with 30 GHz, finding for synchrotron a scaling of
  alpha=0.0073+-0.0011" "for the 143-GHz map." "the posterior on tau is not substantially affected
  by these variations in the handling of synchrotron contamination."
- Other values met on the way (ar5iv): LFI 70 GHz cleaning coefficients "alpha=0.058+-0.004,
  beta=0.0092+-0.0004" (Sect. 2.3.2); tau = 0.0506 +- 0.0086 (lowE, Sect. 2.2.3), 0.063 +- 0.020
  (LFI, Sect. 2.3.4), 0.062 +- 0.012 (WMAP with the 353 GHz dust template, Sect. 2.3.6).

Verdict: CONFIRMED WITH CORRECTION.

Corrected statement: the paper gives the two synchrotron template coefficients for the 100 GHz
map, alpha = (1.83 +- 0.12 [0.18]) x 10^-2 with Planck 30 GHz as tracer (Table 2) and alpha =
0.0119 +- 0.0006 with WMAP K band as tracer (Sect. 2.2.5), and says that tau is not substantially
affected. The ratio 1.54 +- 0.13 and the expectation of 1.96 to 1.98 are not in the paper. No
sentence in the parts read compares the two coefficients or gives an implied spectral index.

(inference) Ratio 0.0183/0.0119 = 1.54. Naive error from the chi^2 errors: 0.13. With the FFP10
simulation error of 0.18 x 10^-2 for the 30 GHz coefficient it is 0.17. Expected for a power law
of index -3.1 between 22.8 and 28.4 GHz: (28.4/22.8)^3.1 = 1.98 in brightness temperature and
1.96 in thermodynamic temperature. The two fits share the 100 GHz map and the dust template, so
the errors are not independent, and the treatment of template noise in the fit affects both
coefficients. With the larger error the difference from 1.96 is about 2.5 sigma instead of 3.2.
If Paper A1 uses the ratio it must present it as its own arithmetic.

## 15. Flags against the repo notes

1. `2026-11-09_science_landscape_and_data.md`, Sect. 1 item 2, states "template coefficient ratio
   30 GHz over K of 1.54 +- 0.13 in Planck 2018 V" as if the paper gave it. Sect. 3.3 of the same
   note marks it as inference. Only the two coefficients are in the paper (item 14).
2. Same note, row 2610.31469, "UNRESOLVED: Table 4 extraction gave +- 0.03": resolved. Both
   numbers are in the paper (item 12a).
3. `2026-11-09_novelty_check_methods.md`, Sect. 6, unresolved "1608.00551 Fig. 11 text against
   caption" and "0811.4280 92 against 96": resolved. Both discrepancies are in the papers
   themselves (items 7 and 6).
4. All three 2026-11-09 notes and `nextgen/LITERATURE.md` call arXiv 2011.08503 "BeyondPlanck XV".
   The published title is "BeyondPlanck XIV" (item 10). The same check should be made for the
   other BeyondPlanck numbers used in the notes (for example 2201.03530, called XVI there), since
   arXiv and A&A numbering differ; not checked here.
5. `2026-11-09_novelty_check_methods.md`, claim 1, describes Stompor's Eq. (10) as having "the
   extra term ln|(A^t N^-1 A)^-1|" without a sign. As printed the term enters -2 ln L with a plus
   sign, which is opposite to the Gaussian integral and to Alonso et al. Eq. (A7) (item 4,
   inference).
6. Same note, claim 2, places the Jeffreys prior of 0709.1058 in "Sect. III.1 and Eq. (36)". Eq.
   (36) is in Sect. 4.2; the formula in Sect. 3.1 is unnumbered and contains the factor G(nu).
7. `2026-11-09_science_landscape_and_data.md`, Sect. 1, says the QUIJOTE XII value is "0.15 to
   0.25 steeper than the same collaboration's map-domain values". The paper itself says its values
   are "in agreement with" and "consistent with" the real-space analyses (abstract, Sects. 5.2.3
   and 6). The note's wording is its own inference and should be marked so.
8. Same note, row 2011.08503, "posterior sd over about 500 samples, 2 chains": "500 ensemble Gibbs
   samples" is from the caption of Fig. 14, which concerns the dust index. For synchrotron the
   text read gives two chains with 10 burn-in samples each and no sample count (HTML-S).
9. `2026-11-09_noisy_mh_and_translation_move_prior_art.md`: every quote from Ball et al., Andrieu
   and Roberts, and Bardenet et al. that I re-read is confirmed word for word. No contradiction.
10. Dunkley et al., Weiland et al. 2018 and 2022, Fuskeland et al. 2021 and de la Hoz et al.:
    the numbers in the landscape note's tables that I re-read agree with the papers. Not
    re-checked: the section number given there for the 0.1 +- 0.2 flattening of 1909.05923 (I
    found it in the abstract), and App. A of 2301.05117.

## 16. What was looked for and not found

- A numerical value of nu_0,s in the text of 2011.08503 or in the first 300 000 characters of
  2011.05609. Found only in the Fig. 6 caption of the former and in Table 2 of 2201.08188.
- A tau number for the widened index prior in 2011.08503 (one qualitative sentence only).
- Any sentence in 1907.12875 that compares the 30 GHz and K-band coefficients.
- An explicit invariant law for MCWM in 0903.5480 (pages read).
- A Gaussian-noise assumption in Sect. 6.3 of 1505.02827 (the text says the opposite).
- QUIJOTE data in 1909.05923 (pages read).
- A statement about calibration uncertainty in the fit of 2610.31469 (pages read; the HTML check
  of the first 100 000 characters found none).
- A quotable sentence from Fuller (1987), Kendall and Stuart, Patefield (1977, 1978) or Neyman and
  Scott (1948): texts not accessible.

Sources tried without success: arXiv PDFs above 10 MB (listed at the top);
`www.aanda.org` full text and PDF of aa43160-22 (HTTP 403); `air.unimi.it/handle/2434/981809`
(HTTP 403); `www.jstor.org/stable/1914288` (page did not load); Project Euclid pages for Gleser
1981 (empty content). Bibliographic checks used the arXiv API (`export.arxiv.org/api/query`),
Crossref (`api.crossref.org/works/<DOI>`), INSPIRE (`inspirehep.net/api/literature`), the JMLR
page `jmlr.org/papers/v18/15-205.html` and RePEc. One web search each was run for a mirror of the
2011.08503 PDF and for Patefield (1977).

## 17. Gaps and suggested follow-up

1. Read 2011.08503 from a local copy of the published PDF and confirm the quotes of item 10, the
   wording and place of the K-band footnote, and the number of synchrotron samples. Decide how
   Paper A1 cites it: "BeyondPlanck XIV, A&A 675, A14" with the arXiv number.
2. Check Eq. (10) of Stompor et al. in the MNRAS version for the sign of the determinant term, and
   read pp. 7 to 14 of the arXiv version.
3. Read Planck 2018 V Sect. 2.2.1 and 2.2.5 in the PDF (two paragraphs and Table 2) to confirm the
   HTML-S quotes and the alpha and beta naming after Eq. (5).
4. Library reads before citing: Patefield (1977, 1978), Gleser (1981), Fuller (1987, the sections
   on the functional model with known error variance), Neyman and Scott (1948) for the page range
   and for its own treatment of the straight line with errors in both variables.
5. Ask the QUIJOTE XII authors, or wait for the published version, about +-0.05 against +-0.03.
6. The PRL version of Ball et al. (equation number and wording) and the JMLR version of Bardenet
   et al. (section number) were not checked against the arXiv versions read.
7. Unread remainders: 0903.5480 pp. 6 to 17 and 21 to 27; 1505.02827 pp. 3 to 19 and beyond 23;
   0709.1058 pp. 13 to 23; 1608.00551 pp. 1 to 3 and 8 to 14; 1801.01226 pp. 11 to 13; 1404.5323
   v1 pp. 10 to 11 and all of v2 except the abstract; 1909.05923 pp. 4 to 6 and 14 to 16;
   2610.31469 pp. 4 to 6 and 15 to 18; 2301.05117 pp. 2 to 3, 12 to 15 and 20 to 25; 1907.12875
   beyond character 200 000; 2011.05609 beyond character 300 000; 2201.08188 beyond character
   100 000.
