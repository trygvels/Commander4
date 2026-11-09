# Prior art: (A) the stationary law of Metropolis with a re-estimated acceptance ratio, (B) the deterministic translation move for linear amplitudes and a nonlinear parameter

Date: 2026-11-09. Author: lit-scout agent. Scope: the two questions only.

Read levels used below:
- DIRECT: I read the stated PDF pages myself (page images of the arXiv PDF). Quotes at this level
  are transcribed by me. Greek letters and typographic dashes are written in ASCII.
- SUM: full text read through an HTML rendering (ar5iv or arXiv HTML) by a summarizer that quotes
  in short pieces. SUM-part: only the stated character range was read. Verify in the PDF before
  citing a quote or a number at this level.
- ABS: abstract only. BIB: bibliographic record only (Crossref or publisher page), content not read.
- (inference) marks my own synthesis. It is not a result of any paper.

## 0. Verdicts

**A. In print, in physics, under the name "stochastic annealing"; close to in print in statistics.**
Ball, Fink and Bowler 2003 (arXiv cond-mat/0301179, Phys. Rev. Lett. 91, 030201), Eq. (7), give for
acceptance on the sign of a Gaussian-noisy energy difference of standard deviation sigma:
T ~= sqrt(pi/8) sigma (1 - 0.018 (Delta E / T)^2 + ...), and say that detailed balance is "not
strictly achieved". That is the same constant, the same mechanism and the same small-step condition
as the E035 law. Bardenet, Doucet and Holmes 2017 (arXiv 1505.02827, Sect. 6.3) give the same
mechanism with a different constant (a global logistic fit with lambda = 1.702, so T = sd / 1.702 =
0.588 sd instead of 0.627 sd, the conversion being my inference). Not found anywhere: the statement
for the "Monte Carlo within Metropolis" or "noisy Metropolis-Hastings" chain in the statistics
papers named in the task (they give existence conditions and total-variation bounds, no explicit
law), the chi-square form T = sqrt(2 pi k) / 4 and the width factor 0.79 k^(1/4), or any use in
CMB or astronomy. The methods paper must cite Ball et al. for the formula and can claim only the
identification of the legacy move as a refreshed pseudo-marginal chain and the k^(1/4) consequence.

**B. The move class is published; the move for spectral indices is not found.** In CMB analysis
Racine et al. 2016 (arXiv 1512.06619, Eqs. 1, 17 and A2 to A8) publish the joint move
s' = s_hat(theta') + F (s - s_hat(theta)) with a deterministic filter F on the fluctuation; the
translation move is the member F = 1 of that family (inference), which the paper does not single
out. In statistics it is a Metropolis update of the parameter in the data-corrected linear
reparametrization X = v^(1/2) X* + m(Theta, Y) with m = E(X | Theta, Y) of Papaspiliopoulos,
Roberts and Skold 2007 (arXiv 0708.3797, Eq. 5.1, "partially noncentered"), location part only
(inference). No CMB or radio component-separation paper was found that shifts or rescales the
amplitudes deterministically when a spectral index is proposed. The closest published things are
the amplitude-marginal step and the "ridge" reset of BeyondPlanck XIII (arXiv 2201.08188) and its
per-proposal refit of band monopoles.

## A. Metropolis with a re-estimated (noisy) acceptance ratio

### A.0 Notation bridge (inference)

E035 has log acceptance ratio -(Delta + Z)/2 with Var(Z) = 4k, so the noise on the log ratio has
sd sigma = sqrt(k), and T = sqrt(2 pi k)/4 = sigma sqrt(pi/8). For sigma much larger than 1 the
Metropolis rule with noise accepts with probability Phi(-Delta_log / sigma) + O(1/sigma), which is
the "accept if the noisy difference is favourable" rule analysed by Ball et al. and by Bardenet et
al. Sect. 6.3. The published results below are for that sign rule, not for the Metropolis rule at
finite sigma. No source was found that treats Metropolis plus Gaussian noise at finite sigma and
gives an effective temperature.

### A.1 "Monte Carlo within Metropolis" (MCWM) against GIMH

- **Andrieu and Roberts 2009**, Ann. Statist. 37, 697, DOI 10.1214/07-AOS574, arXiv 0903.5480.
  DIRECT, pages 1 to 10 and 18 to 29 of 29 (pages 11 to 17, the GIMH convergence proofs, not read).
  - Definition and credit (p. 3): "both Z and Z* are "refreshed" at each iteration independently of
    previously sampled auxiliary variables given theta and theta*. This algorithm is referred to as
    the Monte Carlo within Metropolis (MCWM) in [1] following the terminology of [7]." Their [1] is
    Beaumont 2003 and [7] is O'Neill et al. 2000 (reference list, p. 29).
  - Stationary law (p. 3): "MCWM and the marginal algorithm P are not equivalent. In particular,
    pi(theta) is typically not the invariant distribution density of P~_N^MCWM and therefore will
    not produce samples from pi(theta) even in steady state." It "will asymptotically be
    distributed according to an approximation of pi(theta), which should be all the more precise
    that N is large."
  - GIMH (p. 4): "GIMH can be interpreted as a form of MCWM where Z is "recycled" from the previous
    iteration"; it is a Metropolis-Hastings chain on the extended space with marginal pi(theta).
  - Sect. 6 (p. 18): "the existence of an invariant distribution for P~_N^noisy is not obvious in
    general (it is not a MH update)." Theorem 9 (p. 19) assumes an invariant law exists and bounds
    its total-variation distance to pi, Eq. (6.2): ||pi - pi_N|| <= C epsilon / (1 - rho).
  - No explicit form of the MCWM invariant law, no tempering or temperature language.
- **Beaumont 2003**, Genetics 164, 1139, DOI 10.1093/genetics/164.3.1139. ABS only (publisher
  page; the full text was not accessible: the publisher page is paywalled and PubMed Central
  refused automated access). What is said above about it is Andrieu and Roberts' account. Whether
  Beaumont himself remarks on the width of MCWM posteriors: not found (not read).
- **O'Neill, Balding, Becker, Eerola, Mollison 2000**, J. R. Stat. Soc. C 49, 517, DOI
  10.1111/1467-9876.00210. BIB plus two abstract sentences from Crossref. Its role as the source
  of the term MCWM is Andrieu and Roberts' statement. Content not read.

### A.2 "Noisy Metropolis-Hastings"

- **Alquier, Friel, Everitt, Boland 2016**, Stat. Comput. (DOI 10.1007/s11222-014-9521-x as given
  in the reference list of 1503.07066), arXiv 1403.5496. SUM-part (first 100 000 of 108 986
  characters). Sect. 2.1: "Of course there is no reason for pi to be invariant under P^". Sect.
  2.4: "In MCWM the estimate of the target in the denominator is recomputed at every iteration of
  the MCMC, whereas in GIMH it is reused from the previous iteration." Results are total-variation
  bounds between the exact and the noisy chain (Corollary 2.3, Theorems 2.2 and 3.4); no formula
  for the limiting law, no tempering language.
- **Medina-Aguayo, Lee, Roberts 2016**, "Stability of noisy Metropolis-Hastings", arXiv 1503.07066.
  SUM (both halves, 163 587 characters). Sect. 1.3: the noisy algorithm "is commonly known as
  Monte Carlo within Metropolis (MCWM)"; "the noisy chain is not invariant under pi and it is not
  reversible in general"; it "may not even have an invariant distribution". Transient examples in
  Sects. 2.3 and 2.4; geometric ergodicity for large N (Theorems 3.2, 3.3); convergence of the
  noisy invariant law to pi in total variation with rates (Theorem 4.1, Theorem 4.2, Proposition
  4.1). Sect. 2.2 speaks of "an unknown noisy invariant distribution". No closed form in any
  example, including the homogeneous-weight case of Sect. 2.1 and the log-normal weights of Sect.
  1.4. No tempering language.
- **Nicholls, Fox, Watt 2012**, "Coupled MCMC with a randomized acceptance probability", arXiv
  1205.6857. SUM. "The naive algorithm does not in general target pi, and may not even have an
  equilibrium distribution"; it names "The 'MCWM' algorithm in Andrieu and Roberts 2009" as a
  naive algorithm; Fig. 2 shows naive samplers "slightly over-dispersed". No stationary law. It is
  the only statistics paper found that cites Ball et al. (for an exact "universal rule").
- Also checked: Llorente, Martino, Read, Delgado, survey, arXiv 2108.00490 (SUM-part, first 100 000
  of 179 600 characters): no explicit stationary density for MCWM. Medina-Aguayo, Rudolf,
  Schweizer, arXiv 1809.09547 (ABS): perturbation bounds only.

### A.3 Physics

- **Ball, Fink, Bowler 2003**, "Stochastic Annealing", arXiv cond-mat/0301179 (the arXiv listing
  gives DOI 10.1103/PhysRevLett.91.030201). DIRECT, all 5 pages. This is the decisive source.
  - Setting (p. 1): moves are accepted with probability A(x) "where x is only an estimate of the
    energy change. We will assume that these estimates have statistically independent errors."
    Eq. (6): K(Delta E) = integral f(x | Delta E) A(x) dx.
  - The rule (p. 1): "the crude choice A(x) = 1 for x < 0 and A(x) = 0 otherwise."
  - The decisive sentences (p. 2): "The resemblance between this and the Glauber acceptance
    function (whose symmetry it shares) is striking, showing how the random energy errors make
    the selection look thermal - although in this case the match is not exact, so detailed balance
    is not strictly achieved. The standard deviation sigma of the Gaussian distribution controls
    an approximate effective temperature using this rule, as inverting eq. (3) gives
    T ~= sqrt(pi/8) sigma (1 - 0.018 (Delta E/T)^2 + order (Delta E/T)^4). (7)"
  - "the departure from detailed balance, due to the energy dependent terms at large
    |Delta E/T|, is not a problem" for optimisation (p. 2).
  - sigma is the sd of the estimated difference: "sigma is the standard deviation of x" (p. 4).
  - Their exact Gaussian rule, Eq. (20), A(x) = min(1, exp(-beta (x + beta sigma^2 / 2))), is the
    penalty formula of Ceperley and Dewing (inference: identical at beta = 1); the 14-entry
    reference list cites neither Ceperley and Dewing nor Kennedy and Kuti.
  - (inference) I re-derived Eq. (7) from K = Phi(-Delta E / sigma): ln[K(-D)/K(D)] = 2 sqrt(2/pi)
    (D/sigma) [1 + (2/(3 pi) - 1/6)(D/sigma)^2 + ...], which gives the constant sqrt(pi/8) and the
    coefficient 0.0179. The correction term has the sign of E035 result P4 (wider proposals give a
    slightly narrower law).
- **Ceperley and Dewing 1999**, J. Chem. Phys. 110, 9812, DOI 10.1063/1.478034, arXiv
  physics/9812035. DIRECT pages 1 to 4; SUM for the rest. P. 3: "one could simply ignore the
  effect of fluctuations in the estimate of Delta V(s). If the errors are small then clearly the
  sampled distribution will be changed only a little. If the acceptance ratio as a function of
  Delta V(s) were a linear function there would be no bias, but because it is non-linear,
  fluctuations will bias the asymptotic distribution." Eq. (13): a_P(delta; sigma) = min(1,
  exp(-delta - sigma^2/2)). SUM: Fig. 3 discussion, "ignoring the noise leads to a much smoother
  density than the exact result". No effective-temperature statement and no formula for the bias
  of the uncorrected walk.
- **Kennedy and Kuti 1985**, Phys. Rev. Lett. 54, 2473, DOI 10.1103/PhysRevLett.54.2473. ABS (as
  paraphrased by the summarizer from the publisher page): a method that tolerates large unbiased
  noise in the relative probability. Per Lin, Liu, Sloan 2000 (arXiv hep-lat/9905033, SUM) it is
  a linear acceptance rule that needs an unbiased estimate of the ratio, and "Simple
  exponentiation of this estimator ... inevitably yields a bias". Nothing on temperature found.
- **Lattice QCD and quantum Monte Carlo after 1985.** Lin, Liu, Sloan (hep-lat/9905033, SUM);
  Bakeyev and de Forcrand (hep-lat/0008006, SUM); Pierleoni and Ceperley (physics/0510254,
  SUM-part, first 100 000 of 129 637 characters: uncorrected results "will be biased, the amount
  of bias increasing for increasing noise level"); Doll and Freeman 1994 (chem-ph/9405002, SUM:
  Gaussian noise on the density ratio inside the Metropolis test, "the guarantee of detailed
  balance is lost", effects "surprisingly modest"); Sanderski, Gentinetta, Carleo 2026 (arXiv
  2609.01204, ABS: "a noisy energy biases the sampled distribution"). None states an effective
  temperature. These papers build exact corrections and do not characterize the uncorrected law.
- **Branke, Meisel, Schmidt 2008**, J. Heuristics 14, 627, DOI 10.1007/s10732-007-9058-7. BIB plus
  repository abstract ("an improved version of stochastic annealing that allows for arbitrary
  annealing schedules"). Full text not accessible. Bardenet et al. 2017 credit it with the first
  use of the sign rule; Ball et al. 2003 predates it. Whether it restates Eq. (7): not found.

### A.4 Minibatch ("tall data") Metropolis

- **Bardenet, Doucet, Holmes 2017**, "On Markov chain Monte Carlo methods for tall data", arXiv
  1505.02827. DIRECT, PDF pages 15 to 23 (Sects. 5.2 and 6.1 to 6.3); SUM-part for the rest of the
  first 100 000 characters.
  - Sect. 6.3 (p. 22): "Why, then, not rely on the subsampling noise to guarantee exploration, and
    accept a move if and only if" the subsampled log ratio is positive, Eq. (22). "This idea has
    been first used by Branke et al. (2008)".
  - Sect. 6.3 (p. 23): with |Phi(x) - 1/(1 + exp(-lambda x))| <= C and "C = 0.0095 and lambda =
    1.702" (Bowling et al. 2009): "the acceptance probability of an algorithm that would accept
    the move from theta to theta' by checking whether Lambda*(theta, theta') > 0 is close to the
    acceptance probability of an MCMC algorithm with a Baker acceptance criterion (Robert &
    Casella, 2004, Section 7.8.1) that targets pi^beta with temperature beta = lambda sqrt(t) /
    (n sigma)." And: "a clear link between the variance of the subsampled log likelihood ratios
    and the temperature."
  - (inference) The sd of the estimated total log ratio is s = n sigma / sqrt(t), so this is
    T = s / 1.702 = 0.588 s. The constant differs from sqrt(pi/8) = 0.627 because 1.702 is a
    global fit of a logistic to the normal distribution function, while sqrt(pi/8) matches the
    slope at zero, which is what governs small steps. In temperature the two differ by 6 percent,
    in posterior width by 3 percent.
  - Sect. 6.2.2 (p. 20), on Korattikara's test in a Gaussian example: "the posterior is rightly
    centered, but is slightly too wide. This is a tempering effect due to too small subsamples".
  - Sect. 6.1 ("naive subsampling is nontrivial tempering") concerns an estimate kept with the
    state, whose target is p(theta) E exp(n l^(theta)); it is not the refreshed chain.
- **Korattikara, Chen, Welling 2014**, arXiv 1304.5299. SUM. Sequential t-test; Theorem 1 bounds
  the total-variation distance of the stationary law by Delta_max / (1 - eta). No explicit law,
  no temperature. Does not cite Ceperley and Dewing, Beaumont or Nicholls et al.
- **Bardenet, Doucet, Holmes 2014**, ICML (PMLR 32). ABS: the approximate chain "samples from a
  perturbed version of the target distribution". No tempering in the abstract.
- **Seita, Pan, Chen, Canny 2017**, arXiv 1610.06848. SUM. Barker test written as "Delta(theta,
  theta') + X > 0" with X logistic (Eq. 10); a correction variable makes Gaussian minibatch noise
  logistic, X_log = X_norm + X_corr (Eq. 12), possible only for noise sd below about 1.814. App.
  B: the closest Gaussian has "standard deviation approximately 1.7" with "sup|L(x) - Phi(x/1.7)|
  < 0.01". No statement of what uncorrected Gaussian noise does to the target.
- **Wu, Wang, Wong**, arXiv 1908.02910 (SUM-part, first 100 000 of 140 946 characters), and **Li
  and Wong 2017**, arXiv 1707.09705 (SUM): minibatch MH that "gives rise to approximately a
  tempered stationary distribution", with the temperature designed in by scaling the minibatch
  log likelihood (T = n / c_n, and T = n^(1 - lambda) respectively). A different mechanism from
  noise acting as temperature; no formula in the noise sd.

### A.5 Answer to the four sub-questions

1. MCWM against GIMH: MCWM does not have the target as invariant law; an invariant law need not
   exist; if it does it converges to the target as the number of importance samples grows. No
   explicit form (Andrieu and Roberts 2009, DIRECT).
2. Noisy MH: the same, with rates. No explicit form, no temperature (Alquier et al., Medina-Aguayo
   et al., both SUM).
3. Physics: yes. Ball, Fink, Bowler 2003, Eq. (7): T ~= sqrt(pi/8) sigma for the sign rule with
   Gaussian noise, approximate, with a stated first correction. Kennedy and Kuti and Ceperley and
   Dewing do not state it.
4. Minibatch: Bardenet et al. 2017 Sect. 6.3 state a tempered target pi^beta for the sign rule
   with beta = lambda sqrt(t) / (n sigma), lambda = 1.702. Korattikara et al. and Seita et al.
   give bounds or an exact correction, not a law.

Plain statement: "T = sigma sqrt(pi/8) in the large-noise small-step limit" is **in print** (Ball
et al. 2003, Eq. 7, for the sign rule that the Metropolis rule reduces to at large noise, the
reduction being my inference and E035's own expansion) and **close to in print** in the tall-data
literature (Bardenet et al. 2017, same mechanism, constant 1/1.702). **Not found**: the law
attached to MCWM or noisy MH by name; the chi-square case, T = sqrt(2 pi k)/4 and 0.79 k^(1/4);
the dependence of the coefficient on the number of index steps per refresh (E035 P4); any
occurrence in CMB analysis.

Suggested wording for the paper (inference): the legacy move is Monte Carlo within Metropolis
(O'Neill et al. 2000; Beaumont 2003; Andrieu and Roberts 2009) for the ridge posterior, whose
invariant law is known to differ from the target and is in general not explicit; in the
large-noise limit its acceptance is the sign rule of stochastic annealing, for which Ball, Fink
and Bowler (2003) give the effective temperature sqrt(pi/8) sigma; with sigma = sqrt(k) this
gives the width factor 0.79 k^(1/4).

## B. The deterministic translation move

### B.0 The move in standard terms (inference)

With u = a - a_hat(beta), the move is a Metropolis update of beta at fixed u followed by the map
back to a. The maps T_delta: (beta, a) -> (beta + delta, a + a_hat(beta + delta) - a_hat(beta))
compose as T_delta2 T_delta1 = T_(delta1 + delta2), so they form a one-parameter group acting on
the joint space with unit Jacobian. The acceptance ratio is the amplitude-marginal ratio times
N(u; 0, C(beta')) / N(u; 0, C(beta)), where C is the conditional amplitude covariance: collapsed
sampling up to the change of C along the move. Rescaling u by C(beta')^(1/2) C(beta)^(-1/2) as
well would give the collapsed step exactly.

### B.1 CMB relatives

- **Jewell, Eriksen, Wandelt et al. 2009**, arXiv 0807.0624 (ApJ 697, 258 according to the
  reference list of 1512.06619). DIRECT pages 3 to 5; SUM for the rest.
  - Sect. 3 (p. 3): the aim is "a large variation in C_l in the low signal-to-noise regime" with
    "an associated deterministic change in the CMB map"; "The motivation for a deterministic
    change is that it will avoid introducing ratios of determinants which we cannot compute."
  - General form, Eq. (14): s_(n+1) = F(s_n, C^(n+1), C^(n)), with the detailed-balance condition
    Eq. (15): F^-1(s_(n+1), C^(n+1), C^(n)) = F(s_(n+1), C^(n), C^(n+1)). Their choice, Eq. (16):
    s_(n+1) = (C^(n+1) / C^(n))^(1/2) s_n. Acceptance on the chi-square ratio, Eq. (17).
  - P. 3 to 4: Gibbs sampling in whitened maps x = C^(-1/2) s gives an accept probability
    "numerically identitical to the above"; "there could be other deterministic proposal schemes
    or another change of variables which lead to improvements".
  - P. 4: the step "is not irreducable", "the phases are unchanged in each MCMC step", so it is
    alternated with a standard Gibbs step.
  - A scale move with no mean shift. SUM: no statistics citation for the move apart from Green
    1995; foregrounds appear only as motivation.
  - (inference) The translation move satisfies Eq. (15) with F(a, beta', beta) = a + a_hat(beta')
    - a_hat(beta).
- **Racine, Jewell, Eriksen, Wehus 2016**, arXiv 1512.06619 (ApJ, published 2016 March 15
  according to the paper's header). DIRECT, all 9 pages. The closest CMB relative.
  - Abstract: the low-S/N inefficiency of the rescaling move "is mitigated in the new algorithm by
    subtracting the Wiener filter mean field from the proposed map prior to rescaling".
  - Eqs. (1), (17): s^(i+1) = s_hat^(i+1) + S_(i+1)^(1/2) S_i^(-1/2) (s^i - s_hat^i), where s_hat
    is the Wiener-filter mean evaluated at the current and at the proposed parameters.
  - App. A, Eq. (A3): "general proposals of the form" s^(i+1) = s_hat^(i+1) + F(theta^(i+1),
    theta^i)(s^i - s_hat^i) + beta^(1/2) xi, taken to the "deterministic proposal" limit (A5),
    with the condition F(theta^i, theta^(i+1)) = F^-1(theta^(i+1), theta^i) (A8). Two filters are
    given: S^(1/2)(theta^(i+1)) S^(-1/2)(theta^i) (A6), and (N^-1 + S^-1(theta^(i+1)))^(-1/2)
    (N^-1 + S^-1(theta^i))^(1/2) (A17), for which "everything in the accept probability ratio R
    cancels" when the parameters are proposed from the exact marginal.
  - (inference) F = 1 satisfies (A8) with unit determinant and is exactly the translation move.
    The paper does not mention that choice. Filter (A17) is the full location-and-scale version,
    that is, collapsed sampling.
  - Validation: marginal means within 0.006 sigma and standard deviations within about 3 percent
    of the exact posterior (abstract, Table 1); correlation length about twice that of a direct
    sampler (Sect. 4).
  - Scope: (theta or C_l, s) only. Conclusions: "full physical marginalization over astrophysical
    foreground contamination will be straightforward (Eriksen et al. 2008; Planck Collaboration X
    2015)", which refers to combining with existing foreground sampling, not to applying the move
    to spectral indices. Spectral indices are not mentioned.
- **Later use.** Eskilt et al. 2023, arXiv 2306.15511 (SUM): the Racine move for cosmological
  parameters in Commander3 (Sect. 2.2.2, Eqs. 21 to 26), and "we did not sample any instrumental
  or foreground parameters" (Sect. 2.1). BeyondPlanck XI, arXiv 2208.14276 (SUM-part, first 100 000
  of 161 119 characters), Sect. 2.3: the joint step "has not yet been fully implemented in the
  latest version of the Commander code ... and it is therefore not used in the BeyondPlanck
  processing". Ducrocq, Chopin, Errard, Stompor 2021, arXiv 2111.07664 (SUM-part, first 100 000 of
  134 504 characters): non-centred Gibbs for C_l and interweaving, "also called ASIS in the
  statistical literature"; credits Jewell et al. with the non-centred Metropolis step; CMB and
  noise only. Millea, Anderes, Wandelt 2020, arXiv 2002.00965 (SUM-part, first 100 000 of 155 796
  characters): "mixing" reparametrization of fields before sampling r and A_phi, scalings and a
  lensing remap, no conditional-mean shift, described as "an extension of Racine et al. 2016".
  Morshed et al. 2024, arXiv 2405.18365 (SUM-part, first 100 000 of 136 427 characters): cites
  both as precedents for Metropolis within Gibbs; mixing-matrix elements by Metropolis on a
  profile-type likelihood, no deterministic joint move. A 21 cm paper, arXiv 2604.26885 (SUM-part,
  first 100 000 of 148 521 characters), Sect. 3.2: "Approaches such as deterministically rescaling
  the random draws can help avoid this issue (Racine et al., 2016), but we do not pursue this
  further." arXiv 2604.26890 (SUM-part, first 100 000 of 164 010 characters): fixed foreground
  basis, no spectral parameters sampled.

### B.2 Component separation: what is done when an index is proposed

- BeyondPlanck XIII, arXiv 2201.08188 (SUM-part this session, first 100 000 of 166 791 characters;
  the 2026-10-27 repo note records a full read). Sect. 4.1: Metropolis on the amplitude-marginal
  likelihood, then "we draw one sample from P(a|m_nu, beta) using Eq. (24)"; the ridge option
  "sets a equal to its maximum likelihood value for a given value of beta". Sect. 4.3: "for each
  proposed spectral parameter value, we estimate a new band monopole value conditioned on the
  input amplitude and the proposed parameter value." These are a collapsed step, a deterministic
  reset to the conditional mean, and a deterministic refit of a linear nuisance parameter. None
  keeps the fluctuation and moves the mean. No citation of Jewell 2009 or Racine 2016 there.
- Eriksen et al. 2008, arXiv 0709.1058, and the other Commander-family and parametric papers read
  in `2026-11-09_novelty_check_methods.md` (0811.4280, 0903.4311, 1103.2554, 1502.01588,
  1608.00551, 1409.0858, 1910.14170, 2509.11894): indices conditional on amplitudes, or
  amplitude-marginal, or joint gradient-based sampling. No deterministic amplitude shift (NOTE
  level: taken from that note, not re-read here).
- Not found: any CMB or radio component-separation paper in which a spectral-index proposal is
  accompanied by a deterministic shift or rescaling of the component amplitudes, and any paper
  that applies the Jewell or Racine move to a foreground parameter.

### B.3 Statistics relatives and names

- **Papaspiliopoulos, Roberts, Skold 2007**, Statist. Sci. 22, 59, DOI 10.1214/088342307000000014,
  arXiv 0708.3797. DIRECT, pages 1 to 6 and 10 to 15 of 15. The closest statistical statement.
  - Eq. (2.1): a reparametrization is a pair (X*, Theta) with "a function h such that X = h(X*,
    Theta, Y)", so h may depend on the data.
  - Example 2.1 (p. 2 to 3): "X_i = X_i* + A_i ... for vectors A ... which are allowed to depend
    only on theta. Remarkably, even when A is restricted to depend only linearly on theta, this
    family of reparameterizations is sufficiently flexible to include complete posterior
    orthogonalization of the state space of (X, Theta)."
  - Sect. 5.1.1 (p. 12 to 13): "a linear reparametrization X = v^(1/2)(Theta, Y) X* + m(Theta, Y)
    (5.1)"; "When correcting this parametrization for the presence of data, it is natural to
    consider v(Theta, Y) = Var(X | Theta, Y) and m(Theta, Y) = E(X | Theta, Y). Indeed, when also
    X | Theta, Y is location-scale, this choice will make X* independent of Theta." And: "(5.1)
    can be interpreted as a partially noncentered parametrization".
  - Sect. 4.1 (p. 11): when Theta is updated by Metropolis-Hastings under a reparametrization,
    h "only needs to be computed for the current value ... and a proposed new value Theta'".
  - (inference) The translation move is the Metropolis update of beta in (5.1) with m = E(a |
    beta, d) and v = 1. Example 2.1 covers the case where the conditional covariance does not
    depend on the parameter, in which the translation alone decorrelates completely.
- **Murray and Adams 2010**, arXiv 1006.0868. SUM. "Surrogate data" reparametrization for latent
  Gaussian models: with the whitened variable eta fixed, f' = L_R(theta') eta + m(theta', g)
  (Algorithm 4), with R_theta = (Sigma_theta^-1 + S_theta^-1)^-1 and m = R_theta S_theta^-1 g (Eq.
  9). Shift to the new conditional mean plus rescaling of the fluctuation, the same structure as
  Racine's (A17). A location-only version is not discussed. Cites Christensen et al. 2006.
- **Liu and Sabatti 2000**, Biometrika 87, 353, DOI 10.1093/biomet/87.2.353. ABS only: a
  generalised Gibbs sampler of "conditional moves along the traces of groups" that "encompasses
  a class of recently proposed tricks such as parameter expansion". (inference) By B.0 the
  translation move is a Metropolized move along the orbit of a one-parameter group, so this is
  the natural general reference. The theorem and its Metropolis variant were not read.
- **Liu and Wu 1999**, JASA 94, 1264, DOI 10.1080/01621459.1999.10473879; **Yu and Meng 2011**,
  JCGS 20, 531, DOI 10.1198/jcgs.2011.203main; **van Dyk and Park 2008**, JASA 103, 790, DOI
  10.1198/016214508000000409. BIB only (Crossref; publisher pages refused access). Their relation
  to the move is inference: parameter expansion moves the latent variables with the parameter
  through a working transformation; interweaving alternates a parameter update in the centred
  and in the non-centred coordinates, which is what "translation move, then Gibbs draw of a" does
  with a posterior-centred coordinate in place of the ancillary one; partial collapsing is the
  limit reached when the fluctuation is also rescaled. Read them before citing.

### B.4 Plain statement

A move of this form for spectral indices is **not found in print**. Its general form (deterministic
joint move with an involution condition, Jewell et al. 2009 Eqs. 14 to 15), its mean-recentred
form with a filter on the fluctuation (Racine et al. 2016 Eqs. 17 and A3 to A8), and its
statistical name (data-corrected, partially non-centred linear reparametrization, Papaspiliopoulos
et al. 2007 Eq. 5.1; posterior reparametrization, Murray and Adams 2010) are published. The
pure-translation member (identity filter) is not written out in the CMB papers read, and no
component-separation paper applies any member to (beta, a). A referee will expect Racine et al.
2016 and Papaspiliopoulos et al. 2007 to be cited as the origin of the idea.

## C. Nuances against the repo notes

1. `2026-11-09_novelty_check_methods.md`, Claim 1, lists as not found "a statement of what such a
   step targets". With the E035 identification of the legacy move as a refreshed pseudo-marginal
   chain, the statistics literature does state the qualitative answer under the name MCWM (not
   the target, approximate, possibly no invariant law), and Ball et al. 2003 give the temperature.
   That sentence should be narrowed to "not found for component separation".
2. The same note cites Andrieu and Roberts 2009 as BIB only with inferred relevance. It is now
   read (DIRECT, pages as above) and directly relevant: it names both variants and says what
   distinguishes them.
3. `FINDINGS.md` F1 says the width 0.79 k^(1/4) is "derived". Correct, but the underlying
   T = sigma sqrt(pi/8) has a 2003 antecedent that F1 and Paper A should cite.
4. `LITERATURE.md` line 60 calls 0807.0624 and 1512.06619 "rescaling moves". Racine's move is a
   mean-field recentring plus a rescaling of the fluctuation; only Jewell's is a pure rescaling.
5. `LITERATURE.md` line 61 describes 2111.07664 as auxiliary-variable and overrelaxed Gibbs. It
   also contains non-centred Gibbs and interweaving (ASIS), which is the part that bears on B.
6. No contradiction found with E035 itself.

## D. Searches run, and what was not found

arXiv API (`export.arxiv.org/api/query`), result counts in brackets:
- `abs:"Monte Carlo within Metropolis" OR abs:"noisy Metropolis" OR abs:"noisy Metropolis-Hastings"`
  [10; none states a law or a temperature].
- `(abs:Metropolis OR abs:"Monte Carlo") AND (abs:"noisy energy" OR abs:"noisy energies" OR
  abs:"uncertain energies" OR abs:"noisy acceptance" OR abs:"stochastic annealing")` [3].
- `abs:"spectral index" AND abs:amplitude AND (joint move, joint proposal, joint sampling,
  rescaling, reparameterization) AND (abs:foreground OR abs:CMB)` [1: 0709.1058].
- `(abs:"component separation" OR abs:foreground) AND (non-centered, noncentered, non-centred,
  interweaving, "partially collapsed", "parameter expansion", "joint move", "rescaling move")`
  [1, irrelevant].
- `all:"spectral index" AND all:Gibbs AND all:Metropolis AND all:amplitude AND (joint OR marginal
  OR collapsed)` [0].
- `(abs:"component separation" OR abs:foreground) AND (abs:Gibbs OR abs:Metropolis) AND (spectral
  index, spectral indices, spectral parameters, mixing matrix)` [11, none with a deterministic
  amplitude move].

INSPIRE: citers of Racine et al. 2016 (`refersto:recid:1410785`) [20] and of Jewell et al. 2009
(`refersto:recid:789912`) [42], all titles scanned; opened: 2306.15511, 2208.14276, 2405.18365,
2111.07664, 2002.00965, 2604.26885. Not opened: 1706.08971, 2305.16134, 2210.13260, 1607.00008,
1610.06673, 2509.13435, 2510.26691, 2505.24829, 2610.18036. (The `refersto:arxiv:` form returned
an unfiltered list through the fetch tool and was discarded.)

Semantic Scholar citation list of cond-mat/0301179 [13 returned; one MCMC paper, 1205.6857].

WebSearch (8 queries): the sqrt(pi/8) constant with noisy Metropolis and effective temperature;
"stochastic annealing" with noisy energy and effective temperature; MCWM with tempered posterior
and width inflation; Beaumont 2003 full text; Branke et al. 2008; Liu and Sabatti 2000 (two
queries); Yu and Meng 2011; a deterministic amplitude shift with spectral-index proposals.

Not found in any of the above: (a) an explicit stationary law or temperature for MCWM or noisy
MH in a statistics paper that uses those names; (b) an effective temperature for the Metropolis
rule (as opposed to the sign rule) with Gaussian noise at finite sigma; (c) the chi-square case;
(d) any astronomy or CMB occurrence of either; (e) a deterministic amplitude shift or rescaling
tied to a spectral-index proposal in component separation; (f) the identity-filter member of
Racine's family written out anywhere.

## E. Gaps and suggested follow-up

- Beaumont 2003 and O'Neill et al. 2000 were not read in full. If Paper A quotes them for MCWM,
  read them; in particular check whether Beaumont remarks on MCWM posterior width.
- Branke et al. 2008 full text (first cited use of the sign rule according to Bardenet et al.):
  check whether it restates Ball's Eq. (7) or gives a finite-temperature version. Also Fink's
  thesis and Bowler's 2001 thesis (Ball et al. ref. [14]), and Ball, Branke, Meisel 2018, "Optimal
  Sampling for Simulated Annealing Under Noise" (DOI 10.1287/ijoc.2017.0774, title seen only).
- Liu and Sabatti 2000, Liu and Wu 1999, Yu and Meng 2011, van Dyk and Park 2008, and Christensen,
  Roberts, Skold 2006 (J. Comput. Graph. Statist. 15, 1, as listed in 0708.3797 ref. [6]) need a
  full read through a library before the paper leans on them.
- Unread remainders: 1505.02827 beyond PDF page 23 and before page 15 (SUM-part only); 0903.5480
  pages 11 to 17; 0708.3797 pages 7 to 9; 0807.0624 appendices (SUM only); every SUM-part range
  stated above.
- Wu, Wang, Wong (1908.02910) say that "a connection between tempering and subsample variance was
  also mentioned in [5]"; the reference number was not resolved. Likely Bardenet et al.; check.
- A two-line numerical check worth adding to E035 (inference): the sign rule with Gaussian noise
  should reproduce Ball's correction term, that is, the effective temperature falling as
  1 - 0.018 (Delta E / T)^2 with the proposal size. That would tie E035 P4 at one step per
  refresh to the published formula.
- Citation lists of 1512.06619 and 0807.0624 were scanned by title; nine citers were not opened.

## F. Sources

| Source | What was read | Finding | Bears on |
|---|---|---|---|
| Ball, Fink, Bowler 2003, cond-mat/0301179 | DIRECT, all 5 pages | Sign rule with Gaussian noise: T ~= sqrt(pi/8) sigma (1 - 0.018 (Delta E/T)^2), Eq. 7; detailed balance "not strictly achieved" | A (decisive) |
| Bardenet, Doucet, Holmes 2017, 1505.02827 | DIRECT pp. 15 to 23; SUM-part otherwise | Sign rule is close to Barker targeting pi^beta, beta = lambda sqrt(t)/(n sigma), lambda = 1.702 (Sect. 6.3) | A (close) |
| Andrieu, Roberts 2009, 0903.5480 | DIRECT pp. 1 to 10, 18 to 29 | MCWM defined; target "typically not the invariant distribution"; existence "not obvious"; Theorem 9 bound; no explicit law | A |
| Ceperley, Dewing 1999, physics/9812035 | DIRECT pp. 1 to 4; SUM rest | Uncorrected noise "will bias the asymptotic distribution"; penalty rule Eq. 13; no temperature | A |
| Medina-Aguayo, Lee, Roberts 2016, 1503.07066 | SUM, all | Noisy MH is MCWM; may be transient; invariant law "unknown"; converges to target with rates | A |
| Alquier, Friel, Everitt, Boland 2016, 1403.5496 | SUM-part | Total-variation bounds only; MCWM against GIMH stated | A |
| Nicholls, Fox, Watt 2012, 1205.6857 | SUM | Naive algorithm "does not in general target pi"; "slightly over-dispersed"; cites Ball et al. | A |
| Korattikara, Chen, Welling 2014, 1304.5299 | SUM | Bound on stationary error, no law | A |
| Seita, Pan, Chen, Canny 2017, 1610.06848 | SUM | Barker test with logistic noise; Gaussian sd 1.7 matches the logistic; exact correction only | A |
| Wu, Wang, Wong, 1908.02910; Li, Wong 2017, 1707.09705 | SUM-part; SUM | Tempered stationary law designed in by scaling, not from acceptance noise | A (different mechanism) |
| Bardenet, Doucet, Holmes 2014, PMLR 32 | ABS | "perturbed version of the target" | A |
| Kennedy, Kuti 1985, DOI 10.1103/PhysRevLett.54.2473 | ABS (paraphrase) | Tolerates unbiased noise in the ratio; no temperature found | A |
| Lin, Liu, Sloan, hep-lat/9905033; Bakeyev, de Forcrand, hep-lat/0008006 | SUM | Exact noisy Monte Carlo; uncorrected case not analysed | A |
| Pierleoni, Ceperley, physics/0510254 | SUM-part | Uncorrected results biased, bias grows with noise; no temperature | A |
| Doll, Freeman 1994, chem-ph/9405002 | SUM | Gaussian noise in the Metropolis ratio; no effective law | A |
| Sanderski, Gentinetta, Carleo 2026, 2609.01204 | ABS | "a noisy energy biases the sampled distribution" | A |
| Llorente et al., 2108.00490; Medina-Aguayo, Rudolf, Schweizer, 1809.09547 | SUM-part; ABS | No explicit MCWM law; bounds | A |
| Branke, Meisel, Schmidt 2008, DOI 10.1007/s10732-007-9058-7 | BIB and abstract | Improved stochastic annealing; content not read | A |
| Beaumont 2003, DOI 10.1093/genetics/164.3.1139 | ABS | GIMH origin (per Andrieu and Roberts); not read | A |
| O'Neill et al. 2000, DOI 10.1111/1467-9876.00210 | BIB | MCWM term (per Andrieu and Roberts); not read | A |
| Racine, Jewell, Eriksen, Wehus 2016, 1512.06619 | DIRECT, all 9 pages | s' = s_hat' + F (s - s_hat), deterministic, filter family Eqs. A3 to A8, A17; no foreground use | B (closest CMB) |
| Jewell et al. 2009, 0807.0624 | DIRECT pp. 3 to 5; SUM rest | Deterministic joint move, Eqs. 14 to 16; pure rescaling; change-of-variables remark | B |
| Papaspiliopoulos, Roberts, Skold 2007, 0708.3797 | DIRECT pp. 1 to 6, 10 to 15 | X = v^(1/2) X* + m(Theta, Y), m = E(X | Theta, Y), "partially noncentered" (Eq. 5.1); Example 2.1 | B (closest statistics) |
| Murray, Adams 2010, 1006.0868 | SUM | Surrogate-data reparametrization: shift to the conditional mean and rescale | B |
| Eskilt et al. 2023, 2306.15511 | SUM | Racine move in Commander3; no foreground parameters sampled | B |
| BeyondPlanck XI, 2208.14276 | SUM-part | Joint (a, C_l) step not implemented, not used | B |
| BeyondPlanck XIII, 2201.08188 | SUM-part | Marginal step then amplitude draw; ridge reset; monopoles refit per proposal | B |
| Ducrocq et al. 2021, 2111.07664 | SUM-part | Non-centred Gibbs and ASIS for C_l; no foregrounds | B |
| Millea, Anderes, Wandelt 2020, 2002.00965 | SUM-part | Mixing reparametrization, scalings only | B |
| Morshed et al. 2024, 2405.18365 | SUM-part | Mixing matrix by Metropolis within Gibbs, no joint move | B |
| 2604.26885; 2604.26890 | SUM-part | Rescaling named as an option, not used; no spectral parameters sampled | B |
| Liu, Sabatti 2000, DOI 10.1093/biomet/87.2.353 | ABS | Generalised Gibbs: "conditional moves along the traces of groups" | B |
| Liu, Wu 1999; Yu, Meng 2011; van Dyk, Park 2008 (DOIs in B.3) | BIB | Not read; relation is inference | B |
