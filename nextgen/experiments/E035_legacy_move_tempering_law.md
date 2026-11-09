# E035: the legacy re-solve move as a refreshed pseudo-marginal chain, and its tempering law

- **Date started / finished:** 2026-11-09 22:52 / 2026-11-09 22:57
- **Branch and commit:** `nextgen` at 349759b
- **Open question addressed:** Q1 (closed) and Paper A claim 1; referee objections 8 and the disputed
  k^(1/4) arithmetic of F1 (`../reviews/2026-11-09_referee_report.md`, Sect. 2.1 and 4.1 item 5)
- **Validation level:** L1 (toy chains against a derived stationary law, with exact controls)

## Question

F1 says the legacy `resolve_chisq` move "is not a valid MCMC step" and that its width inflation
"grows roughly as k^(1/4)", an empirical fit that the referee could not reproduce from the record
(14.1 to 16.0 against the quoted 18.5 for E014). Is there a derivation, with a coefficient, of what
the legacy sampler converges to? What changes: Paper A claim 1 becomes a statement of the legacy
chain's actual stationary law, with the real-map number (16.5 and 17.1 at k = 190825) as a test of
it, and the words "not a valid step" are replaced by what is exactly true.

## Target and priors

Setting: flat amplitude priors, one index beta with a symmetric random-walk proposal, k amplitude
degrees of freedom re-solved per proposal. Write chi2(beta, a) = chi2_min(beta) + X, with
X = (a - a_hat)^T A (a - a_hat). For a ~ p(a | beta, d) = N(a_hat, A^-1), X ~ chi^2_k and does not
depend on beta.

1. *The move alone is a pseudo-marginal Metropolis-Hastings kernel* (Andrieu and Roberts 2009). The
   proposal is (beta', a') with a' ~ N(a_hat(beta'), A(beta')^-1) and the acceptance is
   exp(-(chi2' - chi2)/2). That is the Hastings ratio of the extended target
   pi(beta, a) propto exp(-chi2_min(beta)/2) det(A)^(1/2) exp(-(a - a_hat)^T A (a - a_hat)),
   whose beta marginal is the profile (ridge) posterior exp(-chi2_min(beta)/2) and whose amplitude
   conditional is N(a_hat, (2A)^-1), half the posterior covariance. Under it X ~ Gamma(k/2, scale 1),
   mean k/2.
2. *The Gibbs sweep refreshes the amplitudes from the wrong conditional.* Between index updates
   the amplitude group draws a ~ N(a_hat, A^-1), so X is reset to a chi^2_k variate (mean k)
   instead of the Gamma(k/2, 1) the extended target requires. The composite chain leaves neither
   target invariant.
3. *Stationary law of the composite chain with one index step per refresh.* The acceptance
   probability of a move with Delta = chi2_min(beta') - chi2_min(beta), averaged over the fresh X
   and X', is E min(1, exp(-(Delta + Z)/2)) with Z = X' - X, variance 4k. For large k this is
   Phi(-Delta / (2 sqrt k)) + O(k^(-1/2)). For steps with abs(Delta) much smaller than 2 sqrt(k)
   the ratio of forward to backward acceptance is exp(-4 Delta / (2 sqrt(2 pi k))), so the chain is
   reversible with respect to exp(-chi2_min(beta) / (2 T)) with temperature
   T = sqrt(2 pi k) / 4 = 0.627 sqrt(k). The index posterior is the ridge posterior with its width
   multiplied by sqrt(T) = (2 pi)^(1/4) / 2 x k^(1/4) = 0.7916 k^(1/4), and with the ridge mean.

## Setup

Toy: chi2_min(beta) = beta^2 (ridge posterior N(0, 1)), X drawn exactly from chi^2_k. 4000
independent chains, vectorized, 100 burn-in and 300 kept iterations (one iteration = one refresh
plus `numstep` index steps), started at beta = 0. Convergence check: sd over the first and second
halves of the kept iterations. One thread, numpy, seed fixed.

## Reference

**Predictions, written before the first run:**

- P1 (derived): one step per refresh, proposal sd 0.3 times the predicted width: sd(beta) =
  0.7916 k^(1/4) within 3 percent for k = 100, 1000, 1e4, 1e5, 190825, and Gaussian (excess
  kurtosis within 0.1). At k = 10 the large-k expansion is not expected to hold; reported.
- P2 (exact control): no refresh at all (the pure pseudo-marginal chain): sd(beta) = 1 within Monte
  Carlo error and mean X = k/2, at k = 4, 8, 16 where the chain still accepts (acceptance of order
  exp(-k/4)).
- P3 (exact control): refresh from the extended target's own conditional, X ~ Gamma(k/2, 1):
  sd(beta) = 1 at k = 4, 8, 16.
- P4 (no derived prediction): larger proposals (1 and 2.4 predicted widths) and 5 or 20 steps per
  refresh. Expectation: more steps per refresh give a narrower law, because X ratchets down
  between refreshes and the chain moves toward the exact pseudo-marginal kernel. The E014 legacy
  chain had k = 190825 (38165 pixels, five components), 20 steps per refresh, proposals of 2.4
  legacy widths, two indices proposed jointly, and measured 16.5 and 17.1 (ESS 74, so about 8
  percent each). If the toy at those settings, in one and in two dimensions, lands outside 13.5 to
  20, the mechanism as modelled here does not explain the real-map number.

## Command

```bash
cd /Users/trygvels/code/Commander4
OMP_NUM_THREADS=1 .venv/bin/python nextgen/scripts/e035_legacy_tempering.py > nextgen/experiments/E035/toy.txt
```

## Result

Toy run 22:52:55 to 22:54:02 on one thread (`E035/toy.txt`). sd in units of the ridge sd; Monte
Carlo errors 0.2 to 0.5 percent; the two halves of every run agree to 1 percent.

**P1, one step per refresh, small proposals: pass for k of 1000 and above.**

| k | law 0.7916 k^(1/4) | toy sd | ratio | excess kurtosis |
|---|---|---|---|---|
| 10 | 1.408 | 1.577 | 1.120 | -0.01 |
| 100 | 2.503 | 2.584 | 1.032 | -0.02 |
| 1000 | 4.452 | 4.463 | 1.003 | -0.04 |
| 10000 | 7.916 | 7.870 | 0.994 | -0.01 |
| 100000 | 14.077 | 14.027 | 0.996 | -0.02 |
| 190825 | 16.545 | 16.435 | 0.993 | -0.01 |

The law holds to 0.7 percent from k = 1000 up and the stationary law is Gaussian. At k = 100 the
toy is 3.2 percent above the law, 0.2 points outside the 3 percent I wrote down, so P1 fails at
k = 100 as stated; the O(k^(-1/2)) correction is visible there and at k = 10 (12 percent).

**P2 and P3, the exact controls: pass.** Without any refresh the chain returns the ridge posterior
(sd 1.000, 0.999, 1.001 at k = 4, 8, 16) with mean X of 2.000, 3.999, 8.001, which is k/2. With
the refresh drawn from the extended target's own conditional, Gamma(k/2, 1), it returns sd 1.000,
1.000, 1.001. Acceptance falls from 0.28 to 0.10 between k = 4 and 16, as exp(-k/4) predicts.

**P4, more steps per refresh and wider proposals.** Ratio of the toy sd to the one-step law:

| steps per refresh | proposal 0.3 W | 1.0 W | 2.4 W | two indices, 1.0 W | two indices, 2.4 W |
|---|---|---|---|---|---|
| 1 | 0.99 to 1.00 | 0.97 to 0.98 | 0.96 to 0.97 | 0.97 to 0.98 | 0.95 to 0.96 |
| 5 | 0.81 | 0.82 | 0.86 | 0.83 | 0.90 |
| 20 | 0.68 | 0.69 to 0.70 | 0.74 to 0.75 | 0.71 to 0.72 | 0.80 to 0.81 |

The same ratios at k = 1000 and k = 190825, so the k^(1/4) scaling holds at any fixed number of
steps and only the coefficient changes: 0.79 at one step, 0.64 at five, 0.54 at twenty for small
proposals. The expectation (narrower with more steps) holds.

**The toy at the settings of the earlier measurements** (`E035/toy_at_earlier_settings.txt`,
`E035/toy_e014_recorded_proposals.txt`; run after P1 to P4 were read, so not pre-registered):

| measurement | k | steps | toy | measured |
|---|---|---|---|---|
| E002 prototype, 256 pixels, three components, two indices | 768 | 1 | 4.03 to 4.08 | 3.9 |
| E002 prototype, 1024 pixels | 3072 | 1 | 5.66 to 5.72 | 5.7 |
| E005 unit test, 150 pixels, three components | 450 | 4 | 3.04 | 3.2 |
| E005 unit test, 600 pixels | 1800 | 4 | 4.28 | 4.4 |
| E014 real maps, 38165 pixels, five components, recorded proposals | 190825 | 20 | 13.1 and 13.6 | 17.1 and 16.5 |

The five simulated numbers agree with the toy to 5 percent. The real-map number does not: the toy
at E014's settings gives 13.1 to 13.6 (12.5 in one dimension), below the 13.5 to 20 band fixed in
advance in one of the two indices and in one dimension, and 22 to 31 percent below the measured
widths, which carry about 8 percent error each at ESS 74. The toy reproduces the chain's
acceptance (0.067 against 0.073). The measured widths sit at the one-step law (16.5).

## Negative control

P2 and P3 above. Removing the refresh or drawing it from the right conditional restores the ridge
posterior exactly, so the inflation is caused by the refresh from p(a | beta, d) and by nothing
else in the toy.

## Conclusion

1. **The legacy move is a valid kernel for a different target, and the Gibbs sweep breaks it.**
   Taken alone, "propose beta', draw a' from p(a | beta', d), accept on chi-square at the draws"
   is the pseudo-marginal Metropolis-Hastings kernel of an extended target whose index marginal
   is the profile (ridge) posterior and whose amplitudes have half their posterior covariance.
   The amplitude step of the Gibbs sweep then redraws the amplitudes with their full posterior
   covariance, which is not that target's conditional. F1's "not a valid MCMC step" is replaced by
   this statement. It also explains why the legacy means match the ridge posterior.
2. **With one index step per amplitude draw the chain converges to the ridge posterior tempered
   by T = sqrt(2 pi k) / 4**, so the index width is 0.7916 k^(1/4) times the ridge width, with k
   the number of amplitude degrees of freedom re-solved per proposal. Derived above, confirmed to
   0.7 percent for k from 1000 to 190825. More index steps per amplitude draw lower the
   coefficient (0.64 at five, 0.54 at twenty) and leave the k^(1/4) scaling.
3. The correct-looking repair of drawing the refresh from N(a_hat, (2A)^-1) is exact and useless:
   acceptance falls as exp(-k/4). The marginal or ridge acceptance of F2 is the fix.
4. **The real-map width is not explained by the toy.** E014 measured 16.5 and 17.1 where the toy
   at its settings gives 13.1 to 13.6. Candidates, none tested: the legacy chi-square is evaluated
   on the alm-represented amplitudes (F7), which adds acceptance noise the toy does not have and
   would need to multiply the variance of the chi-square difference by 2.4; the width of a chain
   with ESS 74 that a tuned rerun could not traverse; correlated indices. Until one of them is
   shown, Paper A quotes the law for the simulated cases and gives the real-map number as
   measured, 25 percent above it.

What this does not establish: anything for proposals in more than two dimensions, for the CG
amplitude path, where the re-solve is a truncated solve (audit A1), or for other codes.

Changes: F1 rewritten around items 1 and 2 with the old empirical fit kept as history; the
novelty report's sentence that the move "has no stated target" is corrected in the program review.

## Follow-ups

- The legacy mode in Commander4 at small nside with 1 and 20 steps per refresh, with and without
  the alm storage of the per-pixel amplitudes, to test the first candidate of item 4. Needs cores,
  so after E032.
- The O(k^(-1/2)) correction, if Paper A wants the law below k = 1000.

## Negative control

P2 and P3 are the controls: the same code with the refresh removed or corrected must return the
ridge posterior exactly.

## Conclusion

Pending.

## Follow-ups
