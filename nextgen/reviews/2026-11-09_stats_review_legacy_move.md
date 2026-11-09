# Stats review: the legacy `resolve_chisq` move (`_run_resolve_chisq`)

Reviewer: stats-reviewer subagent. Started 2026-11-09 23:04 CET; sections 1 to 6 written
2026-11-09 23:18 CET (small corrections 23:20), before opening `nextgen/experiments/E035*`,
`nextgen/FINDINGS.md` or any other file under `nextgen/reviews/`. Sections 7 to 9 (comparison with
E035, test gaps, verdict) were added 23:36 CET, after reading them.

Code under review: `src/commander4/compsep/mcmc.py:260-348` (`_run_resolve_chisq`), called from
`src/commander4/compsep/processing.py:862-909` (amplitude groups once per iteration, then `numstep`
moves per MCMC group). No repository file other than this report was modified.

Notation: chi2(beta, a) = chi2_min(beta) + X, X = (a - a_hat)^T A (a - a_hat), A = M^T N^-1 M,
k = number of real amplitude degrees of freedom that the re-solve redraws, p(beta) the prior on beta,
q a symmetric random walk. A fresh exact draw a ~ N(a_hat(beta), A(beta)^-1) has X ~ chi2_k for
every beta. "Profile posterior" means P(beta) = p(beta) exp(-chi2_min(beta)/2).

## 1. The legacy move taken alone (question 1)

Proposal density from (beta, a): q(beta'|beta) g(a'|beta'), with
g(a|beta) = |A(beta)|^(1/2) (2 pi)^(-k/2) exp(-X(beta, a)/2). Acceptance in the code
(`mcmc.py:319-326`): min(1, [p(beta')/p(beta)] exp(-(chi2(beta', a') - chi2(beta, a))/2)).

Matching this to the Metropolis-Hastings ratio pi(beta', a') g(a|beta) / (pi(beta, a) g(a'|beta'))
requires pi(beta, a) / g(a|beta) proportional to p(beta) exp(-chi2(beta, a)/2). Hence the move is a
valid MH kernel, reversible with respect to

    pi_K(beta, a)  ∝  p(beta) |A(beta)|^(1/2) exp(-chi2(beta, a) + chi2_min(beta)/2)
                   =  p(beta) exp(-chi2_min(beta)/2) * |A(beta)|^(1/2) exp(-(a - a_hat)^T A (a - a_hat)).

- Marginal of beta: integrating a gives pi^(k/2) |A|^(-1/2), which cancels |A|^(1/2) exactly:
  pi_K(beta) ∝ p(beta) exp(-chi2_min(beta)/2). This is the profile ("ridge") posterior, with no
  ln det A term.
- Conditional of a: N(a_hat(beta), (2 A(beta))^-1). Half the covariance of the true conditional;
  X = chi2_k / 2 (mean k/2, variance k/2) instead of chi2_k.
- So it is not the joint posterior p(beta) exp(-chi2/2) (wrong amplitude scatter by 1/sqrt 2, and
  beta marginal lacks |A|^(-1/2)), but it is a proper target.

Equivalent reading: W = exp(-chi2(beta, a)/2) with a ~ g(.|beta) is an unbiased estimator of
2^(-k/2) exp(-chi2_min(beta)/2), a beta-independent constant times the profile likelihood, and the
move keeps the estimate of the current state. That is the pseudo-marginal construction (Andrieu
and Roberts 2009), which is exact for the profile posterior. Var(ln W) = k/2, so it is exact but
mixes exponentially slowly in k (section 4).

## 2. Composite chain, one legacy move per amplitude redraw (question 2)

The redraw from p(a|beta, d) resets X ~ chi2_k, independent of the past. The redraw preserves laws
of the form mu(beta) N(a_hat, A^-1); the legacy move preserves pi_K with N(a_hat, (2A)^-1). They
share no invariant law, so the composite targets neither.

With one move per redraw the beta chain is Markov by itself, with effective acceptance

    alpha(D) = E[min(1, exp(-(D + Z)/2))],   D = chi2_min(beta') - chi2_min(beta) - 2 ln[p(beta')/p(beta)],
    Z = X' - X,  X, X' iid chi2_k (symmetric, variance 4k).

This is Monte Carlo within Metropolis (both the proposed and the current estimate refreshed), which
is known not to be exact. alpha depends on beta only through D: det A never enters.

Small-step (diffusion) limit. For a step eps with variance s^2, drift = s^2 alpha'(0) dV/dbeta and
diffusion = s^2 alpha(0), with V the D-potential. The stationary density is
exp(2 alpha'(0) V / alpha(0)). With c_k = E[exp(-Z/2); Z > 0]: alpha(0) = 1/2 + c_k and
alpha'(0) = -c_k/2. Therefore

    pi(beta) ∝ [p(beta) exp(-chi2_min(beta)/2)]^(1/T),     T(k) = 1/2 + 1/(4 c_k),

    c_k = 2^-(k+1) Gamma(k+1) / Gamma(k/2+1)^2   (= 2^-(k+1) C(k, k/2); c_0 = 1/2, c_1 = 1/pi, c_2 = 1/4),

    T(k) = sqrt(pi k / 8) (1 + 1/(4k) + ...) + 1/2   ≈  0.6267 sqrt(k) + 0.5.

The profile posterior (prior included) is tempered by T. For a Gaussian profile of width
sigma_prof the stationary width is sigma_prof sqrt(T) ≈ sigma_prof (pi k/8)^(1/4)
= 0.7916 k^(1/4) sigma_prof at large k. T(0) = 1 recovers plain Metropolis. The stationary
acceptance rate is 1/2 + c_k, tending to 1/2 and independent of the (small) step size.

General form, useful when the acceptance noise is not exactly chi2_k: for any symmetric continuous
Z, T = 1/2 + 1/(4 E[exp(-Z/2); Z>0]) ≈ 1/2 + 1/(8 f_Z(0)); for Gaussian Z,
T ≈ 1/2 + sqrt(2 pi) sd(Z)/8 = 1/2 + 0.4431 sd(chi2 at a fresh draw). The width scales as the square
root of the standard deviation of the acceptance noise.

What it rests on:
1. The amplitude draws are exact and the acceptance chi2 is the chi2 of the same linear model, so
   X ~ chi2_k for every beta. Exact in the stated model.
2. The diffusion limit: proposal step small compared with the stationary width. In that limit the
   formula T = 1/2 + 1/(4 c_k) is exact for every k and every shape of chi2_min. The closed form of
   c_k is exact; only the last line (sqrt(pi k/8)) uses large k (Stirling, equivalently the CLT for
   Z). For finite steps the kernel is not reversible with respect to any function of chi2_min
   unless alpha(D)/alpha(-D) is exactly exponential in D, which it is not; see section 3.

## 3. Several moves per redraw, and large steps (question 3)

### 3a. Large steps, one move per redraw
At large k, alpha(D) → Phi(-D / (2 sqrt k)): a probit acceptance. It is close to the Barker rule
1/(1 + exp(D/(2T))) for the tempered target (same slope at D = 0), and Barker is exact for any step.
The probit has lighter tails, ln[alpha(D)/alpha(-D)] = -(D/(2T))(1 + 0.045 x^2 + ...) with
x = D/(2 sqrt k), so far proposals are rejected slightly more often and the law narrows a little and
becomes mildly sub-Gaussian. For chi2_min = beta^2 the scaling limit (u = beta / k^(1/4)) is a
universal function of eps = step / sqrt(T): Var/T falls from 1 (eps → 0) as about 1 - 0.15 eps^2 at
small eps, to a plateau of 0.918 at eps ≥ 4 (0.907 at k = 100), with kurtosis 2.83. So the law of section
2 holds within 8 percent in variance (4 percent in width) for any step size, still ∝ sqrt(k).

### 3b. n moves per redraw
Between redraws the retained X is no longer fresh: an accepted X' is biased low, the next proposal
must beat it. In the scaling limit (k large, xi = (X - k)/sqrt(2k) standard normal, u = beta/k^(1/4))
the move is deterministic: accept iff u'^2 + sqrt2 xi' < u^2 + sqrt2 xi. For small steps the
accepted moves are the lower records of the iid sequence xi_0 (redraw), xi_1, ..., xi_n. The j-th
proposal is accepted with probability 1/(j+1), so the number of accepted moves per iteration is
N_acc(n) = H_(n+1) - 1 (0.5, 0.83, 1.45, 2.02, 3.03, 4.20 for n = 1, 2, 5, 10, 30, 100): it grows
only as ln n and the acceptance fraction N_acc/n is independent of the step size. The drift per
proposal is set by phi(running minimum), which gives

    Var_n(beta) / sigma_prof^2  ≈  sqrt(k) N_acc(n) / (2 sqrt2 N_phi(n)),
    N_phi(n) = sum_{j=1..n} E[phi(min of j iid N(0,1))] = sum_j j ∫ phi(t)^2 (1 - Phi(t))^(j-1) dt.

n = 1 gives sqrt(pi k/8) again. Ratio R_n = T_n/T_1: 1, 0.833 (= 5/6 exactly), 0.744, 0.647, 0.543,
0.467, 0.432, 0.396, 0.357, 0.311, 0.276 for n = 1, 2, 3, 5, 10, 20, 30, 50, 100, 300, 1000. The law
stays a power of the profile posterior (kurtosis 3) and stays ∝ sqrt(k); more moves narrow it only
as about 1/sqrt(ln n). Relative corrections are O(k^(-1/2)) and grow with n (measured: +8 percent
at n = 1 to -8 percent at n = 100 for k = 100; under 1.5 percent at k = 2e5).

As n → ∞ at fixed k the block converges to pi_K and the beta marginal to the profile posterior, but
only after of order exp(0.06 k) moves (section 4), which is never reached at map-level k.

With n > 1 and large steps the narrowing is partly undone: most proposals fail on D alone and do not
consume a record, so the effective n is about n/eps. Measured at k = 1e5, Var/T_1 for n = 10 goes
0.546, 0.548, 0.568, 0.628, 0.710, 0.795 at eps = 0.25, 0.5, 1, 2, 4, 8 (eps = step/sqrt(T_1)), and
tends to the n = 1 plateau 0.918 as eps → ∞. Full table in section 5.

## 4. A redraw that makes the composite exact for the profile posterior (question 4)

Yes. Replace the amplitude redraw by a draw from the conditional of pi_K,
a ~ N(a_hat(beta), (2 A(beta))^-1), that is a = a_hat + (fluctuation)/sqrt 2, while the legacy move
still proposes a' at full covariance. The redraw is then a Gibbs step of pi_K and the legacy move is
MH for pi_K, so any composite of them leaves pi_K invariant, for any numstep and any step size, and
beta has the profile posterior exactly. (The redraw is not even needed: the legacy move alone has
the same target.) The stored amplitudes then have the wrong scatter and need one true-conditional
redraw for output, taken outside the chain state.

Not usable at large k. At stationarity the current X is Gamma(k/2, 1) (mean k/2) and the proposed
X' is chi2_k (mean k). With beta frozen the move is an independence sampler with acceptance rate
2 I_(1/3)(k/2, k/2) (regularized incomplete beta), which decays as (8/9)^(k/2) = exp(-0.059 k):
0.18 at k = 16, 6e-3 at k = 64, 6e-4 at k = 100, 2e-27 at k = 1000, 4e-258 at k = 1e4. Equivalent
statement: the pseudo-marginal estimator has Var(ln W) = k/2, and pseudo-marginal chains need that
to be of order one. Usable for k up to a few tens only. At map-level k the correct fix is to remove
the amplitudes from the acceptance (the ridge/marginal options), not to repair the redraw.

## 5. Numerical check

Own toy, numpy only, single thread: chi2_min = beta^2 (profile variance 1), X ~ chi2_k, flat prior.
Scripts: `/private/tmp/claude-501/-Users-trygvels-code-Commander4/e8b23887-d1a3-4d04-b2ac-228d7b739a81/scratchpad/sr/`
(`sr_lib.py`, `run_selfcheck.py`, `run_A.py` to `run_E.py`); the scratchpad is not permanent, so the
core is listed in the appendix. Two independent references for n = 1:
(i) simulation, 1e4 parallel chains, half started at beta = 0 and half over-dispersed, burn-in of
20 or more autocorrelation times, error bar from the spread of per-chain time averages;
(ii) an exact-kernel reference: alpha(D) by FFT convolution of the chi2_k density (no Gaussian
approximation), transition matrix on a beta grid of 1201 points, stationary vector by direct solve.
The closed form of c_k agrees with the FFT quadrature to 6 digits for k = 10 to 2e5.

### 5a. n = 1, small step (step = 0.1 sqrt(T)), 6000 kept iterations per chain

| k | T = 1/2 + 1/(4 c_k) | sqrt(pi k/8) | exact-kernel grid | simulated Var(beta) | (sim - T)/se | acceptance (1/2 + c_k) |
|---|---|---|---|---|---|---|
| 10 | 2.532 | 1.982 | 2.524 | 2.525 ± 0.008 | -0.8 | 0.621 (0.623) |
| 100 | 6.782 | 6.267 | 6.769 | 6.751 ± 0.024 | -1.3 | 0.538 (0.540) |
| 1000 | 20.32 | 19.82 | 20.29 | 20.26 ± 0.07 | -0.8 | 0.511 (0.513) |
| 1e4 | 63.17 | 62.67 | 63.07 | 63.30 ± 0.23 | +0.6 | 0.503 (0.504) |
| 1e5 | 198.67 | 198.17 | 198.35 | 199.39 ± 0.72 | +1.0 | 0.500 (0.501) |
| 2e5 | 280.75 | 280.25 | 280.31 | 280.22 ± 1.01 | -0.5 | 0.500 (0.501) |

Also at step 0.2 sqrt(T): k = 30: T 3.961, grid 3.926, sim 3.927 ± 0.012; k = 300: T 11.363, grid
11.285, sim 11.319 ± 0.035. Kurtosis 2.98 to 3.01 throughout (Gaussian, as a tempered Gaussian
must be). Under- and over-dispersed starts agree. The +1/2 term is needed at small k (k = 10: sim
is 27 percent above sqrt(pi k/8) and 0.3 percent below T). The grid reference converges to T as the
step shrinks: Var/T = 0.980, 0.994, 0.9984, 0.9996 at step/sqrt(T) = 0.4, 0.2, 0.1, 0.05 (k = 1e5;
same pattern at k = 10, 100, 1000), i.e. the residual 0.2 percent in the table is the finite step.
Negative controls: k = 0 (plain Metropolis, T(0) = 1) run through the same simulator gives
0.9997 ± 0.0002 pooled over 12 seeds; and the half-covariance redraw at k = 1 to 64 gives 1.000
(5e), so the excess is the fresh full-covariance redraw and nothing else.

### 5b. n = 1, step scan (eps = step/sqrt(T)); Var/T from the exact-kernel grid

| eps | 0.25 | 0.5 | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|---|---|
| k = 100 | 0.989 | 0.967 | 0.934 | 0.912 | 0.908 | 0.907 | 0.907 |
| k = 1e4 | 0.991 | 0.973 | 0.943 | 0.922 | 0.918 | 0.918 | 0.917 |
| k = 2e5 | 0.991 | 0.973 | 0.943 | 0.924 | 0.919 | 0.919 | 0.919 |
| acceptance (k = 1e4) | 0.496 | 0.476 | 0.417 | 0.305 | 0.182 | 0.096 | 0.049 |
| kurtosis (k = 1e4) | 2.98 | 2.95 | 2.90 | 2.84 | 2.83 | 2.83 | 2.84 |

Simulation (6000 chains) agrees with the grid in all 21 cells within 2 se (se 0.1 to 0.4 percent).

### 5c. n moves per redraw, small step (0.15 sqrt(V_n)); V_n = sqrt(k) N_acc/(2 sqrt2 N_phi)

| n | R_n predicted | k = 2e5: sim/V_n | k = 1e4: sim/V_n | k = 100: sim/V_n | acceptance fraction, k = 2e5 (N_acc/n) |
|---|---|---|---|---|---|
| 1 | 1.000 | 0.999 ± 0.003 | 1.003 ± 0.003 | 1.081 ± 0.003 | 0.498 (0.500) |
| 2 | 0.833 | 1.005 ± 0.003 | 1.009 ± 0.003 | 1.062 ± 0.003 | 0.416 (0.417) |
| 3 | 0.744 | | 1.007 ± 0.003 | | |
| 5 | 0.647 | 1.005 ± 0.004 | 1.003 ± 0.003 | 1.026 ± 0.004 | 0.290 (0.290) |
| 10 | 0.543 | 1.001 ± 0.004 | 0.998 ± 0.003 | 1.004 ± 0.004 | 0.202 (0.202) |
| 30 | 0.432 | 0.987 ± 0.006 | 1.004 ± 0.005 | 0.953 ± 0.005 | 0.101 (0.101) |
| 100 | 0.357 | 0.994 ± 0.008 | 0.990 ± 0.007 | 0.918 ± 0.006 | 0.042 (0.042) |

The record-process formula holds to about 1 percent at k ≥ 1e4 (it is leading order in k: the n = 1
row at k = 100 is the +1/2 term). Kurtosis 2.98 to 3.05.

### 5d. n by step table, k = 1e5: Var(beta)/T_1 with T_1 = 198.67 (se 0.001 to 0.003), acceptance in brackets

| n | eps = 0.25 | 0.5 | 1 | 2 | 4 | 8 | small-step R_n |
|---|---|---|---|---|---|---|---|
| 1 | 0.994 [0.49] | 0.974 [0.47] | 0.946 [0.42] | 0.923 [0.30] | 0.918 [0.18] | 0.920 [0.10] | 1.000 |
| 2 | 0.828 [0.41] | 0.826 [0.40] | 0.823 [0.35] | 0.851 [0.27] | 0.883 [0.17] | 0.903 [0.09] | 0.833 |
| 5 | 0.643 [0.29] | 0.647 [0.28] | 0.666 [0.25] | 0.723 [0.20] | 0.799 [0.14] | 0.858 [0.08] | 0.647 |
| 10 | 0.546 [0.20] | 0.548 [0.20] | 0.568 [0.18] | 0.628 [0.15] | 0.710 [0.11] | 0.795 [0.07] | 0.543 |
| 30 | 0.431 [0.10] | 0.433 [0.10] | 0.457 [0.09] | 0.510 [0.08] | 0.577 [0.07] | 0.657 [0.05] | 0.432 |
| 100 | 0.350 [0.04] | 0.355 [0.04] | 0.383 [0.04] | 0.427 [0.04] | 0.475 [0.03] | 0.530 [0.02] | 0.357 |

Here eps = step/sqrt(T_1(k)), the same absolute step in every row. Over the whole table the variance
is between 0.35 and 1.0 times sqrt(pi k/8) + 1/2: the k^(1/4) width law is robust, the coefficient
depends on numstep and on the step.

### 5e. Questions 1 and 4 (step 1.5, profile variance 1)

| k | legacy move alone: Var(beta), E[X]/k | half-covariance redraw, n = 1 | half-covariance redraw, n = 5 | full redraw (control) | acceptance | 2 I_(1/3)(k/2,k/2) |
|---|---|---|---|---|---|---|
| 1 | 0.9999 ± 0.0005, 0.500 | 1.0006 ± 0.0008 | 1.0014 ± 0.0008 | 1.219 ± 0.001 | 0.51 | 0.78 |
| 4 | 0.9991 ± 0.0008, 0.500 | 1.0002 ± 0.0008 | 1.0005 ± 0.0008 | 1.694 ± 0.002 | 0.37 | 0.52 |
| 16 | 1.0004 ± 0.0014, 0.500 | 1.0004 ± 0.0011 | 0.9981 ± 0.0011 | 2.850 ± 0.004 | 0.13 | 0.18 |
| 64 | see below | 0.996 ± 0.004 | 1.003 ± 0.004 | 5.269 ± 0.019 | 4.9e-3 | 6.2e-3 |

(The acceptance column includes the beta move, so it sits below the frozen-beta value in the last
column.) A first run of "legacy alone" at k = 1 gave 0.9980 ± 0.0008 (-2.6 se); a repeat with a new
seed and twice the length gave the value in the table, and k = 2 and 8 gave 0.9999 ± 0.0005 and
1.0009 ± 0.0007. I read the first as a fluctuation.
k = 64, legacy alone, 1500 chains, 3e4 iterations of which 1e4 burn-in: started from X ~ chi2_k,
Var = 1.021 ± 0.013 and E[X]/k = 0.512, not yet at 0.5; started from X ~ chi2_k/2 (the predicted
stationary law), E[X]/k stays at 0.4995 and Var = 0.977 ± 0.016. Invariance holds, convergence is
already poor at k = 64 (low-X states hold for 1e5 or more iterations). That is the practical content
of section 4.

## 6. Code observations (`mcmc.py:260-348`, `processing.py:862-909`)

1. Minor, documentation. `processing.py:447-449` says the legacy acceptance "does not leave any
   posterior invariant", and the `run` docstring (`mcmc.py:187-191`) says it does not leave the joint
   posterior invariant. The second is true. The first is true only of the composite sweep; the move
   alone is exact MH for pi_K, whose beta marginal is the profile posterior (section 1). The
   `_run_resolve_chisq` docstring's "widens by a factor that grows with the number of pixels" is
   right; the factor is (pi k/8)^(1/4) times an n- and step-dependent coefficient of 0.6 to 1.0.
2. State handling is correct for the move as described: beta and every component's `_data` are
   restored on rejection (`mcmc.py:339-342`), out-of-bounds proposals are rejected without a re-solve
   (`mcmc.py:308-312`), the prior enters next to the likelihood (`mcmc.py:324-325`), and the proposal,
   uniform draw and decision are made on root and broadcast. The stated assumption that
   `resolve_amplitudes` touches only `comp._data` is not checked by the code.
3. The premise "X is fresh at the first move" holds only if the amplitude groups run by the sweep
   (`processing.py:862-865`) are the same as `update_amplitude_groups`, and if a previous MCMC group
   in the same iteration did not leave retained (record-biased) amplitudes behind
   (`processing.py:871-909` loops over MCMC groups without an amplitude redraw between them). With
   two legacy groups the second starts from the first one's last accepted X, so its effective n is
   larger than its `numstep`.
4. The laws above assume the acceptance chi2 is the chi2 the amplitude draw was conditioned on.
   `global_loglike` sums only `chisq_active` bands and only `chisq_mask` pixels, at band resolution
   after a harmonic round trip of the stored amplitudes. Any mismatch with the solve changes the
   distribution of X (no longer chi2_k, possibly beta dependent) and therefore T through the general
   formula of section 2. Not derived further here; see section 7 for the open item.

## 7. Comparison with E035 (added 2026-11-09 23:36 CET, after sections 1 to 6 were on disk)

Read at this point: `nextgen/experiments/E035_legacy_move_tempering_law.md`, `nextgen/experiments/E035/*.txt`,
`nextgen/scripts/e035_legacy_tempering.py`, the F1 and F7 entries of `nextgen/FINDINGS.md`,
`nextgen/experiments/E014_q1_real_maps_lowres.md`, and (read-only) the E014 legacy run log and inputs.

### 7a. Agreement

- Question 1: identical. Same extended target, profile marginal, N(a_hat, (2A)^-1) conditional,
  X ~ Gamma(k/2, 1).
- Question 2, leading order: identical. sqrt(2 pi k)/4 = sqrt(pi k/8); width 0.7916 k^(1/4).
- Question 4: same modification, same conclusion (exact, useless at large k).
- E035's P1 toy rows are reproduced by my exact-kernel reference at E035's own step, including the two
  rows E035 could not explain:

  | k | E035 law | E035 toy | my sqrt(T), T = 1/2 + 1/(4 c_k) | my exact kernel at step 0.3 W | (toy - exact)/se |
  |---|---|---|---|---|---|
  | 10 | 1.408 | 1.577 ± 0.004 | 1.591 | 1.578 | -0.2 |
  | 100 | 2.503 | 2.584 ± 0.007 | 2.604 | 2.586 | -0.2 |
  | 1000 | 4.452 | 4.463 ± 0.012 | 4.508 | 4.479 | -1.3 |
  | 1e4 | 7.916 | 7.870 ± 0.020 | 7.948 | 7.898 | -1.4 |
  | 1e5 | 14.077 | 14.027 ± 0.036 | 14.095 | 14.008 | +0.5 |
  | 190825 | 16.545 | 16.435 ± 0.043 | 16.560 | 16.459 | -0.5 |

- E035's P4 coefficients (no derivation there) are the record-process formula of section 3b:
  0.7916 sqrt(R_5) = 0.637 and 0.7916 sqrt(R_20) = 0.541 against E035's measured 0.64 and 0.54;
  acceptance N_acc/n = 0.290 and 0.132 against E035's 0.286 to 0.297 and 0.131 to 0.139.
- My own two-index toy at the recorded E014 proposals (5.6 and 93 ridge sd, 20 steps, k = 190825)
  gives 13.11 ± 0.03, 13.55 ± 0.01, acceptance 0.0669; E035 has 13.09, 13.56, 0.0669.

### 7b. Where I go further than E035

1. Finite k: T = 1/2 + 1/(4 c_k) = sqrt(pi k/8) + 1/2 + O(k^(-1/2)), exact in the small-step limit.
   This is E035's open "O(k^(-1/2)) correction": the 12 percent (k = 10) and 3.2 percent (k = 100)
   excesses of P1 are the +1/2. P1 did not fail at k = 100; the stated 3 percent band was too tight
   for a leading-order law.
2. The 0.993 to 0.996 ratios of P1 at large k are the finite step 0.3 W (exact kernel: 0.994), not
   noise. "Reversible with respect to exp(-chi2_min/(2T))" is exact only in the small-step limit;
   for one move per redraw the variance is 0.92 to 1.00 of T depending on the step (section 5b).
3. n moves per redraw: derived (record process), sd ∝ about (ln n)^(-1/4); E035 has measurements only.
4. The law in terms of the acceptance noise rather than k: T ≈ 1/2 + 0.4431 sd(chi2 at a fresh draw).
   This is what resolves the open item (7d).

### 7c. Errors and weak points in E035

1. **Acceptance rate of the exact chain (P2 text, Result, Conclusion 3): "exp(-k/4)" is wrong.**
   The stationary acceptance with beta frozen is 2 I_(1/3)(k/2, k/2), decaying as
   (8/9)^(k/2) = exp(-0.059 k), 4.2 times slower in the exponent. E035's own P2 numbers contradict
   exp(-k/4): measured 0.282, 0.197, 0.104 at k = 4, 8, 16; exp(-k/4) gives 0.368, 0.135, 0.018;
   2 I_(1/3) gives 0.519, 0.347, 0.176, and the measured values are a steady 0.54 to 0.59 of that
   (the 2.4 sd beta proposal). "As exp(-k/4) predicts" should be removed. Severity low: the conclusion
   (unusable at map-level k) stands, but at k = 100 the acceptance is 6e-4, not 1e-11.
2. **The error bar on the real-map width is misread (Result, Conclusion 4, and F1: "8 percent error
   each").** ESS 74 belongs to beta_s. beta_d has ESS 1184 (`E014_q1_real_maps_lowres.md:178`; I get
   1215 from the log), so 16.5 carries 2 percent, and the log holds a second chain: beta_d 16.40 and
   16.59 ridge sd (ESS 1215 and 1294), beta_s 17.04 and 16.43 (ESS 74 and 54). Against the toy's
   13.55 the beta_d excess is 22 percent at more than 8 sigma. The candidate "the width of a chain
   with ESS 74" is excluded for beta_d.
3. **"The toy reproduces the chain's acceptance (0.067 against 0.073)" is too generous.** The two
   chains have 0.0735 and 0.0765, each from about 2350 accepted moves (2 percent). 0.0669 is 9 to 13
   percent low, and that deficit is a second symptom of the missing acceptance noise (7d).
4. **Candidate 1 is right in outcome and wrong in mechanism.** The part of the alm storage that F7
   describes (the white fluctuation loses variance above lmax) lowers the acceptance noise, to 0.74
   of its exact value. The increase comes from the part F7 calls fine, the mean (7d).
5. Reproducibility: `nextgen/scripts/e035_legacy_tempering.py` has no P5 block and no anisotropic
   proposal, so `E035/toy_at_earlier_settings.txt` and `E035/toy_e014_recorded_proposals.txt` cannot
   be regenerated from the committed script (my toy reproduces the second). The Setup paragraph
   (100 burn-in, 300 kept) does not match P1 (300, 600) or P2 (400, 1200, 20 steps). The record as
   I read it at 23:21 ended with leftover template sections ("Conclusion: Pending"); those were
   removed from the working copy by someone else at 23:35, while this review was being written.
6. Minor: P4's "proposals of 2.4 legacy widths" (taken from E014) is 0.33 (beta_s) and 5.6 (beta_d)
   legacy widths; the recorded-proposal run uses the right values. The script's Monte Carlo error uses
   the scatter of within-chain variances and ignores the scatter of chain means; adequate here, a
   little optimistic for P1. All chains start at beta = 0; the half-split check is one-sided. None
   of these changes a number beyond its quoted error.
7. No algebra error in the derivation of items 1 to 3 of "Target and priors". No bug in the toy loop.

### 7d. The open item: the real-map width (22 to 31 percent above the toy)

Resolved, with numbers, by three independent routes. No Commander4 run was needed; nothing was
written outside the scratchpad.

**(i) The chain's own log gives the acceptance noise.** `_run_resolve_chisq` logs dloglike and both
states for every proposal (`mcmc.py:332`). At step 1 of each iteration the current and the proposed
amplitudes are both fresh draws, so dloglike = -(D + X' - X)/2 with no selection. Regressing the
3200 post-burn step-1 values (two chains) on a polynomial difference P(beta') - P(beta) leaves a
Gaussian residual (kurtosis 3.01) of sd 682 ± 9, the same for polynomial order 2, 3, 4, 5, for each
chain (665, 695) and each half. Exact draws would give sqrt(k) = 437. So

    Var(chi2 at a fresh draw) = (2.44 ± 0.06) x 2k,   k = 190825.

**(ii) The toy with that one measured number reproduces the chain.** Two indices, 20 steps, recorded
proposals, acceptance-noise variance scaled by f:

| f | sd beta_s | sd beta_d | acceptance |
|---|---|---|---|
| 1.00 (E035's toy) | 13.11 | 13.55 | 0.0669 |
| 2.00 | 15.28 | 15.86 | 0.0734 |
| 2.44 (measured in i) | 15.95 | 16.57 | 0.0751 |
| 2.69 (measured in iii) | 16.37 | 16.94 | 0.0761 |
| E014 chain 1 / chain 2 | 17.04 ± 1.4 / 16.43 ± 1.6 | 16.40 ± 0.33 / 16.59 ± 0.33 | 0.0735 / 0.0765 |

Widths and acceptance agree at f = 2.44. That the measured widths "sit at the one-step law (16.5)"
is a coincidence of 0.82 (20 steps at these proposals) times 2.44^(1/4) = 1.25.

**(iii) The alm storage produces that noise, computed directly on the E014 inputs.** Own script, no MH
code: bands through Commander4's reader, per-pixel A_p and a_hat at the exact ridge mean, full-sky
draws a_hat + delta, each component map through `map2alm_healpix(niter=3, lmax=191)` and
`alm_to_map` as `solve_compsep_perpix` and `get_sky_at_nu` do, chi2 inside the 38165-pixel mask as
`local_loglike` does. X = chi2 - chi2_min, 430 draws:

| quantity | mean / k | variance / 2k |
|---|---|---|
| pixel-space draw (exact; expect 1, 1) | 1.000 | 1.05 ± 0.07 |
| after the alm round trip, as `global_loglike` sees it | 2.634 | 2.69 ± 0.18 |
| quadratic part (P delta)^T A (P delta) | 0.732 | 0.74 ± 0.05 |
| cross term -2 (P delta)^T A (Q a_hat) | 0.000 | 1.92 ± 0.12 |
| offset B = (Q a_hat)^T A (Q a_hat), deterministic | 1.903 | 0 |

P is the storage round trip, Q = 1 - P. Stored amplitudes are P(a_hat + delta), so
X = |P delta - Q a_hat|^2_A. 2.69 ± 0.18 agrees with the 2.44 ± 0.06 of (i). The loss of fluctuation
variance alone would have narrowed the legacy posterior by 7 percent. The noise is added by the
cross term between the band-limited fluctuation and the round-trip residual of the mean amplitude
maps, which is not small: B = 1.90 k = 3.6e5 in chi2 (9.5 per kept pixel) on chi2_min = 8.37e6.
B is not an iteration artefact (2.35 k, 1.90 k, 1.83 k, 1.76 k at niter 0, 3, 10, 30) and not a
mask-edge effect only (median 2.4 to 4.7 per pixel at every latitude; the top 10 percent of pixels
hold 60 percent). Single and double precision give the same B, so float32 is not involved. The
per-pixel mean maps at nside 64 are simply not band-limited at lmax = 3 nside - 1: the round-trip
residual is 0.5 to 1.6 marginal sd per component. F7's "with common-resolution smoothing the mean
is fine" does not hold at the chi2 level on these maps.

**(iv) The same offset explains the legacy mean shift.** The legacy acceptance sees
chi2_min(beta) + B(beta). Finite differences (double precision) give dB/dbeta = +0.83 and -6.38 per
ridge sd of beta_s and beta_d, hence a predicted legacy-minus-ridge shift of -0.11 sd (beta_s) and
+3.12 sd (beta_d). E014 measured +3.0 ± 0.5 in beta_d (the second chain in the log: +3.0 as well) and
-2.7 ± 2.0 / -0.9 in beta_s. The "6 sigma" dust shift of F1 is the gradient of B, not a property of
the move. (Check on my loader: the Hessian of chi2_min returns 0.995 of both quoted ridge sds and a
ridge correlation of 0.095, so "correlated indices" is not a candidate either.)

Limits of this: (i) uses dloglike as logged (4 significant digits) and assumes D is polynomial over
the proposal range; (iii) and (iv) re-implement the draw and the round trip with the same library
calls, they do not execute `solve_compsep_perpix` or `global_loglike` themselves; the toy in (ii)
carries only the variance of the noise. E035's planned Commander4 run with and without alm storage
remains the end-to-end confirmation, but its outcome is now predicted: without the round trip the
legacy widths should fall to 13.1 and 13.6 ridge sd at these settings and the beta_d shift should
vanish.

Consequences for the records:
1. E035 Conclusion 4 and F1: the real-map width is explained. State the law through the acceptance
   noise, T_1 ≈ 1/2 + 0.4431 sd(chi2 at a fresh draw), times the numstep and step factor; k enters only
   when the stored amplitudes are the drawn amplitudes. "The ridge mean" holds only then as well.
2. F7 needs amending: inside the E014 mask the stored mean costs 1.9 k in chi2 and its gradient moves
   beta_d by 3 ridge sd for any acceptance evaluated through `global_loglike`. That includes the
   `conditional` mode (`mcmc.py:210`), whose Gibbs pair is then (beta given P a, a given beta followed
   by projection), not a Gibbs sampler of the stated joint. Not derived further here; it needs its
   own check. `ridge` and `marginal` accept on the pixel-space kernel and are not affected.
3. The k^(1/4) extrapolations of the legacy width to nside 512 and above ("45 at k = 1e7") assume
   exact stored draws. With alm storage the factor depends on B, hence on lmax, beam and signal to
   noise, and is not predicted by k.

## 8. Test gaps

- `tests/test_mcmc_likelihood.py:330-336` asserts only width ratio > 2.0 for `resolve_chisq`. A
  quantitative check is now available: 150 pixels x 3 components, 4 steps: sqrt(T(450) R_4) = 3.08
  at small steps (E035 toy 3.04 at that test's step, E005 measured 3.2). Assert a two-sided band
  justified by the chain's ESS, and add the negative control (half-covariance redraw returns 1).
- No test covers the storage round trip inside an acceptance: chi2 at stored amplitudes against chi2
  at the drawn pixel amplitudes (mean k and variance 2k inside the mask) would have caught 7d.

## 9. Verdict

- Legacy move alone: correct MH kernel, but for pi_K (profile marginal in beta, half-covariance
  amplitudes), not for the stated joint posterior. Composite sweep as coded: wrong for any stated
  target; stationary law of beta is the profile posterior (prior included) tempered by
  T = 1/2 + 1/(4 c_k) ≈ sqrt(pi k/8) + 1/2 for one small move per redraw, times R_n for n moves
  (0.83, 0.65, 0.54, 0.47, 0.36 at n = 2, 5, 10, 20, 100), within 8 percent of that in variance for
  any step size at n = 1, and with a step-dependent partial return toward n = 1 at n > 1.
- E035: derivation correct and in agreement at leading order; three errors of reading (exp(-k/4),
  the ESS of beta_d, the acceptance match), one of mechanism (which part of the alm storage adds
  noise). Its open item is closed by 7d.

## Appendix: core of the toy (the scratchpad copy is not permanent)

```python
import numpy as np
from math import lgamma, log, exp, sqrt

def c_k(k):      # E[exp(-Z/2); Z > 0], Z = X' - X, X, X' iid chi2_k
    return exp(-(k + 1) * log(2.0) + lgamma(k + 1.0) - 2.0 * lgamma(k / 2.0 + 1.0))

def T_small(k):  # one small move per redraw: stationary law = profile^(1/T)
    return 0.5 + 1.0 / (4.0 * c_k(k))

def simulate(k, n, sig, nch, nit, burn, rng, redraw=1.0):
    """chi2_min = beta^2. redraw = 1.0: X ~ chi2_k each iteration (the sweep); 0.5: X ~ chi2_k / 2
    (exact for the profile posterior); None: no redraw (legacy move alone)."""
    beta = np.zeros(nch); beta[nch // 2:] = 2 * sqrt(T_small(k)) * rng.standard_normal(nch - nch // 2)
    X = rng.chisquare(k, nch); s2 = np.zeros(nch)
    for it in range(nit):
        if redraw is not None:
            X = redraw * rng.chisquare(k, nch)
        for _ in range(n):
            bp = beta + sig * rng.standard_normal(nch); Xp = rng.chisquare(k, nch)
            acc = np.log(rng.random(nch)) < -0.5 * (bp * bp - beta * beta + Xp - X)
            beta = np.where(acc, bp, beta); X = np.where(acc, Xp, X)
        if it >= burn:
            s2 += beta * beta
    v = s2 / (nit - burn)
    return v.mean(), v.std(ddof=1) / sqrt(nch)   # Var(beta) and its standard error

# n moves per redraw, small steps: Var = sqrt(k) N_acc / (2 sqrt2 N_phi),
# N_acc = sum_{j=1..n} 1/(j+1), N_phi = sum_{j=1..n} j * integral phi(t)^2 (1 - Phi(t))^(j-1) dt.
```

Real-chain checks of 7d (scratchpad, read-only on the data): `parse_e014_log.py`,
`analyze_e014_log.py` (log regression), `run_F.py` (two-index toy with scaled noise),
`e014_load.py`, `e014_roundtrip.py`, `e014_pool.py`, `e014_Bgrad.py`, `e014_Bmap.py` (direct round-trip
calculation). Inputs: `nextgen/params/q1/param_q1_n64_resolve_chisq_npipemask.yml` and
`/Users/trygvels/data/c4_chains/q1_n64/resolve_chisq_npipemask/logs/run-20261025T200835.395526+0100-op8e.log`.
The recipe for (i): keep step 1 of every post-burn iteration, least-squares fit of dloglike on the
monomials of degree 1 to 3 of the proposed state minus those of the current state, residual variance
divided by k. For (iii): 430 draws of X as defined there.
