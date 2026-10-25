# E012: Hierarchical index-field prototype (M2.1): exact target, then exact-target checks and SBC

- **Date started / finished:** 2026-10-25 /
- **Branch and commit:** nextgen, 0c0c612 (target v1); v2 after the independent review, 2026-10-25
- **Open question addressed:** Q3 (optimal spatially varying index model), Q4 (sampler), Q2 (prior)
- **Validation level:** L1 (exact target) and L2 (simulation-based calibration)

Version 2. Two stats-reviewer agents reviewed v1 (mathematics; validation design). Their reports are
in `E009_c3_nonlin_review/independent_reviews_2026-10-25.md` (items e012-math and
e012-validation-design). Every change from v1 is marked **[v2]** with the reason.

## Question
Can a spectral-index field with an inferred smoothness prior (mean, amplitude and correlation length
sampled as hyperparameters) be sampled exactly and efficiently, with amplitudes collapsed, on
per-pixel simulations? If yes, this is the core of Paper A and M2.4 implements it in Commander4. If
the funnel (Grumitt et al. 2020) or the hyperparameter mixing cannot be tamed, Q3 falls back to
fixed-hyperparameter fields or coarser pixel ladders.

## Target and priors

### Data model (per-pixel, model-consistent; the simulations are drawn from exactly this)
Bands nu_b, b = 1..B, intensity only in the first prototype. Data pixels p = 1..N at nside_d. For
each p:

    d_p = M_p(beta_{r(p)}) a_p + n_p,    n_p ~ N(0, N_p),   N_p = diag(sigma_{b,p}^2)

with a_p in R^C the component amplitudes (C = 3: CMB, synchrotron, thermal dust), M_p the B x C
mixing matrix in uK_RJ (CMB column fixed; synchrotron column (nu/nu_0s)^beta; dust column a
modified black body with fixed beta_d and T_d in the first prototype), and r(p) the index pixel
that contains data pixel p (nested ordering, nside_beta <= nside_d). beta in R^{N_beta} is the
synchrotron index field. No beams, no smoothing, white noise: the per-pixel likelihood is exact for
this model (F7 and F8 do not enter by construction).

### Amplitude treatment: collapsed, three variants
With A_p = M_p^T N_p^-1 M_p (C x C), b_p = M_p^T N_p^-1 d_p:

- **Proper Gaussian amplitude prior (sampler-correctness target) [v2]:** a_p ~ N(m_p, T_p) with
  known m_p, T_p (diagonal). **[v2.1, 2026-10-25]** Default width T_p^(1/2) = max(0.5 |m_p|, 20)
  uK_RJ for synchrotron and dust and 70 uK_RJ for the CMB (m = 0); the first-build setting
  10 max(|m_p|, 10) uK_RJ stays available as "wide" (it allows negative and 10x amplitudes and
  makes |log pi| ~ 1e8, which spoils finite-difference checks). Exact marginal

      log L_G(beta) = sum_p [ 1/2 b'_p^T (A_p + T_p^-1)^-1 b'_p - 1/2 log det (A_p + T_p^-1) ] + const,
      b'_p = b_p + T_p^-1 m_p,

  constant terms independent of beta dropped. This is the ONLY target for which the SBC is a test
  of the sampler, because the simulations can be drawn from exactly this prior (review finding:
  neither the flat prior nor the profile corresponds to a proper generative amplitude prior).
- **w-family (frequentist coverage study, principle 4):** log L_w(beta) = sum_p [ -1/2 chi2min_p
  - (w/2) log det A_p ], chi2min_p = d_p^T N_p^-1 d_p - b_p^T A_p^-1 b_p. w = 1 is the flat-prior
  marginal; w = 0 is the profile ("ridge") likelihood, **[v2]** which equals the marginal under
  the Jeffreys prior of the amplitude block CONDITIONAL on beta, p(a_p | beta) proportional to
  det A_p(beta)^(1/2), with a flat prior on beta. (v1 said "joint Jeffreys prior on (a_p, beta)";
  the review showed by quadrature that the full joint Jeffreys prior sqrt(det I(a, beta)) carries
  an a-dependent factor and gives a different beta posterior, dlogL +0.052 vs +0.006 in its toy.
  FINDINGS F4 is corrected accordingly.) Default for science runs stays w = 0 (D003). The w-family
  runs are reported as bias and coverage versus S/N with template amplitudes, never as a pass/fail
  test of the sampler.

### Index-field prior (hierarchical Gaussian field, non-centred)
Hyperparameters theta = (mu, log sigma, u), with u the unconstrained rho coordinate below. Field:

    beta = mu 1 + sigma S(rho) z,   z ~ N(0, I_K),

with K the number of basis modes (see below), so beta | theta ~ N(mu 1, sigma^2 S S^T).

Two covariance families, both through a fixed orthonormal basis V (N_beta x K) and a spectrum
s_i(rho) >= 0:

    S(rho) = V diag( s_i(rho) ) / c(rho),   c(rho)^2 = (1/N_beta) sum_i s_i(rho)^2,

so that the pixel-averaged marginal variance of beta - mu 1 is exactly sigma^2.

**[v2] The constant mode is removed from both families** (s_0 = 0, i.e. the lambda_0 = 0
eigenvector of the graph Laplacian and the l = 0 harmonic are dropped and c is normalized over the
remaining modes). Reason: with it in, 85 to 94 percent of the unit variance sat in the constant
mode at rho = 180 deg, exactly degenerate with mu; mu now carries the field mean and sigma the sd
of fluctuations about it, which is the quantity Q3 asks for.

**[v2] One unit system for both families.** rho and kappa are in radians; the practical range
of a Matern field with alpha = 2 (nu = 1 in the SPDE correspondence of Lindgren et al. 2011) is
rho with kappa = sqrt(8 nu) / rho = sqrt(8) / rho.

- **GMRF / SPDE (leading candidate):** V from the eigendecomposition of the HEALPix 8-neighbour
  graph Laplacian G at nside_beta (G = D - W). **[v2]** The eigenvalues are rescaled to the
  continuum Laplace-Beltrami operator, lambda_i -> lambda_i / (3 h^2) with h the mean pixel
  spacing in radians (the review measured lambda / (h^2 l(l+1)) = 2.97 at l = 1 and 2.65 at l = 8
  at nside 8 for the unweighted graph Laplacian; the residual 10 percent departure at l ~ nside is
  recorded as a property of the GMRF family, not corrected further). Spectrum s_i = (kappa^2 +
  lambda_i / (3 h^2))^(-1), alpha = 2. Precision Q = c^2 sigma^-2 V diag(s_i^-2) V^T.
- **alm Gaussian random field (Commander3-style baseline, hyperparameters inferred):** V = real
  spherical harmonics on the index pixels for 1 <= l <= lmax = 2 nside_beta (**[v2]** reduced from
  3 nside - 1 so that the pixel-evaluated harmonics are close to orthogonal; the matrix is
  orthonormalized by QR and K = (lmax + 1)^2 - 1 < N_beta is stated explicitly), s_l = (kappa^2 +
  l(l+1))^(-1). The two families differ in basis, band limit and the discrete-Laplacian spectrum.

Hyperpriors (all proper):
- mu ~ N(mu_0, s_mu^2), mu_0 = -3.0, s_mu = 0.5 (synchrotron; dust would use 1.55, 0.3).
- sigma ~ half-Normal(0, 0.3) on the sd, sampled as log sigma with the Jacobian + log sigma.
- rho: log-Uniform on [rho_min, rho_max], rho_min = h (one pixel spacing), rho_max = pi
  (180 deg). **[v2]** Sampled as the unconstrained u = logit((log rho - log rho_min) / (log rho_max
  - log rho_min)) with Jacobian + log[(log rho_max - log rho_min) s(u)(1 - s(u))], so the target is
  smooth and NUTS never meets a hard wall. The second option is the joint PC prior of Fuglstad et
  al. 2019 on (rho, sigma) with P(rho < rho_0) = 0.05 at rho_0 = 2 h and P(sigma > 0.5) = 0.05,
  on the same u coordinate (rho then has a soft prior inside the same bounds).
- No hard bounds on beta (the Gaussian field prior replaces them). The C++ kernel skips pixels
  whose A_p Cholesky fails; **[v2]** the prototype asserts nskipped == 0 on every evaluation, so a
  skip is an error, not a silent change of target.

### Full posterior sampled
Variables x = (z in R^K, mu, log sigma, u), dimension K + 3:

    log pi(x | d) = log L(beta(x)) + log N(z; 0, I) + log N(mu; mu_0, s_mu^2)
                    + log halfN(sigma; 0.3) + log sigma
                    + log p_rho(rho(u)) + log |d rho / d u|          (log-uniform: constant + logit Jacobian)
                    + const,

with log L one of log L_G (sampler test), log L_0, log L_1. Reported quantities: beta = mu 1 +
sigma S z, theta, and the field mean and spatial sd.

### Gradient
d log L / d beta_r = sum_{p in r} g_p, with for the w-family g_p = x_c g^T N_p^-1 (d_p - M_p x) -
w (A_p^-1 M_p^T N_p^-1 g)_c, g = dM_p[:, c]/dbeta, x = A_p^-1 b_p (the envelope theorem for the
chi2min term; confirmed against finite differences by the review to 2e-7). For log L_G replace
A_p by A_p + T_p^-1, b_p by b'_p and x by (A_p + T_p^-1)^-1 b'_p, with w = 1. The C++ kernel
`marginal_loglike_regions` returns the w-family per-region gradient; the prototype uses a NumPy
implementation and checks the kernel against it. Chain rule to x: d/dz = sigma S^T g, d/dmu =
1^T g, d/dlog sigma = sigma (S z)^T g + 1 - sigma^2/0.09, d/du = sigma z^T (dS/du)^T g + prior and
Jacobian terms, with d(s_i/c)/dlog rho = (ds_i - s_i dc/c)/c, ds_i/dlog rho = alpha kappa^2
s_i^2 (for alpha = 2: 2 kappa^2 s_i^2), dc/dlog rho = (1/(c N_beta)) sum s_i ds_i, and dlog rho/du
= (log rho_max - log rho_min) s(u)(1 - s(u)). (The review checked this chain rule, including the
c(rho) term, to 6e-6 on a ring Laplacian.)

### Sampler
NUTS (`src/commander4/math_utils/nuts.py`) over x jointly, dual averaging and diagonal mass
adaptation in warm-up only, samples kept after warm-up. Divergence count and E-BFMI are recorded
per run (the funnel diagnostic). **[v2]** The centred parameterization (beta, mu, log sigma, u)
with its log-det and -K log sigma terms is implemented as well, as an efficiency comparison, not a
negative control: at strong data it is the better geometry. Second sampler for Q4 later (MAMS).

## Setup
- nside_d = 16 (3072 data pixels), nside_beta in {2, 4, 8}.
- Bands: 23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353 GHz, white noise per pixel at levels
  scaled from E006 (RMS x1 and x5).
- Amplitudes for simulations: for the sampler-correctness SBC, a_p ~ N(m_p, T_p) with the template
  means m_p (CMB 0; synchrotron and dust the plane-brightening templates of
  `aux/validate_mcmc_likelihood.py`) and T_p as above. For the coverage study, a_p = template
  values times log-normal scatter (not drawn from any prior the sampler uses; that is the point).
- Code: `nextgen/prototypes/indexfield/` (NumPy + the existing NUTS; numba where needed). No
  Commander4 changes in M2.1.
- Machine: Mac mini (T0/T1). **[v2]** Compute plan: the full 300-replicate SBC runs for the default
  configuration only (log L_G, GMRF family, nside_beta 4, log-uniform rho, non-centred); the other
  factors (alm family, nside_beta 2 and 8, PC prior, centred) run as 100-replicate smoke tests with
  the reduced power stated. The review estimated the full factorial at 7200 NUTS runs, several
  days per configuration, which does not fit T1.

## Reference
1. **L1 exact-target checks** (every component of the code, not only the whole):
   a. Gradient of log pi against central finite differences at random x, tolerance
      |g_fd - g| < 1e-6 |g| + 1e-4 (**[v2]** absolute floor from roundoff on a log density of
      order 1e4 to 1e5; this check protects efficiency, 1b to 1d protect the target).
   b. **[v2] Exact grids with the real SED:** N_beta = 1 (1-D grid in beta, theta fixed) and
      N_beta = 2 (2-D grid), for log L_G, log L_0 and log L_1; NUTS means within 4 MC standard
      errors (ESS-based) and sd ratios within 0.95 to 1.05 of the grid. No MC error in the
      reference. This is the check that catches a sign or factor error in the chi2min gradient
      or the log-det term.
   c. **[v2] Hyperparameters, independent centred Gibbs reference:** nside_beta = 1 (12 index
      pixels), nside_d = 16. Reference chain written independently of the prototype's log density,
      in centred coordinates: beta | theta, d by random-walk MH (12-D, near Gaussian); mu | sigma,
      rho, beta exactly Gaussian (conjugate); (log sigma, log rho) | beta, mu by an exact draw from
      a fine 2-D grid (the likelihood does not enter this conditional). The review showed this has
      about 500x the ESS of a non-centred random walk in this geometry (beta pinned to 0.02 by 256
      data pixels per index pixel while sigma is 0.3). Report ESS of both chains; compare with the
      combined SE sqrt(SE_ref^2 + SE_nuts^2); require ESS > 800 on log sigma in both chains before
      the comparison counts.
   d. **[v2] Linear-Gaussian control with KNOWN fixed amplitudes:** d_p = m_p + F_p beta_r(p) +
      n_p. Then beta | theta, d ~ N(Sigma (F^T N^-1 (d - m) + Q mu 1), Sigma), Sigma = (F^T N^-1 F
      + Q)^-1, and theta | d proportional to p(theta) N(d - m; F mu 1, sigma^2 F K F^T + N) on a
      3-D grid. Zero MC error in the reference; tests the basis, c(rho), the Jacobians and the
      hyperpriors (it does not exercise the SED or the collapse, which 1b does). v1 claimed a
      closed form for the collapsed likelihood with a linear SED; the review showed that chi2min
      and log det A are rational, not quadratic, in beta even then.
   e. **[v2] Negative controls (wrong targets, must fail 1c or 1d at the stated ESS):** (i) the
      centred parameterization with the -K log sigma - 1/2 log det K(rho) term omitted (the classic
      error); (ii) the log rho Jacobian dropped (a 1/rho tilt over the whole prior range, a factor
      3 at nside_beta 1 and 25 at 8, visible in 1d). The dropped log sigma Jacobian shifts the log
      sigma posterior by only 0.2 sd at K = 12 and would need ESS > 800 in both chains; it is kept
      as a third control with that precondition stated. The reference log density is written
      separately from the prototype's so that a shared error cannot cancel.
2. **L2 SBC** (log L_G only, matching generator): 300 replicates for the default configuration.
   Draw theta from the hyperprior, z ~ N(0, I), a_p ~ N(m_p, T_p), data from the likelihood. Run
   NUTS with a fixed budget; **[v2]** thin by per-replicate ESS to L = 99 draws (100 ranks, 20
   equal bins) with the precondition ESS > 400 on every test quantity; replicates that miss it or
   diverge are rerun with a smaller step size or a longer budget, never dropped. Test quantities:
   mu, log sigma, u, four beta pixels (two bright, two faint), the field mean of beta, z^T z,
   (beta - mu 1)^T K^-1 (beta - mu 1)/sigma^2, the log prior density of beta given theta, the
   joint log likelihood (Modrak et al. 2022), and TARP coverage of the full beta field (Lemos et
   al. 2023) as VALIDATION.md requires. Pass: chi-square on 20 bins with a Bonferroni-corrected
   level over the test quantities of one configuration, plus the uniformity of the p-value
   distribution itself, and the ECDF simultaneous 95 percent band (Sailynoja et al. 2022).
   **[v2] Stated power** (n = 300, 20 bins, alpha 0.01): a posterior sd off by a factor 0.8 is
   detected with probability 0.84, by 0.9 with 0.09, by 1.2 with 0.36; a bias of 0.25 sd with
   0.53, 0.5 sd with 1.0. Ten to twenty percent width errors are below the resolution of this
   SBC and are caught by 1b to 1d instead.
3. **Coverage study** (w = 0 and w = 1, template amplitudes): coverage and bias of beta pixels and
   hyperparameters versus S/N, reported as curves, no pass/fail (principle 4). This is where the
   volume term is quantified for the field model.

## Command
```bash
cd nextgen/prototypes/indexfield
../../../.venv/bin/python -m indexfield.model      # self-tests (run from nextgen/prototypes)
../../../.venv/bin/python checks/check_1a_gradient.py
../../../.venv/bin/python checks/check_1b_grid.py
../../../.venv/bin/python checks/check_1c_hyper.py
../../../.venv/bin/python checks/check_1d_linear.py [--scale 12 --tag _long]
```

## Result
**Build (2026-10-25, commit 49f913a).** Modules written by three agents in parallel against
`SPEC.md`; `reference.py` by an agent that was not allowed to read the other modules (independent
code path, own plain-loop likelihood, own basis construction, own 2-D grid draw). Integration and
checks by a fourth agent; results in `nextgen/prototypes/indexfield/checks/out/`. Numbers below are
from that report; the three runs it had to cut short are being rerun.

Module self-tests: gradient of the collapsed likelihood vs finite differences, relative error below
1e-7 for gauss, w0, w1; mode gauss equals brute-force 2-D quadrature over a_p to 1e-13; one
evaluation at N = 3072, B = 12 takes 3 to 4 ms. Bases orthonormal to 1e-14; Matern derivatives
including the c(rho) term within 1e-6; centred minus non-centred log density equals
-K log sigma - 1/2 log det K to 4e-13. NUTS on a 20-D correlated Gaussian: 3 seeds x 30000 samples,
sd ratios 0.992 to 1.024. Reference: grid_1d gauss equals brute-force 3-D quadrature absolutely
(2e-11, dropped constants verified); 2-D grid conditional draw chi-square p = 0.38; centred Gibbs
reference reproduces the linear-Gaussian closed form (ESS 6600 to 20000, all |z| < 1).

L1 checks:
- **1a PASS.** Per-index-pixel gradients, all modes, |g_fd - g| max 1.7e-8 against the tolerance
  1e-6 |g| + 1e-4; posterior gradients (both parameterizations, both families, three modes) worst
  err/tol 1.9e-4. The C++ kernel `marginal_loglike_regions` agrees with the NumPy w0/w1
  likelihood to 7e-12 and the gradients to 5e-14, nskipped 0.
- **1b PASS (gauss, w0), w1 pending rerun.** Prototype vs reference per-pixel log density pointwise:
  relative spread 3e-15 to 8e-15 (|log L| ~ 1e6). Grid vs NUTS at N_beta 1 and 2 with a fixed
  Gaussian prior: mean |z| <= 1.76, sd ratios 0.97 to 1.022, ESS 970 to 2900, no divergences. For
  w1 at noise x5 the flat-prior marginal target sits at beta = -0.25 +- 0.002, 3750 log units above
  beta = -3 (F3 in the toy: the -1/2 sum log det A term rewards near-degenerate cmb/sync columns);
  the first run's grid did not cover it (a check defect, not a sampler one); NUTS agrees with the
  wide-grid maximum; the fine comparison is in the rerun.
- **1c PROVISIONAL.** Pointwise: centred posterior equals the reference joint density up to a
  constant (spread 6e-15), Jacobian to 6e-11. NUTS vs centred Gibbs reference at nside_beta 1,
  18 quantities: all |z| <= 2.2, sd ratios 0.95 to 1.09. Precondition ESS(log sigma) > 800 met by
  the Gibbs reference (15500) but not by NUTS (centred 659, non-centred 64 within the budget);
  rerun at 2.5x budget in progress. **Result in itself:** at this data strength the non-centred
  parameterization mixes 10x worse than the centred one (ESS/n 2 percent, E-BFMI 0.78, no
  divergences), as the review predicted: centred is the better geometry at strong data.
- **1d PASS (centred), non-centred provisional.** Closed form vs NUTS for the fixed-amplitude linear
  model, three cases (graph Laplacian with independent bases; literal harmonic basis with
  independent bases; constant-projected harmonic basis on the shared basis): pointwise spread 2e-15,
  centred chains |z| <= 2.9 on means and sds with ESS 2300 to 4100 (long run: 45500, |z| <= 3.4,
  sd ratios 0.992 to 1.011); non-centred chains agree (|z| <= 2.1) at ESS 40 to 70, long run
  pending.
- **1c (rerun, 2.5x budget) PASS.** Centred: ESS(log sigma) 840 vs 2661 (reference), max |z|
  1.67; non-centred: ESS 928, max |z| 1.23, sd ratios 0.976 to 1.063. Precondition met in both.
- **1d (long run, 96000 draws).** Case A (graph Laplacian, independent bases): centred ESS 45500,
  |z| <= 3.4, sd ratios 0.992 to 1.011, PASS; non-centred ESS(log sigma) 5300, |z| <= 2.7, PASS.
  Case C (constant-projected harmonic basis, shared basis): PASS both parameterizations (ESS 33800
  and 2100). Case B (literal harmonic basis): the first long run showed a mu offset of 0.0176 (z 764) in
  both NUTS chains. Explanation after the rerun: the reference's harmonic basis had been switched
  to the constant-projected construction while that run was in flight, so a literal prototype
  basis was compared with a projected reference; with both sides literal (rerun _long2) case B
  passes (mu z -0.58 and +1.67). The offset is the constant-mode leakage itself (0.17 of a unit
  constant at nside 1), which is why the SPEC removes it. **Decisive rerun _long2 (96000 draws per
  chain):** case A PASS (centred ESS 45500, non-centred 1490, |z| <= 3.5); case C with INDEPENDENT
  constant-projected bases PASS (centred ESS 44800, non-centred 3170, |z| <= 1.8, sd ratios 0.989
  to 1.009). L1 is complete: 1a to 1e all pass.
- **1b (rerun, 12 cases).** gauss and w0 PASS at noise x1 and x5, N_beta 1 and 2 (|z| <= 3.2,
  sd ratios 0.959 to 1.046, ESS 790 to 2600). w1 PASS in 5 of 6 cases; at N_beta 2, noise x5 the
  flat-prior marginal is MULTIMODAL: the grid mode sits at beta = (-0.256, -0.233) with sd 0.003,
  one NUTS chain started near -3 stayed in a secondary basin at -2.3 (z -1240, sd ratio 14) while
  the other pixel agreed. Rewritten check (four dispersed starts, basin search, local grid
  refinement, mode-mass comparison): all 12 cases PASS with pooled ESS 3450 to 10700, worst |z|
  2.13, sd ratios 0.974 to 1.032; secondary basins of w0/w1 exist (log mass -1660 to -6860 below
  the main one) but carry no measurable mass. The w1 collapse to beta ~ -0.23 at noise x5 stands.
- **1e (negative controls) PASS, i.e. every wrong target failed as expected, positive control
  passed.** (i) centred without -K log sigma - 1/2 log det K: log sigma shifted +3.06 sd (z 182),
  sd ratio 1.21. (ii) logit Jacobian dropped: the wrong target is IMPROPER in u (the likelihood
  saturates at rho_min and rho_max so log pi is flat in both tails), stronger than v2's "1/rho tilt":
  NUTS drifts to |u| ~ 1e8 with rho pinned at a bound. (iii) log sigma Jacobian dropped: shift
  -0.197 sd (z -11.6 centred at ESS 3360; z -7.1 non-centred at ESS 1200 with a 10x budget), as
  v2 predicted (0.2 sd at K = 12).

**Finding for Q2 (field-model version of F3):** with amplitudes integrated under a FLAT prior
(w = 1) the field posterior collapses to beta ~ -0.1 to -0.25 for every index pixel even at
per-index-pixel S/N of 30 (coverage smoke test, nside_beta 2: z(mu) ~ +1400), because the
-1/2 sum_p log det A_p term rewards the near-degenerate CMB/synchrotron columns and is summed over
all pixels. The flat-prior marginal is unusable for spatially varying indices; ridge (w = 0) and
the proper Gaussian prior behave. This confirms D003 and is a Paper A result; it goes into FINDINGS
as F12 once the coverage study quantifies it.

**Independent review of the implementation (stats-reviewer, 2026-10-25): correct with caveats.**
Its own 60-line implementation of the E012 density matches both posteriors up to the declared
constants (1e-12) for all 12 combinations of prior, mode and parameterization; the SBC generator
matches the gauss target term by term; rank binning exact. Caveats being fixed: unclean SBC
replicates must not be excluded (escalating retries instead); posterior must catch the ValueError
from a singular A; check 1b needs dispersed starts; TARP local references from a separate chain.

Properties recorded for the GMRF family: after the 1/(3 h^2) rescaling the graph-Laplacian
eigenvalues over l(l+1) are 0.94 for the l = 1 multiplet at nside 8 (within-multiplet spread 0.77
to 1.10 at l = 2, 3; 0.74 to 0.93 at l = 8) and 0.64 at nside 1, 0.90 at nside 2: "rho in radians"
is approximate at the few-percent level for nside_beta >= 4 and at the 30 percent level at
nside_beta 1. E012 v2's "2.97 at l = 1" (0.99 after rescaling) is softened to "about 1 at low l for
nside >= 4".

Decisions taken after the build (SPEC.md "Decisions"): constant projected out of the harmonic
basis in both code paths; PC prior semantics fixed (exponential on sigma replaces the half-normal);
default Gaussian amplitude prior narrowed to T^(1/2) = max(0.5 |m|, 20) uK_RJ (CMB 70 uK_RJ);
zero-density evaluations must be zero in the sampling phase.

## Negative control
Listed under Reference 1e. The volume-term control for science is w = 1 versus w = 0 at low S/N in
the coverage study.

## Conclusion
(pending)

## Follow-ups
- M2.2 uses this prototype with a scalar index (N_beta = 1) for the full-sky coverage table (Q2).
- Two fields (beta_s and beta_d) jointly, and T_d, after the single-field version passes.
- Heavy-tailed increments (Student-t on z, or horseshoe on neighbour differences) as a variant.
