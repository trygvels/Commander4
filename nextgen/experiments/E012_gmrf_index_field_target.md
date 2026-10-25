# E012: Hierarchical index-field prototype (M2.1): exact target, then exact-target checks and SBC

- **Date started / finished:** 2026-10-25 /
- **Branch and commit:** nextgen, 0c0c612 (target written before any prototype code)
- **Open question addressed:** Q3 (optimal spatially varying index model), Q4 (sampler), Q2 (prior)
- **Validation level:** L1 (exact target) and L2 (simulation-based calibration)

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
that contains data pixel p (nested ordering, nside_beta <= nside_d, so each data pixel has exactly
one index pixel and each index pixel has (nside_d/nside_beta)^2 data pixels). beta in R^{N_beta} is
the synchrotron index field. No beams, no smoothing, white noise: the per-pixel likelihood is exact
for this model (F7 and F8 do not enter by construction).

### Amplitude treatment: collapsed
The amplitudes are integrated or profiled out per pixel. With A_p = M_p^T N_p^-1 M_p (C x C),
b_p = M_p^T N_p^-1 d_p and chi2min_p = d_p^T N_p^-1 d_p - b_p^T A_p^-1 b_p:

    log L_w(beta) = sum_p [ -1/2 chi2min_p(beta) - (w/2) log det A_p(beta) ] + const

- w = 1: flat (improper uniform) prior on every a_p, exact marginal.
- w = 0: profile ("ridge") likelihood, equal to the marginal under the joint Jeffreys prior on
  (a_p, beta) for fixed p (F4). **Default w = 0 (D003).** Both are run in the SBC (Q2); the SBC
  for w = 1 draws a_p from a proper prior that is wide enough that the flat-prior approximation is
  the stated target (see Setup); any failure of w = 1 there is a finding about the volume term.

### Index-field prior (hierarchical Gaussian field, non-centred)
Hyperparameters theta = (mu, log sigma, log rho). The field is

    beta = mu 1 + sigma * S(rho) z,    z ~ N(0, I_{N_beta}),

so that beta | theta ~ N(mu 1, sigma^2 S S^T) with Cov = sigma^2 K(rho), K = S S^T, K_ii = 1 on
average over pixels (K is normalized so that the pixel-averaged marginal variance is 1 and sigma is
the field standard deviation).

Two covariance families, both through one fixed orthonormal basis V (eigenvectors of a fixed
operator on the index pixels) and a spectrum s_i(rho) >= 0:

    S(rho) = V diag( s_i(rho) ) / c(rho),  c(rho)^2 = (1/N_beta) sum_i s_i(rho)^2

- **GMRF / SPDE (leading candidate):** V from the eigendecomposition of the HEALPix neighbour-graph
  Laplacian G at nside_beta (G = D - W, W_ij = 1 for the 7 or 8 HEALPix neighbours), and
  s_i(rho) = (kappa(rho)^2 + lambda_i)^(-alpha/2) with alpha = 2 (Matern nu = 1 on the sphere
  surface in the SPDE correspondence of Lindgren et al. 2011), where lambda_i are the eigenvalues
  of G and kappa is tied to the correlation length rho by kappa = sqrt(8 nu) / rho in units of the
  pixel spacing. The precision is Q = c^2 sigma^-2 V diag((kappa^2 + lambda_i)^alpha) V^T, a
  polynomial in G, hence sparse for integer alpha. The eigendecomposition is done once per run
  (N_beta <= 12288 at nside_beta 32); hyperparameter changes cost O(N_beta^2).
- **alm Gaussian random field (Commander3-style baseline, hyperparameters inferred):** V = real
  spherical harmonics on the index pixels up to lmax = 3 nside_beta - 1 (not exactly orthonormal
  on pixels; orthonormalized once by QR), s_l(rho) = (kappa^2 + l(l+1))^(-alpha/2), the same
  Matern spectrum. Differences from the GMRF are then only the basis and the band limit.

Hyperpriors (all proper):
- mu ~ N(mu_0, s_mu^2), mu_0 = -3.0, s_mu = 0.5 (synchrotron; dust would use 1.55, 0.3).
- sigma ~ half-Normal(0, 0.3) on the standard deviation (density 2 phi(sigma/0.3)/0.3, sigma > 0),
  written in log sigma with the Jacobian sigma.
- rho ~ log-Uniform on [rho_min, rho_max] with rho_min = 1 pixel spacing and rho_max = 180 deg;
  as the second option the joint PC prior of Fuglstad et al. 2019 on (rho, sigma) with
  P(rho < rho_0) = 0.05 at rho_0 = 2 pixels and P(sigma > 0.5) = 0.05. Both are run; the SBC must
  pass for whichever is used.
- Hard bounds on beta are not imposed in the prototype (the Gaussian field prior replaces them). If
  the SED needs them for numerical reasons, they are logged as a change of target.

### Full posterior sampled
Variables x = (z in R^{N_beta}, mu, log sigma, log rho), dimension N_beta + 3:

    log pi(x | d) = log L_w(beta(x)) + log N(z; 0, I) + log N(mu; mu_0, s_mu^2)
                    + log halfN(sigma; 0.3) + log sigma            (Jacobian of log sigma)
                    + log p(rho) + log rho                          (Jacobian of log rho)
                    + const.

Reported quantities are beta = mu + sigma S(rho) z and theta. The marginal posterior of beta under
this target is the quantity Q3 is about; its posterior mean is a Wiener-filtered index map with the
filter inferred from the data.

### Gradient
d log L / d beta_r = sum_{p in r} [ -1/2 d chi2min_p / d beta - (w/2) tr(A_p^-1 dA_p/dbeta) ] with
d chi2min_p / d beta = -2 (d_p - M_p a*_p)^T N_p^-1 (dM_p/dbeta) a*_p at a*_p = A_p^-1 b_p (the
envelope theorem; the term through a*_p vanishes at the minimum). This is what the existing C++
kernel `marginal_loglike_regions` returns per region; the prototype first uses a NumPy
implementation and checks the C++ kernel against it. Chain rule to x: d/dz = sigma S^T g,
d/dmu = 1^T g, d/dlog sigma = sigma (S z)^T g + prior terms, d/dlog rho = sigma z^T (dS/dlog rho)^T g
+ prior terms, with g = d log L / d beta and dS/dlog rho = V diag(d(s_i/c)/dlog rho) analytic.

### Sampler
NUTS (`src/commander4/math_utils/nuts.py`) over x jointly, dual averaging and diagonal mass
adaptation in warm-up only, samples kept after warm-up. Divergence count and E-BFMI are recorded
per run (the funnel diagnostic). Second sampler for Q4 later (MAMS); not part of M2.1's pass
criterion.

## Setup
- nside_d = 16 (3072 data pixels), nside_beta in {2, 4, 8} (48, 192, 768 index pixels).
- Bands: 23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353 GHz (WMAP + LFI + HFI centre
  frequencies), white noise per pixel at levels scaled from E006 (RMS x1 and x5).
- Amplitudes for simulations: CMB ~ N(0, 70 uK) per pixel, synchrotron and dust log-normal around
  a plane-brightening template as in `aux/validate_mcmc_likelihood.py`. For the w = 1 SBC the
  proper amplitude prior is N(0, (10^4 uK)^2) per component, so the implied beta volume term equals
  the flat-prior one to better than 1e-6 in log density.
- Code: `nextgen/prototypes/indexfield/` (NumPy + the existing NUTS; numba where needed). No
  Commander4 changes in M2.1.
- Machine: Mac mini (T0/T1).

## Reference
1. **L1 exact-target checks** (every component of the code, not only the whole):
   a. Gradient of log pi against central finite differences (relative error < 1e-6) at random x.
   b. nside_beta = 1 (12 index pixels) and fixed theta: NUTS beta posterior against a dense
      random-walk MH chain of 10^6 steps on the same log pi (independent code path), means within
      4 MC standard errors by ESS, sd ratios within 0.9 to 1.1.
   c. N_beta = 12, theta free: joint NUTS against the same long MH chain, including the
      hyperparameters.
   d. Linear-Gaussian control: replace the SED by a linear model (log L exactly Gaussian in beta);
      then beta | theta, d and theta | d have closed forms (conjugate Gaussian and a 3-D
      marginal evaluated on a grid). The sampler must reproduce them. This control catches basis,
      normalization and Jacobian errors with no MC error in the reference.
   e. Negative control: run the same checks with the log sigma Jacobian dropped and with the
      centred parameterization at a funnel-prone setting; both must fail c or d visibly.
2. **L2 SBC**: 300 replicates per configuration. Draw theta from the hyperprior, z ~ N(0, I),
   amplitudes from their simulation prior, data from the likelihood. Run the sampler with fixed
   budget, thin to L = 100 effectively independent draws, record ranks of mu, log sigma, log rho,
   four beta pixels (two bright, two faint), and the joint log likelihood (Modrak et al. 2022).
   Pass: chi-square test of rank histograms with 20 bins, p > 0.01 for every quantity, and the ECDF
   difference within the 95% simultaneous band (Sailynoja et al. 2022). Repeat for w = 0 and w = 1
   and for both covariance families.

## Command
```bash
# filled in when the prototype exists
```

## Result
(pending)

## Negative control
Listed under Reference 1e. Also the legacy `resolve_chisq` structure is not applicable here because
amplitudes are collapsed; the control for the volume term is w = 1 versus w = 0 at low S/N.

## Conclusion
(pending)

## Follow-ups
- M2.2 uses this prototype with a scalar index (N_beta = 1) for the full-sky coverage table (Q2).
- Two fields (beta_s and beta_d) jointly, and T_d, after the single-field version passes.
- Heavy-tailed increments (Student-t on z, or horseshoe on neighbour differences) as a variant.
