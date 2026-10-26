# E016: Hierarchical amplitude prior in the collapsed index-field likelihood

- **Date started / finished:** 2026-10-26 / (running: target written, no code yet)
- **Branch and commit:** `nextgen`, 1c25d62
- **Open question addressed:** Q2 (prior for calibrated spectral-parameter posteriors), feeding Q3
- **Validation level:** L2 (SBC of the sampler for the stated target), then the Reference-3 style
  coverage study of E012 on template data (not a pass/fail test, principle 4)

## Question
Does a hierarchical amplitude prior, whose mean scale and width per component are learned from the
data, give calibrated spectral-index posteriors on template data where the fixed Gaussian prior is
biased (E012 study_v1: z +2 to +14) and the flat-prior marginal collapses (z +100 to +300)?

If yes, the Phase 2 index-field sampler keeps the collapsed likelihood with a learned amplitude
prior and the ridge default (D003) is replaced for spatially varying indices. If no, the index
resolution must be set from the S/N map so that the profile likelihood stays calibrated (index-pixel
S/N above about 15, F12), and low-S/N regions are regularised by the field prior alone.

## Target and priors
Per data pixel p (N pixels at nside_d), bands b (B), components c in (cmb, sync, dust):
d_p = M_p(beta_{r(p)}) a_p + n_p, n_p ~ N(0, N_p), N_p diagonal, with the E012 mixing matrix
(sync power law with index beta per index pixel r, dust fixed at beta_d = 1.6, T_d = 20 K).

Amplitude prior (new, stage 1, per-pixel independent so that the collapse stays per pixel):
    a_p,c ~ N( mu_c g_p,c , (s_c g_p,c)^2 )        for c = sync, dust,
    a_p,cmb ~ N( 0, 70^2 )                           (fixed, as in E012),
where g_c is a FIXED positive shape template per component (what a real analysis takes from a
previous component-separation pass or a smoothed single-band map) and (mu_c, s_c) are global
hyperparameters, one scale and one relative width per component:
    mu_c ~ N(1, 0.5^2),      log s_c ~ N(log 0.3, 1^2).
So m_p,c = mu_c g_p,c and T_p,c = (s_c g_p,c)^2 in the notation of `model.collapsed_loglike`
(mode "gauss"), and the collapsed likelihood becomes a function of (beta, mu, s):
    log L(beta, mu, s) = sum_p [ 1/2 b'_p^T A'_p^-1 b'_p - 1/2 log det A'_p ]
                          - 1/2 sum_p m_p^T T_p^-1 m_p - 1/2 sum_p log det T_p + const,
with A'_p = M_p^T N_p^-1 M_p + T_p^-1, b'_p = M_p^T N_p^-1 d_p + T_p^-1 m_p. The two trailing
terms, dropped as constants in E012 because (m, T) were fixed, are now part of the target: the
-1/2 sum_p log det T_p = -sum_p sum_c log(s_c g_p,c) term is what prevents the drift toward the
flat-prior (w1) collapse as s_c grows, because for s -> infinity the full expression is the proper
marginal N(d_p; M_p m_p, M_p T_p M_p^T + N_p), whose own log-determinant penalises amplitude
variance that the data do not show. (The w1 collapse of F12 is the s -> infinity limit of the first
line alone, without this normalisation.)

Index field and hyperpriors exactly as E012 / D009: beta = mu_beta 1 + sigma S(rho) z, centred
parameterisation, mu_beta ~ N(-3, 0.5^2), sigma ~ half-N(0.3), rho log-uniform on [rho_min,
rho_max] through the logit u with its Jacobian, graph-Laplacian basis at nside_beta.

Joint target sampled:
    p(beta, theta_beta, mu, log s | d) proportional to
        L(beta, mu, s) p(beta | theta_beta) p(theta_beta) p(mu) p(log s),
with NUTS on all of (beta or z, mu_beta, log sigma, u, mu_sync, mu_dust, log s_sync, log s_dust);
gradients with respect to (mu, log s) through A' and b' (analytic, to be derived in SPEC.md and
checked by finite differences as in check_1a before anything is run).

Stage 2 (only if stage 1 is calibrated on the well-specified test and still biased on template
data): replace the per-pixel prior by a coupled smoothness prior on each amplitude map,
a_c ~ N(mu_c g_c, S_c(s_c, rho_c)) with a GMRF S_c, where the collapse is global,
-1/2 log det(M^T N^-1 M + S^-1), and needs a sparse Cholesky of a 3N x 3N matrix (9216 at
nside 16). Deferred; its target will be written in a separate record.

## Setup
(to be filled when the code exists) Prototype `nextgen/prototypes/indexfield/`, new mode
"hier" in `model.py` taking (g, mu, log s) with gradients; `posterior.py` extended with the four
hyperparameters; `coverage.py --modes hier`. nside_d 16, nside_beta 2 (48 index pixels, as
study_v1) and 4 (SBC, as the L2 run). Template data as study_v1: a_p,c = m_p,c 10^(0.1 N(0,1)),
which is log-normal and therefore mildly misspecified for the Gaussian hierarchy (relative sd
0.23, skewed). Shape template g: (i) g = m (the true template, "ideal shape"); (ii) g = m times a
smooth multiplicative distortion 1 + 0.3 cos(theta_p) ("wrong shape"), to measure sensitivity to
the template, which a real analysis never has exactly.

## Reference
1. **Sampler correctness (L2):** SBC (Talts et al. 2018, as E012 L2) with the generative model
   equal to the target: draw theta_beta, beta, (mu_c, s_c) from their priors, a_p,c ~ N(mu_c g_p,c,
   (s_c g_p,c)^2), cmb ~ N(0, 70^2), data; sample; rank statistics of all 16 quantities (12 of E012
   plus mu_sync, mu_dust, log s_sync, log s_dust) plus the joint log density; 300 replicates;
   Bonferroni p, ECDF band, TARP, same pre-registered power statement as E012.
2. **Scientific question (coverage, no pass/fail):** `coverage.py` on template data with modes
   hier (g = m), hier (g distorted), gauss (default T, the known-biased comparator), w0, at noise
   0.3 / 1 / 3, 12 replicates; bias z and coverage versus index-pixel S/N exactly as study_v1 so the
   rows are comparable. The answer is "calibrated" if the hier rows have mean z within +-0.5 and
   coverage within binomial error of 0.68 / 0.95 in every S/N bin where w0 is not, i.e. S/N below
   15, at noise 0.3 and 1.
3. **Limit checks (L1, cheap):** s_c fixed at a large value (1e3) must reproduce the w1 column of
   study_v1 up to the normalisation term; s_c fixed at 0.5 and mu_c = 1 must reproduce mode
   "gauss" with the default T exactly (same numbers as `collapsed_loglike(..., "gauss")`).

## Command
```bash
# none yet
```

## Result
(pending)

## Negative control
Mode "gauss" with the fixed default T on the same template replicates: known biased (study_v1,
z +14 at noise 0.3). Mode "w1": known collapsed. Both run through the same analysis.

## Conclusion
(pending)

## Follow-ups
- Stage 2 (coupled smoothness prior on amplitude maps) if stage 1 is insufficient.
- Scalar-index version (M2.2 coverage table) with the same code.
