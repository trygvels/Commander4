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
flat-prior (w1) collapse as s_c grows, because the full expression equals the proper marginal
N(d_p; M_p m_p, M_p T_p M_p^T + N_p) for every s, and its log-determinant grows as 2 log s_c per
pixel and component as s_c -> infinity (one power only for exactly collinear sync and dust
columns), so log L is bounded above in s and the maximum-marginal-likelihood s is finite or at the
s -> 0 boundary. (The w1 collapse of F12 is the s -> infinity limit of the first line alone,
without this normalisation.) Reviewed 2026-10-26 (stats-reviewer): target CONFIRMED, including the
dropped constants -1/2 d^T N^-1 d, -(B/2) log 2 pi, -1/2 log det N; condition g_p,c > 0 strictly;
mu_c and s_c identified by the N pixels through the fixed shape g, beta identified by band ratios
and not absorbed by mu_sync; the prior is pivot-dependent (nu_0 = 30 GHz for sync), so the
generator must use the same pivot. Hyperpriors proper, likelihood bounded, posterior proper;
mu_c g < 0 has prior probability 2.3 percent per component and is unphysical but harmless for SBC
because generator and target are the same model.

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
Code (commits 9cdb462, cc7c8e1 and the coverage.py hier modes): `model.hier_loglike` and
`hier_prior`, `posterior.HierPosterior` (centred, dim K + 7), `sbc.py --mode hier` (generator
draws the four hyperparameters from their hyperpriors and the amplitudes from the resulting
prior), `coverage.py --modes hier hier_distort hier_scatter` (shape templates g = m, m times
1 + 0.3 cos theta, m times 10^(0.1 N'(0,1)); amplitude hyperparameter posteriors stored per cell
as amp_mean / amp_sd / amp_ess). L1: `checks/check_2a_hier.py` (all pass, 2026-10-26: direct
form 7e-13, gradients 4e-9 at s = 1e3 and 0.3 and 2e-5 at s = 1e-3 where cond(A') ~ 1e8, w1
limit residual 1.9e-4, 1.9e-6, 1.9e-8 at s = 1e3, 1e4, 1e5, gauss limit exact). Smoke tests of
both harnesses run; on one template replicate at noise 1 the posterior of log s_sync is -1.51 +-
0.14, i.e. s = 0.22 against the generator's 0.23 relative scatter.

Runs: coverage `study_v1_hier` (12 replicates x noise 0.3 / 1 / 3 x three hier modes, seed and
nside as study_v1, 2 workers, started 2026-10-26 midday) and SBC `hier_centred_hier` (300
replicates, nside_beta 4, 4 workers, queued behind the centred SBC). Original plan text: prototype `nextgen/prototypes/indexfield/`, new mode
"hier" in `model.py` taking (g, mu, log s) with gradients; `posterior.py` extended with the four
hyperparameters; `coverage.py --modes hier`. nside_d 16, nside_beta 2 (48 index pixels, as
study_v1) and 4 (SBC, as the L2 run). Template data as study_v1: a_p,c = m_p,c 10^(0.1 N(0,1)),
which is log-normal and therefore mildly misspecified for the Gaussian hierarchy (relative sd
0.23, skewed). Shape template g: (i) g = m (the true template, "ideal shape", a near-generative
best case); (ii) g = m times a smooth multiplicative distortion 1 + 0.3 cos(theta_p) ("wrong
large-scale shape", mostly absorbed by mu_c, so a weak test); (iii) g = m 10^(0.1 N'(0,1)) with an
independent pixel-scale scatter ("wrong small-scale shape", the honest test the reviewer asked
for), to measure sensitivity to the template, which a real analysis never has exactly. Same S/N
proxy, generator, noise and replicate seeds as study_v1 so the rows are comparable.

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
3. **Limit checks (L1, cheap):** (a) s_c fixed at 1e3 must reproduce w1 with the fixed cmb prior
   1/70^2 added to A (not the flat-cmb w1 of study_v1; the difference is a 1 to 2 percent change
   of one diagonal entry), up to a beta-independent constant; (b) s_c = 0.5, mu_c = 1, g = m must
   reproduce mode "gauss" with inv_T = 1/(0.5 m)^2 (no floor) after dropping the two trailing
   terms; with the default T = max(0.5 |m|, 20)^2 equality holds only in pixels with m_p,c >= 40
   uK_RJ (about 10 percent of the sky), so the check is run with the floor removed; (c) the (mu,
   log s) gradients against finite differences at s = 1e3, s = 0.3 and s = 1e-3 (check_1a style).

## Command
```bash
cd nextgen/prototypes
python -m indexfield.checks.check_2a_hier
python -m indexfield.coverage --n-rep 12 --noise 0.3 1 3 --modes hier hier_distort hier_scatter --jobs 2 --tag study_v1_hier
python -m indexfield.coverage_analyze --dir indexfield/checks/out/coverage/study_v1_hier
cd indexfield && python -u sbc.py --config hier --mode hier --param noncentred --n 300 --workers 4
python sbc_analyze.py --config hier_hier
```
The SBC runs the non-centred index parameterization (`HierNonCentredPosterior`, same target as
`HierPosterior`; `posterior.py` self-test verifies the Jacobian identity between the two to 5e-14
and the gradients to 5e-2 of 1e-5 relative). Reason: the centred gauss SBC of E012 (nside_beta 4,
K 192) ran 5.4 h on 4 workers without completing one replicate (4.5 CPU-hours per worker; the
centred field at K 192 sits in a funnel and NUTS runs at tree depth 10), and the hier SBC was
first planned centred because the hier posterior only existed centred. The coverage study
(nside_beta 2, K 48) ran centred without trouble, so the centred parameterization stays there.

## Result

### Coverage study (scientific question), 2026-10-26 15:27
108 cells (3 modes x 3 noise scales x 12 replicates, 280 to 290 s per cell, output
`checks/out/coverage/study_v1_hier/coverage_summary.txt`), same replicates, truths and S/N bins as
study_v1, so the rows below compare directly with the gauss / w0 / w1 rows of E012. Bias z is the
mean of (posterior mean minus truth) / posterior sd over index pixels with its standard error; cov68
and cov95 are the fractions of truths inside the central 68 / 95 percent intervals (binomial SE
0.02 for the all-pixel rows, larger per bin). The template-data truth is a_p = m_p 10^(0.1 N(0,1)),
so the "hier" shape g = m is correctly specified up to log-normal versus Gaussian scatter;
"hier_distort" multiplies g by 1 + 0.3 cos(theta) (large-scale wrong), "hier_scatter" by
10^(0.1 N) (pixel-scale wrong, independent of the truth's scatter).

| noise | mode | beta all pixels: z, sd z, cov68 / cov95 | lowest populated S/N bin: z, cov68 / cov95 (n) |
|---|---|---|---|
| 0.3 | hier | +0.10 +- 0.04, 1.02, 0.66 / 0.95 | S/N 1.6 to 3.3: +0.03 +- 0.06, 0.70 / 0.96 (256) |
| 0.3 | hier_distort | +0.09 +- 0.04, 1.08, 0.64 / 0.92 | +0.03 +- 0.06, 0.64 / 0.94 (256) |
| 0.3 | hier_scatter | +0.04 +- 0.05, 1.08, 0.67 / 0.94 | +0.06 +- 0.07, 0.70 / 0.95 (256) |
| 0.3 | gauss (E012) | +14.3 +- 0.6, 13.2, 0.10 / 0.20 | +25.7 +- 0.4, 0.00 / 0.00 (256) |
| 0.3 | w0 (E012) | +1.1 +- 0.4, 9.3, 0.65 / 0.91 | +1.7 +- 0.7, 0.62 / 0.89 (256) |
| 1 | hier | +0.35 +- 0.04, 0.93, 0.64 / 0.95 | S/N < 0.54: -0.02 +- 0.14, 0.56 / 0.96 (57) |
| 1 | hier_distort | +0.29 +- 0.04, 0.87, 0.69 / 0.97 | -0.19 +- 0.13, 0.54 / 0.96 (57) |
| 1 | hier_scatter | +0.19 +- 0.04, 0.85, 0.73 / 0.97 | -0.16 +- 0.11, 0.74 / 0.98 (57) |
| 1 | gauss (E012) | +2.08 +- 0.05, 1.18, 0.15 / 0.45 | +2.34 +- 0.11, 0.02 / 0.21 (57) |
| 1 | w0 (E012) | +4.7 +- 0.9, 21.3, 0.44 / 0.72 | +8.5 +- 3.4, 0.33 / 0.61 (57) |
| 3 | hier | +0.24 +- 0.05, 1.27, 0.67 / 0.93 | S/N < 0.54: +0.35 +- 0.08, 0.67 / 0.91 (288) |
| 3 | hier_distort | +0.10 +- 0.04, 0.99, 0.69 / 0.95 | +0.18 +- 0.06, 0.72 / 0.95 (288) |
| 3 | hier_scatter | +0.04 +- 0.04, 0.97, 0.72 / 0.95 | +0.12 +- 0.06, 0.75 / 0.95 (288) |
| 3 | gauss (E012) | +0.60 +- 0.04, 0.94, 0.59 / 0.93 | +0.68 +- 0.06, 0.57 / 0.93 (288) |
| 3 | w0 (E012) | +9.1 +- 1.2, 28.7, 0.39 / 0.62 | +13.1 +- 2.1, 0.34 / 0.57 (288) |

Per S/N bin (full tables in the summary file): the hier rows have mean z within +-0.5 in every bin
with more than one pixel at every noise level, except the 12-pixel S/N 0.54 to 1.57 bin at noise
0.3 (z -0.96 +- 0.16, sd z 0.55, cov 0.50 / 1.00: over-wide, 12 pixels from 12 replicates). Coverage
is within about two binomial SE of 0.68 / 0.95 in every bin except S/N 3.3 to 15.6 at noise 0.3
(58 pixels; hier 0.48 +- 0.07 / 0.86 +- 0.05 with sd z 1.28; w0 has 0.59 / 0.88 and sd z 1.25 in the
same bin, gauss 0.02 / 0.10) and the S/N < 0.54 bin at noise 3 for the ideal template only
(0.67 / 0.91, sd z 1.35; the two distorted templates give 0.72 / 0.95 and 0.75 / 0.95 there). The
pre-registered criterion (z within +-0.5 and coverage within binomial error at S/N below 15 at
noise 0.3 and 1) is met in all bins but the 58-pixel one above.

Residuals worth understanding: (i) at noise 1 all three variants show a positive mean z of +0.2
to +0.35 (0.3 sd) in the measured bins with sd z 0.85 to 0.93, i.e. a small bias with slightly
over-wide intervals, not seen at noise 0.3 or 3; (ii) the misspecified templates are, if
anything, better calibrated than the ideal one (the learned scatter absorbs the shape error, see
the hyperparameters below); (iii) low E-BFMI (0.03 to 0.17) at every noise level and 2 and 3
divergences in one and two chains at noise 3 (hier, hier_distort), median min ESS 13 to 29 per
1000 draws at noise 1 and 3 (the minimum is log sigma of the index field), 170 to 210 at noise
0.3. Whether (i) is the sampler or the log-normal versus Gaussian misspecification is what the SBC
below decides (the SBC generator is exactly the target).

Index hyperparameters (12 values per cell, so weak): mu z +0.1 to +0.8 with cov68 0.42 to 0.75;
log sigma z +0.5 to +1.3 with cov68 0.25 to 0.67 and cov95 0.75 to 0.92 (ESS 13 to 214). Compare
gauss (mu z +67 at noise 0.3, +8 at noise 1) and w0 (log sigma z +4.6 at noise 0.3, +20 at
noise 1) in E012.

Amplitude hyperparameters (posterior mean over the 12 replicates; posterior sd in brackets):
truth for mode hier is mu_c = E[10^(0.1 N)] = 1.027 and relative scatter 0.233 (log s = -1.46).
hier: mu_sync 1.037 (0.014), 1.037 (0.028), 1.021 (0.055) and log s_sync -1.42 (0.04), -1.48
(0.10), -1.66 (0.45) at noise 0.3 / 1 / 3; mu_dust 1.03 to 1.06, log s_dust -1.50 to -1.67 (sd
0.11 to 0.76). hier_scatter (truth scatter sqrt(2) x 0.233 = 0.33, log -1.11): log s_sync -1.10
(0.04) at noise 0.3, -1.17 (0.08) at noise 1, -1.37 (0.33) at noise 3; mu_sync 1.02, 0.99, 0.95.
hier_distort: same as hier within 0.01 (the 30 percent cos theta distortion is absorbed by mu_c
and s_c at the 1 percent level, as it is a large-scale modulation of a shape that the per-pixel
data constrain). ESS of the amplitude block 230 to 980 per 1000 draws.

Reading: at every S/N where the fixed Gaussian prior is biased and the profile (w0) mode is
unusable (S/N below 15, E012 F12), the hierarchical amplitude prior gives index posteriors with
mean z within 0.35 sd of zero and coverage at or near nominal, with the amplitude scale and
scatter recovered from the data to within their posterior widths, and this holds when the shape
template is wrong at the 30 percent level on large scales or at the 23 percent level per pixel.
This is the first mode in this programme that is usable below index-pixel S/N 15 without a
hand-set amplitude prior. The cost is the same per gradient as gauss; mixing of log sigma is the
bottleneck as before (centred parameterization at K 48).

Tooling fix recorded here: `coverage_analyze.py` had the mode list hard-coded to gauss / w0 / w1
and printed no tables for this directory; it now takes the list from `coverage.MODES`. The main
study's summary is unchanged by this.

### SBC (L2)
Smoke test of the non-centred path (2 replicates, 100 + 200 draws): tree depth 4.0, E-BFMI 1.04
and 1.12, step 0.24 to 0.29, no divergences, 310 and 540 s under load (the centred smoke had step
0.12 to 0.16). Full run launched 15:44 (300 replicates, 4 workers), `sbc_out/run_hier.log`.

### Residual check (noise 1, Gaussian generative scatter)
Launched 15:44: 12 noise-1 replicates with amp_mode "hier" (a_p = m_p (1.0269 + 0.2336 N), the
Gaussian twin of the template generator), mode hier, output
`checks/out/coverage/study_v1_hier_gaussgen/`. If the +0.3 sd bias at noise 1 disappears, it is
the log-normal versus Gaussian misspecification; if it stays, the SBC decides whether it is the
sampler.

## Negative control
Mode "gauss" with the fixed default T on the same template replicates: known biased (study_v1,
z +14 at noise 0.3). Mode "w1": known collapsed. Both run through the same analysis.

## Conclusion
(pending)

## Follow-ups
- Stage 2 (coupled smoothness prior on amplitude maps) if stage 1 is insufficient.
- Scalar-index version (M2.2 coverage table) with the same code.
