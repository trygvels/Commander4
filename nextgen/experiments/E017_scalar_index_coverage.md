# E017: Scalar (full-sky) spectral index: coverage of profile, flat-prior marginal, fixed Gaussian and hierarchical amplitude priors on template data (M2.2)

- **Date started / finished:** 2026-10-26 / 2026-10-27
- **Branch and commit:** nextgen, b3113aa
- **Open question addressed:** Q2 (prior calibration), ROADMAP M2.2 "full-sky" row
- **Validation level:** L2 (frequentist coverage against known truth; exact 1D grids for three of
  the four modes)

## Question
For one synchrotron index over the whole sky (the Commander4 scalar case of E014, where the
`ridge` and `marginal` chains sit 54 ridge sd apart on real data), which of the four amplitude
treatments gives posteriors that cover the truth, as a function of full-sky S/N, on the same
template skies as the spatially varying study (E012 study_v1, E016)? If profile (ridge) and hier
are calibrated and the flat-prior marginal is not, the 54 sd gap of E014 is a bias of `marginal`,
and the ridge default for scalar indices (D-entry pending) is justified on simulations; if the
profile is biased at low S/N (as it is below index-pixel S/N 15 in E012), the hier prior is the
candidate default for scalar indices as well.

## Target and priors
Data: N = 3072 pixels (nside_d 16), B = 12 bands (simulate.BANDS_GHZ), three components (cmb,
sync, dust); d_p = M_p(beta) a_p + n_p, n_p ~ N(0, N_p) diagonal with the simulate.RMS_UK_RJ
levels times a noise scale; mixing matrix of E012 (sync power law with pivot 30 GHz, dust fixed
beta_d = 1.6, T_d = 20 K). One index: region[p] = 0 for all p, n_beta = 1.

Index prior (all modes): beta ~ N(-3, 0.2^2), equal to the generator. At scale 3 the posterior sd of
beta is 0.15, so the prior (sd 0.2) contributes there; it is the generator's distribution, so the
coverage test stays well specified; at the other scales it is immaterial.

Amplitude treatments, all through `model.collapsed_loglike` / `model.hier_loglike` with the
conventions stated there (the dropped constants are beta-independent and identical across the
grid points of a mode):
- `w0` (profile, Commander4 `ridge`): log p(beta | d) = 1/2 b^T A^-1 b + log p(beta) + const.
- `w1` (flat amplitude prior, Commander4 `marginal`): 1/2 b^T A^-1 b - 1/2 log det A + log p(beta).
- `gauss` (fixed Gaussian prior a_p ~ N(m_p, T_p), m = smooth templates, T^(1/2) =
  max(0.5 |m|, 20) for sync and dust, 70 for cmb; SPEC "default"): 1/2 b'^T A'^-1 b'
  - 1/2 log det A' + log p(beta). Misspecified for the template generator (as in E012).
- `hier` (E016 target with a scalar index): a_p,c ~ N(mu_c g_p,c, (s_c g_p,c)^2) for sync and
  dust with g = the smooth templates, a_p,cmb ~ N(0, 70^2); mu_c ~ N(1, 0.5^2), log s_c ~
  N(log 0.3, 1^2); joint target p(beta, mu, log s | d) proportional to L(beta, mu, s) p(beta)
  p(mu) p(log s) with L of E016 (reviewed 2026-10-26). Also `hier_distort` (g times
  1 + 0.3 cos theta) as in E016.

For `w0`, `w1`, `gauss` the posterior is one-dimensional and is computed on a grid (exact up to
quadrature, no sampler). For the hier modes it is five-dimensional (beta, mu_sync, mu_dust,
log s_sync, log s_dust) and is sampled with NUTS (`sampler.run_nuts`).

## Setup
Generator (coverage.simulate_template logic with one region): beta_true ~ N(-3, 0.2^2);
a_p = m_p 10^(0.1 N(0,1)) for sync and dust (log-normal scatter, mean 1.027 and relative sd 0.234
around the smooth templates m), a_p,cmb ~ N(0, 70^2); noise scales 0.1, 0.3, 1, 3, giving full-sky
S/N proxies (the snr_index of coverage.py summed over the sky) of about 710, 220, 65 and 23, i.e.
from the top of the E012 study_v1 range (index-pixel S/N 15, where the profile became
calibrated) upward by a factor 50; the posterior sd of beta runs from about 0.15 at scale 3 to a
few 1e-3 at scale 0.1. 24 replicates per noise scale, the
four modes plus hier_distort on the same data per replicate (seeded by replicate and noise scale).
Grids: coarse 3001 points on [-4.5, -1.5], then a fine grid of 1601 points on the coarse maximum
plus or minus 8 coarse-curvature sd, widened until it holds 8 posterior sd on each side. NUTS:
500 warm-up, 1000 samples, target acceptance 0.8, max depth 10, random init from the priors.
Machine: Mac mini, 2 worker processes (the E016 SBC is running on 4).

## Reference
Truth is the generator input. For `w0`, `w1`, `gauss` the posterior is exact (1D grid, spacing
below 1/50 of the posterior sd). For hier, a sampler check before the study: with the amplitude
hyperpriors pinned (mu_c ~ N(mu*, 1e-4^2), log s_c ~ N(log s*, 1e-4^2) at the generator values
mu* = 1.027, s* = 0.234) the hier posterior of beta must equal the gauss-mode grid with
m = mu* g, T = (s* g)^2 (same replicate): mean within 0.1 sd and width within 10 percent at
ESS above 300. Calibration statistics as in E012/E016: z = (mean - truth) / sd, 68 and 95
percent coverage, with replicates independent (every cell has its own truth and noise), so the
plain SE over replicates is the right one.

## Command
```bash
cd nextgen/prototypes && python -m indexfield.scalar_coverage --selftest      # pinned-hier check
cd nextgen/prototypes && nohup python -m indexfield.scalar_coverage --n-rep 24 \
    --noise 0.1 0.3 1 3 --jobs 2 --tag study_v1 > indexfield/checks/out/scalar/study_v1.log 2>&1 &
cd nextgen/prototypes && python -m indexfield.scalar_coverage_analyze --dir indexfield/checks/out/scalar/study_v1
```

## Result
**Self-test, first pass (18:30, 500+1000 NUTS, pinned hyperpriors):** on the noise-3 replicate
(truth -3.151, full-sky S/N 23) the pinned hier posterior gives -3.098 +- 0.123 against the gauss
grid -3.083 +- 0.119: width ratio 1.03, mean difference -0.13 grid sd at ESS 196 (Monte Carlo
error 0.07 sd), inconclusive at the 0.1 sd criterion; rerun with 3000 samples. On the same
replicate the flat-prior marginal `w1` has its maximum at the upper edge of the support (-1.5): the
F12/F15 collapse appears for a scalar index at S/N 23. `grid_posterior` now reports such cells
truncated to the support with an `at_edge` flag instead of failing. The grid modes of the study
were launched at 18:52 (2 workers).

**Self-test, second pass (500+3000 NUTS, 40 min under load 28): PASS.** Pinned hier
-3.0849 +- 0.127 against the gauss grid -3.0827 +- 0.119: mean difference -0.02 grid sd at ESS
618, width ratio 1.07, 0 divergences. On the same replicate: w0 z +0.30, gauss z +0.95, w1 at the
edge of the support (mean -1.5008, z +1704). The hier modes (hier, hier_distort) were launched at
19:38 (2 workers).

**Grid modes (288 cells, 81 min, done 20:13)** and the first 5 hier replicates per noise scale,
24 replicates per grid cell (`scalar_summary.txt` of 20:15, provisional for hier):

| noise (S/N) | w0 z, cov68/95 | w1 z, shift/w0 sd | gauss z, shift/w0 sd | hier z, cov68/95 (n) |
|---|---|---|---|---|
| 0.1 (680) | +0.13 +- 0.17, 0.75/1.00 | +6.0, +5.8 | +5.2, +4.9 | +0.07 +- 0.40, 0.83/1.00 (6) |
| 0.3 (230) | +0.35 +- 0.19, 0.50/1.00 | +17.7, +16.2 | +10.1, +8.6 | -0.14 +- 0.60, 0.00/1.00 (5) |
| 1 (67) | -0.20 +- 0.21, 0.54/0.96 | at edge (all 24), +19.7 | +4.4, +3.5 | one chain broken (see below) |
| 3 (22) | -0.34 +- 0.40, sd z 1.98, 0.21/0.58 | at edge (all 24), +10.0 | +0.4, +0.7; 0.46/0.88 | -0.31 +- 0.30, 0.80/1.00 (5) |

Reading so far: the flat-prior marginal is unusable for a scalar index at every S/N (its bias in sd
grows with S/N from 6 sd at S/N 680 to 18 sd at S/N 230, and below S/N 70 the posterior piles up
at the upper boundary of the support, the F12/F15 collapse); the fixed Gaussian prior is 4 to 10 sd
off wherever the index is well measured; the profile is calibrated at S/N 680 but at S/N 22 its
posterior is two times too narrow (sd z 1.98, coverage 0.21 / 0.58) and at S/N 230 and 67 its 68
percent coverage is 0.50 and 0.54 (+- 0.10); 72 more w0 replicates are running to settle that. The
hier amplitude hyperparameters are recovered at the high-S/N scales (mu 1.027 / 1.027 +- 0.01,
log s -1.44 / -1.39 +- 0.03 against the generator 1.027, -1.45).

**w0 with 96 replicates (extension of 72, done 21:05):** z +0.09 +- 0.09 / +0.08 +- 0.09 /
-0.09 +- 0.12 / +0.02 +- 0.19 at noise 0.1 / 0.3 / 1 / 3, sd z 0.93 / 0.89 / 1.14 / 1.90, coverage
0.72 / 0.98, 0.72 / 0.98, 0.56 / 0.93, 0.33 / 0.66 (binomial SE 0.05). So the profile is calibrated
at full-sky S/N 680 and 230, its posterior is about 15 percent too narrow at S/N 68, and two times
too narrow at S/N 23 (no mean bias at any scale). The 24-replicate cov68 of 0.50 at noise 0.3 was a
fluctuation. Hier at noise 3 on the first 10 replicates: z +0.04 +- 0.29, sd z 0.91, coverage
0.80 / 0.90.

**Broken hier chain (rep 4, noise 1):** 983 divergences, E-BFMI 0.03, step size 1.9e-14, posterior
"sd" 8e-13, log s_dust at -19. A numerical floor, not a property of the target: at that log s the
prior precision 1/(s g)^2 is of order 1e16 and swamps M^T N^-1 M in double precision, so the density
is rounding noise and every step diverges; the chain got there in early warm-up (large trial step)
and could not return. The other four noise-1 chains have log s_dust -1.9 to -2.0 +- 0.6 (dust is
barely detected at noise 1, so its scatter is prior-dominated). Fix (2026-10-26 20:30): a guard
log s_c >= -8 in both hier posteriors (`posterior.HIER_LOGS_MIN`; hyperprior mass below it 5e-12,
so the target is unchanged to that level; states below are reported as zero density); the
analyzer lists chains with divergences or E-BFMI below 0.3; the cell was rerun with the guard. The
E016 SBC (42 replicates by then, all with 0 divergences) runs the unguarded code; the guard would
not have acted on any of them. Rerun of rep 4 with the guard: z +0.32, sd 0.074, ESS min 294,
0 divergences, E-BFMI 0.97, log s_dust -1.39 (113 s).

**Final table (2026-10-27 05:30, all 768 cells: 96 w0, 24 each of the other modes per noise
scale; `scalar_summary.txt`).** Eight hier cells (3 at noise 1, 5 at noise 3, six of them
hier_distort) had collapsed at the unguarded log s floor (log s_dust or log s_sync at -18 to -20,
step sizes 1e-14 to 1e-16, posterior sd 1e-13): the study workers had been started at 19:38, before
the guard of 20:30, so they ran the old code. The collapsed records were moved to
`study_v1/collapsed_unguarded/` and the eight cells rerun with the guarded code on the same seeded
data (same truth to 1e-12); all eight came back healthy (0 or 1 divergence, E-BFMI 0.87 to 1.13,
ESS 176 to 745, log s -1.1 to -2.3). The hier rows below use the reruns.

| noise (S/N) | w0 z, sd z, cov68/95 (n 96) | w1 z, shift/w0 sd | gauss z, shift/w0 sd, cov68/95 | hier z, sd z, cov68/95, shift/w0 sd | hier_distort z, cov68/95 |
|---|---|---|---|---|---|
| 0.1 (680) | +0.09 +- 0.09, 0.93, 0.72/0.98 | +6.0, +5.8 | +5.2, +4.9, 0.00/0.00 | +0.10 +- 0.17, 0.82, 0.75/1.00, -0.02 | +0.08, 0.71/1.00 |
| 0.3 (230) | +0.08 +- 0.09, 0.89, 0.72/0.98 | +17.7, +16.2 | +10.1, +8.6, 0.00/0.00 | +0.31 +- 0.19, 0.92, 0.50/1.00, -0.04 | +0.31, 0.58/1.00 |
| 1 (68) | -0.09 +- 0.12, 1.14, 0.56/0.93 | at edge (all 24), +19.7 | +4.4, +3.5, 0.00/0.00 | -0.12 +- 0.16, 0.78, 0.83/1.00, +0.09 | -0.13, 0.79/1.00 |
| 3 (23) | +0.02 +- 0.19, 1.90, 0.33/0.66 | at edge (all 24), +10.0 | +0.4, +0.7, 0.46/0.88 | -0.07 +- 0.23, 1.13, 0.67/0.88, +0.29 | -0.06, 0.67/0.88 |

Binomial SE of the coverages: 0.05 (w0), 0.10 (others). Hier amplitude hyperparameters, mean over
replicates (generator mu 1.027, log s -1.45): mu_sync 1.03 / 1.03 / 1.02 / 1.03, mu_dust 1.02 /
1.03 / 1.03 / 0.98, log s_sync -1.42 / -1.42 / -1.44 / -1.65, log s_dust -1.44 / -1.42 / -1.84 /
-1.57 at noise 0.1 / 0.3 / 1 / 3; at the two noisy scales dust is barely detected and the log s
posteriors are prior-dominated (median sd 0.6 to 0.8). Cost: hier 83 to 100 s per chain (500 +
1000 NUTS, median ESS 340 to 1200), the grids 1 to 2 s.

## Negative control
`w1` on template data: E012/F12/F15 predict a bias from the flat-prior volume term that does not
shrink with S/N in sd units (E014: 54 sd on real data); `gauss` with the misspecified fixed prior
is the second control.

## Conclusion
For a single full-sky index on template data (F18):
- The flat-prior marginal (`marginal` in Commander4) is biased at every S/N: +6 sd at full-sky S/N
  680, +18 sd at 230, and below S/N about 70 the posterior collapses to the upper edge of the
  support (every replicate at noise 1 and 3), the scalar version of the F12/F15 collapse. It is
  the mechanism behind the 54 sd (pass 2) and 16 to 23 sd (passes 3 and 4) marginal-minus-ridge
  gaps of E014.
- The fixed Gaussian prior with a misspecified template is 4 to 10 sd off wherever the index is
  well measured and only looks acceptable (z +0.4, coverage 0.46 / 0.88) where its width is
  irrelevant.
- The profile (ridge, Commander4's default per D003) is unbiased at every S/N (|z| < 0.1 +- 0.1 to
  0.2) and calibrated at S/N 680 and 230, but its posterior is 15 percent too narrow at S/N 68
  (coverage 0.56 / 0.93) and about two times too narrow at S/N 23 (sd z 1.90, coverage 0.33 /
  0.66): the profile ignores the amplitude uncertainty that feeds the index.
- The hierarchical prior is unbiased at every S/N (|z| <= 0.31 +- 0.2), with sd z 0.8 to 1.1 and
  coverage 0.50 to 0.83 / 0.88 to 1.00 (binomial SE 0.10; the 0.50 at noise 0.3 is a 1.8 sigma low
  and comes with 95 percent coverage 1.00), the same with a 30 percent large-scale template
  distortion, and it recovers the amplitude hyperparameters where they are detected. Its mean is
  within 0.3 w0 sd of the profile mean, so on real data (E014) hier and ridge would agree closely
  while marginal does not.
So for scalar indices the ridge default is justified at high S/N and the hier prior is the
candidate wherever the full-sky S/N drops below about 100 or widths matter; this answers the
full-sky row of M2.2. The prior-dependence of real-data indices (E014, F17) is a marginal-mode
artefact plus the zero-level convention, not a property of the profile. Caveat: one template
family (the E012 smooth templates with log-normal scatter), 24 replicates per cell for the
non-profile modes.

## Follow-ups
- Same table with the real-data configuration of E014 (17 bands, nside 64 mask) once the
  prototype can read Commander4 parameter files, or directly in Commander4 once the hier prior is
  implemented (M2.4).
