# E017: Scalar (full-sky) spectral index: coverage of profile, flat-prior marginal, fixed Gaussian and hierarchical amplitude priors on template data (M2.2)

- **Date started / finished:** 2026-10-26 /
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
(running)

## Negative control
`w1` on template data: E012/F12/F15 predict a bias from the flat-prior volume term that does not
shrink with S/N in sd units (E014: 54 sd on real data); `gauss` with the misspecified fixed prior
is the second control.

## Conclusion
(pending)

## Follow-ups
- Same table with the real-data configuration of E014 (17 bands, nside 64 mask) once the
  prototype can read Commander4 parameter files, or directly in Commander4 once the hier prior is
  implemented (M2.4).
