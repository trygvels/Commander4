# E019: NUTS versus MAMS on the hierarchical index-field target: effective samples per gradient (M2.3, Q4)

- **Date started / finished:** 2026-10-27 / 2026-10-27
- **Branch and commit:** nextgen, d5878b9
- **Open question addressed:** Q4 (which gradient sampler for the index field); ROADMAP M2.3
- **Validation level:** L1 for the new sampler (exact-target self-test), then an efficiency measurement (principle 2: cost in gradient evaluations)

## Question
On the E016 hier target (index field at nside_beta 2 and 4 with the hierarchical amplitude prior,
template skies at noise 1 and 0.3), does the Metropolis-adjusted microcanonical sampler (MAMS,
arXiv 2503.01707) give more effective samples per gradient evaluation than our NUTS, in the
non-centred and centred parameterizations? The answer selects the production sampler for M2.4
stage 2 (the Commander4 group currently runs NUTS); a factor below about 2 is not worth a second
sampler in the C++/MPI path, a factor above 3 is.

## Target and priors
The E016 hier target exactly as in `posterior.HierNonCentredPosterior` / `HierPosterior`
(non-centred and centred index block, hyperpriors mu ~ N(-3, 0.5^2), sigma ~ half-N(0, 0.3), rho
log-uniform on [h, pi], amplitude hyperpriors mu_c ~ N(1, 0.5^2), log s_c ~ N(log 0.3, 1^2), guard
log s >= -8), on one template replicate per (nside_beta, noise) from
`coverage.simulate_template(..., amp_mode="hier")` with seed 0 (generator hyperparameters drawn from
the hyperprior; dim = K + 7 = 54 at nside_beta 2, 198 at nside_beta 4). Both samplers target the
same density; the parameterization changes the geometry, not the target.

## Setup
`nextgen/prototypes/indexfield/bench_samplers.py`, nside_d 16, 500 warm-up + 2000 samples, NUTS
(ours, dual averaging to 0.8, windowed diagonal mass adaptation, max depth 10) and MAMS
(`indexfield/mams.py`: leapfrog, exact sphere kick, dual averaging to 0.9 acceptance, diagonal
preconditioner from a warm-up window, ALBA trajectory length with the 0.3 constant measured in
adjusted transitions, floored at 0.3 sqrt(d)). Mac mini under heavy load (wall times are not
comparable and are not used). Tag `bench_v1`, JSON per run under
`prototypes/indexfield/checks/out/bench/bench_v1/`.

## Reference
Correctness of MAMS before any efficiency claim: `python -m indexfield.mams` self-test: the kick
weight equals the log-Jacobian form (d-1) log[cosh(delta) + (e.u) sinh(delta)] to 1e-10 and keeps
|u| = 1; `rescale` gives E[n] = L/eps; a 20-d correlated Gaussian: all 20 means within 4 MC SE and
variances within 10 percent (PASS, 2026-10-27 02:40: max |z| 1.41, variance ratios 0.974 to 1.024,
0 divergences). The efficiency measure is Geyer's ESS per coordinate divided by the gradient
evaluations of the sampling phase (warm-up gradients reported separately); the figure of merit is
the minimum over the beta field (per index pixel) and over the hyperparameters. Both chains must
recover the truth (z rms about 1) and have no divergences.

## Command
```bash
cd nextgen/prototypes && python -m indexfield.mams
for s in nuts mams; do for nb in 2 4; do for noise in 1 0.3; do for par in noncentred centred; do
  python -m indexfield.bench_samplers --sampler $s --nside-beta $nb --noise $noise \
      --parameterization $par --n-warmup 500 --n-samples 2000 --tag bench_v1
done; done; done; done
```

## Result
**NUTS baseline (02:27 to 03:27):**

| nside_beta (K) | noise | param. | gradients | ESS beta min / median | ESS hyper min | ESS per gradient beta min / median, hyper min | div. |
|---|---|---|---|---|---|---|---|
| 2 (47) | 1 | nc | 29976 | 1115 / 2226 | 497 | 0.037 / 0.074, 0.017 | 0 |
| 2 (47) | 1 | c | 26116 | 16 / 82 | 12 | 0.0006 / 0.0031, 0.0005 | 0 |
| 2 (47) | 0.3 | nc | 30000 | 1262 / 1834 | 543 | 0.042 / 0.061, 0.018 | 0 |
| 2 (47) | 0.3 | c | 12054 | 45 / 117 | 18 | 0.0037 / 0.0097, 0.0015 | 0 |
| 4 (191) | 1 | nc | 30000 | 1031 / 1698 | 524 | 0.034 / 0.057, 0.018 | 0 |
| 4 (191) | 1 | c | 55320 | 88 / 226 | 12 | 0.0016 / 0.0041, 0.0002 | 0 |
| 4 (191) | 0.3 | nc | 30000 | 1392 / 2298 | 1308 | 0.046 / 0.077, 0.044 | 0 |
| 4 (191) | 0.3 | c | 29424 | 41 / 200 | 5 | 0.0014 / 0.0068, 0.0002 | 16 |

Non-centred NUTS: 0.034 to 0.046 ESS per gradient for the worst index pixel, flat in K from 47 to
191 (15 gradients per transition throughout); the centred parameterization is 10 to 60 times
worse on this target (the funnel of the field scale), with divergences at nside_beta 4, noise 0.3.
Smoke run at nside_beta 1, nside_d 8 (500 + 1000): NUTS 0.047, MAMS 0.052 ESS per gradient (beta
min).

**MAMS (03:29 to 05:58, same data and seeds, 500 + 2000 as NUTS; ESS per gradient counts the
sampling-phase gradients only, as for NUTS):**

| nside_beta (K) | noise | param. | gradients warm-up + sampling | ESS beta min / median | ESS hyper min | ESS per gradient beta min / median, hyper min | div. | MAMS / NUTS (beta min) | (beta median) |
|---|---|---|---|---|---|---|---|---|---|
| 1 (11) | 1 | nc | 14081 + 3182 | 165 / 229 | 113 | 0.052 / 0.072, 0.035 | 0 | 1.10 | 1.19 |
| 2 (47) | 1 | nc | 25139 + 5698 | 194 / 370 | 87 | 0.034 / 0.065, 0.015 | 0 | 0.92 | 0.87 |
| 2 (47) | 1 | c | 97239 + 7344 | 57 / 158 | 25 | 0.0077 / 0.022, 0.0035 | 0 | 12.7 | 6.9 |
| 2 (47) | 0.3 | nc | 41846 + 5264 | 209 / 386 | 101 | 0.040 / 0.073, 0.019 | 0 | 0.94 | 1.20 |
| 2 (47) | 0.3 | c | 110865 + 7749 | 27 / 71 | 14 | 0.0035 / 0.0092, 0.0019 | 0 | 0.93 | 0.95 |
| 4 (191) | 1 | nc | 169307 + 5952 | 55 / 327 | 145 | 0.0093 / 0.055, 0.024 | 0 | 0.27 | 0.97 |
| 4 (191) | 1 | c | 263894 + 11964 | 12 / 37 | 4 | 0.0010 / 0.0031, 0.0004 | 0 | 0.64 | 0.76 |
| 4 (191) | 0.3 | nc | 182882 + 8201 | 389 / 619 | 331 | 0.048 / 0.076, 0.040 | 0 | 1.02 | 0.99 |
| 4 (191) | 0.3 | c | 532954 + 23211 | 5 / 5 | 4 | 0.0002 / 0.0002, 0.0002 | 0 | 0.16 | 0.03 |

(The nside_beta 1 NUTS row, from the same driver: 21806 + 14608 gradients, ESS 689 / 886, 0.047 /
0.061, 0.033.) All MAMS runs recover the truth (z rms of the index pixels 0.67 to 0.99 on the
non-centred rows; 1.25 on the broken centred nside_beta 4, noise 0.3 run whose ESS is 5).

Reading: on the non-centred target, which is the one in use (F19, D-entries of M2.4), MAMS and NUTS
have the same efficiency: the median-pixel ratio is 0.87 to 1.20 in all five cells, the worst-pixel
ratio 0.92 to 1.10 in four cells and 0.27 in one (nside_beta 4, noise 1, where a single pixel has
ESS 55 of 2000 samples while the median pixel is at parity; with 2000 samples the minimum over 192
ESS estimates is itself noisy at the 30 to 50 percent level). On the centred funnel both are poor
and MAMS is erratic: 13 times better in one cell, 6 times worse in another, and its ALBA
adaptation did not converge at nside_beta 4, noise 0.3 (ESS 5). MAMS warm-up as implemented
(dual averaging, a 600-transition window for the preconditioner, ALBA, re-adaptation) cost 2 to
20 times the NUTS warm-up in gradients (97000 to 533000 against 9000 to 49000); this is excluded
from the ESS per gradient but matters in a Gibbs sampler that re-adapts, or must carry its
adaptation across iterations.

## Negative control
The centred parameterization is the control: the same target, known from E012 to mix badly in a
funnel; both samplers must show it. They do: 5 to 60 times fewer effective samples per gradient
than non-centred for NUTS, 4 to 240 times for MAMS, in every cell.

## Conclusion
On the hierarchical index-field target as sampled in Commander4 (non-centred field, diagonal
preconditioner, K 11 to 191 plus 7 hyperparameters), MAMS (arXiv 2503.01707, our NumPy port with
leapfrog, exact sphere kick, dual averaging and ALBA) gives the same effective samples per
gradient as our NUTS: ratios 0.9 to 1.2 at the median pixel in every cell, 0.3 to 1.1 at the worst
pixel, with a 2 to 20 times more expensive warm-up. This is far below the factor 2 that would make a
second sampler worth maintaining (the threshold stated before the runs) and nowhere near the 3 to
100 reported in the literature for very high-dimensional field targets. Decision for M2.3: NUTS
stays the production sampler for the index field; Q4 is answered for the dimensions of Phase 2
(nside_beta up to 4 to 8). Caveats recorded for a later revisit at nside_beta 16 and above
(d > 3000): one seed per cell; 2000 samples; our own adaptation schedule rather than the
reference implementation; leapfrog rather than the minimal-norm integrator; diagonal
preconditioning for both; wall time not compared (machine under load). The expected place where
MAMS could still pay is the joint amplitude-plus-index target of Q5 (d of order 1e5 to 1e7), which
is a different experiment.

## Follow-ups
- MAMS with Langevin noise (L_partial = 1.25 L) and the minimal-norm integrator (2 gradients per
  step) if the leapfrog version is competitive.
- Preconditioning beyond the diagonal (Fisher / low-rank, nutpie-style) for both samplers.
