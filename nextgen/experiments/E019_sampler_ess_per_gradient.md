# E019: NUTS versus MAMS on the hierarchical index-field target: effective samples per gradient (M2.3, Q4)

- **Date started / finished:** 2026-10-27 /
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

**MAMS (running, queued after the NUTS runs):** (pending)

## Negative control
The centred parameterization is the control: the same target, known from E012 to mix badly in a
funnel; both samplers must show it.

## Conclusion
(pending)

## Follow-ups
- MAMS with Langevin noise (L_partial = 1.25 L) and the minimal-norm integrator (2 gradients per
  step) if the leapfrog version is competitive.
- Preconditioning beyond the diagonal (Fisher / low-rank, nutpie-style) for both samplers.
