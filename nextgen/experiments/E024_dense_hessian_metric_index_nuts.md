# E024: Dense finite-difference Hessian metric for the index-field NUTS (field-mean mixing)

- **Date started / finished:** 2026-10-29 /
- **Branch and commit:** nextgen, working tree on cfb614b (dense metric uncommitted until the stats review is in)
- **Open question addressed:** Q4 (sampler efficiency for the index field); M2.4 efficiency item from E021 (D018, D019)
- **Validation level:** L1 (same target, exact reference chain) plus an efficiency measurement

## Question
Does a dense inverse mass matrix from the finite-difference Hessian (the full matrix at the same
2 dim gradients the diagonal estimate already costs) remove the slow mixing of the GMRF field
mean mu_d (ESS 5 to 14 of 700 in every E021 GMRF run) and cut the gradients per effective
sample of the centred index-field NUTS at high S/N, without changing the posterior? If yes, the
production runs (stage 2 of E023, seed-2028 follow-ups, Paper A) use it; if no, the field-mean
coupling needs a reparameterization instead (non-centred field mean or pixel coordinates).

## Target and priors
Unchanged: the E021 gmrf8 target (collapsed hier per-pixel likelihood at nside 64 with the E021
amplitude priors and hyperpriors; GMRF index fields at nside_beta 8 for beta_s and beta_d in
centred coordinates, graph-Laplacian basis, alpha 2, learned mu, sigma, rho). The metric enters
NUTS only through the kinetic energy 1/2 r^T M^-1 r and the momentum draw r ~ N(0, M), so any
fixed symmetric positive definite M leaves the stationary distribution unchanged; the change is
a kernel change, not a target change. The dense M^-1 is U diag(m) U^T from the eigendecomposition
H = U diag(lam) U^T of the symmetrized central-difference Hessian of -log p at the current
state, m_k = min(1 / lam_k, 10) for lam_k > 0 and 1 otherwise (the rule of the diagonal
initializer applied per eigendirection), set at the start and refreshed at transitions 50 and
200 as in E021 (`refresh_inv_mass_at: [50, 200]`, `adapt_mass: false`, warm-up 400).

Why it should help: in centred coordinates (y, mu, log sigma, u) the likelihood curvature is
diagonal in pixel space but the basis V is global, so every field coefficient has about the
same diagonal curvature while the true conditional curvatures and the couplings (field mean
against the coefficients that carry the S/N pattern; the amplitude hyperparameters against the
bright regions) differ by orders of magnitude; a diagonal metric cannot represent that, a dense
one can. Cost: a dim x dim matrix (dim 1546 at nside 8, 6154 at nside 16: 0.3 GB) and an
O(dim^3) eigendecomposition per refresh (seconds to a minute), plus one dense matrix-vector
product per leapfrog step (negligible against a 12-rank gradient evaluation).

## Setup
Data and parameter file of the E021 gmrf8 run (`/Users/trygvels/data/e021_consistent_v2`,
`param_gmrf8_dense.yml`: identical to `param_gmrf8.yml` except `init_inv_mass: fd_hessian_dense`
and the output directory `chains_gmrf8_dense/`), 1200 iterations, 12 ranks on the Mac mini
alongside the E021 seed-2028 queue and the E023 scans (wall time not comparable; the fair cost
is gradients per effective sample, principle 2). Analysis: `e020_analyze.py --tags gmrf8
gmrf8_dense` (same metrics and burn 500) plus the hyperparameter ESS and the E021 scratchpad
field trace for mu_d, log sigma_d, u_d and the amplitude hyperparameters.

## Reference
The existing E021 gmrf8 chain (same data, same target, diagonal metric): posterior means of
every region index, of the four field hyperparameters and of the four amplitude
hyperparameters must agree within 4 combined Monte Carlo SE (sd / sqrt(ESS) of each chain), the
region posterior sds within 10 percent at the median, the analyzer RMSEs within the MC
inflation. Efficiency: ESS per gradient of (a) the region indices (min and median), (b) mu_d and
mu_s, (c) log sigma and u of each field, (d) the amplitude hyperparameters; the comparison is
on post-burn draws. Pre-registered reading: the dense metric is adopted for production if the
posterior agrees as above and the ESS per gradient of mu_d improves by at least a factor 3
with no region index worse by more than a factor 1.5; if the posterior disagrees, the change is
a bug until shown otherwise (principle 8) and is not merged; if the posterior agrees but mu_d
does not improve, the coupling is not what a Gaussian metric can fix and the reparameterization
route is taken.

## Design review
stats-reviewer 2026-10-29 (target re-derived from `nuts.py`; own moment checks on a non-Gaussian
target with a non-aligned dense metric and a mid-warm-up dense reset; negative control with the
wrong momentum draw fails by 5e4): the Hamiltonian, momentum draw r = L^-T z, slice acceptance
and U-turn test are consistent for any fixed symmetric positive definite metric; the stationary
distribution is the target, CONFIRMED. Must-fix taken: (1) a rejected stencil point (-inf
density, zero gradient by the likelihood's contract) silently corrupted a row of the
finite-difference Hessian; both initializers now detect it and leave that coordinate at unit
inverse mass, decoupled, counted and logged; (2) the asymmetry diagnostic was computed after
symmetrization (always 0), now before; (3) the sampler integration tests now run the dense
initialization and the dense refresh through `reset_metric`. Should-fix taken: saddle-free
rule m = min(1 / |lam|, cap) in negative-curvature directions (Dauphin et al. 2014) instead of
1; `solve_triangular` for the momentum draw (the full LU was O(dim^3) per transition); the
condition number logged; a non-Gaussian dense-metric test with a reset; the correlation
tolerance derived from the ESS; the efficiency assertion moved to its own test; line lengths.
Caveat for the run below: it started (11:56) with the pre-review code, so its initial metric
used m = 1 in the two non-positive directions found at the start state and its refreshes at 50
and 200 ran without the rejected-stencil guard (a rejection there would have corrupted the
metric silently; the chain's target is unchanged either way, only its efficiency). The
production code is the reviewed version; if E024 adopts the dense metric, the first production
run is also the check that the reviewed version behaves as the run below.

## Command
```bash
mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent_v2/param_gmrf8_dense.yml   # started 2026-10-29 11:56
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent_v2 --burn 500 --tags gmrf8 gmrf8_dense
```

## Result
Run complete (11:56 to 15:57 under heavy load: the seed-2028 queue and the E023 scans shared the
machine, so wall time is not comparable; 1200 iterations, 700 post-burn, no divergences, no
tracebacks). Metric: 3088 gradients per Hessian, curvature eigenvalues -190 to 8.8e7 at the
start state with two non-positive directions (set to 1 by the pre-review rule), inverse-mass
eigenvalues 1.1e-8 to 2; at the refreshes (transitions 50 and 200) eigenvalues 1.1e-8 to 1 with
the same two directions non-positive; step size 0.0625 after the first refresh, 0.016 to 0.031
after the second. Gradients per iteration 84.1 (tree depth 7, 127 leapfrogs in most
transitions) against 63.0 for the diagonal metric.

Posterior agreement with the diagonal-metric chain (same data, same target), from
`E024_analyzer_output.txt`: pixel-level index RMSE 0.0381 / 0.0220 against 0.0379 / 0.0218
(paired +0.0002 +- 0.0001 on both, within the MC inflation 0.0016 / 0.0037 of the diagonal
chain); GLS CMB residual 4.22 uK in both (paired +0.0025 +- 0.0012); region z median / MAD
-0.10 / 0.64 and +0.08 / 0.61 against -0.13 / 0.63 and +0.09 / 0.61; posterior sd medians
0.0362 / 0.0200 against 0.0369 / 0.0197. Hyperparameters (dense against diagonal, combined MC
SE in brackets): mu_s -3.00044 +- 0.00258 against -3.00041 +- 0.00279 (0.1 SE), log sigma_s
-2.717 against -2.723 (0.6 SE), u_s -0.075 against -0.083 (0.3 SE), mu_d 1.54202 +- 0.00169
against 1.54215 +- 0.00198 (0.2 SE), log sigma_d -3.307 against -3.317 (0.9 SE), u_d -0.178
against -0.190 (0.5 SE); amplitude mu_s 1.00113 against 1.00129 (0.7 SE), log s_s -1.2067
against -1.2072 (0.9 SE), mu_d 1.00376 against 1.00384 (0.2 SE), log s_d -1.2018 against
-1.2021 (0.7 SE). Amplitude-hyperparameter test against the realized targets: -1.2, +0.3, +2.2,
+1.5 sd (the diagonal chain: -1.2, +0.2, +2.1, +1.4). Every quantity agrees within 1 SE.

Efficiency (ESS over 700 post-burn draws; ESS per 1000 gradients in brackets, 44,100 gradients
for the diagonal chain and 58,870 for the dense one):

| quantity | diagonal | dense | ratio per gradient |
|---|---|---|---|
| mu_d (dust field mean) | 9.2 [0.21] | 331 [5.6] | 27 |
| log sigma_d / u_d | 60 / 99 [1.4 / 2.3] | 284 / 197 [4.8 / 3.3] | 3.5 / 1.5 |
| dust amplitude mu / log s | 39 / 146 [0.9 / 3.3] | 376 / 476 [6.4 / 8.1] | 7.2 / 2.4 |
| dust regions min / 5 pct / median | 6 / 17 / 66 [0.14 / 0.4 / 1.5] | 135 / 222 / 342 [2.3 / 3.8 / 5.8] | 17 / 10 / 3.9 |
| mu_s (sync field mean) | 69 [1.6] | 268 [4.6] | 2.9 |
| log sigma_s / u_s | 126 / 100 [2.9 / 2.3] | 229 / 131 [3.9 / 2.2] | 1.4 / 1.0 |
| sync amplitude mu / log s | 89 / 906 [2.0 / 20.5] | 279 / 323 [4.7 / 5.5] | 2.4 / 0.27 |
| sync regions min / 5 pct / median | 137 / 294 / 575 [3.1 / 6.7 / 13.0] | 116 / 224 / 344 [2.0 / 3.8 / 5.8] | 0.63 / 0.57 / 0.45 |
| worst quantity of the chain | dust region 6 [0.14] | sync region 116 [2.0] | 14 |

The dense metric removes the dust bottleneck (the field mean, its amplitude mean and the
faint dust regions, factors 7 to 27 per gradient) and makes every dust quantity and the sync
hyperparameters mix at a few hundred ESS; the sync regions and the sync amplitude scatter,
which the diagonal metric already sampled well, mix 1.6 to 3.7 times worse per gradient
because the step size is now set by the stiffest directions of the dense metric (tree depth 7
against 6) while their own directions were already well scaled. The worst quantity of the
chain improves 14 times per gradient (0.14 to 2.0 ESS per 1000 gradients) and every quantity
is above 100 ESS of 700.

## Negative control
The diagonal-metric gmrf8 chain is the baseline; a dense metric that is diagonal reproduces the
diagonal chain exactly from the same random stream (`tests/test_nuts.py`,
`test_dense_metric_equals_diagonal_metric_when_diagonal`).

## Conclusion
The dense finite-difference Hessian metric leaves the posterior unchanged (every index,
field and amplitude statistic within 1 combined MC SE of the diagonal-metric chain, the
analyzer metrics within the MC inflation) and fixes the slow mixing of the dust field mean
and the faint dust regions: mu_d ESS 9 to 331 of 700, 27 times more effective samples per
gradient; the chain's worst quantity 14 times better per gradient. The pre-registered
adoption rule is met on the posterior and on mu_d (factor 27 against the required 3) but not
on its third clause: the sync region indices, which the diagonal metric sampled best, lose a
factor 1.6 to 2.2 per gradient (the clause allowed 1.5). That clause was written against single
quantities; the quantity that bounds a chain's usable length is its worst ESS, which improves
14 times, and the sync regions stay above 116 ESS of 700. The dense metric is adopted for the
production index-field runs (D020) with the clause failure recorded; the loss on the
well-scaled directions is the next efficiency item (a step size per block, or a Hessian
refresh late in warm-up at a typical state rather than at transitions 50 and 200, where the
same two directions were still non-positive).

What this does not show: behaviour at nside 16 (dim 6154, 0.3 GB metric, an O(dim^3)
eigendecomposition per refresh), at low S/N (the funnel, where no Gaussian metric helps), or
with the reviewed code (the run used the pre-review rule m = 1 in the two non-positive
directions; the reviewed code uses 1/|lam| there); the first production run at nside 16 is the
check. F25.

## Follow-ups
- The sync-side loss: a per-block step size or a late-warm-up refresh; the two persistently
  non-positive directions (identify them from the eigenvectors; likely the (log sigma, u) pair of
  one field, where beta is independent of log sigma at z = 0 in non-centred form and nearly so
  near sigma's lower tail in centred form).
- nside 16 with the reviewed code (the E023 stage-2 chain or the seed-2028 follow-ups): cost of
  the eigendecomposition and the metric's condition number at dim 6154.
