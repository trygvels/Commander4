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

## Command
```bash
mpirun -n 12 commander4 -p /Users/trygvels/data/e021_consistent_v2/param_gmrf8_dense.yml   # started 2026-10-29 11:56
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e021_consistent_v2 --burn 500 --tags gmrf8 gmrf8_dense
```

## Result
(pending)

## Negative control
The diagonal-metric gmrf8 chain is the baseline; a dense metric that is diagonal reproduces the
diagonal chain exactly from the same random stream (`tests/test_nuts.py`,
`test_dense_metric_equals_diagonal_metric_when_diagonal`).

## Conclusion
(pending)

## Follow-ups
