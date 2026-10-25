# Spectral-index sampler benchmark

Toy test of whether gradient-based sampling improves on Commander3's spectral-index sampling.

## Model
Per-pixel intensity model in uK_RJ over 8 bands (23 to 353 GHz):
CMB + synchrotron power law (beta_s) + modified blackbody dust (beta_d, T_d = 20 K) + white noise.
Flat amplitude priors, Gaussian priors beta_s ~ N(-3.0, 0.2), beta_d ~ N(1.6, 0.1).
Pixels span faint (high latitude) to bright (plane). No beams, no smoothing, no TOD.

## Samplers
| name | what it is |
|---|---|
| MH chisq k=10 / k=100 | Commander3 default: exact amplitude draw, then per-index random-walk MH conditional on amplitudes, proposal length tuned with Commander's 0.3 to 0.7 acceptance rule (comm_nonlin_smod.f90) |
| MH marginal k=10 | Commander3 'marginal' lnL: each index sampled with its own amplitude marginalized, amplitudes redrawn after each index |
| Gibbs + NUTS | same Gibbs loop, MH replaced by one NUTS step on (beta_s, beta_d) given amplitudes |
| Joint NUTS | NUTS over all 5 parameters, no amplitude solves |
| Collapsed NUTS | all amplitudes marginalized analytically, NUTS on the indices, exact amplitude draw per sample |

Gradients come from ForwardDiff on the same Julia functions that evaluate the likelihood. NUTS is AdvancedHMC.jl.

## Run
```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
NPIX=300 NOISE=1.0 julia --project=. -t auto run.jl
```
`NOISE` scales the noise rms (0.3 = high S/N, 3.0 = low S/N). Output: `summary_noise*.txt`, `results_noise*.csv`.

## Results (300 pixels, 2-core x86 VM, Julia 1.13.2)
Amplitude solves per effective sample (median over pixels, worst of 5 parameters; lower is better).
In Commander3 the amplitude step is the expensive global CG solve, so this is the number that transfers.

| sampler | noise x0.3 | noise x1 | noise x3 |
|---|---|---|---|
| MH chisq k=10 (Commander default) | 41.8 | 13.8 | 3.2 |
| MH chisq k=100 | 41.2 | 14.5 | 3.2 |
| MH marginal k=10 | 20.0 | 7.8 | 3.4 |
| Gibbs + NUTS | 41.4 | 13.2 | 3.1 |
| Collapsed NUTS | 1.9 | 1.5 | 1.3 |
| Joint NUTS | 0 solves, 126 / 55 / 42 gradients per ESS | | |

Effective samples per 1000 iterations: 30 / 91 / 395 for MH chisq versus 517 / 672 / 747 for collapsed NUTS.
All samplers agree on posterior means to within ~0.03 sigma, and 68/95% coverage of the true indices is nominal.

## Caveats
- Pixels are independent here. In Commander3 beams, smoothing and the CMB prior couple pixels, so the
  amplitude solve is a global CG and full collapse is not free. Realistic transfers: collapse per region or
  at low resolution (as the 'marginal' lnL already does), or joint HMC with adjoint-SHT gradients.
- Wall-clock numbers compare a hand-written MH loop to AdvancedHMC, so per-iteration and per-solve
  metrics are the fair ones.
- The baseline is a reimplementation of Commander3's sampler logic, not a run of Commander3 itself.

## Commander4-style move (`c4test.jl`, `c4global_check.jl`)
Commander4 (compsep/mcmc.py, spectral_index.py) samples indices by proposing beta', re-solving the
amplitudes as a constrained realization (`optimize: false`, or the per-pixel solver which always adds
fluctuations), and accepting on the change in data chi^2 plus prior. This toy reproduces that move.

```
julia --project=. -t auto c4test.jl          # per-pixel and full-sky scalar beta
julia --project=. -t auto c4global_check.jl  # 4 chains x 12000 iterations, full-sky beta
```

Results:
- Per-pixel indices (300 pixels): beta posteriors 24% too wide (median sd ratio 1.24), 32 amplitude
  solves per effective sample (collapsed NUTS: 1.5).
- Full-sky scalar indices, per-pixel amplitudes (Commander4's `perpix_betasamp` setup):

| npix | beta_s sd ratio (C4 / exact) | beta_s mean shift |
|---|---|---|
| 256 | 3.9 | -0.43 sd (+-0.09) |
| 1024 | 5.7 | -0.99 sd (+-0.20) |

The shifted means match the profile ("ridge") posterior, i.e. the exact posterior without the
amplitude log-determinant term.

Why: the correct acceptance ratio for "propose beta', then draw a' ~ p(a | beta', d)" is the ratio of
marginal posteriors p(beta' | d) / p(beta | d), which does not depend on the drawn amplitudes. Using
chi^2 at the drawn amplitudes adds the random term (X' - X), with X ~ chi^2_k and k the number of
amplitude degrees of freedom, and drops the log det term. The noise in the acceptance grows as sqrt(k),
so the posterior inflation grows with resolution (roughly as k^(1/4) here).

Independent check (`global_gibbs_check.jl`): Commander3-style Gibbs, which shares no math with the
collapsed formula, reproduces the exact full-sky posterior (mean shifts <= 0.06 sd, sd ratio 0.95 to 1.00
at 64 and 256 pixels), while the Commander4-style move at 256 pixels is 3.9x too wide.

Valid alternatives inside Commander4:
- `update_amplitude_groups: []` (propose against fixed amplitudes, Commander3 'chisq' style): valid, slow.
- With the per-pixel solver, accept on the analytic marginal likelihood
  0.5 b^T M^-1 b - 0.5 log det M per pixel (exact, same cost), or sample it with NUTS ("Collapsed NUTS").
