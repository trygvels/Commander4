# E011: Mac mini build, test suite, and reproduction of E005 and E006 (M1.1, M1.2)

- **Date started / finished:** 2026-10-25 / 2026-10-25
- **Branch and commit:** nextgen, 4137c7f (plus the setup-script fix committed with this record)
- **Open question addressed:** none directly; Phase 1 milestones M1.1 and M1.2 in `../ROADMAP.md`
- **Validation level:** L0 (full suite), L1 (E005 toy vs grid), L3 (E006 end to end vs grid)

## Question
Does the Commander4 `nextgen` branch build and pass on the Mac mini, and do the cloud-session
results E005 and E006 (F1, F2, F3, F8) reproduce on a different machine, compiler and MPI stack?
If not, the findings would rest on one environment and the macOS build notes would be wrong.

## Target and priors
Same as E005 and E006, unchanged:
- E005 toy (`tests/test_mcmc_likelihood.py::PerPixelToy`): per-pixel d = M(beta) a + n, 7 bands,
  three components (flat, fixed power law, power law with the sampled scalar index beta), white
  noise sigma = 1, flat amplitude prior, Gaussian prior beta ~ N(1.5, 0.3^2). `marginal` and
  `conditional` target the flat-amplitude-prior marginal p(beta | d); `ridge` targets the profile
  likelihood (no -1/2 ln det A term); `resolve_chisq` is the legacy move with no valid target.
- E006 (`aux/validate_mcmc_likelihood.py`): nside 16, 7 bands 30 to 353 GHz, common 480 arcmin
  beam, CMB + synchrotron (beta_s, nu_ref 30 GHz) + dust (beta_d, T = 20 K, nu_ref 545 GHz),
  flat priors on beta_s in [-4, -2] and beta_d in [1.0, 2.5], flat amplitude priors; exact targets
  on a bounded grid for both log-det settings.

## Setup
Mac mini, Apple M4 (10 cores), 24 GB, macOS 27.0, Apple clang 17, CPython 3.13.11, MPICH 5.0.2
(pip), see `../COMPUTE.md`. Build via `nextgen/scripts/setup_macos.sh` after fixing it (it assumed
Homebrew-installed tools; none were needed). E006: 7 MPI ranks on 10 cores, not oversubscribed,
250 Gibbs iterations, 80 burn-in, 50 MH steps per iteration per mode, modes ridge, marginal,
conditional, resolve_chisq. Simulation seed 2026 as in E006.

## Reference
E005: grid quadrature of the exact posterior from an independent NumPy implementation
(`PerPixelToy._logpost`). E006: coarse-to-fine bounded grid from an independent NumPy marginal
(`_exact_posterior` in the harness). Cloud-session numbers from `registry.md` are the comparison.

## Command
```bash
bash nextgen/scripts/setup_macos.sh
source .venv/bin/activate
python -m pytest -q -n auto -p no:cacheprovider
python nextgen/experiments/E011_mac_mini_repro/e005_repro.py
python aux/validate_mcmc_likelihood.py make --dir $V
for m in ridge marginal conditional resolve_chisq; do mpirun -n 7 commander4 -p $V/param_$m.yml; done
python aux/validate_mcmc_likelihood.py analyze --dir $V
```

## Result
**Build and tests (M1.1):** build succeeded on the first run of the fixed script. Focused tests 24
passed in 35 s; full suite 859 passed, 4 skipped, 0 failed in 25 s (10 xdist workers).

**E005 (M1.2), toy chains vs exact grid posterior, seed 7, 20% burn-in:**

| npix | likelihood | n_iter | sd ratio | mean shift / sd | cloud (E005) |
|---|---|---|---|---|---|
| 150 | marginal | 1500 | 0.99 | -0.02 | matches exact |
| 150 | ridge | 1500 | 0.99 | -0.03 | matches exact |
| 150 | conditional | 2500 | 1.08 | -0.09 | matches exact |
| 150 | resolve_chisq | 2500 | 3.22 | +0.50 | 3.2x, +0.50 sd |
| 600 | marginal | 1500 | 1.03 | -0.04 | |
| 600 | ridge | 1500 | 1.01 | -0.04 | |
| 600 | resolve_chisq | 2500 | 4.35 | +0.83 | 4.4x, +0.83 sd |

Identical to the cloud numbers to the quoted precision (same seeds, deterministic kernels; the
C++ kernel compiled with Apple clang gives the same chains as the Linux build).

**E006 (M1.2), end to end, 170 post-burn-in samples per mode:**

Exact targets (grid): flat-amplitude-prior posterior beta_sync -2.0080 +- 0.0089 (on the upper
bound), beta_dust 2.4995 +- 0.0006 (on the upper bound); ridge target beta_sync -2.7523 +- 0.1314,
beta_dust 1.1834 +- 0.1128 (1.9 sd and 3.7 sd from the input, F8). Both identical to E006.

| likelihood | sync shift / sd | sync sd ratio | dust shift / sd | dust sd ratio | chain mean +- sd | wall time |
|---|---|---|---|---|---|---|
| ridge | +0.03 | 0.94 | -0.02 | 0.97 | -2.749 +- 0.124, 1.181 +- 0.109 | 1 min 39 s |
| marginal | -165 | 17.5 | +0.23 | 0.62 | -3.474 +- 0.156, 2.4996 +- 0.0004 | 1 min 17 s |
| conditional | -0.23 | 2.01 | -0.04 | 0.74 | -2.010 +- 0.018, 2.4994 +- 0.0004 | 1 min 17 s |
| resolve_chisq | -115 | 36 | -1730 | 522 | -3.027 +- 0.318, 1.531 +- 0.292 | 5 min 7 s |

Shifts and ratios are against each mode's own exact target (ridge target for `ridge`, flat-prior
posterior for the others). Ridge matches E006 exactly (cloud: within 0.03 sd, ratios 0.94 / 0.97).
`marginal` is again stuck in the local mode near beta_s = -3.5 (F3; the exact flat-prior target sits
on the hard bounds). `conditional` finds the bound mode (beta_s = -2.010, beta_dust at 2.5) with a
sd 2x the exact (the chain is clipped against the bound and mixes slowly there); this mode was not
tabulated in E006 and is consistent with F4 (conditional Gibbs carries the same volume pull as the
flat-prior marginal). `resolve_chisq` wanders far from any target (sd 0.32 and 0.29 versus 0.009 and
0.0006), as in E006.

## Negative control
`resolve_chisq` in both E005 and E006: failed as expected (see tables).

## Conclusion
M1.1 and M1.2 are met. The macOS build is verified and the setup script is fixed and documented
(COMPUTE.md). F1, F2, F3 and F8 reproduce on an independent toolchain (Apple clang, MPICH 5 arm64)
and are not artefacts of the cloud environment. Wall time per E006 mode is about 80 to 100 s (5 min for resolve_chisq, which re-solves amplitudes at every step) on 7
unoversubscribed ranks, versus about 15 min on 2 oversubscribed cores in the cloud.

## Follow-ups
- Proceed to the independent `stats-reviewer` pass on E009/E010, then Q14 (STATUS.md).
- M1.3: fetch low-resolution maps within the ~100 GB disk budget.
