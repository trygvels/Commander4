# Experiment registry

One row per experiment, newest at the bottom. Each experiment from E009 on gets its own record
`EXXX_<slug>.md` made from `EXPERIMENT_TEMPLATE.md` (`/new-experiment` scaffolds it). E001 to E008
were run on 2026-10-25 in the exploratory cloud session; their details are summarized below and in
`../FINDINGS.md`.

| ID | Date | Question | Where | Result | Findings |
|---|---|---|---|---|---|
| E001 | 2026-10-25 | Does gradient-based or collapsed sampling beat Commander-style Gibbs + MH for spectral indices? | `prototypes/specidx-bench/run.jl` | Collapsed NUTS needs 1.3 to 1.9 amplitude solves per effective sample vs 3.2 to 42 for Gibbs + MH; NUTS inside the Gibbs loop does not help | F5 |
| E002 | 2026-10-25 | Does Commander4's "re-solve then accept on chi^2" move target the right posterior? | `prototypes/specidx-bench/c4test.jl`, `c4global_check.jl` | No. 1.24x too wide per pixel; 3.9x (256 px) and 5.7x (1024 px) too wide for a full-sky scalar index, means shifted to the ridge posterior | F1 |
| E003 | 2026-10-25 | Independent check of the exact full-sky reference used in E002 | `prototypes/specidx-bench/global_gibbs_check.jl` | Commander3-style Gibbs (no shared math) matches the collapsed exact posterior: mean shifts <= 0.06 sd, sd ratio 0.95 to 1.00 | F1 |
| E004 | 2026-10-25 | Is language a bottleneck for the per-pixel marginal likelihood? | `prototypes/langbench/` | nside 512, 8 bands, 2 cores: Julia 0.070 s, Numba 0.060 s, C++ comparable, NumPy 1.8 s | F6 |
| E005 | 2026-10-25 | Do the new `likelihood` modes reproduce the exact target inside Commander4's MCMC class? | `tests/test_mcmc_likelihood.py` | marginal, ridge, conditional match exact; legacy `resolve_chisq` 3.2x / 4.4x too wide (150 / 600 px), mean +0.50 / +0.83 sd | F1, F2 |
| E006 | 2026-10-25 | End to end: full Commander4 (7 MPI ranks, maps only, nside 16, 7 bands) vs exact grid target | `aux/validate_mcmc_likelihood.py` | Ridge within 0.03 sd of exact, width ratios 0.94 / 0.97. Flat-prior marginal: exact target on the hard bounds, chain stuck in a local mode (barrier ~22). Ridge target 1.9 sd / 3.7 sd from input because smoothed noise is treated as white | F2, F3, F8 |
| E007 | 2026-10-25 | Why does an unsmoothed maps-only run give reduced chi^2 ~794? | first attempt of E006 | Per-pixel amplitudes are stored as alms with lmax = 3 nside - 1; the round trip alone raises reduced chi^2 from 0.57 to ~794 | F7 |
| E008 | 2026-10-25 | Are the region kernel and NUTS sampler correct? | `tests/test_nuts.py`, kernel checks in `tests/test_mcmc_likelihood.py` | Region kernel matches NumPy to 1e-16 and finite-difference gradients to ~1e-6; NUTS passes correlated Gaussian, 20-D scales, log-gamma, boundary, and 12-region joint NUTS vs exact | foundation for Q3 |
| E009 | 2026-10-25 | Does Commander3's nonlinear sampler share F1? Are its per-pixel lnL options exact? | `E009_c3_nonlin_review.md` | BP-era alm/region and local chisq samplers valid (amplitudes fixed). 'marginal' lnL has the log-det sign flipped (used for AME nu_p in BP final). 2024 MH module has the F1 structure | F10 |
| E010 | 2026-10-25 | Are Commander3's polarized index (alm/pixel-region) and bandpass samplers valid as used in BeyondPlanck? | `E010_c3_polarized_and_bandpass_review.md` | Yes: amplitudes, n_corr and gains held fixed; one symmetric proposal per iteration for bandpass. Approximations: prior-drawn frozen regions, diagonal spurious-map chi^2 for relative shifts | F10 |
| E011 | 2026-10-25 | Does nextgen build and pass on the Mac mini, and do E005/E006 reproduce there? | `E011_mac_mini_repro.md` | Yes. 859 passed / 4 skipped; E005 ratios 3.22x / 4.35x and shifts +0.50 / +0.83 sd identical; E006 ridge +0.03 sd, ratios 0.94 / 0.97, marginal stuck in the -3.5 mode, exact targets identical. ~80 s per mode on 7 ranks | F1, F2, F3, F8 confirmed |
| E012 | 2026-10-25 | Exact target for the hierarchical index-field prototype (M2.1): collapsed per-pixel likelihood, non-centred Matern/GMRF field, hyperpriors; L1 checks and SBC | `E012_gmrf_index_field_target.md` | running: L1 passed; L2 SBC default configuration PASS (300 replicates, 12 quantities); coverage study (gauss / w0 / w1) and centred SBC running | D008, D009 |
| E013 | 2026-10-25 | Q14: effect of the sign-flipped `marginal` volume term on BeyondPlanck's AME nu_p, exact per-pixel grid on BP data | `E013_bp_ame_nup_volume_term.md` | Negligible: BP nu_p was one full-sky parameter at high S/N; sign flip shifts it by 0.8 of a 0.01 GHz conditional sd (0.008 GHz); per pixel +0.04 ridge sd median, +0.14 at S/N 3 to 10 | F10 |
| E014 | 2026-10-25 | Q1: legacy resolve_chisq vs ridge/conditional/marginal on real DR1+PR4+Haslam maps at nside 64 | `E014_q1_real_maps_lowres.md` | running: pass 2 (NPIPE 353/857, PR3 100, dx12 mask); ridge matches exact grid within 0.2 sd; resolve_chisq 16 to 17x too wide, beta_d shifted 3 sd; conditional and marginal running | F13 |
| E015 | 2026-10-25 | Is BP's absolute bandpass chi-square reduced over MPI ranks? | `E015_bp_abs_bandpass_chisq_reduce.md` | No: accepted on rank 0's scans only, at BP10 and HEAD, also WMAP/DIRBE/SPIDER/LiteBIRD/QUIET modules | F11 |

## Notes on E006 settings
NSIDE 16, common FWHM 480 arcmin, bands 30/44/70/100/143/217/353 GHz, noise RMS x20, dust amplitudes
`10**(0.5 + 1.5 t)`; start beta_s = -2.9, beta_d = 1.595; proposals 0.08 / 0.06; bounds beta_s in
[-4, -2], beta_d in [1.0, 2.5]; 250 iterations, 80 burn-in, 50 MH steps per iteration. About 15 min
per mode on 2 cores. The exact target is a bounded coarse-to-fine grid computed for both log-det
settings by the analysis script.
