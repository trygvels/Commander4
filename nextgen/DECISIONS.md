# Decision log

Append new decisions at the bottom. Never rewrite an old entry; supersede it with a new one that
references it.

## D001 (2026-10-25): build on Commander4, in Python and C++
Commander4 already provides the engineering gains a rewrite would give (Python front end, C++
kernels via nanobind, ducc0, MPI, YAML config, ~860 tests) and is actively developed. Kernel speed is
at parity between Julia, Numba and C++ (F6). A separate Julia rewrite would duplicate an active
official code. Julia stays for fast prototypes only.

## D002 (2026-10-25): spectral-index MCMC must use an exact acceptance
The original resolve-then-chi^2 acceptance is not a valid MH step (F1). New `likelihood` option on
MCMC groups: marginal, ridge, conditional, resolve_chisq (legacy, warns).

## D003 (2026-10-25): `auto` resolves to ridge where available, else conditional
The flat-amplitude-prior marginal can be dominated by its volume term for full-sky indices (F3).
Ridge is the standard spectral likelihood and corresponds to a joint Jeffreys prior (F4). Revisit
after Q2 (coverage study).

## D004 (2026-10-25): spatially varying indices aim for the finest resolution the S/N supports
Trygve's direction. Hard region boundaries are unphysical and alms ring, but the goal is the
statistically optimal model. Leading candidate: hierarchical GMRF index fields with inferred
hyperparameters (Q3). Alms, clustering and fixed ladders are comparison baselines.

## D005 (2026-10-25): research scaffolding lives in `nextgen/`, work happens on branch `nextgen`
`main` keeps tracking upstream. Upstream pull requests are cut from focused branches, and only
after Trygve decides (Q13).

## D006 (2026-10-25): allow macOS development builds on `nextgen`
Upstream `CMakeLists.txt` stops on any non-Linux system. The long-running session will work on a
Mac mini, so `nextgen` allows Darwin with a warning, picks `-march=native` or `-mcpu=native` by
compiler check, and skips `-s` on Apple. Linux behaviour is unchanged (re-verified). Keep this out of
upstream pull requests unless Cosmoglobe wants it.

## D007 (2026-10-25): long-running research happens on the Mac mini, heavy runs on a cluster
The cloud workspace (2 vCPU, 7.8 GB, jobs stopped after each reply) cannot hold research runs. The
Mac handles T0 to T2 (COMPUTE.md); T3/T4 go to ITA or Sigma2 via Slurm.

## D008 (2026-10-25): sampler correctness is tested with a proper Gaussian amplitude prior; prior choices are a coverage study
Neither the flat amplitude prior (marginal, w = 1) nor the profile likelihood (ridge, w = 0) has a
proper generative distribution, so an SBC run with them cannot distinguish a sampler error from a
prior-volume effect (E012 review). The index-field prototype and every later sampler therefore pass
L1/L2 with the exact marginal under a proper Gaussian amplitude prior N(m, T) (a two-line change to
the collapsed kernel), and ridge / marginal / Jeffreys are compared in a separate frequentist
coverage study reported as bias and coverage versus S/N (principle 4). D003 (ridge default) stands.

## D009 (2026-10-25): index-field prior conventions
The hierarchical field drops the constant mode (the field mean is mu, sigma is the sd of
fluctuations about it), uses radians for the correlation length in both the GMRF and the alm family
with the graph-Laplacian eigenvalues rescaled by 1/(3 h^2), and samples the bounded log-uniform
correlation length through a logit transform with its Jacobian. Reasons in E012 v2.

## D010 (2026-10-26): Q1 is answered at nside 64; Q2 (prior choice) moves to the top of the ranking
E014 settled Q1 on real maps: the legacy move is 16 to 17x too wide and shifts beta_d by 3 sd, and the
exact `ridge` and `marginal` chains both reproduce independent grid targets. The same run showed that
the flat-prior volume term moves beta_s by 54 posterior sd with no change in width or chi-square, so
the prior is the dominant systematic for full-sky indices on real data, far above sampler width. The
ranking in OPEN_QUESTIONS therefore becomes Q2 (coverage study, E012 coverage run in progress, then
a scalar-index coverage table at M2.2), Q3/Q4 (index-field model and sampler, M2.x), then the second
half of Q14. Q1's remaining items (IQU, nside 128 scaling, HFI full set) are follow-ups in E014, not
blockers. The `conditional` Gibbs mode is dropped from real-data comparisons (F14: conditional sd
2e-6 against marginal 3.6e-4; of order 3e4 iterations per independent sample); it stays in the code
as a documented option. No default changes.
