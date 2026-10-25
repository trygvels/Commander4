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
