# Compute: where to run what

## Measured: Claude cloud workspace (2026-10-25)

| Resource | Value |
|---|---|
| CPU | 2 vCPU, Intel Xeon 2.8 GHz |
| Memory | 7.8 GB |
| Free disk | ~26 GB |
| GPU | none |
| Network | package registries, GitHub, some data hosts (julialang, conda) via allowlist |
| Background jobs | stopped within minutes after a reply ends |

What fit here: the full test suite (~1.5 min), toy benchmarks, and one nside-16 maps-only validation
run on 7 oversubscribed MPI ranks (~15 min per 250 Gibbs iterations). Nothing larger. Long MCMC runs
cannot outlive a reply, so this workspace is not suitable for the research phase.

## Requirements by analysis tier (estimates)

| Tier | Example | Memory | CPU | Where |
|---|---|---|---|---|
| T0 prototypes, unit and statistical tests | toy samplers, pytest, SBC at nside 8 to 16 | < 4 GB | 2 to 8 cores | anywhere |
| T1 maps-only, low resolution | nside 32 to 128, 10 to 30 bands, hundreds of SBC replicates | 1 to 4 GB per run | 8 to 12 cores, hours | Mac mini |
| T2 maps-only, science resolution | nside 512, 20 to 35 bands IQU, per-pixel or CG with lmax ~1000 | ~6 to 12 GB (per-pixel), more for CG | 10+ cores, hours to days | Mac mini (24 GB: per-pixel at nside 512 only, with headroom) or cluster for CG |
| T3 maps-only, high resolution | nside 1024 to 2048, full CG | 30 to 100+ GB | many nodes | cluster |
| T4 TOD-level | LFI, WMAP, HFI time-ordered data | hundreds of GB to TB (DR1 chains alone 841 GB each) | 170 to 3000 CPU-h per Gibbs sample (published) | cluster only |

Per-pixel memory at nside 512: one IQU float64 map is ~75 MB; 30 bands of sky plus inverse variance
is ~4.5 GB, and the per-pixel solver gathers all bands on one rank. nside 1024 is 4x that.

## Recommendation

1. **Run the long Claude Code session locally on the Mac mini** (measured below: M4, 10 cores,
   24 GB, ~100 GB free). That covers T0 and T1 comfortably and T2 per-pixel runs up to about
   nside 512 with 20 to 30 bands (the per-pixel solver gathers all bands on one rank, ~4.5 GB at
   30 bands; leave headroom for the OS and the Claude session). T2 CG runs with lmax ~1000 and
   anything at nside 1024 go to a cluster. Disk is the tighter limit: ~100 GB free means low
   resolution maps and simulations only; do not fetch full-resolution DR1 chains here. Keep the
   machine awake (`caffeinate -dimsu`) and run Claude Code inside `tmux` so long jobs survive
   disconnects.
2. **Use a cluster for T3/T4**: ITA owls or Sigma2 Olivia (CPU nodes: 256 cores, 768 GiB), via the
   Cosmoglobe group's allocation if Trygve still has access. Claude Code on the Mac can prepare and
   submit Slurm jobs over SSH, and pull results back.
3. **Use cloud Claude Code sessions only for code-only tasks** (refactors, reviews, writing): the
   measured workspace is far too small for the research runs, and cannot keep jobs alive.

## macOS build notes (verified 2026-10-25)

Build, full test suite and the E005/E006 reproductions all passed on the Mac mini below (E011).
Script: `nextgen/scripts/setup_macos.sh`. What the script does and what had to change from the
unverified draft:

- **No Homebrew needed for the build.** Every Python dependency has an arm64 CPython 3.13 wheel on
  PyPI (checked with `pip download --only-binary=:all:` for mpich 5.0.2, mpi4py 4.1.2, pixell,
  qpoint, camb, ducc0, healpy, numba, h5py, numpy-quaternion, cmake, ninja). The draft script ran
  `brew install cmake ninja python@3.13 gcc tmux` unconditionally; none of those were installed on
  the machine and none were needed. The script now installs `cmake` and `ninja` as pip wheels in the
  venv (scikit-build-core picks them up from PATH), uses any `python3.13` on PATH (here the
  uv-managed CPython 3.13.11; Homebrew's `python3` is 3.14), and only reports missing Homebrew tools
  unless `C4_BREW_INSTALL=1` is set.
- **gfortran is not needed**: pixell, qpoint and camb install from wheels. The script sets `FC` only
  if a gfortran exists (the draft pointed `FC` at a Homebrew gfortran that was not there, which
  would have failed under `set -u`).
- **MPI**: pip `mpich` 5.0.2 plus `mpi4py` 4.1.2 work as on Linux. 7 ranks on 10 cores run
  unoversubscribed.
- **Compiler**: Apple clang 17.0.0 (Xcode command line tools, macOS 27.0). The `CMakeLists.txt`
  change on `nextgen` is what makes this build possible: upstream refuses non-Linux. On this
  toolchain `-march=native` is accepted (no fallback to `-mcpu=native` needed), the `-s` strip flag
  is skipped on Apple, and the macOS warning is printed once at configure time. Build time for the
  nanobind extension plus ducc0 sources: under a minute.
- **tmux and juliaup** were installed with Homebrew after Trygve approved it (tmux 3.8; Julia
  release via juliaup, and `nextgen/prototypes/specidx-bench` instantiated: 107 packages). They are
  for the session and the prototypes, not for the build.
- `pip install -e . --group test` needs pip >= 25.1; the script upgrades pip in the venv first.
- Test timings: focused tests (`test_mcmc_likelihood.py`, `test_nuts.py`) 24 passed in 35 s on the
  first run (numba compilation); full suite 859 passed, 4 skipped in 25 s with `-n auto` (10
  workers).

## Mac mini machine record

| Field | Value |
|---|---|
| CPU | Apple M4 (4 performance + 6 efficiency cores) |
| Cores | 10 |
| Memory | 24 GB |
| Free disk | ~100 GB of 460 GB (76% used) at setup time |
| OS | macOS 27.0 (build 26A5388g), arm64 |
| Python | CPython 3.13.11 (uv-managed) in `.venv`; numpy 2.5.4, scipy 1.18.1, ducc0 0.42.0 |
| MPI | MPICH 5.0.2 (pip), mpi4py 4.1.2 |
| Compiler | Apple clang 17.0.0, `-march=native -O3 -ffast-math`, IPO on |
| Build verified | 2026-10-25, commit 4137c7f (nextgen), full suite 859 passed / 4 skipped |
