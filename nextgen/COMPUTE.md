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
| T2 maps-only, science resolution | nside 512, 20 to 35 bands IQU, per-pixel or CG with lmax ~1000 | ~6 to 12 GB (per-pixel), more for CG | 10+ cores, hours to days | Mac mini (32 GB+) or cluster |
| T3 maps-only, high resolution | nside 1024 to 2048, full CG | 30 to 100+ GB | many nodes | cluster |
| T4 TOD-level | LFI, WMAP, HFI time-ordered data | hundreds of GB to TB (DR1 chains alone 841 GB each) | 170 to 3000 CPU-h per Gibbs sample (published) | cluster only |

Per-pixel memory at nside 512: one IQU float64 map is ~75 MB; 30 bands of sky plus inverse variance
is ~4.5 GB, and the per-pixel solver gathers all bands on one rank. nside 1024 is 4x that.

## Recommendation

1. **Run the long Claude Code session locally on the Mac mini**, provided it has at least 32 GB of
   memory and 10+ cores (an M4 Pro with 48 to 64 GB is comfortable for T2), and 500 GB+ free disk
   for maps and simulations. With 16 GB, stay at T0/T1 and push T2 to a cluster. Keep the machine
   awake (`caffeinate -dimsu`) and run Claude Code inside `tmux` so long jobs survive disconnects.
2. **Use a cluster for T3/T4**: ITA owls or Sigma2 Olivia (CPU nodes: 256 cores, 768 GiB), via the
   Cosmoglobe group's allocation if Trygve still has access. Claude Code on the Mac can prepare and
   submit Slurm jobs over SSH, and pull results back.
3. **Use cloud Claude Code sessions only for code-only tasks** (refactors, reviews, writing): the
   measured workspace is far too small for the research runs, and cannot keep jobs alive.

First thing on the Mac: record `sysctl -n machdep.cpu.brand_string hw.ncpu hw.memsize` and free disk
here, and adjust the tiers above.

## macOS build notes (to be verified)

- Upstream `CMakeLists.txt` refused any non-Linux build. On `nextgen` it now allows macOS with a
  warning, tries `-march=native` then `-mcpu=native`, and skips the `-s` strip flag on Apple. The
  Linux build was re-verified after the change; the macOS build is not yet verified.
- Toolchain: Xcode command line tools, Homebrew `cmake`, `ninja`, `python@3.13`, `gcc` (gfortran).
- MPI: pip `mpich` plus `mpi4py`, as on Linux. Fallback: Homebrew `open-mpi` with
  `MPICC=$(brew --prefix open-mpi)/bin/mpicc pip install --no-binary mpi4py mpi4py`.
- `pixell` and `qpoint` may need `gfortran` (`brew install gcc`) if no arm64 wheel exists.
- ducc0 and the nanobind extension build with Apple clang; ducc0 has its own threading, so no
  OpenMP is required.
- Script: `nextgen/scripts/setup_macos.sh`. Record every deviation here.

## Mac mini machine record

| Field | Value |
|---|---|
| CPU | (fill in) |
| Cores | (fill in) |
| Memory | (fill in) |
| Free disk | (fill in) |
| Build verified | (date, commit) |
