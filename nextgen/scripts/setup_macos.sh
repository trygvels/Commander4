#!/usr/bin/env bash
# Set up a Commander4 development environment on macOS (Apple Silicon).
# Verified 2026-10-25 on an Apple M4 Mac mini (macOS 27.0); see nextgen/COMPUTE.md.
#
# Usage, from the repository root on branch nextgen:
#   bash nextgen/scripts/setup_macos.sh                   # core build + focused tests
#   WITH_JULIA=1 bash nextgen/scripts/setup_macos.sh      # also the Julia prototypes
#   C4_BREW_INSTALL=1 bash nextgen/scripts/setup_macos.sh # let the script brew-install missing tools
#
# Nothing here needs Homebrew for the build itself: every Python dependency (including cmake,
# ninja, mpich, mpi4py, pixell, qpoint, camb, ducc0) ships an arm64 wheel for CPython 3.13. The
# script uses any python3.13 on PATH (uv, Homebrew, python.org), installs cmake and ninja into the
# venv, and only reports missing Homebrew tools unless C4_BREW_INSTALL=1 is set.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

echo "== Machine record (copy into nextgen/COMPUTE.md)"
sysctl -n machdep.cpu.brand_string || true
echo "cores: $(sysctl -n hw.ncpu)  perf cores: $(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo ?)  eff cores: $(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null || echo ?)"
echo "memory GB: $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"
echo "macOS: $(sw_vers -productVersion)"
df -h "$ROOT" | tail -1

echo "== Toolchain"
xcode-select -p >/dev/null 2>&1 || { echo "Run: xcode-select --install"; exit 1; }

# Optional Homebrew tools. tmux keeps the long Claude Code session alive; gfortran is only needed
# if some dependency has no arm64 wheel and must build Fortran from source (none did on 2026-10-25).
missing=()
for tool in tmux; do command -v "$tool" >/dev/null || missing+=("$tool"); done
if (( ${#missing[@]} )); then
  if [[ "${C4_BREW_INSTALL:-0}" == "1" ]]; then
    command -v brew >/dev/null || { echo "Install Homebrew first: https://brew.sh"; exit 1; }
    brew install "${missing[@]}"
  else
    echo "Not installed (not needed for the build): ${missing[*]}"
    echo "Install with: brew install ${missing[*]}   (or rerun with C4_BREW_INSTALL=1)"
  fi
fi

# Python 3.13: any interpreter on PATH, then Homebrew's, then uv-managed.
PY=""
for cand in python3.13 "$(brew --prefix python@3.13 2>/dev/null || true)/bin/python3.13"; do
  [[ -n "$cand" && -x "$(command -v "$cand" 2>/dev/null || true)" ]] && { PY="$(command -v "$cand")"; break; }
done
if [[ -z "$PY" ]] && command -v uv >/dev/null; then
  uv python install 3.13
  PY="$(uv python find 3.13)"
fi
if [[ -z "$PY" ]]; then
  if [[ "${C4_BREW_INSTALL:-0}" == "1" ]]; then
    brew install python@3.13
    PY="$(brew --prefix python@3.13)/bin/python3.13"
  else
    echo "No python3.13 found. Install one (uv python install 3.13, or brew install python@3.13)."
    exit 1
  fi
fi
echo "python: $PY ($("$PY" --version))"

if [[ ! -x .venv/bin/python ]]; then
  "$PY" -m venv .venv
fi
source .venv/bin/activate
python -m pip install --upgrade pip

echo "== Build tools (pip wheels; scikit-build-core uses them from PATH)"
pip install cmake ninja

echo "== MPI (pip-installed MPICH, same as the verified Linux setup)"
pip install mpich mpi4py
python -c "from mpi4py import MPI; print('MPI ok:', MPI.Get_library_version().splitlines()[0])"

echo "== Commander4"
git submodule update --init --recursive
if command -v gfortran >/dev/null; then export FC="$(command -v gfortran)"; fi
pip install -e . --group test

echo "== Tests"
python -m pytest -q -n auto -p no:cacheprovider tests/test_mcmc_likelihood.py tests/test_nuts.py
echo "Focused tests passed. Full suite: python -m pytest -q -n auto -p no:cacheprovider"

if [[ "${WITH_JULIA:-0}" == "1" ]]; then
  echo "== Julia prototypes"
  if ! command -v juliaup >/dev/null; then
    command -v brew >/dev/null || { echo "Install Homebrew first: https://brew.sh"; exit 1; }
    brew install juliaup
  fi
  juliaup add release && juliaup default release
  (cd nextgen/prototypes/specidx-bench && julia --project=. -e 'using Pkg; Pkg.instantiate()')
fi

cat <<'MSG'

Done. Next:
  source .venv/bin/activate
  caffeinate -dimsu &            # keep the Mac awake during long runs
  tmux new -s c4                 # run Claude Code inside tmux
  claude                         # then paste nextgen/KICKOFF.md
MSG
