#!/usr/bin/env bash
# Set up a Commander4 development environment on macOS (Apple Silicon).
# Not yet verified on a real Mac: record every deviation in nextgen/COMPUTE.md.
#
# Usage, from the repository root on branch nextgen:
#   bash nextgen/scripts/setup_macos.sh            # core build + tests
#   WITH_JULIA=1 bash nextgen/scripts/setup_macos.sh   # also the Julia prototypes
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

echo "== Machine record (copy into nextgen/COMPUTE.md)"
sysctl -n machdep.cpu.brand_string || true
echo "cores: $(sysctl -n hw.ncpu)  perf cores: $(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo ?)"
echo "memory GB: $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"
df -h "$ROOT" | tail -1

echo "== Toolchain"
xcode-select -p >/dev/null 2>&1 || { echo "Run: xcode-select --install"; exit 1; }
command -v brew >/dev/null || { echo "Install Homebrew first: https://brew.sh"; exit 1; }
# gcc provides gfortran for any dependency that builds Fortran from source (camb, pixell).
brew install cmake ninja python@3.13 gcc tmux

PY="$(brew --prefix python@3.13)/bin/python3.13"
"$PY" -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip

echo "== MPI (pip-installed MPICH, same as the verified Linux setup)"
pip install mpich mpi4py
python -c "from mpi4py import MPI; print('MPI ok:', MPI.Get_library_version().splitlines()[0])"

echo "== Commander4"
git submodule update --init --recursive
export FC="$(brew --prefix)/bin/gfortran"
pip install -e . --group test

echo "== Tests"
python -m pytest -q -n auto -p no:cacheprovider tests/test_mcmc_likelihood.py tests/test_nuts.py
echo "Focused tests passed. Full suite: python -m pytest -q -n auto -p no:cacheprovider"

if [[ "${WITH_JULIA:-0}" == "1" ]]; then
  echo "== Julia prototypes"
  command -v juliaup >/dev/null || brew install juliaup
  juliaup add release && juliaup default release
  (cd nextgen/prototypes/specidx-bench && julia --project=. -e 'using Pkg; Pkg.instantiate()')
fi

cat <<'EOF'

Done. Next:
  source .venv/bin/activate
  caffeinate -dimsu &            # keep the Mac awake during long runs
  tmux new -s c4                 # run Claude Code inside tmux
  claude                         # then paste nextgen/KICKOFF.md
EOF
