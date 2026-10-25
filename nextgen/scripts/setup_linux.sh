#!/usr/bin/env bash
# Set up a Commander4 development environment on Linux (verified 2026-10-25, Python 3.13, x86_64).
# On a cluster, prefer the site's MPI module over pip MPICH: load it, then `pip install mpi4py`
# with MPICC pointing at the site compiler wrapper, and skip the `pip install mpich` line.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

echo "== Machine record"
lscpu | grep -E 'Model name|^CPU\(s\)' || true
free -g | head -2
df -h "$ROOT" | tail -1

PY="${PYTHON:-python3}"
"$PY" -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip

if [[ "${SITE_MPI:-0}" != "1" ]]; then
  pip install mpich
fi
pip install mpi4py
python -c "from mpi4py import MPI; print('MPI ok:', MPI.Get_library_version().splitlines()[0])"

git submodule update --init --recursive
pip install -e . --group test

python -m pytest -q -n auto -p no:cacheprovider tests/test_mcmc_likelihood.py tests/test_nuts.py
echo "Focused tests passed. Full suite: python -m pytest -q -n auto -p no:cacheprovider"
