# CLAUDE.md: Commander4 next-generation research fork

## Mission

This is Trygve Leithe Svalheim's fork (GitHub `trygvels`) of Cosmoglobe/Commander4, the Bayesian
end-to-end CMB Gibbs sampler. Trygve has a PhD in CMB component separation from UiO and was a
BeyondPlanck author (synchrotron polarization, spectral-index sampling).

Goal: proper end-to-end CMB analysis with robust error propagation and explicit control of
degeneracies and biases, and breaking those degeneracies by better methods, better component or
instrument models, or new data. The bar is publication-level results that go beyond BeyondPlanck
and Cosmoglobe DR1/DR2, whose conclusions we are free to disagree with. Work continues until results
are publication-worthy. Negative results are recorded, never hidden.

## Read these first

| File | What it holds |
|---|---|
| `nextgen/STATUS.md` | Latest status, what is running, what needs Trygve |
| `nextgen/FINDINGS.md` | What has been established, with numbers and experiment IDs |
| `nextgen/OPEN_QUESTIONS.md` | What must be found out next, ranked, with how to answer each |
| `nextgen/ROADMAP.md` | Phases, milestones, target papers, acceptance criteria |
| `nextgen/DECISIONS.md` | Decision log with rationale (append only) |
| `nextgen/VALIDATION.md` | The validation ladder every method must pass |
| `nextgen/COMPUTE.md` | Compute tiers, machine requirements, where to run what |
| `nextgen/DATA.md` | Public datasets, URLs, sizes, how to fetch |
| `nextgen/LITERATURE.md` | Curated literature with arXiv IDs and one-line takeaways |
| `nextgen/experiments/registry.md` | Every experiment run so far |
| `nextgen/research/` | Raw research reports from 2026-10-25 (agent summaries; verify before citing) |

## Branches

- `main`: tracks upstream Cosmoglobe/Commander4. Never commit research work here.
- `fix/mcmc-marginal-acceptance`: PR-ready correctness fix for spectral-index MCMC (likelihood
  options marginal/ridge/conditional/resolve_chisq; `auto` resolves to ridge or conditional).
- `feat/spatially-varying-beta-foundation`: per-region likelihood kernel with analytic gradients,
  and a NUTS sampler. Stacked on the fix.
- `nextgen`: integration branch, everything above plus this research scaffolding. Work here, or on
  short-lived feature branches cut from here and merged back.

## Working principles (non-negotiable)

1. **Correctness before speed.** Every sampler change states its target distribution explicitly
   (including priors) and is tested against an exact or independently computed target (grid,
   analytic, or an independent valid sampler) before any efficiency claim.
2. **Measure efficiency fairly.** Effective samples per unit of the dominant cost (amplitude
   solves, SHTs, gradient evaluations, CPU-hours). Never compare raw wall time of different
   implementations as if it were algorithmic.
3. **Simulations with known truth first.** Real data only after the method passes the validation
   ladder in `nextgen/VALIDATION.md`.
4. **Separate "exact for a stated prior" from "unbiased against truth".** Prior choices (flat
   amplitude priors, Jeffreys, ridge/profile) change answers. When they differ, report both.
5. **Log every experiment** in `nextgen/experiments/` (question, setup, command, result,
   conclusion), including failures. Claims in docs cite an experiment ID or a paper.
6. **Independent review for statistical changes.** Before merging a change to a sampler,
   likelihood or prior, run the `stats-reviewer` subagent, which re-derives the target without
   seeing your reasoning.
7. **Keep upstream compatibility where cheap.** Behaviour changes get a config switch and a
   documented default. Never leave the existing test suite red.
8. **Do not overstate.** If a result could be a bug, it is a bug until shown otherwise.

## Ask Trygve before

- changing a default that alters scientific results;
- opening pull requests or pushing anything to upstream Cosmoglobe repositories;
- contacting Cosmoglobe developers, or publishing anything;
- spending large cluster allocations, or downloads over ~200 GB;
- deleting data or rewriting published branch history.

Write a short update to `nextgen/STATUS.md` at every milestone, and at least once per long
autonomous stretch: what was learned (numbers), what is running, what is next, what needs a
decision. Trygve's writing preference: never use em-dashes in text written for him.

## Build and test

Linux (verified on 2026-10-25, Python 3.13, MPICH from pip):

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install mpi4py mpich
git submodule update --init --recursive        # ducc0
pip install -e . --group test                  # rerun after any C++ change
python -m pytest -q -n auto -p no:cacheprovider   # ~860 tests, ~1.5 min on 2 cores
```

macOS (Apple Silicon, not yet verified): `nextgen/scripts/setup_macos.sh`. Record any fixes needed
in `nextgen/COMPUTE.md`.

Focused tests: `tests/test_mcmc_likelihood.py`, `tests/test_nuts.py`.

End-to-end maps-only validation (7 MPI ranks, about 15 min per likelihood mode on 2 cores):

```bash
python aux/validate_mcmc_likelihood.py make --dir /tmp/c4val
mpirun -n 7 commander4 -p /tmp/c4val/param_ridge.yml
python aux/validate_mcmc_likelihood.py analyze --dir /tmp/c4val
```

Fast prototypes (Julia, outside Commander4): `nextgen/prototypes/specidx-bench/`.

## Code map

| Path | Role |
|---|---|
| `src/commander4/cli.py` | Entry point; two chains pipelined across TOD and CompSep ranks |
| `src/commander4/compsep/processing.py` | Sampling groups (CG, per-pixel, MCMC), configs, likelihood resolution |
| `src/commander4/compsep/mcmc.py` | MH loop and likelihood modes; legacy path `_run_resolve_chisq` |
| `src/commander4/compsep/spectral_index.py` | Scalar spectral-index groups |
| `src/commander4/compsep/perpix_solver.py` | Per-pixel amplitude solver; `PerPixelMarginalLikelihood` |
| `src/commander4/compsep/cg_solver.py`, `preconditioners.py` | CG constrained realizations |
| `src/commander4/sky/diffuse_components.py` | SEDs (`get_sed`), `spatially_varying_MM` hooks |
| `src/lib_cpp/compsep_pymod.cc` | C++ kernels: `solve_perpix`, `marginal_loglike_perpix`, `marginal_loglike_regions` |
| `src/commander4/math_utils/nuts.py` | NUTS with dual averaging and windowed mass adaptation |
| `src/commander4/tod/` | Gain, noise, mapmaking, sidelobes, jumps, glitches |
| `params/param_explanation.yml` | Every parameter-file key, documented |
| `aux/validate_mcmc_likelihood.py` | Maps-only end-to-end validation harness |

Internal units are uK_RJ. Existing code uses a 100-character line limit and docstrings that state
the math; follow that.

## Subagents and commands in this repo

- `.claude/agents/stats-reviewer.md`: independent check of a sampler, likelihood or prior change.
- `.claude/agents/lit-scout.md`: literature search with strict citation discipline.
- `.claude/agents/validation-runner.md`: runs the validation ladder for a branch and reports.
- `/new-experiment <question>`: scaffold an experiment record.
- `/validate`: run tests plus the relevant validation level and summarize.
- `/status`: refresh `nextgen/STATUS.md` from the registry and recent commits.
