# Starting the long-running Claude Code session

## On the Mac mini (once)

```bash
git clone -b nextgen https://github.com/trygvels/commander4.git
cd commander4
git remote add upstream https://github.com/Cosmoglobe/Commander4.git
bash nextgen/scripts/setup_macos.sh      # add WITH_JULIA=1 for the Julia prototypes
source .venv/bin/activate
caffeinate -dimsu &                      # keep the machine awake
tmux new -s c4                           # so the session survives a closed terminal
claude                                   # start Claude Code in the repository root
```

Claude Code reads `CLAUDE.md` automatically, and finds the subagents in `.claude/agents/` and the
commands `/new-experiment`, `/validate`, `/status` in `.claude/commands/`.

For long unattended stretches, pre-approve the routine commands yourself with `/permissions` (for
example `python -m pytest`, `mpirun`, `commander4`, `julia`, and local `git` commands other than
push), and keep pushes, deletions and anything touching `upstream` on ask.

## Prompt to paste

```text
You are taking over a long-running research project on this repository (my fork of Cosmoglobe's
Commander4). Read CLAUDE.md, then nextgen/STATUS.md, FINDINGS.md, OPEN_QUESTIONS.md, ROADMAP.md,
DECISIONS.md and VALIDATION.md before doing anything else.

Goal: end-to-end CMB analysis with robust error propagation and explicit control of degeneracies and
biases, and publication-level results beyond BeyondPlanck and Cosmoglobe DR1/DR2. The first target
is spatially varying spectral parameters at the finest resolution the signal-to-noise supports,
done in a statistically optimal way (ROADMAP Phase 2, OPEN_QUESTIONS Q2 to Q4), leading to the
methods paper (Paper A) and then the synchrotron science paper (Paper B).

Work plan:
1. Phase 1 first: M1.1 and M1.2 (build, full test suite, reproduce E005 and E006 on this machine;
   fill in the machine record in nextgen/COMPUTE.md and adjust the compute tiers to it).
2. Then run M1.3/M1.4 (Q1 on real low-resolution maps) and M2.1 (hierarchical GMRF index-field
   prototype) in parallel. Before writing M2.1 code, write the exact target (likelihood, priors,
   hyperpriors, parameterization) into its experiment record and have the stats-reviewer subagent
   check it.
3. Continue through the roadmap. Re-rank OPEN_QUESTIONS when findings change the picture, and log
   any change of plan in DECISIONS.md.

Rules: follow the working principles in CLAUDE.md strictly. Every claim needs an experiment ID or a
paper. Log every experiment, including failures, with /new-experiment. Use /validate before merging
any method change. Run long jobs with nohup and poll them; never let a run die with the shell.
Commit small, logical steps on nextgen (or short feature branches merged back) and push to origin
regularly. Update nextgen/STATUS.md with /status at every milestone and at least every few hours of
autonomous work. Use subagents for independent work in parallel (literature, validation, review),
but keep the core reasoning yourself. Stop and ask me only for the items under "Ask Trygve before"
in CLAUDE.md; otherwise make the reasonable call, record it, and keep going. Never use em-dashes in
anything you write for me.
```

## Checking in

- `nextgen/STATUS.md` is the summary to read first; `git log --oneline` shows the work.
- To redirect: tell the session in plain words, or edit OPEN_QUESTIONS/ROADMAP and say so.
- For cluster work (T3/T4 in COMPUTE.md), give the session SSH access to the cluster and say which
  allocation it may use, and how many CPU-hours.
