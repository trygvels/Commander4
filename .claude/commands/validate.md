---
description: Run tests and the relevant validation level for the current change, then summarize
argument-hint: [level, e.g. L1 or L3; default: decide from the diff]
---

Validate the current branch.

1. Summarize the change: `git diff --stat main...HEAD` and `git log --oneline main..HEAD | head -20`.
2. Decide the validation levels from `nextgen/VALIDATION.md` (or use $ARGUMENTS if given).
3. Delegate the runs to the `validation-runner` subagent with the change summary, the levels, and the
   experiment record to write into (create one with /new-experiment if none exists).
4. If the change touches a sampler, likelihood or prior, also run the `stats-reviewer` subagent on
   the diff, giving it the stated target but not your reasoning.
5. Report in a few lines: pass or fail per level with key numbers, reviewer verdict, and what must
   happen before merging.
