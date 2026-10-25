---
description: Refresh nextgen/STATUS.md from the registry, recent commits and running jobs
---

Update `nextgen/STATUS.md` with a new entry at the top, dated today:

1. Read the current STATUS.md, `nextgen/experiments/registry.md`, and `git log --oneline -20`.
2. Check for running jobs (`ps`, Slurm `squeue -u $USER` if on a cluster) and note what is running.
3. Write the entry with four parts: Learned (numbers, experiment IDs), Running, Next (milestone IDs
   from ROADMAP.md), Needs Trygve (decisions, access, data). Keep it under 25 lines.
4. If any finding changed, make sure FINDINGS.md, OPEN_QUESTIONS.md and ROADMAP.md agree with it.
5. Commit with message "status: <date>" and push the current branch to `origin`.

Never use em-dashes in this file.
