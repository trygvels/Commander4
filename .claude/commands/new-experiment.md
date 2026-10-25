---
description: Scaffold a new experiment record and registry row
argument-hint: <one-line question>
---

Create a new experiment for the question: $ARGUMENTS

1. Read `nextgen/experiments/registry.md` and take the next free ID (EXXX, zero padded).
2. Pick a short kebab-case slug from the question.
3. Copy `nextgen/experiments/EXPERIMENT_TEMPLATE.md` to `nextgen/experiments/EXXX_<slug>.md`, fill
   in the title, date, current branch and commit (`git rev-parse --short HEAD`), the question, and
   the matching Q number from `nextgen/OPEN_QUESTIONS.md` if there is one.
4. Write the "Target and priors" and "Reference" sections now, before running anything. If the
   target or the reference is unclear, say so in the file and stop to think it through.
5. Add a row to the registry table with Result "running" and Findings empty.
6. Print the path of the new record.
