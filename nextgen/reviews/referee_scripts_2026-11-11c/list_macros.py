"""List every \\N{key} used in the paper with the value defined in gen_numbers.tex (referee read, draft 4)."""
import os, re, glob
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "papers", "paperA1")
defs = {}
for line in open(os.path.join(P, "gen_numbers.tex")):
    m = re.match(r"\\expandafter\\def\\csname num@(.+?)\\endcsname\{(.*)\}\s*$", line)
    if m:
        defs[m.group(1)] = m.group(2)
used = {}
for f in sorted(glob.glob(os.path.join(P, "*.tex"))):
    if os.path.basename(f).startswith("gen_"):
        continue
    for i, line in enumerate(open(f), 1):
        for k in re.findall(r"\\N\{([^}]+)\}", line):
            used.setdefault(k, []).append(f"{os.path.basename(f)}:{i}")
for k in sorted(used):
    print(f"{k:32s} = {defs.get(k, '!!! UNDEFINED'):40s} {' '.join(used[k])}")
print(len(used), "keys used;", len(defs), "defined")
