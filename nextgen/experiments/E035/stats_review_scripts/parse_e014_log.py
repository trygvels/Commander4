"""Read-only parse of the E014 legacy (resolve_chisq) run log. Extract every MH step:
chain, iteration, step, dloglike, accept, current (dust, sync), proposed (dust, sync). Save npz."""
import re, sys
import numpy as np
log = sys.argv[1]
out = sys.argv[2]
step_re = re.compile(r"MCMC step (\d+)/(\d+): dloglike=([-+0-9.eE]+|nan|-?inf) \[(ACCEPT|reject)\] "
                     r"\(\{'dust': ([-0-9.eE]+), 'sync': ([-0-9.eE]+)\} -> \{'dust': ([-0-9.eE]+), 'sync': ([-0-9.eE]+)\}\)")
done_re = re.compile(r"Chain (\d+), iteration (\d+) complete")
rows, block = [], []
nother = 0
with open(log) as f:
    for line in f:
        if "MCMC step" in line:
            m = step_re.search(line)
            if m is None:
                nother += 1
                continue
            block.append((int(m.group(1)), float(m.group(3)), m.group(4) == "ACCEPT",
                          float(m.group(5)), float(m.group(6)), float(m.group(7)), float(m.group(8))))
        elif "iteration" in line and "complete" in line:
            m = done_re.search(line)
            if m:
                c, it = int(m.group(1)), int(m.group(2))
                for b in block:
                    rows.append((c, it) + b)
                block = []
a = np.array(rows, dtype=float)
print("rows", a.shape, "unparsed step lines", nother, "leftover block", len(block))
print("chains", np.unique(a[:, 0]), "iterations", a[:, 1].min(), a[:, 1].max())
np.savez(out, a=a)
