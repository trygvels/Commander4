"""Referee check (draft 3): list keys and shapes of the result files used by the paper. Read only."""
import numpy as np, glob, os, sys
pats = sys.argv[1:]
for pat in pats:
    for f in sorted(glob.glob(pat)):
        z = np.load(f, allow_pickle=True)
        print("==", f)
        for k in z.files:
            a = z[k]
            s = "" 
            if a.size <= 12:
                s = np.array2string(np.asarray(a), precision=6, max_line_width=200)
            print(f"   {k:22s} {str(a.shape):14s} {a.dtype} {s}")
