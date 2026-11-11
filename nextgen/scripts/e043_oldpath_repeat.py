"""E043: repeated Cholesky factorizations on the old path (C-ordered array, lower=True, that is LAPACK uplo = U) at
one matrix order, to be run under the malloc canary. Each repetition uses a new positive definite matrix
(rank 3 plus a random diagonal) and factorizes in place, as `RegionLikelihood.solve` did before the fix.

    DYLD_INSERT_LIBRARIES=~/data/c4_tools/libcanary.dylib python -I nextgen/scripts/e043_oldpath_repeat.py 18432 20 2> err.txt
"""
import sys
import time

import numpy as np
import scipy.linalg as sl

n, reps = int(sys.argv[1]), int(sys.argv[2])
rng = np.random.default_rng(43)
t0 = time.time()
for r in range(reps):
    v = rng.standard_normal((n, 3))
    A = v @ v.T
    A[np.diag_indices(n)] += 5.0 + rng.random(n)
    c = sl.cho_factor(A, lower=True, overwrite_a=True, check_finite=False)
    x = sl.cho_solve(c, np.ones(n), check_finite=False)
    print(f"order {n} repetition {r}: logdet {2 * np.log(np.diag(c[0])).sum():.4f}, sum of the solution {x.sum():.8f} [{time.time() - t0:.0f} s]", flush=True)
    del A, c, v
