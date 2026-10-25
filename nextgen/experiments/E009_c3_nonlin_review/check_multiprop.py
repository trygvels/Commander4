"""Exact check of pi-invariance for the Commander3 sample_bp multi-proposal scheme.

Discrete state space, symmetric proposal matrix Q, target pi.  Three schemes:
  MH     : single proposal (ndelta=2), reference.
  ABS    : ndelta-1 proposals all drawn from Q(x1,.) (commander.f90:405-407),
           sequentially accepted against current with min(1, pi(y)/pi(cur)).
  REL    : proposals chained: x_k ~ Q(x_{k-1},.) (commander.f90:393-397, 419),
           same sequential acceptance.
Transition kernel built by exhaustive enumeration, then ||pi K - pi|| is reported.
"""
import itertools, numpy as np
rng = np.random.default_rng(1)
n = 4
pi = np.array([0.5, 0.25, 0.15, 0.10]); pi /= pi.sum()
# symmetric proposal
Q = rng.random((n, n)); Q = Q + Q.T; Q /= Q.sum(1, keepdims=True)
# make it symmetric row-stochastic: use uniform off-diagonal + diag
Q = np.full((n, n), 1.0/n)
assert np.allclose(Q, Q.T)
acc = lambda c, y: min(1.0, pi[y]/pi[c])

def kernel(nprop, scheme):
    K = np.zeros((n, n))
    for x1 in range(n):
        for props in itertools.product(range(n), repeat=nprop):
            # proposal path probability
            p = 1.0
            prev = x1
            for y in props:
                p *= Q[x1, y] if scheme == 'ABS' else Q[prev, y]
                prev = y
            # sequential accept/reject: enumerate accept decisions
            for dec in itertools.product([0, 1], repeat=nprop):
                cur = x1; pd = p
                for y, d in zip(props, dec):
                    a = acc(cur, y)
                    if d: pd *= a; cur = y
                    else: pd *= (1 - a)
                K[x1, cur] += pd
    assert np.allclose(K.sum(1), 1)
    return K

for nprop in (1, 2, 3):
    for scheme in ('ABS', 'REL'):
        K = kernel(nprop, scheme)
        inv = np.abs(pi @ K - pi).max()
        db = np.abs(pi[:, None]*K - (pi[:, None]*K).T).max()
        stat = np.linalg.matrix_power(K, 2000)[0]
        print(f"nprop={nprop} {scheme}: max|piK-pi|={inv:.3e}  max DB viol={db:.3e}  "
              f"stationary={np.round(stat,4)}  target={np.round(pi,4)}")
