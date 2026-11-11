"""Referee toy (draft 3): Eq. (block) with two coupled indices. Region A fits the model. Region B does not (its block
scores are over-dispersed by a factor kappa in variance). The Hessian couples the two indices with correlation rho.
Question: does the block error of A exceed its model error although A has one index? Scores are simulated directly
(no maps): 40 blocks, 2 parameters, 20000 draws, a few seconds, a few MB.

Also checks, for one parameter and independent blocks: E[CR0] = (1 - sum h^2)/D, E[CR3] = sum h/(1-h) / D, and
CR3 = G/(G-1) x delete-one jackknife variance (linear case)."""
import numpy as np
rng = np.random.default_rng(11)
nb, ndraw = 40, 20000

def shares(neff_like, n):
    w = rng.gamma(shape=neff_like / n, size=n) + 1e-4
    return w / w.sum()

# ---- one parameter, independent blocks
h = shares(6.0, nb); D = 1.0
G = rng.standard_normal((ndraw, nb)) * np.sqrt(h * D)          # Var(G_B) = D_B
delta = -G.sum(1) / D                                          # estimate minus truth
Ghat = G + (h * D)[None, :] * delta[:, None]                   # scores at the optimum
cr0 = (Ghat ** 2).sum(1) / D ** 2
cr3 = ((Ghat / (1 - h)) ** 2).sum(1) / D ** 2
# delete-one estimates in the linear case: delta_(-B) = -(sum_{B' != B} G_B') / (D - D_B)
dl = -(G.sum(1)[:, None] - G) / (D - h * D)[None, :]
jk = (nb - 1) / nb * ((dl - delta[:, None]) ** 2).sum(1)
print(f"one parameter: n_eff {1 / (h ** 2).sum():.1f}, max h {h.max():.2f}")
print(f"  true variance {delta.var():.4f} (1/D = 1)")
print(f"  E[CR0] {cr0.mean():.4f}  predicted 1 - sum h^2 = {1 - (h ** 2).sum():.4f}")
print(f"  E[CR3] {cr3.mean():.4f}  predicted sum h/(1-h) = {(h / (1 - h)).sum():.4f}")
print(f"  CR3 / jackknife variance: {np.median(cr3 / jk):.4f}  predicted G/(G-1) = {nb / (nb - 1):.4f}")

# ---- two coupled parameters; region B misfits
for rho, kappa in ((0.2, 1.0), (0.2, 30.0), (0.2, 70.0), (0.0, 70.0)):
    hA, hB = shares(6.0, nb), shares(3.0, nb)
    # put the heaviest block of A and of B in the same block, as on the maps (block 18)
    iA, iB = hA.argmax(), hB.argmax(); hB[[iA, iB]] = hB[[iB, iA]]
    c = rho * np.sqrt(hA * hB) / np.sqrt(hA * hB).sum()        # block shares of the off-diagonal Hessian element; sum = rho
    c = np.minimum(c, 0.99 * np.sqrt(hA * hB))
    J = np.array([[1.0, c.sum()], [c.sum(), 1.0]]); Ji = np.linalg.inv(J)
    Gs = np.empty((ndraw, nb, 2))
    for b in range(nb):
        cov = np.array([[hA[b], c[b]], [c[b], hB[b]]])
        Lc = np.linalg.cholesky(cov + 1e-12 * np.eye(2))
        Gs[:, b, :] = rng.standard_normal((ndraw, 2)) @ Lc.T
    Gs[:, :, 1] *= np.sqrt(kappa)                               # over-dispersed scores for region B only
    S = Gs.sum(1)
    delta = -S @ Ji                                             # (ndraw, 2)
    Dm = np.stack([np.stack([hA, c], -1), np.stack([c, hB], -1)], -2)     # (nb, 2, 2) block Hessians
    Ghat = Gs + np.einsum("brs,ds->dbr", Dm, delta)
    lev = np.stack([hA, hB], -1)
    Gk = Ghat / (1 - lev)[None]
    M = np.einsum("dbr,dbs->drs", Gk, Gk)
    full = np.einsum("r,drs,s->d", Ji[0], M, Ji[0])             # CR3 variance of A, full matrix form (as in block.py)
    own = Ji[0, 0] ** 2 * M[:, 0, 0]                            # own-score term only
    sand = Ji[0, 0]
    print(f"rho {rho}, kappa {kappa:4.0f}: region A  model variance {sand:.3f} | true variance of the estimate {delta[:, 0].var():.3f} | mean CR3 full {full.mean():.3f}, own-score {own.mean():.3f}"
          f" | median ratio of errors: full {np.median(np.sqrt(full / sand)):.2f}, own {np.median(np.sqrt(own / sand)):.2f}; 97.5th percentile of the own ratio {np.percentile(np.sqrt(own / sand), 97.5):.2f}")
