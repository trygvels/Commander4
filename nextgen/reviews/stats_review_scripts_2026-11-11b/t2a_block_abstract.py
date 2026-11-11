"""Addition 2, claims (c) and (d) in the cleanest setting: independent blocks, one index.

Block B gives an estimate b_B ~ N(beta_0 + u_B, 1 / D_B), u_B ~ N(0, tau^2) (tau = 0: the model is true).
beta_hat = sum_B h_B b_B with h_B = D_B / sum D, the block score at the optimum is G_B = D_B (b_B - beta_hat), and

    CR0 = sum G_B^2 / D^2,   CR2 = sum G_B^2 / (1 - h_B) / D^2,   CR3 = sum G_B^2 / (1 - h_B)^2 / D^2.

Exact expectations (tau = 0): E CR0 = (1 - sum h^2) / D, E CR2 = 1 / D, E CR3 = sum h / (1 - h) / D.
Satterthwaite degrees of freedom of CR2 (tau = 0): nu = 1 / (sum h^2 + sum_{B != B'} h_B^2 h_B'^2 / ((1 - h_B)(1 - h_B'))).
The delete-one-block jackknife is exact here: beta_(-B) - beta_hat = -(G_B / D) / (1 - h_B).
"""
import numpy as np

rng = np.random.default_rng(3)
NMC = 400000


def shares(kind):
    if kind.startswith("equal"):
        G = int(kind[5:])
        return np.full(G, 1.0 / G)
    G, decay = {"conc4": (4, 0.33), "conc15": (15, 0.62), "conc23": (23, 0.78), "conc44": (44, 0.80)}[kind]
    h = decay ** np.arange(G)
    return h / h.sum()


def run(h, tau, clip=0.9):
    G = h.size
    D = h * 1.0                                         # total information 1: model variance 1
    b = rng.standard_normal((NMC, G)) * np.sqrt(tau ** 2 + 1.0 / D)[None, :]
    bh = b @ h
    Gs = D[None, :] * (b - bh[:, None])
    hc = np.clip(h, 0.0, clip)
    v0 = np.sum(Gs ** 2, axis=1)
    v2 = np.sum(Gs ** 2 / (1 - hc)[None, :], axis=1)
    v3 = np.sum(Gs ** 2 / (1 - hc)[None, :] ** 2, axis=1)
    jk = np.sum((Gs / (1 - h)[None, :]) ** 2, axis=1) * (G - 1) / G        # delete-one jackknife variance (exact deletion)
    true_var = tau ** 2 * np.sum(h ** 2) + 1.0
    return bh, v0, v2, v3, jk, true_var


def nu_satt(h, tau=0.0):
    """Satterthwaite 2 E[V]^2 / Var(V) for CR2 with G ~ N(0, D (I - h 1^T) S (I - 1 h^T) D), S = diag(tau^2 + 1/D)."""
    G = h.size
    S = np.diag(tau ** 2 + 1.0 / h)
    P = np.eye(G) - np.outer(np.ones(G), h)             # b - beta_hat = P b
    C = np.diag(h) @ P @ S @ P.T @ np.diag(h)           # covariance of the block scores
    out = []
    for p in (0.0, 1.0, 2.0):
        Lam = np.diag(1.0 / (1 - h) ** p)
        m = np.trace(Lam @ C)
        v = 2 * np.trace(Lam @ C @ Lam @ C)
        out.append((m, 2 * m * m / v))
    return out


print("# t2a: independent blocks, model true (tau = 0). Variances in units of the true variance 1 / D.")
print(f"{'shares':8s} {'G':>3s} {'n_eff':>6s} {'max h':>6s} | {'E CR0':>6s} {'E CR2':>6s} {'E CR3':>6s} (exact {'CR0':>5s} {'CR3':>5s}) | nu Satterthwaite CR0, CR2, CR3 | "
      f"rms z: CR0 CR2 CR3 | sqrt(nu/(nu-2)) with nu = n_eff, and with nu(CR3) and the CR3 bias | 68 pct of |z|: CR0 CR2 CR3 | CR3 / jackknife * (G-1)/G")
for kind in ("equal4", "equal8", "equal16", "equal44", "conc4", "conc15", "conc23", "conc44"):
    h = shares(kind)
    G = h.size
    neff = 1 / np.sum(h ** 2)
    bh, v0, v2, v3, jk, tv = run(h, 0.0)
    sat = nu_satt(h)
    z = [bh / np.sqrt(v) for v in (v0, v2, v3)]
    f = lambda nu: np.sqrt(nu / (nu - 2)) if nu > 2 else np.inf
    print(f"{kind:8s} {G:3d} {neff:6.2f} {h.max():6.2f} | {v0.mean():6.3f} {v2.mean():6.3f} {v3.mean():6.3f} (exact {1 - np.sum(h ** 2):5.3f} {np.sum(h / (1 - h)):5.3f}) | "
          f"{sat[0][1]:5.2f} {sat[1][1]:5.2f} {sat[2][1]:5.2f} | " + " ".join(f"{np.sqrt(np.mean(x ** 2)):5.2f}" for x in z)
          + f" | {f(neff):5.2f} {f(sat[2][1]) / np.sqrt(sat[2][0]):5.2f} | " + " ".join(f"{np.percentile(np.abs(x), 68.27):5.2f}" for x in z)
          + f" | {np.median(v3 / jk) * (G - 1) / G:.4f}")
print("# the last column is 1 when no leverage is clipped: CR3 = G / (G - 1) times the delete-one-block jackknife variance, exactly, in this linear setting")
print()
print("# random effect: true index of block B is beta_0 + tau u_B; tau^2 in units of the model variance of the region estimate times n_eff")
print("# (tau^2 sum h^2 = k / D: the spread adds k times the model variance). rms about beta_0 against the quoted errors.")
print(f"{'shares':8s} {'k':>4s} | true sd | sqrt(E CR0) sqrt(E CR2) sqrt(E CR3) | nu Satterthwaite CR3 | rms z: CR0 CR2 CR3 | 68 pct of |z|: CR0 CR2 CR3 | coverage of +-1 CR3 sd")
for kind in ("equal8", "conc4", "conc15", "conc23", "conc44"):
    h = shares(kind)
    for k in (1.0, 9.0):
        tau = np.sqrt(k / np.sum(h ** 2))
        bh, v0, v2, v3, jk, tv = run(h, tau)
        sat = nu_satt(h, tau)
        z = [bh / np.sqrt(v) for v in (v0, v2, v3)]
        print(f"{kind:8s} {k:4.0f} | {np.sqrt(tv):7.3f} | {np.sqrt(v0.mean()):11.3f} {np.sqrt(v2.mean()):11.3f} {np.sqrt(v3.mean()):11.3f} | {sat[2][1]:5.2f} | "
              + " ".join(f"{np.sqrt(np.mean(x ** 2)):6.2f}" for x in z) + " | " + " ".join(f"{np.percentile(np.abs(x), 68.27):5.2f}" for x in z)
              + f" | {np.mean(np.abs(z[2]) < 1):.3f}")
print()
print("# clipping of the leverage at 0.9: one block with share h_1, the rest equal (G = 6), model true")
for h1 in (0.5, 0.8, 0.9, 0.95, 0.99):
    h = np.concatenate([[h1], np.full(5, (1 - h1) / 5)])
    bh, v0, v2, v3, jk, tv = run(h, 0.0)
    print(f"   h_1 = {h1:.2f}: n_eff {1 / np.sum(h ** 2):.2f}; E CR0 {v0.mean():.3f}, E CR2 {v2.mean():.3f}, E CR3 (clipped) {v3.mean():.3f}; exact unclipped CR3 {np.sum(h / (1 - h)):.2f}; "
          f"68 pct of |z| with CR3: {np.percentile(np.abs(bh / np.sqrt(v3)), 68.27):.2f}; coverage of +-1 CR3 sd {np.mean(np.abs(bh / np.sqrt(v3)) < 1):.3f}")
