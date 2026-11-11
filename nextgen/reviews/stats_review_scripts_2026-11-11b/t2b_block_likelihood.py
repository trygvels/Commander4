"""Addition 2 on the real code path: `block_scores` and `block_errors` of e040f_block_errors_v2.py on a small
RegionLikelihood (96 pixels x (Q, U) = 192 elements per band, five bands, 24 blocks of 4 pixels, two regions, one
block straddling the boundary), with the information concentrated in a few blocks.

Arms:  white   diagonal noise, white CMB prior, CMB drawn white             (block scores independent under the model)
       dense   noise correlated over two block lengths in four bands, dense CMB prior, CMB drawn from it
       modes   a sky coherent over many blocks and eight large-scale noise modes of 50 times the white variance
each with a spread tau of the true index between blocks of 0 (model true) or 0.1.

    python -I t2b_block_likelihood.py <arm> <tau> <nreal> [njack] [correlation length of the offsets in blocks]
"""
import sys
import time

import numpy as np

sys.path.insert(0, "nextgen")
sys.path.insert(0, "nextgen/scripts")
sys.path.insert(0, "nextgen/reviews/stats_review_scripts_2026-11-11b")
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
import e040f_block_errors_v2 as be  # noqa: E402
import rv2_common as rv  # noqa: E402

arm, tau, nreal = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
njack = int(sys.argv[4]) if len(sys.argv) > 4 else 0
corr = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0          # correlation length of the index offsets, in blocks (0: independent)
NP, NBLK = 96, 24
n = 2 * NP
pix = np.arange(NP)
block = np.tile(pix // 4, 2)
region = np.tile((pix >= 66).astype(int), 2)                       # block 16 (pixels 64..67) straddles the boundary
beta_reg = np.array([-3.10, -2.95])
freqs = rv.FREQS
rng0 = np.random.default_rng(17)
scale = np.array([4, 4, 3.5, 3.5, 3, 3, 2.5, 2.5, 2, 2, 2, 1.5, 1.5, 1.5, 1, 1, 2, 5, 2, 3.5, 1.5, 2.5, 1.5, 1.0])   # synchrotron rms per block: a few carry the information
a30 = rng0.standard_normal(n) * scale[block] * 1.2
a_d = rng0.standard_normal(n) * 4.0
fd = [rv.dust_coef(f) for f in freqs]
dist = np.abs(pix[:, None] - pix[None, :])
dist = np.minimum(dist, NP - dist)
K1 = np.exp(-dist / 8.0)
K = np.block([[K1, 0 * K1], [0 * K1, K1]])
sig = 0.2
if arm == "white":
    covs = [sig ** 2 * np.eye(n) for _ in freqs]
    Cc = 0.14 ** 2 * np.eye(n)
elif arm == "modes":
    # a sky that is coherent over many blocks, and eight large-scale noise modes with 50 times the white variance in
    # the four low bands (the generalized least squares fit gives them no weight); white CMB prior
    a30 = (1.0 + 0.3 * rng0.standard_normal(n)) * scale[block] * 1.2 * np.concatenate([np.ones(NP), -np.ones(NP)])
    um = []
    for k_ in (1, 2):
        for fn in (np.cos, np.sin):
            v = fn(2 * np.pi * k_ * pix / NP) * np.sqrt(2.0 / NP)
            um += [np.concatenate([v, 0 * v]), np.concatenate([0 * v, v])]
    um = np.array(um).T
    covs = [sig ** 2 * (np.eye(n) + 50.0 * um @ um.T) if i < 4 else sig ** 2 * np.eye(n) for i in range(len(freqs))]
    Cc = 0.14 ** 2 * np.eye(n)
else:
    covs = [sig ** 2 * (0.4 * np.eye(n) + 0.6 * K) if i < 4 else sig ** 2 * np.eye(n) for i in range(len(freqs))]
    Kc = np.block([[np.exp(-dist / 16.0), 0 * K1], [0 * K1, np.exp(-dist / 16.0)]])
    Cc = 0.14 ** 2 * (0.1 * np.eye(n) + 0.9 * Kc)
chol = [np.linalg.cholesky(c) for c in covs]
ninv = [np.linalg.inv(c) for c in covs]
ninv = [0.5 * (m + m.T) for m in ninv]
Lc = np.linalg.cholesky(Cc)
prec = np.linalg.inv(Cc)
prec = 0.5 * (prec + prec.T)


bd = np.abs(np.arange(NBLK)[:, None] - np.arange(NBLK)[None, :])
bd = np.minimum(bd, NBLK - bd)
OFFCHOL = np.linalg.cholesky(np.exp(-bd / corr) + 1e-12 * np.eye(NBLK)) if corr > 0 else np.eye(NBLK)


def make_like(dv, keep=None):
    if keep is None:
        bands = [BandData(f"{int(f)}", d, ni, rv.sync_fun(f), fdi) for f, d, ni, fdi in zip(freqs, dv, ninv, fd)]
        return RegionLikelihood(bands, region, cmb_prior_prec=prec)
    # delete data: the likelihood of the remaining elements (sub-blocks of the covariances, not of the precisions)
    bands = []
    for f, d, c, fdi in zip(freqs, dv, covs, fd):
        ni = np.linalg.inv(c[np.ix_(keep, keep)])
        bands.append(BandData(f"{int(f)}", d[keep], 0.5 * (ni + ni.T), rv.sync_fun(f), fdi))
    pk = np.linalg.inv(Cc[np.ix_(keep, keep)])
    return RegionLikelihood(bands, region[keep], cmb_prior_prec=0.5 * (pk + pk.T))


def simulate(r):
    rng = np.random.default_rng([99, int(round(100 * tau)), r])
    off = tau * (OFFCHOL @ rng.standard_normal(NBLK))
    beta_el = beta_reg[region] + off[block]
    a_s = a30 * (23.0 / 30.0) ** beta_el
    a_c = Lc @ rng.standard_normal(n)
    out = []
    for i, f in enumerate(freqs):
        out.append(np.exp(beta_el * np.log(f / rv.NU_REF)) * a_s + fd[i] * a_d + a_c + chol[i] @ rng.standard_normal(n))
    return out, off


def full_leverage(L, bh, h=2e-3):
    """Block errors with the full Jacobian D_B = dG_B/dbeta^T (nreg x nreg per block) and the observed Hessian
    H = sum_B D_B as the bread: CR0 = H^-1 sum G G^T H^-T, CR3 with (1 - D_B H^-1)^-1 G_B, no clipping."""
    G = be.block_scores(L, bh, block, NBLK)
    nr = L.nreg
    Dfull = np.zeros((NBLK, nr, nr))
    for s in range(nr):
        e = np.zeros(nr)
        e[s] = h
        Dfull[:, :, s] = (be.block_scores(L, bh + e, block, NBLK) - be.block_scores(L, bh - e, block, NBLK)) / (2 * h)
    H = Dfull.sum(0)
    Hi = np.linalg.inv(H)
    G3 = np.array([np.linalg.solve(np.eye(nr) - Dfull[b] @ Hi, G[b]) for b in range(NBLK)])
    cr0 = Hi @ (G.T @ G) @ Hi.T
    cr3 = Hi @ (G3.T @ G3) @ Hi.T
    off = np.abs(Dfull[:, 0, 1]).sum() + np.abs(Dfull[:, 1, 0]).sum()
    return np.sqrt(np.diag(cr0)), np.sqrt(np.diag(cr3)), H, G, Dfull, off / (np.abs(Dfull[:, 0, 0]).sum() + np.abs(Dfull[:, 1, 1]).sum())


print(f"# t2b: arm {arm}, spread of the true index between blocks {tau} (correlation length {corr} blocks), {nreal} realizations; {n} elements, {NBLK} blocks, regions of "
      f"{np.bincount(region).tolist()} elements, truth {beta_reg}", flush=True)
# exact covariance of the element scores at the truth under the model (tau = 0), from the reviewer's reference:
#   Cov(g_p, g_q) = a_p X_pq a_q + [At^-1]_ss,pq X_pq,  X = Dm^T Wt Dm.
# Independent blocks would make the block-to-block part vanish: sum_B Var(G_B) = Var(sum_B G_B).
ref = rv.Reference(freqs, covs, region, cmb_cov=Cc)
M_, Dm_, At_, Ati_, Wt_ = ref.parts(beta_reg)
X_ = Dm_.T @ Wt_ @ Dm_
a_true = a30 * (23.0 / 30.0) ** beta_reg[region]
Kel = a_true[:, None] * X_ * a_true[None, :] + Ati_[:n, :n] * X_
line = []
for k in range(2):
    sel = region == k
    Kr = Kel[np.ix_(sel, sel)]
    bl_ = block[sel]
    same = bl_[:, None] == bl_[None, :]
    Db = np.array([Kr[bl_ == B].sum() for B in np.unique(bl_)])                 # Cov(G_B, sum G): expected block information
    Vb = np.array([Kr[np.ix_(bl_ == B, bl_ == B)].sum() for B in np.unique(bl_)])
    line.append(f"region {k}: sum_B Var(G_B) / Var(sum_B G_B) = {Kr[same].sum() / Kr.sum():.3f}; blocks with Cov(G_B, sum G) < 0: {int(np.sum(Db < 0))} of {Db.size}; "
                f"largest Var(G_B) / Cov(G_B, sum G) among blocks with more than 5 percent of the information: {np.max((Vb / Db)[Db > 0.05 * Db.sum()]):.2f}")
print("## exact, model true: " + " | ".join(line), flush=True)
if nreal == 0:
    sys.exit(0)
t0 = time.time()
rows, jrows, chk = [], [], []
for r in range(nreal):
    dv, off = simulate(r)
    L = make_like(dv)
    bh, res = fit_profile(L, beta_reg)
    s_model, s_blk, neff, levmax, dsum = be.block_errors(L, bh, block, NBLK)
    f0, f3, H, G, Dfull, offrel = full_leverage(L, bh)
    Jp, T, _ = L.information(bh)
    Jc = Jp - T
    if r < 3:
        # claim (a): the sum of the element scores is the gradient of chi2_min / 2, and it vanishes at the optimum
        g_code = L.profile(bh, grad=True)[1]
        g_fd = rv.fd_grad(L.profile, bh, h=1e-5)
        b1 = bh + np.array([0.03, -0.02])
        G1 = be.block_scores(L, b1, block, NBLK)
        Hfd = np.array([(L.profile(bh + e, grad=True)[1] - L.profile(bh - e, grad=True)[1]) / 4e-3 for e in (np.array([2e-3, 0]), np.array([0, 2e-3]))]).T
        chk.append((np.abs(G.sum(0)).max() / np.abs(G).max(), np.abs(G1.sum(0) - rv.fd_grad(L.profile, b1, h=1e-5)).max() / np.abs(G1.sum(0)).max(),
                    np.abs(np.diag(H) - np.diag(Hfd)).max() / np.abs(np.diag(Hfd)).max(), np.diag(H) / np.diag(Jc)))
    # realized information-weighted mean of the block indices, first order
    Dd = np.array([Dfull[:, k, k] for k in range(2)]).T
    wmean = beta_reg + (Dd * off[:, None]).sum(0) / Dd.sum(0)
    rows.append(np.concatenate([bh, s_model, s_blk["CR0"], s_blk["CR2"], s_blk["CR3"], f0, f3, neff, levmax, dsum, np.diag(H) / np.diag(Jc), [offrel],
                                wmean, [float(np.any(Dd < 0))]]))
    if r < njack:
        # the true delete-one-block jackknife: refit without the data of block B
        d = np.zeros((NBLK, 2))
        for B in range(NBLK):
            keep = np.where(block != B)[0]
            Lk = make_like(dv, keep)
            d[B] = fit_profile(Lk, bh)[0]
        # blocks that hold no element of a region leave that region's other blocks; all 12 deletions are used
        vj = (NBLK - 1) / NBLK * np.sum((d - d.mean(0)[None, :]) ** 2, axis=0)
        vj0 = np.sum((d - bh[None, :]) ** 2, axis=0)                              # about the full-data estimate, no factor
        jrows.append(np.concatenate([np.sqrt(vj), np.sqrt(vj0), s_blk["CR3"], f3]))
    if r % 500 == 0:
        print(f"  realization {r}: index {bh.round(4)}, sandwich {s_model.round(4)}, CR3 {s_blk['CR3'].round(4)}, n_eff {neff.round(1)} [{time.time() - t0:.0f} s]", flush=True)
a = np.array(rows)
est, sm = a[:, 0:2], a[:, 2:4]
c0, c2, c3, f0, f3 = a[:, 4:6], a[:, 6:8], a[:, 8:10], a[:, 10:12], a[:, 12:14]
neff, levmax, dsum, hsum, offrel, wmean, negD = a[:, 14:16], a[:, 16:18], a[:, 18:20], a[:, 20:22], a[:, 22], a[:, 23:25], a[:, 25]
print("## claim (a), first three realizations: |sum_B G_B| / max |G_B| at the optimum; |sum G - finite-difference gradient| / |gradient| away from it; "
      "|sum_B D_B - finite-difference Hessian| / Hessian; sum_B D_B / J_c per region")
for c in chk:
    print(f"   {c[0]:.1e}   {c[1]:.1e}   {c[2]:.1e}   " + " ".join(f"{x:.3f}" for x in c[3]))
print(f"## claim (b): sum_B D_B / J_c over the realizations (code's printed check): median {np.median(dsum, axis=0).round(3)}, 2.5 and 97.5 percentiles "
      f"{np.percentile(dsum, 2.5, axis=0).round(3)} {np.percentile(dsum, 97.5, axis=0).round(3)}")
print(f"## n_eff median {np.median(neff, axis=0).round(2)}; largest leverage median {np.median(levmax, axis=0).round(2)}; realizations with a negative D_B {int(negD.sum())}; "
      f"with a leverage at the 0.9 clip {int(np.sum(levmax >= 0.9 - 1e-12))}; cross-region share of |dG_B,r/dbeta_s| {np.median(offrel):.4f}")
truth = beta_reg
scat = est.std(0, ddof=1)
rmsd = np.sqrt(np.mean((est - truth) ** 2, axis=0))
print(f"## estimate: mean - region truth {(est - truth).mean(0).round(5)} +- {(scat / np.sqrt(len(a))).round(5)}; rms about the region truth {rmsd.round(5)}; "
      f"rms about the realized information-weighted mean of the block indices {np.sqrt(np.mean((est - wmean) ** 2, axis=0)).round(5)}")
print(f"{'error':42s} {'sqrt(mean var)/rms':>20s} {'median/rms':>14s} {'rms of z':>14s} {'68 pct of |z|':>14s} {'cover +-1':>12s} {'cover +-2':>12s}")
for name, s in (("sandwich (model)", sm), ("CR0 (code)", c0), ("CR2 (code)", c2), ("CR3 (code)", c3), ("CR0, observed Hessian, full Jacobian", f0), ("CR3, observed Hessian, full Jacobian", f3)):
    z = (est - truth) / s
    print(f"{name:42s} " + " ".join(f"{x:9.3f}" for x in np.sqrt(np.mean(s ** 2, axis=0)) / rmsd) + "  " + " ".join(f"{x:6.3f}" for x in np.median(s, axis=0) / rmsd)
          + "  " + " ".join(f"{x:6.2f}" for x in np.sqrt(np.mean(z ** 2, axis=0))) + "  " + " ".join(f"{x:6.2f}" for x in np.percentile(np.abs(z), 68.27, axis=0))
          + "  " + " ".join(f"{x:5.3f}" for x in np.mean(np.abs(z) < 1, axis=0)) + "  " + " ".join(f"{x:5.3f}" for x in np.mean(np.abs(z) < 2, axis=0)))
f = lambda nu: np.sqrt(nu / (nu - 2)) if nu > 2 else np.inf
print(f"   sqrt(nu / (nu - 2)) at nu = median n_eff: " + " ".join(f"{f(x):.2f}" for x in np.median(neff, axis=0)) + "; at n_eff - 1: " + " ".join(f"{f(x - 1):.2f}" for x in np.median(neff, axis=0))
      + "; a t distribution has 68 percent of |z| below 1.00 (nu large), 1.14 (nu 4), 1.20 (nu 3), 1.32 (nu 2)")
print(f"   CR3 (code) / CR3 (observed Hessian, full Jacobian): median {np.median(c3 / f3, axis=0).round(4)}, 2.5 to 97.5 percent {np.percentile(c3 / f3, 2.5, axis=0).round(3)} {np.percentile(c3 / f3, 97.5, axis=0).round(3)}")
if jrows:
    j = np.array(jrows)
    print(f"## true delete-one-block jackknife (data of the block removed, refit), {len(j)} realizations, {NBLK} blocks:")
    for nm, col in (("CR3 (code)", j[:, 4:6]), ("CR3 (observed Hessian, full Jacobian)", j[:, 6:8])):
        q = col / j[:, 2:4]
        print(f"   {nm} / [sum_B (beta_(-B) - beta_hat)^2]^(1/2): median " + " ".join(f"{x:.3f}" for x in np.median(q, axis=0)) + "; 16 to 84 percent "
              + " ".join(f"{x:.3f}" for x in np.percentile(q, 16, axis=0)) + " to " + " ".join(f"{x:.3f}" for x in np.percentile(q, 84, axis=0)))
    q = j[:, 4:6] / j[:, 0:2]
    print("   CR3 (code) / jackknife sd with the factor (G - 1) / G: median " + " ".join(f"{x:.3f}" for x in np.median(q, axis=0)) + f"; sqrt(G / (G - 1)) = {np.sqrt(NBLK / (NBLK - 1)):.3f}")
np.savez(f"nextgen/reviews/stats_review_scripts_2026-11-11b/out_t2b_{arm}_tau{tau:g}_corr{corr:g}.npz", rows=a, jack=np.array(jrows))
print(f"# done, {time.time() - t0:.0f} s")
