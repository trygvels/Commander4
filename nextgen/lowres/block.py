"""Block-robust errors of the region indices from the scatter of the profile score over blocks of sky.

At the profile optimum the gradient of chi2_min / 2 with respect to the index of region r is a sum over
elements (envelope theorem),

    g_p = -sum_b (N_b^-1 r_b)_p (ds_b / dbeta)_p a_s,p ,        sum_{p in r} g_p = 0 ,

and under the model, at the truth, its covariance is J + T. Summing g over blocks B gives G_{B,r}. With
D_{B,r} = dG_{B,r} / dbeta_r (amplitudes re-solved; sum_B D_B is the observed profile Hessian) and the
leverage h_B = D_B / sum_B D_B, three covariances J_c^-1 (sum_B Gk_B Gk_B^T) J_c^-1 are returned:

    CR0   Gk = G                      E = (1 - sum h^2) / D   for independent blocks under the model: too small
    CR2   Gk = G / sqrt(1 - h)        E = 1 / D               when Var(G_B) = D_B (1 - h_B), blocks independent
    CR3   Gk = G / (1 - h)            E = sum h / (1 - h) / D ; G / (G - 1) times the delete-one-block jackknife
                                      variance when the covariances do not couple blocks; conservative

What the block error is and is not (independent statistical review, 2026-11-11,
`nextgen/reviews/2026-11-11_stats_review_prior_block_lowpass.md`):

- It is the error of the estimate as a measurement of the mean of a population of exchangeable, independent
  blocks. For the information-weighted mean index of the region as realized on this sky, the sandwich error is
  the right one.
- If index offsets are correlated between neighbouring blocks it is too small (0.6 to 0.7 of the true rms
  for a correlation length of two blocks).
- If the noise covariance or a prior couples blocks strongly (large-scale modes with tens of times the
  white variance) the block scores are negatively correlated and all three forms are too large under the
  model, by up to a factor of two. The ratio to the sandwich error under the model therefore has to be
  taken from simulations with the covariances and the prior of the fit.
- The variance estimate has fewer degrees of freedom than n_eff = (sum D)^2 / sum D^2 once a spread dominates.
- No block error is returned for a region whose largest leverage exceeds `h_max` (CR2 and CR3 are not
  defined as h approaches one) or whose debiased information is not positive.
"""
import numpy as np

from .fit import sandwich_plugin

KINDS = ("CR0", "CR2", "CR3")


def block_scores(like, beta, block, nblk):
    """G (nblk, nreg): the gradient of chi2_min / 2 with respect to each region index, split by block.
    `block` (n,) gives the block of each element. Columns sum to the gradient returned by `like.profile`."""
    sol = like.solve(beta)
    n = like.n
    g = np.zeros(n)
    for i, (b, r) in enumerate(zip(like.bands, like.residuals(sol))):
        g -= b.apply(r) * sol["ds"][i] * sol["x"][:n]
    G = np.zeros((nblk, like.nreg))
    np.add.at(G, (block, like.region), g)
    return G


def block_errors(like, beta_hat, block, nblk, h=2e-3, h_max=0.9):
    """Block-robust errors at the profile optimum `beta_hat`.

    Returns a dict with
      sd_model (nreg,)       debiased sandwich error
      CR0, CR2, CR3 (nreg,)  block errors (nan where not defined)
      neff (nreg,)           effective number of blocks, (sum D)^2 / sum D^2 over the positive D
      lev_max (nreg,)        largest leverage
      n_neg (nreg,)          number of blocks with negative information D_B (set to zero in the leverages)
      d_over_jc (nreg,)      sum_B D_B over the debiased information J_c (a check: 1 where T / J is small)
      cross (float)          share of the off-diagonal terms |dG_{B,r} / dbeta_s|, r != s, in the block Jacobian
      G (nblk, nreg), D (nblk, nreg), ok (nreg,)
    """
    _, cs, Jp, T, _, ok = sandwich_plugin(like, beta_hat)
    nreg = like.nreg
    G = block_scores(like, beta_hat, block, nblk)
    Jac = np.zeros((nblk, nreg, nreg))                       # dG_{B,r} / dbeta_s
    for s in range(nreg):
        e = np.zeros(nreg)
        e[s] = h
        Jac[:, :, s] = (block_scores(like, beta_hat + e, block, nblk) - block_scores(like, beta_hat - e, block, nblk)) / (2 * h)
    D = np.einsum("brr->br", Jac).copy()
    tot = np.abs(Jac).sum()
    cross = float((tot - np.abs(D).sum()) / tot) if tot > 0 else 0.0
    n_neg = (D < 0).sum(0)
    Dpos = np.clip(D, 0.0, None)
    dsum = Dpos.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        lev = np.where(dsum > 0, Dpos / dsum, 0.0)
        neff = np.where(dsum > 0, dsum ** 2 / (Dpos ** 2).sum(0), 0.0)
    lev_max = lev.max(0)
    Jc = Jp - T
    good = np.asarray(ok, dtype=bool) & (np.diag(Jc) > 0)
    out = {"sd_model": np.sqrt(np.diag(cs)), "neff": neff, "lev_max": lev_max, "n_neg": n_neg, "cross": cross,
           "d_over_jc": np.where(np.diag(Jc) > 0, D.sum(0) / np.where(np.diag(Jc) > 0, np.diag(Jc), 1.0), np.nan),
           "G": G, "D": D, "ok": good}
    Ji = np.zeros((nreg, nreg))
    if good.any():
        idx = np.flatnonzero(good)
        Ji[np.ix_(idx, idx)] = np.linalg.inv(Jc[np.ix_(idx, idx)])
    levc = np.clip(lev, 0.0, h_max)                           # only to keep the arithmetic finite; masked below
    for kind, f in zip(KINDS, (np.ones_like(levc), 1 / np.sqrt(1 - levc), 1 / (1 - levc))):
        Gk = G * f
        sd = np.sqrt(np.clip(np.diag(Ji @ (Gk.T @ Gk) @ Ji), 0.0, None))
        sd[~good] = np.nan
        if kind != "CR0":
            sd[lev_max > h_max] = np.nan
        out[kind] = sd
    return out
