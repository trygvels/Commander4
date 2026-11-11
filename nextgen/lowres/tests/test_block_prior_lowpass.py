"""Exact and Monte Carlo checks of the three additions of 2026-11-11: block-robust errors (`lowres.block`), a
dense prior on the CMB amplitudes (`cmb_prior_prec`, `lowres.window.lcdm_cmb_precision`), and the low-pass
construction used to read the index from large scales. The independent review that asked for these tests is
`nextgen/reviews/2026-11-11_stats_review_prior_block_lowpass.md`. Run: pytest -q nextgen/lowres/tests"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))
from lowres.block import block_errors, block_scores  # noqa: E402
from lowres.fit import fit_profile  # noqa: E402
from lowres.fit import sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band  # noqa: E402
from test_likelihood import fd, make  # noqa: E402

NUS = [22.8, 28.4, 33.0, 40.7, 60.8, 353.0]


def setup(n, nreg, seed, noise=2.0, amp=6.0, dense=()):
    """Regions, true amplitudes and the bands (spectrum object, dust response, noise covariance). The data are
    drawn by `draw`, so that a test controls the CMB and the index of every element."""
    rng = np.random.default_rng(seed)
    region = np.arange(n) % nreg
    a_s = amp * rng.standard_normal(n) * (1.0 + 3.0 * rng.random(n))
    a_d = 0.5 * rng.standard_normal(n)
    specs = []
    for i, nu in enumerate(NUS):
        B = Band.delta(f"b{i}", nu)
        sig = noise * (0.2 if nu > 300 else rng.uniform(1.0, 2.5))
        if i in dense:
            X = rng.standard_normal((n, n))
            Cv = sig ** 2 * (np.eye(n) + 0.5 * X @ X.T / n)
            specs.append((B, B.mbb(1.55, 19.6, 353.0), np.linalg.inv(Cv), np.linalg.cholesky(Cv)))
        else:
            v = sig ** 2 * rng.uniform(0.6, 1.6, n)
            specs.append((B, B.mbb(1.55, 19.6, 353.0), 1.0 / v, np.sqrt(v)))
    return region, a_s, a_d, specs


def draw(specs, a_s, a_d, a_c, beta_el, rng):
    """One data set: synchrotron with index beta_el per element, dust, CMB a_c, and a new noise draw."""
    out = []
    for B, dust, ninv, root in specs:
        xi = rng.standard_normal(a_s.size)
        noise = root @ xi if root.ndim == 2 else root * xi
        d = B.power_law(beta_el, 23.0)[0] * a_s + dust * a_d + a_c + noise
        out.append(BandData(B.name, d, ninv, (lambda b, B=B: B.power_law(b, 23.0)), dust))
    return out


def test_block_scores_sum_to_gradient_and_block_information_sums_to_hessian():
    build, region, _ = make(n=96, nreg=3, dense=(0, 2), seed=1)
    L = RegionLikelihood(build(), region, cmb_prior_sd=0.3)
    block = np.arange(96) // 8                                  # 12 blocks, each holding all three regions
    beta = np.array([-3.0, -3.2, -2.95])
    G = block_scores(L, beta, block, 12)
    _, grad = L.profile(beta, grad=True)
    assert np.allclose(G.sum(0), grad, rtol=1e-9, atol=1e-9), (G.sum(0), grad)
    bh, _ = fit_profile(L, np.full(3, -3.1))
    out = block_errors(L, bh, block, 12)
    assert np.abs(out["G"].sum(0)).max() < 1e-4 * np.abs(out["G"]).max()        # the score vanishes at the optimum
    h = 2e-3
    for r in range(3):
        e = np.zeros(3)
        e[r] = h
        hess_rr = (L.profile(bh + e, grad=True)[1][r] - L.profile(bh - e, grad=True)[1][r]) / (2 * h)
        assert np.isclose(out["D"].sum(0)[r], hess_rr, rtol=1e-8), (r, out["D"].sum(0)[r], hess_rr)
    assert np.all(np.abs(out["d_over_jc"] - 1.0) < 0.15), out["d_over_jc"]      # observed Hessian against J(a_hat) - T
    assert 0.0 <= out["cross"] < 0.2


def _one_region(n=240, seed=4):
    """Diagonal noise and a white CMB prior: the blocks are independent. One region; blocks of unequal size."""
    build, region, _ = make(n=n, nreg=1, dense=(), noise=2.0, amp=6.0, seed=seed)
    edges = np.array([0, 72, 120, 150, 174, 192, 204, 216, 224, 232, 240])          # 10 blocks, the first holds 30 percent
    block = np.searchsorted(edges, np.arange(n), side="right") - 1
    return build, region, block


def test_block_errors_expectations_for_independent_blocks():
    """Under the model with independent blocks: E[CR2] is the variance of the estimate, E[CR0] is smaller by
    1 - sum h^2, E[CR3] is larger, and the debiased sandwich is right."""
    build, region, block = _one_region()
    rng = np.random.default_rng(21)
    est, v = [], {k: [] for k in ("model", "CR0", "CR2", "CR3")}
    lev2 = []
    for _ in range(300):
        L = RegionLikelihood(build(rng), region, cmb_prior_sd=0.3)
        bh, _ = fit_profile(L, np.array([-3.1]))
        out = block_errors(L, bh, block, 10)
        est.append(bh[0])
        v["model"].append(out["sd_model"][0] ** 2)
        for k in ("CR0", "CR2", "CR3"):
            v[k].append(out[k][0] ** 2)
        lev2.append(1.0 / out["neff"][0])
    var = np.var(est, ddof=1)
    m = {k: np.mean(x) for k, x in v.items()}
    sum_h2 = np.mean(lev2)
    assert abs(m["model"] / var - 1) < 0.2, (m, var)
    assert abs(m["CR2"] / var - 1) < 0.2, (m, var)
    assert abs(m["CR0"] / var - (1 - sum_h2)) < 0.2, (m["CR0"] / var, 1 - sum_h2)
    assert m["CR0"] < m["CR2"] < m["CR3"]
    assert abs(np.mean(est) + 3.1) < 4 * np.sqrt(var / len(est)) + 0.01


def test_block_errors_follow_a_true_spread_between_blocks():
    """A random index offset per block (new in every realization): the scatter of the estimate about the region
    mean is far above the sandwich error, and the block error follows it."""
    n, tau = 240, 0.15
    region, a_s, a_d, specs = setup(n, 1, seed=4, noise=0.5)         # high signal to noise: the spread dominates
    edges = np.array([0, 72, 120, 150, 174, 192, 204, 216, 224, 232, 240])
    block = np.searchsorted(edges, np.arange(n), side="right") - 1
    rng = np.random.default_rng(22)
    est, cr3, mod = [], [], []
    for _ in range(200):
        beta_el = -3.1 + tau * rng.standard_normal(10)[block]
        L = RegionLikelihood(draw(specs, a_s, a_d, 0.2 * rng.standard_normal(n), beta_el, rng), region, cmb_prior_sd=0.2)
        bh, _ = fit_profile(L, np.array([-3.1]))
        out = block_errors(L, bh, block, 10)
        est.append(bh[0])
        cr3.append(out["CR3"][0])
        mod.append(out["sd_model"][0])
    rms = np.sqrt(np.mean((np.array(est) + 3.1) ** 2))
    assert rms > 1.5 * np.median(mod), (rms, np.median(mod))                 # the model error is far too small
    ratio = np.sqrt(np.mean(np.square(cr3))) / rms
    assert 0.8 < ratio < 1.6, ratio


def test_no_leverage_corrected_error_when_one_block_carries_the_region():
    build, region, _ = make(n=120, nreg=1, dense=(), seed=6)
    block = np.where(np.arange(120) < 116, 0, 1 + (np.arange(120) - 116))    # block 0 holds 97 percent of the elements
    L = RegionLikelihood(build(), region, cmb_prior_sd=0.3)
    bh, _ = fit_profile(L, np.array([-3.1]))
    out = block_errors(L, bh, block, 5)
    assert out["lev_max"][0] > 0.9
    assert np.isfinite(out["CR0"][0]) and np.isnan(out["CR2"][0]) and np.isnan(out["CR3"][0])


def test_dense_cmb_prior_reduces_to_white_and_has_right_gradients():
    build, region, _ = make(n=64, nreg=4, dense=(0, 2), seed=2)
    bands = build()
    beta = np.array([-3.0, -3.2, -2.9, -3.15])
    sd = 0.3
    Lw = RegionLikelihood(bands, region, cmb_prior_sd=sd)
    Ld = RegionLikelihood(bands, region, cmb_prior_prec=np.eye(64) / sd ** 2)
    assert np.isclose(Lw.profile(beta), Ld.profile(beta), rtol=1e-12)
    assert np.isclose(Lw.marginal(beta), Ld.marginal(beta), rtol=1e-12)
    Jw, Tw, _ = Lw.information(beta)
    Jd, Td, _ = Ld.information(beta)
    assert np.allclose(Jw, Jd, rtol=1e-10) and np.allclose(Tw, Td, rtol=1e-10)
    rng = np.random.default_rng(3)
    X = rng.standard_normal((64, 64))
    P = X @ X.T / 64 + 0.05 * np.eye(64)                                      # dense prior precision, condition number of order 100
    L = RegionLikelihood(bands, region, cmb_prior_prec=4.0 * P)
    f, g = L.profile(beta, grad=True)
    assert np.allclose(g, fd(L.profile, beta), rtol=2e-5, atol=1e-6)
    f, g = L.marginal(beta, grad=True)
    assert np.allclose(g, fd(L.marginal, beta), rtol=2e-5, atol=1e-6)


def test_profile_is_calibrated_when_the_cmb_is_drawn_from_a_dense_prior():
    n = 96
    region, a_s, a_d, specs = setup(n, 2, seed=7, dense=(0,))
    rng = np.random.default_rng(8)
    X = rng.standard_normal((n, 6))
    cov = 0.25 * X @ X.T / 6 + 0.01 * np.eye(n)                              # six strong modes and a floor
    prec, root = np.linalg.inv(cov), np.linalg.cholesky(cov)
    truth = np.array([-3.1, -3.0])
    z = []
    for _ in range(250):
        L = RegionLikelihood(draw(specs, a_s, a_d, root @ rng.standard_normal(n), truth[region], rng), region, cmb_prior_prec=prec)
        bh, _ = fit_profile(L, truth)
        _, cs, _, _, _, ok = sandwich_plugin(L, bh)
        z.append((bh - truth) / np.sqrt(np.diag(cs)))
    z = np.array(z)
    assert np.all(np.abs(z.mean(0)) < 4 / np.sqrt(len(z))), z.mean(0)
    assert np.all(np.abs(z.std(0, ddof=1) - 1.0) < 0.2), z.std(0, ddof=1)


def test_lcdm_covariance_is_consistent():
    pytest.importorskip("camb")
    from lowres.window import WindowMatcher, lcdm_cmb_precision
    WM = WindowMatcher(2, target_arcmin=60.0)
    prec, cov = lcdm_cmb_precision(WM)
    n = cov.shape[0]
    assert np.allclose(cov, cov.T) and np.linalg.eigvalsh(cov).min() > 0
    assert np.abs(prec @ cov - np.eye(n)).max() < 1e-8
    # the same matrix from the basis and the spectra, written out here
    import camb
    import healpy as hp
    pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=300)
    cls = camb.get_results(pars).get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]
    lmax = WM.lmax
    w2 = (hp.gauss_beam(np.radians(1.0), lmax=lmax, pol=True)[:, 1] * hp.pixwin(2, pol=True, lmax=lmax)[1]) ** 2
    omega = 4 * np.pi / hp.nside2npix(2)
    nE = WM.Y.shape[1] // 2
    assert np.all(WM.ell[:nE] == WM.ell[nE:])                                # E block then B block, same multipoles
    cl = np.concatenate([(cls[:lmax + 1, 1] * w2)[WM.ell[:nE]], (cls[:lmax + 1, 2] * w2)[WM.ell[nE:]]])
    ref = (WM.Y * (cl / omega)) @ WM.Y.T + 0.01 ** 2 * np.eye(n)
    assert np.allclose(cov, ref, rtol=1e-10, atol=1e-14)
    # E modes dominate: the variance carried by the E columns is far larger than that of the B columns
    assert cl[:nE].sum() > 50 * cl[nE:].sum()


def test_lowpass_projection_identities():
    """The low-pass band (d' = Y_k Y_k^+ d, N'^-1 = Y_k^{+T} (Y_k^+ N Y_k^{+T})^-1 Y_k^+, with Y_k^+ the rows of the
    full pseudo-inverse) is exactly the Gaussian likelihood of the kept coefficients, blind to the other columns
    and to everything outside the span of the basis."""
    rng = np.random.default_rng(5)
    n, m, k = 60, 20, 8
    Y = rng.standard_normal((n, m)) + 0.3 * rng.standard_normal((n, 1))      # full column rank, columns not orthogonal
    Yp = np.linalg.solve(Y.T @ Y, Y.T)
    Yk, Ypk = Y[:, :k], Yp[:k]
    assert np.abs(Ypk @ Yk - np.eye(k)).max() < 1e-10
    assert np.abs(Ypk @ Y[:, k:]).max() < 1e-10
    out = rng.standard_normal(n)
    out -= Y @ (Yp @ out)                                                    # a vector outside the span
    assert np.abs(Ypk @ out).max() < 1e-10
    X = rng.standard_normal((n, n))
    N = X @ X.T / n + 0.5 * np.eye(n)
    C = Ypk @ N @ Ypk.T
    Nlow = Ypk.T @ np.linalg.inv(C) @ Ypk
    assert np.abs(Nlow @ N @ Nlow - Nlow).max() < 1e-9 * np.abs(Nlow).max()
    d, mean = rng.standard_normal(n), rng.standard_normal(n)
    dlow = Yk @ (Ypk @ d)
    q1 = (dlow - mean) @ Nlow @ (dlow - mean)
    c = Ypk @ d - Ypk @ mean
    q2 = c @ np.linalg.solve(C, c)
    assert np.isclose(q1, q2, rtol=1e-10)
    # the pseudo-inverse of the kept columns alone is not blind to the columns above the cut
    assert np.abs(np.linalg.pinv(Yk) @ Y[:, k:]).max() > 1e-3


def test_masking_is_the_schur_complement_and_matches_the_kept_block():
    from lowres.mask import marginalize, restrict_prior
    rng = np.random.default_rng(31)
    n = 40
    X = rng.standard_normal((n, n))
    N = X @ X.T / n + np.eye(n)
    keep = np.ones(n, dtype=bool)
    keep[[2, 5, 17, 30, 31]] = False
    d = rng.standard_normal(n)
    dk, nk = marginalize(d, np.linalg.inv(N), keep)
    assert np.allclose(dk, d[keep]) and np.allclose(nk, np.linalg.inv(N[np.ix_(keep, keep)]), rtol=1e-9, atol=1e-12)
    assert np.allclose(restrict_prior(N, keep), nk, rtol=1e-9, atol=1e-12)
    # a singular inverse covariance (a band known only through 8 linear combinations): the pseudo-inverse form leaves
    # a positive semidefinite matrix whose range lies in the kept elements' part of those combinations
    L = rng.standard_normal((8, n))
    low = L.T @ np.linalg.inv(L @ N @ L.T) @ L
    _, lk = marginalize(d, low, keep, pinv=True)
    w = np.linalg.eigvalsh(lk)
    assert w.min() > -1e-9 * w.max() and (w > 1e-9 * w.max()).sum() <= 8


def test_own_score_block_error_ignores_a_misfitting_neighbour():
    """Two regions coupled through dense noise; region 1 is given a strong index spread between blocks. The joint
    block error of region 0 inherits it through the off-diagonal of J_c^-1, the own-score block error does not."""
    n = 192
    region, a_s, a_d, specs = setup(n, 2, seed=11, noise=0.6, dense=(0, 1, 2, 3, 4))
    block = np.arange(n) // 16                                    # 12 blocks, each holding both regions
    rng = np.random.default_rng(12)
    joint, own, mod = [], [], []
    for _ in range(40):
        off = 0.5 * rng.standard_normal(12)[block] * (region == 1)
        L = RegionLikelihood(draw(specs, a_s, a_d, 0.2 * rng.standard_normal(n), -3.1 + off, rng), region, cmb_prior_sd=0.2)
        bh, _ = fit_profile(L, np.full(2, -3.1))
        out = block_errors(L, bh, block, 12)
        joint.append(out["CR3"]); own.append(out["CR3_own"]); mod.append(out["sd_model"])
    joint, own, mod = np.median(joint, axis=0), np.median(own, axis=0), np.median(mod, axis=0)
    assert own[1] > 3 * mod[1]                                    # the region with the spread is far from its model error
    assert own[0] <= joint[0] + 1e-12                             # and its neighbour's own-score error is the smaller one
    assert own[0] < 2.5 * mod[0], (own, joint, mod)


def test_cmb_prior_mismatch_predicts_bias_and_scatter():
    """A white prior fitted to skies whose CMB has a dense covariance with strong modes: the first-order bias and
    covariance of `lowres.mismatch.cmb_prior_mismatch` against 500 simulated fits; and no change when the prior
    is the truth."""
    from lowres.mismatch import cmb_prior_mismatch
    n, sd = 96, 0.3
    region, a_s, a_d, specs = setup(n, 2, seed=13, noise=1.0, amp=6.0, dense=(0,))
    rng = np.random.default_rng(14)
    X = rng.standard_normal((n, 5))
    S_true = 0.6 * X @ X.T / 5 + 0.02 * np.eye(n)                 # five strong modes, total variance well above the prior's in them
    root = np.linalg.cholesky(S_true)
    truth = np.array([-3.1, -3.0])
    est = []
    for _ in range(500):
        L = RegionLikelihood(draw(specs, a_s, a_d, root @ rng.standard_normal(n), truth[region], rng), region, cmb_prior_sd=sd)
        est.append(fit_profile(L, truth)[0])
    est = np.array(est)
    # prediction at the truth, from a noise-free data set with no CMB (so that the fitted amplitudes are the true ones)
    clean = [BandData(B.name, B.power_law(truth[region], 23.0)[0] * a_s + dust * a_d, ninv, (lambda b, B=B: B.power_law(b, 23.0)), dust)
             for B, dust, ninv, _ in specs]
    Lc = RegionLikelihood(clean, region, cmb_prior_sd=sd)
    out = cmb_prior_mismatch(Lc, truth, 1.0 / sd ** 2, S_true - sd ** 2 * np.eye(n))
    sd_true, sd_matched = np.sqrt(np.diag(out["cov_true"])), np.sqrt(np.diag(out["cov_matched"]))
    emp = est.std(axis=0, ddof=1)
    assert np.max(sd_true / sd_matched) > 1.2, (sd_true, sd_matched)            # the configuration makes the effect visible
    assert np.all(np.abs(emp / sd_true - 1.0) < 0.15), (emp, sd_true, sd_matched)
    assert np.all(np.abs(est.mean(axis=0) - truth - out["bias"]) < 4 * emp / np.sqrt(len(est)) + 0.15 * np.abs(out["bias"])), (est.mean(0) - truth, out["bias"])
    same = cmb_prior_mismatch(Lc, truth, 1.0 / sd ** 2, np.zeros((n, n)))
    assert np.allclose(same["bias"], 0.0) and np.allclose(same["cov_true"], np.linalg.inv(same_J(Lc, truth)), rtol=0.2)


def same_J(like, beta):
    J, T, _ = like.information(beta)
    return J
