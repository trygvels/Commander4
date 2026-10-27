"""Direct tests of the compiled per-pixel component-separation solver in ``commander4.backend``."""
import numpy as np
import pytest
from numpy.typing import NDArray

from commander4.backend import compsep as cpp_compsep

NBAND, NCOMP, NPIX = 5, 3, 40


def _inputs(seed: int) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Random sky and inverse variance (nband, npix), mixing matrix M (nband, ncomp), and eta
    (npix, nband)."""
    rng = np.random.default_rng(seed)
    return (rng.normal(size=(NBAND, NPIX)), rng.uniform(0.25, 4.0, (NBAND, NPIX)),
            rng.uniform(0.5, 2.0, (NBAND, NCOMP)), rng.normal(size=(NPIX, NBAND)))


@pytest.mark.parametrize("nthreads", [1, 4])
def test_solve_perpix_matches_numpy(nthreads: int) -> None:
    """Each pixel solves (M^T N^-1 M) x = M^T (N^-1 d + N^-1/2 eta), N^-1 the inverse variance."""
    map_sky, map_inv_var, M, rand = _inputs(seed=1)
    w = map_inv_var                                                    # (nband, npix)
    A = np.einsum("bi,bp,bj->pij", M, w, M)                            # (npix, ncomp, ncomp)
    rhs = np.einsum("bi,bp->pi", M, w*map_sky + np.sqrt(w)*rand.T)     # (npix, ncomp)
    expected = np.linalg.solve(A, rhs[..., None])[..., 0].T            # (ncomp, npix)
    comp_maps = cpp_compsep.solve_perpix(map_sky, map_inv_var, M, rand, nthreads)
    np.testing.assert_allclose(comp_maps, expected, rtol=1e-10, atol=1e-12)


def test_solve_perpix_zeroes_unsolvable_pixels() -> None:
    """A pixel with zero inverse variance in every band carries no information: amplitudes 0."""
    map_sky, map_inv_var, M, rand = _inputs(seed=2)
    map_inv_var[:, 7] = 0.0
    comp_maps = cpp_compsep.solve_perpix(map_sky, map_inv_var, M, rand)
    np.testing.assert_array_equal(comp_maps[:, 7], 0.0)
    assert np.all(comp_maps[:, :7] != 0.0)


def test_solve_perpix_rejects_wrong_shapes() -> None:
    map_sky, map_inv_var, M, rand = _inputs(seed=3)
    with pytest.raises(RuntimeError):  # rand must be (npix, nband), not (nband, npix)
        cpp_compsep.solve_perpix(map_sky, map_inv_var, M, np.ascontiguousarray(rand.T))
    with pytest.raises(RuntimeError):  # M must have one row per band
        cpp_compsep.solve_perpix(map_sky, map_inv_var, M[:4].copy(), rand)


# --- per-region solver and the extras output of the region likelihood (M2.4 stage 1) ---------

NREG = 4


def _region_inputs(seed: int) -> tuple[NDArray, NDArray, NDArray, NDArray, NDArray]:
    """Random sky and inverse variance (nband, npix), region map (npix,), per-region mixing
    matrices (nreg, nband, ncomp) and eta (npix, nband)."""
    rng = np.random.default_rng(seed)
    region = rng.integers(0, NREG, NPIX).astype(np.int64)
    return (rng.normal(size=(NBAND, NPIX)), rng.uniform(0.25, 4.0, (NBAND, NPIX)), region,
            rng.uniform(0.5, 2.0, (NREG, NBAND, NCOMP)), rng.normal(size=(NPIX, NBAND)))


def _normal_eqs(map_sky: NDArray, w: NDArray, M_pix: NDArray) -> tuple[NDArray, NDArray]:
    """A_p = M_p^T N^-1 M_p (npix, ncomp, ncomp) and b_p = M_p^T N^-1 d_p (npix, ncomp) for a
    per-pixel mixing matrix M_pix (npix, nband, ncomp)."""
    A = np.einsum("pbi,bp,pbj->pij", M_pix, w, M_pix)
    b = np.einsum("pbi,bp->pi", M_pix, w * map_sky)
    return A, b


@pytest.mark.parametrize("nthreads", [1, 3])
def test_solve_perpix_regions_matches_numpy(nthreads: int) -> None:
    """Each pixel solves (M_r^T N^-1 M_r) x = M_r^T (N^-1 d + N^-1/2 eta) with its region's M_r."""
    map_sky, w, region, M_reg, rand = _region_inputs(seed=11)
    A = np.einsum("pbi,bp,pbj->pij", M_reg[region], w, M_reg[region])
    rhs = np.einsum("pbi,bp->pi", M_reg[region], w * map_sky + np.sqrt(w) * rand.T)
    expected = np.linalg.solve(A, rhs[..., None])[..., 0].T
    got = cpp_compsep.solve_perpix_regions(map_sky, w, region, M_reg, rand, nthreads)
    np.testing.assert_allclose(got, expected, rtol=1e-10, atol=1e-12)


def test_solve_perpix_regions_masked_pixels_are_zero() -> None:
    map_sky, w, region, M_reg, rand = _region_inputs(seed=12)
    region[[3, 17]] = -1
    got = cpp_compsep.solve_perpix_regions(map_sky, w, region, M_reg, rand)
    np.testing.assert_array_equal(got[:, [3, 17]], 0.0)
    assert np.all(got[:, region >= 0] != 0.0)


def test_solve_perpix_regions_single_region_is_bitwise_solve_perpix() -> None:
    """One region reproduces the scalar kernel exactly (shared body, same operation order)."""
    map_sky, map_inv_var, M, rand = _inputs(seed=13)
    map_inv_var[:, 5] = 0.0
    scalar = cpp_compsep.solve_perpix(map_sky, map_inv_var, M, rand, 2)
    regions = cpp_compsep.solve_perpix_regions(map_sky, map_inv_var,
                                               np.zeros(NPIX, dtype=np.int64), M[None], rand, 2)
    np.testing.assert_array_equal(regions, scalar)


def test_solve_perpix_regions_rejects_bad_region_index() -> None:
    map_sky, w, region, M_reg, rand = _region_inputs(seed=14)
    region[0] = NREG
    with pytest.raises(RuntimeError):
        cpp_compsep.solve_perpix_regions(map_sky, w, region, M_reg, rand)


def test_marginal_loglike_regions_extras_match_numpy() -> None:
    """want_extras returns x = A^-1 b, diag(A^-1) and the used mask; pixels with region -1 or
    fewer than ncomp observed bands are unused and zero; the 3-tuple is unchanged without it."""
    map_sky, w, region, M_reg, _ = _region_inputs(seed=15)
    region[2] = -1
    w[:, 9] = 0.0            # unobserved
    w[:NBAND - NCOMP + 1, 20] = 0.0   # ncomp - 1 observed bands: under-determined
    dM = np.zeros_like(M_reg)
    plain = cpp_compsep.marginal_loglike_regions(map_sky, w, region, M_reg, dM, 1)
    full = cpp_compsep.marginal_loglike_regions(map_sky, w, region, M_reg, dM, 1,
                                                want_extras=True)
    assert len(plain) == 3 and len(full) == 6
    np.testing.assert_array_equal(full[0], plain[0])
    assert full[2] == plain[2] == 0
    x, diag, used = full[3:]
    assert x.shape == diag.shape == (NCOMP, NPIX) and used.shape == (NPIX,)
    assert used.dtype == np.bool_
    expected_used = region >= 0
    expected_used[[9, 20]] = False
    np.testing.assert_array_equal(used, expected_used)
    A, b = _normal_eqs(map_sky, w, M_reg[region])
    Ainv = np.linalg.inv(A[used])
    np.testing.assert_allclose(x[:, used], np.einsum("pij,pj->ip", Ainv, b[used]),
                               rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(diag[:, used], np.einsum("pii->ip", Ainv), rtol=1e-10, atol=1e-12)
    np.testing.assert_array_equal(x[:, ~used], 0.0)
    np.testing.assert_array_equal(diag[:, ~used], 0.0)


def _augment(map_sky: NDArray, w: NDArray, M_reg: NDArray, m: NDArray, t: NDArray
             ) -> tuple[NDArray, NDArray, NDArray]:
    """Append one pseudo-band per component: map m_c (ncomp, npix), inverse variance 1/t_c^2
    (0 for a flat component, t = inf), mixing row e_c in every region."""
    with np.errstate(divide="ignore"):
        inv_t2 = np.where(np.isfinite(t), 1.0 / t**2, 0.0)
    sky_aug = np.concatenate([map_sky, m], axis=0)
    w_aug = np.concatenate([w, inv_t2], axis=0)
    rows = np.broadcast_to(np.eye(NCOMP), (M_reg.shape[0], NCOMP, NCOMP))
    return sky_aug, w_aug, np.ascontiguousarray(np.concatenate([M_reg, rows], axis=1))


def test_pseudo_band_augmentation_gives_gaussian_prior_marginal() -> None:
    """With the prior a_c ~ N(m_c, t_c^2) as pseudo-bands the kernel's (A, b) are (A', b') =
    (A + T^-1, b + T^-1 m): the per-region sum of 1/2 b'^T A'^-1 b' - 1/2 log det A' to 1e-12,
    including a flat component (t = inf) that contributes nothing."""
    map_sky, w, region, M_reg, _ = _region_inputs(seed=16)
    rng = np.random.default_rng(17)
    m = rng.normal(size=(NCOMP, NPIX))
    t = np.broadcast_to(np.array([0.7, np.inf, 2.5])[:, None], (NCOMP, NPIX)).copy()
    sky_aug, w_aug, M_aug = _augment(map_sky, w, M_reg, m, t)
    ll, _, nskip = cpp_compsep.marginal_loglike_regions(sky_aug, w_aug, region, M_aug,
                                                        np.zeros_like(M_aug), 1)
    assert nskip == 0
    A, b = _normal_eqs(map_sky, w, M_reg[region])
    inv_t2 = w_aug[NBAND:].T                                  # (npix, ncomp)
    A_p = A + np.einsum("pi,ij->pij", inv_t2, np.eye(NCOMP))
    b_p = b + inv_t2 * m.T
    ll_pix = 0.5 * np.einsum("pi,pi->p", b_p, np.linalg.solve(A_p, b_p[..., None])[..., 0]) \
        - 0.5 * np.linalg.slogdet(A_p)[1]
    expected = np.bincount(region, weights=ll_pix, minlength=NREG)
    np.testing.assert_allclose(ll, expected, rtol=1e-12, atol=1e-10)


def test_pseudo_band_draw_has_prior_posterior_moments() -> None:
    """solve_perpix_regions on the augmented stream draws a_p ~ N(A'^-1 b', A'^-1): sample mean
    and covariance of 20000 draws in two pixels (one component flat, two with proper priors)
    within 4 standard errors. The draws are one kernel call on the pixels repeated 20000 times."""
    map_sky, w, region, M_reg, _ = _region_inputs(seed=18)
    rng = np.random.default_rng(19)
    m = rng.normal(size=(NCOMP, NPIX))
    t = np.broadcast_to(np.array([0.5, 1.5, np.inf])[:, None], (NCOMP, NPIX)).copy()
    sky_aug, w_aug, M_aug = _augment(map_sky, w, M_reg, m, t)
    A, b = _normal_eqs(map_sky, w, M_reg[region])
    inv_t2 = w_aug[NBAND:].T
    A_p = A + np.einsum("pi,ij->pij", inv_t2, np.eye(NCOMP))
    cov = np.linalg.inv(A_p)
    mean = np.einsum("pij,pj->pi", cov, b + inv_t2 * m.T)
    ndraw, pixels = 20000, [0, NPIX - 1]
    rep = np.repeat(pixels, ndraw)                                   # (2 ndraw,)
    eta = rng.normal(size=(rep.size, NBAND + NCOMP))
    draws = cpp_compsep.solve_perpix_regions(np.ascontiguousarray(sky_aug[:, rep]),
                                             np.ascontiguousarray(w_aug[:, rep]), region[rep],
                                             M_aug, eta, 2)          # (ncomp, 2 ndraw)
    for k, p in enumerate(pixels):
        d = draws[:, k * ndraw:(k + 1) * ndraw]
        sd = np.sqrt(np.diag(cov[p]))
        assert np.all(np.abs(d.mean(axis=1) - mean[p]) <= 4 * sd / np.sqrt(ndraw))
        se_cov = np.sqrt((np.outer(sd, sd)**2 + cov[p]**2) / (ndraw - 1))
        assert np.all(np.abs(np.cov(d) - cov[p]) <= 4 * se_cov)
