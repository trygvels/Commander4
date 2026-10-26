"""Parity of commander4.math_utils.gmrf_field with the E012 prototype
(nextgen/prototypes/indexfield/field.py and posterior.py, loaded by file path), basis
properties, the centred/non-centred Jacobian identity and finite-difference gradients."""
import importlib
import sys
from pathlib import Path

import numpy as np
import pytest

from commander4.math_utils.gmrf_field import (NSIDE_BETA_MAX, GraphLaplacianBasis, HyperPrior,
                                              IndexFieldPrior, MaternField, RhoTransform)

_PROTO_DIR = Path(__file__).resolve().parents[1] / "nextgen" / "prototypes"


@pytest.fixture(scope="module")
def proto():
    """The prototype package (field, posterior) imported as `indexfield` from nextgen/prototypes."""
    if not (_PROTO_DIR / "indexfield" / "field.py").is_file():
        pytest.skip("prototype nextgen/prototypes/indexfield not present")
    if str(_PROTO_DIR) not in sys.path:
        sys.path.insert(0, str(_PROTO_DIR))
    return importlib.import_module("indexfield.field"), importlib.import_module(
        "indexfield.posterior")


NSIDES = (2, 4)
MU0, S_MU, SIGMA_SCALE = -3.0, 0.5, 0.3
# (mu, sigma, u) triples; rho follows from u through the transform (rho_min = h, rho_max = pi).
THETAS = [(-3.1, 0.05, -3.0), (-2.8, 0.3, -0.5), (-3.0, 1.0, 0.0), (-2.9, 0.12, 1.5),
          (-3.3, 0.6, 4.0)]


def _build(nside):
    basis = GraphLaplacianBasis(nside)
    field = MaternField(basis, alpha=2.0)
    hyper = HyperPrior(mu0=MU0, s_mu=S_MU, sigma_scale=SIGMA_SCALE)
    return basis, field, hyper


def _build_proto(proto, nside):
    pf, pp = proto
    basis = pf.GraphLaplacianBasis(nside)
    field = pf.MaternField(basis, alpha=2.0)
    hyper = pp.HyperPrior(mu0=MU0, s_mu=S_MU, sigma_scale=SIGMA_SCALE, rho_prior="loguniform")
    return basis, field, hyper


def _gauss_loglike(nb, rng):
    """Gaussian stand-in for the collapsed likelihood, not part of the package."""
    prec = 10.0 ** rng.uniform(1.0, 3.0, size=nb)
    b0 = -3.0 + 0.2 * rng.normal(size=nb)

    def loglike(beta):
        r = beta - b0
        return -0.5 * float(np.sum(prec * r * r)), -prec * r
    return loglike


def _fd4_grad(f, x, h):
    """Fourth-order central-difference gradient of a scalar function of a vector."""
    g = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h
        g[i] = (8.0 * (f(x + e) - f(x - e)) - (f(x + 2 * e) - f(x - 2 * e))) / (12.0 * h)
    return g


def _assert_fd(grad, g_fd):
    assert np.all(np.abs(g_fd - grad) <= 1e-6 * np.abs(grad) + 1e-8)


# ------------------------------------------------------------------------------- basis properties
@pytest.mark.parametrize("nside", NSIDES)
def test_basis_orthonormal_and_constant_mode_absent(nside):
    basis = _build(nside)[0]
    assert basis.K == basis.n_beta - 1
    assert basis.orthonormality_error() < 1e-14
    assert np.max(np.abs(basis.V.T @ np.ones(basis.n_beta))) < 1e-12        # S^T 1 = 0
    assert basis.h == pytest.approx(np.sqrt(4 * np.pi / basis.n_beta), rel=0, abs=0)
    assert np.allclose(basis.lam, basis.graph_eigenvalues / (3 * basis.h ** 2), rtol=1e-14,
                       atol=0)
    assert np.all(basis.lam > 0)


def test_nside_cap_and_validation():
    with pytest.raises(ValueError, match="NSIDE_BETA_MAX"):
        GraphLaplacianBasis(2 * NSIDE_BETA_MAX)
    with pytest.raises(ValueError):
        GraphLaplacianBasis(3)
    with pytest.raises(ValueError):
        GraphLaplacianBasis(0)


def test_matern_normalization():
    basis, field, _ = _build(2)
    for rho in (basis.h, 0.3, np.pi):
        assert abs(field.pixel_variances(rho).mean() - 1.0) < 1e-12     # trace(S S^T) = Nb
        z = np.random.default_rng(1).normal(size=basis.K)
        q = field.Kinv_quadform(rho, field.S_times(rho, z))
        assert abs(q - z @ z) < 1e-10 * (z @ z)


# ------------------------------------------------------------------------- parity with prototype
@pytest.mark.parametrize("nside", NSIDES)
def test_field_operators_match_prototype(proto, nside):
    basis, field, hyper = _build(nside)
    pbasis, pfield, _ = _build_proto(proto, nside)
    rng = np.random.default_rng(nside)
    assert np.array_equal(basis.V, pbasis.V)
    assert np.array_equal(basis.lam, pbasis.lam)
    assert basis.h == pbasis.h
    tr = RhoTransform(basis.h, np.pi)
    for _, _, u in THETAS:
        rho = float(tr.rho(u))
        z = rng.normal(size=basis.K)
        v = rng.normal(size=basis.n_beta)
        assert np.max(np.abs(field.S_times(rho, z) - pfield.S_times(rho, z))) < 1e-12
        assert np.max(np.abs(field.ST_times(rho, v) - pfield.ST_times(rho, v))) < 1e-12
        assert np.max(np.abs(field.dS_dlogrho_times(rho, z)
                             - pfield.dS_dlogrho_times(rho, z))) < 1e-12
        assert abs(field.logdet_K(rho) - pfield.logdet_K(rho)) < 1e-12
        w, dw = field.weights(rho)
        pw, pdw = pfield.weights(rho)
        assert np.max(np.abs(w - pw)) < 1e-12 and np.max(np.abs(dw - pdw)) < 1e-12
        # <v, S z> = <S^T v, z>
        assert abs(v @ field.S_times(rho, z) - field.ST_times(rho, v) @ z) < 1e-9


@pytest.mark.parametrize("nside", NSIDES)
def test_rho_transform_matches_prototype(proto, nside):
    pf = proto[0]
    h = GraphLaplacianBasis(nside).h
    tr, ptr = RhoTransform(h, np.pi), pf.RhoTransform(h, np.pi)
    assert tr.L == ptr.L and tr.log_rho_min == ptr.log_rho_min
    assert tr.L == pytest.approx(np.log(np.pi / h), rel=0, abs=1e-15)
    for _, _, u in THETAS:
        assert abs(tr.rho(u) - ptr.rho(u)) < 1e-12
        assert abs(tr.log_jacobian(u) - ptr.log_jacobian(u)) < 1e-12
        assert abs(tr.dlogrho_du(u) - ptr.dlogrho_du(u)) < 1e-12
        assert abs(tr.dlog_jacobian_du(u) - ptr.dlog_jacobian_du(u)) < 1e-12
        s = 1 / (1 + np.exp(-u))
        assert abs(tr.log_jacobian(u) - (np.log(tr.L) + np.log(s) + np.log1p(-s))) < 1e-12
        assert abs(tr.u(tr.rho(u)) - u) < 1e-9


@pytest.mark.parametrize("nside", NSIDES)
def test_hyper_terms_match_prototype(proto, nside):
    _, pp = proto
    basis, field, hyper = _build(nside)
    pbasis, pfield, phyper = _build_proto(proto, nside)
    prior = IndexFieldPrior(field, hyper)
    ppost = pp.NonCentredPosterior(pfield, phyper, loglike=lambda b: (0.0, np.zeros_like(b)))
    assert prior.hyper.rho_min == ppost.hyper.rho_min == basis.h
    for mu, sigma, u in THETAS:
        ours = prior.hyper_terms(mu, np.log(sigma), u)
        theirs = ppost._hyper_terms(mu, np.log(sigma), u)
        assert np.max(np.abs(np.array(ours) - np.array(theirs))) < 1e-12


@pytest.mark.parametrize("nside", NSIDES)
@pytest.mark.parametrize("centred", [False, True])
def test_logp_and_grad_match_prototype(proto, nside, centred):
    _, pp = proto
    _, field, hyper = _build(nside)
    _, pfield, phyper = _build_proto(proto, nside)
    rng = np.random.default_rng(10 * nside + centred)
    loglike = _gauss_loglike(field.n_beta, rng)
    prior = IndexFieldPrior(field, hyper, centred=centred)
    Post = pp.CentredPosterior if centred else pp.NonCentredPosterior
    ppost = Post(pfield, phyper, loglike=loglike)
    for mu, sigma, u in THETAS:
        block = rng.normal(size=field.K)
        if centred:
            block *= sigma * field.weights(float(prior.transform.rho(u)))[0]
        x = prior.pack(block, mu, np.log(sigma), u)
        lp, grad = prior.logp_and_grad(x, loglike)
        plp, pgrad = ppost.logp_and_grad(x)
        assert abs(lp - plp) <= 1e-12 * max(1.0, abs(plp))
        assert np.max(np.abs(grad - pgrad)) <= 1e-12 * max(1.0, np.max(np.abs(pgrad)))
        assert np.max(np.abs(prior.beta(x) - ppost.beta(x))) < 1e-12
        # Prior-only part equals the prototype with a zero likelihood.
        lp0, grad0 = prior.prior_logp_and_grad(x)
        zero = Post(pfield, phyper, loglike=lambda b: (0.0, np.zeros_like(b)))
        plp0, pgrad0 = zero.logp_and_grad(x)
        assert abs(lp0 - plp0) <= 1e-12 * max(1.0, abs(plp0))
        assert np.max(np.abs(grad0 - pgrad0)) <= 1e-12 * max(1.0, np.max(np.abs(pgrad0)))


# ------------------------------------------------------------------------ Jacobian identity
@pytest.mark.parametrize("nside", NSIDES)
def test_centred_minus_noncentred_is_jacobian(nside):
    _, field, hyper = _build(nside)
    rng = np.random.default_rng(7)
    loglike = _gauss_loglike(field.n_beta, rng)
    nc = IndexFieldPrior(field, hyper, centred=False)
    c = IndexFieldPrior(field, hyper, centred=True)
    for mu, sigma, u in THETAS:
        x_nc = nc.pack(rng.normal(size=field.K), mu, np.log(sigma), u)
        x_c = c.from_noncentred(x_nc)
        assert np.max(np.abs(c.to_noncentred(x_c) - x_nc)) < 1e-12
        assert np.max(np.abs(nc.beta(x_nc) - c.beta(x_c))) < 1e-12
        rho = float(nc.transform.rho(u))
        expected = -field.K * np.log(sigma) - 0.5 * field.logdet_K(rho)
        assert abs(c.centred_log_jacobian(x_c) - expected) < 1e-12
        for ll in (None, loglike):
            lp_nc = nc.logp_and_grad(x_nc, ll)[0]
            lp_c = c.logp_and_grad(x_c, ll)[0]
            assert abs((lp_c - lp_nc) - expected) <= 1e-10 * max(1.0, abs(expected))


# ------------------------------------------------------------------- finite-difference gradients
@pytest.mark.parametrize("nside", NSIDES)
def test_spectrum_derivatives_fd(nside):
    _, field, _ = _build(nside)
    rng = np.random.default_rng(3)
    z = rng.normal(size=field.K)
    for rho in (field.basis.h, 0.3, 2.0):
        lr, h = np.log(rho), 1e-3

        def fd(f):
            return (8.0 * (f(lr + h) - f(lr - h)) - (f(lr + 2 * h) - f(lr - 2 * h))) / (12.0 * h)

        s, c, ds, dc = field.spectrum(rho)
        _assert_fd(ds, fd(lambda t: field.spectrum(np.exp(t))[0]))
        _assert_fd(np.array([dc]), np.array([fd(lambda t: field.spectrum(np.exp(t))[1])]))
        _assert_fd(field.weights(rho)[1], fd(lambda t: field.weights(np.exp(t))[0]))
        _assert_fd(field.dS_dlogrho_times(rho, z), fd(lambda t: field.S_times(np.exp(t), z)))
        _assert_fd(np.array([2 * np.sum(field.weights(rho)[1] / field.weights(rho)[0])]),
                   np.array([fd(lambda t: field.logdet_K(np.exp(t)))]))


def test_rho_transform_fd():
    tr = RhoTransform(GraphLaplacianBasis(2).h, np.pi)
    for u in (-6.0, -1.0, 0.0, 2.5, 7.0):
        t = np.array([u])
        _assert_fd(np.array([tr.dlogrho_du(u)]), _fd4_grad(lambda v: np.log(tr.rho(v[0])), t, 1e-3))
        _assert_fd(np.array([tr.dlog_jacobian_du(u)]),
                   _fd4_grad(lambda v: tr.log_jacobian(v[0]), t, 1e-3))


@pytest.mark.parametrize("nside", NSIDES)
def test_hyper_terms_fd(nside):
    _, field, hyper = _build(nside)
    prior = IndexFieldPrior(field, hyper)
    for mu, sigma, u in THETAS:
        t = np.array([mu, np.log(sigma), u])
        _, *grad = prior.hyper_terms(*t)
        _assert_fd(np.array(grad), _fd4_grad(lambda v: prior.hyper_terms(*v)[0], t, 1e-4))


@pytest.mark.parametrize("nside", NSIDES)
@pytest.mark.parametrize("centred", [False, True])
def test_logp_and_grad_fd(nside, centred):
    _, field, hyper = _build(nside)
    rng = np.random.default_rng(100 + nside + centred)
    loglike = _gauss_loglike(field.n_beta, rng)
    prior = IndexFieldPrior(field, hyper, centred=centred)
    for mu, sigma, u in THETAS[:3]:
        block = rng.normal(size=field.K)
        if centred:
            block *= sigma * field.weights(float(prior.transform.rho(u)))[0]
        x = prior.pack(block, mu, np.log(sigma), u)
        for ll in (None, loglike):
            lp, grad = prior.logp_and_grad(x, ll)
            assert np.isfinite(lp)
            _assert_fd(grad, _fd4_grad(lambda v: prior.logp_and_grad(v, ll)[0], x, 1e-3))


def test_zero_density_contract():
    _, field, hyper = _build(2)
    prior = IndexFieldPrior(field, hyper)
    x = prior.pack(np.zeros(field.K), MU0, np.log(0.3), 0.0)

    def failing(beta):
        raise np.linalg.LinAlgError("not positive definite")

    for ll in (failing, lambda b: (np.nan, np.zeros_like(b))):
        lp, grad = prior.logp_and_grad(x, ll)
        assert lp == -np.inf and np.all(grad == 0.0) and grad.shape == (prior.dim,)
    with pytest.raises(ValueError):
        IndexFieldPrior(field, HyperPrior(mu0=MU0, s_mu=S_MU, rho_prior="pc"))
