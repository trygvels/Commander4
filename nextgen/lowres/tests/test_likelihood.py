"""Exact-target checks of the dense region likelihood (L0/L1). Run: pytest -q nextgen/lowres/tests"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from lowres.fit import fit_profile, sandwich, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band, g_rj, mbb_rj  # noqa: E402

NUS = [22.8, 28.4, 33.0, 40.7, 60.8, 353.0]
NU_REF = 23.0


def make(n=64, nreg=4, dense=(0, 2), seed=0, noise=1.0, beta_true=-3.1, amp=8.0, with_dust=True):
    rng = np.random.default_rng(seed)
    region = np.arange(n) % nreg
    bands, cov = [], []
    a_s = amp * rng.standard_normal(n) * (1.0 + 3.0 * rng.random(n))
    a_d = (0.5 * rng.standard_normal(n)) if with_dust else np.zeros(n)
    a_c = 0.2 * rng.standard_normal(n)
    for i, nu in enumerate(NUS):
        B = Band.delta(f"b{i}", nu)
        dust = B.mbb(1.55, 19.6, 353.0)
        sig = noise * (0.2 if nu > 300 else rng.uniform(1.0, 2.5))
        if i in dense:
            X = rng.standard_normal((n, n))
            Cv = sig ** 2 * (np.eye(n) + 0.5 * X @ X.T / n)
        else:
            Cv = np.diag(sig ** 2 * rng.uniform(0.6, 1.6, n))
        cov.append(Cv)
        s0, _ = B.power_law(np.full(nreg, beta_true), NU_REF)
        mean = s0[region] * a_s + dust * a_d + a_c
        d = mean + np.linalg.cholesky(Cv) @ rng.standard_normal(n)
        ninv = np.linalg.inv(Cv) if i in dense else 1.0 / np.diag(Cv)
        bands.append((B, dust, d, ninv, mean))
    def build(noise_draw=None, clean=False):
        out = []
        for k, (B, dust, d, ninv, mean) in enumerate(bands):
            if clean:
                dd = mean
            else:
                dd = d if noise_draw is None else mean + np.linalg.cholesky(cov[k]) @ noise_draw.standard_normal(n)
            out.append(BandData(B.name, dd, ninv, lambda b, B=B: B.power_law(b, NU_REF), dust))
        return out
    return build, region, a_s


def fd(f, x, h=1e-5):
    g = np.zeros_like(x)
    for i in range(x.size):
        e = np.zeros_like(x); e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g


def test_gradients():
    build, region, _ = make()
    beta = np.array([-3.0, -3.2, -2.9, -3.15])
    for kw in ({}, {"cmb_prior_sd": 0.3}, {"use_cmb": False}, {"use_dust": False, "use_cmb": False}):
        L = RegionLikelihood(build(), region, **kw)
        f, g = L.profile(beta, grad=True)
        assert np.allclose(g, fd(L.profile, beta), rtol=2e-5, atol=1e-6), (kw, g, fd(L.profile, beta))
        f, g = L.marginal(beta, grad=True)
        assert np.allclose(g, fd(L.marginal, beta), rtol=2e-5, atol=1e-6), kw
    L = RegionLikelihood(build(), region, cmb_prior_sd=0.3, sync_prior=True)
    lt = np.log(np.array([8.0, 12.0, 20.0, 5.0]))
    f, gb, gt = L.marginal(beta, lt, grad=True)
    assert np.allclose(gb, fd(lambda b: L.marginal(b, lt), beta), rtol=2e-5, atol=1e-6)
    assert np.allclose(gt, fd(lambda t: L.marginal(beta, t), lt), rtol=2e-5, atol=1e-6)


def test_diagonal_noise_matches_per_element_formulas():
    """Synchrotron only, diagonal noise: the dense code must equal the per-element closed forms."""
    build, region, _ = make(dense=(), with_dust=False)
    bands = build()
    L = RegionLikelihood(bands, region, use_dust=False, use_cmb=False)
    beta = np.array([-3.05, -3.2, -2.95, -3.1])
    S = sum(b.ninv * b.d ** 2 for b in bands)
    s = np.array([b.sync(beta)[0][region] for b in bands])
    bb = sum(s[i] * b.ninv * b.d for i, b in enumerate(bands))
    A = sum(s[i] ** 2 * b.ninv for i, b in enumerate(bands))
    assert np.isclose(L.profile(beta), 0.5 * np.sum(S - bb ** 2 / A), rtol=1e-12)
    assert np.isclose(L.marginal(beta), 0.5 * np.sum(S - bb ** 2 / A) + 0.5 * np.sum(np.log(A)), rtol=1e-12)
    # zero-mean prior
    Lz = RegionLikelihood(bands, region, use_dust=False, use_cmb=False, sync_prior=True)
    lt = np.log(np.array([6.0, 9.0, 14.0, 4.0])); t2 = np.exp(2 * lt)[region]
    ref = 0.5 * np.sum(S - bb ** 2 * t2 / (1 + t2 * A)) + 0.5 * np.sum(np.log1p(t2 * A))
    # the dense objective carries ln det(A + 1/tau^2) + ln tau^2 = ln(1 + tau^2 A): identical
    assert np.isclose(Lz.marginal(beta, lt), ref, rtol=1e-11)


def test_information_matches_single_component_identity():
    """One component, diagonal noise: J = sum u V and T = sum V per region (u = a^2 A)."""
    build, region, _ = make(dense=(), with_dust=False)
    bands = build()
    L = RegionLikelihood(bands, region, use_dust=False, use_cmb=False)
    beta = np.full(4, -3.1)
    J, T, sol = L.information(beta)
    s = np.array([b.sync(beta)[0][region] for b in bands]); ds = np.array([b.sync(beta)[1][region] for b in bands])
    w = np.array([b.ninv for b in bands])
    A = (w * s ** 2).sum(0)
    V = ((w * ds ** 2).sum(0) - (w * s * ds).sum(0) ** 2 / A) / A
    a = sol["x"][:L.n]
    assert np.allclose(np.diag(J), np.bincount(region, weights=a ** 2 * A * V), rtol=1e-9)
    assert np.allclose(np.diag(T), np.bincount(region, weights=V), rtol=1e-9)
    assert np.allclose(J - np.diag(np.diag(J)), 0, atol=1e-9 * np.abs(J).max())


def test_sandwich_matches_monte_carlo():
    """Dense noise, three components, low S/N: the scatter of the profile optimum over noise draws equals
    the sandwich covariance at the truth, and exceeds the naive J^-1."""
    build, region, a_s = make(n=96, nreg=3, dense=(0, 1, 2), noise=3.0, amp=6.0, seed=3)
    rng = np.random.default_rng(11)
    est = []
    for _ in range(300):
        L = RegionLikelihood(build(rng), region, cmb_prior_sd=0.3)
        bh, res = fit_profile(L, np.full(3, -3.1))
        est.append(bh)
    est = np.array(est)
    # information at the truth: a noiseless data vector, so that the fitted amplitudes are the truth
    Ltruth = RegionLikelihood(build(clean=True), region, cmb_prior_sd=0.3)
    Ji, Sw, J, T, sol = sandwich(Ltruth, np.full(3, -3.1))
    assert np.allclose(sol["x"][:96], a_s, atol=0.05 * np.abs(a_s).max())
    emp = est.std(axis=0, ddof=1)
    assert np.all(np.abs(emp / np.sqrt(np.diag(Sw)) - 1.0) < 0.12), (emp, np.sqrt(np.diag(Sw)), np.sqrt(np.diag(Ji)))
    assert np.all(np.sqrt(np.diag(Sw)) > np.sqrt(np.diag(Ji)))
    assert np.all(np.abs(est.mean(axis=0) + 3.1) < 4 * emp / np.sqrt(len(est)) + 0.02)


def test_plugin_information_is_biased_by_T_and_debiased_sandwich_covers():
    """On noisy data J(a_hat) averages to J + T, J(a_hat) - T to J, and intervals from the debiased
    plug-in sandwich cover the truth at the nominal rate where the plain plug-in is too narrow."""
    kw = dict(n=384, nreg=3, dense=(0, 1, 2), noise=5.0, amp=6.0, seed=3)
    build, region, a_s = make(**kw)
    truth = np.full(3, -3.1)
    _, Sw0, J0, T0, _ = sandwich(RegionLikelihood(build(clean=True), region, cmb_prior_sd=0.3), truth)
    t = np.diag(T0) / np.diag(J0)
    assert t.min() > 0.1, t                      # the configuration must make the bias visible
    rng = np.random.default_rng(12)
    Jp_sum, est, z_plain, z_deb, nok = np.zeros((3, 3)), [], [], [], 0
    nrep = 300
    for _ in range(nrep):
        L = RegionLikelihood(build(rng), region, cmb_prior_sd=0.3)
        bh, _ = fit_profile(L, truth)
        est.append(bh)
        Jp_sum += L.information(truth)[0]
        Jh, Th, _ = L.information(bh)
        Ji = np.linalg.inv(Jh)
        z_plain.append((bh - truth) / np.sqrt(np.diag(Ji @ (Jh + Th) @ Ji)))
        _, cs, _, _, _, ok = sandwich_plugin(L, bh)
        if ok.all():
            nok += 1
            z_deb.append((bh - truth) / np.sqrt(np.diag(cs)))
    Jp_mean = Jp_sum / nrep
    assert np.allclose(np.diag(Jp_mean), np.diag(J0 + T0), rtol=0.05), (np.diag(Jp_mean), np.diag(J0 + T0), np.diag(J0))
    est, z_plain, z_deb = np.array(est), np.array(z_plain), np.array(z_deb)
    emp = est.std(axis=0, ddof=1)
    assert np.all(np.abs(emp / np.sqrt(np.diag(Sw0)) - 1.0) < 0.12), (emp, np.sqrt(np.diag(Sw0)))
    assert nok > 0.95 * nrep
    assert z_plain.std() > 1.08, (z_plain.std(), t)             # too narrow, by about (1 + t)^(3/2) / (1 + 2t)^(1/2)
    assert abs(z_deb.std() - 1.0) < 0.08, (z_deb.std(), z_plain.std(), t)
    assert abs(np.mean(np.abs(z_deb) < 1) - 0.683) < 0.05, np.mean(np.abs(z_deb) < 1)


def test_reference_frequency_shift_of_flat_marginal():
    """Moving the reference adds beta ln(ref'/ref) per synchrotron amplitude to (1/2) ln det A."""
    build, region, _ = make(dense=(0,), with_dust=True)
    bands = build()
    beta = np.array([-3.0, -3.2, -2.9, -3.15])
    L23 = RegionLikelihood(bands, region)
    b30 = [BandData(b.name, b.d, b.ninv, (lambda be, nu=nu: Band.delta("x", nu).power_law(be, 30.0)), b.dust)
           for b, nu in zip(bands, NUS)]
    L30 = RegionLikelihood(b30, region)
    nper = np.bincount(region)
    # ln det A(ref') = ln det A(ref) - 2 sum_p beta_p ln(ref'/ref)  ->  (1/2) ln det shifts by -sum beta ln(30/23)
    assert np.isclose(L30.marginal(beta) - L23.marginal(beta), -np.sum(nper * beta) * np.log(30.0 / 23.0), rtol=1e-9)
    assert np.isclose(L30.profile(beta), L23.profile(beta), rtol=1e-11)
