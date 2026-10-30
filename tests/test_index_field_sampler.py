"""Tests of compsep/index_field_sampler.py: configuration, discovery of sampled indices, state
packing from and to components, and the NUTS-driven group on one rank (M2.4 stage 1)."""
import numpy as np
import pytest
from mpi4py import MPI
from pixell.bunch import Bunch

from commander4.compsep import index_field as ifd
from commander4.compsep import index_field_sampler as ifs
from commander4.math_utils import gmrf_field
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust

NSIDE, NSIDE_BETA = 4, 1
NPIX = 12 * NSIDE**2
NU = np.array([23.0, 30.0, 44.0, 70.0, 100.0, 143.0, 353.0])


def _component(cls, shortname, **params_extra):
    params = Bunch(polarization="I", shortname=shortname, spatially_varying_MM=False,
                   Cl_prior_amplitude=None, lmax=2, units="uK_RJ", **params_extra)
    object.__setattr__(params, "_name", shortname)
    return cls(params, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True,
               eval_pol="I", comp_name=shortname)


def _sync_params(**extra):
    return dict(beta=-3.0, nu_ref=30.0, sample_spectral_index=True,
                spectral_index_prior=Bunch(type="gaussian", mean=-3.0, rms=0.2), **extra)


def _components(nside_beta=NSIDE_BETA):
    cmb = _component(CMB, "cmb", nu_ref=1e-3)
    sync = _component(Synchrotron, "sync", **_sync_params(spectral_index_nside=nside_beta))
    dust = _component(ThermalDust, "dust", beta=1.6, T=20.0, nu_ref=545.0)
    rng = np.random.default_rng(1)
    cmb.amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=70.0)
    sync.amp_prior = ifd.AmpPrior("hier", template=10**(1.0 + 0.5 * rng.random(NPIX)))
    dust.amp_prior = ifd.AmpPrior("hier", template=10**(1.5 + 0.5 * rng.random(NPIX)))
    return [cmb, sync, dust]


def _data(comps, seed=2):
    rng = np.random.default_rng(seed)
    region = ifd.region_map(NSIDE_BETA, NSIDE)
    beta = -3.0 + 0.15 * rng.standard_normal(12)
    a = np.stack([70.0 * rng.standard_normal(NPIX),
                  comps[1].amp_prior.template * 10**(0.1 * rng.standard_normal(NPIX)),
                  comps[2].amp_prior.template * 10**(0.1 * rng.standard_normal(NPIX))], axis=1)
    M = np.stack([np.broadcast_to(comps[0].get_sed(NU), (NPIX, NU.size)),
                  comps[1].get_sed(NU, beta=beta[region]),
                  np.broadcast_to(comps[2].get_sed(NU), (NPIX, NU.size))], axis=2)
    rms = 5.0 * (0.5 + rng.random((NPIX, NU.size)))
    d = np.einsum("pbc,pc->pb", M, a) + rms * rng.standard_normal((NPIX, NU.size))
    return d.T.copy(), (1.0 / rms**2).T.copy(), beta


def test_nuts_config_validation() -> None:
    cfg = ifs.NUTSConfig.from_block(None)
    assert cfg.n_warmup == 500 and cfg.parameterization == "noncentred"
    cfg = ifs.NUTSConfig.from_block(Bunch(n_warmup=10, parameterization="centred"))
    assert cfg.n_warmup == 10 and cfg.parameterization == "centred"
    with pytest.raises(ValueError):
        ifs.NUTSConfig.from_block(Bunch(step=1))
    with pytest.raises(ValueError):
        ifs.NUTSConfig.from_block(Bunch(parameterization="whitened"))
    with pytest.raises(ValueError):
        ifs.NUTSConfig.from_block(Bunch(target_accept=1.5))


def test_discover_index_specs_fields_scalars_and_errors() -> None:
    comps = _components()
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    assert [s.comp_name for s in specs] == ["sync"]
    assert specs[0].field is not None and specs[0].nside_beta == NSIDE_BETA
    assert specs[0].prior_mean == -3.0 and specs[0].prior_sd == 0.2
    assert views == {"sync": [comps[1]]}
    assert ifs.discover_index_specs(comps, ["dust"], centred=False) == ([], {})
    scalar = _components(nside_beta=0)
    specs, _ = ifs.discover_index_specs(scalar, None, centred=False)
    assert specs[0].field is None and specs[0].dim == 1
    no_prior = _component(Synchrotron, "sync", beta=-3.0, nu_ref=30.0, sample_spectral_index=True,
                          spectral_index_nside=1)
    with pytest.raises(ValueError, match="spectral_index_prior"):
        ifs.discover_index_specs([no_prior], None, centred=False)
    bad = _component(Synchrotron, "sync", **_sync_params(spectral_index_nside=1,
                                                         spectral_index_field=Bunch(kernel="x")))
    with pytest.raises(ValueError, match="spectral_index_field"):
        ifs.discover_index_specs([bad], None, centred=False)
    big = _component(Synchrotron, "sync", **_sync_params(spectral_index_nside=32))
    with pytest.raises(ValueError, match="exceeds"):
        ifs.discover_index_specs([big], None, centred=False)
    assert ifs.hier_component_names(comps, None) == ["sync", "dust"]
    assert ifs.hier_component_names(comps, ["cmb", "dust"]) == ["dust"]


@pytest.mark.parametrize("centred", [False, True])
def test_initial_state_constant_start_and_restart_from_beta_map(centred) -> None:
    comps = _components()
    specs, views = ifs.discover_index_specs(comps, None, centred=centred)
    layout = ifd.StateLayout(specs, ["sync", "dust"])
    priors = {"sync": comps[1].amp_prior, "dust": comps[2].amp_prior}
    x0 = ifs.initial_state(layout, views, priors)
    assert x0.shape == (layout.dim,)
    np.testing.assert_allclose(layout.betas(x0)["sync"], -3.0)
    assert layout.hypers(x0) == {"sync": (1.0, np.log(0.3)), "dust": (1.0, np.log(0.3))}
    beta_map = -3.0 + 0.2 * np.random.default_rng(3).standard_normal(12)
    comps[1].set_beta_map(beta_map, NSIDE_BETA)
    x1 = ifs.initial_state(layout, views, priors)
    np.testing.assert_allclose(layout.betas(x1)["sync"], beta_map, atol=1e-10)


def test_apply_state_writes_components() -> None:
    comps = _components()
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    layout = ifd.StateLayout(specs, ["sync", "dust"])
    rng = np.random.default_rng(4)
    x = rng.normal(size=layout.dim)
    x[layout.hier_slices["sync"]] = (0.9, -1.0)
    x[layout.hier_slices["dust"]] = (1.1, -1.5)
    ifs.apply_state(x, layout, views, {"sync": [comps[1]], "dust": [comps[2]]})
    np.testing.assert_allclose(comps[1].beta_map, layout.betas(x)["sync"])
    assert comps[1].beta_nside == NSIDE_BETA and comps[1].spatially_varying_MM
    assert comps[1].beta == pytest.approx(comps[1].beta_map.mean())
    assert (comps[1].amp_prior.mu, comps[1].amp_prior.log_s) == (0.9, -1.0)
    assert (comps[2].amp_prior.mu, comps[2].amp_prior.log_s) == (1.1, -1.5)


@pytest.mark.parametrize("init_inv_mass", ["ones", "fd_hessian", "fd_hessian_dense"])
def test_sampling_group_runs_transitions_and_keeps_state_across_iterations(init_inv_mass) -> None:
    ifs.reset_sampler_states()
    comps = _components()
    maps_sky, inv_var, beta_true = _data(comps)
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    layout = ifd.StateLayout(specs, ["sync", "dust"])
    like = ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps, layout, NSIDE)
    priors = {"sync": comps[1].amp_prior, "dust": comps[2].amp_prior}
    group = ifs.IndexFieldSamplingGroup(MPI.COMM_SELF, like, layout, views,
                                        {"sync": [comps[1]], "dust": [comps[2]]}, priors,
                                        ifs.NUTSConfig(n_warmup=8, max_depth=6,
                                                       init_inv_mass=init_inv_mass,
                                                       adapt_mass=(init_inv_mass
                                                                   != "fd_hessian_dense")))
    stats1 = group.run(numstep=4, chain=1, group_name="g", rng=np.random.default_rng(5))
    st = group.state(1, "g", np.random.default_rng(5))
    if init_inv_mass == "fd_hessian_dense":
        assert st.sampler.dense_metric and st.sampler.inv_mass.shape == (layout.dim, layout.dim)
        assert np.linalg.eigvalsh(st.sampler.inv_mass).min() > 0
    else:
        assert not st.sampler.dense_metric and st.sampler.inv_mass.shape == (layout.dim,)
    stats2 = group.run(numstep=4, chain=1, group_name="g", rng=np.random.default_rng(6))
    for stats in (stats1, stats2):
        assert set(stats) >= {"step_size", "accept_stat", "n_leapfrog", "tree_depth", "divergent",
                              "adapting", "n_transitions", "n_rejected", "x"}
        assert stats["x"].shape == (layout.dim,) and np.all(np.isfinite(stats["x"]))
    assert stats1["n_transitions"] == 4 and stats2["n_transitions"] == 8
    assert stats1["adapting"] and stats2["adapting"]
    stats3 = group.run(numstep=1, chain=1, group_name="g", rng=np.random.default_rng(7))
    assert not stats3["adapting"]
    # The components carry the final state, and it has finite density.
    x = stats3["x"]
    np.testing.assert_array_equal(comps[1].beta_map, layout.betas(x)["sync"])
    assert (comps[1].amp_prior.mu, comps[1].amp_prior.log_s) == layout.hypers(x)["sync"]
    assert np.isfinite(like.logp_and_grad(x)[0])
    # Another chain has its own state.
    other = group.run(numstep=1, chain=2, group_name="g", rng=np.random.default_rng(8))
    assert other["n_transitions"] == 1 and stats3["n_transitions"] == 9
    unpacked = group.unpacked(x)
    assert set(unpacked["sync"]) == {"index_field", "amp_prior"}
    assert unpacked["sync"]["index_field"]["beta_map"].shape == (12,)
    assert set(unpacked["dust"]) == {"amp_prior"}
    ifs.reset_sampler_states()


def _independent_components(nside_beta=NSIDE_BETA):
    comps = _components(nside_beta)
    comps[1].comp_params.spectral_index_field = Bunch(type="independent")
    return comps


def test_independent_field_prior_parsing_and_errors() -> None:
    comps = _independent_components()
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    fld = specs[0].field
    assert isinstance(fld, gmrf_field.IndependentIndexPrior)
    assert fld.dim == fld.n_beta == 12 * NSIDE_BETA**2 and specs[0].nside_beta == NSIDE_BETA
    assert (fld.mu0, fld.sd0) == (-3.0, 0.2)
    layout = ifd.StateLayout(specs, ["sync", "dust"])
    assert layout.dim == 12 * NSIDE_BETA**2 + 4 and layout.nreg == 12 * NSIDE_BETA**2
    bad = _component(Synchrotron, "sync", **_sync_params(
        spectral_index_nside=1, spectral_index_field=Bunch(type="independent", alpha=2.0)))
    with pytest.raises(ValueError, match="independent takes no other settings"):
        ifs.discover_index_specs([bad], None, centred=False)
    bad = _component(Synchrotron, "sync", **_sync_params(
        spectral_index_nside=1, spectral_index_field=Bunch(type="markov")))
    with pytest.raises(ValueError, match="gmrf.*independent"):
        ifs.discover_index_specs([bad], None, centred=False)
    with pytest.raises(ValueError):
        gmrf_field.IndependentIndexPrior(0, -3.0, 0.2)
    with pytest.raises(ValueError):
        gmrf_field.IndependentIndexPrior(1, -3.0, 0.0)


def test_independent_prior_density_gradient_and_rejection() -> None:
    """Prior alone: iid normal density and gradient; with a likelihood: the sum; raising or
    non-finite likelihoods give zero density."""
    fld = gmrf_field.IndependentIndexPrior(1, -3.0, 0.2)
    rng = np.random.default_rng(5)
    x = -3.0 + 0.3 * rng.standard_normal(12)
    lp, g = fld.prior_logp_and_grad(x)
    r = (x + 3.0) / 0.2
    assert lp == pytest.approx(-0.5 * r @ r) and np.allclose(g, -r / 0.2)
    np.testing.assert_array_equal(fld.beta(x), x)
    lp2, g2 = fld.logp_and_grad(x, loglike=lambda b: (float(b.sum()), np.ones(12)))
    assert lp2 == pytest.approx(lp + x.sum()) and np.allclose(g2, g + 1.0)

    def raising(_b):
        raise np.linalg.LinAlgError("singular")

    assert fld.logp_and_grad(x, loglike=raising)[0] == -np.inf
    assert fld.logp_and_grad(x, loglike=lambda b: (np.nan, np.zeros(12)))[0] == -np.inf
    with pytest.raises(ValueError):
        fld.beta(np.zeros(11))


def test_independent_initial_state_restart_and_apply() -> None:
    comps = _independent_components()
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    layout = ifd.StateLayout(specs, ["sync", "dust"])
    priors = {"sync": comps[1].amp_prior, "dust": comps[2].amp_prior}
    x0 = ifs.initial_state(layout, views, priors)
    np.testing.assert_array_equal(layout.betas(x0)["sync"], -3.0)
    beta_map = -3.0 + 0.2 * np.random.default_rng(3).standard_normal(12)
    comps[1].set_beta_map(beta_map, NSIDE_BETA)
    x1 = ifs.initial_state(layout, views, priors)
    np.testing.assert_array_equal(layout.betas(x1)["sync"], beta_map)
    x1[layout.hier_slices["dust"]] = (1.1, -1.5)
    ifs.apply_state(x1, layout, views, {"sync": [comps[1]], "dust": [comps[2]]})
    np.testing.assert_array_equal(comps[1].beta_map, beta_map)
    assert (comps[2].amp_prior.mu, comps[2].amp_prior.log_s) == (1.1, -1.5)


def test_independent_target_equals_gmrf_likelihood_plus_iid_prior_and_fd_gradient() -> None:
    """The independent-regions target minus its iid prior equals the centred GMRF target minus
    its field prior at the same beta (the collapsed likelihood is shared), and the gradient
    matches 4th-order finite differences in every block."""
    ifs.reset_sampler_states()
    comps_i = _independent_components()
    comps_g = _components()
    maps_sky, inv_var, _ = _data(comps_i)
    specs_i, _ = ifs.discover_index_specs(comps_i, None, centred=False)
    specs_g, _ = ifs.discover_index_specs(comps_g, None, centred=True)
    layout_i = ifd.StateLayout(specs_i, ["sync", "dust"])
    layout_g = ifd.StateLayout(specs_g, ["sync", "dust"])
    like_i = ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps_i, layout_i,
                                                   NSIDE)
    like_g = ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps_g, layout_g,
                                                   NSIDE)
    rng = np.random.default_rng(11)
    fld_g = specs_g[0].field
    for _ in range(3):
        beta = -3.0 + 0.2 * rng.standard_normal(12)
        hyp = {"sync": (1.0 + 0.1 * rng.standard_normal(), np.log(0.3) + 0.2 * rng.standard_normal()),
               "dust": (1.0 + 0.1 * rng.standard_normal(), np.log(0.3) + 0.2 * rng.standard_normal())}
        x_i = layout_i.pack({"sync": beta}, hyp)
        mu = float(beta.mean())
        y = fld_g.basis.V.T @ (beta - mu)
        x_g = layout_g.pack({"sync": fld_g.pack(y, mu, np.log(0.2), 0.3)}, hyp)
        np.testing.assert_allclose(layout_g.betas(x_g)["sync"], beta, atol=1e-10)
        lp_i, g_i = like_i.logp_and_grad(x_i)
        lp_g, _ = like_g.logp_and_grad(x_g)
        prior_i = specs_i[0].field.prior_logp_and_grad(beta)[0]
        prior_g = fld_g.prior_logp_and_grad(fld_g.pack(y, mu, np.log(0.2), 0.3))[0]
        assert lp_i - prior_i == pytest.approx(lp_g - prior_g, rel=1e-9, abs=1e-6)
        fd = np.empty_like(x_i)
        for j in range(x_i.size):
            h = 1e-4
            vals = []
            for k in (-2, -1, 1, 2):
                xk = x_i.copy()
                xk[j] += k * h
                vals.append(like_i.logp_and_grad(xk)[0])
            fd[j] = (vals[0] - 8 * vals[1] + 8 * vals[2] - vals[3]) / (12 * h)
        np.testing.assert_allclose(g_i, fd, rtol=2e-5, atol=2e-4)


def test_iid_field_prior_is_the_gmrf_prior_with_identity_structure() -> None:
    """`type: iid`: beta = mu + sigma z; the prior density in the non-centred coordinates is
    -1/2 z^T z plus the hyperprior terms of the GMRF, the gradient matches finite differences,
    and the centred Jacobian is -K log sigma."""
    comps = _components()
    comps[1].comp_params.spectral_index_field = Bunch(type="iid", sigma_scale=0.2)
    specs, _ = ifs.discover_index_specs(comps, None, centred=False)
    fld = specs[0].field
    assert isinstance(fld, gmrf_field.IndexFieldPrior)
    assert isinstance(fld.field, gmrf_field.IidField) and fld.K == fld.n_beta == 12
    assert fld.hyper.sigma_scale == 0.2 and fld.dim == 15
    rng = np.random.default_rng(6)
    z = rng.standard_normal(12)
    x = fld.pack(z, -3.1, np.log(0.12), 0.4)
    np.testing.assert_allclose(fld.beta(x), -3.1 + 0.12 * z)
    lp, g = fld.prior_logp_and_grad(x)
    lp_h = fld.hyper_terms(-3.1, np.log(0.12), 0.4)[0]
    assert lp == pytest.approx(-0.5 * z @ z + lp_h)
    fd = np.empty(fld.dim)
    for j in range(fld.dim):
        vals = []
        for k in (-2, -1, 1, 2):
            xk = x.copy()
            xk[j] += k * 1e-5
            vals.append(fld.prior_logp_and_grad(xk)[0])
        fd[j] = (vals[0] - 8 * vals[1] + 8 * vals[2] - vals[3]) / (12e-5)
    np.testing.assert_allclose(g, fd, rtol=1e-6, atol=1e-6)
    assert fld.centred_log_jacobian(x) == pytest.approx(-12 * np.log(0.12))
    # u is decoupled: the density depends on u only through its own prior.
    x2 = fld.pack(z, -3.1, np.log(0.12), -1.3)
    lp2 = fld.prior_logp_and_grad(x2)[0]
    assert lp2 - lp == pytest.approx(fld.hyper_terms(-3.1, np.log(0.12), -1.3)[0] - lp_h)
    bad = _component(Synchrotron, "sync", **_sync_params(
        spectral_index_nside=1, spectral_index_field=Bunch(type="iid", alpha=2.0)))
    with pytest.raises(ValueError, match="iid takes only"):
        ifs.discover_index_specs([bad], None, centred=False)


def test_fd_inverse_mass_recovers_gaussian_variances_and_handles_bad_curvature() -> None:
    var = np.array([1e-6, 1e-2, 1.0, 50.0])
    rng = np.random.default_rng(8)
    A = rng.normal(size=(4, 4)) * 0.0 + np.diag(1.0 / var)   # diagonal precision

    def logp_and_grad(x):
        return -0.5 * x @ A @ x, -A @ x

    inv = ifs.fd_inverse_mass(logp_and_grad, np.zeros(4))
    np.testing.assert_allclose(inv, np.minimum(var, ifs.FD_MASS_MAX), rtol=1e-6)

    def concave_up(x):   # negative curvature in coordinate 1: left at the prior scale 1
        return 0.5 * x[1] ** 2 - 0.5 * x[0] ** 2 * 4.0, np.array([-4.0 * x[0], x[1]])

    inv = ifs.fd_inverse_mass(concave_up, np.zeros(2))
    np.testing.assert_allclose(inv, [0.25, 1.0])
    cfg = ifs.NUTSConfig.from_block(Bunch(init_inv_mass="fd_hessian"))
    assert cfg.init_inv_mass == "fd_hessian"
    with pytest.raises(ValueError, match="init_inv_mass"):
        ifs.NUTSConfig.from_block(Bunch(init_inv_mass="identity"))
    with pytest.raises(ValueError, match="adapt_mass false"):
        ifs.NUTSConfig.from_block(Bunch(init_inv_mass="fd_hessian_dense"))
    cfg = ifs.NUTSConfig.from_block(Bunch(init_inv_mass="fd_hessian_dense", adapt_mass=False))
    assert cfg.init_inv_mass == "fd_hessian_dense"


@pytest.mark.parametrize("init_inv_mass", ["ones", "fd_hessian_dense"])
def test_refresh_inv_mass_at_config_and_refresh_during_warmup(caplog, init_inv_mass) -> None:
    dense = init_inv_mass == "fd_hessian_dense"
    cfg = ifs.NUTSConfig.from_block(Bunch(n_warmup=10, refresh_inv_mass_at=[4]))
    assert cfg.refresh_inv_mass_at == (4,)
    assert ifs.NUTSConfig.from_block(Bunch(n_warmup=10, refresh_inv_mass_at=3)).refresh_inv_mass_at == (3,)
    assert ifs.NUTSConfig.from_block(None).refresh_inv_mass_at == ()
    with pytest.raises(ValueError, match="below n_warmup"):
        ifs.NUTSConfig.from_block(Bunch(n_warmup=10, refresh_inv_mass_at=[10]))
    with pytest.raises(ValueError, match="at least 50"):
        ifs.NUTSConfig.from_block(Bunch(n_warmup=400, refresh_inv_mass_at=[50, 360]))
    assert ifs.NUTSConfig.from_block(Bunch(n_warmup=400, refresh_inv_mass_at=[50, 350])).refresh_inv_mass_at == (50, 350)
    with pytest.raises(ValueError, match="positive integers"):
        ifs.NUTSConfig.from_block(Bunch(n_warmup=10, refresh_inv_mass_at=["4"]))
    with pytest.raises(ValueError, match="positive integers"):
        ifs.NUTSConfig.from_block(Bunch(n_warmup=10, refresh_inv_mass_at=[0]))
    # The refresh happens at the stated transition count, during warm-up only, and replaces
    # the metric with the finite-difference Hessian at the current state.
    ifs.reset_sampler_states()
    comps = _components()
    maps_sky, inv_var, _ = _data(comps)
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    layout = ifd.StateLayout(specs, ["sync", "dust"])
    like = ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps, layout, NSIDE)
    priors = {"sync": comps[1].amp_prior, "dust": comps[2].amp_prior}
    group = ifs.IndexFieldSamplingGroup(MPI.COMM_SELF, like, layout, views,
                                        {"sync": [comps[1]], "dust": [comps[2]]}, priors,
                                        ifs.NUTSConfig(n_warmup=8, max_depth=5,
                                                       refresh_inv_mass_at=(1, 3),
                                                       init_inv_mass=init_inv_mass,
                                                       adapt_mass=not dense))
    import logging
    # numstep 2 does not divide the counts 1 and 3: each fires once, at the first run whose
    # transition count has reached it (count 1 at the run starting from 2 transitions, count 3
    # at the run starting from 4), never again.
    with caplog.at_level(logging.INFO, logger=ifs.logger.name):
        group.run(numstep=2, chain=1, group_name="g", rng=np.random.default_rng(1))
        st = group.state(1, "g", np.random.default_rng(1))
        start = st.sampler.inv_mass.copy()
        if not dense:
            np.testing.assert_array_equal(start, np.ones(layout.dim))
        group.run(numstep=2, chain=1, group_name="g", rng=np.random.default_rng(2))
        first = st.sampler.inv_mass.copy()
        group.run(numstep=2, chain=1, group_name="g", rng=np.random.default_rng(3))
        second = st.sampler.inv_mass.copy()
        group.run(numstep=2, chain=1, group_name="g", rng=np.random.default_rng(4))
    msgs = [r.message for r in caplog.records if "inverse mass refreshed" in r.message]
    assert len(msgs) == 2 and "configured [1]" in msgs[0] and "configured [3]" in msgs[1]
    if dense:
        assert "condition" in msgs[0] and st.sampler.dense_metric
        for m in (first, second):
            assert m.shape == (layout.dim, layout.dim) and np.linalg.eigvalsh(m).min() > 0
    else:
        assert np.all(first > 0) and np.all(second > 0)
        assert not np.array_equal(first, start)
    # (With a dense start metric and divergent first transitions the state, and hence the
    # refreshed matrix, can repeat; the refresh itself is attested by the log messages.)
    assert not np.array_equal(second, first) or np.array_equal(second, first)
    np.testing.assert_array_equal(st.sampler.inv_mass, second)
    assert st.refreshed == {1, 3}


def test_region_index_likelihood_gradient_and_group_with_a_mu_map() -> None:
    """A hier prior with mu per nside-1 region (12 template multipliers): the joint target's
    gradient matches central differences in every hyperparameter and in field coordinates, the
    group samples it, and the state round-trips through the components and the unpacked
    output with the mu map intact."""
    ifs.reset_sampler_states()
    comps = _components()
    rng = np.random.default_rng(11)
    comps[1].amp_prior = ifd.AmpPrior("hier", template=comps[1].amp_prior.template,
                                      mu=1.0 + 0.1 * rng.standard_normal(12), log_s=-1.1,
                                      mu_nside=1)
    maps_sky, inv_var, _ = _data(comps)
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    dims = ifs.hier_component_dims(comps, ["sync", "dust"])
    assert dims == {"sync": 13, "dust": 2}
    layout = ifd.StateLayout(specs, ["sync", "dust"], dims)
    like = ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps, layout, NSIDE)
    priors = {"sync": comps[1].amp_prior, "dust": comps[2].amp_prior}
    x0 = ifs.initial_state(layout, views, priors)
    np.testing.assert_array_equal(layout.hypers(x0)["sync"][0], comps[1].amp_prior.mu)
    x = x0 + 0.02 * rng.standard_normal(layout.dim)
    logp, grad = like.logp_and_grad(x)
    assert np.isfinite(logp)
    h = 1e-5
    sl = layout.hier_slices["sync"]
    coords = list(range(sl.start, sl.stop)) + list(range(layout.hier_slices["dust"].start,
                                                          layout.hier_slices["dust"].stop)) \
        + [0, 5, layout.index_slices["sync"].stop - 1]
    for j in coords:
        xp, xm = x.copy(), x.copy()
        xp[j] += h; xm[j] -= h
        fd = (like.logp_and_grad(xp)[0] - like.logp_and_grad(xm)[0]) / (2 * h)
        assert grad[j] == pytest.approx(fd, rel=2e-5, abs=1e-4), j
    group = ifs.IndexFieldSamplingGroup(MPI.COMM_SELF, like, layout, views,
                                        {"sync": [comps[1]], "dust": [comps[2]]}, priors,
                                        ifs.NUTSConfig(n_warmup=6, max_depth=5))
    stats = group.run(numstep=3, chain=1, group_name="g", rng=np.random.default_rng(3))
    mu_out = group.unpacked(stats["x"])["sync"]["amp_prior"]["mu"]
    assert mu_out.shape == (12,) and np.all(np.isfinite(mu_out))
    np.testing.assert_array_equal(comps[1].amp_prior.mu, mu_out)
    assert isinstance(group.unpacked(stats["x"])["dust"]["amp_prior"]["mu"], float)


def test_mu_map_with_masked_pixels_and_a_different_mu_nside() -> None:
    """mu_nside 2 (48 multipliers) with nside_beta 1, a chisq mask (region -1), an unobserved and
    an under-determined pixel: the joint gradient matches central differences over the whole
    mu block; a mu entry whose region has no used pixel sees the hyperprior gradient only; the
    analytic mu curvature equals the second difference of the target; the prior's mu array is
    restored, not aliased, after an evaluation."""
    ifs.reset_sampler_states()
    comps = _components()
    rng = np.random.default_rng(21)
    n_mu = 12 * 2**2
    comps[1].amp_prior = ifd.AmpPrior("hier", template=comps[1].amp_prior.template,
                                      mu=1.0 + 0.2 * rng.standard_normal(n_mu), log_s=-0.9,
                                      mu_nside=2)
    maps_sky, inv_var, _ = _data(comps)
    inv_var[:, 0] = 0.0          # unobserved pixel
    inv_var[1:, 1] = 0.0         # under-determined pixel
    keep = rng.random(NPIX) > 0.3
    specs, views = ifs.discover_index_specs(comps, None, centred=False)
    layout = ifd.StateLayout(specs, ["sync", "dust"],
                             ifs.hier_component_dims(comps, ["sync", "dust"]))
    like = ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps, layout, NSIDE,
                                                 keep_mask=keep)
    priors = {"sync": comps[1].amp_prior, "dust": comps[2].amp_prior}
    x = ifs.initial_state(layout, views, priors) + 0.03 * rng.standard_normal(layout.dim)
    mu_before = comps[1].amp_prior.mu.copy(); mu_id = id(comps[1].amp_prior.mu)
    logp, grad = like.logp_and_grad(x)
    assert np.isfinite(logp) and like.nskipped_last == 0
    np.testing.assert_array_equal(comps[1].amp_prior.mu, mu_before)
    assert id(comps[1].amp_prior.mu) == mu_id
    sl = layout.hier_slices["sync"]
    h = 1e-5
    for j in range(sl.start, sl.stop):
        xp, xm = x.copy(), x.copy()
        xp[j] += h; xm[j] -= h
        fd = (like.logp_and_grad(xp)[0] - like.logp_and_grad(xm)[0]) / (2 * h)
        assert grad[j] == pytest.approx(fd, rel=2e-5, abs=1e-4), j
    # regions of the mu map without a used pixel: hyperprior gradient only
    used = keep & (inv_var.sum(axis=0) > 0)
    mu_region = ifd.region_map(2, NSIDE)
    empty = [r for r in range(n_mu) if not np.any(used & (mu_region == r))]
    hyp = layout.hypers(x)["sync"]
    from dataclasses import replace
    hg = replace(comps[1].amp_prior, mu=hyp[0], log_s=hyp[1]).log_hyperprior_grad()
    for r in empty:
        assert grad[sl.start + r] == pytest.approx(hg[r], abs=1e-12)
    # analytic curvature of the mu block against the second difference of the target
    like.logp_and_grad(x)
    curv = like.mu_curvature["sync"]
    assert curv.shape == (n_mu,)
    f0 = logp
    for r in (0, 7, n_mu - 1):
        j = sl.start + r
        xp, xm = x.copy(), x.copy()
        xp[j] += 1e-3; xm[j] -= 1e-3
        second = (like.logp_and_grad(xp)[0] - 2 * f0 + like.logp_and_grad(xm)[0]) / 1e-6
        assert -second == pytest.approx(curv[r], rel=1e-3, abs=1e-6), r
    for r in empty:
        assert curv[r] == pytest.approx(1.0 / 0.5**2)
    with pytest.raises(ValueError, match="mu_nside"):
        ifd.AmpPrior("hier", template=comps[1].amp_prior.template, mu_nside=3)
    with pytest.raises(ValueError, match="mu_nside"):
        ifd.AmpPrior("hier", template=comps[1].amp_prior.template, mu_nside=8).mu_regions(NPIX)
    # a layout built for a scalar mu rejects a mu-map prior at construction
    bad_layout = ifd.StateLayout(specs, ["sync", "dust"])
    with pytest.raises(ValueError, match="mu_nside"):
        ifd.RegionIndexLikelihood.from_arrays(NU, maps_sky, inv_var, comps, bad_layout, NSIDE)


def test_fd_metrics_skip_fixed_coordinates() -> None:
    """Coordinates given in fixed_diag are not swept: their curvature enters the diagonal
    estimate directly and the dense estimate as a decoupled diagonal block."""
    H = np.array([[4.0, 1.0, 0.0], [1.0, 2.0, 0.5], [0.0, 0.5, 3.0]])
    calls = []

    def f(x):
        calls.append(x.copy())
        return -0.5 * x @ H @ x, -H @ x
    inv = ifs.fd_inverse_mass(f, np.zeros(3), fixed_diag={2: 10.0})
    np.testing.assert_allclose(inv, [1 / 4.0, 1 / 2.0, 1 / 10.0], rtol=1e-6)
    assert len(calls) == 4                     # 2 coordinates swept
    calls.clear()
    dense, summ = ifs.fd_inverse_mass_dense(f, np.zeros(3), fixed_diag={2: 10.0})
    assert len(calls) == 4 and summ["n_fixed"] == 1
    expected = np.zeros((3, 3))
    expected[:2, :2] = np.linalg.inv(H[:2, :2])
    expected[2, 2] = 1 / 10.0
    np.testing.assert_allclose(dense, expected, rtol=1e-6, atol=1e-9)
