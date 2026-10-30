"""Configuration, validation, chain IO and one in-process Gibbs loop of the `likelihood: hier`
index-field group (M2.4 stage 1, commits 6 to 8)."""
import os

import h5py
import healpy as hp
import numpy as np
import pytest
from mpi4py import MPI
from pixell.bunch import Bunch

from commander4.compsep import index_field as ifd
from commander4.compsep import index_field_sampler as ifs
from commander4.compsep.processing import (MCMCSamplingGroupConfig, PerPixelSamplingGroupConfig,
                                           _resolve_mcmc_likelihoods, _validate_index_fields,
                                           process_compsep)
from commander4.data_models.detector_map import DetectorMap
from commander4.file_io import paths
from commander4.file_io.chain_writer import write_compsep_chain_to_file
from commander4.sky.comp_list import CompList

NSIDE = 4
NPIX = 12 * NSIDE**2


def _comp_cfg(name, cls, **params):
    cfg = Bunch(enabled=True, component_class=cls,
                params=Bunch(polarization="I", shortname=name, spatially_varying_MM=False,
                             Cl_prior_amplitude=None, lmax=3, units="uK_RJ", **params))
    object.__setattr__(cfg, "_name", name)
    return cfg


def _comp_list(template_path=None, nside_beta=1, likelihood_hier=True):
    sync_extra = {}
    if template_path is not None:
        sync_extra["amp_prior"] = Bunch(type="hier", template=template_path,
                                        mu_prior=[1.0, 0.5], log_s_prior=[-1.2, 1.0])
    comps = Bunch({
        "cmb": _comp_cfg("cmb", "CMB", nu_ref=1e-3,
                         amp_prior=Bunch(type="gaussian", mean=0.0, rms=70.0)),
        "sync": _comp_cfg("sync", "Synchrotron", beta=-3.0, nu_ref=30.0,
                          sample_spectral_index=True, spectral_index_nside=nside_beta,
                          spectral_index_prior=Bunch(type="gaussian", mean=-3.0, rms=0.2),
                          **sync_extra),
        "dust": _comp_cfg("dust", "ThermalDust", beta=1.6, T=20.0, nu_ref=545.0),
    })
    params = Bunch(compsep=Bunch(nside=NSIDE, double_precision=True, bands=Bunch()))
    return CompList.init_from_params(comps, params)


def _template(tmp_path, seed=0):
    path = os.path.join(tmp_path, "sync_template.fits")
    g = 10**(1.0 + 0.5 * np.random.default_rng(seed).random(NPIX))
    hp.write_map(path, g, overwrite=True)
    return path, g


def test_mcmc_group_config_nuts_block() -> None:
    group = MCMCSamplingGroupConfig.from_block(
        "g", Bunch(likelihood="hier", nuts=Bunch(n_warmup=20, parameterization="centred")))
    assert group.nuts_config.n_warmup == 20 and group.nuts_config.parameterization == "centred"
    assert MCMCSamplingGroupConfig.from_block("g", Bunch(likelihood="hier")).nuts_config.n_warmup \
        == 500
    with pytest.raises(ValueError, match="nuts"):
        MCMCSamplingGroupConfig.from_block("g", Bunch(likelihood="ridge", nuts=Bunch(n_warmup=1)))
    with pytest.raises(ValueError, match="parameterization"):
        MCMCSamplingGroupConfig.from_block(
            "g", Bunch(likelihood="hier", nuts=Bunch(parameterization="x")))


def test_hier_likelihood_resolution_and_field_validation() -> None:
    comp_list = _comp_list()
    amp = {"amps": PerPixelSamplingGroupConfig(name="amps", comps=("cmb", "sync", "dust"))}
    hier = MCMCSamplingGroupConfig(name="g", comps=("sync",), update_amplitude_groups=("amps",),
                                   likelihood="hier")
    resolved = _resolve_mcmc_likelihoods({"g": hier}, amp, "per_pixel", comp_list, report=False)
    assert resolved == {"g": "hier"}
    _validate_index_fields({"g": hier}, resolved, comp_list)
    with pytest.raises(ValueError, match="unavailable"):   # CG amplitudes: no per-pixel marginal
        _resolve_mcmc_likelihoods({"g": hier}, amp, "cg", comp_list, report=False)
    ridge = MCMCSamplingGroupConfig(name="g", comps=("sync",), update_amplitude_groups=("amps",),
                                    likelihood="ridge")
    with pytest.raises(ValueError, match="likelihood: hier"):
        _validate_index_fields({"g": ridge}, {"g": "ridge"}, comp_list)
    with pytest.raises(ValueError, match="no enabled MCMC group"):
        _validate_index_fields({}, {}, comp_list)


def test_load_amp_priors_from_params(tmp_path) -> None:
    path, g = _template(tmp_path)
    comp_list = _comp_list(template_path=path)
    ifs.load_amp_priors(comp_list, "per_pixel")
    by_name = {c.comp_name: c for c in comp_list}
    assert by_name["cmb"].amp_prior.type == "gaussian" and by_name["cmb"].amp_prior.rms == 70.0
    prior = by_name["sync"].amp_prior
    assert prior.type == "hier" and prior.mu_prior == (1.0, 0.5) and prior.log_s_prior == (-1.2, 1.0)
    np.testing.assert_allclose(prior.template, np.maximum(g, 1e-3 * np.mean(g)))
    assert prior.mu == 1.0 and prior.log_s == -1.2
    assert by_name["dust"].amp_prior is None
    with pytest.raises(ValueError, match="per-pixel"):
        ifs.load_amp_priors(_comp_list(template_path=path), "cg")
    bad = _comp_list()
    bad[0].comp_params.amp_prior = Bunch(type="flat")
    with pytest.raises(ValueError, match="amp_prior type"):
        ifs.load_amp_priors(bad, "per_pixel")


def test_chain_round_trip_of_index_field_and_hier_hyperparameters(tmp_path) -> None:
    path, _ = _template(tmp_path)
    params = Bunch(output=Bunch(dir=str(tmp_path), chains=Bunch(write=[1], maps_nside=2)))
    params.parameter_file_as_string = "parameters"
    paths.create_output_dirs(params.output)
    comp_list = _comp_list(template_path=path)
    ifs.load_amp_priors(comp_list, "per_pixel")
    rng = np.random.default_rng(1)
    for comp in comp_list:
        comp.allocate_empty_alms()
    sync = [c for c in comp_list if c.comp_name == "sync"][0]
    beta_map = -3.0 + 0.1 * rng.standard_normal(12)
    sync.set_beta_map(beta_map, 1)
    sync.amp_prior.mu, sync.amp_prior.log_s = 1.07, -1.33
    write_compsep_chain_to_file(comp_list, params, chain=1, iter=1, diagnostics={},
                                band_frequencies={"b": 30.0})
    file = os.path.join(paths.subdir(params, paths.CHAINS_COMPSEP), "chain01_iter0001.h5")
    with h5py.File(file) as f:
        np.testing.assert_array_equal(f["comps/sync/index_field/beta_map"][()], beta_map)
        assert int(f["comps/sync/index_field/beta_nside"][()]) == 1
        assert float(f["comps/sync/sed/beta"][()]) == pytest.approx(beta_map.mean())
        assert float(f["comps/sync/amp_prior/mu"][()]) == 1.07
        assert "index_field" not in f["comps/dust"] and "amp_prior" not in f["comps/cmb"]
    fresh = _comp_list(template_path=path)
    ifs.load_amp_priors(fresh, "per_pixel")
    fresh.load_initial_state(file, apply_init_from=False)
    sync2 = [c for c in fresh if c.comp_name == "sync"][0]
    np.testing.assert_array_equal(sync2.beta_map, beta_map)
    assert sync2.beta_nside == 1 and sync2.spatially_varying_MM
    assert (sync2.amp_prior.mu, sync2.amp_prior.log_s) == (1.07, -1.33)
    # A run at another field nside starts constant from the restored mean.
    other = _comp_list(template_path=path, nside_beta=2)
    ifs.load_amp_priors(other, "per_pixel")
    other.load_initial_state(file, apply_init_from=False)
    sync3 = [c for c in other if c.comp_name == "sync"][0]
    assert sync3.beta_map is None and sync3.beta == pytest.approx(beta_map.mean())


def test_process_compsep_runs_a_hier_group_end_to_end(tmp_path, monkeypatch) -> None:
    """One intensity band on COMM_SELF: per-pixel amplitude group plus a hier group with a sync
    field at nside 1; three Gibbs iterations run, the field and hyperparameters move, the
    NUTS state persists, and the diagnostics reach the chain writer."""
    ifs.reset_sampler_states()
    path, g = _template(tmp_path)
    comp_list = _comp_list(template_path=path)
    ifs.load_amp_priors(comp_list, "per_pixel")
    for comp in comp_list:
        comp.allocate_empty_alms()
    rng = np.random.default_rng(2)
    sky = 70.0 * rng.standard_normal(NPIX) + g * (30.0 / 30.0)**(-3.0) + 20.0
    inv_var = np.full(NPIX, 1.0 / 5.0**2)
    detector_data = DetectorMap(sky[None], inv_var[None], nu=30.0, fwhm=60.0, nside=NSIDE,
                                double_precision=True)
    mpi_info = Bunch(compsep=Bunch(comm=MPI.COMM_SELF, subcomm=MPI.COMM_SELF, rank=0, master=0,
                                   size=1))
    amp = PerPixelSamplingGroupConfig(name="amps", comps=("cmb", "sync", "dust"))
    hier = MCMCSamplingGroupConfig(name="g", comps=("sync",), update_amplitude_groups=("amps",),
                                   likelihood="hier", numstep=2, nuts=Bunch(n_warmup=4, max_depth=4))
    state = Bunch(band_name="b", band_identifier="b_I", target_pol="I", amplitude_method="per_pixel",
                  amplitude_groups={"amps": amp}, mcmc_groups={"g": hier}, double_precision=True,
                  num_threads=1, chisq_masks={}, nside_chisq=1, include_chisq_map=False,
                  include_residual_maps=False, mcmc_likelihoods={"g": "hier"})
    written = []
    monkeypatch.setattr("commander4.compsep.processing.write_compsep_chain_to_file",
                        lambda comp_list, params, chain, iter, diagnostics, band_frequencies:
                        written.append(diagnostics))
    sync = [c for c in comp_list if c.comp_name == "sync"][0]
    for it in range(1, 4):
        process_compsep(mpi_info, state, detector_data, iter=it, chain=1, params=Bunch(),
                        comp_list=comp_list, rng=np.random.default_rng(10 + it))
    assert sync.beta_map is not None and sync.beta_map.shape == (12,)
    assert not np.allclose(sync.beta_map, -3.0)
    assert sync.amp_prior.mu != 1.0
    stats = written[-1]["mcmc"]["g"]
    assert stats["nuts"]["n_transitions"] == 6
    np.testing.assert_array_equal(stats["state"]["sync"]["index_field"]["beta_map"], sync.beta_map)
    assert stats["state"]["sync"]["amp_prior"]["mu"] == sync.amp_prior.mu
    assert all(np.isfinite(c.alms).all() for c in comp_list)
    ifs.reset_sampler_states()


def test_chain_round_trip_of_a_mu_map(tmp_path) -> None:
    """A hier prior with mu_nside 1 writes its 12 multipliers to the chain and reads them back;
    a scalar run reads their mean, a run at another mu_nside keeps its initial mu, and a
    scalar chain read into a mu-map run is broadcast."""
    path, _ = _template(tmp_path)
    params = Bunch(output=Bunch(dir=str(tmp_path), chains=Bunch(write=[1], maps_nside=2)))
    params.parameter_file_as_string = "parameters"
    paths.create_output_dirs(params.output)

    def with_mu_nside(comp_list, mu_nside):
        ifs.load_amp_priors(comp_list, "per_pixel")
        for comp in comp_list:
            if comp.comp_name == "sync":
                comp.amp_prior = ifd.AmpPrior("hier", template=comp.amp_prior.template,
                                              mu_nside=mu_nside)
        return comp_list

    comp_list = with_mu_nside(_comp_list(template_path=path), 1)
    for comp in comp_list:
        comp.allocate_empty_alms()
    sync = [c for c in comp_list if c.comp_name == "sync"][0]
    mu_map = 1.0 + 0.1 * np.arange(12)
    sync.amp_prior.mu = mu_map.copy()
    sync.amp_prior.log_s = -1.33
    write_compsep_chain_to_file(comp_list, params, chain=1, iter=1, diagnostics={},
                                band_frequencies={"b": 30.0})
    file = os.path.join(paths.subdir(params, paths.CHAINS_COMPSEP), "chain01_iter0001.h5")
    with h5py.File(file) as f:
        np.testing.assert_array_equal(f["comps/sync/amp_prior/mu"][()], mu_map)
    for mu_nside, expected in ((1, mu_map), (0, mu_map.mean()), (2, np.ones(48))):
        fresh = with_mu_nside(_comp_list(template_path=path), mu_nside)
        fresh.load_initial_state(file, apply_init_from=False)
        got = [c for c in fresh if c.comp_name == "sync"][0].amp_prior
        np.testing.assert_allclose(got.mu, expected)
        assert got.log_s == -1.33
    # scalar chain into a mu-map run: broadcast
    scalar_list = _comp_list(template_path=path)
    ifs.load_amp_priors(scalar_list, "per_pixel")
    for comp in scalar_list:
        comp.allocate_empty_alms()
    s = [c for c in scalar_list if c.comp_name == "sync"][0]
    s.amp_prior.mu, s.amp_prior.log_s = 1.07, -1.0
    write_compsep_chain_to_file(scalar_list, params, chain=1, iter=2, diagnostics={},
                                band_frequencies={"b": 30.0})
    file2 = os.path.join(paths.subdir(params, paths.CHAINS_COMPSEP), "chain01_iter0002.h5")
    fresh = with_mu_nside(_comp_list(template_path=path), 1)
    fresh.load_initial_state(file2, apply_init_from=False)
    np.testing.assert_allclose([c for c in fresh if c.comp_name == "sync"][0].amp_prior.mu,
                               np.full(12, 1.07))
