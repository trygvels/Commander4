"""End-to-end check of the MCMC spectral-index likelihoods on a small maps-only simulation.

Simulates per-pixel CMB + synchrotron + thermal dust intensity maps in 7 bands (white noise, no
beams beyond a common smoothing), runs Commander4 once per `likelihood` setting of the MCMC group, and compares the beta
chains with the exact posterior, evaluated on a grid by an independent NumPy implementation of the
per-pixel marginal likelihood.

    python aux/validate_mcmc_likelihood.py make    --dir /tmp/c4val
    for m in marginal conditional resolve_chisq; do
        mpirun -n 7 commander4 -p /tmp/c4val/param_$m.yml
    done
    python aux/validate_mcmc_likelihood.py analyze --dir /tmp/c4val

Mode `hier_field` (E018, M2.4 stage 1): a synchrotron index field at nside_beta 2 with a
hierarchical amplitude prior, sampled with `likelihood: hier` (NUTS), compared with the E016
prototype chain on the same data (nextgen/prototypes/indexfield, imported by path):

    python aux/validate_mcmc_likelihood.py make    --mode hier_field --dir /tmp/c4val_hier
    mpirun -n 12 commander4 -p /tmp/c4val_hier/param_hier_field.yml
    mpirun -n 12 commander4 -p /tmp/c4val_hier/param_hier_scalar.yml
    python aux/validate_mcmc_likelihood.py analyze --mode hier_field --dir /tmp/c4val_hier
"""
import argparse
import glob
import json
import os
import sys

import h5py
import healpy as hp
import numpy as np
from pixell.bunch import Bunch

from commander4.sky.comp_list import CompList

FREQS = [30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0]
RMS = [60.0, 60.0, 70.0, 24.0, 16.0, 24.0, 80.0]  # uK_RJ per pixel, before smoothing
TRUTH = {"sync": -3.0, "dust": 1.6}
START = {"sync": -2.9, "dust": 1.595}      # chains start several sd away from the truth
PROPOSAL = {"sync": 0.08, "dust": 0.06}  # random-walk step lengths, about the ridge target sd
BOUNDS = {"sync": [-4.0, -2.0], "dust": [1.0, 2.5]}  # hard limits (flat prior inside)
NSIDE = 16
FWHM_ARCMIN = 480.0                        # common beam, about two pixels
MODES = {"ridge": 50, "marginal": 50, "conditional": 50, "resolve_chisq": 50}  # -> numstep
NITER = 250
BURN = 80


def _components(beta_sync: float, beta_dust: float) -> CompList:
    params = Bunch(compsep=Bunch(nside=NSIDE, double_precision=True))
    blocks = {
        "CMB": dict(component_class="CMB", params=dict(shortname="cmb")),
        "Synchrotron": dict(component_class="Synchrotron",
                            params=dict(shortname="sync", beta=beta_sync, nu_ref=30.0)),
        "ThermalDust": dict(component_class="ThermalDust",
                            params=dict(shortname="dust", beta=beta_dust, T=20.0, nu_ref=545.0)),
    }
    comps = Bunch()
    for name, block in blocks.items():
        cfg = Bunch(enabled=True, component_class=block["component_class"], params=Bunch(
            lmax=3 * NSIDE - 1, polarization="I", spatially_varying_MM=False,
            Cl_prior_amplitude=None, **block["params"]))
        object.__setattr__(cfg, "_name", name)
        comps[name] = cfg
    return CompList.init_from_params(comps, params)


_COMPS = None


def mixing(beta_sync: float, beta_dust: float) -> np.ndarray:
    """(nband, 3) mixing matrix from Commander4's own SED implementations."""
    global _COMPS
    nu = np.array(FREQS)
    if _COMPS is None:
        comps = list(_components(TRUTH["sync"], TRUTH["dust"]).split_for_eval_pol("I"))
        _COMPS = (comps, comps[0].get_sed(nu))  # the CMB column does not depend on beta
    comps, cmb_sed = _COMPS
    comps[1].beta, comps[2].beta = beta_sync, beta_dust
    return np.stack([cmb_sed, comps[1].get_sed(nu), comps[2].get_sed(nu)], axis=1)


def make(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(2026)
    npix = hp.nside2npix(NSIDE)
    lat = np.abs(90 - np.degrees(hp.pix2ang(NSIDE, np.arange(npix))[0]))
    t = np.exp(-lat / 15.0)                     # bright towards the plane
    amps = np.stack([rng.normal(0, 70, npix),
                     10**(0.3 + 2.0 * t + 0.1 * rng.normal(size=npix)),
                     10**(0.5 + 1.5 * t + 0.1 * rng.normal(size=npix))])
    data = mixing(TRUTH["sync"], TRUTH["dust"]) @ amps
    # Every band is smoothed to one common beam, as in real per-pixel runs (compsep.common_res_fwhm);
    # Commander4 stores per-pixel amplitudes as alms, which cannot hold unsmoothed pixel-scale
    # structure. The rms map is the per-pixel rms of the smoothed noise.
    fwhm = np.radians(FWHM_ARCMIN / 60.0)
    lmax = 3 * NSIDE - 1
    for i, (nu, rms) in enumerate(zip(FREQS, RMS)):
        noise = hp.smoothing(rng.normal(0, rms, npix), fwhm=fwhm, lmax=lmax)
        sky = hp.smoothing(data[i], fwhm=fwhm, lmax=lmax) + noise
        hp.write_map(f"{outdir}/band{int(nu):03d}_map.fits", sky, column_names=["TEMPERATURE"],
                     column_units="uK_RJ", overwrite=True, dtype=np.float64)
        hp.write_map(f"{outdir}/band{int(nu):03d}_rms.fits", np.full(npix, noise.std()),
                     column_names=["TEMPERATURE"], column_units="uK_RJ", overwrite=True,
                     dtype=np.float64)
    for mode, numstep in MODES.items():
        with open(f"{outdir}/param_{mode}.yml", "w") as f:
            f.write(_param_file(outdir, mode, numstep))
    print(f"Simulation and parameter files written to {outdir}.")


def _param_file(outdir: str, mode: str, numstep: int) -> str:
    bands = ""
    for nu in FREQS:
        name = f"Band{int(nu):03d}"
        bands += f"""
    {name}:
      enabled: true
      get_from: "file"
      eval_nside: {NSIDE}
      signal_map:
        path: "{outdir}/band{int(nu):03d}_map.fits"
        dataset_names: ["TEMPERATURE", "", ""]
      rms_map:
        type: "rms"
        path: "{outdir}/band{int(nu):03d}_rms.fits"
        dataset_names: ["TEMPERATURE", "", ""]
      fwhm: {FWHM_ARCMIN}
      freq: {nu}
      polarization: "I"
      units: "uK_RJ"
"""
    return f"""gibbs:
  num_iterations: {NITER}
  seed: 7
resources:
  tod:
    num_threads: 1
  compsep:
    num_threads: 1
output:
  dir: "{outdir}/chains_{mode}/"
  logging:
    console:
      level: 'warning'
    file:
      level: 'info'
  profiling: False
  chains:
    write: [1,]
    interval:
      bands: 1
      compsep: 1
      maps: 1000
    maps_nside: "native"
    include:
      sky_model_maps: False
tod_processing: {{}}
compsep:
  double_precision: true
  per_pixel_sampling_groups:
    sample_amps_perpix:
      enabled: true
  mcmc_sampling_groups:
    sample_spectral_indices:
      enabled: true
      parameters: spectral_indices
      comps: [Synchrotron, ThermalDust]
      update_amplitude_groups: [sample_amps_perpix]
      likelihood: {mode}
      numstep: {numstep}
  bands:{bands}
components:
  CMB:
    enabled: true
    component_class: "CMB"
    params:
      polarization: "I"
      shortname: "cmb"
      Cl_prior_amplitude: null
      lmax: {3 * NSIDE - 1}
      spatially_varying_MM: False
  Synchrotron:
    enabled: true
    component_class: "Synchrotron"
    params:
      polarization: "I"
      shortname: "sync"
      Cl_prior_amplitude: null
      beta: {START["sync"]}
      nu_ref: 30
      lmax: {3 * NSIDE - 1}
      spatially_varying_MM: False
      sample_spectral_index: True
      spectral_index_proposal_sigma: {PROPOSAL["sync"]}
      spectral_index_bounds: {BOUNDS["sync"]}
  ThermalDust:
    enabled: true
    component_class: "ThermalDust"
    params:
      polarization: "I"
      shortname: "dust"
      Cl_prior_amplitude: null
      beta: {START["dust"]}
      T: 20
      nu_ref: 545
      lmax: {3 * NSIDE - 1}
      spatially_varying_MM: False
      sample_spectral_index: True
      spectral_index_proposal_sigma: {PROPOSAL["dust"]}
      spectral_index_bounds: {BOUNDS["dust"]}
experiments: {{}}
"""


def _exact_posterior(outdir: str, use_logdet: bool = True):
    """Exact 2D posterior of (beta_sync, beta_dust) on a grid, flat prior inside the bounds.

    use_logdet=True is the flat-amplitude-prior posterior ('marginal', 'conditional'); False is
    the ridge target ('ridge').
    """
    d = np.array([hp.read_map(f"{outdir}/band{int(nu):03d}_map.fits") for nu in FREQS])
    w = np.array([1.0 / hp.read_map(f"{outdir}/band{int(nu):03d}_rms.fits")**2 for nu in FREQS])

    def logp(bs, bd):
        M = mixing(bs, bd)
        A = np.einsum("bp,bi,bj->pij", w, M, M)
        b = np.einsum("bp,bi->pi", w * d, M)
        quad = np.einsum("pi,pi->p", b, np.linalg.solve(A, b[..., None])[..., 0])
        logdet = np.linalg.slogdet(A)[1] if use_logdet else 0.0
        return np.sum(0.5 * quad - 0.5 * logdet)

    # Flat prior inside the hard bounds of the parameter file. Coarse grid over the whole box, then
    # two refinements around the mass (a flat-amplitude-prior posterior can pile up at a bound).
    lo, hi = np.array(BOUNDS["sync"] + BOUNDS["dust"]).reshape(2, 2).T

    def grid_moments(lo_, hi_, n):
        gs, gd = np.linspace(lo_[0], hi_[0], n), np.linspace(lo_[1], hi_[1], n)
        L = np.array([[logp(s_, t_) for t_ in gd] for s_ in gs])
        p = np.exp(L - L.max())
        p /= p.sum()
        ps, pd = p.sum(1), p.sum(0)
        mean = np.array([np.sum(ps * gs), np.sum(pd * gd)])
        sd = np.sqrt([np.sum(ps * (gs - mean[0])**2), np.sum(pd * (gd - mean[1])**2)])
        step = np.array([gs[1] - gs[0], gd[1] - gd[0]])
        return mean, sd, step

    mean, sd, step = grid_moments(lo, hi, 41)
    for _ in range(2):
        half = np.maximum(8 * sd, 2 * step)
        mean, sd, step = grid_moments(np.maximum(lo, mean - half), np.minimum(hi, mean + half), 61)
    return mean, sd


def analyze(outdir: str) -> None:
    targets = {True: _exact_posterior(outdir, True), False: _exact_posterior(outdir, False)}
    print(f"input: beta_sync {TRUTH['sync']}, beta_dust {TRUTH['dust']}")
    for use_logdet, label in ((True, "flat-amplitude-prior posterior"), (False, "ridge target")):
        mean, sd = targets[use_logdet]
        print(f"exact {label:31s} beta_sync {mean[0]:.4f} +- {sd[0]:.4f}, "
              f"beta_dust {mean[1]:.4f} +- {sd[1]:.4f}")
    print(f"{'likelihood':14s} {'sync shift':>11s} {'sync sd ratio':>14s} "
          f"{'dust shift':>11s} {'dust sd ratio':>14s}  (vs its own exact target, in target sd)")
    for mode in MODES:
        mean, sd = targets[mode != "ridge"]
        files = sorted(glob.glob(f"{outdir}/chains_{mode}/chains_compsep/chain01_iter*.h5"))
        if len(files) <= BURN:
            print(f"{mode:14s} (not run yet, {len(files)} chain files)")
            continue
        b = []
        for fn in files[BURN:]:
            with h5py.File(fn, "r") as f:
                b.append([float(np.ravel(f["comps/sync/sed/beta"][()])[0]),
                          float(np.ravel(f["comps/dust/sed/beta"][()])[0])])
        b = np.array(b)
        shift = (b.mean(0) - mean) / sd
        ratio = b.std(0) / sd
        print(f"{mode:14s} {shift[0]:+11.2f} {ratio[0]:14.2f} {shift[1]:+11.2f} {ratio[1]:14.2f}"
              f"   chain: {b[:, 0].mean():.4f} +- {b[:, 0].std():.4f}, "
              f"{b[:, 1].mean():.5f} +- {b[:, 1].std():.5f} ({len(b)} samples)")


# ----------------------------------------------------------------------------------------------
# Mode hier_field (E018): index field + hierarchical amplitude prior, likelihood: hier
# ----------------------------------------------------------------------------------------------

HF_NSIDE = 16
HF_NSIDE_BETA = 2
HF_FREQS = [23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 94.0, 100.0, 143.0, 217.0, 353.0]
HF_RMS = [70.0, 60.0, 60.0, 55.0, 60.0, 60.0, 70.0, 40.0, 24.0, 16.0, 24.0, 80.0]  # uK_RJ, white
HF_TRUTH = {"mu": -3.0, "sigma": 0.15, "beta_dust": 1.6, "T_dust": 20.0}
HF_AMP = {"mu": 1.0269, "s": 0.2336, "cmb_sd": 70.0}   # Gaussian twin of the template scatter
HF_INDEX_PRIOR = {"mean": -3.0, "rms": 0.3}
HF_NITER = 3000
HF_WARMUP = 500
HF_BURN = HF_WARMUP + 100
HF_GROUP = "sample_index_fields"
HF_RUNS = {"hier_field": HF_NSIDE_BETA, "hier_scalar": 0}
_PROTO_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "nextgen",
                          "prototypes")


def _hf_templates(nside: int):
    """Smooth plane-brightening templates g (npix, 3): cmb 0, sync at 30 GHz, dust at 545 GHz
    (the forms of `make`, without pixel scatter), uK_RJ."""
    npix = hp.nside2npix(nside)
    lat = np.abs(90 - np.degrees(hp.pix2ang(nside, np.arange(npix))[0]))
    t = np.exp(-lat / 15.0)
    return np.stack([np.zeros(npix), 10**(0.3 + 2.0 * t), 10**(0.5 + 1.5 * t)], axis=1)


def _hf_components(nside: int) -> list:
    """Commander4 components with the prototype's mixing columns (cmb nu_ref 1e-3 GHz so that
    its column is cmb_to_rj(nu); sync pivot 30; dust 1.6, 20 K, pivot 545)."""
    params = Bunch(compsep=Bunch(nside=nside, double_precision=True))
    blocks = {
        "CMB": dict(component_class="CMB", params=dict(shortname="cmb", nu_ref=1e-3)),
        "Synchrotron": dict(component_class="Synchrotron",
                            params=dict(shortname="sync", beta=HF_TRUTH["mu"], nu_ref=30.0)),
        "ThermalDust": dict(component_class="ThermalDust",
                            params=dict(shortname="dust", beta=HF_TRUTH["beta_dust"],
                                        T=HF_TRUTH["T_dust"], nu_ref=545.0)),
    }
    comps = Bunch()
    for name, block in blocks.items():
        cfg = Bunch(enabled=True, component_class=block["component_class"], params=Bunch(
            lmax=3 * nside - 1, polarization="I", spatially_varying_MM=False,
            Cl_prior_amplitude=None, **block["params"]))
        object.__setattr__(cfg, "_name", name)
        comps[name] = cfg
    return list(CompList.init_from_params(comps, params).split_for_eval_pol("I"))


def make_hier_field(outdir: str) -> None:
    """Simulate the E018 sky: sync index field drawn from the stated GMRF prior at fixed
    hyperparameters, hier-twin amplitudes around the smooth templates, white noise; write band
    maps, rms maps, the two prior templates, truth.npz and the two parameter files."""
    from commander4.compsep import index_field as ifd
    from commander4.math_utils import gmrf_field
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(2026)
    npix = hp.nside2npix(HF_NSIDE)
    nu = np.array(HF_FREQS)
    field = gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(HF_NSIDE_BETA))
    h = field.basis.h
    rho = float(np.sqrt(h * np.pi))                      # middle of the log-uniform range
    z = rng.standard_normal(field.K)
    beta = HF_TRUTH["mu"] + HF_TRUTH["sigma"] * field.S_times(rho, z)
    region = ifd.region_map(HF_NSIDE_BETA, HF_NSIDE)
    g = _hf_templates(HF_NSIDE)
    a = g * (HF_AMP["mu"] + HF_AMP["s"] * rng.standard_normal(g.shape))
    a[:, 0] = HF_AMP["cmb_sd"] * rng.standard_normal(npix)
    comps = _hf_components(HF_NSIDE)
    M = np.stack([np.broadcast_to(comps[0].get_sed(nu), (npix, nu.size)),
                  comps[1].get_sed(nu, beta=beta[region]),
                  np.broadcast_to(comps[2].get_sed(nu), (npix, nu.size))], axis=2)
    sky = np.einsum("pbc,pc->pb", M, a)
    for i, (f, rms) in enumerate(zip(HF_FREQS, HF_RMS)):
        d = sky[:, i] + rms * rng.standard_normal(npix)
        hp.write_map(f"{outdir}/band{int(f):03d}_map.fits", d, column_names=["TEMPERATURE"],
                     column_units="uK_RJ", overwrite=True, dtype=np.float64)
        hp.write_map(f"{outdir}/band{int(f):03d}_rms.fits", np.full(npix, rms),
                     column_names=["TEMPERATURE"], column_units="uK_RJ", overwrite=True,
                     dtype=np.float64)
    hp.write_map(f"{outdir}/template_sync_30GHz.fits", g[:, 1], overwrite=True, dtype=np.float64)
    hp.write_map(f"{outdir}/template_dust_545GHz.fits", g[:, 2], overwrite=True,
                 dtype=np.float64)
    np.savez(f"{outdir}/truth.npz", beta=beta, z=z, mu=HF_TRUTH["mu"], sigma=HF_TRUTH["sigma"],
             rho=rho, a=a, g=g, region=region, nu=nu, rms=np.array(HF_RMS),
             amp_mu=HF_AMP["mu"], amp_s=HF_AMP["s"])
    for tag, nside_beta in HF_RUNS.items():
        with open(f"{outdir}/param_{tag}.yml", "w") as f:
            f.write(_hf_param_file(outdir, tag, nside_beta))
    print(f"hier_field simulation and parameter files written to {outdir} (truth: mu "
          f"{HF_TRUTH['mu']}, sigma {HF_TRUTH['sigma']}, rho {rho:.3f}, K {field.K}).")


def _hf_param_file(outdir: str, tag: str, nside_beta: int) -> str:
    bands = ""
    for f in HF_FREQS:
        name = f"Band{int(f):03d}"
        bands += f"""
    {name}:
      enabled: true
      get_from: "file"
      eval_nside: {HF_NSIDE}
      signal_map:
        path: "{outdir}/band{int(f):03d}_map.fits"
        dataset_names: ["TEMPERATURE", "", ""]
      rms_map:
        type: "rms"
        path: "{outdir}/band{int(f):03d}_rms.fits"
        dataset_names: ["TEMPERATURE", "", ""]
      fwhm: {FWHM_ARCMIN}
      freq: {f}
      polarization: "I"
      units: "uK_RJ"
"""
    field_block = ""
    if nside_beta > 0:
        field_block = f"""
      spectral_index_nside: {nside_beta}
      spectral_index_field:
        basis: graph
        alpha: 2.0
        sigma_scale: 0.3
        rho_prior: loguniform"""
    return f"""gibbs:
  num_iterations: {HF_NITER}
  seed: 7
resources:
  tod:
    num_threads: 1
  compsep:
    num_threads: 1
output:
  dir: "{outdir}/chains_{tag}/"
  logging:
    console:
      level: 'warning'
    file:
      level: 'info'
  profiling: False
  chains:
    write: [1,]
    interval:
      bands: 1
      compsep: 1
      maps: 100000
    maps_nside: "native"
    include:
      sky_model_maps: False
tod_processing: {{}}
compsep:
  double_precision: true
  per_pixel_sampling_groups:
    sample_amps_perpix:
      enabled: true
  mcmc_sampling_groups:
    {HF_GROUP}:
      enabled: true
      parameters: spectral_indices
      comps: [Synchrotron]
      update_amplitude_groups: [sample_amps_perpix]
      likelihood: hier
      numstep: 1
      nuts:
        n_warmup: {HF_WARMUP}
        target_accept: 0.8
        max_depth: 10
        parameterization: noncentred
  bands:{bands}
components:
  CMB:
    enabled: true
    component_class: "CMB"
    params:
      polarization: "I"
      shortname: "cmb"
      nu_ref: 0.001
      Cl_prior_amplitude: null
      lmax: {3 * HF_NSIDE - 1}
      spatially_varying_MM: False
      amp_prior:
        type: gaussian
        mean: 0.0
        rms: {HF_AMP["cmb_sd"]}
  Synchrotron:
    enabled: true
    component_class: "Synchrotron"
    params:
      polarization: "I"
      shortname: "sync"
      Cl_prior_amplitude: null
      beta: {HF_TRUTH["mu"]}
      nu_ref: 30
      lmax: {3 * HF_NSIDE - 1}
      spatially_varying_MM: False
      sample_spectral_index: True
      spectral_index_prior:
        type: gaussian
        mean: {HF_INDEX_PRIOR["mean"]}
        rms: {HF_INDEX_PRIOR["rms"]}{field_block}
      amp_prior:
        type: hier
        template: "{outdir}/template_sync_30GHz.fits"
        template_floor: 1.0e-3
        mu_prior: [1.0, 0.5]
        log_s_prior: [{float(np.log(0.3))}, 1.0]
  ThermalDust:
    enabled: true
    component_class: "ThermalDust"
    params:
      polarization: "I"
      shortname: "dust"
      Cl_prior_amplitude: null
      beta: {HF_TRUTH["beta_dust"]}
      T: {HF_TRUTH["T_dust"]}
      nu_ref: 545
      lmax: {3 * HF_NSIDE - 1}
      spatially_varying_MM: False
      amp_prior:
        type: hier
        template: "{outdir}/template_dust_545GHz.fits"
        template_floor: 1.0e-3
        mu_prior: [1.0, 0.5]
        log_s_prior: [{float(np.log(0.3))}, 1.0]
experiments: {{}}
"""


def _hf_read_chain(outdir: str, tag: str) -> dict | None:
    """Per-iteration samples of the hier run `tag` after HF_BURN: beta per index pixel (n, Nb) or
    (n, 1), field hyperparameters (field run), amplitude hyperparameters, NUTS diagnostics."""
    files = sorted(glob.glob(f"{outdir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))
    if len(files) <= HF_BURN:
        print(f"{tag}: not run yet or too short ({len(files)} chain files, burn {HF_BURN})")
        return None
    out = {k: [] for k in ("beta", "mu", "log_sigma", "u", "mu_s", "log_s_s", "mu_d", "log_s_d",
                           "divergent", "n_leapfrog", "step_size", "accept_stat")}
    for fn in files[HF_BURN:]:
        with h5py.File(fn, "r") as f:
            st = f[f"mcmc/{HF_GROUP}/state"]
            nuts = f[f"mcmc/{HF_GROUP}/nuts"]
            if "index_field" in st["Synchrotron"]:
                fld = st["Synchrotron/index_field"]
                out["beta"].append(np.asarray(fld["beta_map"][()]))
                for k in ("mu", "log_sigma", "u"):
                    out[k].append(float(fld[k][()]))
            else:
                out["beta"].append(np.atleast_1d(np.asarray(st["Synchrotron/beta"][()])))
            out["mu_s"].append(float(st["Synchrotron/amp_prior/mu"][()]))
            out["log_s_s"].append(float(st["Synchrotron/amp_prior/log_s"][()]))
            out["mu_d"].append(float(st["ThermalDust/amp_prior/mu"][()]))
            out["log_s_d"].append(float(st["ThermalDust/amp_prior/log_s"][()]))
            for k in ("divergent", "n_leapfrog", "step_size", "accept_stat"):
                out[k].append(float(nuts[k][()]))
    return {k: np.array(v) for k, v in out.items() if len(v)}


def _hf_prototype_chain(outdir: str, nside_beta: int, rng: np.random.Generator) -> dict:
    """The E016 prototype (HierNonCentredPosterior or ScalarHierPosterior) sampled with
    sampler.run_nuts on the same data; returns arrays as _hf_read_chain."""
    sys.path.insert(0, _PROTO_DIR)
    from indexfield import field as pfield, posterior, sampler, scalar_coverage
    from indexfield.model import PixelData
    nu = np.array(HF_FREQS)
    d = np.array([hp.read_map(f"{outdir}/band{int(f):03d}_map.fits") for f in HF_FREQS]).T
    rms = np.array([hp.read_map(f"{outdir}/band{int(f):03d}_rms.fits") for f in HF_FREQS]).T
    g = np.stack([np.ones(d.shape[0]), hp.read_map(f"{outdir}/template_sync_30GHz.fits"),
                  hp.read_map(f"{outdir}/template_dust_545GHz.fits")], axis=1)
    g[:, 1:] = np.maximum(g[:, 1:], 1e-3)
    if nside_beta > 0:
        from commander4.compsep import index_field as ifd
        region = ifd.region_map(nside_beta, HF_NSIDE)
        n_beta = 12 * nside_beta**2
    else:
        region = np.zeros(d.shape[0], dtype=np.int64)
        n_beta = 1
    data = PixelData(nu=nu, d=d, inv_var=1.0 / rms**2, region=region, n_beta=n_beta,
                     beta_d=HF_TRUTH["beta_dust"], T_d=HF_TRUTH["T_dust"])
    out = {}
    if nside_beta > 0:
        fld = pfield.MaternField(pfield.GraphLaplacianBasis(nside_beta))
        hyper = posterior.HyperPrior(mu0=HF_INDEX_PRIOR["mean"], s_mu=HF_INDEX_PRIOR["rms"])
        post = posterior.HierNonCentredPosterior(fld, hyper, data, g)
        res = sampler.run_nuts(post, rng, HF_WARMUP, HF_NITER - HF_BURN, max_depth=10)
        K = post.K
        x = res.samples
        out["beta"] = np.array([post.beta(xi[:K + 3]) for xi in x])
        out["mu"], out["log_sigma"], out["u"] = x[:, K], x[:, K + 1], x[:, K + 2]
        out["mu_s"], out["mu_d"] = x[:, K + 3], x[:, K + 4]
        out["log_s_s"], out["log_s_d"] = x[:, K + 5], x[:, K + 6]
    else:
        post = scalar_coverage.ScalarHierPosterior(data, g)
        res = sampler.run_nuts(post, rng, HF_WARMUP, HF_NITER - HF_BURN, max_depth=10)
        x = res.samples
        out["beta"] = x[:, :1]
        out["mu_s"], out["mu_d"], out["log_s_s"], out["log_s_d"] = x[:, 1], x[:, 2], x[:, 3], x[:, 4]
    out["n_divergent"] = res.n_divergent
    out["n_leapfrog"] = res.n_leapfrog
    return out


def _ess(x: np.ndarray) -> float:
    """Geyer initial positive sequence ESS of a 1-D chain."""
    x = np.asarray(x, dtype=float) - np.mean(x)
    n = x.size
    var = np.mean(x * x)
    if var == 0:
        return float(n)
    acf = np.correlate(x, x, "full")[n - 1:] / (np.arange(n, 0, -1) * var)
    s = 0.0
    for k in range(0, n - 1, 2):
        pair = acf[k] + (acf[k + 1] if k + 1 < n else 0.0)
        if pair <= 0 and k > 0:
            break
        s += pair
    return float(n / max(2 * s - 1.0, 1.0))


def _hf_compare(name: str, ours: np.ndarray, theirs: np.ndarray, truth=None) -> dict:
    """Mean difference in pooled MC SE, sd ratio, ESS of both chains, and z against the truth."""
    ours = np.atleast_2d(ours.T).T
    theirs = np.atleast_2d(theirs.T).T
    rows = []
    for j in range(ours.shape[1]):
        e1, e2 = _ess(ours[:, j]), _ess(theirs[:, j])
        se = np.sqrt(ours[:, j].var() / e1 + theirs[:, j].var() / e2)
        rows.append(dict(diff_se=float((ours[:, j].mean() - theirs[:, j].mean()) / se),
                         sd_ratio=float(ours[:, j].std() / theirs[:, j].std()),
                         ess_c4=float(e1), ess_proto=float(e2),
                         mean_c4=float(ours[:, j].mean()), sd_c4=float(ours[:, j].std()),
                         mean_proto=float(theirs[:, j].mean()),
                         z_truth=(float((ours[:, j].mean() - truth[j]) / ours[:, j].std())
                                  if truth is not None else None)))
    ds = np.array([r["diff_se"] for r in rows])
    sr = np.array([r["sd_ratio"] for r in rows])
    ok = bool(np.all(np.abs(ds) < 4) and np.all((sr > 0.85) & (sr < 1.15)))
    print(f"  {name:10s} n {len(rows):3d}  |diff|/se max {np.abs(ds).max():5.2f} rms "
          f"{np.sqrt(np.mean(ds**2)):5.2f}  sd ratio {sr.min():.3f} to {sr.max():.3f}  ESS c4 min "
          f"{min(r['ess_c4'] for r in rows):5.0f} proto min {min(r['ess_proto'] for r in rows):5.0f}"
          + (f"  z(truth) rms {np.sqrt(np.mean([r['z_truth']**2 for r in rows])):.2f}"
             if truth is not None else "") + ("  PASS" if ok else "  FAIL"))
    return dict(rows=rows, pass_=ok)


def analyze_hier_field(outdir: str) -> None:
    truth = np.load(f"{outdir}/truth.npz")
    report = {}
    for tag, nside_beta in HF_RUNS.items():
        c4 = _hf_read_chain(outdir, tag)
        if c4 is None:
            continue
        n = c4["beta"].shape[0]
        print(f"\n{tag}: {n} Commander4 samples after burn {HF_BURN}; divergences "
              f"{int(c4['divergent'].sum())} ({100 * c4['divergent'].mean():.2f} percent), mean "
              f"gradients per iteration {c4['n_leapfrog'].mean():.1f}, step size "
              f"{c4['step_size'][-1]:.3g}, mean accept {c4['accept_stat'].mean():.3f}")
        proto = _hf_prototype_chain(outdir, nside_beta, np.random.default_rng(11))
        print(f"  prototype chain: {proto['beta'].shape[0]} samples, divergences "
              f"{proto['n_divergent']}, gradients {proto['n_leapfrog']}")
        res = {}
        beta_truth = truth["beta"] if nside_beta > 0 else np.array([truth["mu"]])
        res["beta"] = _hf_compare("beta", c4["beta"], proto["beta"], beta_truth)
        if nside_beta > 0:
            res["mu"] = _hf_compare("mu", c4["mu"], proto["mu"], [float(truth["mu"])])
            res["log_sigma"] = _hf_compare("log_sigma", c4["log_sigma"], proto["log_sigma"],
                                           [float(np.log(truth["sigma"]))])
            res["u"] = _hf_compare("u", c4["u"], proto["u"])
        for k, tv in (("mu_s", truth["amp_mu"]), ("mu_d", truth["amp_mu"]),
                      ("log_s_s", np.log(truth["amp_s"])), ("log_s_d", np.log(truth["amp_s"]))):
            res[k] = _hf_compare(k, c4[k], proto[k], [float(tv)])
        res["divergent_fraction"] = float(c4["divergent"].mean())
        res["pass_"] = bool(all(v["pass_"] for v in res.values() if isinstance(v, dict))
                            and res["divergent_fraction"] < 0.01)
        print(f"  {tag}: {'PASS' if res['pass_'] else 'FAIL'} (all |diff| < 4 se, sd ratios within "
              "0.85 to 1.15, divergences < 1 percent)")
        report[tag] = res
    with open(f"{outdir}/analysis_hier_field.json", "w") as f:
        json.dump(report, f, indent=1, default=float)
    print(f"\nwritten {outdir}/analysis_hier_field.json")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["make", "analyze"])
    parser.add_argument("--dir", required=True)
    parser.add_argument("--mode", choices=["mh", "hier_field"], default="mh",
                        help="mh: the four MH likelihood modes (default); hier_field: E018")
    args = parser.parse_args()
    if args.mode == "hier_field":
        make_hier_field(args.dir) if args.action == "make" else analyze_hier_field(args.dir)
    else:
        make(args.dir) if args.action == "make" else analyze(args.dir)
