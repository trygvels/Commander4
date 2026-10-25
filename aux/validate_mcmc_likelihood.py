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
"""
import argparse
import glob
import os

import h5py
import healpy as hp
import numpy as np
from pixell.bunch import Bunch

from commander4.sky.comp_list import CompList

FREQS = [30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0]
RMS = [60.0, 60.0, 70.0, 24.0, 16.0, 24.0, 80.0]  # uK_RJ per pixel, before smoothing
TRUTH = {"sync": -3.0, "dust": 1.6}
START = {"sync": -2.9, "dust": 1.595}      # chains start several sd away from the truth
PROPOSAL = {"sync": 0.02, "dust": 0.0007}  # about the exact posterior sd of each index
NSIDE = 16
FWHM_ARCMIN = 480.0                        # common beam, about two pixels
MODES = {"marginal": 20, "ridge": 20, "conditional": 20, "resolve_chisq": 20}  # -> numstep
NITER = 300
BURN = 100


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


def mixing(beta_sync: float, beta_dust: float) -> np.ndarray:
    """(nband, 3) mixing matrix from Commander4's own SED implementations."""
    comps = list(_components(beta_sync, beta_dust).split_for_eval_pol("I"))
    nu = np.array(FREQS)
    return np.stack([c.get_sed(nu) for c in comps], axis=1)


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
      spectral_index_bounds: [-4.0, -2.0]
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
      spectral_index_bounds: [1.0, 2.5]
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

    centre = np.array([TRUTH["sync"], TRUTH["dust"]])
    half = np.array([0.4, 0.15])
    for _ in range(3):  # zoom in twice around the posterior
        gs = np.linspace(centre[0] - half[0], centre[0] + half[0], 61)
        gd = np.linspace(centre[1] - half[1], centre[1] + half[1], 61)
        L = np.array([[logp(s, t) for t in gd] for s in gs])
        p = np.exp(L - L.max())
        p /= p.sum()
        ps, pd = p.sum(1), p.sum(0)
        mean = np.array([np.sum(ps * gs), np.sum(pd * gd)])
        sd = np.sqrt([np.sum(ps * (gs - mean[0])**2), np.sum(pd * (gd - mean[1])**2)])
        centre, half = mean, 8 * sd
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
        if not files:
            print(f"{mode:14s} (no chain files)")
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
              f"   ({len(b)} samples)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("action", choices=["make", "analyze"])
    parser.add_argument("--dir", required=True)
    args = parser.parse_args()
    make(args.dir) if args.action == "make" else analyze(args.dir)
