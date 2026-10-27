"""E020 simulation: PySM d1 + s1 + c1 intensity sky at nside 64 in the 12 E018 bands, white noise,
pseudo-templates for the hier amplitude prior built from the data, truth maps, and one Commander4
parameter file per index-prior option (gmrf / independent at nside_beta 4, 8, 16; scalar).

Truth stored in truth.npz: the PySM parameter maps degraded to nside (pl_index, mbb_index,
mbb_temperature), the c1 CMB map in uK_CMB, the noiseless band maps, the noise rms per band.

    python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e020_d1s1 --nside 64 \
        --noise-scale 1
"""
from __future__ import annotations

import argparse
import os

import healpy as hp
import numpy as np

FREQS = [23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 94.0, 100.0, 143.0, 217.0, 353.0]
RMS16 = [70.0, 60.0, 60.0, 55.0, 60.0, 60.0, 70.0, 40.0, 24.0, 16.0, 24.0, 80.0]  # uK_RJ at nside 16
INDEX_PRIOR = {"sync": (-3.0, 0.3), "dust": (1.55, 0.3)}
CMB_SD = 70.0
RUNS = {"gmrf4": ("gmrf", 4), "gmrf8": ("gmrf", 8), "gmrf16": ("gmrf", 16),
        "indep4": ("independent", 4), "indep8": ("independent", 8),
        "indep16": ("independent", 16), "scalar": (None, 0)}
H_OVER_K = 0.0479924466  # K / GHz


def _bb_ratio(nu, nu0, T):
    """B(nu, T) / B(nu0, T) up to the (nu/nu0)^3 factor, in the RJ convention used by
    Commander4's ThermalDust: (nu/nu0)^(beta + 1) expm1(x0) / expm1(x)."""
    x, x0 = H_OVER_K * nu / T, H_OVER_K * nu0 / T
    return np.expm1(x0) / np.expm1(x)


def simulate(outdir: str, nside: int, noise_scale: float, seed: int, n_iter: int,
             n_warmup: int) -> None:
    import pysm3
    import pysm3.units as u
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    npix = hp.nside2npix(nside)
    fg = pysm3.Sky(nside=nside, preset_strings=["d1", "s1"])
    cmb = pysm3.Sky(nside=nside, preset_strings=["c1"])
    dust, sync = fg.components
    beta_d = np.asarray(getattr(dust.mbb_index, "value", dust.mbb_index), dtype=float)
    T_d = np.asarray(dust.mbb_temperature.to(u.K).value, dtype=float)
    beta_s = np.asarray(getattr(sync.pl_index, "value", sync.pl_index), dtype=float)
    cmb_map = cmb.get_emission(100.0 * u.GHz)[0].to(
        u.uK_CMB, equivalencies=u.cmb_equivalencies(100.0 * u.GHz)).value
    T_fixed = float(T_d.mean())
    rms = np.array(RMS16) * (nside / 16.0) * noise_scale
    sky = np.empty((npix, len(FREQS)))
    for i, f in enumerate(FREQS):
        nu = f * u.GHz
        tot = fg.get_emission(nu)[0] + cmb.get_emission(nu)[0]
        sky[:, i] = tot.to(u.uK_RJ, equivalencies=u.cmb_equivalencies(nu)).value
    data = sky + rms[None, :] * rng.standard_normal(sky.shape)
    for i, f in enumerate(FREQS):
        hp.write_map(f"{outdir}/band{int(f):03d}_map.fits", data[:, i],
                     column_names=["TEMPERATURE"], column_units="uK_RJ", overwrite=True,
                     dtype=np.float64)
        hp.write_map(f"{outdir}/band{int(f):03d}_rms.fits", np.full(npix, rms[i]),
                     column_names=["TEMPERATURE"], column_units="uK_RJ", overwrite=True,
                     dtype=np.float64)
    # Pseudo-templates from the data: 23 GHz scaled to 30 GHz with beta -3; 353 GHz scaled to
    # 545 GHz with beta 1.55 and T_fixed. They carry CMB, noise, the other foreground and the
    # wrong constant index; the hier prior learns their scale and scatter.
    g_sync = data[:, FREQS.index(23.0)] * (30.0 / 23.0) ** INDEX_PRIOR["sync"][0]
    g_dust = data[:, FREQS.index(353.0)] * (545.0 / 353.0) ** (INDEX_PRIOR["dust"][0] + 1.0) \
        * _bb_ratio(545.0, 353.0, T_fixed)
    for name, g in (("template_sync_30GHz", g_sync), ("template_dust_545GHz", g_dust)):
        hp.write_map(f"{outdir}/{name}.fits", g, overwrite=True, dtype=np.float64)
    np.savez(f"{outdir}/truth.npz", beta_s=beta_s, beta_d=beta_d, T_d=T_d, T_fixed=T_fixed,
             cmb=cmb_map, sky=sky, rms=rms, nu=np.array(FREQS), nside=nside,
             noise_scale=noise_scale, seed=seed, g_sync=g_sync, g_dust=g_dust)
    for tag, (kind, nside_beta) in RUNS.items():
        with open(f"{outdir}/param_{tag}.yml", "w") as f:
            f.write(param_file(outdir, tag, kind, nside_beta, nside, T_fixed, n_iter, n_warmup))
    print(f"E020 sky written to {outdir}: nside {nside}, noise scale {noise_scale}, T_fixed "
          f"{T_fixed:.3f} K; truth sd beta_s {beta_s.std():.4f}, beta_d {beta_d.std():.4f}, "
          f"cmb rms {cmb_map.std():.1f} uK_CMB; {len(RUNS)} parameter files.")


def _field_block(kind: str | None, nside_beta: int) -> str:
    if kind is None:
        return ""
    if kind == "independent":
        return f"""
      spectral_index_nside: {nside_beta}
      spectral_index_field:
        type: independent"""
    return f"""
      spectral_index_nside: {nside_beta}
      spectral_index_field:
        type: gmrf
        basis: graph
        alpha: 2.0
        sigma_scale: 0.3
        rho_prior: loguniform"""


def param_file(outdir: str, tag: str, kind: str | None, nside_beta: int, nside: int,
               T_fixed: float, n_iter: int, n_warmup: int) -> str:
    fwhm = 2.0 * hp.nside2resol(nside, arcmin=True)
    bands = ""
    for f in FREQS:
        bands += f"""
    Band{int(f):03d}:
      enabled: true
      get_from: "file"
      eval_nside: {nside}
      signal_map:
        path: "{outdir}/band{int(f):03d}_map.fits"
        dataset_names: ["TEMPERATURE", "", ""]
      rms_map:
        type: "rms"
        path: "{outdir}/band{int(f):03d}_rms.fits"
        dataset_names: ["TEMPERATURE", "", ""]
      fwhm: {fwhm:.1f}
      freq: {f}
      polarization: "I"
      units: "uK_RJ"
"""
    lmax = 3 * nside - 1
    log_s0 = float(np.log(0.3))
    return f"""gibbs:
  num_iterations: {n_iter}
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
    sample_index_fields:
      enabled: true
      parameters: spectral_indices
      comps: [Synchrotron, ThermalDust]
      update_amplitude_groups: [sample_amps_perpix]
      likelihood: hier
      numstep: 1
      nuts:
        n_warmup: {n_warmup}
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
      lmax: {lmax}
      spatially_varying_MM: False
      amp_prior:
        type: gaussian
        mean: 0.0
        rms: {CMB_SD}
  Synchrotron:
    enabled: true
    component_class: "Synchrotron"
    params:
      polarization: "I"
      shortname: "sync"
      Cl_prior_amplitude: null
      beta: {INDEX_PRIOR["sync"][0]}
      nu_ref: 30
      lmax: {lmax}
      spatially_varying_MM: False
      sample_spectral_index: True
      spectral_index_prior:
        type: gaussian
        mean: {INDEX_PRIOR["sync"][0]}
        rms: {INDEX_PRIOR["sync"][1]}{_field_block(kind, nside_beta)}
      amp_prior:
        type: hier
        template: "{outdir}/template_sync_30GHz.fits"
        template_floor: 1.0e-3
        mu_prior: [1.0, 0.5]
        log_s_prior: [{log_s0}, 1.0]
  ThermalDust:
    enabled: true
    component_class: "ThermalDust"
    params:
      polarization: "I"
      shortname: "dust"
      Cl_prior_amplitude: null
      beta: {INDEX_PRIOR["dust"][0]}
      T: {T_fixed:.4f}
      nu_ref: 545
      lmax: {lmax}
      spatially_varying_MM: False
      sample_spectral_index: True
      spectral_index_prior:
        type: gaussian
        mean: {INDEX_PRIOR["dust"][0]}
        rms: {INDEX_PRIOR["dust"][1]}{_field_block(kind, nside_beta)}
      amp_prior:
        type: hier
        template: "{outdir}/template_dust_545GHz.fits"
        template_floor: 1.0e-3
        mu_prior: [1.0, 0.5]
        log_s_prior: [{log_s0}, 1.0]
experiments: {{}}
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--nside", type=int, default=64)
    ap.add_argument("--noise-scale", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=2027)
    ap.add_argument("--n-iter", type=int, default=2000)
    ap.add_argument("--n-warmup", type=int, default=500)
    a = ap.parse_args(argv)
    simulate(a.dir, a.nside, a.noise_scale, a.seed, a.n_iter, a.n_warmup)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
