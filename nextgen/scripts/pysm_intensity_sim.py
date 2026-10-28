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
# Index prior means offset from the PySM sky means (-3.000, 1.540) by more than one truth sd
# (0.063, 0.037), so that shrinkage towards mu0 cannot reproduce the truth by construction
# (stats review of E020, item 3); sd 0.3 covers both truth ranges.
INDEX_PRIOR = {"sync": (-3.10, 0.3), "dust": (1.60, 0.3)}
INDEP_SD = 1.0            # per-region prior sd of the 'independent' option (flat where data speak)
CMB_SD = 100.0            # a_cmb ~ N(0, CMB_SD^2); the c1 map rms at nside 64 is 84 uK
TEMPLATE_FWHM_DEG = 5.0   # smoothing of the data-built pseudo-templates
TEMPLATE_FLOOR_PCT = 1.0  # floor at this percentile of the smoothed template
RUNS = {"gmrf4": ("gmrf", 4), "gmrf8": ("gmrf", 8), "gmrf16": ("gmrf", 16),
        "iid4": ("iid", 4), "iid8": ("iid", 8), "iid16": ("iid", 16),
        "indep4": ("independent", 4), "indep8": ("independent", 8),
        "indep16": ("independent", 16), "scalar": (None, 0)}
H_OVER_K = 0.0479924466  # K / GHz


def _dust_sed(nu, beta, T):
    """(..., nband) ThermalDust mixing column at pivot 545 GHz for beta (...,) and scalar T."""
    nu = np.asarray(nu, dtype=float)
    beta = np.asarray(beta, dtype=float)
    x, x0 = H_OVER_K * nu / T, H_OVER_K * 545.0 / T
    return (nu / 545.0) ** (beta[..., None] + 1.0) * (np.expm1(x0) / np.expm1(x))


def fisher_sd(a_true, dsed, rms, region, nreg):
    """Per-region Fisher sd of an index with the amplitudes known: 1/sqrt(sum_{p in r} sum_b
    (a_p dsed_pb / rms_b)^2)."""
    f = np.sum((a_true[:, None] * dsed / rms[None, :]) ** 2, axis=1)
    F = np.bincount(region, f, minlength=nreg)
    with np.errstate(divide="ignore"):
        return 1.0 / np.sqrt(F)


def fisher_sd_profiled(M, dM, rms, region, nreg):
    """Per-region Fisher sd of (beta_s, beta_d) with the three per-pixel amplitudes profiled
    out: per pixel the 5 x 5 Fisher of (a_cmb, a_s, a_d, beta_s, beta_d) at the truth, Schur
    complement onto the two indices, summed over the region, inverted. M (npix, B, 3) is the
    mixing at the truth, dM (npix, B, 2) the derivative columns a_c dsed_c/dbeta_c."""
    w = 1.0 / rms ** 2
    J = np.concatenate([M, dM], axis=2)                       # (npix, B, 5)
    F = np.einsum("pbi,b,pbj->pij", J, w, J)                  # (npix, 5, 5)
    Faa, Fab, Fbb = F[:, :3, :3], F[:, :3, 3:], F[:, 3:, 3:]
    S = Fbb - np.einsum("pij,pjk->pik", np.transpose(Fab, (0, 2, 1)), np.linalg.solve(Faa, Fab))
    Sr = np.zeros((nreg, 2, 2))
    np.add.at(Sr, region, S)
    out = np.full((nreg, 2), np.inf)
    ok = np.linalg.det(Sr) > 0
    out[ok] = np.sqrt(np.diagonal(np.linalg.inv(Sr[ok]), axis1=1, axis2=2))
    return out[:, 0], out[:, 1]


def beta_eff_fixed_T(dust_bands, nu, rms, T_fixed, grid=np.linspace(1.0, 2.2, 601)):
    """Noise-free best-fit dust index per pixel with T fixed (rms-weighted least squares over the
    bands, dust amplitude profiled, dust emission alone): the two-parameter effective index."""
    npix = dust_bands.shape[0]
    out = np.empty(npix)
    w = 1.0 / rms ** 2
    for lo in range(0, npix, 4096):
        d = dust_bands[lo:lo + 4096]                          # (n, B)
        sed = _dust_sed(nu, grid, T_fixed)                    # (G, B)
        num = d @ (sed * w[None, :]).T                        # (n, G): sum_b w d s
        den = np.sum(sed ** 2 * w[None, :], axis=1)           # (G,)
        chi2 = np.sum(d ** 2 * w[None, :], axis=1)[:, None] - num ** 2 / den[None, :]
        out[lo:lo + 4096] = grid[np.argmin(chi2, axis=1)]
    return out


def beta_eff_fixed_T_3comp(sky, nu, rms, beta_s, T_fixed, cmb_col,
                           grid=np.linspace(1.0, 2.2, 601)):
    """Noise-free best-fit dust index per pixel with T fixed and the cmb, sync and dust
    amplitudes profiled (three-column GLS at the true beta_s on the full noiseless sky): the
    effective index the model fit actually converges to (E020 review 2, item 4)."""
    npix = sky.shape[0]
    w = 1.0 / rms ** 2
    sync_col = (nu[None, :] / 30.0) ** beta_s[:, None]
    best = np.full(npix, np.inf)
    out = np.empty(npix)
    cmb_b = np.broadcast_to(cmb_col, (npix, nu.size))
    for b in grid:
        dust_col = np.broadcast_to(_dust_sed(nu, np.array([b]), T_fixed), (npix, nu.size))
        M = np.stack([cmb_b, sync_col, dust_col], axis=2)
        A = np.einsum("pbi,b,pbj->pij", M, w, M)
        rhs = np.einsum("pbi,b,pb->pi", M, w, sky)
        a = np.linalg.solve(A, rhs[..., None])[..., 0]
        resid = sky - np.einsum("pbi,pi->pb", M, a)
        chi2 = np.sum(resid ** 2 * w[None, :], axis=1)
        better = chi2 < best
        best[better] = chi2[better]
        out[better] = b
    return out


def cmb_sed(nu):
    """uK_CMB to uK_RJ conversion per band (the CMB mixing column)."""
    x = H_OVER_K * np.asarray(nu, dtype=float) / 2.7255
    return x ** 2 * np.exp(x) / np.expm1(x) ** 2


def _bb_ratio(nu, nu0, T):
    """B(nu, T) / B(nu0, T) up to the (nu/nu0)^3 factor, in the RJ convention used by
    Commander4's ThermalDust: (nu/nu0)^(beta + 1) expm1(x0) / expm1(x)."""
    x, x0 = H_OVER_K * nu / T, H_OVER_K * nu0 / T
    return np.expm1(x0) / np.expm1(x)


def simulate(outdir: str, nside: int, noise_scale: float, seed: int, n_iter: int,
             n_warmup: int, sky_mode: str = "pysm", nside_truth: int = 16,
             template_scatter: float = 0.3, adapt_mass: bool = True) -> None:
    """`sky_mode` "pysm": the PySM d1 s1 c1 sky as it is (E020; T_d varies, indices vary per nside
    pixel, templates built from the data). "consistent" (E021): a sky drawn from the model class
    Commander4 fits: indices constant per nside_truth region (the PySM region means), one modified
    blackbody at T_fixed, PySM amplitudes at the pivots and the c1 CMB, and templates equal to the
    true amplitudes divided by (1 + template_scatter * eps), eps iid N(0, 1), the ratio floored at
    0.2, so that a_true = g (1 + template_scatter eps) matches the hier prior a ~ N(mu g, (s g)^2)
    with mu = 1, s = template_scatter up to the floor (hit by 0.4 percent of pixels)."""
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
    nu_arr = np.array(FREQS)
    # Truth amplitudes at the model pivots (uK_RJ) and the noise-free dust emission per band.
    a_sync = sync.get_emission(30.0 * u.GHz)[0].to(
        u.uK_RJ, equivalencies=u.cmb_equivalencies(30.0 * u.GHz)).value
    a_dust = dust.get_emission(545.0 * u.GHz)[0].to(
        u.uK_RJ, equivalencies=u.cmb_equivalencies(545.0 * u.GHz)).value
    sky = np.empty((npix, len(FREQS)))
    dust_bands = np.empty((npix, len(FREQS)))
    if sky_mode == "consistent":
        # Region-constant indices (plain mean of the PySM map over each nside_truth region), one
        # MBB at T_fixed; the amplitudes at the pivots are the PySM ones, so the sky differs from
        # PySM only through the SED shapes.
        beta_s = hp.ud_grade(hp.ud_grade(beta_s, nside_truth), nside)
        beta_d = hp.ud_grade(hp.ud_grade(beta_d, nside_truth), nside)
        T_d = np.full(npix, T_fixed)
        dust_bands[:] = a_dust[:, None] * _dust_sed(nu_arr, beta_d, T_fixed)
        sky[:] = (cmb_sed(nu_arr)[None, :] * cmb_map[:, None]
                  + a_sync[:, None] * (nu_arr[None, :] / 30.0) ** beta_s[:, None] + dust_bands)
    elif sky_mode == "pysm":
        for i, f in enumerate(FREQS):
            nu = f * u.GHz
            tot = fg.get_emission(nu)[0] + cmb.get_emission(nu)[0]
            sky[:, i] = tot.to(u.uK_RJ, equivalencies=u.cmb_equivalencies(nu)).value
            dust_bands[:, i] = dust.get_emission(nu)[0].to(
                u.uK_RJ, equivalencies=u.cmb_equivalencies(nu)).value
    else:
        raise ValueError(f"sky_mode must be 'pysm' or 'consistent', got {sky_mode!r}")
    data = sky + rms[None, :] * rng.standard_normal(sky.shape)
    for i, f in enumerate(FREQS):
        hp.write_map(f"{outdir}/band{int(f):03d}_map.fits", data[:, i],
                     column_names=["TEMPERATURE"], column_units="uK_RJ", overwrite=True,
                     dtype=np.float64)
        hp.write_map(f"{outdir}/band{int(f):03d}_rms.fits", np.full(npix, rms[i]),
                     column_names=["TEMPERATURE"], column_units="uK_RJ", overwrite=True,
                     dtype=np.float64)
    # Pseudo-templates from the data (not the truth): 23 GHz scaled to 30 GHz and 353 GHz scaled
    # to 545 GHz with the prior-mean indices and T_fixed, smoothed to TEMPLATE_FWHM_DEG and
    # floored at the TEMPLATE_FLOOR_PCT percentile of the smoothed map. They carry CMB, noise,
    # the other foreground and a wrong constant index; the hier prior learns scale and scatter.
    lmax = 3 * nside - 1
    g_sync = data[:, FREQS.index(23.0)] * (30.0 / 23.0) ** INDEX_PRIOR["sync"][0]
    g_dust = data[:, FREQS.index(353.0)] * (545.0 / 353.0) ** (INDEX_PRIOR["dust"][0] + 1.0) \
        * _bb_ratio(545.0, 353.0, T_fixed)
    lat = hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]
    hi = np.abs(lat) > 20
    report = []
    templates = {}
    for name, g, a_true in (("sync", g_sync, a_sync), ("dust", g_dust, a_dust)):
        if sky_mode == "consistent":
            # Template = truth / (1 + s eps) with the ratio floored at 0.2: a_true = g (1 + s eps).
            ratio_ta = np.maximum(1.0 + template_scatter * rng.standard_normal(npix), 0.2)
            gs = np.asarray(a_true, dtype=float) / ratio_ta
            floor = 1e-3 * float(np.mean(np.abs(gs)))
            gf = np.maximum(gs, floor)
            templates[name] = gf
            ratio = gf[hi] / a_true[hi]
            report.append(f"template {name} (truth with iid scatter {template_scatter}): ratio floor "
                          f"0.2 hit by {np.mean(ratio_ta <= 0.2) * 100:.2f} percent of pixels, "
                          f"positivity floor {floor:.3g} hits {np.mean(gs < floor) * 100:.2f} percent; "
                          f"g/true at |b| > 20 16/50/84 pct {np.percentile(ratio, 16):.3f}/"
                          f"{np.percentile(ratio, 50):.3f}/{np.percentile(ratio, 84):.3f}")
            continue
        gs = hp.smoothing(g, fwhm=np.radians(TEMPLATE_FWHM_DEG), lmax=lmax, iter=3)
        floor = float(np.percentile(gs, TEMPLATE_FLOOR_PCT))
        floor = max(floor, 1e-3 * float(np.mean(np.abs(gs))))
        gf = np.maximum(gs, floor)
        templates[name] = gf
        ratio = gf[hi] / a_true[hi]
        report.append(f"template {name}: floor {floor:.3g} hits {np.mean(gs < floor) * 100:.1f} "
                      f"percent of pixels ({np.mean(gs[hi] < floor) * 100:.1f} at |b| > 20); "
                      f"g/true at |b| > 20 16/50/84 pct {np.percentile(ratio, 16):.3f}/"
                      f"{np.percentile(ratio, 50):.3f}/{np.percentile(ratio, 84):.3f}")
    hp.write_map(f"{outdir}/template_sync_30GHz.fits", templates["sync"], overwrite=True,
                 dtype=np.float64)
    hp.write_map(f"{outdir}/template_dust_545GHz.fits", templates["dust"], overwrite=True,
                 dtype=np.float64)
    if sky_mode == "consistent":
        beta_eff = beta_d.copy()        # the dust SED is one MBB at T_fixed: no effective index
        beta_eff3 = beta_d.copy()
    else:
        beta_eff = beta_eff_fixed_T(dust_bands, nu_arr, rms, T_fixed)
        beta_eff3 = beta_eff_fixed_T_3comp(sky, nu_arr, rms, beta_s, T_fixed, cmb_sed(nu_arr))
    # Fisher sd of the region indices: amplitudes known (lower bound) and amplitudes profiled
    # (the 5-parameter per-pixel Fisher reduced to the two indices by Schur complement).
    sed_s = (nu_arr[None, :] / 30.0) ** beta_s[:, None]
    sed_d = _dust_sed(nu_arr, beta_d, T_fixed)
    dsed_s = sed_s * np.log(nu_arr / 30.0)[None, :]
    dsed_d = sed_d * np.log(nu_arr / 545.0)[None, :]
    M_true = np.stack([np.broadcast_to(cmb_sed(nu_arr), sed_s.shape), sed_s, sed_d], axis=2)
    dM_true = np.stack([a_sync[:, None] * dsed_s, a_dust[:, None] * dsed_d], axis=2)
    fisher = {}
    for nb in (4, 8, 16):
        region = hp.ud_grade(np.arange(12 * nb * nb).astype(float), nside, order_in="RING",
                             order_out="RING").astype(int)
        prof = fisher_sd_profiled(M_true, dM_true, rms, region, 12 * nb * nb)
        for name, a_true, dsed, sdp in (("sync", a_sync, dsed_s, prof[0]),
                                        ("dust", a_dust, dsed_d, prof[1])):
            sd = fisher_sd(a_true, dsed, rms, region, 12 * nb * nb)
            fisher[f"{name}_{nb}"] = sd
            fisher[f"{name}_{nb}_profiled"] = sdp
            tsd = (beta_s if name == "sync" else beta_d).std()
            report.append(f"Fisher sd {name} nside_beta {nb}: amplitudes known median "
                          f"{np.median(sd):.4f} (16/84 pct {np.percentile(sd, 16):.4f}/"
                          f"{np.percentile(sd, 84):.4f}); amplitudes profiled median "
                          f"{np.median(sdp):.4f} (16/84 pct {np.percentile(sdp, 16):.4f}/"
                          f"{np.percentile(sdp, 84):.4f}), fraction above the truth sd "
                          f"{tsd:.3f}: {np.mean(sdp > tsd) * 100:.0f} percent")
    for name, a_true, dsed in (("sync", a_sync, dsed_s), ("dust", a_dust, dsed_d)):
        report.append(f"Fisher sd {name} full sky: "
                      f"{fisher_sd(a_true, dsed, rms, np.zeros(npix, int), 1)[0]:.5f}")
    report.append(f"beta_eff (fixed T {T_fixed:.3f} K, dust only) minus beta_d: rms "
                  f"{np.std(beta_eff - beta_d):.4f}, mean {np.mean(beta_eff - beta_d):+.4f}, range "
                  f"{np.min(beta_eff - beta_d):+.3f} to {np.max(beta_eff - beta_d):+.3f}; "
                  f"three-component beta_eff3 minus beta_d: rms {np.std(beta_eff3 - beta_d):.4f}, "
                  f"beta_eff3 minus beta_eff rms {np.std(beta_eff3 - beta_eff):.4f}")
    for name, m in (("beta_s", beta_s), ("beta_d", beta_d), ("beta_eff", beta_eff),
                    ("beta_eff3", beta_eff3), ("T_d", T_d)):
        report.append(f"truth {name}: mean {m.mean():.4f} sd {m.std():.4f} (|b| > 20: {m[hi].mean():.4f}"
                      f" / {m[hi].std():.4f}), range {m.min():.3f} to {m.max():.3f}")
    np.savez(f"{outdir}/truth.npz", beta_s=beta_s, beta_d=beta_d, T_d=T_d, T_fixed=T_fixed,
             beta_eff=beta_eff, beta_eff3=beta_eff3, cmb=cmb_map, sky=sky, dust_bands=dust_bands, rms=rms, nu=nu_arr,
             nside=nside, noise_scale=noise_scale, seed=seed, g_sync=templates["sync"],
             g_dust=templates["dust"], a_sync=a_sync, a_dust=a_dust, index_prior_sync=INDEX_PRIOR["sync"],
             index_prior_dust=INDEX_PRIOR["dust"], indep_sd=INDEP_SD, cmb_sd=CMB_SD,
             sky_mode=sky_mode, nside_truth=nside_truth if sky_mode == "consistent" else 0,
             template_scatter=template_scatter if sky_mode == "consistent" else 0.0,
             **{f"fisher_sd_{k}": v for k, v in fisher.items()})
    with open(f"{outdir}/sim_report.txt", "w") as f:
        f.write("\n".join(report) + "\n")
    print("\n".join(report))
    for tag, (kind, nside_beta) in RUNS.items():
        with open(f"{outdir}/param_{tag}.yml", "w") as f:
            f.write(param_file(outdir, tag, kind, nside_beta, nside, T_fixed, n_iter, n_warmup,
                               adapt_mass))
    print(f"{sky_mode} sky written to {outdir}: nside {nside}, noise scale {noise_scale} (rms per pixel "
          f"{np.round(rms, 2).tolist()} uK_RJ), T_fixed {T_fixed:.3f} K, cmb rms "
          f"{cmb_map.std():.1f} uK_CMB; {len(RUNS)} parameter files.")


def _field_block(kind: str | None, nside_beta: int) -> str:
    if kind is None:
        return ""
    if kind == "independent":
        return f"""
      spectral_index_nside: {nside_beta}
      spectral_index_field:
        type: independent"""
    if kind == "iid":
        return f"""
      spectral_index_nside: {nside_beta}
      spectral_index_field:
        type: iid
        sigma_scale: 0.3"""
    return f"""
      spectral_index_nside: {nside_beta}
      spectral_index_field:
        type: gmrf
        basis: graph
        alpha: 2.0
        sigma_scale: 0.3
        rho_prior: loguniform"""


def param_file(outdir: str, tag: str, kind: str | None, nside_beta: int, nside: int,
               T_fixed: float, n_iter: int, n_warmup: int, adapt_mass: bool = True) -> str:
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
    prior_sd = INDEP_SD if kind == "independent" else INDEX_PRIOR["sync"][1]
    return f"""gibbs:
  num_iterations: {n_iter}
  seed: 7
resources:
  tod:
    num_threads: 1
  compsep:
    num_threads: 4
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
        parameterization: centred
        init_inv_mass: fd_hessian
        adapt_mass: {str(adapt_mass).lower()}
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
        rms: {prior_sd}{_field_block(kind, nside_beta)}
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
        rms: {prior_sd}{_field_block(kind, nside_beta)}
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
    ap.add_argument("--sky", choices=["pysm", "consistent"], default="pysm",
                    help="pysm: the PySM sky as is (E020); consistent: a sky from the model class "
                         "(E021), see simulate()")
    ap.add_argument("--nside-truth", type=int, default=16,
                    help="consistent sky: nside of the region-constant true indices")
    ap.add_argument("--template-scatter", type=float, default=0.3,
                    help="consistent sky: iid relative scatter of the true amplitudes about the templates")
    ap.add_argument("--no-adapt-mass", action="store_true",
                    help="write nuts.adapt_mass: false (keep the Hessian metric for the whole run)")
    a = ap.parse_args(argv)
    simulate(a.dir, a.nside, a.noise_scale, a.seed, a.n_iter, a.n_warmup, sky_mode=a.sky,
             nside_truth=a.nside_truth, template_scatter=a.template_scatter,
             adapt_mass=not a.no_adapt_mass)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
