"""E020 analysis: index-map error, region-level calibration, CMB recovery, fit and cost of every
run in a simulation directory written by pysm_intensity_sim.py (metrics as defined in
nextgen/experiments/E020_index_prior_comparison_pysm_d1s1.md before the runs).

    python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e020_d1s1 [--burn 600]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from datetime import datetime

import h5py
import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "aux"))
from validate_mcmc_likelihood import _ess  # noqa: E402  Geyer ESS

GROUP = "sample_index_fields"
COMPS = {"sync": "Synchrotron", "dust": "ThermalDust"}
H_OVER_K = 0.0479924466
T_CMB = 2.7255
N_CMB_DRAWS = 60      # thinned index draws for the GLS CMB estimate
SNR_BINS = 4


def cmb_sed(nu: np.ndarray) -> np.ndarray:
    """uK_CMB to uK_RJ conversion factor per band (the CMB mixing column)."""
    x = H_OVER_K * np.asarray(nu, dtype=float) / T_CMB
    return x ** 2 * np.exp(x) / np.expm1(x) ** 2


def gls_cmb(d: np.ndarray, rms: np.ndarray, nu: np.ndarray, beta_s: np.ndarray,
            beta_d: np.ndarray, T: float) -> np.ndarray:
    """Flat-prior generalized-least-squares CMB amplitude per pixel at the given indices:
    the cmb row of (M^T N^-1 M)^-1 M^T N^-1 d, the standard parametric-fit CMB map."""
    M = np.stack([np.broadcast_to(cmb_sed(nu), d.shape), sed("sync", nu, beta_s, T),
                  sed("dust", nu, beta_d, T)], axis=2)              # (npix, B, 3)
    w = 1.0 / rms ** 2
    A = np.einsum("pbi,b,pbj->pij", M, w, M)
    b = np.einsum("pbi,b,pb->pi", M, w, d)
    return np.linalg.solve(A, b[..., None])[:, 0, 0]


def sed(comp: str, nu: np.ndarray, beta: np.ndarray, T: float) -> np.ndarray:
    """(npix, nband) mixing column of the component at its pivot (30 or 545 GHz)."""
    if comp == "sync":
        return (nu[None, :] / 30.0) ** beta[:, None]
    x, x0 = H_OVER_K * nu / T, H_OVER_K * 545.0 / T
    return (nu[None, :] / 545.0) ** (beta[:, None] + 1.0) * (np.expm1(x0) / np.expm1(x))[None, :]


def region_map(nside_beta: int, nside: int) -> np.ndarray:
    if nside_beta == 0:
        return np.zeros(hp.nside2npix(nside), dtype=int)
    return hp.ud_grade(np.arange(12 * nside_beta**2).astype(float), nside, order_in="RING",
                       order_out="RING").astype(int)


def read_run(outdir: str, tag: str, burn: int) -> dict | None:
    files = sorted(glob.glob(f"{outdir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))
    if len(files) <= burn:
        return None
    out = {"beta": {c: [] for c in COMPS}, "hyper": {c: [] for c in COMPS},
           "amp": {c: [] for c in COMPS}, "cmb_alms": [], "chi2": [], "n_leapfrog": [],
           "divergent": [], "times": [], "kind": None, "nside_beta": 0}
    for fn in files[burn:]:
        with h5py.File(fn, "r") as f:
            st = f[f"mcmc/{GROUP}/state"]
            for c, cls in COMPS.items():
                g = st[cls]
                if "index_field" in g:
                    fld = g["index_field"]
                    out["beta"][c].append(np.asarray(fld["beta_map"][()]))
                    out["nside_beta"] = int(fld["beta_nside"][()])
                    kind = fld["parameterization"][()]
                    kind = kind.decode() if isinstance(kind, bytes) else str(kind)
                    out["kind"] = "independent" if kind == "independent" else "gmrf"
                    if "mu" in fld:
                        out["hyper"][c].append([float(fld["mu"][()]), float(fld["log_sigma"][()]),
                                                float(fld["u"][()])])
                else:
                    out["beta"][c].append(np.atleast_1d(np.asarray(g["beta"][()])))
                    out["kind"] = "scalar"
                out["amp"][c].append([float(g["amp_prior/mu"][()]), float(g["amp_prior/log_s"][()])])
            out["cmb_alms"].append(np.asarray(f["comps/cmb/alms"][()])[0])
            out["chi2"].append(float(f["chi2/reduced"][()]))
            nuts = f[f"mcmc/{GROUP}/nuts"]
            out["n_leapfrog"].append(float(nuts["n_leapfrog"][()]))
            out["divergent"].append(int(nuts["divergent"][()]))
            out["times"].append(f["metadata/datetime"][()])
    for c in COMPS:
        out["beta"][c] = np.array(out["beta"][c])
        out["hyper"][c] = np.array(out["hyper"][c]) if out["hyper"][c] else None
        out["amp"][c] = np.array(out["amp"][c])
    out["n_files"] = len(files)
    return out


def seconds_per_iteration(times) -> float:
    def parse(t):
        t = t.decode() if isinstance(t, bytes) else str(t)
        for fmt in ("%Y-%m-%d %H:%M:%S.%f", "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S.%f",
                    "%Y-%m-%dT%H:%M:%S"):
            try:
                return datetime.strptime(t, fmt)
            except ValueError:
                pass
        return None
    ts = [parse(t) for t in times]
    ts = [t for t in ts if t is not None]
    if len(ts) < 2:
        return float("nan")
    return (ts[-1] - ts[0]).total_seconds() / (len(ts) - 1)


def analyze_run(run: dict, truth: dict, data: np.ndarray, nside: int, hi_mask: np.ndarray,
                snr_bin_edges: dict) -> dict:
    nu = truth["nu"]
    rms = truth["rms"]
    T_fixed = float(truth["T_fixed"])
    nside_beta = run["nside_beta"]
    region = region_map(nside_beta, nside)
    nreg = region.max() + 1
    res = {"kind": run["kind"], "nside_beta": nside_beta, "n_samples": len(run["chi2"]),
           "chi2_reduced_mean": float(np.mean(run["chi2"])),
           "divergences": int(np.sum(run["divergent"])),
           "gradients_per_iteration": float(np.mean(run["n_leapfrog"])),
           "seconds_per_iteration": seconds_per_iteration(run["times"])}
    truths = {"sync": {"pysm": truth["beta_s"]},
              "dust": {"pysm": truth["beta_d"], "fixedT": truth["beta_eff"]}}
    for c in COMPS:
        a_true = truth["a_sync"] if c == "sync" else truth["a_dust"]
        bt = truths[c]["pysm"]
        chain = run["beta"][c]                       # (n, nreg)
        mean_r, sd_r = chain.mean(axis=0), chain.std(axis=0)
        beta_hat = mean_r[region]
        snr_w = np.sum((a_true[:, None] * sed(c, nu, bt, T_fixed) / rms[None, :]) ** 2, axis=1)
        counts = np.bincount(region, minlength=nreg)
        tm_unw = np.bincount(region, bt, minlength=nreg) / counts
        tm_w = np.bincount(region, bt * snr_w, minlength=nreg) / np.bincount(region, snr_w,
                                                                               minlength=nreg)
        z_unw = (mean_r - tm_unw) / sd_r
        z_w = (mean_r - tm_w) / sd_r
        ess = np.array([_ess(chain[:, r]) for r in range(nreg)])
        # RMSE of the posterior mean, and the posterior distribution of the RMSE of a draw.
        rmse_draws = np.sqrt(np.mean((chain[:, region] - bt[None, :]) ** 2, axis=1))
        r = {"post_sd_region_median": float(np.median(sd_r)),
             "z_unweighted_median": float(np.median(z_unw)),
             "z_unweighted_mad": float(np.median(np.abs(z_unw - np.median(z_unw)))),
             "z_weighted_median": float(np.median(z_w)),
             "z_weighted_mad": float(np.median(np.abs(z_w - np.median(z_w)))),
             "frac_abs_z_unw_lt1": float(np.mean(np.abs(z_unw) < 1)),
             "frac_abs_z_unw_lt2": float(np.mean(np.abs(z_unw) < 2)),
             "frac_abs_z_w_lt1": float(np.mean(np.abs(z_w) < 1)),
             "frac_abs_z_w_lt2": float(np.mean(np.abs(z_w) < 2)),
             "ess_min": float(ess.min()), "ess_median": float(np.median(ess)),
             "amp_mu_mean": float(run["amp"][c][:, 0].mean()),
             "amp_log_s_mean": float(run["amp"][c][:, 1].mean()),
             "rmse_draw_mean": float(rmse_draws.mean()), "rmse_draw_sd": float(rmse_draws.std()),
             "region_floor_unweighted": float(np.sqrt(np.mean((tm_unw[region] - bt) ** 2))),
             "region_floor_weighted": float(np.sqrt(np.mean((tm_w[region] - bt) ** 2)))}
        for tname, tmap in truths[c].items():
            err = beta_hat - tmap
            r[f"rmse_{tname}_all"] = float(np.sqrt(np.mean(err ** 2)))
            r[f"rmse_{tname}_hi"] = float(np.sqrt(np.mean(err[hi_mask] ** 2)))
            r[f"bias_{tname}_all"] = float(np.mean(err))
            r[f"truth_{tname}_sd_all"] = float(tmap.std())
            r[f"truth_{tname}_sd_hi"] = float(tmap[hi_mask].std())
        edges = snr_bin_edges[c]
        which = np.clip(np.digitize(snr_w, edges) - 1, 0, SNR_BINS - 1)
        r["rmse_pysm_by_snr_bin"] = [float(np.sqrt(np.mean((beta_hat - bt)[which == k] ** 2)))
                                     for k in range(SNR_BINS)]
        r["truth_sd_by_snr_bin"] = [float(bt[which == k].std()) for k in range(SNR_BINS)]
        if run["hyper"][c] is not None:
            h = run["hyper"][c]
            hgrid = np.sqrt(4 * np.pi / (12 * nside_beta ** 2))
            L = np.log(np.pi / hgrid)
            rho = np.exp(np.log(hgrid) + L / (1 + np.exp(-h[:, 2])))
            r["field_mu"] = [float(h[:, 0].mean()), float(h[:, 0].std())]
            r["field_sigma"] = [float(np.exp(h[:, 1]).mean()), float(np.exp(h[:, 1]).std())]
            r["field_rho_deg"] = [float(np.degrees(rho).mean()), float(np.degrees(rho).std())]
            r["ess_hyper_min"] = float(min(_ess(h[:, j]) for j in range(3)))
        res[c] = r
    # CMB, primary: flat-prior GLS estimate at thinned index draws, averaged, against the input.
    n = len(run["chi2"])
    idx = np.unique(np.linspace(0, n - 1, min(N_CMB_DRAWS, n)).astype(int))
    cmb_gls = np.mean([gls_cmb(data, rms, nu, run["beta"]["sync"][k][region],
                               run["beta"]["dust"][k][region], T_fixed) for k in idx], axis=0)
    cmb_in = truth["cmb"]
    res["cmb"] = cmb_metrics(cmb_gls - cmb_in, cmb_in, hi_mask)
    res["cmb"]["estimator"] = f"GLS at {idx.size} thinned index draws, averaged"
    # CMB, secondary: the chain's own (prior-shrunk) amplitude alms.
    alms = np.mean(np.array(run["cmb_alms"]), axis=0)
    lmax = hp.Alm.getlmax(alms.size)
    cmb_hat = hp.alm2map(np.ascontiguousarray(alms, dtype=np.complex128), nside, lmax=lmax)
    cmb_in_bl = hp.alm2map(hp.map2alm(cmb_in, lmax=lmax, iter=3), nside, lmax=lmax)
    res["cmb_chain_alms"] = cmb_metrics(cmb_hat - cmb_in_bl, cmb_in_bl, hi_mask)
    res["cmb_chain_alms"]["note"] = "posterior-mean chain alms (shrunk by the N(0, cmb_sd^2) " \
        "prior); input band-limited to the chain lmax"
    return res


def cmb_metrics(resid: np.ndarray, cmb_in: np.ndarray, hi_mask: np.ndarray) -> dict:
    """RMS of the residual, residual power over input power and the regression coefficient of
    the estimate on the input (1 + cross/in) in l 2 to 30 and 30 to 100."""
    lmax = 100
    cl_res = hp.anafast(resid, lmax=lmax)
    cl_in = hp.anafast(cmb_in, lmax=lmax)
    cl_x = hp.anafast(resid, cmb_in, lmax=lmax)
    out = {"resid_rms_all": float(resid.std()), "resid_rms_hi": float(resid[hi_mask].std())}
    for name, sl in (("l2_30", slice(2, 31)), ("l30_100", slice(30, 101))):
        out[f"resid_power_frac_{name}"] = float(cl_res[sl].sum() / cl_in[sl].sum())
        out[f"regression_coeff_{name}"] = float(1.0 + cl_x[sl].sum() / cl_in[sl].sum())
    return out


def snr_bin_edges(truth: dict) -> dict:
    """Quartile edges of the per-pixel index S/N proxy per component (shared by all runs)."""
    nu, rms, T = truth["nu"], truth["rms"], float(truth["T_fixed"])
    edges = {}
    for c, a, bt in (("sync", truth["a_sync"], truth["beta_s"]),
                     ("dust", truth["a_dust"], truth["beta_d"])):
        w = np.sum((a[:, None] * sed(c, nu, bt, T) / rms[None, :]) ** 2, axis=1)
        edges[c] = np.quantile(w, np.linspace(0, 1, SNR_BINS + 1))
    return edges


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--burn", type=int, default=600)
    ap.add_argument("--tags", nargs="*", default=None)
    a = ap.parse_args(argv)
    truth = dict(np.load(f"{a.dir}/truth.npz"))
    nside = int(truth["nside"])
    nu = truth["nu"]
    data = np.array([hp.read_map(f"{a.dir}/band{int(f):03d}_map.fits") for f in nu]).T
    lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)[1]
    hi = np.abs(lat) > 20
    edges = snr_bin_edges(truth)
    tags = a.tags or sorted(os.path.basename(d)[7:] for d in glob.glob(f"{a.dir}/chains_*")
                            if "timing" not in d)
    # Reference: GLS CMB at the true indices (the noise floor no index model can beat) and at
    # the fixed-T effective dust index.
    ref = {"gls_at_truth": cmb_metrics(gls_cmb(data, truth["rms"], nu, truth["beta_s"],
                                               truth["beta_d"], float(truth["T_fixed"]))
                                       - truth["cmb"], truth["cmb"], hi),
           "gls_at_fixedT_truth": cmb_metrics(gls_cmb(data, truth["rms"], nu, truth["beta_s"],
                                                      truth["beta_eff"], float(truth["T_fixed"]))
                                              - truth["cmb"], truth["cmb"], hi)}
    summary = {"reference": ref, "runs": {}}
    lines = [f"E020 analysis of {a.dir} (burn {a.burn}; nside {nside}; truth sd beta_s "
             f"{truth['beta_s'].std():.4f}, beta_d {truth['beta_d'].std():.4f}, beta_eff "
             f"{truth['beta_eff'].std():.4f}; beta_eff - beta_d rms "
             f"{np.std(truth['beta_eff'] - truth['beta_d']):.4f})",
             f"reference GLS CMB at the true indices: rms {ref['gls_at_truth']['resid_rms_all']:.2f}"
             f" / {ref['gls_at_truth']['resid_rms_hi']:.2f} uK, residual power fraction "
             f"{ref['gls_at_truth']['resid_power_frac_l2_30']:.4f} (l 2-30) "
             f"{ref['gls_at_truth']['resid_power_frac_l30_100']:.4f} (l 30-100); at beta_eff: "
             f"rms {ref['gls_at_fixedT_truth']['resid_rms_all']:.2f}, fractions "
             f"{ref['gls_at_fixedT_truth']['resid_power_frac_l2_30']:.4f} / "
             f"{ref['gls_at_fixedT_truth']['resid_power_frac_l30_100']:.4f}", "",
             f"{'run':<8s} {'n':>5s} {'rmse_s all/hi':>15s} {'rmse_d all/hi':>15s} "
             f"{'rmse_d(fixT)':>12s} {'floor_s/d':>13s} {'rmse draw s/d':>19s} "
             f"{'sd_s/sd_d':>13s} {'z_s med/mad':>12s} {'z_d med/mad':>12s} "
             f"{'cmb rms all/hi':>15s} {'P<30 P<100':>13s} {'coef<30':>7s} {'chi2':>6s} "
             f"{'grad':>5s} {'s/it':>5s} {'ESS_s/d':>10s} {'div':>4s}"]
    for tag in tags:
        run = read_run(a.dir, tag, a.burn)
        if run is None:
            lines.append(f"{tag:<8s} (not enough samples)")
            continue
        r = analyze_run(run, truth, data, nside, hi, edges)
        summary["runs"][tag] = r
        s, d, c = r["sync"], r["dust"], r["cmb"]
        lines.append(
            f"{tag:<8s} {r['n_samples']:>5d} {s['rmse_pysm_all']:.4f}/{s['rmse_pysm_hi']:.4f} "
            f"{d['rmse_pysm_all']:.4f}/{d['rmse_pysm_hi']:.4f} {d['rmse_fixedT_all']:.4f}/"
            f"{d['rmse_fixedT_hi']:.4f} {s['region_floor_unweighted']:.4f}/"
            f"{d['region_floor_unweighted']:.4f} {s['rmse_draw_mean']:.4f}+-{s['rmse_draw_sd']:.4f}"
            f"/{d['rmse_draw_mean']:.4f}+-{d['rmse_draw_sd']:.4f} {s['post_sd_region_median']:.4f}/"
            f"{d['post_sd_region_median']:.4f} {s['z_weighted_median']:+5.2f}/"
            f"{s['z_weighted_mad']:4.2f} {d['z_weighted_median']:+5.2f}/{d['z_weighted_mad']:4.2f} "
            f"{c['resid_rms_all']:7.2f}/{c['resid_rms_hi']:6.2f} {c['resid_power_frac_l2_30']:.4f} "
            f"{c['resid_power_frac_l30_100']:.4f} {c['regression_coeff_l2_30']:7.3f} "
            f"{r['chi2_reduced_mean']:6.3f} {r['gradients_per_iteration']:5.1f} "
            f"{r['seconds_per_iteration']:5.1f} {s['ess_min']:4.0f}/{d['ess_min']:4.0f} "
            f"{r['divergences']:4d}")
        lines.append(f"    rmse by S/N quartile (low to high): sync "
                     f"{' '.join(f'{v:.4f}' for v in s['rmse_pysm_by_snr_bin'])} (truth sd "
                     f"{' '.join(f'{v:.4f}' for v in s['truth_sd_by_snr_bin'])}); dust "
                     f"{' '.join(f'{v:.4f}' for v in d['rmse_pysm_by_snr_bin'])} (truth sd "
                     f"{' '.join(f'{v:.4f}' for v in d['truth_sd_by_snr_bin'])}); |z_w|<1,<2: sync "
                     f"{s['frac_abs_z_w_lt1']:.2f}/{s['frac_abs_z_w_lt2']:.2f}, dust "
                     f"{d['frac_abs_z_w_lt1']:.2f}/{d['frac_abs_z_w_lt2']:.2f}")
        for cname in ("sync", "dust"):
            if "field_mu" in r[cname]:
                f = r[cname]
                lines.append(f"    {cname} field: mu {f['field_mu'][0]:.4f} +- {f['field_mu'][1]:.4f},"
                             f" sigma {f['field_sigma'][0]:.4f} +- {f['field_sigma'][1]:.4f}, rho "
                             f"{f['field_rho_deg'][0]:.1f} +- {f['field_rho_deg'][1]:.1f} deg; amp "
                             f"mu {f['amp_mu_mean']:.3f} log s {f['amp_log_s_mean']:.2f}; ESS hyper"
                             f" min {f['ess_hyper_min']:.0f}")
    lines += ["", "columns: rmse = posterior-mean index at the data pixels minus the PySM truth, "
              "all sky / |b| > 20; rmse_d(fixT) = dust against the fixed-T effective index "
              "beta_eff; floor = RMSE of the truth's own region means (unbeatable by a region "
              "model at that nside); rmse draw = mean +- sd over posterior draws of the RMSE of a "
              "draw (the uncertainty scale for comparing runs); sd = median posterior sd of a "
              "region index; z = (region mean - S/N-weighted truth region mean) / sd, median and "
              "MAD over regions; cmb = flat-prior GLS CMB at thinned index draws minus the input, "
              "rms in uK_CMB, residual power over input power at l 2-30 and 30-100, regression "
              "coefficient of estimate on input at l 2-30; grad = gradients per Gibbs iteration; "
              "s/it = wall seconds per iteration (not comparable across load); ESS = min Geyer ESS "
              "over regions."]
    text = "\n".join(lines)
    print(text)
    with open(f"{a.dir}/e020_summary.txt", "w") as f:
        f.write(text + "\n")
    with open(f"{a.dir}/e020_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
