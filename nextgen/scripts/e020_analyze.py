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


def analyze_run(run: dict, truth: dict, nside: int, hi_mask: np.ndarray) -> dict:
    nu = truth["nu"]
    rms = truth["rms"]
    nside_beta = run["nside_beta"]
    region = region_map(nside_beta, nside)
    nreg = region.max() + 1
    res = {"kind": run["kind"], "nside_beta": nside_beta, "n_samples": len(run["chi2"]),
           "chi2_reduced_mean": float(np.mean(run["chi2"])),
           "divergences": int(np.sum(run["divergent"])),
           "gradients_per_iteration": float(np.mean(run["n_leapfrog"])),
           "seconds_per_iteration": seconds_per_iteration(run["times"])}
    for c in COMPS:
        bt = truth["beta_s"] if c == "sync" else truth["beta_d"]
        a_true = truth["a_sync"] if c == "sync" else truth["a_dust"]
        chain = run["beta"][c]                       # (n, nreg)
        mean_r, sd_r = chain.mean(axis=0), chain.std(axis=0)
        beta_hat = mean_r[region]
        err = beta_hat - bt
        snr_w = np.sum((a_true[:, None] * sed(c, nu, bt, float(truth["T_fixed"])) / rms[None, :])
                       ** 2, axis=1)
        tm_unw = np.bincount(region, bt, minlength=nreg) / np.bincount(region, minlength=nreg)
        tm_w = np.bincount(region, bt * snr_w, minlength=nreg) / np.bincount(region, snr_w,
                                                                               minlength=nreg)
        z_unw = (mean_r - tm_unw) / sd_r
        z_w = (mean_r - tm_w) / sd_r
        ess = np.array([_ess(chain[:, r]) for r in range(nreg)])
        r = {"rmse_all": float(np.sqrt(np.mean(err**2))),
             "rmse_hi": float(np.sqrt(np.mean(err[hi_mask]**2))),
             "bias_all": float(np.mean(err)),
             "truth_sd_all": float(bt.std()), "truth_sd_hi": float(bt[hi_mask].std()),
             "region_floor_rmse": float(np.sqrt(np.mean((tm_unw[region] - bt)**2))),
             "post_sd_region_median": float(np.median(sd_r)),
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
             "amp_log_s_mean": float(run["amp"][c][:, 1].mean())}
        if run["hyper"][c] is not None:
            h = run["hyper"][c]
            hgrid = np.sqrt(4 * np.pi / (12 * nside_beta**2))
            L = np.log(np.pi / hgrid)
            rho = np.exp(np.log(hgrid) + L / (1 + np.exp(-h[:, 2])))
            r["field_mu"] = [float(h[:, 0].mean()), float(h[:, 0].std())]
            r["field_sigma"] = [float(np.exp(h[:, 1]).mean()), float(np.exp(h[:, 1]).std())]
            r["field_rho_deg"] = [float(np.degrees(rho).mean()), float(np.degrees(rho).std())]
            r["ess_hyper_min"] = float(min(_ess(h[:, j]) for j in range(3)))
        res[c] = r
    # CMB: posterior-mean alms to a map at the data nside against the input c1 map.
    alms = np.mean(np.array(run["cmb_alms"]), axis=0)
    lmax = hp.Alm.getlmax(alms.size)
    cmb_hat = hp.alm2map(np.ascontiguousarray(alms, dtype=np.complex128), nside, lmax=lmax)
    cmb_in = truth["cmb"]
    cmb_in_bl = hp.alm2map(hp.map2alm(cmb_in, lmax=lmax, iter=3), nside, lmax=lmax)
    resid = cmb_hat - cmb_in_bl
    cl_res = hp.anafast(resid, lmax=min(100, lmax))
    cl_in = hp.anafast(cmb_in_bl, lmax=min(100, lmax))
    res["cmb"] = {"resid_rms_all": float(resid.std()), "resid_rms_hi": float(resid[hi_mask].std()),
                  "resid_power_frac_l2_30": float(cl_res[2:31].sum() / cl_in[2:31].sum()),
                  "resid_power_frac_l30_100": float(cl_res[30:101].sum() / cl_in[30:101].sum()),
                  "note": "input CMB band-limited to the chain lmax before differencing"}
    return res


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--burn", type=int, default=600)
    ap.add_argument("--tags", nargs="*", default=None)
    a = ap.parse_args(argv)
    truth = dict(np.load(f"{a.dir}/truth.npz"))
    nside = int(truth["nside"])
    lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)[1]
    hi = np.abs(lat) > 20
    tags = a.tags or sorted(os.path.basename(d)[7:] for d in glob.glob(f"{a.dir}/chains_*")
                            if "timing" not in d)
    summary = {}
    lines = [f"E020 analysis of {a.dir} (burn {a.burn}; nside {nside}; truth sd beta_s "
             f"{truth['beta_s'].std():.4f}, beta_d {truth['beta_d'].std():.4f})", "",
             f"{'run':<9s} {'n':>5s} {'rmse_s all/hi':>15s} {'rmse_d all/hi':>15s} "
             f"{'floor_s/d':>13s} {'sd_s/sd_d':>13s} {'z_s med/mad':>12s} {'z_d med/mad':>12s} "
             f"{'cmb rms all/hi':>15s} {'P<30 P<100':>12s} {'chi2':>6s} {'grad':>5s} {'s/it':>5s} "
             f"{'ESS_s/d':>10s} {'div':>4s}"]
    for tag in tags:
        run = read_run(a.dir, tag, a.burn)
        if run is None:
            lines.append(f"{tag:<9s} (not enough samples)")
            continue
        r = analyze_run(run, truth, nside, hi)
        summary[tag] = r
        s, d, c = r["sync"], r["dust"], r["cmb"]
        lines.append(
            f"{tag:<9s} {r['n_samples']:>5d} {s['rmse_all']:.4f}/{s['rmse_hi']:.4f} "
            f"{d['rmse_all']:.4f}/{d['rmse_hi']:.4f} {s['region_floor_rmse']:.4f}/"
            f"{d['region_floor_rmse']:.4f} {s['post_sd_region_median']:.4f}/"
            f"{d['post_sd_region_median']:.4f} {s['z_weighted_median']:+5.2f}/"
            f"{s['z_weighted_mad']:4.2f} {d['z_weighted_median']:+5.2f}/{d['z_weighted_mad']:4.2f} "
            f"{c['resid_rms_all']:7.2f}/{c['resid_rms_hi']:6.2f} {c['resid_power_frac_l2_30']:.4f} "
            f"{c['resid_power_frac_l30_100']:.4f} {r['chi2_reduced_mean']:6.3f} "
            f"{r['gradients_per_iteration']:5.1f} {r['seconds_per_iteration']:5.1f} "
            f"{s['ess_min']:4.0f}/{d['ess_min']:4.0f} {r['divergences']:4d}")
        for cname in ("sync", "dust"):
            if "field_mu" in r[cname]:
                f = r[cname]
                lines.append(f"    {cname} field: mu {f['field_mu'][0]:.4f} +- {f['field_mu'][1]:.4f}, "
                             f"sigma {f['field_sigma'][0]:.4f} +- {f['field_sigma'][1]:.4f}, rho "
                             f"{f['field_rho_deg'][0]:.1f} +- {f['field_rho_deg'][1]:.1f} deg; amp "
                             f"mu {f['amp_mu_mean']:.3f} log s {f['amp_log_s_mean']:.2f}; ESS hyper "
                             f"min {f['ess_hyper_min']:.0f}")
    lines += ["", "columns: rmse = posterior-mean index at the data pixels minus truth, all sky / "
              "|b| > 20; floor = RMSE of the truth's own region means (unbeatable by a region "
              "model at that nside); sd = median posterior sd of a region index; z = (region mean "
              "- S/N-weighted truth region mean) / sd, median and MAD over regions; cmb rms = "
              "posterior-mean CMB map minus input (band-limited), uK_CMB; P<30, P<100 = residual "
              "power over input CMB power at l 2 to 30 and 30 to 100; grad = gradients per Gibbs "
              "iteration; s/it = wall seconds per iteration (not comparable across load); ESS = "
              "min Geyer ESS over regions."]
    text = "\n".join(lines)
    print(text)
    with open(f"{a.dir}/e020_summary.txt", "w") as f:
        f.write(text + "\n")
    with open(f"{a.dir}/e020_summary.json", "w") as f:
        json.dump(summary, f, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
