"""E032 item (viii): posterior-mean template multiplier maps against the realized template error.

The cs template is g = a_true exp(f), so the multiplier of a region should track the template-weighted
region mean of a_true / g. Per S/N quartile: correlation of the posterior mean with that ratio, and
the rms difference over the sky. Same construction as the E030 chain check sets.

    python nextgen/scripts/e032_mu_maps.py --dir D --tags TAG ... --refs REFS.npz [--burn 500]
"""
import argparse
import glob

import h5py
import healpy as hp
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True); ap.add_argument("--tags", nargs="+", required=True)
ap.add_argument("--refs", required=True); ap.add_argument("--burn", type=int, default=500)
ap.add_argument("--sync-template", default="template_sync_cs_30GHz.fits")
ap.add_argument("--dust-template", default="template_dust_cs_545GHz.fits")
a = ap.parse_args()
tr = np.load(f"{a.dir}/truth.npz", allow_pickle=True); r = np.load(a.refs)
q = {"Synchrotron": r["qs16"], "ThermalDust": r["qd16"]}
reg = hp.ud_grade(np.arange(3072.0), 64).astype(int)


def ratio(comp, tmpl):
    g = hp.read_map(f"{a.dir}/{tmpl}", field=0)
    amp = np.asarray(tr[{"ThermalDust": "a_dust", "Synchrotron": "a_sync"}[comp]], float)
    w = np.abs(g)
    return np.bincount(reg, w * amp / np.maximum(g, 1e-3), 3072) / np.maximum(np.bincount(reg, w, 3072), 1e-30)


print("posterior-mean template multiplier maps, post-burn draws; quartiles Q1 bright to Q4 faint")
for tag in a.tags:
    files = sorted(glob.glob(f"{a.dir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
    for comp, tm in (("Synchrotron", a.sync_template), ("ThermalDust", a.dust_template)):
        mu = np.array([h5py.File(f, "r")[f"mcmc/sample_index_fields/state/{comp}/amp_prior/mu"][()] for f in files])
        pm, sd = mu.mean(0), mu.std(0); rr = ratio(comp, tm)
        print(f"{tag:40s} {comp:12s} mean {pm.mean():.4f} (realized {rr.mean():.4f}) | corr per quartile "
              + " ".join(f"{np.corrcoef(pm[q[comp] == k], rr[q[comp] == k])[0, 1]:+.3f}" for k in range(4))
              + f" | all {np.corrcoef(pm, rr)[0, 1]:+.3f} | rms(pm - a/g) {np.sqrt(np.mean((pm - rr) ** 2)):.4f} per quartile "
              + " ".join(f"{np.sqrt(np.mean((pm - rr)[q[comp] == k] ** 2)):.3f}" for k in range(4))
              + " | median post sd " + " ".join(f"{np.median(sd[q[comp] == k]):.3f}" for k in range(4)))
