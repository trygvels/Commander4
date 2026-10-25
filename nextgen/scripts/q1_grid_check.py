"""E014 reference: exact ridge and flat-prior-marginal posteriors of (beta_s, beta_d) on a grid.

Reads the bands of a Commander4 Q1 parameter file through Commander4's own map reader (same
degrade, unit conversion and noise-weighted smoothing as the chain saw), builds the per-pixel
mixing matrix from the same component classes, restricts to the group's chisq mask, and evaluates
on a grid

    log p_ridge(beta)    = -chi2_min(beta)/2 + log prior(beta),
    log p_marginal(beta) = log p_ridge(beta) - 1/2 sum_p ln det A_p(beta),

with A_p = M^T N_p^-1 M and chi2_min the per-pixel generalized-least-squares minimum (flat amplitude
priors). The marginal is also the beta marginal of the `conditional` chain's joint target. Reports
grid means and sds (coarse grid over the prior box, then a refined grid around the mode) and writes
them to <chain_dir>/grid_check.json. Independent of the MH code in `compsep/mcmc.py`; shares the
data path and SEDs with the chain by construction (that is what is being checked against).

    python -I nextgen/scripts/q1_grid_check.py nextgen/params/q1/param_q1_n64_ridge_npipemask.yml
"""
import json
import os
import sys

import healpy as hp
import numpy as np

from commander4.compsep.mcmc import resolve_chisq_mask
from commander4.file_io.map_reader import read_data_map_from_file
from commander4.parameters.parse import load_params
from commander4.sky.comp_list import CompList


def load(param_file):
    params, _, _ = load_params(param_file)
    bands = params.compsep.bands
    d, ivar, freqs, names = [], [], [], []
    for name in bands:
        b = bands[name]
        if not b.enabled:
            continue
        if "_name" not in b:
            b._name = name
        dm = read_data_map_from_file(b, params)
        d.append(dm.map_sky[0]); ivar.append(dm.map_inv_var[0]); freqs.append(float(b.freq))
        names.append(name)
    d, ivar, freqs = np.array(d), np.array(ivar), np.array(freqs)
    nside = hp.npix2nside(d.shape[1])
    comps = CompList.init_from_params(params.components, params)
    comps = list(comps.split_for_eval_pol("I"))
    group = params.compsep.mcmc_sampling_groups.sample_spectral_indices
    keep = None
    if "chisq_mask" in group and group.chisq_mask:
        keep = resolve_chisq_mask(hp.read_map(group.chisq_mask, field=0, dtype=np.float64), nside)
    if keep is None:
        keep = np.ones(d.shape[1], bool)
    return params, d[:, keep], ivar[:, keep], freqs, names, comps


def find(comps, cls):
    out = [c for c in comps if type(c).__name__ == cls]
    assert out, cls
    return out


def loglikes(d, ivar, M):
    """Per-pixel chi2_min and ln det A summed over pixels for mixing matrix M (nband, ncomp)."""
    # A_p = M^T diag(ivar_p) M ; b_p = M^T (ivar_p * d_p)
    A = np.einsum("bi,bp,bj->pij", M, ivar, M)
    b = np.einsum("bi,bp->pi", M, ivar * d)
    x = np.linalg.solve(A, b[..., None])[..., 0]
    chi2 = np.sum(ivar * d * d) - np.sum(b * x)
    sign, logdet = np.linalg.slogdet(A)
    assert np.all(sign > 0)
    return chi2, float(np.sum(logdet))


def main():
    pf = sys.argv[1]
    params, d, ivar, freqs, names, comps = load(pf)
    sync = find(comps, "Synchrotron")
    dust = find(comps, "ThermalDust")
    dp = params.components.ThermalDust.params
    sp = params.components.Synchrotron.params
    lo = np.array([sp.spectral_index_bounds[0], dp.spectral_index_bounds[0]])
    hi = np.array([sp.spectral_index_bounds[1], dp.spectral_index_bounds[1]])
    pm, pr = float(dp.spectral_index_prior.mean), float(dp.spectral_index_prior.rms)
    print(f"{d.shape[0]} bands, {d.shape[1]} pixels kept, comps "
          f"{[type(c).__name__ for c in comps]}")

    def logpost(bs, bd):
        for c in sync: c.beta = bs
        for c in dust: c.beta = bd
        M = np.stack([c.get_sed(freqs) for c in comps], axis=1)
        chi2, logdet = loglikes(d, ivar, M)
        lp = -0.5 * ((bd - pm) / pr) ** 2
        return -0.5 * chi2 + lp, -0.5 * chi2 + lp - 0.5 * logdet

    def moments(lo_, hi_, n):
        gs = np.linspace(lo_[0], hi_[0], n); gd = np.linspace(lo_[1], hi_[1], n)
        L = np.array([[logpost(a, b) for b in gd] for a in gs])  # (ns, nd, 2)
        out = {}
        for k, key in enumerate(["ridge", "marginal"]):
            w = np.exp(L[..., k] - L[..., k].max()); w /= w.sum()
            ms = np.sum(w * gs[:, None]); md = np.sum(w * gd[None, :])
            ss = np.sqrt(np.sum(w * (gs[:, None] - ms) ** 2)); sd = np.sqrt(np.sum(w * (gd[None, :] - md) ** 2))
            out[key] = dict(mean_sync=ms, sd_sync=ss, mean_dust=md, sd_dust=sd,
                            step=[gs[1] - gs[0], gd[1] - gd[0]])
        return out

    res = {"coarse": moments(lo, hi, 41)}
    print("coarse", json.dumps(res["coarse"], indent=1))
    final = {}
    for key in ["ridge", "marginal"]:
        m = res["coarse"][key]
        c = np.array([m["mean_sync"], m["mean_dust"]])
        half = np.maximum(5 * np.array([m["sd_sync"], m["sd_dust"]]), 2 * np.array(m["step"]))
        for it in range(3):
            r = moments(np.maximum(lo, c - half), np.minimum(hi, c + half), 61)[key]
            c = np.array([r["mean_sync"], r["mean_dust"]])
            half = 6 * np.array([r["sd_sync"], r["sd_dust"]])
            print(key, "refine", it, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
        final[key] = r
    res["final"] = final
    out = os.path.join(params.output.dir if "dir" in params.output else ".", "grid_check.json")
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    print("written", out)


if __name__ == "__main__":
    main()
