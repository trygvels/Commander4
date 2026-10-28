"""E022 figure: per-region 1-D posteriors and Gibbs traces for one variant (result_<tag>.json)."""
import argparse, json, os, sys
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from bpxv_gibbs_coupling import Surrogate, REGIONS, PRIOR_MEAN, PRIOR_SD, moments
ap = argparse.ArgumentParser(); ap.add_argument("--dir", default=os.path.expanduser("~/data/bpxv_surrogate"))
ap.add_argument("--out", default="nextgen/experiments/E022_bp_xv_polarized_index_gibbs_coupling"); ap.add_argument("--tag", default="bp_cmb0")
a = ap.parse_args()
res = json.load(open(os.path.join(a.out, f"result_{a.tag}.json"))); v = res["variant"]
sur = Surrogate(a.dir, lfi_noise=v["lfi_noise"], cmb_prior_sd=v["cmb_prior_sd"]); grid = np.linspace(-4.5, -1.5, 2001)
rng = np.random.default_rng(5)
beta_star = np.array([PRIOR_MEAN, res["regions"]["Spur"]["marg_prior"]["mean"], PRIOR_MEAN, res["regions"]["Plane"]["marg_prior"]["mean"]])
amp = sur.draw_amp(beta_star, rng)
fig, ax = plt.subplots(2, 4, figsize=(16, 6.5))
for r, name in enumerate(REGIONS):
    for lp, lab, ls in ((sur.collapsed_logp(grid, r, prior=True), "collapsed (T2) + prior", "-"), (sur.collapsed_logp(grid, r, prior=True, logdet=False), "profile + prior", "--"), (sur.cond_logp(grid, amp, r, prior=True), "conditional at one amplitude draw (T1)", ":")):
        m, sd, p = moments(grid, lp); ax[0, r].plot(grid, p / p.max(), ls, label=f"{lab}: {m:.3f}+-{sd:.3f}")
    ax[0, r].plot(grid, np.exp(-0.5 * ((grid - PRIOR_MEAN) / PRIOR_SD) ** 2), color="gray", alpha=0.5, label="prior N(-3.3, 0.1)")
    ax[0, r].set_xlim(-3.6, -2.4); ax[0, r].set_title(f"{name} ({res['regions'][name]['npix']} px)"); ax[0, r].legend(fontsize=6); ax[0, r].set_xlabel("beta_s")
    for ctag, lab, al in (("bp", "Gibbs, BP scheme", 1.0), ("all", "Gibbs, all regions from likelihood", 0.6)):
        fn = os.path.join(a.out, f"gibbs_chain_{ctag}_{a.tag}.npy")
        if os.path.exists(fn): ax[1, r].plot(np.load(fn)[:1500, r], lw=0.5, alpha=al, label=lab)
    g = res["gibbs_bp_scheme"].get(name); ga = res["gibbs_all_regions"].get(name)
    ax[1, r].set_title((f"tau {g['tau']:.1f}+-{g['tau_err']:.1f} (BP scheme), " if g else "") + (f"{ga['tau']:.1f}+-{ga['tau_err']:.1f} (all)" if ga else ""), fontsize=9)
    ax[1, r].set_xlabel("Gibbs iteration"); ax[1, r].legend(fontsize=6)
fig.suptitle(f"E022 variant {a.tag}: LFI noise {v['lfi_noise']}, CMB prior sd {v['cmb_prior_sd']} uK", fontsize=10)
plt.tight_layout(); plt.savefig(os.path.join(a.out, f"posteriors_and_chains_{a.tag}.png"), dpi=90); print("saved", a.tag)
