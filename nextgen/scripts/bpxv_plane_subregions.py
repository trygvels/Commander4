"""E022 follow-up: does the Plane's resolution dependence come from fitting one index to a
region whose index varies? Split the Plane into nside-4 sub-regions and compute the collapsed
target (zero-mean data-scaled priors on all amplitudes, with the N(-3.3, 0.1) prior) per
sub-region at nside 16, 32 and 64; a resolution-stable set of sub-region values with a moving
single-region value supports the hypothesis."""
import os, sys, numpy as np, healpy as hp
sys.path.insert(0, os.path.dirname(__file__))
from bpxv_gibbs_coupling import Surrogate, moments
D = os.path.expanduser("~/data"); grid = np.linspace(-4.0, -2.0, 801); NSUB = 4
cmb = {16: 0.187, 32: 0.321, 64: 0.439}
out = {}
for ns in (16, 32, 64):
    sur = Surrogate(f"{D}/bpxv_surrogate_n{ns}_nowmap", lfi_noise="white", cmb_prior_sd=cmb[ns], synch_prior_sd="rms", dust_prior_sd="rms")
    th, ph = hp.pix2ang(ns, np.arange(sur.npix)); sub = hp.ang2pix(NSUB, th, ph)
    plane = sur.region == 3
    m, sd, _ = moments(grid, sur.collapsed_logp(grid, 3, prior=True)); out[(ns, "all")] = (m, sd, int(plane.sum()))
    labels = np.where(plane, 100 + sub, sur.region); sur.region = labels
    for sp in np.unique(sub[plane]):
        n = int((labels == 100 + sp).sum())
        if n < 4: continue
        m, sd, _ = moments(grid, sur.collapsed_logp(grid, 100 + sp, prior=True)); out[(ns, int(sp))] = (m, sd, n)
    print(f"nside {ns} done", flush=True)
subs = sorted({k[1] for k in out if k[1] != "all"})
print("| Plane part | " + " | ".join(f"nside {ns}: mean +- sd (npix)" for ns in (16, 32, 64)) + " |")
print("|---|---|---|---|")
for sp in ["all"] + subs:
    cells = []
    for ns in (16, 32, 64):
        v = out.get((ns, sp)); cells.append(f"{v[0]:.3f} +- {v[1]:.3f} ({v[2]})" if v else "-")
    print(f"| {'whole Plane' if sp == 'all' else f'nside-4 pixel {sp}'} | " + " | ".join(cells) + " |")
# pixel-weighted mean of the sub-region values at each nside, against the whole-Plane fit
for ns in (16, 32, 64):
    vals = [(out[(ns, sp)][0], out[(ns, sp)][2]) for sp in subs if (ns, sp) in out]
    w = np.array([v[1] for v in vals], float); mu = np.array([v[0] for v in vals])
    print(f"nside {ns}: pixel-weighted mean of sub-region values {np.sum(w*mu)/w.sum():.3f}, spread (weighted sd) {np.sqrt(np.sum(w*(mu-np.sum(w*mu)/w.sum())**2)/w.sum()):.3f}, whole-Plane fit {out[(ns,'all')][0]:.3f}")
