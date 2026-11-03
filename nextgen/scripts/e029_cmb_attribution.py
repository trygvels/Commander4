#!/usr/bin/env python
"""E029: attribute the GLS CMB residual difference between two chains to the sync and dust
index fields by evaluating the flat-prior GLS CMB at hybrid posterior-mean fields (sync from one
chain, dust from the other), and at the chain's fields with the dust index replaced by a
reference in a subset of regions.
usage: e029_cmb_attribution.py --dir DIR --refs REFS.npz --a TAG_A --b TAG_B [--burn 500]"""
import argparse, glob, sys
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); ap.add_argument("--refs", required=True)
ap.add_argument("--a", required=True); ap.add_argument("--b", required=True); ap.add_argument("--burn", type=int, default=500)
a = ap.parse_args(); sys.argv = [sys.argv[0]]
import numpy as np, healpy as hp, h5py
sys.path.insert(0, "nextgen/scripts")
from e020_analyze import gls_cmb, region_map
truth = dict(np.load(f"{a.dir}/truth.npz")); nside = int(truth["nside"]); nu = truth["nu"]
data = np.array([hp.read_map(f"{a.dir}/band{int(f):03d}_map.fits") for f in nu]).T
lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)[1]; hi = np.abs(lat) > 20
T = float(truth["T_fixed"]); reg = region_map(16, nside); r = np.load(a.refs); qd = r["qd16"]; qs = r["qs16"]

def pm(tag, comp):
    files = sorted(glob.glob(f"{a.dir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
    return np.mean([h5py.File(f)[f"mcmc/sample_index_fields/state/{comp}/index_field/beta_map"][()] for f in files], 0)

def resid(bs_r, bd_r, label):
    m = gls_cmb(data, truth["rms"], nu, bs_r[reg], bd_r[reg], T) - truth["cmb"]
    print(f"{label:60s} rms all {np.sqrt(np.mean(m**2)):7.3f}  |b|>20 {np.sqrt(np.mean(m[hi]**2)):7.3f}")

S = {t: pm(t, "Synchrotron") for t in (a.a, a.b)}; Dd = {t: pm(t, "ThermalDust") for t in (a.a, a.b)}
print(f"GLS CMB residual (uK) at posterior-mean index fields; A = {a.a}, B = {a.b}")
resid(S[a.a], Dd[a.a], "A sync, A dust"); resid(S[a.b], Dd[a.b], "B sync, B dust")
resid(S[a.a], Dd[a.b], "A sync, B dust"); resid(S[a.b], Dd[a.a], "B sync, A dust")
for name in ("bd_fixT", "bd_trueT"):
    resid(S[a.a], r[name], f"A sync, dust = {name} everywhere")
    for k in range(4):
        d = Dd[a.a].copy(); d[qd == k] = r[name][qd == k]
        resid(S[a.a], d, f"A sync, A dust with Q{k+1} replaced by {name}")
    d = Dd[a.a].copy(); d[qd >= 2] = r[name][qd >= 2]; resid(S[a.a], d, f"A sync, A dust with faint half replaced by {name}")
resid(r["bs_fixT"], Dd[a.a], "sync = bs_fixT, A dust"); resid(r["bs_fixT"], r["bd_fixT"], "sync = bs_fixT, dust = bd_fixT")
