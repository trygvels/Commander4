"""E040 and E039 at the configuration of the map fit, in stages that each fit the machine's memory.

Simulations at nside 16 with the window-matched WMAP and 353 GHz covariances, a LCDM CMB and one true index per
region (the model is true). `e040h_sim_n16.py` did all of this in one process and needed 21 GB; here the data are
drawn once and saved, the window-matched covariances are cached on disk, and each kind of fit is a separate run.

    draw      draw and save the data vectors of realizations 0 to nreal - 1 (and cache the covariances)
    blocknull LCDM-prior profile fit, sandwich error and block errors for each realization
    white     white-prior (0.2 uK) profile fit and sandwich error
    lowpass   LCDM-prior fit with Ka, Q, V through multipoles 2 to L only, and its sandwich error

    python -I nextgen/scripts/e040j_sim_n16_stages.py draw --nreal 12 --dir ~/data/c4_tools/sim_n16
    python -I nextgen/scripts/e040j_sim_n16_stages.py blocknull --nreal 12 --dir ~/data/c4_tools/sim_n16 --out FILE.npz
"""
import argparse
import gc
import os
import sys
import time

import healpy as hp
import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.block import KINDS, block_errors  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import g_rj, hfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher, lcdm_cmb_precision  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
BANDS = ("K", "Ka", "Q", "V", "353")
BETA_TRUE = np.array([-3.15, -3.25, -2.85, -2.95])


def band_objects():
    out = []
    for b in BANDS:
        B = hfi_band(353, psb=True) if b == "353" else wmap_band(b, D.WMAP_BANDPASS)
        out.append(B)
    return out


def geometry(ns):
    region = D.bpxv_regions()
    if ns != 16:
        rr = region[:region.size // 2].reshape(-1, (16 // ns) ** 2)
        region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in rr]), 2)
    n = region.size
    block = np.tile(np.arange(n // 2) // (ns // 2) ** 2, 2)
    return region, block


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def draw(a):
    """Window-matched inverse covariances to disk, then the data of every realization to one file."""
    ns = a.nside
    os.makedirs(a.dir, exist_ok=True)
    WM = WindowMatcher(ns, target_arcmin=60.0)
    Bs = band_objects()
    region, _ = geometry(ns)
    n = region.size
    d353 = None
    roots = []
    for b in BANDS:
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, ns)
            d353 = d.copy()
        else:
            d, ninv = D.load_wmap(b)
            if ns != 16:
                d, ninv = D.degrade(d, ninv, 16, ns)
        _, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        np.save(os.path.join(a.dir, f"ninv_{b}.npy"), ninv)
        roots.append(np.linalg.cholesky(ninv))               # N^-1 = L L^T; noise = L^-T xi
        del ninv
        gc.collect()
    prec, cov_c = lcdm_cmb_precision(WM)
    np.save(os.path.join(a.dir, "cmb_prec.npy"), prec)
    Lc = np.linalg.cholesky(cov_c)
    del prec, cov_c
    a_s = truth_sync_30(ns) * (23.0 / 30.0) ** BETA_TRUE[region]
    a_d = d353 * g_rj(353.0)
    fd = [B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF) for B in Bs]
    s_true = [B.power_law(BETA_TRUE, 23.0)[0][region] for B in Bs]
    data = np.zeros((a.nreal, len(BANDS), n))
    for r_ in range(a.nreal):
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = Lc @ np.random.default_rng([a.seed, 11, r_]).standard_normal(n)
        for i, (s0, f, L) in enumerate(zip(s_true, fd, roots)):
            data[r_, i] = s0 * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False)
    keep = WM.ell <= a.lmax
    Yk, Ypk = WM.Y[:, keep], WM.Yp[keep]
    err_id = float(np.abs(Ypk @ Yk - np.eye(Yk.shape[1])).max())
    assert err_id < 1e-8, err_id
    np.savez(os.path.join(a.dir, "sim.npz"), data=data, fd=np.array(fd), Yk=Yk, Ypk=Ypk, lmax=a.lmax, nside=ns, seed=a.seed)
    print(f"# drew {a.nreal} realizations at nside {ns} into {a.dir}; low-pass cut {a.lmax} ({int(keep.sum())} of {WM.Y.shape[1]} modes, "
          f"|Y_k^+ Y_k - 1| {err_id:.1e})", flush=True)


def load(a, lowpass=False):
    z = np.load(os.path.join(a.dir, "sim.npz"))
    ns = int(z["nside"])
    region, block = geometry(ns)
    Bs = band_objects()
    ninv = [np.load(os.path.join(a.dir, f"ninv_{b}.npy")) for b in BANDS]
    low = {}
    if lowpass:
        Yk, Ypk = z["Yk"], z["Ypk"]
        for i in (1, 2, 3):
            Cc = Ypk @ np.linalg.inv(ninv[i]) @ Ypk.T
            nl = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
            low[i] = 0.5 * (nl + nl.T)
            ninv[i] = None                                    # the full covariance of a low-passed band is not needed
            del Cc, nl
            gc.collect()
    return z, ns, region, block, Bs, ninv, low


def bands_for(z, Bs, ninv, low, r_):
    out = []
    for i, B in enumerate(Bs):
        d_ = z["data"][r_, i]
        if i in low:
            out.append(BandData(B.name, z["Yk"] @ (z["Ypk"] @ d_), low[i], (lambda be, B=B: B.power_law(be, 23.0)), float(z["fd"][i])))
        else:
            out.append(BandData(B.name, d_, ninv[i], (lambda be, B=B: B.power_law(be, 23.0)), float(z["fd"][i])))
    return out


def fits(a):
    t0 = time.time()
    z, ns, region, block, Bs, ninv, low = load(a, lowpass=(a.mode == "lowpass"))
    if a.mode == "white":
        kw = dict(cmb_prior_sd=0.2)
    else:
        kw = dict(cmb_prior_sd=None, cmb_prior_prec=np.load(os.path.join(a.dir, "cmb_prec.npy")))
    nreal = min(a.nreal, z["data"].shape[0])
    print(f"# E040 nside-{ns} simulation, stage {a.mode}: truth {BETA_TRUE}; {nreal} realizations; setup {time.time() - t0:.0f} s", flush=True)
    keys = ["b", "s"] + (["neff", "levmax", "nneg", "cross", "dj"] + list(KINDS) if a.mode == "blocknull" else [])
    store = {k: [] for k in keys}
    for r_ in range(a.first, nreal):
        L = RegionLikelihood(bands_for(z, Bs, ninv, low, r_), region, **kw)
        bh, res = fit_profile(L, BETA_TRUE)
        if a.mode == "blocknull":
            out = block_errors(L, bh, block, 48)
            sd = out["sd_model"]
            store["neff"].append(out["neff"]); store["levmax"].append(out["lev_max"]); store["nneg"].append(out["n_neg"])
            store["cross"].append(out["cross"]); store["dj"].append(out["d_over_jc"])
            for k in KINDS:
                store[k].append(out[k])
            extra = (" | CR3/sd " + " ".join(f"{x:.2f}" for x in out["CR3"] / sd) + " | CR0/sd " + " ".join(f"{x:.2f}" for x in out["CR0"] / sd)
                     + " | n_eff " + " ".join(f"{x:.1f}" for x in out["neff"]) + f" | neg D {int(out['n_neg'].sum())} cross {out['cross']:.3f} sumD/Jc "
                     + " ".join(f"{x:.2f}" for x in out["d_over_jc"]))
            del out
        else:
            res_s = sandwich_plugin(L, bh)
            sd = np.sqrt(np.diag(res_s[1]))
            del res_s
            extra = ""
        del L
        gc.collect()
        store["b"].append(bh); store["s"].append(sd)
        print(f"real {r_} {a.mode} " + " ".join(f"{x:+.4f}" for x in bh) + " sd " + " ".join(f"{x:.4f}" for x in sd) + extra
              + f" | conv {res.success} [{time.time() - t0:.0f} s]", flush=True)
        if a.out:
            np.savez(a.out, truth=BETA_TRUE, first=a.first, **{k: np.array(v) for k, v in store.items()})
    S = {k: np.array(v) for k, v in store.items()}
    zz = (S["b"] - BETA_TRUE) / S["s"]
    print(f"## {a.mode}, {len(S['b'])} realizations: mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (S['b'] - BETA_TRUE).mean(0))
          + f"; pooled mean z {zz.mean():+.3f} +- {1 / np.sqrt(zz.size):.3f}, sd z {zz.std(ddof=1):.3f}, coverage of 1 sd {np.mean(np.abs(zz) < 1):.3f}; "
          + "median error " + " ".join(f"{x:.4f}" for x in np.median(S["s"], axis=0)))
    if a.mode == "blocknull":
        print("## block error over sandwich error under the model, per region: median, maximum")
        for k in KINDS:
            q = S[k] / S["s"]
            print(f"   {k}: " + "; ".join(f"{NAMES[j]} {np.nanmedian(q[:, j]):.2f} {np.nanmax(q[:, j]):.2f}" for j in range(4)))
        print("   effective number of blocks, median: " + " ".join(f"{x:.1f}" for x in np.median(S["neff"], axis=0))
              + f"; blocks with negative information, total {int(S['nneg'].sum())}; cross-region share, largest {S['cross'].max():.3f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["draw", "blocknull", "white", "lowpass"])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--nreal", type=int, default=12)
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--lmax", type=int, default=20)
    ap.add_argument("--seed", type=int, default=39)
    ap.add_argument("--nside", type=int, default=16, choices=[8, 16], help="8 is for an end-to-end test of this script only")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    a.dir = os.path.expanduser(a.dir)
    draw(a) if a.mode == "draw" else fits(a)


if __name__ == "__main__":
    main()
