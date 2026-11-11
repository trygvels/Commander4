"""E040 and E039 at the configuration of the map fit: simulations at nside 16 with the window-matched WMAP and
353 GHz covariances, a LCDM CMB, and one true index per region (the model is true).

For each realization, on the same simulated maps:
  1. the all-scales profile fit with the LCDM prior, its debiased sandwich error, and the block errors of
     `lowres.block` over the 48 nside-2 pixels (the null distribution of block error over sandwich error that the
     independent review of 2026-11-11 asked for);
  and for the first `--full` realizations also
  2. the same fit with the white 0.2 uK prior (the paired difference between the two priors, E040 C4 at nside 16);
  3. the fit that reads the index from multipoles 2 to L only (K band and 353 GHz at all scales; E039 addendum B
     at nside 16).

    python -I nextgen/scripts/e040h_sim_n16.py --nreal 12 --full 6 --lmax 20 --out FILE.npz
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


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nreal", type=int, default=12)
    ap.add_argument("--full", type=int, default=6, help="realizations that also get the white-prior and the low-pass fits")
    ap.add_argument("--first", type=int, default=0)
    ap.add_argument("--lmax", type=int, default=20)
    ap.add_argument("--seed", type=int, default=39)
    ap.add_argument("--out", default="")
    ap.add_argument("--nside", type=int, default=16, choices=[8, 16], help="8 is for a quick end-to-end test of this script only")
    a = ap.parse_args()
    t0 = time.time()
    ns = a.nside
    WM = WindowMatcher(ns, target_arcmin=60.0)
    raw = []                                                    # (name, Band, inverse covariance as in the map fit)
    d353 = None
    for b in ("K", "Ka", "Q", "V", "353"):
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, ns)
            B = hfi_band(353, psb=True)
            d353 = d.copy()
        else:
            d, ninv = D.load_wmap(b)
            if ns != 16:
                d, ninv = D.degrade(d, ninv, 16, ns)
            B = wmap_band(b, D.WMAP_BANDPASS)
        _, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        raw.append((B.name, B, ninv))
    region = D.bpxv_regions()
    if ns != 16:
        rr = region[:region.size // 2].reshape(-1, (16 // ns) ** 2)
        region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in rr]), 2)
    n = region.size
    block = np.tile(np.arange(n // 2) // (ns // 2) ** 2, 2)       # nside-2 pixel of each pixel (NESTED), Q and U together
    beta_true = np.array([-3.15, -3.25, -2.85, -2.95])
    a_s = truth_sync_30(ns) * (23.0 / 30.0) ** beta_true[region]
    a_d = d353 * g_rj(353.0)
    fd = [B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF) for _, B, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]       # N^-1 = L L^T; noise = L^-T xi
    prec, cov_c = lcdm_cmb_precision(WM)
    Lc = np.linalg.cholesky(cov_c)
    del cov_c
    keep = WM.ell <= a.lmax
    Yk, Ypk = WM.Y[:, keep], WM.Yp[keep]
    err_id = float(np.abs(Ypk @ Yk - np.eye(Yk.shape[1])).max())
    assert err_id < 1e-8, err_id
    low_ninv = {}
    if a.full > 0:
        for (nm, B, ninv) in raw[1:4]:
            Cc = Ypk @ np.linalg.inv(ninv) @ Ypk.T
            nl = Ypk.T @ np.linalg.inv(0.5 * (Cc + Cc.T)) @ Ypk
            low_ninv[nm] = 0.5 * (nl + nl.T)
            del Cc, nl
            gc.collect()
    s_true = [B.power_law(beta_true, 23.0)[0][region] for _, B, _ in raw]
    print(f"# E040 nside-16 simulation at the configuration of the map fit: window-matched covariances, LCDM CMB, truth {beta_true}; "
          f"low-pass cut {a.lmax} ({int(keep.sum())} of {WM.Y.shape[1]} modes, |Y_k^+ Y_k - 1| {err_id:.1e}); setup {time.time() - t0:.0f} s", flush=True)

    def simulate(r_):
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = Lc @ np.random.default_rng([a.seed, 11, r_]).standard_normal(n)
        return [s0 * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False)
                for s0, f, L in zip(s_true, fd, chol)]

    def bands_for(dv, lowpass=False):
        out = []
        for (nm, B, ninv), d_, f in zip(raw, dv, fd):
            if lowpass and nm in low_ninv:
                out.append(BandData(nm, Yk @ (Ypk @ d_), low_ninv[nm], (lambda be, B=B: B.power_law(be, 23.0)), f))
            else:
                out.append(BandData(nm, d_, ninv, (lambda be, B=B: B.power_law(be, 23.0)), f))
        return out

    keys = ["b_l", "s_l", "neff", "levmax", "nneg", "cross", "dj"] + list(KINDS) + ["b_w", "s_w", "b_p", "s_p"]
    store = {k: [] for k in keys}
    for r_ in range(a.first, a.first + a.nreal):
        dv = simulate(r_)
        L = RegionLikelihood(bands_for(dv), region, cmb_prior_sd=None, cmb_prior_prec=prec)
        bh, res = fit_profile(L, beta_true)
        out = block_errors(L, bh, block, 48)
        del L
        gc.collect()
        store["b_l"].append(bh); store["s_l"].append(out["sd_model"]); store["neff"].append(out["neff"])
        store["levmax"].append(out["lev_max"]); store["nneg"].append(out["n_neg"]); store["cross"].append(out["cross"])
        store["dj"].append(out["d_over_jc"])
        for k in KINDS:
            store[k].append(out[k])
        line = (f"real {r_} lcdm " + " ".join(f"{x:+.4f}" for x in bh) + " sd " + " ".join(f"{x:.4f}" for x in out["sd_model"])
                + " | CR3/sd " + " ".join(f"{x:.2f}" for x in out["CR3"] / out["sd_model"]) + " | CR0/sd " + " ".join(f"{x:.2f}" for x in out["CR0"] / out["sd_model"])
                + " | n_eff " + " ".join(f"{x:.1f}" for x in out["neff"]) + f" | neg D {int(out['n_neg'].sum())} cross {out['cross']:.3f} sumD/Jc "
                + " ".join(f"{x:.2f}" for x in out["d_over_jc"]) + f" | conv {res.success}")
        if r_ - a.first < a.full:
            L = RegionLikelihood(bands_for(dv), region, cmb_prior_sd=0.2)
            bw, resw = fit_profile(L, beta_true)
            _, cw, _, _, _, _ = sandwich_plugin(L, bw)
            del L
            gc.collect()
            L = RegionLikelihood(bands_for(dv, lowpass=True), region, cmb_prior_sd=None, cmb_prior_prec=prec)
            bp, resp = fit_profile(L, beta_true)
            _, cp, _, _, _, _ = sandwich_plugin(L, bp)
            del L
            gc.collect()
            store["b_w"].append(bw); store["s_w"].append(np.sqrt(np.diag(cw))); store["b_p"].append(bp); store["s_p"].append(np.sqrt(np.diag(cp)))
            line += (" | white " + " ".join(f"{x:+.4f}" for x in bw) + " sd " + " ".join(f"{x:.4f}" for x in np.sqrt(np.diag(cw))) + f" conv {resw.success}"
                     + f" | low-pass " + " ".join(f"{x:+.4f}" for x in bp) + " sd " + " ".join(f"{x:.4f}" for x in np.sqrt(np.diag(cp))) + f" conv {resp.success}")
        print(line + f" [{time.time() - t0:.0f} s]", flush=True)
        if a.out:
            np.savez(a.out, truth=beta_true, **{k: np.array(v) for k, v in store.items()})
    S = {k: np.array(v) for k, v in store.items()}
    m = len(S["b_l"])
    z = (S["b_l"] - beta_true) / S["s_l"]
    print(f"## LCDM prior, all scales, {m} realizations: mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (S['b_l'] - beta_true).mean(0))
          + f"; pooled mean z {z.mean():+.3f} +- {1 / np.sqrt(z.size):.3f}, sd z {z.std(ddof=1):.3f}, coverage of 1 sd {np.mean(np.abs(z) < 1):.3f}")
    print("## block error over sandwich error under the model, per region: median, maximum")
    for k in KINDS:
        q = S[k] / S["s_l"]
        print(f"   {k}: " + "; ".join(f"{NAMES[j]} {np.nanmedian(q[:, j]):.2f} {np.nanmax(q[:, j]):.2f}" for j in range(4)))
    print("   effective number of blocks, median: " + " ".join(f"{x:.1f}" for x in np.median(S["neff"], axis=0))
          + f"; blocks with negative information, total {int(S['nneg'].sum())}; cross-region share, largest {S['cross'].max():.3f}")
    if len(S["b_w"]):
        k_ = len(S["b_w"])
        dlt = S["b_w"] - S["b_l"][:k_]
        zw = (S["b_w"] - beta_true) / S["s_w"]
        print(f"## white prior on the same {k_} realizations: pooled mean z {zw.mean():+.3f}, sd z {zw.std(ddof=1):.3f}; white minus LCDM, paired: mean "
              + " ".join(f"{x:+.4f}" for x in dlt.mean(0)) + "; sd " + " ".join(f"{x:.4f}" for x in dlt.std(0, ddof=1)))
        zp = (S["b_p"] - beta_true) / S["s_p"]
        dp = S["b_p"] - S["b_l"][:k_]
        print(f"## index from l <= {a.lmax} on the same {k_} realizations: mean estimate - truth " + " ".join(f"{x:+.4f}" for x in (S['b_p'] - beta_true).mean(0))
              + f"; pooled mean z {zp.mean():+.3f} +- {1 / np.sqrt(zp.size):.3f}, sd z {zp.std(ddof=1):.3f}, coverage of 1 sd {np.mean(np.abs(zp) < 1):.3f}; "
              + "median error " + " ".join(f"{x:.4f}" for x in np.median(S["s_p"], axis=0)) + " against all scales " + " ".join(f"{x:.4f}" for x in np.median(S["s_l"][:k_], axis=0))
              + "; low-pass minus all scales, paired: mean " + " ".join(f"{x:+.4f}" for x in dp.mean(0)) + "; sd " + " ".join(f"{x:.4f}" for x in dp.std(0, ddof=1))
              + "; expected sd sqrt(sd_low^2 - sd_all^2) " + " ".join(f"{x:.4f}" for x in np.sqrt(np.clip(np.median(S["s_p"], axis=0) ** 2 - np.median(S["s_l"][:k_], axis=0) ** 2, 0, None))))


if __name__ == "__main__":
    main()
