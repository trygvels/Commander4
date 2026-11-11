"""E040 stage 5: block-robust errors of the region index with leverage corrections, on simulations and on the maps.

At the profile optimum the score of region r is a sum over elements, g_p = -sum_b (N_b^-1 r_b)_p (ds_b a_s)_p,
which is summed over blocks of sky B (the 48 nside-2 pixels) to G_{B,r} (g is the gradient of chi2_min / 2, so
its derivative is the information). With D_{B,r} = dG_{B,r}/dbeta_r (the information that block B carries about
the index of its region, by a central difference with the amplitudes re-solved; sum_B D_B / J_c is printed as a check) and the leverage h_B = D_B / sum_B D_B:

    CR0   cov = J_c^-1 (sum_B G_B G_B^T) J_c^-1                      (biased low when few blocks carry the information)
    CR2   the same with G_B / sqrt(1 - h_B)                          (unbiased under the model for equal-variance scores)
    CR3   the same with G_B / (1 - h_B)                              (the delete-one-block jackknife to first order; conservative)

n_eff = (sum_B D_B)^2 / sum_B D_B^2 is the effective number of blocks. Leverages are clipped to [0, 0.9].

    python -I nextgen/scripts/e040f_block_errors_v2.py sim --nreal 60 [--spread 0.1]
    python -I nextgen/scripts/e040f_block_errors_v2.py maps [--cmb-prior-lcdm] [--out FILE.npz]
"""
import argparse
import os
import sys
import time

import healpy as hp
import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D_  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import Band, g_rj, hfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher  # noqa: E402

DR1 = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
KINDS = ("CR0", "CR2", "CR3")


def block_scores(L, beta, block, nblk):
    sol = L.solve(beta)
    n = L.n
    g = np.zeros(n)
    for i, (b, r) in enumerate(zip(L.bands, L.residuals(sol))):
        g -= b.apply(r) * sol["ds"][i] * sol["x"][:n]
    G = np.zeros((nblk, L.nreg))
    np.add.at(G, (block, L.region), g)
    return G


def block_errors(L, bh, block, nblk, h=2e-3):
    """Returns sd_model (nreg,), {kind: sd (nreg,)}, n_eff (nreg,), max leverage (nreg,), sum of D over J_c."""
    _, cs, Jp, T, _, ok = sandwich_plugin(L, bh)
    G = block_scores(L, bh, block, nblk)
    Dm = np.zeros_like(G)
    for r in range(L.nreg):
        e = np.zeros(L.nreg)
        e[r] = h
        Dm[:, r] = (block_scores(L, bh + e, block, nblk)[:, r] - block_scores(L, bh - e, block, nblk)[:, r]) / (2 * h)
    Dpos = np.clip(Dm, 0.0, None)
    lev = np.clip(Dpos / Dpos.sum(0)[None, :], 0.0, 0.9)
    neff = Dpos.sum(0) ** 2 / (Dpos ** 2).sum(0)
    Ji = np.linalg.inv(Jp - T)
    out = {}
    for kind, f in zip(KINDS, (np.ones_like(lev), 1 / np.sqrt(1 - lev), 1 / (1 - lev))):
        Gk = G * f
        out[kind] = np.sqrt(np.diag(Ji @ (Gk.T @ Gk) @ Ji))
    return np.sqrt(np.diag(cs)), out, neff, lev.max(0), Dm.sum(0) / np.diag(Jp - T)


def truth_sync_30(nside):
    m = np.array(hp.read_map(DR1, field=(0, 1, 2)), dtype=float)
    m = hp.smoothing(hp.ud_grade(m, 64), fwhm=np.radians(5.0), pol=True)
    qu = hp.ud_grade(m[1:], nside)
    return np.concatenate([hp.reorder(qu[0], r2n=True), hp.reorder(qu[1], r2n=True)])


def run_sim(a):
    ns = 8
    n = 2 * hp.nside2npix(ns)
    raw = []
    for b in ("K", "Ka", "Q", "V"):
        d, ninv = D_.load_wmap(b)
        d, ninv = D_.degrade(d, ninv, 16, ns)
        raw.append((f"WMAP_{b}", wmap_band(b, D_.WMAP_BANDPASS), d, ninv))
    d3, n3 = D_.load_planck_qu(D_.PR4_353, ns)
    raw.append(("Planck_353", Band.delta("Planck_353", 353.0), d3, n3))
    reg16 = D_.bpxv_regions()
    r = reg16[:reg16.size // 2].reshape(-1, 4)
    region = np.tile(np.array([np.bincount(row, minlength=4).argmax() for row in r]), 2)
    block = np.tile(np.arange(n // 2) // 16, 2)
    nblk = np.array([np.unique(block[region == k]).size for k in range(4)])
    beta_reg = np.array([-3.15, -3.25, -2.85, -2.95])
    a30 = truth_sync_30(ns)
    a_d = d3 * g_rj(353.0)
    fd = [B.mbb(D_.DUST_BETA, D_.DUST_T, D_.DUST_REF) for _, B, _, _ in raw]
    chol = [np.linalg.cholesky(ninv) for *_, ninv in raw]

    def simulate(r_):
        off = a.spread * np.random.default_rng([a.seed, 3, r_]).standard_normal(48)
        beta_el = beta_reg[region] + off[block]
        a_s = a30 * (23.0 / 30.0) ** beta_el
        rng = np.random.default_rng([a.seed, 7, r_])
        a_c = 0.2 * rng.standard_normal(n)
        return [B.power_law(beta_el, 23.0)[0] * a_s + f * a_d + a_c + sl.solve_triangular(L.T, rng.standard_normal(n), lower=False)
                for (_, B, _, _), f, L in zip(raw, fd, chol)]

    print(f"# E040 stage 5, simulation check: nside 8, blocks per region {nblk.tolist()}, true index spread between blocks {a.spread} "
          "(a new offset per block in every realization)", flush=True)
    t0 = time.time()
    est, sm, sb, ne = [], [], {k: [] for k in KINDS}, []
    for r_ in range(a.nreal):
        dv = simulate(r_)
        bands = [BandData(nm, d_, ninv, (lambda be, B=B: B.power_law(be, 23.0)), f) for (nm, B, _, ninv), d_, f in zip(raw, dv, fd)]
        L = RegionLikelihood(bands, region, cmb_prior_sd=0.2)
        bh, _ = fit_profile(L, beta_reg)
        s_model, s_blk, neff, levmax, dsum = block_errors(L, bh, block, 48)
        est.append(bh); sm.append(s_model); ne.append(neff)
        for k in KINDS:
            sb[k].append(s_blk[k])
        if r_ % 10 == 0:
            print("real", r_, "index", " ".join(f"{x:+.4f}" for x in bh), "| CR0/CR2/CR3 over sandwich, High Latitude",
                  " ".join(f"{s_blk[k][0] / s_model[0]:.2f}" for k in KINDS), "Plane", " ".join(f"{s_blk[k][3] / s_model[3]:.2f}" for k in KINDS),
                  "| n_eff", " ".join(f"{x:.1f}" for x in neff), "| sum D / J_c", " ".join(f"{x:.2f}" for x in dsum), f"[{time.time() - t0:.0f} s]", flush=True)
    est, sm, ne = np.array(est), np.array(sm), np.array(ne)
    sb = {k: np.array(v) for k, v in sb.items()}
    scat = np.sqrt(np.mean((est - beta_reg) ** 2, axis=0))            # rms about the region truth
    print(f"## over {len(est)} realizations: rms of the estimate about the region truth, against the median quoted errors")
    print(f"{'region':16s} {'blocks':>6s} {'n_eff':>6s} {'rms':>8s} {'sandwich':>9s} " + " ".join(f"{k:>8s}" for k in KINDS)
          + "   rms of z with: sandwich " + " ".join(f"{k:>6s}" for k in KINDS))
    for k_ in range(4):
        zs = [np.sqrt(np.mean(((est[:, k_] - beta_reg[k_]) / s[:, k_]) ** 2)) for s in (sm, sb["CR0"], sb["CR2"], sb["CR3"])]
        print(f"{NAMES[k_]:16s} {nblk[k_]:6d} {np.median(ne[:, k_]):6.1f} {scat[k_]:8.4f} {np.median(sm[:, k_]):9.4f} "
              + " ".join(f"{np.median(sb[k][:, k_]):8.4f}" for k in KINDS) + "                    " + " ".join(f"{z:6.2f}" for z in zs))


def run_maps(a):
    t0 = time.time()
    WM = WindowMatcher(16, target_arcmin=60.0)
    bands = []
    for b in ("K", "Ka", "Q", "V", "353"):
        if b == "353":
            d, ninv = D_.load_planck_qu(D_.PR4_353, 16)
            B = hfi_band(353, psb=True)
        else:
            d, ninv = D_.load_wmap(b)
            B = wmap_band(b, D_.WMAP_BANDPASS)
        d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        bands.append(BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, 23.0)), B.mbb(D_.DUST_BETA, D_.DUST_T, D_.DUST_REF)))
    region = D_.bpxv_regions()
    n = region.size
    kw = dict(cmb_prior_sd=0.2)
    if a.cmb_prior_lcdm:
        from lowres.window import lcdm_cmb_precision
        prec, _ = lcdm_cmb_precision(WM)
        kw = dict(cmb_prior_sd=None, cmb_prior_prec=prec)
    L = RegionLikelihood(bands, region, use_dust=True, **kw)
    bh, res = fit_profile(L, np.full(4, -3.1))
    block = np.tile(np.arange(n // 2) // 64, 2)
    nblk = np.array([np.unique(block[region == k]).size for k in range(4)])
    s_model, s_blk, neff, levmax, dsum = block_errors(L, bh, block, 48)
    print(f"# E040 stage 5, maps: K, Ka, Q, V, 353, window-matched, nside 16, CMB prior {'LCDM' if a.cmb_prior_lcdm else 'white 0.2 uK'}; {time.time() - t0:.0f} s")
    print(f"{'region':16s} {'blocks':>6s} {'n_eff':>6s} {'max h':>6s} {'sum D/J_c':>9s} {'index':>9s} {'sandwich':>9s} " + " ".join(f"{k:>8s}" for k in KINDS)
          + "  ratio to sandwich: " + " ".join(f"{k:>6s}" for k in KINDS))
    for k_ in range(4):
        print(f"{NAMES[k_]:16s} {nblk[k_]:6d} {neff[k_]:6.1f} {levmax[k_]:6.2f} {dsum[k_]:9.2f} {bh[k_]:+9.4f} {s_model[k_]:9.4f} "
              + " ".join(f"{s_blk[k][k_]:8.4f}" for k in KINDS) + "                     " + " ".join(f"{s_blk[k][k_] / s_model[k_]:6.2f}" for k in KINDS))
    if a.out:
        np.savez(a.out, beta=bh, sd_model=s_model, neff=neff, **{k: s_blk[k] for k in KINDS})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["sim", "maps"])
    ap.add_argument("--nreal", type=int, default=60)
    ap.add_argument("--seed", type=int, default=39)
    ap.add_argument("--spread", type=float, default=0.0)
    ap.add_argument("--cmb-prior-lcdm", action="store_true")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    run_sim(a) if a.mode == "sim" else run_maps(a)


if __name__ == "__main__":
    main()
