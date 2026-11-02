"""E027: noise-free region-fit references and template-S/N quartiles at a given index nside.

Writes refs_nside{N}.npz with bs_fixT, bd_fixT (region fits at fixed T), bs_trueT, bd_trueT,
qs16, qd16 (quartile labels at that nside, named as e026_analyze.py expects), region Fisher sds
and the within-region sd of the truth (representation error of a constant index per region).
usage: e027_refs.py --dir DIR --nside 8 4 [--out nextgen/experiments/E027]
"""
import argparse, os, sys, time
import numpy as np, healpy as hp
sys.path.insert(0, os.path.dirname(__file__))
from e023_exact_scan import FREQS, Sky, region_fit, fisher_sd_regions, quartiles   # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True); ap.add_argument("--nside", type=int, nargs="+", required=True)
    ap.add_argument("--out", default="nextgen/experiments/E027")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t = np.load(f"{a.dir}/truth.npz"); nside = int(t["nside"])
    d = np.array([hp.read_map(f"{a.dir}/band{int(f):03d}_map.fits", field=0) for f in FREQS]).T
    rms = np.array([hp.read_map(f"{a.dir}/band{int(f):03d}_rms.fits", field=0) for f in FREQS]).T
    skymap = np.asarray(t["sky"], float)
    if skymap.shape != d.shape:
        skymap = skymap.T
    sky = Sky(d, skymap, rms)
    T_fixed = float(t["T_fixed"]); T_true = np.asarray(t["T_d"], float)
    a_true = {"sync": np.asarray(t["a_sync"], float), "dust": np.asarray(t["a_dust"], float)}
    bs, bd, be3 = (np.asarray(t[k], float) for k in ("beta_s", "beta_d", "beta_eff3"))
    for nb in a.nside:
        nreg = hp.nside2npix(nb); reg = hp.ud_grade(np.arange(nreg, dtype=float), nside).astype(int)
        out = {}
        for sed in ("fixT", "trueT"):
            T = np.full(d.shape[0], T_fixed) if sed == "fixT" else T_true
            bs_mean = hp.ud_grade(bs, nb); bd_mean = hp.ud_grade(be3 if sed == "fixT" else bd, nb)
            t0 = time.time()
            out[f"bs_{sed}"], out[f"bd_{sed}"] = region_fit(skymap, rms, T, np.asarray(sky.cm), reg, nreg, bs_mean, bd_mean)
            print(f"nside {nb} {sed}: region fit in {time.time() - t0:.0f} s; fit minus region mean rms sync "
                  f"{np.sqrt(np.mean((out[f'bs_{sed}'] - bs_mean) ** 2)):.4f} dust {np.sqrt(np.mean((out[f'bd_{sed}'] - bd_mean) ** 2)):.4f}")
        # quartiles by the template (true-amplitude) region Fisher sd, as the E023 scans at nside 16
        fs = fisher_sd_regions("sync", a_true["sync"], bs, rms, None, reg, nreg)
        fd = fisher_sd_regions("dust", a_true["dust"], be3, rms, np.full(d.shape[0], T_fixed), reg, nreg)
        out["qs16"], out["qd16"] = quartiles(fs), quartiles(fd)
        out["fisher_sd_sync"], out["fisher_sd_dust"] = fs, fd
        # within-region sd of the truth (representation error of one index per region)
        for name, m in (("sync", bs), ("dust", be3)):
            mean_r = hp.ud_grade(m, nb)[reg]
            out[f"within_sd_{name}"] = np.sqrt(np.bincount(reg, (m - mean_r) ** 2, nreg) / np.bincount(reg, minlength=nreg))
        for k in range(4):
            print(f"  nside {nb} Q{k+1}: sync Fisher sd median {np.median(fs[out['qs16'] == k]):.3f} within-region truth sd median "
                  f"{np.median(out['within_sd_sync'][out['qs16'] == k]):.4f}; dust {np.median(fd[out['qd16'] == k]):.3f} / "
                  f"{np.median(out['within_sd_dust'][out['qd16'] == k]):.4f}")
        np.savez(f"{a.out}/refs_nside{nb}.npz", **out)
        print(f"  wrote {a.out}/refs_nside{nb}.npz")


if __name__ == "__main__":
    main()
