"""E044: the flat-prior marginal along the index of one region, at the configuration of the map fit.

For the fit of `e044_fit.py` (named regions, mask, K, Ka, Q, V and 353 GHz, nside 16) with a flat or the LCDM prior on
the CMB: the profile objective chi2_min / 2 and the flat-prior marginal objective chi2_min / 2 + (1/2) ln det A along
the index of region `--region`, the other indices held (a) at the profile optimum and (b) at the flat-prior optimum
for 23 GHz. For the reference at 30 GHz the marginal follows from the one at 23 GHz by the identity of the paper,
ln L(30) = ln L(23) + sum_r n_r beta_r ln(30 / 23); two direct evaluations at 30 GHz check that.

    python -I nextgen/scripts/e044_flat_scan.py --prior flat --fit nextgen/experiments/E044/flat.npz --out FILE.npz
"""
import argparse
import gc
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.mask import marginalize, restrict_prior  # noqa: E402
from lowres.regions import NAMES, named_regions  # noqa: E402
from lowres.spectra import hfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher, lcdm_cmb_precision  # noqa: E402

BANDS = ["K", "Ka", "Q", "V", "353"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", required=True, choices=["lcdm", "flat"])
    ap.add_argument("--fit", required=True, help="npz of the e044_fit.py run with this prior and the flat estimators")
    ap.add_argument("--region", type=int, default=0)
    ap.add_argument("--lo", type=float, default=-4.6)
    ap.add_argument("--hi", type=float, default=-1.6)
    ap.add_argument("--step", type=float, default=0.15)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    fit = np.load(a.fit, allow_pickle=True)
    region, keep, _ = named_regions(16, mask=True, bright_fraction=0.05)
    WM = WindowMatcher(16, target_arcmin=60.0)
    raw = []
    for b in BANDS:
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=True)
        else:
            d, ninv = D.load_wmap(b)
            B = wmap_band(b, D.WMAP_BANDPASS, drift=True)
        d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        d, ninv = marginalize(d, ninv, keep)
        raw.append((B, d, ninv, B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF)))
        gc.collect()
    region = region[keep]
    n, nreg = region.size, int(region.max()) + 1
    nper = np.bincount(region, minlength=nreg)
    assert np.array_equal(nper, fit["nper"]) and np.array_equal(keep, fit["keep"])
    if a.prior == "lcdm":
        _, cov_full = lcdm_cmb_precision(WM, floor_uk=0.01)
        kw = dict(cmb_prior_sd=None, cmb_prior_prec=restrict_prior(cov_full[np.ix_(keep, keep)], np.ones(n, dtype=bool)))
        del cov_full
    else:
        kw = dict(cmb_prior_sd=None)
    del WM
    gc.collect()

    def like(nu_ref):
        return RegionLikelihood([BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, nu_ref)), f) for B, d, ninv, f in raw], region, **kw)

    def objectives(L, beta):
        sol = L.solve(np.asarray(beta, float))
        out = 0.5 * float(sol["chi2"]), 0.5 * float(sol["chi2"]) + 0.5 * float(sol["logdet"])
        del sol
        gc.collect()
        return out

    r = a.region
    bp, b23 = np.array(fit["beta"], float), np.array(fit["flat23"], float)
    grid = np.round(np.arange(a.lo, a.hi + 1e-9, a.step), 6)
    print(f"# E044 flat-prior marginal along the index of {NAMES[r]} ({nper[r]} elements), CMB prior {a.prior}; profile optimum "
          + " ".join(f"{x:.4f}" for x in bp) + "; flat-prior optimum at 23 GHz " + " ".join(f"{x:.4f}" for x in b23)
          + "; at 30 GHz " + " ".join(f"{x:.4f}" for x in fit["flat30"]) + f"; setup {time.time() - t0:.0f} s", flush=True)
    L = like(23.0)
    p0, m0 = objectives(L, bp)
    tilt = lambda be: float(np.sum(nper * be) * np.log(30.0 / 23.0))          # ln L(30) - ln L(23)
    print(f"# at the profile optimum: profile objective {p0:.3f}, marginal objective at 23 GHz {m0:.3f}", flush=True)
    rows = []
    for tag, others in (("others at the profile optimum", bp), ("others at the flat-prior optimum for 23 GHz", b23)):
        print(f"# slice: {tag}. Columns: index; ln L_profile, ln L_flat(23 GHz), ln L_flat(30 GHz), each minus its value at the profile optimum", flush=True)
        for x in grid:
            be = others.copy()
            be[r] = x
            p, m = objectives(L, be)
            lp_, l23, l30 = -(p - p0), -(m - m0), -(m - m0) + tilt(be) - tilt(bp)
            rows.append((0 if others is bp else 1, x, lp_, l23, l30))
            print(f"{x:+.2f} {lp_:+12.3f} {l23:+12.3f} {l30:+12.3f}   [{time.time() - t0:.0f} s]", flush=True)
            if a.out:
                np.savez(a.out, rows=np.array(rows), profile_opt=bp, flat23=b23, flat30=fit["flat30"], nper=nper, region=r, prior=a.prior)
    del L
    gc.collect()
    L30 = like(30.0)
    for x in (grid[len(grid) // 4], grid[3 * len(grid) // 4]):
        be = bp.copy()
        be[r] = x
        _, m30 = objectives(L30, be)
        _, m30_0 = objectives(L30, bp)
        direct = -(m30 - m30_0)
        ident = [q[4] for q in rows if q[0] == 0 and abs(q[1] - x) < 1e-9][0]
        print(f"# identity check at {x:+.2f}: ln L_flat(30 GHz) direct {direct:+.4f}, from 23 GHz by the identity {ident:+.4f}, difference {direct - ident:+.2e}", flush=True)
    print(f"# done [{time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
