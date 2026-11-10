"""E042: anatomy of the disagreement between WMAP K band and Planck 30 GHz in polarization.

    python -I nextgen/scripts/e042_k30_anatomy.py --cases base cut30 cutK worstK ctrlKa amp_global amp_region angle
        [--lcuts 3 5 8 12 20] [--out FILE.npz]

All fits: four BP XV regions, nside 16, window-matched to 60 arcmin, profile likelihood, CMB prior 0.2 uK,
dust from PR4 353 GHz. See the E042 record for what each case changes and for the expectations.
"""
import argparse
import gc
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.fit import fit_profile, sandwich_plugin  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import hfi_band, lfi_band, wmap_band  # noqa: E402
from lowres.window import FWHM_ARCMIN, WindowMatcher  # noqa: E402

NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
T0 = time.time()


def project_out(ninv, V):
    """N^-1 with the columns of V given infinite variance: N^-1 - N^-1 V (V^T N^-1 V)^-1 V^T N^-1."""
    NV = ninv @ V
    out = ninv - NV @ np.linalg.solve(V.T @ NV, NV.T)
    return 0.5 * (out + out.T)


def fit(bands, region, label, store, extra=""):
    L = RegionLikelihood(bands, region, cmb_prior_sd=0.2, use_dust=True)
    bh, res = fit_profile(L, np.full(4, -3.1))
    _, Sw, J, T, sol, ok = sandwich_plugin(L, bh)
    sd = np.sqrt(np.diag(Sw))
    rb = L.residuals(sol)
    n = region.size
    perband = "  ".join(f"{b.name} {float(r @ b.apply(r)) / n:.3f}" for b, r in zip(L.bands, rb))
    print(f"{label:28s} " + "  ".join(f"{b:+.4f} +- {s:.4f}" for b, s in zip(bh, sd)) + f" | chi2 {2 * res.fun:10.1f} | {perband} {extra}| "
          f"{time.time() - T0:.0f} s", flush=True)
    store[label] = np.concatenate([bh, sd, [2 * res.fun]])
    del L, sol
    gc.collect()
    return bh, sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="+", default=["base", "cut30", "cutK", "worstK", "ctrlKa", "amp_global", "amp_region", "angle"])
    ap.add_argument("--lcuts", nargs="+", type=int, default=[3, 5, 8, 12, 20])
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    WM = WindowMatcher(16, target_arcmin=60.0)
    Y, ell = WM.Y, WM.ell
    raw = {}
    for b in ("K", "Ka", "30", "353"):
        if b == "353":
            d, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=True)
        elif b == "30":
            d, ninv = D.load_lfi(30, 16)
            B = lfi_band(30)
        else:
            d, ninv = D.load_wmap(b)
            B = wmap_band(b, D.WMAP_BANDPASS)
        d, ninv = WM.apply(d, ninv, FWHM_ARCMIN[b])
        raw[b] = (B, d, ninv)
    region = D.bpxv_regions()
    n = region.size

    def band(b, d=None, ninv=None):
        B, d0, n0 = raw[b]
        return BandData(B.name, d0 if d is None else d, n0 if ninv is None else ninv, (lambda be, B=B: B.power_law(be, 23.0)),
                        B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF))

    store = {}
    print("# E042: profile index and debiased sandwich error per region (High Latitude, Spur, Galactic Center, Plane); chi2; chi2 per band and element")
    if "base" in a.cases:
        fit([band("K"), band("Ka"), band("353")], region, "K Ka 353", store)
        fit([band("K"), band("30"), band("353")], region, "K 30 353", store)
        fit([band("K"), band("Ka"), band("30"), band("353")], region, "K Ka 30 353", store)
    if "cut30" in a.cases:
        for L in a.lcuts:
            V = Y[:, ell <= L]
            fit([band("K"), band("30", ninv=project_out(raw["30"][2], V)), band("353")], region, f"K 30 353, 30 cut l<={L}", store,
                extra=f"modes {V.shape[1]} ")
    if "cutK" in a.cases:
        for L in a.lcuts:
            V = Y[:, ell <= L]
            fit([band("K", ninv=project_out(raw["K"][2], V)), band("30"), band("353")], region, f"K 30 353, K cut l<={L}", store,
                extra=f"modes {V.shape[1]} ")
    if "ctrlKa" in a.cases:
        for L in (8, 20):
            V = Y[:, ell <= L]
            fit([band("K"), band("Ka", ninv=project_out(raw["Ka"][2], V)), band("353")], region, f"K Ka 353, Ka cut l<={L}", store)
            fit([band("K", ninv=project_out(raw["K"][2], V)), band("Ka"), band("353")], region, f"K Ka 353, K cut l<={L}", store)
    if "worstK" in a.cases:
        cov = np.linalg.inv(raw["K"][2])
        w, U = np.linalg.eigh(0.5 * (cov + cov.T))
        sel = w > 10 * np.median(w)
        print(f"# K band: {int(sel.sum())} eigenmodes of the covariance with variance above 10 times the median (largest {w.max() / np.median(w):.0f} times)")
        V = U[:, sel]
        del cov, U
        fit([band("K", ninv=project_out(raw["K"][2], V)), band("30"), band("353")], region, "K 30 353, worst K modes out", store,
            extra=f"modes {V.shape[1]} ")
        fit([band("K", ninv=project_out(raw["K"][2], V)), band("Ka"), band("353")], region, "K Ka 353, worst K modes out", store)

    def joint(kind):
        """Bands K, Ka, 30, 353 with a nuisance on the 30 GHz map: an amplitude factor (global or per region)
        or a polarization angle. Joint L-BFGS in (beta, nuisance) with analytic gradients."""
        B30, d30, n30 = raw["30"]
        npar = {"amp_global": 1, "amp_region": 4, "angle": 1}[kind]
        hq = n // 2

        def transformed(z):
            if kind == "angle":
                c, s = np.cos(2 * z[0]), np.sin(2 * z[0])
                dq, du = d30[:hq], d30[hq:]
                dp = np.concatenate([c * dq + s * du, -s * dq + c * du])
                R = np.zeros((n, n))
                i = np.arange(hq)
                R[i, i], R[i, hq + i], R[hq + i, i], R[hq + i, hq + i] = c, s, -s, c
                return dp, R @ n30 @ R.T, None
            cvec = (np.full(n, z[0]) if kind == "amp_global" else z[region])
            return d30 / cvec, (cvec[:, None] * n30) * cvec[None, :], cvec

        def fun(x):
            beta, z = x[:4], x[4:]
            dp, npv, cvec = transformed(z)
            b30 = band("30", d=dp, ninv=npv)
            L = RegionLikelihood([band("K"), band("Ka"), b30, band("353")], region, cmb_prior_sd=0.2, use_dust=True)
            f, gb = L.profile(beta, grad=True)
            sol = L.solve(beta)
            r30 = L.residuals(sol)[2]
            w = npv @ r30                                   # N'^-1 r'
            m = dp - r30                                    # fitted model at 30 GHz
            if kind == "angle":
                Gm = 2 * np.concatenate([m[hq:], -m[:hq]])
                gz = np.array([float(w @ Gm)])
            else:
                # F = (1/2)(d - C m)^T N^-1 (d - C m): dF/dc_p = -m_p [N^-1 (d - C m)]_p = -m_p w_p / c_p
                per = -(w * m) / cvec
                gz = np.array([per.sum()]) if kind == "amp_global" else np.bincount(region, weights=per, minlength=4)
            del L, sol
            gc.collect()
            return f, np.concatenate([gb, gz])

        x0 = np.concatenate([np.full(4, -3.07), np.zeros(1) if kind == "angle" else np.ones(npar)])
        # self-check of the analytic gradient in the first nuisance parameter and the first index, at the start point
        f0, g0 = fun(x0)
        for j, hstep in ((4, 2e-3), (0, 2e-3)):
            e = np.zeros(x0.size)
            e[j] = hstep
            num = (fun(x0 + e)[0] - fun(x0 - e)[0]) / (2 * hstep)
            print(f"# {kind}: gradient check, parameter {j}: analytic {g0[j]:+.4f}, numerical {num:+.4f}", flush=True)
            if abs(g0[j] - num) > 2e-3 * max(1.0, abs(num)) + 0.05:
                raise SystemExit(f"gradient check failed for {kind}, parameter {j}")
        bnds = [(-4.6, -1.6)] * 4 + ([(-0.5, 0.5)] if kind == "angle" else [(0.5, 1.5)] * npar)
        res = minimize(fun, x0, jac=True, method="L-BFGS-B", bounds=bnds, options={"ftol": 1e-11, "gtol": 1e-7, "maxiter": 60})
        # curvature of the nuisance by central differences of its gradient with beta re-profiled is too costly here:
        # report the joint Hessian diagonal from finite differences of the gradient in the nuisance only (beta fixed)
        h = 2e-3
        H = np.zeros((npar, npar))
        for i in range(npar):
            e = np.zeros(4 + npar)
            e[4 + i] = h
            gp, gm = fun(res.x + e)[1], fun(res.x - e)[1]
            H[:, i] = (gp[4:] - gm[4:]) / (2 * h)
        sdz = np.sqrt(np.diag(np.linalg.inv(0.5 * (H + H.T))))
        z = res.x[4:]
        unit = (lambda v: np.degrees(v)) if kind == "angle" else (lambda v: v)
        print(f"{'K Ka 30 353, ' + kind:28s} " + "  ".join(f"{b:+.4f}" for b in res.x[:4]) + f" | chi2 {2 * res.fun:10.1f} | nuisance "
              + "  ".join(f"{unit(v):+.4f} +- {unit(s):.4f}" for v, s in zip(z, sdz)) + f" (sd at fixed indices) | {res.nit} it, converged {res.success} | "
              f"{time.time() - T0:.0f} s", flush=True)
        store["K Ka 30 353, " + kind] = np.concatenate([res.x, sdz, [2 * res.fun]])
        return res.x

    for kind in ("amp_global", "amp_region", "angle"):
        if kind in a.cases:
            x = joint(kind)
            if kind == "angle":
                # the K, 30 index with the angle fixed at its fitted value
                c, s = np.cos(2 * x[4]), np.sin(2 * x[4])
                hq = n // 2
                d30, n30 = raw["30"][1], raw["30"][2]
                dp = np.concatenate([c * d30[:hq] + s * d30[hq:], -s * d30[:hq] + c * d30[hq:]])
                R = np.zeros((n, n))
                i = np.arange(hq)
                R[i, i], R[i, hq + i], R[hq + i, i], R[hq + i, hq + i] = c, s, -s, c
                fit([band("K"), band("30", d=dp, ninv=R @ n30 @ R.T), band("353")], region, "K 30 353, angle fixed", store)
    if a.out:
        np.savez(a.out, **{k.replace(" ", "_").replace(",", "").replace("<=", "le"): v for k, v in store.items()})


if __name__ == "__main__":
    main()
