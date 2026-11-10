"""Check 6: spectra.py conventions and quadrature. Run with `python -I` (reads downloaded bandpass files).

(a) delta band gives f(nu0)/g(nu0) for the power law and the modified blackbody (independent g and constants).
(b) power_law derivative vs finite differences; shapes.
(c) np.gradient weights vs trapezoid vs a converged reference, smooth bandpass on a non-uniform grid, with the
    response (i) vanishing at the ends, (ii) truncated at 5 and 30 percent of the peak.
(d) the real WMAP / LFI / HFI bandpasses: end-point response, and the change of the band-averaged spectra when
    np.gradient weights are replaced by trapezoid weights.
(e) effective_frequency: defining property and resolution.
(f) wmap_band: unit-area per detector (code) vs raw sum (what the docstring says).
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/Users/trygvels/code/Commander4/nextgen")
from lowres import spectra as S  # noqa: E402
from rv_common import HK, TCMB, dust_ref, g_ref, sync_ref  # noqa: E402

np.set_printoptions(precision=6, linewidth=160)
print("== constants ==")
print(f"H_OVER_K code {S.H_OVER_K!r} vs CODATA 2018 {HK!r}: rel diff {S.H_OVER_K / HK - 1:.1e}; g(353) rel diff "
      f"{S.g_rj(353.0) / g_ref(353.0) - 1:.1e}")

print("\n== (a) delta band ==")
for nu0 in (22.8, 44.1, 70.4, 353.0):
    B = S.Band.delta("x", nu0)
    s, ds = B.power_law(np.array([-3.1, -2.7]), 23.0)
    sr, dsr = sync_ref(nu0, np.array([-3.1, -2.7]), 23.0)
    print(f"nu0 {nu0:6.1f}: power law rel diff {np.max(np.abs(s / sr - 1)):.1e}, derivative rel diff {np.max(np.abs(ds / dsr - 1)):.1e}, "
          f"mbb rel diff {abs(B.mbb(1.55, 19.6, 353.0) / dust_ref(nu0) - 1):.1e}, cmb column (f=g) {B.average(S.g_rj(B.nu)):.15f}")

print("\n== (b) broad band: derivative and shapes ==")
nu = np.linspace(18.0, 28.0, 201)
resp = np.exp(-0.5 * ((nu - 23.0) / 2.0) ** 2)
B = S.Band("g", nu, resp)
beta = np.array([-3.3, -3.1, -2.6])
s, ds = B.power_law(beta, 23.0)
h = 1e-6
fdd = (B.power_law(beta + h, 23.0)[0] - B.power_law(beta - h, 23.0)[0]) / (2 * h)
print("shapes", s.shape, ds.shape, " derivative vs FD rel diff", np.max(np.abs(ds / fdd - 1)))
print("power_law(scalar) shape", B.power_law(-3.1, 23.0)[0].shape, " CMB column", B.average(S.g_rj(B.nu)))


def band_ratio(nu, resp, f, kind):
    if kind == "gradient":
        w = resp * np.gradient(nu)
        return np.sum(w * f(nu)) / np.sum(w * S.g_rj(nu))
    return np.trapezoid(resp * f(nu), nu) / np.trapezoid(resp * S.g_rj(nu), nu)


print("\n== (c) quadrature: np.gradient weights vs trapezoid, smooth bandpass, non-uniform grid ==")
rng = np.random.default_rng(0)
f_pl = lambda x: (x / 23.0) ** -3.1                                     # noqa: E731
f_d = lambda x: S.mbb_rj(x, 1.55, 19.6, 353.0)                           # noqa: E731
print(f"{'case':46s} {'npts':>5s} {'gradient':>12s} {'trapezoid':>12s}   (relative error of the band-averaged power law, beta=-3.1)")
for label, lo, hi in (("ends at 4 sigma (response ~3e-4 of peak)", 23 - 8.0, 23 + 8.0),
                      ("truncated at 5% of peak", 23 - 4.9, 23 + 4.9),
                      ("truncated at 30% of peak", 23 - 3.1, 23 + 3.1)):
    rfun = lambda x: np.exp(-0.5 * ((x - 23.0) / 2.0) ** 2) * (1 + 0.1 * np.sin(x))   # noqa: E731
    fine = np.linspace(lo, hi, 400001)
    exact = np.trapezoid(rfun(fine) * f_pl(fine), fine) / np.trapezoid(rfun(fine) * S.g_rj(fine), fine)
    for npts in (20, 60, 200):
        # non-uniform grid: jittered spacing with a 3:1 spread of step sizes
        steps = rng.uniform(0.5, 1.5, npts - 1)
        g = lo + (hi - lo) * np.concatenate([[0.0], np.cumsum(steps) / steps.sum()])
        eg = band_ratio(g, rfun(g), f_pl, "gradient") / exact - 1
        et = band_ratio(g, rfun(g), f_pl, "trapz") / exact - 1
        print(f"{label:46s} {npts:5d} {eg:+12.2e} {et:+12.2e}")
print("interior weights identical to trapezoid on a non-uniform grid:",
      np.allclose(np.gradient(g)[1:-1], (g[2:] - g[:-2]) / 2), "; end weights are", np.gradient(g)[0] / ((g[1] - g[0]) / 2), "x trapezoid")

print("\n== (d) real bandpasses ==")
bands = []
try:
    for b in ("K", "Ka", "Q", "V"):
        bands.append((f"WMAP {b}", S.wmap_band(b, "/Users/trygvels/data/wmap9_bandpass/extracted")))
except Exception as e:  # noqa: BLE001
    print("WMAP bandpasses not read:", type(e).__name__, e)
try:
    for f in (30, 44, 70):
        bands.append((f"LFI {f}", S.lfi_band(f)))
except Exception as e:  # noqa: BLE001
    print("LFI bandpasses not read:", type(e).__name__, e)
try:
    bands.append(("HFI 353 psb", S.hfi_band(353, psb=True)))
    bands.append(("HFI 353 full", S.hfi_band(353, psb=False)))
except Exception as e:  # noqa: BLE001
    print("HFI bandpasses not read:", type(e).__name__, e)
print(f"{'band':14s} {'npts':>5s} {'nu range':>18s} {'resp end/peak':>22s} {'grid':>10s} {'d s(-3.1)':>11s} {'d ds/s':>11s} {'d dust':>11s} {'s(-3.1) @23':>12s}")
for name, B in bands:
    nu = B.nu
    resp = B.w / np.gradient(nu)
    wt = np.empty_like(nu)
    wt[1:-1] = (nu[2:] - nu[:-2]) / 2
    wt[0], wt[-1] = (nu[1] - nu[0]) / 2, (nu[-1] - nu[-2]) / 2
    wt = wt * resp
    s, ds = B.power_law(-3.1, 23.0)
    lr = np.log(nu / 23.0)
    st = np.sum(wt * (nu / 23.0) ** -3.1) / np.sum(wt * S.g_rj(nu))
    dst = np.sum(wt * (nu / 23.0) ** -3.1 * lr) / np.sum(wt * S.g_rj(nu))
    dt = np.sum(wt * f_d(nu)) / np.sum(wt * S.g_rj(nu))
    steps = np.diff(nu)
    uni = "uniform" if np.allclose(steps, steps[0], rtol=1e-6) else f"{steps.min() / steps.max():.2f}"
    print(f"{name:14s} {nu.size:5d} {nu.min():8.2f}-{nu.max():8.2f} {resp[0] / resp.max():10.2e} {resp[-1] / resp.max():10.2e} {uni:>10s} "
          f"{s[0] / st - 1:+11.2e} {(ds[0] / s[0]) - (dst / st):+11.2e} {B.mbb(1.55, 19.6, 353.0) / dt - 1:+11.2e} {s[0]:12.5e}")
print("(d s, d dust: relative change code/trapezoid - 1; d ds/s: absolute change of the mean log frequency ratio)")

print("\n== (d1) docstring claim: sum(w f / g) / sum(w) against sum(w f) / sum(w g), WMAP K to V ==")
for name, B in bands[:4]:
    nu = B.nu
    out = []
    for be in (-3.1, -2.5):
        f = (nu / 23.0) ** be
        out.append((np.sum(B.w * f / S.g_rj(nu)) / np.sum(B.w)) / (np.sum(B.w * f) / np.sum(B.w * S.g_rj(nu))) - 1)
    print(f"{name:10s} relative difference of the two conventions: beta=-3.1 {out[0]:+.2e}, beta=-2.5 {out[1]:+.2e}")

print("\n== (d2) hfi_band: samples below the cut are dropped, then np.gradient spans the gaps ==")
try:
    from astropy.io import fits
    for psb in (True, False):
        with fits.open(S.HFI_RIMO) as h:
            t = h[f"BANDPASS_{'PSB' if psb else 'F'}353"].data
            nu0, tau0 = np.asarray(t["WAVENUMBER"], float) * S.C_CM_GHZ, np.asarray(t["TRANSMISSION"], float)
        sel = (nu0 > 0.5 * 353) & (nu0 < 2.0 * 353)
        nu1, tau1 = nu0[sel], tau0[sel]
        keep = tau1 > 1e-5 * tau0.max()
        resp_full = np.where(keep, tau1, 0.0) * nu1 ** 2            # same cut, but zeros kept in place
        B = S.hfi_band(353, psb=psb)
        steps = np.diff(B.nu)
        med = np.median(steps)
        norm = np.trapezoid(resp_full * S.g_rj(nu1), nu1)
        s_ref = np.trapezoid(resp_full * (nu1 / 23.0) ** -3.1, nu1) / norm
        d_ref = np.trapezoid(resp_full * f_d(nu1), nu1) / norm
        # in-band only (300 to 420 GHz) for scale
        inb = (nu1 > 300) & (nu1 < 420)
        d_inb = np.trapezoid((resp_full * f_d(nu1))[inb], nu1[inb]) / np.trapezoid((resp_full * S.g_rj(nu1))[inb], nu1[inb])
        print(f"{'PSB' if psb else 'all'}: native samples in 176 to 706 GHz {nu1.size}, kept {B.nu.size}; native step {np.median(np.diff(nu1)):.3f} GHz; "
              f"kept-grid steps above 3x median: {(steps > 3 * med).sum()}, largest {steps.max():.2f} GHz")
        print(f"     power law s(-3.1): hfi_band / zeros-in-place trapezoid - 1 = {B.power_law(-3.1, 23.0)[0][0] / s_ref - 1:+.2e};  "
              f"dust mbb: {B.mbb(1.55, 19.6, 353.0) / d_ref - 1:+.2e};  dust mbb, 300 to 420 GHz only vs full: {d_inb / d_ref - 1:+.2e}")
except Exception as e:  # noqa: BLE001
    print("skipped:", type(e).__name__, e)

print("\n== (e) effective_frequency ==")
for name, B in bands[:1] + bands[4:5] + bands[-1:]:
    for be in (-3.1, -2.5):
        ne = B.effective_frequency(be)
        s, _ = B.power_law(be, 1.0)
        print(f"{name:14s} beta {be}: nu_eff {ne:.4f}; delta-band response at nu_eff / band average - 1 = {ne ** be / S.g_rj(ne) / s[0] - 1:+.1e}; "
              f"grid step {(B.nu.max() - B.nu.min()) / 4000:.4f} GHz")
print("delta band:", S.Band.delta("x", 30.0).effective_frequency(-3.1))
# non-monotonic case: nu^beta / g(nu) has a minimum where d ln g / d ln nu = beta (near 287 GHz for beta = -3.1);
# a band across it has two roots and argmin returns whichever grid point is closest
nu = np.linspace(220.0, 420.0, 400)
B = S.Band("wide", nu, np.ones_like(nu))
grid = np.linspace(220, 420, 4001)
y = grid ** -3.1 / S.g_rj(grid)
print(f"flat band 220 to 420 GHz, beta=-3.1: nu^beta/g has its minimum at {grid[np.argmin(y)]:.1f} GHz (inside the band, not monotonic); "
      f"effective_frequency returns {B.effective_frequency(-3.1):.2f}; number of roots of (y - s): "
      f"{int(np.sum(np.diff(np.sign(y - B.power_law(-3.1, 1.0)[0][0])) != 0))}")
for name, BB in bands:
    g2 = np.linspace(BB.nu.min(), BB.nu.max(), 2001)
    y2 = g2 ** -3.1 / S.g_rj(g2)
    print(f"   {name:14s} nu^-3.1/g monotonic over the band: {bool(np.all(np.diff(y2) < 0) or np.all(np.diff(y2) > 0))}")

print("\n== (f) wmap_band detector weighting ==")
try:
    dln, lnu = [], []
    for band in ("K", "Ka", "Q", "V"):
        nus, r_unit, r_raw, areas_all, same = None, None, None, [], True
        for rad in S.WMAP_RADIOMETERS[band]:
            a = np.loadtxt(os.path.join("/Users/trygvels/data/wmap9_bandpass/extracted", f"wmap_bandpass_{rad}_v5.cbp"))
            if nus is None:
                nus, r_unit, r_raw = a[:, 0], np.zeros(a.shape[0]), np.zeros(a.shape[0])
            areas = [np.trapezoid(a[:, 1], a[:, 0]), np.trapezoid(a[:, 2], a[:, 0])]
            areas_all += areas
            same &= (a.shape[0] == nus.size and bool(np.allclose(a[:, 0], nus)))
            r_unit += np.interp(nus, a[:, 0], a[:, 1] / areas[0], left=0, right=0) + np.interp(nus, a[:, 0], a[:, 2] / areas[1], left=0, right=0)
            r_raw += np.interp(nus, a[:, 0], a[:, 1], left=0, right=0) + np.interp(nus, a[:, 0], a[:, 2], left=0, right=0)
        su = S.Band("u", nus, np.clip(r_unit, 0, None)).power_law(-3.1, 23.0)[0][0]
        sr = S.Band("r", nus, np.clip(r_raw, 0, None)).power_law(-3.1, 23.0)[0][0]
        code = S.wmap_band(band, "/Users/trygvels/data/wmap9_bandpass/extracted").power_law(-3.1, 23.0)[0][0]
        dln.append(np.log(su / sr))
        lnu.append(np.log(S.WMAP_NOMINAL[band]))
        print(f"WMAP {band:2s}: detector areas min {min(areas_all):.3f} max {max(areas_all):.3f}; all radiometers on one grid: {same}; "
              f"s(-3.1): wmap_band {code:.6e}, unit area per detector {su:.6e}, raw sum (area weighted) {sr:.6e}; unit/raw - 1 = {su / sr - 1:+.2e}")
    slope = np.polyfit(lnu, dln, 1)[0]
    print(f"change of ln s between the two weightings, regressed on ln nu over K, Ka, Q, V: slope {slope:+.4f} "
          "(the size of the shift in a fitted index if the other weighting is the right one)")
except Exception as e:  # noqa: BLE001
    print("skipped:", type(e).__name__, e)
