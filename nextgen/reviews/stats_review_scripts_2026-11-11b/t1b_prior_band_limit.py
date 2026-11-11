"""Addition 1, claim (a), second part: what the band limit lmax = 3 nside - 1 and the 0.01 uK floor leave out.

(i) Exact covariance of a pixel-averaged field at nside 4 (average of 64 sub-pixel values at nside 32), built mode by
    mode to l = 64, against the model of `lcdm_cmb_precision` (pointwise harmonics times the pixel window, l <= 11,
    plus the floor). Gives the effective window w_eff^2(l) = (variance of the pixel average) / (pointwise variance),
    exact for l <= 64 at nside 4, i.e. to l / nside = 16.
(ii) Variance per element of the LCDM polarization at nside 16 (60 arcmin window) and nside 8, in band and above the
    band limit, as one-dimensional sums over l with the tabulated HEALPix window where it exists (l <= 4 nside) and
    w_eff(l / nside) from (i) above it (HEALPix pixels are self-similar; the scaling is checked on the tabulated
    range).
No map above nside 4 is used as a working grid; nside 32 is only the quadrature grid inside the nside-4 pixels.
"""
import sys
import time

import camb
import healpy as hp
import numpy as np

sys.path.insert(0, "nextgen")
from lowres.window import WindowMatcher, lcdm_cmb_precision  # noqa: E402

t0 = time.time()
NLO, NHI, LTOP = 4, 32, 64
npl, nph = hp.nside2npix(NLO), hp.nside2npix(NHI)
k = nph // npl
size = hp.Alm.getsize(LTOP)
zero = np.zeros(size, dtype=complex)
omega = 4 * np.pi / npl
PY, ells, kinds = [], [], []
for which in (1, 2):
    for l in range(2, LTOP + 1):
        for m in range(0, l + 1):
            for part in ((1.0,) if m == 0 else (1.0, 1.0j)):
                alm = np.zeros(size, dtype=complex)
                alm[hp.Alm.getidx(LTOP, l, m)] = part
                _, q, u = hp.alm2map([zero, alm, zero] if which == 1 else [zero, zero, alm], NHI, lmax=LTOP, pol=True)
                q = hp.reorder(q, r2n=True).reshape(npl, k).mean(1)             # pixel average, NESTED
                u = hp.reorder(u, r2n=True).reshape(npl, k).mean(1)
                PY.append(np.concatenate([q, u]) / (1.0 if m == 0 else np.sqrt(2.0)))   # unit: field per unit real coefficient
                ells.append(l)
                kinds.append(which)
PY, ells, kinds = np.array(PY).T, np.array(ells), np.array(kinds)
print(f"# pixel-averaged modes built: {PY.shape}, {time.time() - t0:.0f} s")
# effective window: mean over elements of sum_modes(l) (P y)^2 divided by the pointwise value (2l+1)/(8 pi) per E or B
weff2 = np.zeros(LTOP + 1)
for l in range(2, LTOP + 1):
    sel = ells == l
    weff2[l] = np.mean(np.sum(PY[:, sel] ** 2, axis=1)) / (2 * (2 * l + 1) / (8 * np.pi))
pw4 = hp.pixwin(NLO, pol=True, lmax=4 * NLO)[1]
print("# (i) effective window of the pixel average at nside 4 against the tabulated polarization window squared")
print("   l        : " + " ".join(f"{l:6d}" for l in (2, 4, 8, 11, 12, 16, 24, 32, 48, 64)))
print("   exact    : " + " ".join(f"{weff2[l]:6.4f}" for l in (2, 4, 8, 11, 12, 16, 24, 32, 48, 64)))
print("   tabulated: " + " ".join(f"{pw4[l] ** 2:6.4f}" if l <= 4 * NLO else "     -" for l in (2, 4, 8, 11, 12, 16, 24, 32, 48, 64)))

pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=2000, lens_potential_accuracy=1)
cls = camb.get_results(pars).get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]
pars0 = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=300)
cls0 = camb.get_results(pars0).get_cmb_power_spectra(pars0, CMB_unit="muK", raw_cl=True)["total"]
print(f"# CAMB as called by the code (lmax 300) against lmax 2000 with accurate lensing: EE ratio at l = 2, 10, 47: "
      + " ".join(f"{cls0[l, 1] / cls[l, 1]:.4f}" for l in (2, 10, 47)) + "; BB ratio: " + " ".join(f"{cls0[l, 2] / cls[l, 2]:.3f}" for l in (2, 10, 47)))


def weff_scaled(l, nside):
    """w_eff^2 at multipole l for `nside`, from the nside-4 curve at the same l / nside."""
    return np.interp(l * NLO / nside, np.arange(2, LTOP + 1), weff2[2:])


print("# check of the self-similar scaling on the tabulated range: tabulated w^2(nside 16, l) against w_eff^2(nside 4, l / 4)")
pw16 = hp.pixwin(16, pol=True, lmax=64)[1]
print("   l = 16, 32, 47, 64: tabulated " + " ".join(f"{pw16[l] ** 2:.4f}" for l in (16, 32, 47, 64)) + "; scaled "
      + " ".join(f"{weff_scaled(l, 16):.4f}" for l in (16, 32, 47, 64)))

print("# (ii) variance per element of pixel-averaged LCDM Q or U [uK^2]: sum_l (2l+1)/(8 pi) (C_EE + C_BB) b_l^2 w_l^2")
for nside, target in ((16, 60.0), (8, 60.0), (4, 60.0)):
    lmax = 3 * nside - 1
    ltop = min(16 * nside, 2000)
    l = np.arange(ltop + 1)
    bl = hp.gauss_beam(np.radians(target / 60.0), lmax=ltop, pol=True)[:, 1]
    pwt = hp.pixwin(nside, pol=True, lmax=4 * nside)[1]
    w2 = np.where(l <= 4 * nside, np.concatenate([pwt, np.zeros(ltop - 4 * nside)]) ** 2, weff_scaled(l, nside))
    w2s = weff_scaled(np.maximum(l, 2), nside)
    term = (2 * l + 1) / (8 * np.pi) * (cls[:ltop + 1, 1] + cls[:ltop + 1, 2]) * bl ** 2
    v_in = np.sum((term * w2)[2:lmax + 1])
    v_mid = np.sum((term * w2)[lmax + 1:4 * nside + 1])
    v_hi = np.sum((term * w2)[4 * nside + 1:])
    v_out_s = np.sum((term * w2s)[lmax + 1:])
    print(f"   nside {nside:2d}, window {target:.0f} arcmin: l <= {lmax}: {v_in:.5f} (rms {np.sqrt(v_in):.3f} uK); "
          f"{lmax + 1} <= l <= {4 * nside} (tabulated window): {v_mid:.5f}; {4 * nside} < l <= {ltop} (scaled window): {v_hi:.5f}; "
          f"above the band limit in total: {v_mid + v_hi:.5f} (rms {np.sqrt(v_mid + v_hi):.3f} uK; with the exact-average window throughout {v_out_s:.5f}) "
          f"= {(v_mid + v_hi) / v_in:.2f} of the in-band variance and {(v_mid + v_hi) / 1e-4:.0f} times the floor variance 1e-4")

# (iii) exact structure at nside 4: where does the power above the band limit go, relative to span(Y)?
WM = WindowMatcher(NLO, target_arcmin=60.0)
Y, Yp = WM.Y, WM.Yp
Pspan = Y @ Yp
n = Y.shape[0]
hi = ells > WM.lmax
print(f"# (iii) nside 4, unit white spectrum (C_l = 1) above the band limit, exact pixel average: share of that variance inside span(Y) "
      f"({Y.shape[1]} of {n} dimensions) and outside")
for lo_, hi_ in ((12, 16), (17, 32), (33, 64)):
    sel = (ells >= lo_) & (ells <= hi_)
    C = PY[:, sel] @ PY[:, sel].T
    tin, tout = np.trace(Pspan @ C @ Pspan.T), np.trace((np.eye(n) - Pspan) @ C @ (np.eye(n) - Pspan).T)
    low = WM.ell <= 4
    tl = np.trace(Yp[low] @ C @ Yp[low].T)
    print(f"   l {lo_:2d} to {hi_:2d}: inside {tin / np.trace(C):.3f}, outside {tout / np.trace(C):.3f}; "
          f"variance per coefficient at l <= 4 (leak into the largest scales) / variance per element: {tl / low.sum() / (np.trace(C) / n):.4f}")
# model check for l <= lmax: exact pixel-average covariance against Y diag(w^2) Y^T / Omega for a white unit spectrum
sel = ells <= WM.lmax
C_exact = PY[:, sel] @ PY[:, sel].T
C_model = (Y * (pw4[WM.ell] ** 2 / omega)[None, :]) @ Y.T
ev = np.linalg.eigvalsh(np.linalg.solve(np.linalg.cholesky(C_model + 1e-6 * np.eye(n)), np.linalg.solve(np.linalg.cholesky(C_model + 1e-6 * np.eye(n)), C_exact).T))
print(f"# in band (l <= 11, unit spectrum), exact pixel average against the pixel-window model: max |difference| / max |exact| "
      f"{np.abs(C_exact - C_model).max() / np.abs(C_exact).max():.3f}; diagonal ratio model/exact {np.mean(np.diag(C_model)) / np.mean(np.diag(C_exact)):.4f} "
      f"(element to element {np.min(np.diag(C_model) / np.diag(C_exact)):.3f} to {np.max(np.diag(C_model) / np.diag(C_exact)):.3f})")
np.savez("nextgen/reviews/stats_review_scripts_2026-11-11b/pixavg_modes_n4.npz", PY=PY, ells=ells, kinds=kinds, weff2=weff2)
print(f"# done, {time.time() - t0:.0f} s")
