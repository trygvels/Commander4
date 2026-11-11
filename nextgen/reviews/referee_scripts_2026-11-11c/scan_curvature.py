"""Referee check of claim (a): curvature and slope of the flat-prior marginal along the High Latitude index, from the
saved scans (E044/scan_flat.npz, scan_lcdm.npz). Objective = -ln L; H = -d2 ln L / dbeta2. Reads only."""
import numpy as np
np.set_printoptions(precision=2, suppress=True, linewidth=170)
for k in ("flat", "lcdm"):
    z = np.load(f"nextgen/experiments/E044/scan_{k}.npz", allow_pickle=True)
    fit = np.load(f"nextgen/experiments/E044/{'flat' if k == 'flat' else 'lcdm_all'}.npz", allow_pickle=True)
    print(f"=== CMB prior {k}: profile optimum {z['profile_opt']}, flat23 optimum {z['flat23']}, flat30 {z['flat30']}")
    print("    stored Hmarg diag (evaluated in e044_fit.py at the flat23 optimum):", np.diag(fit["Hmarg"]), " J diag", np.diag(fit["J"]), " Jc diag", np.diag(fit["J"] - fit["T"]))
    for sl in (0, 1):
        r = z["rows"][z["rows"][:, 0] == sl]
        x, lp, l23, l30 = r[:, 1], r[:, 2], r[:, 3], r[:, 4]
        h = x[1] - x[0]
        print(f"  slice {sl} ({'others at profile optimum' if sl == 0 else 'others at flat23 optimum'}): n = {x.size}, step {h:.2f}")
        print("    monotone increasing 23:", bool(np.all(np.diff(l23) > 0)), " 30:", bool(np.all(np.diff(l30) > 0)), " argmax 23:", x[np.argmax(l23)], " 30:", x[np.argmax(l30)])
        xm = x[1:-1]
        Hp = -(lp[2:] - 2 * lp[1:-1] + lp[:-2]) / h ** 2
        H23 = -(l23[2:] - 2 * l23[1:-1] + l23[:-2]) / h ** 2
        H30 = -(l30[2:] - 2 * l30[1:-1] + l30[:-2]) / h ** 2
        s23 = (l23[2:] - l23[:-2]) / (2 * h)
        s30 = (l30[2:] - l30[:-2]) / (2 * h)
        print("    x          ", xm)
        print("    H profile  ", Hp)
        print("    H flat23   ", H23)
        print("    H flat30   ", H30, "(identical to flat23: the tilt is linear)")
        print("    volume curv", H23 - Hp)
        print("    slope lnL23", s23)
        print("    slope lnL30", s30)
        j = np.argmin(np.abs(xm - z["profile_opt"][0]))
        print(f"    nearest grid point to the profile optimum: {xm[j]:.2f}: H_flat23 = {H23[j]:+.0f}, H_profile = {Hp[j]:+.0f}, slope of lnL_flat23 = {s23[j]:+.0f}, "
              f"first-order distance to an optimum slope/H = {s23[j] / H23[j]:+.2f} (range ends {(-1.6 - xm[j]):.2f} above)")
        # one-sided estimate at the upper edge by a quadratic through the last three points
        c = np.polyfit(x[-4:], l23[-4:], 3)
        print(f"    cubic through the last 4 points: H at -1.6 = {-np.polyval(np.polyder(c, 2), -1.6):+.0f}, slope at -1.6 = {np.polyval(np.polyder(c, 1), -1.6):+.0f}")
        # quartic local fit around profile optimum
        m = np.abs(x - z["profile_opt"][0]) < 0.5
        c = np.polyfit(x[m] - z["profile_opt"][0], l23[m], 4)
        print(f"    quartic fit within 0.5 of the profile optimum: H_flat23 at the profile optimum = {-2 * c[-3]:+.0f}, slope = {c[-2]:+.0f}")
