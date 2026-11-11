"""Referee check (draft 3): how much of the Plane's excess chi-square sits in the pixel that contains Tau A and its
neighbour. The band covariances are not loaded; the noise per element is taken as the rms residual of the Plane
elements in the band divided by the square root of the printed chi-square per element of that band (a diagonal
proxy). Read only."""
import numpy as np
c = np.load("experiments/E040/s4_bpxv_n16_wm_lcdm.npz", allow_pickle=True)
res, reg = c["resid"], c["region"]
n = reg.size // 2
chi_band = np.array([0.141, 1.007, 1.057, 0.997])        # printed chi2 per band and element (whole sky), K, Ka, Q, V
pl = np.flatnonzero(reg == 3)
rms = np.sqrt((res[:4][:, pl] ** 2).mean(1))
sig = rms / np.sqrt(chi_band)                             # crude noise per element
print("Plane residual rms by band", rms.round(2), " noise proxy", sig.round(2))
tot_unw = (res[:4][:, pl] ** 2).sum()
for name, p in (("Tau A pixel 1591", 1591), ("neighbour 1597", 1597)):
    el = np.array([p, p + n])
    unw = (res[:4][:, el] ** 2).sum()
    chi = ((res[:4][:, el] / sig[:, None]) ** 2).sum()
    print(f"{name}: share of the Plane's unweighted residual sum of squares {unw / tot_unw:.3f}; chi-square under the proxy {chi:.0f}")
nel = pl.size
print(f"Plane: {nel} elements; printed chi2 per element and band 0.770 against 0.6: total {0.770 * nel * 5:.0f}, expected {0.6 * nel * 5:.0f}, excess {0.170 * nel * 5:.0f}")
# the same for the brightest pixels of the other regions: share of each region's unweighted residual power in its 5 brightest pixels
a = c["a_s"]; P = np.hypot(a[:n], a[n:])
for r, nm in enumerate(["High Latitude", "Spur", "Galactic Center", "Plane"]):
    s = np.flatnonzero(reg[:n] == r)
    top = s[np.argsort(-P[s])[:5]]
    el = np.r_[top, top + n]; allel = np.r_[s, s + n]
    print(f"{nm:16s}: 5 brightest of {s.size} pixels hold {(res[:4][:, el] ** 2).sum() / (res[:4][:, allel] ** 2).sum():.2f} of the unweighted residual sum of squares")
