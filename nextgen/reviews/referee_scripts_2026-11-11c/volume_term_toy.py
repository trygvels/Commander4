"""Referee: why the flat-prior marginal with a free CMB amplitude has no optimum below -1.6. One element, delta bands at
the effective frequencies of K, Ka, Q, V, equal noise in thermodynamic units, synchrotron (reference 23 GHz) and CMB
amplitudes free (dust is pinned by 353 GHz and left out). Volume term v(beta) = -(1/2) ln det A, A = [[s.s, s.1],[s.1, 1.1]].
Compare its slope per amplitude at the profile optimum with the slope measured in the scan of the real fit
(+396 per unit index for 2974 amplitudes = 0.133 per amplitude; curvature -175 / 2974 = -0.059)."""
import numpy as np
nu = np.array([22.41, 32.66, 40.30, 59.99])
x = nu / 56.78
g = x ** 2 * np.exp(x) / np.expm1(x) ** 2                       # dT_RJ / dT_thermo
def v(beta, ref=23.0, w=np.ones(4)):
    s = (nu / ref) ** beta / g
    A = np.array([[np.sum(w * s * s), np.sum(w * s)], [np.sum(w * s), np.sum(w)]])
    return -0.5 * np.log(np.linalg.det(A))
for w, lab in ((np.ones(4), "equal weights"), (1 / np.array([2.0, 2.0, 2.0, 2.2]) ** 2 * g ** 2, "noise 2.0, 2.0, 2.0, 2.2 uK_RJ")):
    b = np.arange(-4.6, 0.51, 0.3)
    vv = np.array([v(t, w=w) for t in b])
    h = 1e-3
    for b0 in (-3.04,):
        s1 = (v(b0 + h, w=w) - v(b0 - h, w=w)) / (2 * h)
        s2 = (v(b0 + h, w=w) - 2 * v(b0, w=w) + v(b0 - h, w=w)) / h ** 2
        print(f"{lab}: slope of the volume term per amplitude at {b0}: {s1:+.3f}; second derivative {s2:+.3f} (curvature of the objective {-s2:+.3f}); implied pivot {23 * np.exp(-s1):.1f} GHz")
    j = int(np.argmax(vv))
    fine = np.arange(-1.0, 0.6, 0.01); vf = np.array([v(t, w=w) for t in fine])
    print("   v(beta) - v(-3.04):", " ".join(f"{t:+.1f}:{y - v(-3.04, w=w):+.2f}" for t, y in zip(b, vv)))
    print(f"   maximum of the volume term at beta = {fine[np.argmax(vf)]:+.2f} (synchrotron nearly degenerate with the CMB), value {vf.max() - v(-3.04, w=w):+.2f} per amplitude above its value at -3.04")
