"""Check 5: lfi_noise.py against brute-force dense algebra.

(a) NoiseModel.neg2loglike vs z^T R^-1 z + ln det R with dense R, random small basis Y (not orthonormal),
    more and fewer modes than data points, extreme r/A ratios.
(b) the same with the real spin-2 basis from spin2_basis at nside 2 and 4; orthonormality of Y.
(c) sqrt_blocks vs scipy.linalg.sqrtm on random 2x2 SPD blocks, including nearly singular ones.
(d) lfi_ninv and whitened_half_difference vs dense matrix square roots (the file loader is replaced by an
    in-memory stand-in inside this script only; nothing under review is edited).
(e) score of neg2loglike has zero mean at the truth (Monte Carlo), i.e. quad and logdet are mutually consistent.
"""
import os
import sys

import numpy as np
import scipy.linalg as sl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/Users/trygvels/code/Commander4/nextgen")
import lowres.lfi_noise as LN  # noqa: E402

np.set_printoptions(precision=5, linewidth=160, suppress=True)
rng = np.random.default_rng(3)


def fake_model(n, k, ortho=False):
    m = LN.NoiseModel.__new__(LN.NoiseModel)
    Y = rng.standard_normal((n, k)) / np.sqrt(n)
    if ortho:
        Y = np.linalg.qr(Y)[0]
    m.Y, m.ell = Y, rng.integers(2, 48, k).astype(float)
    m.G = Y.T @ Y
    m.nside, m.lmax = None, None
    return m


def brute(m, theta, z):
    R = m.R(theta)
    sign, ld = np.linalg.slogdet(R)
    assert sign > 0
    return z @ np.linalg.solve(R, z) + ld


print("== (a) random basis ==")
print(f"{'n':>4s} {'k':>4s} {'theta':>22s} {'code':>18s} {'brute':>18s} {'abs diff':>10s}")
for n, k in ((60, 20), (60, 59), (40, 40)):
    m = fake_model(n, k)
    z = rng.standard_normal(n) * 1.3
    for theta in ((0.2, 0.0, 1.0), (-1.0, 2.0, 2.5), (1.5, -3.0, 0.3), (-6.0, 6.0, 3.0), (3.0, -12.0, 0.0)):
        c, b = m.neg2loglike(np.array(theta), z), brute(m, np.array(theta), z)
        print(f"{n:4d} {k:4d} {str(theta):>22s} {c:18.9f} {b:18.9f} {abs(c - b):10.1e}")
# R model consistency: R(theta) itself against an elementwise construction
m = fake_model(30, 7)
th = np.array([0.3, 0.7, 1.4])
R2 = np.exp(th[0]) * np.eye(30)
for j in range(7):
    R2 += np.exp(th[1]) * (m.ell[j] / 10.0) ** (-th[2]) * np.outer(m.Y[:, j], m.Y[:, j])
print("R(theta) vs elementwise sum: max abs diff", np.max(np.abs(m.R(th) - R2)))

print("\n== (b) real spin-2 basis ==")
for nside in (2, 4):
    m = LN.NoiseModel(nside=nside)
    n, k = m.Y.shape
    off = m.G - np.eye(k)
    w = np.linalg.eigvalsh(m.G)
    print(f"nside {nside}: Y {m.Y.shape}, lmax {m.lmax}; max|G - I| {np.max(np.abs(off)):.3f}; eig(G) min {w.min():.2e} max {w.max():.3f}; "
          f"modes with eig < 0.5: {(w < 0.5).sum()}")
    z = rng.standard_normal(n)
    for theta in ((0.2, 0.0, 1.0), (-0.5, 1.0, 2.0)):
        c, b = m.neg2loglike(np.array(theta), z), brute(m, np.array(theta), z)
        print(f"   theta {theta}: code {c:.9f} brute {b:.9f} diff {abs(c - b):.1e}")
    # per-l orthonormality restricted to l <= 2 nside (where the pixelization supports the modes)
    lowl = m.ell <= 2 * nside
    print(f"   max|G - I| restricted to l <= 2 nside: {np.max(np.abs(off[np.ix_(lowl, lowl)])):.3f}")

print("\n== (c) sqrt_blocks ==")
k = 2000
X = rng.standard_normal((k, 2, 2))
S = X @ np.transpose(X, (0, 2, 1)) + 1e-3 * np.eye(2)
S[:5] = [[[1.0, 0.0], [0.0, 4.0]], [[2.0, 1.999999], [1.999999, 2.0]], [[1e-8, 0.0], [0.0, 1e8]], [[3.0, -1.0], [-1.0, 0.5]],
         [[1.0, 0.0], [0.0, 1.0]]]
(sa, sb, sc), (ia, ib, ic) = LN.sqrt_blocks(S[:, 0, 0], S[:, 1, 1], S[:, 0, 1])
e1 = e2 = e3 = 0.0
for i in range(k):
    Rt = sl.sqrtm(S[i]).real
    Sq = np.array([[sa[i], sc[i]], [sc[i], sb[i]]])
    Iq = np.array([[ia[i], ic[i]], [ic[i], ib[i]]])
    scale = np.max(np.abs(Rt))
    e1 = max(e1, np.max(np.abs(Sq - Rt)) / scale)
    e2 = max(e2, np.max(np.abs(Sq @ Sq - S[i])) / np.max(np.abs(S[i])))
    e3 = max(e3, np.max(np.abs(Iq @ Sq - np.eye(2))))
print(f"max rel |sqrt - sqrtm| {e1:.1e};  max rel |sqrt^2 - S| {e2:.1e};  max |isqrt sqrt - I| {e3:.1e}  ({k} blocks, incl. edge cases)")

print("\n== (d) lfi_ninv and whitened_half_difference with an in-memory loader ==")
npix = 24
n2 = 2 * npix


def block_ninv(scale):
    Xb = rng.standard_normal((npix, 2, 2))
    Sb = scale * (Xb @ np.transpose(Xb, (0, 2, 1)) + 0.2 * np.eye(2))      # covariance blocks
    W = np.zeros((n2, n2))
    i = np.arange(npix)
    W[i, i], W[npix + i, npix + i], W[i, npix + i], W[npix + i, i] = Sb[:, 0, 0], Sb[:, 1, 1], Sb[:, 0, 1], Sb[:, 0, 1]
    return W, np.linalg.inv(W)


store = {}
for split in ("full", "full-ringhalf-1", "full-ringhalf-2"):
    W, Wi = block_ninv(1.0 if split == "full" else 2.0)
    store[split] = (rng.standard_normal(n2), Wi, W)
orig = LN.D.load_planck_qu
LN.D.load_planck_qu = lambda path, nside=16, scale=1e6: store[[s for s in store if path.endswith(s + ".fits")][0]][:2]
m = fake_model(n2, 15)
theta = np.array([0.1, 0.8, 1.7])
d, ninv = LN.lfi_ninv(30, theta, m)
Wfull = store["full"][2]
Wh = sl.sqrtm(Wfull).real
Ndense = Wh @ m.R(theta) @ Wh
ref = np.linalg.inv(Ndense)
print("lfi_ninv vs inv(W^1/2 R W^1/2): max abs diff / max abs", np.max(np.abs(ninv - ref)) / np.max(np.abs(ref)))
print("lfi_ninv returns the full-mission map:", np.array_equal(d, store["full"][0]))
z = LN.whitened_half_difference(30)
Wavg = (store["full-ringhalf-1"][2] + store["full-ringhalf-2"][2]) / 4
zref = np.linalg.inv(sl.sqrtm(Wavg).real) @ (store["full-ringhalf-1"][0] - store["full-ringhalf-2"][0]) / 2
print("whitened_half_difference vs dense: max abs diff", np.max(np.abs(z - zref)))
LN.D.load_planck_qu = orig

print("\n== (e) Monte Carlo: mean gradient of neg2loglike at the truth (should be zero) ==")
m = fake_model(80, 25)
theta = np.array([0.2, 0.5, 1.2])
Lc = np.linalg.cholesky(m.R(theta))
g = []
for _ in range(4000):
    z = Lc @ rng.standard_normal(80)
    gi = np.zeros(3)
    for j in range(3):
        e = np.zeros(3); e[j] = 1e-4
        gi[j] = (m.neg2loglike(theta + e, z) - m.neg2loglike(theta - e, z)) / 2e-4
    g.append(gi)
g = np.array(g)
print("mean gradient", g.mean(0), "+-", g.std(0, ddof=1) / np.sqrt(len(g)), " z =", g.mean(0) / (g.std(0, ddof=1) / np.sqrt(len(g))))
# negative control: drop the log-determinant, the gradient must then be biased
print("negative control (quad only, by subtracting ln det R): mean gradient in ln r",
      np.mean([(lambda z: ((m.neg2loglike(theta + [1e-4, 0, 0], z) - np.linalg.slogdet(m.R(theta + [1e-4, 0, 0]))[1])
                           - (m.neg2loglike(theta - [1e-4, 0, 0], z) - np.linalg.slogdet(m.R(theta - [1e-4, 0, 0]))[1])) / 2e-4)(
                               Lc @ rng.standard_normal(80)) for _ in range(500)]))
