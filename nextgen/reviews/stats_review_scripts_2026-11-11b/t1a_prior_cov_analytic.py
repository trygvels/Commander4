"""Addition 1, claim (a): the Q/U covariance from `spin2_basis` and `lcdm_cmb_precision` against an analytic
two-point function that uses no harmonic transform.

Reference: for two directions separated by theta, with Q_r, U_r defined relative to the great circle joining them,
    <Q_r Q_r> = (xi+ + xi-) / 2,   <U_r U_r> = (xi+ - xi-) / 2,   <Q_r U_r> = 0,
    xi+ = sum_l (2l+1)/(4 pi) (C_E + C_B) d^l_{2,2}(theta),  xi- = sum_l (2l+1)/(4 pi) (C_E - C_B) d^l_{2,-2}(theta),
    d^l_{2,2} = cos^4(theta/2) P_{l-2}^{(0,4)}(cos theta),  d^l_{2,-2} = sin^4(theta/2) P_{l-2}^{(4,0)}(cos theta).
Local Q, U (HEALPix COSMO convention: angle from e_theta, south, toward e_phi, east):
    Q + iU = exp(2 i alpha) (Q_r + i U_r),  alpha = angle of the great-circle tangent from e_theta toward e_phi.
Pixel positions from hp.pix2vec(nest=True); no alm2map, no reorder.
"""
import sys

import healpy as hp
import numpy as np
from scipy.special import eval_jacobi

sys.path.insert(0, "nextgen")
from lowres.lfi_noise import spin2_basis  # noqa: E402
from lowres.window import WindowMatcher, lcdm_cmb_precision  # noqa: E402


def analytic_cov(nside, clE, clB, handed=+1.0):
    """(2 npix, 2 npix) covariance, NESTED, Q then U. clE, clB indexed by l (pointwise field)."""
    npix = hp.nside2npix(nside)
    ip = np.arange(npix)
    th, ph = hp.pix2ang(nside, ip, nest=True)
    r = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)]).T
    eth = np.array([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)]).T
    eph = np.array([-np.sin(ph), np.cos(ph), np.zeros(npix)]).T
    cosb = np.clip(r @ r.T, -1.0, 1.0)
    sinb = np.sqrt(np.clip(1.0 - cosb ** 2, 0.0, None))
    # tangent at i toward j
    t_i = r[None, :, :] - cosb[:, :, None] * r[:, None, :]
    small = sinb < 1e-6                 # identical or antipodal directions (r.r rounds to 1 - 1e-16, sin to 1e-8)
    nrm = np.where(small, 1.0, sinb)
    t_i = t_i / nrm[:, :, None]
    t_i[small] = np.broadcast_to(eth[:, None, :], t_i.shape)[small]            # any tangent when the circle is undefined
    # the same great circle continued to j: t_j = -sin(theta) r_i + cos(theta) t_i
    t_j = -sinb[:, :, None] * r[:, None, :] + cosb[:, :, None] * t_i
    a_i = np.arctan2(np.einsum("ijk,ik->ij", t_i, eph), np.einsum("ijk,ik->ij", t_i, eth)) * handed
    a_j = np.arctan2(np.einsum("ijk,jk->ij", t_j, eph), np.einsum("ijk,jk->ij", t_j, eth)) * handed
    xp, xm = np.zeros_like(cosb), np.zeros_like(cosb)
    ch4, sh4 = ((1 + cosb) / 2) ** 2, ((1 - cosb) / 2) ** 2
    for l in range(2, len(clE)):
        f = (2 * l + 1) / (4 * np.pi)
        xp += f * (clE[l] + clB[l]) * ch4 * eval_jacobi(l - 2, 0, 4, cosb)
        xm += f * (clE[l] - clB[l]) * sh4 * eval_jacobi(l - 2, 4, 0, cosb)
    qq, uu = 0.5 * (xp + xm), 0.5 * (xp - xm)
    c1, s1, c2, s2 = np.cos(2 * a_i), np.sin(2 * a_i), np.cos(2 * a_j), np.sin(2 * a_j)
    QQ = c1 * c2 * qq + s1 * s2 * uu
    QU = c1 * s2 * qq - s1 * c2 * uu
    UQ = s1 * c2 * qq - c1 * s2 * uu
    UU = s1 * s2 * qq + c1 * c2 * uu
    return np.block([[QQ, QU], [UQ, UU]])


def code_cov(nside, lmax, clE, clB):
    Y, ell = spin2_basis(nside, lmax)
    nE = Y.shape[1] // 2
    omega = 4 * np.pi / hp.nside2npix(nside)
    cl = np.concatenate([clE[ell[:nE]], clB[ell[nE:]]])
    return (Y * (cl / omega)[None, :]) @ Y.T, Y, ell


def rel(a, b):
    return np.abs(a - b).max() / np.abs(b).max()


print("# t1a: basis covariance Y diag(C_l / Omega) Y^T against the analytic two-point function")
rng = np.random.default_rng(1)
for nside in (2, 4):
    lmax = 3 * nside - 1
    clE = np.zeros(lmax + 1)
    clB = np.zeros(lmax + 1)
    clE[2:] = rng.uniform(0.5, 2.0, lmax - 1)
    clB[2:] = rng.uniform(0.05, 0.4, lmax - 1)
    Cc, Y, ell = code_cov(nside, lmax, clE, clB)
    Ca = analytic_cov(nside, clE, clB)
    n = Cc.shape[0] // 2
    print(f"nside {nside}, lmax {lmax}, random C_l^EE in [0.5, 2], C_l^BB in [0.05, 0.4]:")
    print(f"   max |code - analytic| / max |analytic|: all {rel(Cc, Ca):.2e}; QQ {rel(Cc[:n, :n], Ca[:n, :n]):.2e}; "
          f"QU {np.abs(Cc[:n, n:] - Ca[:n, n:]).max() / np.abs(Ca).max():.2e}; UU {rel(Cc[n:, n:], Ca[n:, n:]):.2e}")
    print(f"   diagonal: code min {np.diag(Cc).min():.6f} max {np.diag(Cc).max():.6f}; analytic sum (2l+1)/(8 pi)(C_E + C_B) = "
          f"{np.sum((2 * np.arange(lmax + 1) + 1) / (8 * np.pi) * (clE + clB)):.6f}")
    # negative controls: what a wrong convention would give
    print(f"   controls (same statistic): E and B swapped {rel(Cc, analytic_cov(nside, clB, clE)):.2e}; "
          f"U sign flipped (IAU) {rel(Cc, analytic_cov(nside, clE, clB, handed=-1.0)):.2e}; "
          f"RING order instead of NESTED {rel(Cc, (lambda P: Ca[np.ix_(P, P)])(np.concatenate([hp.ring2nest(nside, np.arange(n)), n + hp.ring2nest(nside, np.arange(n))]))):.2e}; "
          f"m > 0 columns not divided by sqrt 2 {rel(Cc, Ca * 2):.2e}")
    G = Y.T @ Y
    w = np.linalg.eigvalsh(G)
    print(f"   basis: {Y.shape[1]} columns on {Y.shape[0]} elements; eigenvalues of Y^T Y from {w.min():.4f} to {w.max():.4f} "
          f"(condition number {w.max() / w.min():.1f}); max |Y^T Y - 1| {np.abs(G - np.eye(G.shape[0])).max():.3f}")

print("# lcdm_cmb_precision itself (CAMB spectra), nside 4")
for target in (60.0, 600.0):
    WM = WindowMatcher(4, target_arcmin=target)
    prec, cov = lcdm_cmb_precision(WM)
    import camb
    pars = camb.set_params(H0=67.36, ombh2=0.02237, omch2=0.1200, tau=0.0544, As=2.1e-9, ns=0.9649, lmax=max(300, WM.lmax))
    cls = camb.get_results(pars).get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=True)["total"]
    bl = hp.gauss_beam(np.radians(target / 60.0), lmax=WM.lmax, pol=True)[:, 1]
    pw = hp.pixwin(4, pol=True, lmax=WM.lmax)[1]
    w2 = (bl * pw) ** 2
    Ca = analytic_cov(4, cls[:WM.lmax + 1, 1] * w2, cls[:WM.lmax + 1, 2] * w2)
    n = cov.shape[0]
    print(f"target {target:.0f} arcmin: max |cov - floor - analytic| / max |analytic| {rel(cov - 0.01 ** 2 * np.eye(n), Ca):.2e}; "
          f"rms per element {np.sqrt(np.mean(np.diag(cov))):.4f} uK; prec * cov - 1: {np.abs(prec @ cov - np.eye(n)).max():.1e}; "
          f"eigenvalues of cov {np.linalg.eigvalsh(cov).min():.2e} to {np.linalg.eigvalsh(cov).max():.2e} (floor 1e-4)")
    print(f"   controls: EE and BB columns swapped {rel(cov - 1e-4 * np.eye(n), analytic_cov(4, cls[:WM.lmax + 1, 2] * w2, cls[:WM.lmax + 1, 1] * w2)):.2e}; "
          f"window not squared {rel(cov - 1e-4 * np.eye(n), analytic_cov(4, cls[:WM.lmax + 1, 1] * bl * pw, cls[:WM.lmax + 1, 2] * bl * pw)):.2e}; "
          f"temperature pixel window {rel(cov - 1e-4 * np.eye(n), analytic_cov(4, cls[:WM.lmax + 1, 1] * (bl * hp.pixwin(4, pol=True, lmax=WM.lmax)[0]) ** 2, cls[:WM.lmax + 1, 2] * (bl * hp.pixwin(4, pol=True, lmax=WM.lmax)[0]) ** 2)):.2e}")
print("# CAMB spectra used: D_l^EE = l(l+1) C_l / 2 pi at l = 2, 4, 6, 10, 20, 47 [uK^2]:",
      " ".join(f"{l * (l + 1) * cls[l, 1] / (2 * np.pi):.4f}" for l in (2, 4, 6, 10, 20, 47)))
print("#                   C_l^BB at l = 2, 10, 47 [uK^2]:", " ".join(f"{cls[l, 2]:.3e}" for l in (2, 10, 47)),
      "; white-floor equivalent C = floor^2 Omega_pix at nside 16:", f"{1e-4 * 4 * np.pi / 3072:.3e}")
