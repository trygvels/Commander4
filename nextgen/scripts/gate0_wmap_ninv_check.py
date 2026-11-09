"""Gate 0(a) of the polarized track: sanity of a WMAP nine-year res-4 Q/U inverse noise covariance.

Checks symmetry, the spectrum (null space, conditioning), the pixel rms implied by the inverse, how far
the covariance is from white (the poorly determined large-scale modes), and the map's chi-square.

    python -I nextgen/scripts/gate0_wmap_ninv_check.py --ninv N.fits --map M.fits --nu 22.8 [--label K]
"""
import argparse

import healpy as hp
import numpy as np
from astropy.io import fits

H_K, TCMB = 0.0479924466, 2.7255


def cmb_to_rj(nu):
    x = H_K * nu / TCMB
    return x ** 2 * np.exp(x) / np.expm1(x) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ninv", required=True); ap.add_argument("--map", required=True)
    ap.add_argument("--nu", type=float, required=True); ap.add_argument("--label", default="")
    a = ap.parse_args()
    with fits.open(a.ninv) as h:
        N = np.asarray(h[0].data if h[0].data is not None else h[1].data, dtype=np.float64)
        order = h[0].header.get("ORDERING", h[1].header.get("ORDERING", "unknown") if len(h) > 1 else "unknown")
    with fits.open(a.map) as h:
        cols = h[1].columns.names
        q = np.asarray(h[1].data[cols[0]], dtype=np.float64).ravel()
        u = np.asarray(h[1].data[cols[1]], dtype=np.float64).ravel()
        units = [h[1].header.get(f"TUNIT{i}") for i in (1, 2)]
    n = N.shape[0]; npix = n // 2; nside = hp.npix2nside(npix)
    print(f"# {a.label} N^-1 {N.shape}, ordering {order}; map columns {cols[:4]} units {units}; nside {nside}")
    asym = np.abs(N - N.T).max() / np.abs(N).max()
    print(f"symmetry: max |N - N^T| / max |N| = {asym:.2e}")
    N = 0.5 * (N + N.T)
    ev, V = np.linalg.eigh(N)
    nnull = int(np.sum(ev < 1e-10 * ev.max()))
    print(f"eigenvalues: min {ev.min():.3e} max {ev.max():.3e}; below 1e-10 max: {nnull}; negative: {int(np.sum(ev < 0))}; "
          f"condition (positive part) {ev.max() / ev[ev > 1e-10 * ev.max()].min():.3e}")
    zero_rows = int(np.sum(np.abs(N).sum(axis=1) == 0))
    print(f"all-zero rows (masked elements): {zero_rows}")
    pos = ev > 1e-10 * ev.max()
    # pixel rms from the pseudo-inverse, in uK_RJ (maps and N^-1 are in mK thermodynamic)
    Cdiag = np.einsum("ij,j,ij->i", V[:, pos], 1.0 / ev[pos], V[:, pos])
    rms_uK_rj = np.sqrt(Cdiag) * 1e3 * cmb_to_rj(a.nu)
    ok = rms_uK_rj > 0
    print(f"pixel rms from diag(C), uK_RJ: median {np.median(rms_uK_rj[ok]):.3f}, 16/84 pct "
          f"{np.percentile(rms_uK_rj[ok], 16):.3f}/{np.percentile(rms_uK_rj[ok], 84):.3f}")
    white = 1.0 / np.sqrt(np.diag(N)[np.diag(N) > 0]) * 1e3 * cmb_to_rj(a.nu)
    print(f"white-noise rms 1/sqrt(diag N^-1), uK_RJ: median {np.median(white):.3f}; ratio median diag(C)^1/2 over it "
          f"{np.median(rms_uK_rj[ok]) / np.median(white):.3f}")
    # how far from white: variance of the worst modes relative to the median mode
    var = 1.0 / ev[pos]
    vs = np.sort(var)[::-1]
    print(f"mode variances relative to the median mode: worst {vs[0] / np.median(var):.1f}, 2nd {vs[1] / np.median(var):.1f}, "
          f"5th {vs[4] / np.median(var):.1f}, 10th {vs[9] / np.median(var):.1f}, 50th {vs[49] / np.median(var):.1f}, "
          f"100th {vs[99] / np.median(var):.2f}; fraction of total variance in the worst 10 modes "
          f"{vs[:10].sum() / var.sum():.3f}, worst 100 {vs[:100].sum() / var.sum():.3f}")
    # multipole content of the worst modes (Q part as a scalar map, crude)
    nest = str(order).upper().startswith("NEST")
    for k in range(3):
        vec = V[:, pos][:, np.argsort(var)[::-1][k]]
        mq = hp.reorder(vec[:npix], n2r=True) if nest else vec[:npix]
        cl = hp.anafast(mq, lmax=3 * nside - 1)
        ell = np.arange(cl.size)
        print(f"  worst mode {k + 1}: mean multipole of its Q part {np.sum(ell * cl * (2 * ell + 1)) / np.sum(cl * (2 * ell + 1)):.1f}")
    d = np.concatenate([q, u])
    chi2 = float(d @ N @ d)
    print(f"map chi-square d^T N^-1 d = {chi2:.1f} for {int(pos.sum())} modes (signal plus noise; above the mode count "
          f"when the sky is detected)")


if __name__ == "__main__":
    main()
