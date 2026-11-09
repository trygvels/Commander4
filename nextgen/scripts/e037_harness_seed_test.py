"""E037: the exact ridge and marginal targets of the validation simulation over many seeds.

Reuses the SEDs, bands, noise levels and amplitude model of `aux/validate_mcmc_likelihood.py`
(imported by path) and evaluates the exact 2-d posterior of (beta_s, beta_d) on a grid for three
treatments of the same sky and noise draw: the harness (smoothed to a common beam, rms = sd of the
smoothed noise, pixels independent), a white-noise control (no smoothing), and the harness with the
noise variance multiplied by the number of pixels per smoothing kernel.
"""
import argparse
import os
import sys

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "aux"))
import validate_mcmc_likelihood as V  # noqa: E402


def posterior(d, w, use_logdet):
    """Mean and sd of (beta_s, beta_d) on a refined grid; d, w are (nband, npix)."""
    def logp(bs, bd):
        M = V.mixing(bs, bd)
        A = np.einsum("bp,bi,bj->pij", w, M, M)
        b = np.einsum("bp,bi->pi", w * d, M)
        quad = np.einsum("pi,pi->p", b, np.linalg.solve(A, b[..., None])[..., 0])
        logdet = np.linalg.slogdet(A)[1] if use_logdet else 0.0
        return np.sum(0.5 * quad - 0.5 * logdet)
    lo, hi = np.array(V.BOUNDS["sync"] + V.BOUNDS["dust"]).reshape(2, 2).T

    def grid(lo_, hi_, n):
        gs, gd = np.linspace(lo_[0], hi_[0], n), np.linspace(lo_[1], hi_[1], n)
        L = np.array([[logp(s_, t_) for t_ in gd] for s_ in gs])
        p = np.exp(L - L.max()); p /= p.sum()
        ps, pd = p.sum(1), p.sum(0)
        mean = np.array([np.sum(ps * gs), np.sum(pd * gd)])
        sd = np.sqrt([np.sum(ps * (gs - mean[0]) ** 2), np.sum(pd * (gd - mean[1]) ** 2)])
        return mean, sd, np.array([gs[1] - gs[0], gd[1] - gd[0]])
    mean, sd, step = grid(lo, hi, 21)
    for _ in range(3):
        half = np.maximum(6 * sd, 2 * step)
        mean, sd, step = grid(np.maximum(lo, mean - half), np.minimum(hi, mean + half), 25)
    return mean, sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nseed", type=int, default=40)
    a = ap.parse_args()
    nside, npix = V.NSIDE, hp.nside2npix(V.NSIDE)
    fwhm, lmax = np.radians(V.FWHM_ARCMIN / 60.0), 3 * V.NSIDE - 1
    sig_b = fwhm / np.sqrt(8.0 * np.log(2.0))
    neff = 4.0 * np.pi * sig_b ** 2 / hp.nside2pixarea(nside)
    truth = np.array([V.TRUTH["sync"], V.TRUTH["dust"]])
    lat = np.abs(90 - np.degrees(hp.pix2ang(nside, np.arange(npix))[0]))
    t = np.exp(-lat / 15.0)
    M0 = V.mixing(*truth)
    print(f"# E037: nside {nside}, common beam {V.FWHM_ARCMIN} arcmin, pixels per smoothing kernel {neff:.2f}, {a.nseed} seeds")
    print("# z = (posterior mean - input) / posterior sd for beta_s, beta_d; variants h (harness), w (white control), "
          "c (harness, noise variance times N_eff); targets ridge and flat-prior marginal")
    Z = {}
    for seed in range(a.nseed):
        rng = np.random.default_rng([37, seed])
        amps = np.stack([rng.normal(0, 70, npix), 10 ** (0.3 + 2.0 * t + 0.1 * rng.normal(size=npix)),
                         10 ** (0.5 + 1.5 * t + 0.1 * rng.normal(size=npix))])
        data = M0 @ amps
        white = np.array([rng.normal(0, rms, npix) for rms in V.RMS])
        dh, wh = np.empty_like(data), np.empty_like(data)
        for i in range(len(V.FREQS)):
            n = hp.smoothing(white[i], fwhm=fwhm, lmax=lmax)
            dh[i] = hp.smoothing(data[i], fwhm=fwhm, lmax=lmax) + n
            wh[i] = 1.0 / n.std() ** 2
        ww = np.array([np.full(npix, 1.0 / rms ** 2) for rms in V.RMS])
        cases = {"h": (dh, wh), "w": (data + white, ww), "c": (dh, wh / neff)}
        line = [f"{seed:3d}"]
        for cn, (d, w) in cases.items():
            for tn, ld in (("ridge", False), ("marg", True)):
                if tn == "marg" and cn == "c":
                    continue
                m, s = posterior(d, w, ld)
                z = (m - truth) / s
                Z.setdefault((cn, tn), []).append(np.concatenate([z, s]))
                line.append(f"{cn}/{tn} {z[0]:+6.2f} {z[1]:+6.2f}")
        print("  ".join(line), flush=True)
    print("\n# summary: variant/target, mean z (se) and sd of z for beta_s and beta_d, median posterior sd")
    for (cn, tn), v in Z.items():
        v = np.array(v); n = len(v)
        print(f"{cn}/{tn:5s} beta_s mean z {v[:, 0].mean():+6.2f} ({v[:, 0].std(ddof=1) / np.sqrt(n):.2f}) sd z "
              f"{v[:, 0].std(ddof=1):5.2f} | beta_d mean z {v[:, 1].mean():+6.2f} ({v[:, 1].std(ddof=1) / np.sqrt(n):.2f}) "
              f"sd z {v[:, 1].std(ddof=1):5.2f} | median sd {np.median(v[:, 2]):.4f} {np.median(v[:, 3]):.4f}")


if __name__ == "__main__":
    main()
