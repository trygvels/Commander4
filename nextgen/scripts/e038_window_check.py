"""E038 (c): bias of region indices from the different beams of the bands, noise off.

    python -I nextgen/scripts/e038_window_check.py [--bands K Ka Q V 353] [--regions bpxv] [--common-beam 60]

Truth: DR1 polarized synchrotron (60 arcmin) plus a synthetic small-scale field continuing its E and B
spectra as power laws to lmax, modulated by the local polarized intensity; index -3.1 everywhere. Dust: DR1
polarized dust (10 arcmin), modified blackbody (1.53, 19.6 K). Each band sees the sky through its own
Gaussian beam, is averaged to nside 16 and fitted with the real inverse covariance as weights.
"""
import argparse
import os
import sys
import time

import healpy as hp
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from lowres import data as D  # noqa: E402
from lowres.fit import fit_profile  # noqa: E402
from lowres.likelihood import BandData, RegionLikelihood  # noqa: E402
from lowres.spectra import g_rj, hfi_band, lfi_band, mbb_rj, wmap_band  # noqa: E402

SYNCH = "/Users/trygvels/data/cosmoglobe_DR1/CG_synch_IQU_n1024_v1.fits"
DUST = "/Users/trygvels/data/cosmoglobe_DR1/CG_dust_IQU_n1024_v1.fits"
FWHM = {"K": 52.8, "Ka": 39.6, "Q": 30.6, "V": 21.0, "30": 32.3, "44": 27.9, "70": 13.2, "353": 4.9}   # arcmin
BETA = -3.1


def beam(fwhm_arcmin, lmax):
    return hp.gauss_beam(np.radians(fwhm_arcmin / 60.0), lmax=lmax, pol=True)[:, 1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bands", nargs="+", default=["K", "Ka", "Q", "V", "353"])
    ap.add_argument("--regions", default="bpxv", choices=["bpxv", "hp2"])
    ap.add_argument("--common-beam", type=float, default=0.0, help="arcmin; if set, all bands get this beam (control)")
    ap.add_argument("--no-small-scales", action="store_true")
    ap.add_argument("--nside-hi", type=int, default=512)
    ap.add_argument("--seed", type=int, default=3803)
    a = ap.parse_args()
    t0 = time.time()
    nh, lmax = a.nside_hi, 2 * a.nside_hi
    # synchrotron truth at 30 GHz, RJ: map (60 arcmin) as alms
    q, u = hp.read_map(SYNCH, field=(1, 2), dtype=np.float64)
    almE, almB = hp.map2alm([np.zeros_like(q), q, u], lmax=lmax, pol=True, iter=1)[1:]
    b60 = beam(60.0, lmax)
    cE, cB = hp.alm2cl(almE), hp.alm2cl(almB)
    ell = np.arange(lmax + 1)
    fitsl = slice(20, 81)
    extra = []
    for c in (cE, cB):
        ctrue = c[fitsl] / b60[fitsl] ** 2
        al, ln0 = np.polyfit(np.log(ell[fitsl]), np.log(ctrue), 1)
        cfit = np.zeros(lmax + 1)
        cfit[2:] = np.exp(ln0) * ell[2:] ** al
        extra.append((al, cfit))
    print(f"# power-law fits to the deconvolved DR1 synchrotron spectra over l 20 to 80: EE index {extra[0][0]:.2f}, BB index {extra[1][0]:.2f}")
    # deconvolved large scales (l <= 100) plus a synthetic field carrying the power the 60 arcmin beam removed
    lcut = 100
    taper = np.where(ell <= lcut, 1.0 / b60, 0.0)
    tE, tB = hp.almxfl(almE, taper), hp.almxfl(almB, taper)
    if not a.no_small_scales:
        rng = np.random.default_rng(a.seed)
        np.random.seed(a.seed)
        cl_hi = [np.where(ell > lcut, extra[i][1], 0.0) for i in (0, 1)]
        gE, gB = hp.synalm(cl_hi[0], lmax=lmax), hp.synalm(cl_hi[1], lmax=lmax)
        # modulate by the local polarized intensity (5 degrees), unit mean square over the sky
        P = np.sqrt(q ** 2 + u ** 2)
        Ps = hp.ud_grade(hp.smoothing(P, fwhm=np.radians(5.0), lmax=128), nh)
        mod = Ps / np.sqrt(np.mean(Ps ** 2))
        _, gq, gu = hp.alm2map([np.zeros_like(gE), gE, gB], nh, lmax=lmax, pol=True)
        sE, sB = hp.map2alm([np.zeros_like(gq), gq * mod, gu * mod], lmax=lmax, pol=True, iter=1)[1:]
        tE, tB = tE + sE, tB + sB
    # dust truth, RJ at 353 GHz (the DR1 polarized dust map is referenced at 70 GHz in RJ units, F13)
    dq, du = hp.read_map(DUST, field=(1, 2), dtype=np.float64)
    fac = mbb_rj(353.0, D.DUST_BETA, D.DUST_T, 70.0)
    dE, dB = hp.map2alm([np.zeros_like(dq), dq * fac, du * fac], lmax=lmax, pol=True, iter=1)[1:]
    b10 = beam(10.0, lmax)

    def band_map(almE_, almB_, fw_true, fw_have):
        """Q, U at nside 16 NESTED of a field with beam fw_have already applied, seen through fw_true."""
        bl = beam(fw_true, lmax) / (beam(fw_have, lmax) if fw_have > 0 else 1.0)
        bl = np.where(np.isfinite(bl) & (bl < 50), bl, 0.0)
        _, mq, mu = hp.alm2map([np.zeros_like(almE_), hp.almxfl(almE_, bl), hp.almxfl(almB_, bl)], nh, lmax=lmax, pol=True)
        return np.concatenate([hp.reorder(hp.ud_grade(mq, 16), r2n=True), hp.reorder(hp.ud_grade(mu, 16), r2n=True)])

    bands = []
    for b in a.bands:
        if b == "353":
            _, ninv = D.load_planck_qu(D.PR4_353, 16)
            B = hfi_band(353, psb=True)
        elif b in ("30", "44", "70"):
            _, ninv = D.load_lfi(int(b), 16)
            B = lfi_band(int(b))
        else:
            _, ninv = D.load_wmap(b)
            B = wmap_band(b, D.WMAP_BANDPASS)
        fw = a.common_beam if a.common_beam > 0 else FWHM[b]
        s0 = B.power_law(BETA, 30.0)[0][0]                    # thermodynamic response per uK_RJ at 30 GHz
        fd = B.mbb(D.DUST_BETA, D.DUST_T, D.DUST_REF)
        d = s0 * band_map(tE, tB, fw, 0.0) + fd * band_map(dE, dB, max(fw, 10.0), 10.0)
        bands.append(BandData(B.name, d, ninv, (lambda be, B=B: B.power_law(be, 23.0)), fd))
        print(f"# {B.name}: beam {fw:.1f} arcmin, rms of the nside-16 signal {np.sqrt(np.mean(d ** 2)):.2f} uK; {time.time() - t0:.0f} s", flush=True)
    region = D.bpxv_regions() if a.regions == "bpxv" else D.healpix_regions(2, 16)
    nreg = int(region.max()) + 1
    L = RegionLikelihood(bands, region, cmb_prior_sd=0.2, use_dust=("353" in a.bands))
    bh, res = fit_profile(L, np.full(nreg, BETA))
    print(f"# bands {a.bands}, regions {a.regions}, common beam {a.common_beam}, small scales {not a.no_small_scales}; "
          f"{res.nit} iterations, chi2_min {2 * res.fun:.2f}; {time.time() - t0:.0f} s")
    print("bias per region (fitted index + 3.1): " + " ".join(f"{x + 3.1:+.4f}" for x in bh))
    print(f"largest abs bias {np.max(np.abs(bh + 3.1)):.4f}; mean {np.mean(bh + 3.1):+.4f}")


if __name__ == "__main__":
    main()
