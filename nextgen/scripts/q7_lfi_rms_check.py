"""Q7 data point: BP's nside-16, 60 arcmin 'rms_true' (from noise simulations with correlated
noise) against white-noise propagation of the BP Q_RMS/U_RMS maps through the same smoothing
and degrading, for LFI 30, 44, 70."""
import numpy as np, healpy as hp, sys
sys.path.insert(0, "nextgen/scripts")
from prep_bpxv_surrogate import col, smooth_pol, degrade, BP
rng = np.random.default_rng(7)
for name, nside, fwhm in (("030", 512, 32.29), ("044", 512, 27.0), ("070", 1024, 13.25)):
    fn = f"{BP}/BP_{name}_IQU_n{nside:04d}_v2.fits"
    rq, ru = col(fn, "Q_RMS"), col(fn, "U_RMS")
    # white noise at native nside, pixel average to 256 first (exact for white noise), then smooth
    nmc = 256; fac = (nside // nmc) ** 2
    rq256 = np.sqrt(hp.ud_grade(rq**2, nmc) / fac); ru256 = np.sqrt(hp.ud_grade(ru**2, nmc) / fac)
    nreal = 40; acc = np.zeros((2, nreal, hp.nside2npix(16)))
    for i in range(nreal):
        q = rq256 * rng.standard_normal(rq256.size); u = ru256 * rng.standard_normal(ru256.size)
        qs, us = smooth_pol(q, u, fwhm)  # kernel sqrt(60^2 - fwhm^2): the map's own smoothing step
        acc[0, i], acc[1, i] = degrade(qs), degrade(us)
    white = acc.std(axis=1, ddof=1)
    tq, tu = col(f"{BP}/BP7_{name}_rms_n0016_60arcmin_uK_true.fits", "Q_POLARISATION"), col(f"{BP}/BP7_{name}_rms_n0016_60arcmin_uK_true.fits", "U_POLARISATION")
    for s, (t, w) in enumerate(((tq, white[0]), (tu, white[1]))):
        r = t / w
        print(f"LFI {name} {'QU'[s]}: rms_true median {np.median(t):.3f} uK, white-propagated median {np.median(w):.3f} uK, "
              f"ratio median {np.median(r):.2f}, 10/90 pct {np.percentile(r,10):.2f}/{np.percentile(r,90):.2f}, max {r.max():.2f}")
