"""Loaders for the low-resolution polarization data. Everything is returned in NESTED ordering, as
(Q elements, then U elements), in uK thermodynamic, at nside 16 unless degraded.

Run any script that uses these on downloaded files with `python -I`.
"""
import numpy as np
import healpy as hp
from astropy.io import fits

from .likelihood import BandData
from .spectra import Band, wmap_band, lfi_band, hfi_band

WMAP_DIRS = {"K": "/Users/trygvels/data/wmap9_lowres_K", "Ka": "/Users/trygvels/data/wmap9_lowres",
             "Q": "/Users/trygvels/data/wmap9_lowres", "V": "/Users/trygvels/data/wmap9_lowres"}
WMAP_BANDPASS = "/Users/trygvels/data/wmap9_bandpass/extracted"
PR4_353 = "/Users/trygvels/data/planck_pr4/HFI_SkyMap_353_2048_R4.00_full.fits"
PR4_LFI = "/Users/trygvels/data/planck_pr4_lfi/LFI_SkyMap_{freq:03d}_1024_R4.00_{split}.fits"
DUST_BETA, DUST_T, DUST_REF = 1.53, 19.6, 353.0


def load_wmap(band):
    """WMAP nine-year res 4 Q/U map and dense inverse covariance: d (6144,) in uK, Ninv in uK^-2."""
    base = WMAP_DIRS[band]
    with fits.open(f"{base}/wmap_band_qumap_r4_9yr_{band}_v5.fits") as h:
        t = h[1].data
        d = np.concatenate([np.asarray(t["Q_POLARISATION"], float).ravel(),
                            np.asarray(t["U_POLARISATION"], float).ravel()]) * 1e3
        assert str(h[1].header.get("ORDERING", "")).upper().startswith("NEST")
    with fits.open(f"{base}/wmap_band_quninv_r4_9yr_{band}_v5.fits") as h:
        ninv = np.asarray(h[0].data, dtype=np.float64) * 1e-6
    return d, 0.5 * (ninv + ninv.T)


def load_planck_qu(path, nside_out=16, scale=1e6):
    """Planck Q/U map with its per-pixel Q/U covariance, brought to nside_out by an inverse-covariance
    weighted average over sub-pixels. Returns d (2n,) in uK and the dense Ninv (2n, 2n), which holds the
    Q/U cross term of each coarse pixel."""
    q, u, qq, qu, uu = [np.asarray(m, dtype=np.float64) for m in
                        hp.read_map(path, field=(1, 2, 7, 8, 9), nest=True, dtype=np.float64)]
    bad = (q < -1e20) | (qq <= 0) | (uu <= 0)
    q, u = q * scale, u * scale
    qq, qu, uu = qq * scale ** 2, qu * scale ** 2, uu * scale ** 2
    det = qq * uu - qu ** 2
    wqq, wuu, wqu = uu / det, qq / det, -qu / det
    for w in (wqq, wuu, wqu):
        w[bad] = 0.0
    k = (hp.npix2nside(q.size) // nside_out) ** 2
    n = hp.nside2npix(nside_out)
    red = lambda x: x.reshape(n, k).sum(1)
    Wqq, Wuu, Wqu = red(wqq), red(wuu), red(wqu)
    bq, bu = red(wqq * q + wqu * u), red(wqu * q + wuu * u)
    D = Wqq * Wuu - Wqu ** 2
    dq, du = (Wuu * bq - Wqu * bu) / D, (Wqq * bu - Wqu * bq) / D
    ninv = np.zeros((2 * n, 2 * n))
    i = np.arange(n)
    ninv[i, i], ninv[n + i, n + i], ninv[i, n + i], ninv[n + i, i] = Wqq, Wuu, Wqu, Wqu
    return np.concatenate([dq, du]), ninv


def degrade(d, ninv, nside_in, nside_out):
    """Degrade a NESTED (Q, U) data vector and dense inverse covariance, assuming the signal is constant
    inside a coarse pixel: Ninv_c = P^T Ninv P, d_c = Ninv_c^-1 P^T Ninv d."""
    k = (nside_in // nside_out) ** 2
    n_in, n_out = hp.nside2npix(nside_in), hp.nside2npix(nside_out)
    idx = np.concatenate([np.arange(n_in) // k, n_out + np.arange(n_in) // k])
    P = np.zeros((2 * n_in, 2 * n_out))
    P[np.arange(2 * n_in), idx] = 1.0
    nc = P.T @ ninv @ P
    return np.linalg.solve(nc, P.T @ (ninv @ d)), nc


def wmap_banddata(band, nu_ref, nside=16, dust_beta=DUST_BETA):
    d, ninv = load_wmap(band)
    if nside != 16:
        d, ninv = degrade(d, ninv, 16, nside)
    B = wmap_band(band, WMAP_BANDPASS)
    return BandData(f"WMAP_{band}", d, ninv, lambda b, B=B: B.power_law(b, nu_ref), B.mbb(dust_beta, DUST_T, DUST_REF))


def planck353_banddata(nu_ref, nside=16, path=PR4_353, bandpass="psb", dust_beta=DUST_BETA):
    """PR4 353 GHz Q/U. bandpass: "psb" (the polarization-sensitive bolometers, which carry the
    polarization), "full" (all bolometers) or "delta" (353.0 GHz; misstates the dust response by 12 to
    14 percent and is kept as a sensitivity arm only)."""
    d, ninv = load_planck_qu(path, nside)
    B = Band.delta("Planck_353", 353.0) if bandpass == "delta" else hfi_band(353, psb=(bandpass == "psb"))
    return BandData("Planck_353", d, ninv, lambda b, B=B: B.power_law(b, nu_ref), B.mbb(dust_beta, DUST_T, DUST_REF))


def planck_lfi_banddata(freq, nu_ref, nside=16, split="full", dust_beta=DUST_BETA):
    """PR4 LFI Q/U at 30, 44 or 70 GHz with the white-noise covariance of the map file, band-averaged
    spectra from the RIMO bandpass. The covariance has no correlated-noise term: see the ring-half
    check before using it at large scales."""
    d, ninv = load_planck_qu(PR4_LFI.format(freq=int(freq), split=split), nside)
    B = lfi_band(freq)
    return BandData(f"Planck_{int(freq):03d}", d, ninv, lambda b, B=B: B.power_law(b, nu_ref),
                    B.mbb(dust_beta, DUST_T, DUST_REF))


def healpix_regions(nside_region, nside=16):
    """Region index per (Q, U) element: the nside_region pixel that contains each pixel (NESTED)."""
    k = (nside // nside_region) ** 2
    return np.tile(np.arange(hp.nside2npix(nside)) // k, 2)


def latitude_regions(edges, nside=16):
    """Regions by absolute Galactic latitude of the pixel centre; edges in degrees, e.g. (0, 20, 45, 90)."""
    lat = np.abs(hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest=True, lonlat=True)[1])
    return np.tile(np.clip(np.digitize(lat, edges[1:-1]), 0, len(edges) - 2), 2)


def label_regions(path, field=1, nside=16):
    """Regions from a label map (NESTED, nside at or above the target), by majority vote among the
    labelled sub-pixels. Unlabelled pixels (NaN or UNSEEN) form one extra region with the last index.
    Returns the region of each (Q, U) element and the number of regions."""
    lab = np.asarray(hp.read_map(path, field=field, nest=True, dtype=np.float64))
    ok = np.isfinite(lab) & (lab > -1e20)
    ids = np.unique(lab[ok])
    code = np.full(lab.size, ids.size, dtype=int)
    code[ok] = np.searchsorted(ids, lab[ok])
    k = (hp.npix2nside(lab.size) // nside) ** 2
    code = code.reshape(hp.nside2npix(nside), k)
    out = np.array([np.bincount(row[row < ids.size], minlength=ids.size).argmax() if np.any(row < ids.size)
                    else ids.size for row in code])
    nreg = int(out.max()) + 1
    return np.tile(out, 2), nreg


def bpxv_regions(path="/Users/trygvels/data/bpxv_surrogate/comps.npz"):
    """The four BeyondPlanck XV polarization regions at nside 16 (0 High Latitude, 1 Spur, 2 Galactic
    Center, 3 Plane), as stored in RING order by the E022 surrogate, returned NESTED per (Q, U) element."""
    reg_ring = np.load(path, allow_pickle=True)["region"].astype(int)
    return np.tile(hp.reorder(reg_ring, r2n=True), 2)
