"""Build the E022 low-resolution surrogate of the BeyondPlanck XV polarized beta_s data.

Output `surrogate.npz` in --out with, at HEALPix nside 16 (RING), Stokes Q and U, uK_RJ:
  d[band, stokes, pix], sigma[band, stokes, pix], freqs[band] (GHz), band names,
  amp[comp, stokes, pix] (synch at 30 GHz uK_RJ, dust at 353 GHz uK_RJ, cmb uK_RJ-equivalent
  per band is applied in the analysis; stored here in uK_CMB), T_d[pix], beta_d (scalar),
  region[pix] (0 High Latitude, 1 Spur, 2 Galactic Center, 3 Plane), beta_bp[pix], sd_bp[pix].
Every map is smoothed to 60 arcmin from its native beam and degraded by pixel averaging; WMAP
nside-16 maps are used as provided. Per-band intermediate files make reruns cheap. The WMAP
stage runs when the LAMBDA files are present (see E022 Setup).
"""
import argparse, os, sys
import numpy as np
import healpy as hp
from astropy.io import fits

T_CMB = 2.7255
H_OVER_K = 0.0479924  # h/k in K/GHz
NSIDE_OUT = 16  # overridden by --nside
FWHM_OUT = 60.0  # arcmin
LMAX = 600
BP = os.path.expanduser("~/data/beyondplanck")
WMAP = os.path.expanduser("~/data/wmap9_lowres")


def cmb_to_rj(nu_ghz):
    x = H_OVER_K * nu_ghz / T_CMB
    return x**2 * np.exp(x) / np.expm1(x) ** 2


def col(fn, name, hdu=1):
    with fits.open(fn) as h:
        d = np.asarray(h[hdu].data[name], dtype=np.float64).ravel()
        order = h[hdu].header.get("ORDERING", "RING")
    if order.upper().startswith("NEST"):
        d = hp.reorder(d, n2r=True)
    return d


def smooth_pol(q, u, fwhm_in_arcmin, lmax=LMAX):
    """Spin-2 Gaussian smoothing from fwhm_in to FWHM_OUT (arcmin); no-op if already there."""
    extra = FWHM_OUT**2 - fwhm_in_arcmin**2
    if extra <= 1.0:
        return q, u
    fwhm = np.radians(np.sqrt(extra) / 60.0)
    _, qs, us = hp.smoothing([np.zeros_like(q), q, u], fwhm=fwhm, lmax=lmax, pol=True, iter=1)
    return qs, us


def degrade(m):
    return hp.ud_grade(m, NSIDE_OUT, order_in="RING", order_out="RING")


def stage_band(out, name, q, u, fwhm_in, nu, to_rj, sq=None, su=None, mc_rms=None, seed=0):
    """Smooth, degrade and convert one band; noise sd either given at nside 16 (sq, su, uK_CMB)
    or Monte-Carlo propagated from a native-resolution white-noise rms map (mc_rms = (rms_q,
    rms_u, nside_mc, nreal))."""
    fn = os.path.join(out, f"band_{name}.npz")
    if os.path.exists(fn):
        print(f"  {name}: cached"); return
    qs, us = smooth_pol(q, u, fwhm_in)
    q16, u16 = degrade(qs), degrade(us)
    if mc_rms is not None:
        rq, ru, nside_mc, nreal = mc_rms[:4]; fwhm_noise = mc_rms[4] if len(mc_rms) > 4 else 0.0
        rng = np.random.default_rng(seed)
        acc_q = np.zeros((nreal, hp.nside2npix(NSIDE_OUT))); acc_u = acc_q.copy()
        for i in range(nreal):
            nq = rq * rng.standard_normal(rq.size); nu_ = ru * rng.standard_normal(ru.size)
            nqs, nus = smooth_pol(nq, nu_, fwhm_noise)  # the kernel the map itself received
            acc_q[i], acc_u[i] = degrade(nqs), degrade(nus)
        sq, su = acc_q.std(axis=0, ddof=1), acc_u.std(axis=0, ddof=1)
    np.savez(fn, q=q16 * to_rj, u=u16 * to_rj, sq=sq * to_rj, su=su * to_rj, nu=nu)
    print(f"  {name}: nu {nu} GHz, uK_CMB->uK_RJ {to_rj:.4f}, |Q| median {np.median(np.abs(q16*to_rj)):.3g}, "
          f"sigma_Q median {np.median(sq*to_rj):.3g} uK_RJ")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.expanduser("~/data/bpxv_surrogate"))
    ap.add_argument("--nreal", type=int, default=60)
    ap.add_argument("--nside", type=int, default=16, help="output nside (WMAP nside-16 maps only at 16)")
    ap.add_argument("--no-wmap", action="store_true")
    a = ap.parse_args()
    global NSIDE_OUT
    NSIDE_OUT = a.nside
    os.makedirs(a.out, exist_ok=True)

    print("LFI bands (BP release maps, BP nside-16 60 arcmin rms)")
    for name, nside, fwhm, nu in (("030", 512, 32.29, 28.4), ("044", 512, 27.0, 44.1), ("070", 1024, 13.25, 70.1)):
        fn = f"{BP}/BP_{name}_IQU_n{nside:04d}_v2.fits"
        q, u = col(fn, "Q_MEAN"), col(fn, "U_MEAN")
        if NSIDE_OUT == 16:
            rf = f"{BP}/BP7_{name}_rms_n0016_60arcmin_uK_true.fits"
            sq, su = col(rf, "Q_POLARISATION"), col(rf, "U_POLARISATION")
            stage_band(a.out, name, q, u, fwhm, nu, cmb_to_rj(nu), sq=sq, su=su)
        else:  # no BP simulated rms at this nside: white propagation (Q7 factor 1.6 applies)
            rq, ru = col(fn, "Q_RMS"), col(fn, "U_RMS")
            nside_mc = 256; fac = (nside // nside_mc) ** 2
            rq = np.sqrt(hp.ud_grade(rq**2, nside_mc) / fac); ru = np.sqrt(hp.ud_grade(ru**2, nside_mc) / fac)
            # the white-noise realization is smoothed with the map's own kernel, so the MC step uses fwhm 0
            # on noise already at the map's resolution: emulate by passing mc_rms and fwhm_in through stage_band
            stage_band(a.out, name, q, u, fwhm, nu, cmb_to_rj(nu), mc_rms=(rq, ru, nside_mc, a.nreal, fwhm), seed=int(name))

    print("LFI white-noise variant (Monte-Carlo propagation of the release Q_RMS/U_RMS maps)")
    for name, nside, fwhm, nu in (() if NSIDE_OUT != 16 else (("030", 512, 32.29, 28.4), ("044", 512, 27.0, 44.1), ("070", 1024, 13.25, 70.1))):
        fn = os.path.join(a.out, f"band_{name}_white.npz")
        if os.path.exists(fn):
            print(f"  {name} white: cached"); continue
        mf = f"{BP}/BP_{name}_IQU_n{nside:04d}_v2.fits"
        rq, ru = col(mf, "Q_RMS"), col(mf, "U_RMS")
        nside_mc = 256; fac = (nside // nside_mc) ** 2
        rq = np.sqrt(hp.ud_grade(rq**2, nside_mc) / fac); ru = np.sqrt(hp.ud_grade(ru**2, nside_mc) / fac)
        rng = np.random.default_rng(int(name))
        acc = np.zeros((2, a.nreal, hp.nside2npix(NSIDE_OUT)))
        for i in range(a.nreal):
            nq = rq * rng.standard_normal(rq.size); nu_ = ru * rng.standard_normal(ru.size)
            nqs, nus = smooth_pol(nq, nu_, fwhm)  # the map's own smoothing step
            acc[0, i], acc[1, i] = degrade(nqs), degrade(nus)
        z = np.load(os.path.join(a.out, f"band_{name}.npz"))
        np.savez(fn, q=z["q"], u=z["u"], sq=acc[0].std(axis=0, ddof=1) * cmb_to_rj(nu), su=acc[1].std(axis=0, ddof=1) * cmb_to_rj(nu), nu=nu)
        print(f"  {name} white: sigma_Q median {np.median(acc[0].std(axis=0, ddof=1) * cmb_to_rj(nu)):.3g} uK_RJ (BP rms_true {np.median(z['sq']):.3g})")

    print("HFI 353 (NPIPE A split, Monte-Carlo white-noise rms)")
    if not os.path.exists(os.path.join(a.out, "band_353.npz")):
        fn = f"{BP}/npipe/npipe6v20A_353_map_uK.fits"
        q, u = col(fn, "Q-POLARISATION"), col(fn, "U-POLARISATION")
        rf = f"{BP}/npipe/npipe6v20A_353_rms_uK_mcscaled.fits"
        rq, ru = col(rf, "Q-POLARISATION"), col(rf, "U-POLARISATION")
        nside_mc = 256; fac = (2048 // nside_mc) ** 2
        q, u = hp.ud_grade(q, nside_mc), hp.ud_grade(u, nside_mc)
        rq = np.sqrt(hp.ud_grade(rq**2, nside_mc) / fac); ru = np.sqrt(hp.ud_grade(ru**2, nside_mc) / fac)
        stage_band(a.out, "353", q, u, 4.94, 353.0, cmb_to_rj(353.0), mc_rms=(rq, ru, nside_mc, a.nreal), seed=353)
    else:
        print("  353: cached")

    print("WMAP Ka, Q, V (nine-year nside-16 band maps, diagonal of N from the dense N^-1)")
    for name, nu in (() if (a.no_wmap or NSIDE_OUT != 16) else (("Ka", 33.0), ("Q", 40.6), ("V", 60.8))):
        fn = os.path.join(a.out, f"band_W{name}.npz")
        if os.path.exists(fn):
            print(f"  W{name}: cached"); continue
        mf = f"{WMAP}/wmap_band_qumap_r4_9yr_{name}_v5.fits"; nf = f"{WMAP}/wmap_band_quninv_r4_9yr_{name}_v5.fits"
        if not (os.path.exists(mf) and os.path.exists(nf)):
            print(f"  W{name}: files missing, skipped"); continue
        with fits.open(mf) as h:
            names = h[1].columns.names; hdr = h[1].header
            print(f"  W{name} map columns {names}, ordering {hdr.get('ORDERING')}, units {[hdr.get(f'TUNIT{i+1}') for i in range(len(names))]}")
            qn = [n for n in names if n.upper().startswith("Q")][0]; un = [n for n in names if n.upper().startswith("U")][0]
        q, u = col(mf, qn), col(mf, un)
        with fits.open(nf) as h:
            ninv = np.asarray(h[0].data if h[0].data is not None else h[1].data, dtype=np.float64)
            nhdr = h[0].header if h[0].data is not None else h[1].header
            print(f"  W{name} ninv shape {ninv.shape}, ordering {nhdr.get('ORDERING', 'unknown')}")
        npix = hp.nside2npix(NSIDE_OUT); assert ninv.shape == (2 * npix, 2 * npix)
        cov = np.linalg.inv(ninv); dvar = np.diag(cov)
        nest = str(nhdr.get("ORDERING", "NESTED")).upper().startswith("NEST")
        sq, su = np.sqrt(dvar[:npix]), np.sqrt(dvar[npix:])
        if nest:
            sq, su = hp.reorder(sq, n2r=True), hp.reorder(su, n2r=True)
        # WMAP maps are in mK_CMB
        np.savez(fn, q=q * 1e3 * cmb_to_rj(nu), u=u * 1e3 * cmb_to_rj(nu), sq=sq * 1e3 * cmb_to_rj(nu),
                 su=su * 1e3 * cmb_to_rj(nu), nu=nu, ninv_file=nf, ninv_nested=nest)
        print(f"  W{name}: |Q| median {np.median(np.abs(q))*1e3*cmb_to_rj(nu):.3g}, sigma_Q median {np.median(sq)*1e3*cmb_to_rj(nu):.3g} uK_RJ")

    print("Components and regions (BP release maps)")
    fn = os.path.join(a.out, "comps.npz")
    if not os.path.exists(fn):
        sf = f"{BP}/BP_synch_IQU_n1024_v2.fits"; df = f"{BP}/BP_dust_IQU_n1024_v2.fits"; cf = f"{BP}/BP_cmb_IQU_n1024_v2.fits"
        sq, su = degrade(col(sf, "Q_MEAN")), degrade(col(sf, "U_MEAN"))
        dq, du = smooth_pol(col(df, "Q_MEAN"), col(df, "U_MEAN"), 10.0); dq, du = degrade(dq), degrade(du)
        cq, cu = smooth_pol(col(cf, "Q_MEAN"), col(cf, "U_MEAN"), 14.0); cq, cu = degrade(cq), degrade(cu)
        td = degrade(col(df, "QU_T_MEAN")); beta_d = float(np.mean(col(df, "QU_BETA_MEAN")))
        beta_bp = degrade(col(sf, "QU_BETA_MEAN")); sd_bp = degrade(col(sf, "QU_BETA_STDDEV"))
        # regions: nearest plateau in (beta, sd), then GC split by distance from the centre
        plateaus = np.array([[-3.30, 0.101], [-3.17, 0.063], [-3.00, 0.051]])  # prior-only, Spur, Plane
        z = np.stack([beta_bp / 0.1, sd_bp / 0.03], axis=1); pz = plateaus / np.array([0.1, 0.03])
        lab = np.argmin(((z[:, None, :] - pz[None, :, :]) ** 2).sum(-1), axis=1)
        region = np.where(lab == 1, 1, np.where(lab == 2, 3, 0))
        th, ph = hp.pix2ang(NSIDE_OUT, np.arange(hp.nside2npix(NSIDE_OUT)))
        ang = np.degrees(np.arccos(np.clip(np.sin(np.pi / 2 - th) * 0 + np.cos(th - np.pi / 2) * np.cos(ph), -1, 1)))
        region[(region == 0) & (ang < 25.0)] = 2
        # at nside above 16 the labels come from the same plateau rule on the finer BP map
        np.savez(fn, amp=np.array([[sq, su], [dq, du], [cq, cu]]), T_d=td, beta_d=beta_d, region=region,
                 beta_bp=beta_bp, sd_bp=sd_bp)
        for r, nm in enumerate(("High Latitude", "Spur", "Galactic Center", "Plane")):
            m = region == r
            print(f"  region {r} {nm:16s}: {m.sum():4d} pixels, BP beta mean {beta_bp[m].mean():.3f}, sd map mean {sd_bp[m].mean():.3f}")
        print(f"  beta_d {beta_d:.4f}, T_d median {np.median(td):.2f} K; synch |Q| median {np.median(np.abs(sq)):.3g} uK_RJ@30, dust {np.median(np.abs(dq)):.3g} uK_RJ@353, cmb {np.median(np.abs(cq)):.3g} uK_CMB")
    else:
        print("  comps: cached")

    # assemble
    bands = [b for b in ("030", "044", "070", "WKa", "WQ", "WV", "353") if os.path.exists(os.path.join(a.out, f"band_{b}.npz"))]
    d = []; s = []; nus = []
    for b in bands:
        z = np.load(os.path.join(a.out, f"band_{b}.npz")); d.append([z["q"], z["u"]]); s.append([z["sq"], z["su"]]); nus.append(float(z["nu"]))
    c = np.load(fn)
    s_white = []
    for b in bands:
        wf = os.path.join(a.out, f"band_{b}_white.npz")
        z = np.load(wf if os.path.exists(wf) else os.path.join(a.out, f"band_{b}.npz")); s_white.append([z["sq"], z["su"]])
    np.savez(os.path.join(a.out, "surrogate.npz"), d=np.array(d), sigma=np.array(s), sigma_white=np.array(s_white), freqs=np.array(nus), bands=np.array(bands),
             amp=c["amp"], T_d=c["T_d"], beta_d=c["beta_d"], region=c["region"], beta_bp=c["beta_bp"], sd_bp=c["sd_bp"], nside=NSIDE_OUT)
    print("assembled surrogate.npz with bands", bands)


if __name__ == "__main__":
    main()
