"""E034: checks of the map ingest path raised by the 2026-11-09 pipeline audit (items B1/D5, B14/D9).

Two independent tests, selected with --test:

  floor   For every enabled band of a Commander4 parameter file, read the band through
          Commander4's own reader (`read_data_map_from_file`: unit conversion, optional
          `add_signal_fraction_to_rms`, ud_grade to eval_nside, noise-weighted smoothing with the
          variance propagated) and report the per-pixel rms 1/sqrt(inv_var) that reaches the
          likelihood at |b| > 30 deg, (a) as configured, (b) with the floor removed, (c) the
          degraded rms before the smoothing with the nominal floor added in quadrature after the
          degradation, which is what the "rms (with floor)" column of the E014 high-latitude table
          is suspected to hold. Also the nominal floor (fraction times the full-sky mean of the
          absolute native map).

  order   A Gaussian CMB realization at nside 512 (lmax 1535) observed by bands with different
          Gaussian beams, brought to nside 64 and a 120 arcmin beam (a) as the code does it:
          ud_grade then Gaussian smoothing by sqrt(120^2 - fwhm^2) on the coarse map (the
          noise-weighted smoothing reduces to this for uniform weights), (b) reference: the same
          alms with the 120 arcmin beam synthesized at nside 64 pixel centres, (c) reference with
          the nside-64 pixel window as well, (d) smooth at nside 512 then ud_grade. Reports the
          rms of (code - reference) per band and of the between-band differences, which is the
          band-dependent part that no CMB amplitude can absorb.

    python -I nextgen/scripts/e034_ingest_checks.py --test floor --param FILE.yml
    python -I nextgen/scripts/e034_ingest_checks.py --test order [--cl CAMB_OR_NONE] [--seed 1]
"""
import argparse
import copy

import healpy as hp
import numpy as np


def gal_lat_mask(nside, bcut):
    _, lat = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    return np.abs(lat) > bcut


def test_floor(a):
    from commander4.file_io.map_reader import read_data_map_from_file, retrieve_map_from_fits_file
    from commander4.parameters.parse import load_params

    params, _, _ = load_params(a.param)
    bands = params.compsep.bands
    common = float(params.compsep.common_res_fwhm)
    print(f"# parameter file {a.param}; common_res_fwhm {common} arcmin; rms in uK_RJ at |b| > 30 deg")
    print("# cols: band, freq, fwhm, native nside, floor fraction, nominal floor (frac * mean|map|), then the"
          " median over |b|>30 of:")
    print("#   A  rms the likelihood gets (reader as configured)")
    print("#   B  the same without add_signal_fraction_to_rms")
    print("#   C  degraded rms before smoothing (no floor), i.e. native white noise per coarse pixel")
    print("#   D  sqrt(C^2 + nominal floor^2): degraded rms with the nominal floor, no smoothing reduction")
    print("#   A/D and the effective floor that reaches the likelihood, sqrt(A^2 - B^2)")
    print(f"{'band':28s} {'GHz':>7s} {'fwhm':>6s} {'nside':>5s} {'frac':>5s} {'nominal':>10s} "
          f"{'A':>10s} {'B':>10s} {'C':>10s} {'D':>10s} {'A/D':>6s} {'eff floor':>10s}")
    for name in bands:
        b = bands[name]
        if not b.enabled:
            continue
        if "_name" not in b:
            b._name = name
        frac = float(b.add_signal_fraction_to_rms) if "add_signal_fraction_to_rms" in b else 0.0
        dm = read_data_map_from_file(b, params)
        keep = gal_lat_mask(dm.nside, 30.0)
        A = np.median(1.0 / np.sqrt(dm.map_inv_var[0][keep]))
        b2 = copy.deepcopy(b)
        if "add_signal_fraction_to_rms" in b2:
            del b2["add_signal_fraction_to_rms"]
        dm2 = read_data_map_from_file(b2, params)
        B = np.median(1.0 / np.sqrt(dm2.map_inv_var[0][keep]))
        # degraded rms without smoothing and without floor, straight from the files
        sky = retrieve_map_from_fits_file(b, "I", "signal")
        rms = np.asarray(retrieve_map_from_fits_file(b, "I", "rms"), dtype=np.float64)
        nside_in = hp.npix2nside(sky.size)
        nominal = frac * float(np.nanmean(np.abs(sky)))
        ivar = np.divide(1.0, rms ** 2, out=np.zeros_like(rms), where=rms > 0)
        ivar64 = hp.ud_grade(ivar, dm.nside, power=-2) if nside_in != dm.nside else ivar
        C = np.median(1.0 / np.sqrt(ivar64[keep]))
        D = np.sqrt(C ** 2 + nominal ** 2)
        eff = np.sqrt(max(A ** 2 - B ** 2, 0.0))
        print(f"{name:28s} {float(b.freq):7.3f} {float(b.fwhm):6.1f} {nside_in:5d} {frac:5.2f} {nominal:10.4g} "
              f"{A:10.4g} {B:10.4g} {C:10.4g} {D:10.4g} {A / D:6.3f} {eff:10.4g}", flush=True)


def test_order(a):
    rng = np.random.default_rng(a.seed)
    nside_in, nside_out, lmax = 512, 64, 3 * 512 - 1
    ell = np.arange(lmax + 1)
    if a.cl:
        cl = np.loadtxt(a.cl)[: lmax + 1]
    else:
        # Smooth stand-in for the TT spectrum, uK_CMB^2: D_l = 1000 uK^2 flat to l = 30, a first
        # peak of 5500 uK^2 at l = 220 and an exponential damping tail. Only the amount of power
        # around l = 100 to 600 matters for this test.
        D = 1000.0 + 4500.0 * np.exp(-0.5 * ((ell - 220.0) / 90.0) ** 2) \
            + 1500.0 * np.exp(-0.5 * ((ell - 540.0) / 120.0) ** 2)
        D *= np.exp(-(ell / 1400.0) ** 2)
        cl = np.zeros(lmax + 1)
        cl[2:] = 2.0 * np.pi * D[2:] / (ell[2:] * (ell[2:] + 1.0))
    np.random.seed(a.seed)
    alm = hp.synalm(cl, lmax=lmax, new=True)
    target = np.deg2rad(a.target / 60.0)
    lmax_out = 3 * nside_out - 1
    keep_l = hp.Alm.getlm(lmax)[0] <= lmax_out
    idx_out = hp.Alm.getidx(lmax_out, *hp.Alm.getlm(lmax_out))
    alm_lo = np.zeros(hp.Alm.getsize(lmax_out), dtype=complex)
    l_all, m_all = hp.Alm.getlm(lmax)
    alm_lo[hp.Alm.getidx(lmax_out, l_all[keep_l], m_all[keep_l])] = alm[keep_l]
    b_t = hp.gauss_beam(target, lmax=lmax_out)
    pw = hp.pixwin(nside_out, lmax=lmax_out)
    ref = hp.alm2map(hp.almxfl(alm_lo, b_t), nside_out, lmax=lmax_out)             # (b)
    ref_pw = hp.alm2map(hp.almxfl(alm_lo, b_t * pw), nside_out, lmax=lmax_out)     # (c)
    print(f"# CMB-only, nside {nside_in} -> {nside_out}, target beam {a.target} arcmin, seed {a.seed}; uK_CMB")
    print(f"# rms of the reference map (b): {ref.std():.3f}")
    print("# per band: rms of code path minus reference (b: no pixel window, c: with the nside-64 pixel window),"
          " and of smooth-then-degrade minus the same references")
    print(f"{'fwhm':>6s} {'code-b':>9s} {'code-c':>9s} {'s.t.d.-b':>9s} {'s.t.d.-c':>9s} {'unsmoothed alias rms':>22s}")
    code, std = {}, {}
    for fwhm in a.fwhm:
        bl = hp.gauss_beam(np.deg2rad(fwhm / 60.0), lmax=lmax)
        native = hp.alm2map(hp.almxfl(alm, bl * hp.pixwin(nside_in, lmax=lmax)), nside_in, lmax=lmax)
        coarse = hp.ud_grade(native, nside_out)
        extra = np.deg2rad(np.sqrt(a.target ** 2 - fwhm ** 2) / 60.0)
        code[fwhm] = hp.smoothing(coarse, fwhm=extra)                               # (a) the code
        fine = hp.smoothing(native, fwhm=extra, lmax=lmax)
        std[fwhm] = hp.ud_grade(fine, nside_out)                                    # (d)
        # aliasing content of the coarse map before the smoothing: coarse minus the band-limited
        # (l <= 191) version of the pixel-averaged, beam-smoothed sky
        bl_lo = hp.gauss_beam(np.deg2rad(fwhm / 60.0), lmax=lmax_out)
        bandlim = hp.alm2map(hp.almxfl(alm_lo, bl_lo * pw), nside_out, lmax=lmax_out)
        print(f"{fwhm:6.1f} {np.std(code[fwhm] - ref):9.4f} {np.std(code[fwhm] - ref_pw):9.4f} "
              f"{np.std(std[fwhm] - ref):9.4f} {np.std(std[fwhm] - ref_pw):9.4f} "
              f"{np.std(coarse - bandlim):22.4f}", flush=True)
    print("# between-band differences (the part no common CMB amplitude absorbs): rms of band_i - band_j")
    print(f"{'pair':>14s} {'code':>9s} {'s.t.d.':>9s}")
    f = list(a.fwhm)
    for i in range(len(f)):
        for j in range(i + 1, len(f)):
            print(f"{f[i]:6.1f}-{f[j]:<6.1f} {np.std(code[f[i]] - code[f[j]]):9.4f} "
                  f"{np.std(std[f[i]] - std[f[j]]):9.4f}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--test", required=True, choices=["floor", "order"])
    p.add_argument("--param")
    p.add_argument("--cl")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--target", type=float, default=120.0)
    p.add_argument("--fwhm", type=float, nargs="+", default=[5.0, 13.2, 32.3, 52.8, 60.0])
    a = p.parse_args()
    {"floor": test_floor, "order": test_order}[a.test](a)
