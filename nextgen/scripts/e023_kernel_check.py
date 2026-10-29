#!/usr/bin/env python
"""E023 code control: the numpy marginal of e023_exact_scan.py against Commander4's own collapsed
hier likelihood (`marginal_loglike_regions` plus the hier trailing terms, as RegionIndexLikelihood
evaluates it), at three index states and two hyperparameter settings; realization mode,
data-built templates, fixed T, global layout. Differences between states must agree to 1e-8
relative; the absolute offset is the state-independent constant -1/2 d^T N^-1 d - 1/2 log det N.
"""
import argparse, sys
ap = argparse.ArgumentParser(); ap.add_argument("--dir", required=True); a = ap.parse_args()
sys.argv = [sys.argv[0]]
import numpy as np, healpy as hp
from pixell.bunch import Bunch
from commander4.compsep import index_field as ifd
from commander4.compsep import index_field_sampler as ifs
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust
from commander4.backend import compsep as cpp
sys.path.insert(0, "nextgen/scripts")
import e023_exact_scan as E

D = a.dir; NSIDE, NB = 64, 16
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; freqs = np.array(bands, float)
maps = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands])
rms = np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands])
inv_var = 1.0 / rms**2
t = np.load(f"{D}/truth.npz"); T_fixed = float(t["T_fixed"])

def comp(cls, name, short, **sed):
    p = Bunch(polarization="I", shortname=short, spatially_varying_MM=False, Cl_prior_amplitude=None,
              lmax=191, units="uK_RJ", sample_spectral_index=True, spectral_index_nside=NB,
              spectral_index_prior=Bunch(type="gaussian", mean=sed.get("beta", 0.0), rms=0.3),
              spectral_index_field=Bunch(type="gmrf", basis="graph", alpha=2.0, sigma_scale=0.3,
                                         rho_prior="loguniform"), **sed)
    object.__setattr__(p, "_name", name)
    return cls(p, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True, eval_pol="I", comp_name=name)

cmb = comp(CMB, "CMB", "cmb", nu_ref=1e-3)
sync = comp(Synchrotron, "Synchrotron", "sync", beta=-3.1, nu_ref=30.0)
dust = comp(ThermalDust, "ThermalDust", "dust", beta=1.6, T=T_fixed, nu_ref=545.0)
comps = [cmb, sync, dust]
cmb.amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=100.0)
gmaps = {}
for c, tpl, key in ((sync, "template_sync_30GHz.fits", "sync"), (dust, "template_dust_545GHz.fits", "dust")):
    g = np.maximum(hp.read_map(f"{D}/{tpl}", field=0), 1e-3); gmaps[key] = g
    c.amp_prior = ifd.AmpPrior("hier", template=g, mu=1.0, log_s=float(np.log(0.3)))
fs = ifs.field_prior_from_params(sync.comp_params, NB, True); fd = ifs.field_prior_from_params(dust.comp_params, NB, True)
layout = ifd.StateLayout([ifd.IndexSpec("Synchrotron", fs, -3.1, 0.3), ifd.IndexSpec("ThermalDust", fd, 1.6, 0.3)], ["Synchrotron", "ThermalDust"])
like = ifd.RegionIndexLikelihood.from_arrays(freqs, maps, inv_var, comps, layout, NSIDE, nthreads=2)
reg = hp.ud_grade(np.arange(hp.nside2npix(NB), dtype=float), NSIDE).astype(int)
assert np.array_equal(like.region, reg), "region map differs"
print("Commander4 likelihood built; pixels", like.region.size, "regions", NB**2 * 12)

def kernel_total(bs16, bd16, hyp_s, hyp_d):
    sync.amp_prior.mu, sync.amp_prior.log_s = hyp_s; dust.amp_prior.mu, dust.amp_prior.log_s = hyp_d
    fr, sky, w = like.streams[0]
    M = np.empty((reg.max() + 1, 12, 3)); dM = np.zeros_like(M)
    M[:, :, 0] = cmb.get_sed(fr)[None, :]; M[:, :, 1] = sync.get_sed(fr, beta=bs16); M[:, :, 2] = dust.get_sed(fr, beta=bd16)
    sky_aug, w_aug = ifd.augment_stream(sky, w, comps); M_aug, dM_aug = ifd.augment_mixing(M, dM)
    out = cpp.marginal_loglike_regions(sky_aug, w_aug, like.region, M_aug, dM_aug, 2, True, True, True)
    ms, Ps = sync.amp_prior.pseudo_band(reg.size); md, Pd = dust.amp_prior.pseudo_band(reg.size)
    trail = -0.5 * (ms**2 * Ps + md**2 * Pd) + 0.5 * (np.log(Ps) + np.log(Pd) + np.log(1 / 100.0**2))
    return float(np.sum(out[0]) + np.sum(trail))

d = maps.T; r = rms.T
S = E.Sky(d, d, r); S.gs, S.gd = gmaps["sync"], gmaps["dust"]
H = E.Hyper("global", np.zeros(reg.size, int), np.zeros(reg.size, int))
Tf = np.full(reg.size, T_fixed)
def numpy_total(bs16, bd16, hyp_s, hyp_d):
    th = np.array([hyp_s[0], hyp_s[1], hyp_d[0], hyp_d[1]])
    chi, ld, tr, _ = S.terms(E.sync_sed(bs16[reg]), E.dust_sed(bd16[reg], Tf), th, H, expected=False)
    return float(chi + ld)
const = float(-0.5 * np.sum(d**2 / r**2) - 0.5 * np.sum(np.log(r**2)))

bs0 = hp.ud_grade(np.asarray(t["beta_s"], float), NB); bd0 = hp.ud_grade(np.asarray(t["beta_eff3"], float), NB)
fs16 = np.asarray(t["fisher_sd_sync_16_profiled"]); fd16 = np.asarray(t["fisher_sd_dust_16_profiled"])
qs = np.argsort(fs16)[-768:]; qd = np.argsort(fd16)[:768]
h0 = ((1.0, float(np.log(0.3))), (1.0, float(np.log(0.3)))); h1 = ((1.1, float(np.log(0.5))), (0.9, float(np.log(0.27))))
bsB = bs0.copy(); bsB[qs] += 0.5; bdC = bd0.copy(); bdC[qd] += 0.05
states = [("A truth, hyp0", bs0, bd0, h0), ("B sync faint +0.5, hyp1", bsB, bd0, h1), ("C dust bright +0.05, hyp0", bs0, bdC, h0),
          ("D truth, hyp1", bs0, bd0, h1)]
vals = []
for lab, bs, bd, (hs, hd) in states:
    k = kernel_total(bs, bd, hs, hd); n = numpy_total(bs, bd, hs, hd)
    vals.append((lab, k, n)); print(f"{lab:28s} kernel+trail {k:18.6f}  numpy {n:18.6f}  numpy-kernel {n-k:16.6f}  (const {const:.6f}, numpy-kernel-const {n-k-const:.3e})")
k0, n0 = vals[0][1], vals[0][2]
ok = True
for lab, k, n in vals[1:]:
    dk, dn = k - k0, n - n0; rel = abs(dn - dk) / max(abs(dk), 1.0)
    ok &= rel < 1e-8
    print(f"difference to A: {lab:28s} kernel {dk:16.6f} numpy {dn:16.6f} relative mismatch {rel:.2e}")
print("KERNEL CHECK", "PASS" if ok else "FAIL")
