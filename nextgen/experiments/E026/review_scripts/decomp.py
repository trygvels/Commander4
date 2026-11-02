"""Independent numpy re-implementation of the collapsed hier target's data term, verified against
Commander4's RegionIndexLikelihood at one draw, then decompositions along the sync level."""
import os, sys, glob, time
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np, healpy as hp, h5py
D = "/Users/trygvels/data/e023_pysm_n1"
bands = [23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353]; freqs = np.array(bands, float)
maps = np.array([hp.read_map(f"{D}/band{b:03d}_map.fits", field=0) for b in bands])
rmsm = np.array([hp.read_map(f"{D}/band{b:03d}_rms.fits", field=0) for b in bands])
inv_var = 1.0 / rmsm**2
t = np.load(f"{D}/truth.npz"); T_fixed = float(t["T_fixed"])
sky = t["sky"].T  # (12, npix) noiseless?
print("rms maps constant per band:", np.allclose(rmsm.min(1), rmsm.max(1)), rmsm[:, 0])
print("map - sky std per band:", np.round((maps - sky).std(1), 1), "truth rms:", t["rms"])
print("map - sky mean per band (in units of se):", np.round((maps - sky).mean(1) / (t["rms"] / np.sqrt(49152)), 1))
g_s = np.maximum(hp.read_map(f"{D}/template_sync_exact_30GHz.fits", field=0), 1e-3)
g_d = np.maximum(hp.read_map(f"{D}/template_dust_exact_545GHz.fits", field=0), 1e-3)
print("template == a_sync?", np.allclose(g_s, np.maximum(t["a_sync"], 1e-3)), "dust:", np.allclose(g_d, np.maximum(t["a_dust"], 1e-3)))
NSIDE, NB = 64, 16
labels = np.arange(12 * NB * NB, dtype=float)
reg = np.rint(hp.ud_grade(labels, NSIDE)).astype(int)
H_K = 4.799243073366221e-11; T_CMB = 2.7255

from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust
from pixell.bunch import Bunch
def comp(cls, name, short, beta0, **sed):
    p = Bunch(polarization="I", shortname=short, spatially_varying_MM=False, Cl_prior_amplitude=None,
              lmax=191, units="uK_RJ", sample_spectral_index=True, spectral_index_nside=NB,
              spectral_index_prior=Bunch(type="gaussian", mean=beta0, rms=0.3),
              spectral_index_field=Bunch(type="gmrf", basis="graph", alpha=2.0, sigma_scale=0.3, rho_prior="loguniform"), beta=beta0, **sed)
    object.__setattr__(p, "_name", name)
    return cls(p, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True, eval_pol="I", comp_name=name)
cmb = comp(CMB, "CMB", "cmb", 0.0, nu_ref=1e-3); sync = comp(Synchrotron, "Synchrotron", "sync", -3.1, nu_ref=30.0)
dust = comp(ThermalDust, "ThermalDust", "dust", 1.6, T=T_fixed, nu_ref=545.0)
sed_cmb = np.asarray(cmb.get_sed(freqs), float)          # (12,)

def mixing(bs_reg, bd_reg):
    """(nreg, 12, 3) mixing per region."""
    M = np.empty((bs_reg.size, 12, 3))
    M[:, :, 0] = sed_cmb[None, :]
    M[:, :, 1] = sync.get_sed(freqs, beta=bs_reg)
    M[:, :, 2] = dust.get_sed(freqs, beta=bd_reg)
    return M

def data_term(bs_reg, bd_reg, mu_s, ls_s, mu_d, ls_d, d=maps, expected=False, parts=False, pix_sel=None):
    """sum_p [1/2 b'^T A'^-1 b' - 1/2 log det A'] + trailing hier terms, using per-pixel 3x3.
    expected: replace d by the noiseless sky and add 1/2 tr(A'^-1 M^T N^-1 M) (noise average)."""
    M = mixing(bs_reg, bd_reg)[reg]                       # (npix, 12, 3)
    P = np.zeros((M.shape[0], 3)); m = np.zeros((M.shape[0], 3))
    P[:, 0] = 1e-8
    P[:, 1] = 1.0 / (np.exp(ls_s) * g_s)**2; m[:, 1] = mu_s[reg] * g_s
    P[:, 2] = 1.0 / (np.exp(ls_d) * g_d)**2; m[:, 2] = mu_d[reg] * g_d
    W = inv_var.T                                          # (npix, 12)
    A = np.einsum("pbi,pb,pbj->pij", M, W, M)
    Ap = A + np.einsum("pi,ij->pij", P, np.eye(3))
    dd = (sky if expected else d).T
    b = np.einsum("pbi,pb,pb->pi", M, W, dd) + P * m
    Ainv = np.linalg.inv(Ap)
    quad = 0.5 * np.einsum("pi,pij,pj->p", b, Ainv, b)
    logdet = -0.5 * np.linalg.slogdet(Ap)[1]
    trail = -0.5 * np.sum(P[:, 1:] * m[:, 1:]**2, 1) + 0.5 * np.sum(np.log(P[:, 1:]), 1)
    tr = 0.5 * np.einsum("pij,pji->p", Ainv, A) if expected else np.zeros(quad.size)
    if pix_sel is not None:
        quad, logdet, trail, tr = quad[pix_sel], logdet[pix_sel], trail[pix_sel], tr[pix_sel]
    if parts:
        return quad.sum(), logdet.sum(), trail.sum(), tr.sum()
    return quad.sum() + logdet.sum() + trail.sum() + tr.sum()

S = "mcmc/sample_index_fields/state/"
def state_of(fn):
    st = {}
    with h5py.File(fn) as f:
        for c in ["Synchrotron", "ThermalDust"]:
            g = f[S + c]; fi = g["index_field"]
            st[c] = dict(z=fi["z"][()], mu=float(fi["mu"][()]), ls=float(fi["log_sigma"][()]), u=float(fi["u"][()]),
                         amu=g["amp_prior/mu"][()], als=float(g["amp_prior/log_s"][()]), beta=fi["beta_map"][()])
    return st
files = sorted(glob.glob(f"{D}/chains_gmrf16_exact_mumap_nc/chains_compsep/chain01_iter*.h5"))[500:]

# ---- verification against Commander4 at draw 0 (iter0501) ----
if "--verify" in sys.argv:
    from commander4.compsep import index_field as ifd, index_field_sampler as ifs
    cmb.amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=1e4)
    sync.amp_prior = ifd.AmpPrior("hier", template=g_s, mu_nside=16); dust.amp_prior = ifd.AmpPrior("hier", template=g_d, mu_nside=16)
    comps = [cmb, sync, dust]
    t0 = time.time(); specs, views = ifs.discover_index_specs(comps, None, centred=False); print(f"basis {time.time()-t0:.0f} s")
    names = ["Synchrotron", "ThermalDust"]
    layout = ifd.StateLayout(specs, names, ifs.hier_component_dims(comps, names))
    like = ifd.RegionIndexLikelihood.from_arrays(freqs, maps, inv_var, comps, layout, NSIDE, nthreads=1)
    fields = {s.comp_name: s.field for s in specs}
    st = state_of(files[0])
    bl = {c: fields[c].pack(st[c]["z"], st[c]["mu"], st[c]["ls"], st[c]["u"]) for c in names}
    x = layout.pack(bl, {c: (st[c]["amu"], st[c]["als"]) for c in names})
    lp, _ = like.logp_and_grad(x)
    pri = sum(fields[c].prior_logp_and_grad(bl[c])[0] for c in names)
    pri += sum(ifd.AmpPrior("hier", template=g_s, mu_nside=16, mu=st[c]["amu"], log_s=st[c]["als"]).log_hyperprior() for c in names)
    beta_s = fields["Synchrotron"].beta(bl["Synchrotron"]); beta_d = fields["ThermalDust"].beta(bl["ThermalDust"])
    print("beta from z equals stored beta_map:", np.abs(beta_s - st["Synchrotron"]["beta"]).max(), np.abs(beta_d - st["ThermalDust"]["beta"]).max())
    mine = data_term(beta_s, beta_d, st["Synchrotron"]["amu"], st["Synchrotron"]["als"], st["ThermalDust"]["amu"], st["ThermalDust"]["als"])
    print(f"Commander4 log target {lp:.2f} (script reported 1275513.32); minus priors {lp - pri:.2f}; mine {mine:.2f}; diff {lp - pri - mine:.4f}")
    print("prior terms:", pri)
    sys.exit(0)
