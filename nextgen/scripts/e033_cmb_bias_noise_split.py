#!/usr/bin/env python
"""E033: split the flat-prior GLS CMB residual at given index fields into its foreground-bias and
noise terms.

The estimator is linear in the data. Per pixel c_hat = w(theta)^T d with

    w^T = e_c^T (M^T N^-1 M)^-1 M^T N^-1,

M the (bands, 3) mixing matrix of CMB, power-law synchrotron and modified-blackbody dust at the
assumed (beta_s, beta_d, T_d), and N = diag(rms^2). With d = sky + n the residual is exactly

    c_hat - c = [w^T sky - c] + w^T n,

a foreground-bias term that vanishes at the true per-pixel SEDs, and a noise term with variance
sum_b w_b^2 rms_b^2 = [(M^T N^-1 M)^-1]_cc, which depends on the assumed SEDs and not on the sky.

usage: e033_cmb_bias_noise_split.py --dir DIR --refs REFS.npz --t-map TMAP.fits
           [--tags TAG ...] [--tags-fixed-t TAG ...] [--burn 500]
"""
import argparse, glob, sys
ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True); ap.add_argument("--refs", required=True)
ap.add_argument("--t-map", required=True, help="per-region dust T map (FITS) of the oracle-T chains")
ap.add_argument("--tags", nargs="*", default=[], help="chains run with the T map")
ap.add_argument("--tags-fixed-t", nargs="*", default=[], help="chains run with the scalar T_fixed")
ap.add_argument("--burn", type=int, default=500)
a = ap.parse_args(); sys.argv = [sys.argv[0]]
import numpy as np, healpy as hp, h5py
sys.path.insert(0, "nextgen/scripts")
from e020_analyze import cmb_sed, gls_cmb, region_map, sed

truth = np.load(f"{a.dir}/truth.npz", allow_pickle=True)
nside = int(truth["nside"]); nu = truth["nu"]; rms = truth["rms"]; npix = hp.nside2npix(nside)
data = np.array([hp.read_map(f"{a.dir}/band{int(f):03d}_map.fits") for f in nu]).T
sky = truth["sky"]; noise = data - sky; cmb = truth["cmb"]
lat = hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]; hi = np.abs(lat) > 20
reg = region_map(16, nside); r = np.load(a.refs); qs = r["qs16"][reg]
Tmap = hp.read_map(a.t_map, field=0)[reg]; Tfix = float(truth["T_fixed"])
ones = np.ones(npix)


def weights(beta_s, beta_d, T, dust_col=None):
    """(npix, B) CMB row of (M^T N^-1 M)^-1 M^T N^-1 at the assumed SEDs."""
    dcol = sed("dust", nu, beta_d, T) if dust_col is None else dust_col
    M = np.stack([np.broadcast_to(cmb_sed(nu), dcol.shape), sed("sync", nu, beta_s, T), dcol], axis=2)
    Ni = 1.0 / rms ** 2
    A = np.einsum("pbi,b,pbj->pij", M, Ni, M)
    return np.linalg.solve(A, np.transpose(M, (0, 2, 1)) * Ni[None, None, :])[:, 0, :]


def rmsf(x, m=None):
    return float(np.sqrt(np.mean((x if m is None else x[m]) ** 2)))


def row(label, beta_s, beta_d, T, dust_col=None, check=False):
    w = weights(beta_s, beta_d, T, dust_col)
    bias = np.einsum("pb,pb->p", w, sky) - cmb
    nz = np.einsum("pb,pb->p", w, noise)
    var = np.einsum("pb,b->p", w ** 2, rms ** 2)
    tot = bias + nz
    if check:   # the weights reproduce gls_cmb, the estimator behind every earlier number
        ref = gls_cmb(data, rms, nu, beta_s, beta_d, T, dust_col) - cmb
        print(f"   [check] max |weights estimate - gls_cmb| = {np.max(np.abs(ref - tot)):.2e} uK")
    q = " ".join(f"{rmsf(bias, qs == k):7.2f}" for k in range(4))
    clos = rmsf(tot) - np.sqrt(rmsf(bias) ** 2 + rmsf(nz) ** 2)
    print(f"{label:58s} {rmsf(bias):7.2f} {rmsf(bias, hi):7.2f}  {q}  {np.sqrt(var.mean()):7.2f} {rmsf(nz):7.2f} "
          f"{np.sqrt(var[hi].mean()):7.2f}  {rmsf(tot):7.2f} {rmsf(tot, hi):7.2f} {clos:+6.2f}  {np.mean(bias):+6.2f}")
    ROWS.append((label, bias, nz))
    return bias


def spectra():
    """Exploratory block (not pre-registered): band powers of the bias and of the realized noise term,
    sqrt(mean D_l) in uK with D_l = l(l+1)C_l/2pi, full sky and at |b|>20 (pseudo-C_l over f_sky, no mode
    decoupling). The bias does not depend on the overall noise scale and the noise term scales with it, so the
    last columns give the bias-to-noise power ratio at this noise scale and at 0.03 times it."""
    bins = ((2, 10), (11, 30), (31, 100), (101, 191)); ell = np.arange(3 * nside)
    fac = ell * (ell + 1) / (2 * np.pi)

    def bp(m, mask):
        cl = hp.anafast(np.where(mask, m, 0.0), lmax=3 * nside - 1) / mask.mean()
        return np.array([np.mean((fac * cl)[lo:hi + 1]) for lo, hi in bins])

    print("\nband powers sqrt(mean D_l) [uK] in l = 2-10, 11-30, 31-100, 101-191: bias, then bias/noise power ratio "
          "at this noise scale, then the same ratio at 0.03 times this noise")
    for name, mask in (("full sky", np.ones(npix, bool)), ("|b|>20", hi)):
        print(f"-- {name}")
        for label, bias, nz in ROWS:
            b, n = bp(bias, mask), bp(nz, mask)
            print(f"{label:58s} " + " ".join(f"{v:7.3f}" for v in np.sqrt(b)) + "   "
                  + " ".join(f"{v:9.2e}" for v in b / n) + "   " + " ".join(f"{v:9.2e}" for v in b / n / 0.03 ** 2))


ROWS = []


def pm(tag, comp):
    files = sorted(glob.glob(f"{a.dir}/chains_{tag}/chains_compsep/chain01_iter*.h5"))[a.burn:]
    if not files:
        return None
    return np.mean([h5py.File(f, "r")[f"mcmc/sample_index_fields/state/{comp}/index_field/beta_map"][()]
                    for f in files], 0)


def regmean(x):
    return (np.bincount(reg, weights=x) / np.bincount(reg))[reg]


bs, bd, Td = truth["beta_s"], truth["beta_d"], truth["T_d"]
print(f"E033 split of the GLS CMB residual, {a.dir}, noise scale {float(truth['noise_scale'])}, rms in uK_CMB")
print("columns: bias = w^T sky - cmb (all sky, |b|>20, then per sync S/N quartile Q1 bright to Q4 faint); noise = "
      "analytic sqrt(mean var), realized rms of w^T n, analytic at |b|>20; total = bias + noise (all sky, |b|>20); "
      "closure = total - sqrt(bias^2 + realized noise^2); mean bias")
print(f"{'fields':58s} {'bias':>7s} {'b>20':>7s}  {'Q1':>7s} {'Q2':>7s} {'Q3':>7s} {'Q4':>7s}  {'n_an':>7s} {'n_real':>7s} "
      f"{'n_an>20':>7s}  {'total':>7s} {'tot>20':>7s} {'clos':>6s}  {'mean':>6s}")
row("true per-pixel beta_s and true dust band column (control)", bs, bd, Td, dust_col=truth["dust_bands"], check=True)
row("true per-pixel beta_s, beta_d, T_d (parametric control)", bs, bd, Td)
row("true per-pixel indices, region-mean T map", bs, bd, Tmap)
row("truth nside-16 region means, T map", regmean(bs), regmean(bd), Tmap)
row("constant at truth means, T map", bs.mean() * ones, bd.mean() * ones, Tmap)
row("constant at truth means, T fixed", bs.mean() * ones, bd.mean() * ones, Tfix)
row("constant at prior means (-3.1, 1.6), T map", -3.1 * ones, 1.6 * ones, Tmap)
for d in (-0.30, -0.164, +0.10, +0.30):
    row(f"truth region means, sync shifted {d:+.3f}, T map", regmean(bs) + d, regmean(bd), Tmap)
row("region-fit references bs_trueT, bd_trueT, T map", r["bs_trueT"][reg], r["bd_trueT"][reg], Tmap)
row("region-fit references bs_fixT, bd_fixT, T fixed", r["bs_fixT"][reg], r["bd_fixT"][reg], Tfix)
for tags, T, tl in ((a.tags, Tmap, "T map"), (a.tags_fixed_t, Tfix, "T fixed")):
    for tag in tags:
        S, Dd = pm(tag, "Synchrotron"), pm(tag, "ThermalDust")
        if S is None:
            print(f"{tag}: no post-burn files"); continue
        print(f"-- {tag} ({tl}): field means {S.mean():+.4f} {Dd.mean():+.4f}")
        row(f"   chain sync, chain dust", S[reg], Dd[reg], T)
        row(f"   chain sync, true dust column", S[reg], bd, Td, dust_col=truth["dust_bands"])
        row(f"   true sync, chain dust", bs, Dd[reg], T)
        row(f"   constant at the chain's field means", S.mean() * ones, Dd.mean() * ones, T)
spectra()
