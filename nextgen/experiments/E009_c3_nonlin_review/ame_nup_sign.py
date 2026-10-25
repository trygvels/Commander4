"""Toy check of Commander3's per-pixel 'marginal' lnL (comp_lnL_max_chisq_diagonal, use_det=.true.).

Model per pixel p, band b: d_pb = M_b(nu_p) a_p + n_pb, flat prior on a_p, Gaussian prior on nu_p.
Exact flat-prior marginal:   lnL = -chi2_min/2 - 0.5 ln A_p,  A_p = sum_b M_b^2 / sigma_pb^2
Ridge (profile):             lnL = -chi2_min/2
Commander3 code ('marginal'): lnL = -chi2_min/2 - 0.5 ln(1/A_p) = -chi2_min/2 + 0.5 ln A_p
Setup loosely follows BeyondPlanck's AME nu_p step: one full-sky nu_p, intensity, bands 15-100 GHz,
nside 16 low-resolution pixels in an AME-bright mask, prior nu_p ~ N(21, 3) GHz.
"""
import numpy as np
import healpy as hp

rng = np.random.default_rng(1)
nside = 16
lat = 90 - np.degrees(hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))[0])
mask = np.abs(lat) < 25
npx = mask.sum()
bands = np.array([22.8, 28.4, 33.0, 40.7, 44.1, 60.8, 70.4, 93.5])
nu_ref, width, nup_true = 22.8, 0.6, 22.0
prior_mu, prior_sd = 21.0, 3.0

def sed(nup):
    s = lambda nu: np.exp(-0.5 * (np.log(nu / nup) / width) ** 2) / nu**2   # flux -> RJ
    return s(bands) / s(nu_ref)

amp = 300.0 * np.exp(-np.abs(lat[mask]) / 10.0) * rng.lognormal(0, 0.5, npx) + 10.0
grid = np.linspace(17, 30, 2601)
Mg = np.array([sed(g) for g in grid])                       # (ngrid, nband)

def posts(d, sig):
    w = 1 / sig**2                                          # (npx, nband)
    A = (Mg[:, None, :] ** 2 * w[None]).sum(-1)             # (ngrid, npx)
    b = (Mg[:, None, :] * (d * w)[None]).sum(-1)
    ridge = 0.5 * (b**2 / A).sum(1)
    vol = 0.5 * np.log(A).sum(1)
    lp = -0.5 * ((grid - prior_mu) / prior_sd) ** 2
    out = {}
    for name, ll in [("ridge", ridge), ("marginal (correct)", ridge - vol), ("C3 code 'marginal'", ridge + vol)]:
        ll = ll + lp
        p = np.exp(ll - ll.max()); p /= p.sum()
        m = (grid * p).sum(); s = np.sqrt(((grid - m) ** 2 * p).sum())
        out[name] = (m, s)
    return out

print(f"pixels in mask: {npx}, truth nu_p = {nup_true} GHz")
for noise in [2.0, 20.0, 100.0]:
    res = {}
    for r in range(20):
        sig = np.full((npx, bands.size), noise)
        d = amp[:, None] * sed(nup_true)[None] + sig * rng.standard_normal(sig.shape)
        for k, (m, s) in posts(d, sig).items():
            res.setdefault(k, []).append(((m - nup_true) / s, m, s))
    print(f"\nnoise {noise:6.1f} uK per band per pixel (20 realizations)")
    for k, v in res.items():
        v = np.array(v)
        print(f"  {k:22s} mean {v[:,1].mean():7.3f}  sd {v[:,2].mean():.4f}  bias/sd {v[:,0].mean():+7.2f}")
