"""E016 L1 checks of the hierarchical amplitude prior (model.hier_loglike, posterior.HierPosterior).

1. Direct form: per pixel, F_p equals -1/2 r^T Sigma^-1 r - 1/2 log det Sigma + const with
   Sigma = N + M T M^T, r = d - M m, the constant being -1/2 log det N + 1/2 d^T N^-1 d
   (the core drops -1/2 d^T N^-1 d; both forms drop -(B/2) log 2 pi); checked over (beta, mu, s)
   including s = 1e3 and negative mu (reviewer's check, reproduced here).
2. Gradients of the full log posterior (K + 7 coordinates) against 4-point finite differences at
   s = 1e3, 0.3 and 1e-3.
3. Limit (a): s = 1e3 reproduces w1 with the cmb prior 1/70^2 added to A up to a beta-independent
   constant (difference of differences across beta, and beta gradients).
4. Limit (b): s = 0.5, mu = 1, g = m reproduces mode "gauss" with inv_T = 1/(0.5 m)^2 (no floor)
   once the two trailing terms are subtracted.
Run from nextgen/prototypes:  python -m indexfield.checks.check_2a_hier
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from indexfield import coverage, field, posterior, simulate               # noqa: E402
from indexfield.model import (PixelData, collapsed_loglike, hier_loglike, hier_prior,  # noqa
                              mixing, CMB_SD)

OUT = os.path.join(HERE, "out", "check_2a_hier.json")
results = {}
all_ok = True


def report(name, value, thresh):
    global all_ok
    ok = bool(value < thresh)
    all_ok &= ok
    results[name] = dict(value=float(value), thresh=float(thresh), ok=ok)
    print(f"  {'PASS' if ok else 'FAIL'} {name:48s} {value:.3e} < {thresh:.0e}")


def fd4(f, x, h):
    g = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x); e[i] = h
        g[i] = (-f(x + 2 * e) + 8 * f(x + e) - 8 * f(x - e) + f(x - 2 * e)) / (12 * h)
    return g


rng = np.random.default_rng(20261026)
nside_d, nside_beta = 8, 1
basis = field.GraphLaplacianBasis(nside_beta)
fld = field.MaternField(basis)
hyper = posterior.HyperPrior(rho_prior="loguniform").resolved(basis.h)
data, truth = coverage.simulate_template(nside_d, nside_beta, fld, hyper, rng, 1.0, None,
                                         "default")
g = simulate.templates(nside_d)
N = data.d.shape[0]
print(f"check_2a_hier: nside_d {nside_d}, nside_beta {nside_beta}, N {N}, K {fld.K}")

# 1. direct form, per pixel
print("1. direct form")
worst = 0.0
for p in rng.choice(N, 6, replace=False):
    one = PixelData(nu=data.nu, d=data.d[p:p + 1], inv_var=data.inv_var[p:p + 1],
                    region=np.zeros(1, int), n_beta=1, beta_d=data.beta_d, T_d=data.T_d)
    gp = g[p:p + 1]
    consts = []
    for beta, mu, ls in [(-3.0, (1.0, 1.0), (np.log(0.3),) * 2), (-2.5, (-0.3, 1.4), (0.0, -2.0)),
                         (-3.4, (0.8, 0.9), (np.log(1e3),) * 2), (-2.8, (1.2, 0.5), (-4.0, 1.0))]:
        F = hier_loglike(np.array([beta]), one, gp, mu, ls)[0]
        m, inv_T = hier_prior(gp, np.asarray(mu), np.asarray(ls))
        M, _ = mixing(data.nu, np.array([beta]), beta_d=data.beta_d, T_d=data.T_d)
        M = M[0]
        Nn = np.diag(1.0 / one.inv_var[0])
        Sigma = Nn + M @ np.diag(1.0 / inv_T[0]) @ M.T
        r = one.d[0] - M @ m[0]
        direct = -0.5 * r @ np.linalg.solve(Sigma, r) - 0.5 * np.linalg.slogdet(Sigma)[1]
        const = -0.5 * np.linalg.slogdet(Nn)[1] + 0.5 * one.d[0] @ (one.inv_var[0] * one.d[0])
        consts.append(F - direct - const)
    worst = max(worst, np.ptp(consts) / max(1.0, abs(np.mean(consts))))
report("direct form: residual constant across (beta, mu, s), rel", worst, 1e-8)

# 2. gradients of the full posterior
print("2. gradients (fd4)")
post = posterior.HierPosterior(fld, hyper, data, g)
x0 = post.random_init(rng)
for s_val in (1e3, 0.3, 1e-3):
    x = x0.copy(); x[post.K + 5:] = np.log(s_val)
    lp, grad = post.logp_and_grad(x)
    f = lambda z: post.logp_and_grad(z)[0]       # noqa: E731
    # The log posterior is of order N mu^2 / s^2 (1e9 at s = 1e-3), so the rounding error of a
    # difference quotient scales as 1e-16 |lp| / h: take the best of three steps per coordinate.
    rel = np.full(post.dim, np.inf)
    for h in (1e-3, 1e-4, 1e-5):
        num = fd4(f, x, h)
        rel = np.minimum(rel, np.abs(num - grad) / np.maximum(np.abs(grad), 1e-3 * np.max(np.abs(grad))))
    # At s = 1e-3 (5.7 hyperprior sd below the median, prior sd 0.1 percent of the template) A' is
    # dominated by T^-1 and cond(A') ~ 1e8, so x = A'^-1 b' carries ~1e-8 relative error, which
    # the field-block gradient (through x) inherits at the 1e-5 level: threshold 1e-4 there.
    thr = 1e-4 if s_val < 1e-2 else 1e-5
    print(f"     worst coordinate {int(np.argmax(rel))} of {post.dim} (field block 0..{post.K - 1})")
    report(f"grad vs fd4 at s = {s_val:g} (|lp| {abs(lp):.1e}), max rel", rel.max(), thr)
    report(f"grad vs fd4 at s = {s_val:g}, amp block max rel", rel[post.K + 3:].max(), 1e-5)

# 3. limit (a): s -> infinity vs w1 with cmb prior
print("3. limit (a)")
m0 = np.zeros_like(g); inv0 = np.zeros_like(g); inv0[:, 0] = 1.0 / CMB_SD ** 2
nb = data.n_beta
w1c = PixelData(nu=data.nu, d=data.d, inv_var=data.inv_var, region=data.region, n_beta=nb,
                beta_d=data.beta_d, T_d=data.T_d, m=m0, inv_T=inv0)
betas = (-3.3, -3.0, -2.7)
ptp_w1 = np.ptp([collapsed_loglike(np.full(nb, b), w1c, "gauss")[0] for b in betas])
# The residual T^-1 = 1/(s g)^2 relative to A must shrink as 1/s^2: report s = 1e3, 1e4, 1e5
# and test at 1e5 (the s = 1e3 value is the "1 to 2 percent of one diagonal entry" of the review).
for s_val in (1e3, 1e4, 1e5):
    diffs, gdiff = [], 0.0
    for beta in betas:
        bv = np.full(nb, beta)
        Fh, gh, _, _ = hier_loglike(bv, data, g, (0.7, 1.3), (np.log(s_val),) * 2)
        Fw, gw = collapsed_loglike(bv, w1c, "gauss")
        diffs.append(Fh - Fw); gdiff = max(gdiff, np.max(np.abs(gh - gw) / np.abs(gw)))
    thr = 1e-6 if s_val >= 1e5 else 1.0
    report(f"s = {s_val:g} vs w1 + cmb prior: ptp(diff) / ptp(w1)", np.ptp(diffs) / ptp_w1, thr)
    report(f"s = {s_val:g} vs w1 + cmb prior: beta gradient rel", gdiff, 10 * thr)

# 4. limit (b): s = 0.5, mu = 1 vs gauss with inv_T = 1/(0.5 m)^2, trailing terms removed
print("4. limit (b)")
m1, invT1 = hier_prior(g, np.array([1.0, 1.0]), np.log([0.5, 0.5]))
gauss = PixelData(nu=data.nu, d=data.d, inv_var=data.inv_var, region=data.region, n_beta=nb,
                  beta_d=data.beta_d, T_d=data.T_d, m=m1, inv_T=invT1)
worst = 0.0
for beta in (-3.3, -3.0):
    bv = np.full(nb, beta)
    Fh, gh, _, _ = hier_loglike(bv, data, g, (1.0, 1.0), np.log([0.5, 0.5]))
    trailing = -0.5 * np.sum(m1 * m1 * invT1) + 0.5 * np.sum(np.log(invT1))
    Fg, gg = collapsed_loglike(bv, gauss, "gauss")
    worst = max(worst, abs(Fh - trailing - Fg) / abs(Fg), np.max(np.abs(gh - gg) / np.abs(gg)))
report("s = 0.5, mu = 1 vs gauss (no floor), rel", worst, 1e-10)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(dict(all_ok=all_ok, results=results), open(OUT, "w"), indent=1)
print("ALL PASS" if all_ok else "SOME FAIL", "->", OUT)
sys.exit(0 if all_ok else 1)
