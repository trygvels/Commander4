"""Toy check of Sky._marginal_terms: one region, 20 pixels; the formula against a 2-d quadrature
over (mu_s, mu_d) of exp(chi(mu) + log prior) using the non-marginal perregmu law; and the
expected mode against the average of realization-mode marginals over noise draws."""
import sys, numpy as np
sys.argv = ["x"]; sys.path.insert(0, "nextgen/scripts")
import e023_exact_scan as E
rng = np.random.default_rng(0); n = 20
gs = 10**(1 + 0.5*rng.random(n)); gd = 10**(1.5 + 0.5*rng.random(n))
bs = -3.0 + 0.1*rng.standard_normal(n); bd = 1.55 + 0.05*rng.standard_normal(n); T = np.full(n, 20.0)
Ms, Md = E.sync_sed(bs), E.dust_sed(bd, T)
a = np.stack([30*rng.standard_normal(n), gs*np.exp(0.05*rng.standard_normal(n)), gd*np.exp(0.05*rng.standard_normal(n))], 1)
M = np.stack([E.cmb_mixing(E.FREQS)[None, :].repeat(n, 0), Ms, Md], -1)
sky = np.einsum("nij,nj->ni", M, a); rms = np.full((n, 12), 3.0) * (0.5 + rng.random((n, 12)))
d = sky + rms*rng.standard_normal((n, 12))
S = E.Sky(d, sky, rms); S.gs, S.gd = gs, gd
reg = np.zeros(n, int); q = np.zeros(n, int)
hm = E.Hyper("perregmarg", q, q, reg=reg); hp_ = E.Hyper("perregmu", q, q, reg=reg)
ls = np.log([0.05]*4); th_m = np.concatenate([ls, ls])
for expected in (False, True):
    chi_m, ld_m, tr_m, _ = S.terms(Ms, Md, th_m, hm, expected)
    # quadrature over (mu_s, mu_d) with the explicit-mu law (one region): chi(mu) + prior
    def chi_of(mus, mud):
        th = np.concatenate([[mus], ls, [mud], ls]); c, ld, tr, _ = S.terms(Ms, Md, th, hp_, expected); return c + tr
    g = np.linspace(0.6, 1.4, 161); dg = g[1]-g[0]
    vals = np.array([[chi_of(x, y) - 0.5*((x-1)/0.5)**2 - 0.5*((y-1)/0.5)**2 for y in g] for x in g])
    m = vals.max(); quad = m + np.log(np.sum(np.exp(vals - m)) * dg*dg) - np.log(2*np.pi*0.25)
    if expected:
        # the expected-mode quadrature above has E[chi(mu)] inside the integral (Jensen): compare instead
        # the marginal formula with the average over noise draws of realization-mode marginals
        vals_r = []
        for k in range(400):
            S2 = E.Sky(sky + rms*rng.standard_normal((n, 12)), sky, rms); S2.gs, S2.gd = gs, gd
            vals_r.append(S2.terms(Ms, Md, th_m, hm, False)[0])
        print(f"expected mode: formula chi+tr {chi_m + tr_m:.4f}; mean over 400 noise draws of realization marginals {np.mean(vals_r):.4f} +- {np.std(vals_r)/20:.4f}")
    else:
        print(f"realization mode: formula chi {chi_m:.6f}; quadrature {quad:.6f}; difference {chi_m - quad:.2e}")
