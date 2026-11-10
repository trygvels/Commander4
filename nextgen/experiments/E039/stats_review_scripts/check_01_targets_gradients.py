"""Check 1: profile and marginal against brute force; analytic gradients against finite differences.

(a) profile and (1/2) ln det A against a QR solution of the whitened least-squares problem (no normal equations).
(b) marginal against the explicit Gaussian density of the data, -ln N(d; 0, N + M Sigma M^T), when every
    amplitude prior is proper (sync prior on, CMB prior on, no dust), absolute value including constants.
(c) flat priors as the limit of wide proper priors (differences between two beta, constants cancel).
(d) tiny problem: marginal by direct numerical quadrature over the amplitudes (no determinant identities).
(e) gradients by central finite differences, all prior configurations.
"""
import numpy as np
from scipy import integrate

from rv_common import Problem, RegionLikelihood, fd_grad, ref_solve

np.set_printoptions(precision=6, linewidth=160)
prob = Problem(n=48, nreg=3, dense=(0, 2, 5), seed=1)
rng = np.random.default_rng(5)
d = prob.mean() + prob.noise(rng)
bands = prob.banddata(d)
beta = np.array([-3.0, -3.25, -2.85])
lt = np.log(np.array([8.0, 15.0, 4.0]))
n = prob.n

print("== (a) profile and logdet vs QR reference (n=48, 6 bands, 3 dense) ==")
configs = [("dust+cmb flat", dict(), dict()),
           ("dust+cmb prior 0.3", dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3)),
           ("dust only", dict(use_cmb=False), dict(use_cmb=False)),
           ("sync only", dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False)),
           ("cmb only prior", dict(use_dust=False, cmb_prior_sd=0.3), dict(use_dust=False, cmb_sd=0.3))]
print(f"{'config':22s} {'profile(code)':>16s} {'rel.diff':>10s} {'marginal(code)':>16s} {'rel.diff':>10s} {'max|x-xref|':>12s}")
for name, kw, rkw in configs:
    L = RegionLikelihood(bands, prob.region, **kw)
    ref = ref_solve(prob, d, beta, **rkw)
    fp, fm = L.profile(beta), L.marginal(beta)
    sol = L.solve(beta)
    rp = 0.5 * ref["chi2"]
    rm = rp + 0.5 * ref["logdetA"]
    print(f"{name:22s} {fp:16.8f} {abs(fp - rp) / abs(rp):10.1e} {fm:16.8f} {abs(fm - rm) / abs(rm):10.1e} "
          f"{np.max(np.abs(sol['x'] - ref['x'])):12.1e}")
name = "sync prior + cmb prior"
L = RegionLikelihood(bands, prob.region, cmb_prior_sd=0.3, sync_prior=True)
ref = ref_solve(prob, d, beta, cmb_sd=0.3, log_tau=lt)
fm = L.marginal(beta, lt)
rm = 0.5 * ref["chi2"] + 0.5 * ref["logdetA"] + np.sum(lt[prob.region])
print(f"{name:22s} {'':16s} {'':10s} {fm:16.8f} {abs(fm - rm) / abs(rm):10.1e}")

print("\n== (b) proper priors on everything (sync prior, CMB prior, no dust): absolute -ln p(d) ==")
L = RegionLikelihood(bands, prob.region, use_dust=False, cmb_prior_sd=0.3, sync_prior=True)
N = prob.N()
_, ldN = np.linalg.slogdet(N)
dv = d.ravel()
for bt, ltt in ((beta, lt), (np.array([-3.4, -2.7, -3.1]), np.log(np.array([3.0, 30.0, 11.0])))):
    M = prob.M(bt, use_dust=False, use_cmb=True)
    Sig = np.concatenate([np.exp(2 * ltt)[prob.region], np.full(n, 0.3 ** 2)])
    C = N + (M * Sig[None, :]) @ M.T
    _, ldC = np.linalg.slogdet(C)
    brute = 0.5 * dv @ np.linalg.solve(C, dv) + 0.5 * ldC          # -ln p(d) without the (2 pi) constant
    code = L.marginal(bt, ltt) + 0.5 * ldN + n * np.log(0.3)        # constants the code drops: ln det N, n ln sigma_c
    print(f"brute {brute:.8f}   code+constants {code:.8f}   diff {code - brute:+.2e}")

print("\n== (c) flat priors as the wide-prior limit: f(beta1) - f(beta2), code flat vs explicit Gaussian density ==")
b1, b2 = beta, np.array([-3.3, -2.9, -3.05])
for name, kw, use_dust, use_cmb, cmb_sd in (("dust+cmb flat", dict(), True, True, None),
                                            ("dust flat, cmb prior", dict(cmb_prior_sd=0.3), True, True, 0.3),
                                            ("sync only flat", dict(use_dust=False, use_cmb=False), False, False, None)):
    L = RegionLikelihood(bands, prob.region, **kw)
    code = L.marginal(b1) - L.marginal(b2)
    out = []
    for big in (1e3, 1e4):
        vals = []
        for bt in (b1, b2):
            M = prob.M(bt, use_dust, use_cmb)
            sig = [np.full(n, big ** 2)]
            if use_dust:
                sig.append(np.full(n, big ** 2))
            if use_cmb:
                sig.append(np.full(n, cmb_sd ** 2 if cmb_sd else big ** 2))
            Sig = np.concatenate(sig)
            C = N + (M * Sig[None, :]) @ M.T
            _, ldC = np.linalg.slogdet(C)
            vals.append(0.5 * dv @ np.linalg.solve(C, dv) + 0.5 * ldC)
        out.append(vals[0] - vals[1])
    print(f"{name:22s} code {code:+.6f}   wide prior sd=1e3 {out[0]:+.6f}   sd=1e4 {out[1]:+.6f}")

print("\n== (d) direct quadrature, 2 elements, sync only, 4 bands, band 0 dense: -ln int exp(-Q) da ==")
p2 = Problem(n=2, nreg=2, nus=(22.8, 33.0, 44.1, 70.4), dense=(0,), seed=3, amp=3.0, random_regions=False)
d2 = p2.mean(a_d=np.zeros(2), a_c=np.zeros(2)) + p2.noise(np.random.default_rng(1))
bb = np.array([-3.0, -3.3])
W2 = p2.W()
for name, kw, ltt in (("flat", dict(), None), ("sync prior", dict(sync_prior=True), np.log(np.array([2.0, 5.0])))):
    L = RegionLikelihood(p2.banddata(d2), p2.region, use_dust=False, use_cmb=False, **kw)
    M = p2.M(bb, use_dust=False, use_cmb=False)
    lam = np.zeros(2) if ltt is None else np.exp(-2 * ltt)[p2.region]

    def integrand(a1, a0):
        r = d2.ravel() - M @ np.array([a0, a1])
        return np.exp(-0.5 * r @ W2 @ r - 0.5 * (lam[0] * a0 ** 2 + lam[1] * a1 ** 2))
    val, err = integrate.dblquad(integrand, -60, 60, -60, 60, epsabs=0, epsrel=1e-10)
    brute = -np.log(val) + np.log(2 * np.pi)                 # code drops (2 pi)^(dim/2), dim = 2
    if ltt is not None:
        brute += np.sum(ltt[p2.region])                      # prior normalization 1/tau per element (2 pi dropped)
    code = L.marginal(bb, ltt) if ltt is not None else L.marginal(bb)
    print(f"{name:12s} quadrature {brute:.8f}   code {code:.8f}   diff {code - brute:+.1e}  (quad rel err {err / val:.0e})")

print("\n== (e) gradients vs central finite differences (h=1e-5), max |analytic - FD| / max|analytic| ==")
for name, kw, _ in configs:
    L = RegionLikelihood(bands, prob.region, **kw)
    _, g = L.profile(beta, grad=True)
    e1 = np.max(np.abs(g - fd_grad(L.profile, beta))) / np.max(np.abs(g))
    _, g = L.marginal(beta, grad=True)
    e2 = np.max(np.abs(g - fd_grad(L.marginal, beta))) / np.max(np.abs(g))
    print(f"{name:22s} profile {e1:.1e}   marginal {e2:.1e}")
for name, kw in (("sync+cmb prior, dust", dict(cmb_prior_sd=0.3, sync_prior=True)),
                 ("sync prior only", dict(use_dust=False, use_cmb=False, sync_prior=True))):
    L = RegionLikelihood(bands, prob.region, **kw)
    _, gb, gt = L.marginal(beta, lt, grad=True)
    e1 = np.max(np.abs(gb - fd_grad(lambda b: L.marginal(b, lt), beta))) / np.max(np.abs(gb))
    e2 = np.max(np.abs(gt - fd_grad(lambda t: L.marginal(beta, t), lt))) / np.max(np.abs(gt))
    print(f"{name:22s} d/dbeta {e1:.1e}   d/dlog_tau {e2:.1e}")

print("\n== (f) profile() and information() with sync_prior=True (no log_tau argument on profile) ==")
L = RegionLikelihood(bands, prob.region, sync_prior=True)
for nm, fn in (("profile(beta)", lambda: L.profile(beta)), ("information(beta)", lambda: L.information(beta))):
    try:
        fn()
        print(nm, "ran")
    except Exception as e:  # noqa: BLE001
        print(nm, "raises", type(e).__name__, str(e)[:80])
