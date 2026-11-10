"""Check 8: does fit_profile depend on the start value at low signal-to-noise?
All tests in the repository (and checks 3 above) start the optimizer at the truth. Here: three starts per data set."""
import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info
from lowres.fit import fit_profile

np.set_printoptions(precision=4, linewidth=170, suppress=True)
for (label, pk, lk, rk, sync_only) in [
        ("sync only, n=320, T/J~0.9", dict(n=320, nreg=2, amp=0.35, noise_slope=2.5), dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False), True),
        ("sync only, n=96, T/J~0.3-0.5", dict(n=96, nreg=2, amp=0.7, noise_slope=2.5), dict(use_dust=False, use_cmb=False), dict(use_dust=False, use_cmb=False), True),
        ("3 comp, cmb prior, n=384, T/J~0.46", dict(n=384, nreg=2, amp=0.6, noise_slope=2.5), dict(cmb_prior_sd=0.3), dict(cmb_sd=0.3), False)]:
    p = Problem(dense=(0, 2, 5), seed=7, noise=1.0, amp_c=0.3, **pk)
    info = ref_info(p, p.beta_true, p.a_s, **rk)
    Ji = np.linalg.inv(info["J"])
    sd = np.sqrt(np.diag(Ji @ (info["J"] + info["T"]) @ Ji))
    rng = np.random.default_rng(17)
    zero = np.zeros(p.n)
    nrep, ndiff, maxdiff, nb, worse = 150, 0, 0.0, 0, 0
    for _ in range(nrep):
        mean = p.mean(a_d=zero, a_c=zero) if sync_only else p.mean(a_c=0.3 * rng.standard_normal(p.n))
        L = RegionLikelihood(p.banddata(mean + p.noise(rng)), p.region, **lk)
        fits = [fit_profile(L, np.full(2, s)) for s in (-3.1, -2.0, -4.2)]
        b = np.array([f[0] for f in fits])
        fv = np.array([f[1].fun for f in fits])
        dmax = np.max(np.abs(b - b[0]) / sd)
        maxdiff = max(maxdiff, dmax)
        ndiff += int(dmax > 0.02)
        worse += int(np.any(fv[1:] > fv[0] + 1e-3))
        nb += int(np.any((b < -4.6 + 1e-6) | (b > -1.6 - 1e-6)))
    print(f"{label}: sandwich sd {sd}; {nrep} data sets, starts at -3.1 (truth), -2.0, -4.2: optima differing by more than 0.02 sd: {ndiff}; "
          f"largest difference {maxdiff:.3f} sd; a far start ending at a higher objective (by > 1e-3): {worse}; any fit at a bound: {nb}")
