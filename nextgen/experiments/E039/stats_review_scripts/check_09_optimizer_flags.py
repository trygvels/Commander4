"""Check 9: what the 'not converged' flags of fit_profile mean (3 to 7 of 800 fits in the three-component runs)."""
import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info
from lowres.fit import fit_profile

np.set_printoptions(precision=5, linewidth=170, suppress=True)
p = Problem(dense=(0, 2, 5), seed=7, noise=1.0, amp_c=1.0, n=384, nreg=2, amp=0.6, noise_slope=2.5)
info = ref_info(p, p.beta_true, p.a_s, cmb_sd=0.3)
Ji = np.linalg.inv(info["J"])
sd = np.sqrt(np.diag(Ji @ (info["J"] + info["T"]) @ Ji))
rng = np.random.default_rng(99)
nfail = 0
for it in range(400):
    L = RegionLikelihood(p.banddata(p.mean(a_c=0.3 * rng.standard_normal(p.n)) + p.noise(rng)), p.region, cmb_prior_sd=0.3)
    bh, res = fit_profile(L, p.beta_true)
    if not res.success:
        nfail += 1
        _, g = L.profile(bh, grad=True)
        step = np.linalg.solve(info["J"], g)          # Newton step that remains, with the expected curvature
        b2, res2 = fit_profile(L, bh)                 # restart from the returned point
        print(f"fit {it}: message '{res.message}', nit {res.nit}, nfev {res.nfev}; remaining Newton step / sd {step / sd}; "
              f"restart moves by {(b2 - bh) / sd} sd, restart success {res2.success}")
print(f"{nfail} of 400 fits returned success=False")
