"""Scratch: configurations for the coverage runs."""
import time
import numpy as np
from rv_common import Problem, ref_info, RegionLikelihood
from lowres.fit import fit_profile
np.set_printoptions(precision=4, linewidth=170, suppress=True)
print("sync only, flat prior")
for (n, amp, slope) in [(320, 0.5, 2.5), (320, 0.35, 2.5), (320, 0.25, 2.5), (320, 0.6, 1.5), (320, 0.4, 1.5)]:
    p = Problem(n=n, nreg=2, dense=(0, 2, 5), seed=7, noise=1.0, amp=amp, noise_slope=slope)
    i = ref_info(p, p.beta_true, p.a_s, use_dust=False, use_cmb=False)
    J, T = i["J"], i["T"]
    Ji = np.linalg.inv(J)
    print(n, amp, slope, " J", np.diag(J), " T/J", np.diag(T) / np.diag(J), " sd_sw", np.sqrt(np.diag(Ji @ (J + T) @ Ji)),
          " sd_naive", np.sqrt(np.diag(Ji)))
print("full, cmb prior")
for (n, amp, slope, sc) in [(384, 1.2, 2.5, 0.3), (384, 0.8, 2.5, 0.3), (384, 0.6, 2.5, 0.3), (384, 0.8, 2.5, 1.0)]:
    p = Problem(n=n, nreg=2, dense=(0, 2, 5), seed=7, noise=1.0, amp=amp, amp_c=sc, noise_slope=slope)
    i = ref_info(p, p.beta_true, p.a_s, cmb_sd=sc)
    J, T = i["J"], i["T"]
    Ji = np.linalg.inv(J)
    print(n, amp, slope, sc, " J", np.diag(J), " T/J", np.diag(T) / np.diag(J), " sd_sw", np.sqrt(np.diag(Ji @ (J + T) @ Ji)),
          " sd_naive", np.sqrt(np.diag(Ji)))
    rng = np.random.default_rng(1)
    t0 = time.time()
    for _ in range(3):
        d = p.mean(a_c=sc * rng.standard_normal(n)) + p.noise(rng)
        L = RegionLikelihood(p.banddata(d), p.region, cmb_prior_sd=sc)
        bh, res = fit_profile(L, p.beta_true)
    print("   time per fit %.2f s, nfev %d, beta_hat" % ((time.time() - t0) / 3, res.nfev), bh)
