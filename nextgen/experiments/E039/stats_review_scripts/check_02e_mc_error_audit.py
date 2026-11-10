"""Check 2e: audit of the Monte Carlo errors quoted in check_02 (reviewer's own check of the reviewer's numbers).

Two of about 36 standardized deviations in check_02 (three seeds) were beyond 3 sigma. Either chance, an
underestimated Monte Carlo error, or a real effect of a few tenths of a percent. Here: case B (CMB prior, CMB
redrawn from the prior), 200000 realizations in 20 batches; the scatter between batch means is compared with the
within-batch standard error, and the tails of the score are inspected. Arrays are saved.
"""
import sys

import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info

np.set_printoptions(precision=4, linewidth=170, suppress=True)
NREP, NB, SEED = 200000, 20, int(sys.argv[1]) if len(sys.argv) > 1 else 90210
SC = 2.0
prob = Problem(n=40, nreg=2, dense=(0, 2, 5), seed=7, noise=1.0, amp=1.6, amp_d=0.5, amp_c=SC)
beta0, n = prob.beta_true, prob.n
info = ref_info(prob, beta0, prob.a_s, cmb_sd=SC)
J, T = info["J"], info["T"]
rng = np.random.default_rng(SEED)
g = np.empty((NREP, 2))
Jh = np.empty((NREP, 2))
for it in range(NREP):
    d = prob.mean(a_c=SC * rng.standard_normal(n)) + prob.noise(rng)
    L = RegionLikelihood(prob.banddata(d), prob.region, cmb_prior_sd=SC)
    g[it] = L.profile(beta0, grad=True)[1]
    Jh[it] = np.diag(L.information(beta0)[0])
np.savez("/Users/trygvels/code/Commander4/nextgen/experiments/E039/stats_review_scripts/out_02e_arrays.npz", g=g, Jh=Jh, J=J, T=T)
K = np.diag(J + T)
print(f"case B, seed {SEED}, {NREP} realizations, {NB} batches of {NREP // NB}")
print("score mean            ", g.mean(0), "+-", g.std(0, ddof=1) / np.sqrt(NREP), " z", g.mean(0) / (g.std(0, ddof=1) / np.sqrt(NREP)))
v = (g - g.mean(0)) ** 2
print("score var / (J+T)     ", v.mean(0) / K, "+-", v.std(0, ddof=1) / np.sqrt(NREP) / K, " z", (v.mean(0) / K - 1) / (v.std(0, ddof=1) / np.sqrt(NREP) / K))
print("E[J(a_hat)] / (J+T)   ", Jh.mean(0) / K, "+-", Jh.std(0, ddof=1) / np.sqrt(NREP) / K, " z", (Jh.mean(0) / K - 1) / (Jh.std(0, ddof=1) / np.sqrt(NREP) / K))
gb = g.reshape(NB, -1, 2)
bm = gb.mean(1)
se_within = gb.std(1, ddof=1).mean(0) / np.sqrt(NREP // NB)
print("batch means: sd between batches", bm.std(0, ddof=1), " expected from within-batch sd", se_within,
      " ratio", bm.std(0, ddof=1) / se_within, "+-", 1 / np.sqrt(2 * (NB - 1)))
vb = (gb ** 2).mean(1)
se_v = (gb ** 2).std(1, ddof=1).mean(0) / np.sqrt(NREP // NB)
print("batch variances: sd between batches", vb.std(0, ddof=1), " expected", se_v, " ratio", vb.std(0, ddof=1) / se_v)
z = g / np.sqrt(K)
print("score tails: max |g|/sd", np.abs(z).max(0), " excess kurtosis", (z ** 4).mean(0) / (z ** 2).mean(0) ** 2 - 3,
      " skewness", (z ** 3).mean(0) / (z ** 2).mean(0) ** 1.5)
print("lag-1 autocorrelation of the score sequence", [float(np.corrcoef(g[:-1, r], g[1:, r])[0, 1]) for r in range(2)],
      " (expected 0 +- %.4f)" % (1 / np.sqrt(NREP)))
