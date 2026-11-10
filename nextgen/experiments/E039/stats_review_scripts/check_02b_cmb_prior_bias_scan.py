"""Check 2b: closed-form size of the profile-score bias when the CMB carries a Gaussian prior.

Validated against Monte Carlo in check_02 (cases C, D, E). With x the true amplitudes, Lam the prior precision,
B = A^-1 M^T W:
    E[g_r | x] = -x^T M_r^T B^T Lam x + x^T Lam A^-1 M_r^T B^T Lam x - t_r,     t_r = tr(Lam A^-1 M_r^T B^T).
If the CMB sky is drawn from N(0, kappa^2 sigma_c^2 I) while the code assumes sigma_c:
    E[g_r] = (kappa^2 - 1) t_r          (zero only if kappa = 1; kappa = 0 is "no CMB in the data").
For one fixed sky drawn with kappa = 1 the conditional mean has zero mean over skies and an rms printed below.
The first-order bias of beta_hat is -J^-1 E[g]; in units of the sandwich sd it is about -E[g_r] / sqrt(J_rr + T_rr).

Toy: 8 delta bands (23, 30, 33, 41, 44, 61, 70, 353 GHz), noise 1 to 2.5 per element in every band below 100 GHz
(flat in frequency), 0.5 at 353, dense covariance in three bands, one region.
"""
import numpy as np

from rv_common import Problem, ref_info

np.set_printoptions(precision=4, linewidth=170, suppress=True)
NUS = (23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0)
print(f"{'n':>5s} {'sigma_c':>8s} {'T/J':>6s} {'sd_sw':>7s} {'t/sqrt(J+T)':>12s} {'bias/sd no CMB':>15s} {'bias/sd 2x sd':>14s} "
      f"{'rms cond. mean/sd, fixed matched sky':>38s} {'shrink':>7s}")
for n in (48, 192, 768):
    for sc in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 100.0):
        p = Problem(n=n, nreg=1, nus=NUS, dense=(0, 3, 7), seed=11, noise=1.0, noise353=0.5, amp=2.5, amp_d=3.0, amp_c=sc)
        i = ref_info(p, p.beta_true, p.a_s, cmb_sd=sc)
        J, T, Ai, B, Lam, Mr = i["J"][0, 0], i["T"][0, 0], i["Ai"], i["B"], np.diag(i["lam"]), i["dM"][0]
        G = Lam @ Ai @ Mr.T @ B.T
        t = np.trace(G)
        sd_score = np.sqrt(J + T)
        # fixed matched sky: g_bar(a_c) = -x^T M_r^T B^T Lam x + x^T G Lam x - t; linear part dominates the variance
        xs = np.concatenate([p.a_s, np.zeros(2 * n)])
        lin = -(Mr @ xs) @ B.T @ Lam                      # coefficient vector multiplying x (only the CMB block is nonzero)
        c = slice(2 * n, 3 * n)
        var_lin = sc ** 2 * np.sum(lin[c] ** 2)
        Gs = 0.5 * ((G @ Lam)[c, c] + (G @ Lam)[c, c].T)
        var_quad = 2 * sc ** 4 * np.sum(Gs * Gs)
        shrink = 1.0 - np.mean(np.diag(Ai)[c]) / sc ** 2      # fraction of the CMB variance the data determine
        print(f"{n:5d} {sc:8.2f} {T / J:6.3f} {np.sqrt((J + T)) / J:7.3f} {t / sd_score:+12.4f} {t / sd_score:+15.4f} {-3 * t / sd_score:+14.4f} "
              f"{np.sqrt(var_lin + var_quad) / sd_score:38.4f} {shrink:7.3f}")
    print()
print("bias/sd columns: first-order (beta_hat - truth) / sandwich sd = -E[g] / sqrt(J+T); 'no CMB' is kappa = 0, '2x sd' is kappa = 2.")
print("'shrink' = 1 - posterior variance / prior variance of a CMB element (0: prior dominated, 1: data dominated).")
