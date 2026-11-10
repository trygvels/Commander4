"""Check 2: Monte Carlo of the profile score, its covariance, the expected Hessian and E[J(a_hat)],
all at the true beta, against J and T from (i) the code and (ii) the reviewer's explicit matrices.

Cases
  A  flat priors on synchrotron, dust and CMB; truth fixed; noise redrawn.
  B  Gaussian CMB prior sd sc; CMB redrawn from N(0, sc^2 I) for every realization (prior is the truth).
  C  Gaussian CMB prior sd sc; ONE fixed CMB sky with rms sc; noise redrawn (the real-data situation).
  D  Gaussian CMB prior sd sc; CMB redrawn with sd 2 sc (prior variance wrong by 4).
  E  Gaussian CMB prior sd sc; ONE fixed CMB sky with rms 0 (no CMB in the data at all).
"""
import sys
import time

import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info

np.set_printoptions(precision=4, linewidth=170, suppress=True)
NREP = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
NHESS = NREP // 4
SC = 2.0
prob = Problem(n=40, nreg=2, dense=(0, 2, 5), seed=7, noise=1.0, amp=1.6, amp_d=0.5, amp_c=SC)
beta0 = prob.beta_true
n, nreg = prob.n, prob.nreg


def run(case, rng):
    kw = {} if case == "A" else {"cmb_prior_sd": SC}
    rkw = {} if case == "A" else {"cmb_sd": SC}
    info = ref_info(prob, beta0, prob.a_s, **rkw)
    a_c_fixed = prob.a_c * (SC / prob.a_c.std())             # fixed sky with rms exactly SC
    g_all, Jhat_sum, H_sum, H_n, H_sq = [], np.zeros((nreg, nreg)), np.zeros((nreg, nreg)), 0, np.zeros(nreg)
    Jhat_all = []
    T_code = None
    for it in range(NREP):
        if case in ("A", "C"):
            a_c = a_c_fixed
        elif case == "B":
            a_c = SC * rng.standard_normal(n)
        elif case == "D":
            a_c = 2 * SC * rng.standard_normal(n)
        else:
            a_c = np.zeros(n)
        d = prob.mean(a_c=a_c) + prob.noise(rng)
        L = RegionLikelihood(prob.banddata(d), prob.region, **kw)
        _, g = L.profile(beta0, grad=True)
        g_all.append(g)
        Jh, Th, _ = L.information(beta0)
        Jhat_all.append(np.diag(Jh))
        Jhat_sum += Jh
        T_code = Th
        if it < NHESS:
            h = 1e-4
            H = np.zeros((nreg, nreg))
            for r in range(nreg):
                e = np.zeros(nreg); e[r] = h
                H[r] = (L.profile(beta0 + e, grad=True)[1] - L.profile(beta0 - e, grad=True)[1]) / (2 * h)
            H_sum += 0.5 * (H + H.T)
            H_sq += np.diag(H) ** 2
            H_n += 1
    g_all, Jhat_all = np.array(g_all), np.array(Jhat_all)
    gm, gse = g_all.mean(0), g_all.std(0, ddof=1) / np.sqrt(NREP)
    gc = g_all - gm
    K = gc.T @ gc / (NREP - 1)
    Kse = np.array([[np.std(gc[:, i] * gc[:, j], ddof=1) / np.sqrt(NREP) for j in range(nreg)] for i in range(nreg)])
    # J "at the truth" as the code would give it: noiseless data, fitted amplitudes
    Lc = RegionLikelihood(prob.banddata(prob.mean(a_c=a_c_fixed if case in ("A", "C") else np.zeros(n))), prob.region, **kw)
    Jc, Tc, solc = Lc.information(beta0)
    J, T = info["J"], info["T"]
    print(f"\n--- case {case}: {NREP} realizations ---")
    print("T/J (diag, reviewer matrices)      ", np.diag(T) / np.diag(J))
    print("max |T_code - T_ref| / max|T_ref|  ", np.max(np.abs(T_code - T)) / np.max(np.abs(T)))
    print("J_ref diag                         ", np.diag(J), "  J_code(noiseless data) diag", np.diag(Jc))
    print("score mean                         ", gm, " +- ", gse, "  -> z =", gm / gse)
    print("score sd expected sqrt(diag(J+T))  ", np.sqrt(np.diag(J + T)))
    print("score cov diag / (J+T) diag        ", np.diag(K) / np.diag(J + T), " +- ", np.diag(Kse) / np.diag(J + T))
    print("score cov off-diag, MC vs J+T      ", K[0, 1], "+-", Kse[0, 1], " vs ", (J + T)[0, 1])
    Hm = H_sum / H_n
    Hse = np.sqrt(np.maximum(H_sq / H_n - np.diag(Hm) ** 2, 0.0) / H_n)
    print(f"E[Hessian] diag / J diag ({H_n} reps)", np.diag(Hm) / np.diag(J), " +- ", Hse / np.diag(J), "   [J+T would give", np.diag(J + T) / np.diag(J), "]")
    Jm, Jse = Jhat_all.mean(0), Jhat_all.std(0, ddof=1) / np.sqrt(NREP)
    print("E[J(a_hat)] diag / (J+T) diag      ", Jm / np.diag(J + T), " +- ", Jse / np.diag(J + T))
    print("E[J(a_hat) - T] diag / J diag      ", (Jm - np.diag(T)) / np.diag(J), " +- ", Jse / np.diag(J))
    return info, gm, gse


t0 = time.time()
CASES = tuple(sys.argv[2]) if len(sys.argv) > 2 else ("A", "B", "C", "D", "E")
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 2026
rng = np.random.default_rng(SEED)
print(f"seed {SEED}, cases {CASES}")
for case in CASES:
    info, gm, gse = run(case, rng)
    if case in ("C", "E"):
        # reviewer's closed form for the conditional mean of the score given a fixed amplitude vector x:
        #   E[g_r | x] = -x^T M_r^T B^T Lam x + x^T Lam A^-1 M_r^T B^T Lam x - tr(Lam A^-1 M_r^T B^T)
        a_c = prob.a_c * (SC / prob.a_c.std()) if case == "C" else np.zeros(n)
        x = np.concatenate([prob.a_s, prob.a_d, a_c])
        Lam, Ai, B = np.diag(info["lam"]), info["Ai"], info["B"]
        pred = np.array([-x @ Mr.T @ B.T @ Lam @ x + x @ Lam @ Ai @ Mr.T @ B.T @ Lam @ x - np.trace(Lam @ Ai @ Mr.T @ B.T)
                         for Mr in info["dM"]])
        print("conditional score mean, closed form", pred, "  (MC above)", "  in units of score sd:",
              pred / np.sqrt(np.diag(info["J"] + info["T"])))
        print("implied first-order bias of beta_hat, -J^-1 E[g]:", -np.linalg.solve(info["J"], pred),
              "  sandwich sd:", np.sqrt(np.diag(np.linalg.inv(info["J"]) @ (info["J"] + info["T"]) @ np.linalg.inv(info["J"]))))
print(f"\nelapsed {time.time() - t0:.0f} s")
