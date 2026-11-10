"""Check 2f: closed-form prediction of the bias of the profile optimum for the configurations F5 and F6 of check_03
(CMB prior sd 0.3; truth sd 0.6, or no CMB). E[g_r] = (kappa^2 - 1) t_r, first-order bias = -J^-1 E[g].
Also t_r written with the quantities that RegionLikelihood.information already forms (C and A^-1), to show how the
code under review could report it:  t_r = sigma_c^-2 sum_{p in r} sum_k (C A^-1)[p, cmb k] (A^-1)[cmb k, p]."""
import numpy as np

from rv_common import Problem, RegionLikelihood, ref_info

np.set_printoptions(precision=4, linewidth=170, suppress=True)
p = Problem(dense=(0, 2, 5), seed=7, noise=1.0, amp_c=1.0, n=384, nreg=2, amp=2.5, amp_d=3.0, noise353=0.5,
            nus=(23.0, 30.0, 33.0, 41.0, 44.0, 61.0, 70.0, 353.0))
sc = 0.3
i = ref_info(p, p.beta_true, p.a_s, cmb_sd=sc)
J, T, Ai, B, Lam = i["J"], i["T"], i["Ai"], i["B"], np.diag(i["lam"])
t = np.array([np.trace(Lam @ Ai @ Mr.T @ B.T) for Mr in i["dM"]])
Ji = np.linalg.inv(J)
sd = np.sqrt(np.diag(Ji @ (J + T) @ Ji))
print("T/J", np.diag(T) / np.diag(J), " sandwich sd", sd, " t_r", t)
for label, k2 in (("F5: truth sd 0.6 (kappa = 2)", 4.0), ("F6: no CMB (kappa = 0)", 0.0), ("F7: matched (kappa = 1)", 1.0)):
    bias = -Ji @ ((k2 - 1.0) * t)
    print(f"{label}: predicted first-order bias {bias}, in sandwich sd {bias / sd}")

# the same t_r from the arrays the code under review forms in information()
n = p.n
L = RegionLikelihood(p.banddata(p.mean(a_c=np.zeros(n))), p.region, cmb_prior_sd=sc)
sol = L.solve(p.beta_true, need_inverse=True)
Ainv = sol["Ainv"]
C = np.zeros((n, 3 * n))
C[:, :n] = sum(b.scaled(sol["ds"][k], sol["s"][k]) for k, b in enumerate(L.bands))
one = np.ones(n)
for (j, cj) in L._nuisance_columns():
    C[:, j * n:(j + 1) * n] = sum(b.scaled(sol["ds"][k], one * getattr(b, cj)) for k, b in enumerate(L.bands))
CA = C @ Ainv
c = slice(L.ic * n, (L.ic + 1) * n)
t_code = np.bincount(p.region, weights=np.einsum("pk,kp->p", CA[:, c], Ainv[c, :n]) / sc ** 2, minlength=2)
print("t_r from the code's own arrays", t_code, " rel diff to explicit matrices", t_code / t - 1)
