"""Scratch: find configurations with T/J in 0.2 to 1 and a sandwich sd small enough for the linear regime."""
import numpy as np
from rv_common import Problem, ref_info
np.set_printoptions(precision=4, linewidth=170, suppress=True)
for (n, nreg, amp, noise, slope, sc, cmb_sd) in [(96, 2, 1.6, 1.0, 0.0, 0.3, 0.3), (96, 2, 3.0, 1.0, 2.0, 0.3, 0.3),
                                                 (96, 2, 2.0, 1.0, 2.0, 0.3, 0.3), (96, 2, 2.0, 1.0, 2.5, 0.3, 0.3),
                                                 (96, 2, 1.2, 1.0, 2.5, 0.3, 0.3), (96, 2, 1.2, 1.0, 2.5, 0.3, None),
                                                 (96, 2, 2.0, 1.0, 2.5, 0.3, None), (120, 2, 1.5, 1.0, 3.0, 0.3, 0.3)]:
    p = Problem(n=n, nreg=nreg, dense=(0, 2, 5), seed=7, noise=noise, amp=amp, amp_c=sc, noise_slope=slope)
    kw = {"cmb_sd": cmb_sd} if cmb_sd else {}
    i = ref_info(p, p.beta_true, p.a_s, **kw)
    J, T = i["J"], i["T"]
    Ji = np.linalg.inv(J)
    print(n, nreg, amp, noise, slope, cmb_sd, " J", np.diag(J), " T/J", np.diag(T) / np.diag(J), " sd_sw", np.sqrt(np.diag(Ji @ (J + T) @ Ji)),
          " sd_naive", np.sqrt(np.diag(Ji)), " cond(A) %.1e" % np.linalg.cond(i["A"]))
