"""Addition 1: exact first-order bias of the profile index when the true CMB covariance is the prior covariance plus
a white part (the power that a pixel-averaged sky has above the band limit), nside-4 toy, reviewer's reference.

    E[g_r] = -tr(B_r U dC U^T),  B_r = Wt D_r [At^-1]_s. M^T Sigma^-1 (all built with the prior),  bias = -J^-1 E[g],
    and the true covariance of the score is the exact quadratic-form value under the true Sigma.

Two noise regimes: the toy of t1d (CMB rms / noise per element 0.6) and one with the ratio of the nside-16 maps (0.08),
with the synchrotron and dust amplitudes scaled by the same factor as the noise so that the index errors are equal.
"""
import sys

import numpy as np

sys.path.insert(0, "nextgen")
sys.path.insert(0, "nextgen/reviews/stats_review_scripts_2026-11-11b")
from lowres.window import lcdm_cmb_precision  # noqa: E402
import rv2_common as rv  # noqa: E402
import toy_n4 as tn  # noqa: E402

toy = tn.Toy()
prec, cov = lcdm_cmb_precision(toy.WM)
n = tn.N
U = np.tile(np.eye(n), (len(tn.FREQS), 1))
print("# t1e: exact first-order bias from CMB power missing in the prior (white, as a share of the prior variance per element)")
for label, fac in (("toy of t1d: noise 0.2 uK, CMB / noise 0.6", 1.0), ("noise and foregrounds times 8.5: noise 1.7 uK, CMB / noise 0.08", 8.5)):
    covs = [c * fac ** 2 for c in toy.covs]
    ref = rv.Reference(tn.FREQS, covs, toy.region, cmb_cov=cov)
    a_s, a_d = toy.a_s * fac, toy.a_d * fac
    M, Dm, At, Ati, Wt = ref.parts(tn.BETA_TRUE)
    J, T = ref.info(tn.BETA_TRUE, a_s)
    Ji = np.linalg.inv(J)
    sd = np.sqrt(np.diag(Ji @ (J + T) @ Ji))
    G = Ati[:n] @ M.T @ ref.Si
    m = M @ np.concatenate([a_s, a_d])
    Bs = []
    for r in range(ref.nreg):
        B = Wt @ (Dm * (toy.region == r)[None, :]) @ G
        Bs.append(B)
    print(f"## {label}: sandwich sd {sd.round(4)}, T/J {(np.diag(T) / np.diag(J)).round(4)}")
    for share in (0.13, 0.44, 1.0):
        dC = share * np.mean(np.diag(cov)) * np.eye(n)
        Sx = U @ dC @ U.T
        Eg = np.array([-np.trace(B @ Sx) for B in Bs])
        bias = -Ji @ Eg
        St = ref.Sigma + Sx
        Sy = [0.5 * (B + B.T) for B in Bs]
        SyS = [S_ @ St for S_ in Sy]
        K = np.array([[2 * np.sum(SyS[r] * SyS[q].T) + 4 * m @ Sy[r] @ St @ Sy[q] @ m for q in range(ref.nreg)] for r in range(ref.nreg)])
        sd_true = np.sqrt(np.diag(Ji @ K @ Ji))
        print(f"   missing white share {share:.2f}: bias {bias.round(5)} = " + " ".join(f"{x:+.3f}" for x in bias / sd) + " quoted sd; true sd / quoted sd "
              + " ".join(f"{x:.3f}" for x in sd_true / sd))
