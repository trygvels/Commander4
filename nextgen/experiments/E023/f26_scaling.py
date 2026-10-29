import numpy as np, sys
sys.argv = ["x"]
exec(open("/private/tmp/claude-501/-Users-trygvels-code-Commander4/e8b23887-d1a3-4d04-b2ac-228d7b739a81/scratchpad/f26_check.py").read().split("deltas = np.arange")[0])
# sigma_c: conditional noise sd of the per-pixel CMB amplitude given the other columns at their priors
mus, ss, mud, sd = THETA[100.]
M = np.stack([np.broadcast_to(cm, (49152, 12)), sync(bs0), dust(bd0, T)], -1)
for name, Tv in (("noise only", None), ("sync+dust at profiled prior", np.stack([(ss*a_s)**2, (sd*a_d)**2], -1))):
    if Tv is None:
        Cm = rms[:, None, :]**2*np.eye(12)[None]
    else:
        Cm = np.einsum("nij,nj,nkj->nik", M[:, :, 1:], Tv, M[:, :, 1:]) + rms[:, None, :]**2*np.eye(12)[None]
    sc = 1/np.sqrt(np.einsum("i,nij,j->n", cm, np.linalg.inv(Cm), cm))
    print(name, "sigma_c median per dust quartile:", [round(float(np.median(sc[qd == k])), 1) for k in range(4)],
          "leak fraction sigma_c^2/(100^2+sigma_c^2):", [round(float(np.median(sc[qd == k]**2/(1e4+sc[qd == k]**2))), 3) for k in range(4)])
flipped = sky - 2*cmb[:, None]*cm[None]
for ns, dl in ((1.0, 0.02), (0.3, 0.004), (0.1, 0.001)):
    deltas = np.arange(-5*dl, 5*dl+1e-9, dl)
    a = scan(f"d_Q2 sky noise {ns}", qd == 1, sky, ns*rms, 100., "d", deltas, THETA[100.])
    b = scan(f"d_Q2 flipped noise {ns}", qd == 1, flipped, ns*rms, 100., "d", deltas, THETA[100.])
    # linear (CMB-induced) slope at delta 0 from the antisymmetric part
    anti = (a - b)/2; slope = (anti[6] - anti[4])/(2*dl)
    sym = (a + b)/2; curv = (sym[6] - 2*sym[5] + sym[4])/dl**2
    print(f"   noise {ns}: CMB-induced slope {slope:+.1f} nats per unit delta, curvature {curv:.0f}, implied shift {-slope/curv:+.4f}")
