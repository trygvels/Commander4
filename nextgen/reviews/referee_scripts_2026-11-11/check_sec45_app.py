"""Referee check of Sect. 4, Sect. 5 and the appendix of paper A1.

Small dense problems, written from the equations in the paper and nothing else.
  (A) tilt identity, Eq. (tilt): multi-component, dense noise, band-averaged spectra, two regions.
  (B) Hessian of (1/2) ln det A: equals 2T for one component with diagonal noise; what it is with a
      second (index-independent) component and with dense noise.
  (C) first-order shift H^-1 g against the exact flat-prior optimum, with H = J + 2T and with J.
  (D) Monte Carlo: E[s] = 0, Var(s) = J + T, E[curvature] = J, E[J(a_hat)] = J + T, sandwich.
  (E) zero-mean prior with a maximized scale, common noise pattern: ln L = q/2 - (n/2) ln q.
  (F) appendix: expected score when a Gaussian prior on an index-independent component has the
      wrong variance; sign and size of t_r.
A few seconds, well under 1 GB.
"""
import numpy as np
from scipy import optimize

rng = np.random.default_rng(1111)
nu = np.array([22.4, 28.1, 32.7, 40.3, 43.8, 60.0])   # GHz, illustrative
nb = len(nu)


def mix_sync(beta, nuref, width=0.0):
    """Band-averaged power law (top-hat band of fractional width), any reference."""
    if width == 0.0:
        return (nu / nuref) ** beta
    out = np.empty(nb)
    for i, n0 in enumerate(nu):
        x = np.linspace(n0 * (1 - width / 2), n0 * (1 + width / 2), 201)
        out[i] = np.mean((x / nuref) ** beta)
    return out


def build_M(betas, region, nuref, comps=("sync",), width=0.0):
    """M of shape (npix*nb, npix*ncomp); data ordered pixel-major."""
    npix = len(region)
    nc = len(comps)
    M = np.zeros((npix * nb, npix * nc))
    for p in range(npix):
        for c, name in enumerate(comps):
            if name == "sync":
                col = mix_sync(betas[region[p]], nuref, width)
            elif name == "cmb":
                col = np.ones(nb)
            elif name == "dust":
                col = (nu / 353.0) ** 1.5
            M[p * nb:(p + 1) * nb, p * nc + c] = col
    return M


def pieces(M, Ninv, d, Lam=None):
    A = M.T @ Ninv @ M
    if Lam is not None:
        A = A + Lam
    Ainv = np.linalg.inv(A)
    ahat = Ainv @ (M.T @ Ninv @ d)
    P = Ninv - Ninv @ M @ Ainv @ M.T @ Ninv
    chimin = d @ P @ d
    sign, logdet = np.linalg.slogdet(A)
    return A, Ainv, ahat, P, chimin, logdet


def random_dense_N(n, cond=30.0):
    Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
    ev = np.exp(rng.uniform(0, np.log(cond), n))
    N = (Q * ev) @ Q.T
    return N, np.linalg.inv(N)


# ---------------------------------------------------------------- (A)
print("(A) tilt identity, Eq. (tilt)")
npix = 12
region = np.array([0] * 7 + [1] * 5)
comps = ("sync", "cmb", "dust")
N, Ninv = random_dense_N(npix * nb)
d = rng.normal(size=npix * nb) * 3
for betas in [np.array([-3.1, -2.8]), np.array([-2.5, -3.6])]:
    out = []
    for nuref in (23.0, 30.0):
        M = build_M(betas, region, nuref, comps, width=0.2)
        A, Ainv, ahat, P, chimin, logdet = pieces(M, Ninv, d)
        out.append((chimin, -0.5 * chimin - 0.5 * logdet))
    pred = np.sum(np.bincount(region) * betas) * np.log(30.0 / 23.0)
    print(f"  betas={betas}: d chimin={out[1][0]-out[0][0]:+.2e}  d lnL_flat={out[1][1]-out[0][1]:+.10f}"
          f"  sum n_r beta_r ln(30/23)={pred:+.10f}")

# ---------------------------------------------------------------- (B)
print("(B) Hessian of (1/2) ln det A against 2T")


def JT(betas, region, nuref, comps, Ninv, a, h=1e-5, Lam=None, width=0.0):
    nr = len(betas)
    M = build_M(betas, region, nuref, comps, width)
    A = M.T @ Ninv @ M + (0 if Lam is None else Lam)
    Ainv = np.linalg.inv(A)
    P = Ninv - Ninv @ M @ Ainv @ M.T @ Ninv
    Mr = []
    for r in range(nr):
        e = np.zeros(nr); e[r] = h
        Mr.append((build_M(betas + e, region, nuref, comps, width) - build_M(betas - e, region, nuref, comps, width)) / (2 * h))
    J = np.array([[(Mr[r] @ a) @ P @ (Mr[s] @ a) for s in range(nr)] for r in range(nr)])
    T = np.array([[np.trace(Ainv @ Mr[r].T @ P @ Mr[s]) for s in range(nr)] for r in range(nr)])
    return J, T, M, Mr, P, Ainv


def half_logdet(betas, region, nuref, comps, Ninv, Lam=None):
    M = build_M(betas, region, nuref, comps)
    A = M.T @ Ninv @ M + (0 if Lam is None else Lam)
    return 0.5 * np.linalg.slogdet(A)[1]


def num_hess(f, x, h=1e-3):
    n = len(x); Hm = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            ei = np.zeros(n); ej = np.zeros(n); ei[i] = h; ej[j] = h
            Hm[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
    return Hm


def num_grad(f, x, h=1e-4):
    n = len(x); g = np.zeros(n)
    for i in range(n):
        e = np.zeros(n); e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g


npix = 40
region1 = np.zeros(npix, dtype=int)
sig = np.exp(rng.normal(0, 0.5, size=(npix, nb)))            # per-pixel, per-band noise, varying pattern
Ninv_diag = np.diag(1.0 / sig.ravel() ** 2)
b0 = np.array([-3.1])
for label, comps_, Ninv_ in [("one component, diagonal noise", ("sync",), Ninv_diag),
                             ("sync + free CMB, diagonal noise", ("sync", "cmb"), Ninv_diag),
                             ("one component, dense noise", ("sync",), random_dense_N(npix * nb)[1])]:
    a = rng.normal(size=npix * len(comps_)) * 5
    J, T, *_ = JT(b0, region1, 23.0, comps_, Ninv_, a)
    Hd = num_hess(lambda x: half_logdet(x, region1, 23.0, comps_, Ninv_), b0)
    print(f"  {label:34s}: Hess(1/2 ln det A)={Hd[0,0]:9.4f}  T={T[0,0]:9.4f}  ratio Hess/T={Hd[0,0]/T[0,0]:.4f}")

# ---------------------------------------------------------------- (C)
print("(C) first-order shift of the flat-prior optimum from the profile optimum")
npix = 400
region1 = np.zeros(npix, dtype=int)
w_common = 1.0 / np.array([1.0, 1.4, 1.2, 1.3, 2.0, 1.8]) ** 2
for snr_scale in (6.0, 2.0, 1.0):
    Ninv_ = np.diag(np.tile(w_common, npix))
    a = rng.standard_t(4, size=npix) * snr_scale
    Mtrue = build_M(b0, region1, 23.0, ("sync",))
    # noise-free data so that the comparison isolates the response formula (profile optimum = truth).
    # The sky is fixed (amplitude a at 23 GHz); J and T below are evaluated with the amplitude
    # expressed at the reference in use, a (nuref/23)^beta, so that they describe the same sky.
    dat = Mtrue @ a
    a23 = a.copy()
    for nuref in (23.0, 26.0, 30.0):
        a = a23 * (nuref / 23.0) ** b0[0]
        def nl_prof(x):
            return 0.5 * pieces(build_M(x, region1, nuref, ("sync",)), Ninv_, dat)[4]
        def nl_flat(x):
            pc = pieces(build_M(x, region1, nuref, ("sync",)), Ninv_, dat)
            return 0.5 * pc[4] + 0.5 * pc[5]
        bp = optimize.minimize_scalar(lambda x: nl_prof(np.array([x])), bounds=(-4.6, -1.6), method="bounded", options={"xatol": 1e-10}).x
        bf = optimize.minimize_scalar(lambda x: nl_flat(np.array([x])), bounds=(-4.6, -1.6), method="bounded", options={"xatol": 1e-10}).x
        J, T, *_ = JT(b0, region1, nuref, ("sync",), Ninv_, a)
        g = -num_grad(lambda x: half_logdet(x, region1, nuref, ("sync",), Ninv_), b0)[0]
        m2 = w_common * (nu / nuref) ** (2 * b0[0])
        nupiv = np.exp(np.sum(m2 * np.log(nu)) / np.sum(m2))
        print(f"  T/J={T[0,0]/J[0,0]:6.3f} nuref={nuref:.0f} nupiv={nupiv:.2f}: exact shift={bf-bp:+.5f}  g/(J+2T)={g/(J[0,0]+2*T[0,0]):+.5f}"
              f"  g/J={g/J[0,0]:+.5f}  n ln(nuref/nupiv)/(J+2T)={npix*np.log(nuref/nupiv)/(J[0,0]+2*T[0,0]):+.5f}")

# ---------------------------------------------------------------- (D)
print("(D) Monte Carlo of score, curvature and plug-in information (2 regions, sync+CMB+dust, dense noise)")
npix = 10
region2 = np.array([0] * 6 + [1] * 4)
comps = ("sync", "cmb", "dust")
N, Ninv = random_dense_N(npix * nb, cond=10.0)
L = np.linalg.cholesky(N)
bt = np.array([-3.1, -2.9])
a = rng.normal(size=npix * 3) * 1.5
J, T, M, Mr, P, Ainv = JT(bt, region2, 23.0, comps, Ninv, a)
h = 1e-4


def chimin_of(x, dat):
    return pieces(build_M(x, region2, 23.0, comps), Ninv, dat)[4]


Ms = {}
for r in range(2):
    for sgn in (+1, -1):
        e = np.zeros(2); e[r] = sgn * h
        Mx = build_M(bt + e, region2, 23.0, comps)
        Ax = np.linalg.inv(Mx.T @ Ninv @ Mx)
        Ms[(r, sgn)] = Ninv - Ninv @ Mx @ Ax @ Mx.T @ Ninv
Pr = [(Ms[(r, 1)] - Ms[(r, -1)]) / (2 * h) for r in range(2)]
Prr0 = (Ms[(0, 1)] - 2 * P + Ms[(0, -1)]) / h**2
nsim = 40000
S = np.zeros((nsim, 2)); C00 = np.zeros(nsim); Jhat = np.zeros((nsim, 2, 2))
G = Ainv @ M.T @ Ninv
MrP = [[Mr[r].T @ P @ Mr[s] for s in range(2)] for r in range(2)]
for i in range(nsim):
    dat = M @ a + L @ rng.normal(size=npix * nb)
    S[i] = [-0.5 * dat @ Pr[r] @ dat for r in range(2)]
    C00[i] = 0.5 * dat @ Prr0 @ dat
    ah = G @ dat
    Jhat[i] = [[ah @ MrP[r][s] @ ah for s in range(2)] for r in range(2)]
covS = np.cov(S.T)
print(f"  mean score = {S.mean(0)} +- {S.std(0)/np.sqrt(nsim)}")
print(f"  Var(s)/(J+T): [0,0] {covS[0,0]/(J+T)[0,0]:.4f}  [1,1] {covS[1,1]/(J+T)[1,1]:.4f}  [0,1] cov {covS[0,1]:.4f} vs {(J+T)[0,1]:.4f}")
print(f"  E[curvature of chimin/2]_00 = {C00.mean():.4f} +- {C00.std()/np.sqrt(nsim):.4f};  J_00 = {J[0,0]:.4f};  (J+T)_00 = {(J+T)[0,0]:.4f}")
print(f"  E[J(a_hat)]_00 = {Jhat[:,0,0].mean():.4f} +- {Jhat[:,0,0].std()/np.sqrt(nsim):.4f};  (J+T)_00 = {(J+T)[0,0]:.4f};  T/J_00 = {T[0,0]/J[0,0]:.3f}")
sand = np.linalg.inv(J) @ (J + T) @ np.linalg.inv(J)
lin = S @ np.linalg.inv(J).T          # first-order estimator error
print(f"  first-order estimator variance / sandwich: {np.var(lin[:,0])/sand[0,0]:.4f}, {np.var(lin[:,1])/sand[1,1]:.4f}")

# ---------------------------------------------------------------- (E)
print("(E) zero-mean prior, maximized scale, common noise pattern")
npix = 300
a = rng.standard_t(3, size=npix) * 1.2
m_true = (nu / 23.0) ** (-3.1)
dat = a[:, None] * m_true[None, :] + rng.normal(size=(npix, nb)) / np.sqrt(w_common)[None, :]


def q_of(beta):
    m = (nu / 23.0) ** beta
    A = np.sum(w_common * m * m)
    b = dat @ (w_common * m)
    return np.sum(b * b) / A, A, b


def lnL_tau_max(beta):
    q, A, b = q_of(beta)
    f = lambda lt: -(0.5 * q * (np.exp(lt) * A) / (1 + np.exp(lt) * A) - 0.5 * npix * np.log1p(np.exp(lt) * A))
    res = optimize.minimize_scalar(f, bounds=(-20, 20), method="bounded", options={"xatol": 1e-12})
    return -res.fun


bb = np.linspace(-3.4, -2.8, 13)
num = np.array([lnL_tau_max(x) for x in bb])
qq = np.array([q_of(x)[0] for x in bb])
form = 0.5 * qq - 0.5 * npix * np.log(qq)
dd = (num - form); dd -= dd.mean()
print(f"  q/n range {qq.min()/npix:.2f} to {qq.max()/npix:.2f}; max |numerical - (q/2 - n/2 ln q) - const| = {np.abs(dd).max():.2e}")
bq = optimize.minimize_scalar(lambda x: -q_of(x)[0], bounds=(-4, -2), method="bounded", options={"xatol": 1e-11}).x
bz = optimize.minimize_scalar(lambda x: -lnL_tau_max(x), bounds=(-4, -2), method="bounded", options={"xatol": 1e-11}).x
hh = 2e-3
cq = -(0.5 * q_of(bq + hh)[0] - q_of(bq)[0] + 0.5 * q_of(bq - hh)[0]) / hh**2
cz = -(lnL_tau_max(bz + hh) - 2 * lnL_tau_max(bz) + lnL_tau_max(bz - hh)) / hh**2
qopt = q_of(bq)[0]
u = a * a * np.sum(w_common * m_true**2)
print(f"  optimum: profile {bq:.6f}, zero-mean {bz:.6f}; curvature ratio {cz/cq:.4f}; 1-n/q = {1-npix/qopt:.4f}; "
      f"1/(1+T/J) at truth = {1/(1+npix/np.sum(u)):.4f}")

# ---------------------------------------------------------------- (F)
print("(F) appendix: expected score with a Gaussian prior of the wrong variance on the CMB")
npix = 30
region1 = np.zeros(npix, dtype=int)
comps = ("sync", "cmb")
sig = np.exp(rng.normal(0, 0.3, size=(npix, nb)))
Ninv_ = np.diag(1.0 / sig.ravel() ** 2)
sig_c = 0.7
Lam = np.zeros((npix * 2, npix * 2))
Lam[1::2, 1::2] = np.eye(npix) / sig_c**2
b0 = np.array([-3.1])
a_s = rng.normal(size=npix) * 3
hh = 1e-4


def Ptil(x):
    M = build_M(x, region1, 23.0, comps)
    A = M.T @ Ninv_ @ M + Lam
    return Ninv_ - Ninv_ @ M @ np.linalg.inv(A) @ M.T @ Ninv_, M, A


P0, M0, A0 = Ptil(b0)
Pp, _, _ = Ptil(b0 + hh)
Pm, _, _ = Ptil(b0 - hh)
Pr = (Pp - Pm) / (2 * hh)
Mr0 = (build_M(b0 + hh, region1, 23.0, comps) - build_M(b0 - hh, region1, 23.0, comps)) / (2 * hh)
A0inv = np.linalg.inv(A0)
t_paper = np.trace(Lam @ A0inv @ Mr0.T @ (A0inv @ M0.T @ Ninv_).T)
C = M0[:, 1::2]
for kappa2 in (1.0, 0.25, 4.0):
    # exact expectation of s = -1/2 d^T P_r d over noise and CMB with variance kappa2 * sig_c^2
    mean_d = M0[:, 0::2] @ a_s
    Cov_d = np.diag(sig.ravel() ** 2) + kappa2 * sig_c**2 * C @ C.T
    Es = -0.5 * (mean_d @ Pr @ mean_d + np.trace(Pr @ Cov_d))
    print(f"  kappa^2={kappa2:4.2f}: exact E[s]={Es:+.5f}   (kappa^2-1) t={(kappa2-1)*t_paper:+.5f}   -(kappa^2-1) t={-(kappa2-1)*t_paper:+.5f}")
