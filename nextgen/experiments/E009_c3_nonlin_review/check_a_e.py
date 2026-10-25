# Independent check of E012 points (a) and (e). Bands x components model d = M(beta) a + n.
import numpy as np
from scipy.integrate import nquad
rng = np.random.default_rng(1)
nu = np.array([23., 30., 44., 70., 100., 143.]); nu0 = 23.
def M(b): return np.column_stack([np.ones_like(nu), (nu/nu0)**b])          # C = 2: CMB + sync
def dM(b): return np.column_stack([np.zeros_like(nu), (nu/nu0)**b*np.log(nu/nu0)])
sig = np.array([3., 2., 2., 1.5, 1., 1.]); Ninv = np.diag(1/sig**2)
a_true = np.array([50., 30.]); beta_true = -3.0
d = M(beta_true) @ a_true + rng.normal(0, sig)
def parts(b):
    Mb = M(b); A = Mb.T @ Ninv @ Mb; bb = Mb.T @ Ninv @ d
    astar = np.linalg.solve(A, bb); chi2min = d @ Ninv @ d - bb @ astar
    return Mb, A, bb, astar, chi2min
def logL(b, w):
    _, A, _, _, c2 = parts(b); return -0.5*c2 - 0.5*w*np.linalg.slogdet(A)[1]
# (a) brute force: integrate exp(-chi2/2) over a (flat) and exp(-chi2/2)*sqrt(det A) (a-independent
# Jeffreys-for-a) and exp(-chi2/2)*sqrt(det I_joint(a,beta)) (full joint Jeffreys); compare beta
# dependence with w=1, w=0 formulas via log-ratio between two beta values.
def joint_jeffreys(a, b):
    Mb = M(b); J = np.column_stack([Mb, dM(b) @ a]); return np.sqrt(np.linalg.det(J.T @ Ninv @ J))
def integ(b, prior):
    Mb, A, bb, astar, c2 = parts(b)
    f = lambda a0, a1: np.exp(-0.5*((d - Mb@[a0, a1]) @ Ninv @ (d - Mb@[a0, a1]) - c2))*prior(np.array([a0, a1]), b)
    sd = np.sqrt(np.diag(np.linalg.inv(A))); r = [[astar[i]-8*sd[i], astar[i]+8*sd[i]] for i in range(2)]
    return np.log(nquad(f, r)[0]) - 0.5*c2
b1, b2 = -3.1, -2.9
for name, prior, w in [("flat -> w=1", lambda a, b: 1.0, 1),
                       ("sqrt(det A(beta)) -> w=0", lambda a, b: np.sqrt(np.linalg.det(M(b).T@Ninv@M(b))), 0),
                       ("joint Jeffreys sqrt(det I(a,beta))", joint_jeffreys, 0)]:
    print(f"(a) {name:38s} brute dlogL={integ(b2, prior)-integ(b1, prior):+.6f}  formula w={w}: {logL(b2, w)-logL(b1, w):+.6f}")
# (e) envelope theorem + trace term vs finite differences, matching kernel formulas
def grad_formula(b, w):
    Mb, A, bb, astar, c2 = parts(b); g = dM(b)[:, 1]; r = d - Mb @ astar
    return astar[1] * g @ Ninv @ r - w * np.linalg.solve(A, Mb.T @ Ninv @ g)[1]
for w in (0, 1):
    h = 1e-5; fd = (logL(-3.0+h, w) - logL(-3.0-h, w))/(2*h)
    print(f"(e) w={w}: analytic {grad_formula(-3.0, w):+.8f}  FD {fd:+.8f}  rel {abs(fd/grad_formula(-3.0,w)-1):.1e}")
# (e) chain rule through non-centred field with Matern spectrum and c(rho): ring graph Laplacian, N=12
N = 12; G = 2*np.eye(N) - np.roll(np.eye(N), 1, 1) - np.roll(np.eye(N), -1, 1)
lam, V = np.linalg.eigh(G); alpha = 2.0
def S(logrho):
    kap2 = 8.0/np.exp(logrho)**2; s = (kap2 + lam)**(-alpha/2); c = np.sqrt(np.mean(s**2)); return V @ np.diag(s/c)
def dS_dlogrho(logrho):
    kap2 = 8.0/np.exp(logrho)**2; s = (kap2 + lam)**(-alpha/2); ds = alpha*kap2*s/(kap2 + lam)
    c = np.sqrt(np.mean(s**2)); dc = np.mean(s*ds)/c; return V @ np.diag(ds/c - s*dc/c**2)
D = [rng.normal(0, 1, nu.size) + M(-3.0 + 0.2*rng.normal()) @ a_true for _ in range(N)]  # one data pixel per index pixel
def logpost(x, w=1):
    z, mu, ls, lr = x[:N], x[N], x[N+1], x[N+2]; beta = mu + np.exp(ls)*S(lr) @ z
    ll = 0.0
    for p in range(N):
        global d; d = D[p]; ll += logL(beta[p], w)
    return ll - 0.5*z@z - 0.5*(mu+3)**2/0.25 - 0.5*np.exp(ls)**2/0.09 + ls   # log-uniform rho: const
def gradpost(x, w=1):
    z, mu, ls, lr = x[:N], x[N], x[N+1], x[N+2]; sg = np.exp(ls); Sm = S(lr); beta = mu + sg*Sm @ z
    g = np.zeros(N)
    for p in range(N):
        global d; d = D[p]; g[p] = grad_formula(beta[p], w)
    return np.concatenate([sg*Sm.T@g - z, [g.sum() - (mu+3)/0.25, sg*(Sm@z)@g - sg**2/0.09 + 1, sg*(dS_dlogrho(lr)@z)@g]])
x0 = np.concatenate([rng.normal(0, 1, N), [-3.05, np.log(0.2), np.log(3.0)]])
an = gradpost(x0); fd = np.array([(logpost(x0+h*e) - logpost(x0-h*e))/(2*h) for e in np.eye(N+3)])
print("(e) chain rule max rel err over z,mu,log sigma,log rho:", np.max(np.abs(an-fd)/np.maximum(np.abs(fd),1e-8)))
print("    log rho component analytic/FD:", an[-1], fd[-1])
