import numpy as np, healpy as hp
nside = 8; N = hp.nside2npix(nside); h = np.sqrt(4*np.pi/N)       # pixel spacing in rad
W = np.zeros((N, N))
for p in range(N):
    for q in hp.get_all_neighbours(nside, p, nest=True):
        if q >= 0: W[p, q] = 1
W = np.maximum(W, W.T); G = np.diag(W.sum(1)) - W
lam, V = np.linalg.eigh(G)
print("(c) nside", nside, "h[deg]", np.degrees(h), "degrees per pixel:", np.unique(W.sum(1)))
for l in (1, 2, 3, 4, 6, 8):
    i0 = l*l; blk = lam[i0:i0+2*l+1]
    print(f"    l={l}: lambda mean {blk.mean():.4f} spread {blk.std():.4f}; ratio lambda/(h^2 l(l+1)) = {blk.mean()/(h*h*l*(l+1)):.3f}")
print("    lambda_max", lam.max(), " (continuum 3h^2 l(l+1) at l=3nside-1:", 3*h*h*23*24, ")")
# (b) K_ii uniformity and zero-mode share for the normalized S at several rho (rho in pixel units)
alpha = 2
for rho_pix in (1, 2, 5, 20, 180/np.degrees(h)):
    kap2 = 8.0/rho_pix**2; s = (kap2 + lam)**(-alpha/2); c = np.sqrt(np.mean(s**2)); S = V*(s/c)
    Kii = (S**2).sum(1)
    print(f"(b) rho={rho_pix:6.1f} pix: mean K_ii {Kii.mean():.3f} min/max {Kii.min():.3f}/{Kii.max():.3f}; zero-mode share of variance {(s[0]/c)**2/N:.3f}")
# w=1 SBC: proper prior N(0, tau^-1) vs flat, beta-dependence of the correction and induced beta bias
nu = np.array([23,30,33,41,44,61,70,94,100,143,217,353.]); nu0 = 23.; tau = 1e-8
sig = np.full(nu.size, 10.0); Ninv = np.diag(1/sig**2)
def Mb(b): return np.column_stack([np.ones_like(nu), (nu/nu0)**b, (nu/353.)**1.55*np.exp(1e9*353*6.626e-34/(1.38e-23*19.6))/ (np.exp(1e9*nu*6.626e-34/(1.38e-23*19.6))-1)*(np.exp(1e9*353*6.626e-34/(1.38e-23*19.6))-1)/np.exp(1e9*353*6.626e-34/(1.38e-23*19.6))])
rng = np.random.default_rng(0)
for amp in (np.array([70., 300., 300.]), rng.normal(0, 1e4, 3)):
    d = Mb(-3.0) @ amp + rng.normal(0, sig)
    def corr(b):
        A = Mb(b).T@Ninv@Mb(b); bb = Mb(b).T@Ninv@d
        return 0.5*bb@np.linalg.solve(A+tau*np.eye(3), bb) - 0.5*np.linalg.slogdet(A+tau*np.eye(3))[1] - (0.5*bb@np.linalg.solve(A, bb) - 0.5*np.linalg.slogdet(A)[1])
    def ll1(b):
        A = Mb(b).T@Ninv@Mb(b); bb = Mb(b).T@Ninv@d; return 0.5*bb@np.linalg.solve(A, bb) - 0.5*np.linalg.slogdet(A)[1]
    hh = 1e-4; curv = -(ll1(-3+hh) - 2*ll1(-3) + ll1(-3-hh))/hh**2; dcorr = (corr(-3+hh)-corr(-3-hh))/(2*hh)
    print(f"(f) amps {np.round(amp)}: correction {corr(-3.0):+.3e} (log dens), d corr/d beta {dcorr:+.3e}, post sd {curv**-0.5:.2e}, bias/sd {dcorr/np.sqrt(curv):.2e}")
