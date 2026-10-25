"""Independent check of the BP10 'marginal' lnL (comp_lnL_max_chisq_diagonal with use_det).

Toy: one pixel, nb bands, single component with amplitude a and spectral parameter theta.
Model d_k = M_k(theta) a + n_k, n_k ~ N(0, sigma_k^2). Flat prior on a.
True marginal: ln int exp(-chi2/2) da = -chi2_min/2 - 0.5 ln(M^T N^-1 M) + 0.5 ln(2 pi).
Code: -chi2_min/2 - 0.5 ln(invMNM) = -chi2_min/2 + 0.5 ln(M^T N^-1 M).
We compare theta-dependence (differences across theta), since constants drop out of MH.
"""
import numpy as np
rng = np.random.default_rng(1)
nu = np.array([30., 44., 70., 100.])
sig = np.array([3., 2., 1.5, 1.0])
a_true, th_true = 50.0, -3.0
def M(th): return (nu/30.)**th
d = M(th_true)*a_true + rng.normal(0, sig)
invN = 1/sig**2

def code_max_chisq(th, use_det):
    m = M(th); MN = m*invN; MNd = (MN*d).sum(); MNM = (MN*m).sum()
    invMNM = 1/MNM; amp = MNd*invMNM
    chisq = (((d - m*amp)**2)*invN).sum()
    out = -0.5*chisq
    if use_det: out -= 0.5*np.log(invMNM)   # exact transcription of f90 line 3513
    return out

def code_marginal_diag(th, use_det):
    m = M(th); MN = m*invN; MNd = (MN*d).sum(); MNM = (MN*m).sum()
    invMNM = 1/MNM
    out = 0.5*MNd*invMNM*MNd
    if use_det: out -= 0.5*np.log(invMNM)   # f90 line 3426
    return out

def quad_marginal(th):
    m = M(th)
    MNM = (m*m*invN).sum(); MNd = (m*d*invN).sum(); ahat = MNd/MNM
    w = 20/np.sqrt(MNM)
    a = np.linspace(ahat-w, ahat+w, 200001)
    chisq = (((d[:,None] - m[:,None]*a[None,:])**2)*invN[:,None]).sum(0)
    return np.log(np.trapezoid(np.exp(-0.5*(chisq-chisq.min())), a)) - 0.5*chisq.min()

ths = np.linspace(-3.6, -2.4, 7)
q = np.array([quad_marginal(t) for t in ths])
cm = np.array([code_max_chisq(t, True) for t in ths])
cr = np.array([code_max_chisq(t, False) for t in ths])
cd = np.array([code_marginal_diag(t, True) for t in ths])
flip = np.array([code_max_chisq(t, False) - 0.5*np.log((M(t)**2*invN).sum()) for t in ths])
print("theta       quad_marg   code_marg   ridge   code_marg_diag  ridge-0.5lnMNM  (all minus value at theta[0])")
for i,t in enumerate(ths):
    print(f"{t:7.3f} {q[i]-q[0]:12.5f} {cm[i]-cm[0]:11.5f} {cr[i]-cr[0]:9.5f} {cd[i]-cd[0]:12.5f} {flip[i]-flip[0]:12.5f}")
print("max |quad - code_marginal| (theta-dependent part):", np.abs((q-q[0])-(cm-cm[0])).max())
print("max |quad - (ridge - 0.5 ln MNM)| :", np.abs((q-q[0])-(flip-flip[0])).max())
print("max |code_marg_diag - code_max_chisq| (should be const):", np.abs((cd-cd[0])-(cm-cm[0])).max())
# posterior mean of theta under each
pr = np.linspace(-4.5,-1.5,3001)
def pm(f):
    l = np.array([f(t) for t in pr]); w = np.exp(l-l.max()); return (pr*w).sum()/w.sum(), np.sqrt(((pr**2)*w).sum()/w.sum()-((pr*w).sum()/w.sum())**2)
print("post mean/std theta: quad   ", pm(quad_marginal))
print("post mean/std theta: code   ", pm(lambda t: code_max_chisq(t, True)))
print("post mean/std theta: ridge  ", pm(lambda t: code_max_chisq(t, False)))
