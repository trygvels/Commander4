"""Independent check of the Commander3 comm_mh_specind_mod move.

Model: d_{p,nu} = f(nu;beta) a_p + n, n ~ N(0, sig^2). k pixels, nb bands, one beta.
Flat priors on a and beta. Compare:
  exact_marg : grid of p(beta|d) = int p(d|beta,a) da  ~ exp(-chi2min/2) |M|^-1/2
  exact_prof : grid of exp(-chi2min/2)
  gibbs_none : UPDATE_CG=none inside Gibbs: a~p(a|beta,d); MH on beta given a   (valid Gibbs)
  cg_numstep1: UPDATE_CG=group, NUMSTEP=1 inside Gibbs (a redrawn at entry each iteration)
  cg_joint   : UPDATE_CG=group, many steps without external a redraw (joint kernel)
  optimize   : OPERATION=optimize (a = ahat(beta')), deterministic -> profile
"""
import numpy as np, sys
rng = np.random.default_rng(1)
nu = np.array([30., 44., 70.]); nu0 = 30.
beta_true = -3.0

def run(k, nstep=200000, sig=5.0, step=None):
    a_true = 100 + 20*rng.standard_normal(k)
    f = lambda b: (nu/nu0)**b
    d = np.outer(a_true, f(beta_true)) + sig*rng.standard_normal((k, len(nu)))
    def Mscalar(b): return np.sum(f(b)**2)/sig**2
    def ahat(b): return d @ f(b) / np.sum(f(b)**2)
    def chi2(b, a): return np.sum((d - np.outer(a, f(b)))**2)/sig**2
    def chi2min(b): return chi2(b, ahat(b))
    # exact grid
    grid = np.linspace(beta_true-0.5, beta_true+0.5, 4001)
    lp_prof = np.array([-0.5*chi2min(b) for b in grid])
    lp_marg = lp_prof - 0.5*k*np.log(np.array([Mscalar(b) for b in grid]))
    def stats(lp):
        w = np.exp(lp-lp.max()); w/=w.sum()
        m = np.sum(w*grid); return m, np.sqrt(np.sum(w*(grid-m)**2))
    m_marg, s_marg = stats(lp_marg); m_prof, s_prof = stats(lp_prof)
    if step is None: step = 1.5*s_marg
    def draw_a(b): return ahat(b) + rng.standard_normal(k)/np.sqrt(Mscalar(b))

    out = {}
    # gibbs_none
    b = beta_true; acc=0; ch=[]
    for i in range(nstep):
        a = draw_a(b); c_old = chi2(b, a)
        bp = b + step*rng.standard_normal(); c_new = chi2(bp, a)
        if np.log(rng.uniform()) < -(c_new-c_old)/2: b = bp; acc+=1
        ch.append(b)
    out['gibbs_none'] = (np.mean(ch), np.std(ch), acc/nstep)
    # cg_numstep1 (composed kernel, as in mbbtab or specind NUMSTEP=1 in Gibbs loop)
    b = beta_true; acc=0; ch=[]
    for i in range(nstep):
        a = draw_a(b); c_old = chi2(b, a)
        bp = b + step*rng.standard_normal(); ap = draw_a(bp); c_new = chi2(bp, ap)
        if np.log(rng.uniform()) < -(c_new-c_old)/2: b = bp; acc+=1
        ch.append(b)
    out['cg_numstep1'] = (np.mean(ch), np.std(ch), acc/nstep)
    # cg_joint (joint (beta,a) kernel, no external redraw)
    b = beta_true; a = draw_a(b); acc=0; ch=[]
    for i in range(nstep):
        c_old = chi2(b, a)
        bp = b + step*rng.standard_normal(); ap = draw_a(bp); c_new = chi2(bp, ap)
        if np.log(rng.uniform()) < -(c_new-c_old)/2: b, a = bp, ap; acc+=1
        ch.append(b)
    out['cg_joint'] = (np.mean(ch), np.std(ch), acc/nstep)
    # optimize
    b = beta_true; acc=0; ch=[]
    for i in range(nstep):
        c_old = chi2min(b)
        bp = b + step*rng.standard_normal(); c_new = chi2min(bp)
        if np.log(rng.uniform()) < -(c_new-c_old)/2: b = bp; acc+=1
        ch.append(b)
    out['optimize'] = (np.mean(ch), np.std(ch), acc/nstep)
    print(f"\nk={k}  (MC se on mean ~ std/sqrt(ESS); nstep={nstep}, step={step:.4f})")
    print(f"  exact marginal : mean={m_marg:+.5f} std={s_marg:.5f}")
    print(f"  exact profile  : mean={m_prof:+.5f} std={s_prof:.5f}")
    for kk,v in out.items():
        print(f"  {kk:15s}: mean={v[0]:+.5f} std={v[1]:.5f} acc={v[2]:.3f}  std/std_marg={v[1]/s_marg:.2f}")

for k in [1, 4, 16, 64]:
    run(k, nstep=int(sys.argv[1]) if len(sys.argv)>1 else 100000)
