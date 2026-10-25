import numpy as np, sys
rng = np.random.default_rng(2)
nu = np.array([30., 44., 70.]); nu0 = 30.; beta_true=-3.0
def setup(k, sig=5.0):
    a_true = 100 + 20*rng.standard_normal(k)
    f = lambda b: (nu/nu0)**b
    d = np.outer(a_true, f(beta_true)) + sig*rng.standard_normal((k, len(nu)))
    M = lambda b: np.sum(f(b)**2)/sig**2
    ahat = lambda b: d @ f(b) / np.sum(f(b)**2)
    chi2 = lambda b, a: np.sum((d - np.outer(a, f(b)))**2)/sig**2
    chi2min = lambda b: chi2(b, ahat(b))
    grid = np.linspace(beta_true-3, beta_true+3, 30001)
    lp_prof = np.array([-0.5*chi2min(b) for b in grid])
    lp_marg = lp_prof - 0.5*k*np.log(np.array([M(b) for b in grid]))
    def stats(lp):
        w = np.exp(lp-lp.max()); w/=w.sum(); m=np.sum(w*grid); return m, np.sqrt(np.sum(w*(grid-m)**2))
    draw_a = lambda b: ahat(b) + rng.standard_normal(k)/np.sqrt(M(b))
    return chi2, chi2min, draw_a, stats(lp_marg), stats(lp_prof)

mode = sys.argv[1]
if mode == 'composed':
    # UPDATE_CG=group with a redrawn from p(a|beta,d) before each MH step (mbbtab always; specind NUMSTEP=1)
    for k in [4, 64, 256, 1024]:
        chi2, chi2min, draw_a, (mm, sm), (mp, sp) = setup(k)
        nstep = 100000; step = 1.5*sm; b=beta_true; ch=[]; acc=0
        for i in range(nstep):
            a = draw_a(b); c_old = chi2(b,a)
            bp = b + step*rng.standard_normal(); ap = draw_a(bp); c_new = chi2(bp,ap)
            if np.log(rng.uniform()) < -(c_new-c_old)/2: b=bp; acc+=1
            ch.append(b)
        ch=np.array(ch)
        print(f"k={k:5d} exact marg mean={mm:+.4f} std={sm:.4f} | prof mean={mp:+.4f} std={sp:.4f} | "
              f"composed mean={ch.mean():+.4f} std={ch.std():.4f} ratio={ch.std()/sm:.2f} acc={acc/nstep:.3f}  "
              f"sqrt(2k)={np.sqrt(2*k):.1f}")
elif mode == 'joint':
    # UPDATE_CG=group, many NUMSTEP without external redraw; prediction: theta-marginal = profile
    k = 8
    chi2, chi2min, draw_a, (mm, sm), (mp, sp) = setup(k)
    nstep = 1500000; step = 1.5*sm; b=beta_true; a=draw_a(b); ch=np.empty(nstep); acc=0
    c_old = chi2(b,a)
    for i in range(nstep):
        bp = b + step*rng.standard_normal(); ap = draw_a(bp); c_new = chi2(bp,ap)
        if np.log(rng.uniform()) < -(c_new-c_old)/2: b,a,c_old = bp,ap,c_new; acc+=1
        ch[i]=b
    # batch-means standard error
    nb=50; bm = ch.reshape(nb,-1).mean(1); se = bm.std(ddof=1)/np.sqrt(nb)
    print(f"k={k} exact marg mean={mm:+.5f} std={sm:.5f} | prof mean={mp:+.5f} std={sp:.5f}")
    print(f"joint kernel: mean={ch.mean():+.5f} (se {se:.5f}) std={ch.std():.5f} acc={acc/nstep:.3f}")
    print(f"  (mean - prof)/se = {(ch.mean()-mp)/se:+.2f}   (mean - marg)/se = {(ch.mean()-mm)/se:+.2f}")
elif mode == 'scale':
    # sample_template_mh: x' = (1 + N(0,s)) x, accept min(1, L'/L), no Hastings term.
    x0, sx = 10.0, 2.0   # likelihood N(x; x0, sx)
    for s in [0.05, 0.2, 0.5]:
        nstep=1000000; x=x0; ch=np.empty(nstep); acc=0
        for i in range(nstep):
            xp = x*(1+s*rng.standard_normal())
            if np.log(rng.uniform()) < -((xp-x0)**2-(x-x0)**2)/(2*sx**2): x=xp; acc+=1
            ch[i]=x
        nb=50; bm=ch.reshape(nb,-1).mean(1); se=bm.std(ddof=1)/np.sqrt(nb)
        # leading-order prediction: stationary ~ |x| N(x;x0,sx) -> mean = x0 + sx^2/x0 (for x0>>sx)
        print(f"s={s:.2f}: mean={ch.mean():.4f} (se {se:.4f}) std={ch.std():.4f} acc={acc/nstep:.3f}; "
              f"exact x0={x0}, predicted |x|L mean={x0+sx**2/x0:.4f}")
