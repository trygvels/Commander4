import numpy as np
rng = np.random.default_rng(3)
# sample_template_mh proposal x' = x(1+N(0,s)), accept min(1,L'/L); L = N(x; x0, sx)
for x0, sx, s, nstep in [(10,2,0.05,5000),(10,2,0.05,50000),(10,0.5,0.05,50000),(10,0.5,0.02,50000),(10,2,0.2,5000)]:
    means=[]; collapsed=0
    for rep in range(40):
        x=x0; ch=np.empty(nstep)
        for i in range(nstep):
            xp = x*(1+s*rng.standard_normal())
            if np.log(rng.uniform()) < -((xp-x0)**2-(x-x0)**2)/(2*sx**2): x=xp
            ch[i]=x
        means.append(ch.mean()); collapsed += (ch[-1] < x0-5*sx)
    means=np.array(means)
    print(f"x0={x0} sx={sx} s={s} nstep={nstep}: chain mean over 40 reps = {means.mean():.4f} +- {means.std()/np.sqrt(40):.4f} "
          f"(exact {x0}); reps ending collapsed (x<x0-5sx): {collapsed}/40")
