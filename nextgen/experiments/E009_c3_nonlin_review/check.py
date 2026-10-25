import numpy as np
from math import erf, sqrt, exp, log, pi
rng = np.random.default_rng(1)

# (a) Redraw-until-inside proposal with symmetric-ratio acceptance.
# Derived claim: stationary density is pi(x)*Z(x), Z(x)=P(x+eps in B).
def Phi(z): return 0.5*(1+erf(z/sqrt(2)))
lo = -1.0; sp = 1.0            # bound and proposal std (exaggerated to make effect visible)
def logpi(x): return -0.5*x*x # N(0,1) truncated at x>lo
N = 2_000_000
x = 0.0; lp = logpi(x); xs = np.empty(N)
eps = rng.normal(0, sp, size=4*N); u = rng.random(N); k = 0
for i in range(N):
    while True:
        y = x + eps[k]; k += 1
        if y > lo: break
    lpy = logpi(y)
    if lpy >= lp or u[i] < exp(lpy-lp): x, lp = y, lpy
    xs[i] = x
xs = xs[N//10:]
grid = np.linspace(lo, 3.5, 46); c = 0.5*(grid[1:]+grid[:-1]); dx = grid[1]-grid[0]
h, _ = np.histogram(xs, grid, density=True)
pi_t = np.exp(-0.5*c*c); pi_t /= pi_t.sum()*dx
Z = np.array([Phi((v-lo)/sp) for v in c]); piZ = pi_t*Z; piZ /= piZ.sum()*dx
print("boundary redraw: max|h-pi|/max(pi) = %.3f   max|h-piZ|/max(pi) = %.3f" %
      (np.max(np.abs(h-pi_t))/pi_t.max(), np.max(np.abs(h-piZ))/pi_t.max()))
print("  mean(x): chain %.4f  pi %.4f  piZ %.4f" % (xs.mean(), (c*pi_t).sum()*dx, (c*piZ).sum()*dx))

# BP10 synch numbers: bounds [-4.5,-1.5], posterior near -3.1, step ~ steplen*L ~ 0.01..0.1
for s in (0.01, 0.05, 0.1, 0.3):
    Zmin = Phi((-1.5-(-3.1+0.3))/s) - Phi((-4.5-(-3.1-0.3))/s)
    print("  BP10 synch: step %.2f  min Z over beta in [-3.4,-2.8] = %.3e  (log-bias %.1e)" % (s, Zmin, -log(Zmin)))

# (b) low-res chisq weighting: code uses Nsub*mean(1/sig^2); exact for udgraded white residual is Nsub/mean(sig^2)
sig = rng.lognormal(0, 0.5, size=(10000, 4096))   # 4096 sub-pixels per nside-16 pixel from nside 1024
code  = 4096*np.mean(1/sig**2, axis=1)
exact = 4096/np.mean(sig**2, axis=1)
print("lowres invN ratio code/exact for lognormal(0.5) rms scatter: mean %.3f" % np.mean(code/exact))
sig = rng.lognormal(0, 0.1, size=(10000, 4096))
print("lowres invN ratio code/exact for lognormal(0.1) rms scatter: mean %.3f" % np.mean(4096*np.mean(1/sig**2,1)/(4096/np.mean(sig**2,1))))

# (c) band limit: 600 arcmin Gaussian b_l at l=100 (lmax_ind) and at l=48 (3*nside16)
sg = 600/60*pi/180/sqrt(8*log(2))
for l in (48, 100): print("b_l(600') at l=%d: %.2e" % (l, exp(-0.5*l*(l+1)*sg**2)))
