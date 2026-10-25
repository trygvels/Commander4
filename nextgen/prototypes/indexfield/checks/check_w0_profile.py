"""Profile scan of the collapsed log-likelihoods (w0, w1, gauss) in one region's beta, others at
truth, for the coverage-probe replicate (seed 20261025, rep 0, noise 1). Shows whether the ridge
posterior's inflated field variance comes from the target rather than the sampler."""
import sys, numpy as np, healpy as hp
sys.path.insert(0, "/Users/trygvels/code/Commander4/nextgen/prototypes")
from indexfield import field, posterior, coverage, model
nside_d, nside_beta, seed, rep, noise = 16, 2, 20261025, 0, 1.0
basis = field.GraphLaplacianBasis(nside_beta); fld = field.MaternField(basis)
hyper = posterior.HyperPrior(rho_prior="loguniform").resolved(basis.h)
rng = np.random.default_rng([seed, rep, int(round(1000 * noise))])
data, truth = coverage.simulate_template(nside_d, nside_beta, fld, hyper, rng, noise, None, "default")
snr_index, snr_pix = coverage.snr_proxy(data, truth["beta"])
beta_t = truth["beta"].copy()
grid = np.linspace(-7.0, 1.0, 801)
print(f"truth mu {truth['theta']['mu']:.3f} sigma {truth['theta']['sigma']:.3f}; beta range "
      f"{beta_t.min():.2f} {beta_t.max():.2f}")
ml = {m: np.zeros(len(beta_t)) for m in ("w0", "w1", "gauss")}
width = {m: np.zeros(len(beta_t)) for m in ("w0", "w1", "gauss")}
for r in range(len(beta_t)):
    for m in ml:
        ll = np.empty(grid.size)
        for i, b in enumerate(grid):
            bb = beta_t.copy(); bb[r] = b
            ll[i] = model.collapsed_loglike(bb, data, m)[0]
        ll -= ll.max(); p = np.exp(ll); p /= p.sum()
        mean = (grid * p).sum(); sd = np.sqrt(((grid - mean) ** 2 * p).sum())
        ml[m][r] = grid[np.argmax(ll)]; width[m][r] = sd
for r in np.argsort(snr_index)[[0, 1, 12, 20, 24, 47]]:
    print(f"region {r:2d} S/N_index {snr_index[r]:5.2f} truth {beta_t[r]:6.2f} | "
          + " | ".join(f"{m}: ML {ml[m][r]:6.2f} sd {width[m][r]:5.2f}" for m in ml))
for m in ml:
    print(f"{m:5s}: ML scatter about truth, all regions: rms {np.sqrt(np.mean((ml[m]-beta_t)**2)):.3f}, "
          f"low S/N (<1) rms {np.sqrt(np.mean((ml[m]-beta_t)[snr_index<1]**2)):.3f}, "
          f"median conditional sd {np.median(width[m]):.3f}, in low S/N {np.median(width[m][snr_index<1]):.3f}")
