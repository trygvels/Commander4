# Independent stats review of the E026 sync-level claims (stats-reviewer, 2026-11-02, 17:50)

Scripts in `review_scripts/` (run from that folder; `decomp.py` rebuilds the target in numpy and
matches Commander4 to 1e-4 nats at draw iter0501, 1280908.95 both; the full target reproduces the
conditional-check log value 1275513.32).

1. Conditional-CDF (Gibbs-consistency) check: valid as a necessary condition (the probability
   integral transform of p(mu_c | rest) at the draw is U(0,1) under the target); implementation
   matches; cannot detect an error in the target itself, errors in unscanned coordinates, or
   D < 0.3 with K = 20. Sync: PIT mean 0.428 +- 0.065, KS p 0.15, optimum minus draw +0.0045 +-
   0.0032 (se): consistent. Dust: PIT mean 0.401, KS p 0.055; optimum minus draw +0.0007 with se
   0.0003 (2.3 sigma, irrelevant in absolute terms): borderline, step 2.1 sd too coarse. Extra
   PITs along log sigma_s (p 0.07), z_0 (0.91), z_50 (0.88), z_1000 (0.25), u (0.47): consistent.
   Recommendation: K = 100, dust step at most 1 sd, a negative control (corrupted target).
2. Roughness coupling: NOT supported. The block-shift optimum is linear in the z scale (draw 0:
   +0.076, +0.038, +0.019 at scales 1, 0.5, 0.25), the first-order response of the data-pinned
   bright regions moving off their optima when z is shrunk (z x 0 is the scalar model, optimum
   -2.93 to -2.99 = the scalar run). With the sync multipliers integrated exactly (Gaussian in
   mu_r, checked against brute force to 1e-8) the faint-region data level at the chain state is
   -2.79 +- 0.065 for a constant, the reference, the chain's, a shuffled and an iid N(0, 0.15)
   pattern alike (-2.784 to -2.796): roughness changes nothing at the 0.02 level. Analytic
   smoothing shift: +0.002 / -0.005. The cause of +0.13: the marginalization over the 3072
   per-region multipliers with hyperprior N(1, 0.5^2) when the truth is 1: the volume term
   -1/2 sum_r log(-Q_r) has slope +130 nats per unit beta over the faint regions (+251 whole sky;
   66 percent of the sync amplitude information is at 23 GHz, a steeper index gives each mu_r
   more leverage). Net of the multipliers' noise fitting, the exact marginal puts the faint level
   +0.19 above the reference in this realization and +0.15 +- 0.05 over 6 simulated noise
   realizations (all positive): a model bias. The chain's faint regions at +0.13 are held back
   by the field prior and the mu hyperprior; 0.25 x (-2.90) + 0.75 x (-2.857) = -2.868 = mu.
   Hyperprior sd 0.2 or 0.1: Occam slope +30 / +8, faint bias +0.06 / zero within noise.
   With the multipliers integrated the whole-sky optimum at the posterior-mean field is -0.003
   +- 0.025: the chain state is the target optimum. The record's scans profiled 3072 multipliers
   (Neyman-Scott bias toward steep, -0.2 to -0.4 here) or used one global multiplier (a
   different model, -0.109 / 7.3 nats, faint level -2.997): neither is a stand-in for the
   chain's marginal; conclusions from those cells need rechecking.
3. Scalar run: the pseudo-true value under the hierarchical model (exact 1-d marginal with the
   multipliers integrated, 8 noise realizations, dust at the scalar run's state) is -2.956 +-
   0.003, the exact marginal for this realization -2.951 +- 0.025 (dust draws -2.927 to -2.962),
   the chain -2.944 +- 0.031 (MC se 0.003): the chain agrees with an exact marginal. The right
   yardstick: realization scatter 0.009 plus dust-state dependence 0.015 plus MC 0.003; the
   record's "+0.039 (1.3 sd)" used the wrong reference (-2.983, amplitude-known Fisher weights)
   and the wrong yardstick. Caveat: with the dust at the truth (region-mean beta_d, mu_d 1) the
   scalar pseudo-true is -2.83: the sync index is sensitive at the 0.1 level to the fixed-T dust
   state.
