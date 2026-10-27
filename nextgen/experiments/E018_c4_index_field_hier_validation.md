# E018: Commander4 index field with a hierarchical amplitude prior (likelihood: hier) against the prototype and independent references (M2.4 stage 1, L1 and L3)

- **Date started / finished:** 2026-10-27 /
- **Branch and commit:** nextgen, 2e3fa1b (commits 1 to 8 of the M2.4 stage-1 plan)
- **Open question addressed:** Q2 (prior calibration), Q3 (spatially varying indices in Commander4); ROADMAP M2.4
- **Validation level:** L1 (exact-target tests in process) and L3 (end-to-end 12-rank Commander4 run on a model-consistent simulation)

## Question
Does the Commander4 implementation of the hierarchical index-field target (`RegionIndexLikelihood`
plus `IndexFieldSamplingGroup`, `likelihood: hier`) sample the distribution stated in
`nextgen/plans/M2.4_commander4_index_field.md` section 1, as checked (L1) against references
that share no code with it and (L3) against the E016 prototype chain on the same data through the
full Gibbs loop (NUTS over MPI, amplitude redraw, chain file)? If yes, M2.4a is done up to the L2
gate (prototype hier SBC) and the method can be run on E014-style real inputs; if the end-to-end
chain differs from the in-process target, the MPI driver, the amplitude redraw or the chain IO is
wrong, not the target.

## Target and priors
Exactly section 1 of the plan (stats-reviewed twice, 2026-10-26 and 2026-10-27). Per pixel p at
nside_d = 16 (3072 pixels), 12 bands (prototype `simulate.BANDS_GHZ`: 23 to 353 GHz), components
cmb, sync, dust in uK_RJ; d_p = M_p a_p + n_p, n_p ~ N(0, diag(rms_b^2)) white (no smoothing,
all bands declared at one common beam so no post-smoothing takes place, F7/F8 do not enter).
- sync index field at nside_beta = 2 (Nb = 48, K = 47): beta = mu 1 + sigma S(rho) z, graph
  Laplacian basis, Matern alpha 2, mu ~ N(-3.0, 0.3^2), sigma ~ half-N(0, 0.3), rho log-uniform on
  [h, pi], h = sqrt(4 pi / 48). Non-centred parameterization.
- dust index fixed at 1.6 (T 20 K, pivot 545 GHz), sync pivot 30 GHz, cmb column cmb_to_rj(nu)
  (nu_ref 1e-3 GHz so that the column equals the prototype's).
- amplitude priors: cmb a ~ N(0, 70^2) fixed; sync and dust hier, a_p,c ~ N(mu_c g_p,c,
  (s_c g_p,c)^2) with g the smooth plane-brightening templates of `simulate.templates` (sync at
  30 GHz, dust at 545 GHz), mu_c ~ N(1, 0.5^2), log s_c ~ N(log 0.3, 1^2), guard log s >= -8.
- collapsed likelihood with log det kept; joint target over (z, mu, log sigma, u, mu_s, log s_s,
  mu_d, log s_d), dim 54.
L1 targets (tests/test_index_field_l1.py, nside_d 4, 7 bands): (a) one scalar sync index with a
fixed Gaussian amplitude prior a ~ N(m, T) (m = 1.03 g, T = (0.23 g)^2, cmb N(0, 70^2)) and
index prior N(-3, 0.2^2): a one-dimensional posterior; (b) the sync field at nside_beta 1
(K = 11) with the same fixed amplitude prior and the field hyperpriors above: the E012 v2 field
target in centred coordinates.

## Setup
Simulation (`aux/validate_mcmc_likelihood.py make --mode hier_field --dir <dir>`, seed 2026):
truth drawn from the stated model with the field hyperparameters fixed at mu = -3.0, sigma =
0.15, rho = sqrt(h pi) (the middle of the log range); amplitudes a_p,c = g_p,c (1.0269 + 0.2336
N(0,1)) for sync and dust (the Gaussian twin of the template scatter, E016), cmb N(0, 70^2);
white noise `simulate.RMS_UK_RJ`. Written as FITS maps per band (map and rms), the two
templates, and `truth.npz`. Commander4: 12 CompSep ranks (one per band), per-pixel amplitude
group over all three components, MCMC group `sample_index_fields` with `likelihood: hier`,
`numstep: 1`, `nuts: {n_warmup: 500, target_accept: 0.8, max_depth: 10, parameterization:
noncentred}`, 3000 Gibbs iterations, chain written every iteration (`param_hier_field.yml`). A
second run `param_hier_scalar.yml` with `spectral_index_nside: 0` (one sync index plus the hier
hyperparameters, the E017 scalar target) on the same data. Machine: Mac mini (M-series), 12
ranks, 1 thread each.

## Reference
L1: (a) `reference.grid_1d(data, "gauss", grid, theta_fixed)` of the prototype (plain NumPy
loops, no shared code with the kernel): posterior mean and sd on a fine grid; the Commander4
chain (RegionIndexLikelihood through NUTSSampler, 500 warm-up + 2000) must have mean within 4
MC SE (sd / sqrt(ESS)) and sd within 10 percent. (b) `reference.centred_gibbs_reference(data,
fspec, hyper, "gauss", ...)`: MH-within-Gibbs in centred coordinates with exact inverse-CDF
draws of (log sigma, u) on a grid, the E012 reference 1c, which shares neither the basis
construction (`reference.field_spec`) nor the likelihood code; per index pixel |mean difference|
within 4 pooled MC SE and sd ratio within 0.8 to 1.25, and the same for mu, log sigma, log rho.
L3: the prototype `HierNonCentredPosterior` (E016, SBC-validated at K 192 in E012/E016) sampled
with `sampler.run_nuts` (500 + 3000) on the PixelData rebuilt from the FITS maps: per index pixel
and per hyperparameter, |mean difference| within 4 pooled MC SE (Geyer ESS on both chains), sd
ratio within 0.85 to 1.15, Commander4 divergences below 1 percent after warm-up and min ESS
above 100 per 1000 iterations; plus z-scores against the truth as a sanity check (the field is a
single realization, so only gross miscalibration would show). The scalar run is compared with
the prototype `ScalarHierPosterior` chain the same way.

## Command
```bash
source .venv/bin/activate
python aux/validate_mcmc_likelihood.py make --mode hier_field --dir /Users/trygvels/data/c4val_hier
mpirun -n 12 commander4 -p /Users/trygvels/data/c4val_hier/param_hier_field.yml
mpirun -n 12 commander4 -p /Users/trygvels/data/c4val_hier/param_hier_scalar.yml
python aux/validate_mcmc_likelihood.py analyze --mode hier_field --dir /Users/trygvels/data/c4val_hier
python -m pytest -q -p no:cacheprovider tests/test_index_field_l1.py
```

## Result
**L1 (tests/test_index_field_l1.py, 2026-10-27 02:05, C4_SLOW_TESTS=1, 141 s): PASS.** (a) The
scalar chain (300 + 2000) against the reference grid: mean within 4 MC SE and width within 10
percent (ESS > 200, no divergences). (b) The nside-1 field (K 11) with the synchrotron 30 times
brighter than the template (per-pixel index error 0.2 against a field scatter of 0.15): all 12
beta means, mu, log sigma and log rho within 4 pooled SE of the centred Gibbs reference (20000
iterations, 5000 burn), sd ratios within 0.8 to 1.25. Diagnostic chains behind the test (same
data, 500 + 6000 NUTS, Gibbs 20000): log sigma -1.675 +- 0.014 (non-centred NUTS, ESS 1485),
-1.653 +- 0.029 (centred NUTS, ESS 263), -1.660 +- 0.035 (reference Gibbs, ESS 198); mu -3.0736
+- 0.0010 / -3.0730 / -3.0695 +- 0.0025; log rho +0.581 / +0.571 / +0.582 (+- 0.004); beta_0
-3.0758 / -3.0702 / -3.0706 (+- 0.004). At the template brightness (sync 30 times fainter) the
same comparison fails on log sigma only: the posterior of log sigma has a long tail towards small
sigma (sd 1.1 in the log), the non-centred NUTS chain gives -1.83 +- 0.02 while the centred
reference gives -1.31 (6000 iterations, ESS 84) drifting to -1.52 (20000, ESS 110), i.e. the
centred MH-within-Gibbs reference does not reach the funnel; beta, mu and log rho agree there
too. The centred target itself was checked against the reference's density pieces
(`loglike_pixel_loops` + `hyper_conditional_logdensity` + `_lp_mu`) at five states: constant
offset to 6e-12. So (b) is a sampler-mixing limitation of the reference at low S/N, not a
discrepancy of the target; the test is run where the reference mixes.

**L3:** chains started 02:14 (hier_field, 12 ranks, about 5 s per Gibbs iteration under load;
hier_scalar follows). (pending)

## Negative control
The same hier_field analysis applied to a chain run with `likelihood: ridge` on the scalar
configuration (profile likelihood, no amplitude prior) must disagree with the hier prototype
chain by many sd (different target; E017 shows the profile is 0.1 to 1 sd off the hier mean at
these S/N and has no hyperparameters), and the field run started from a constant field must leave
the constant start (sigma moves away from its 0.15 initial value only if the data demand it: the
truth has sigma 0.15, so the posterior of log sigma must cover log 0.15).

## Conclusion
(pending)

## Follow-ups
- L2 gate: the prototype hier SBC (E016, 300 replicates) must pass before M2.4a is declared done.
- Real-data configuration of E014 with a sync field at nside 2 to 4 (needs Q15 resolved first).
