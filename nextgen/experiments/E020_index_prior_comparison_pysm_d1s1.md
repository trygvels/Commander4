# E020: Index-field prior comparison on a PySM d1s1 sky: hierarchical GMRF with inferred smoothness against independent regions at a fixed nside and a single index (M2.5 first cut, Q3)

- **Date started / finished:** 2026-10-27 /
- **Branch and commit:** nextgen, a07f538
- **Open question addressed:** Q3 (optimal spatially varying index model); ROADMAP M2.5
- **Validation level:** L4 (misspecified realistic sky with known truth; the samplers themselves
  are validated at L1 to L3 in E018)

## Question
On a realistic, model-misspecified sky (PySM d1 dust with spatially varying beta_d and T_d, s1
synchrotron with spatially varying beta_s, c1 CMB) at nside 64 in intensity with Planck-like
white noise, does the hierarchical GMRF index field with inferred mean, scale and correlation
length recover the index maps and the CMB better than (a) independent indices per region at a
fixed nside (the fgbuster-style ladder, here with the same collapsed likelihood and amplitude
prior) and (b) a single full-sky index, and is the GMRF answer stable when nside_beta is refined
from 8 to 16 (that is, does the inferred smoothness make nside_beta a nuisance choice)? If yes,
the GMRF field is the Phase 2 index model and the remaining Q3 variants (alm GRF, clustering,
heavy-tailed increments) are compared against it later; if the fixed ladder is as good, the
field prior buys nothing at this S/N and the methods paper says so.

## Target and priors
Revised 2026-10-27 after the stats review of the first design (its seven "must fix" findings are
listed under "Design review" below; the first design is kept in git history, commit 21634bd).

Data: 12 bands at nu_b = 23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353 GHz (delta bandpasses),
intensity, nside 64 (N = 49152 pixels), d_p = sum_c M_pc(beta) a_pc + n_p, n_p ~ N(0, N_p)
diagonal (white, known rms per band: 8.4, 7.2, 7.2, 6.6, 7.2, 7.2, 8.4, 4.8, 2.9, 1.9, 2.9, 9.6
uK_RJ per nside-64 pixel, within a factor 2 of WMAP/Planck at this resolution). Model components
(the data are not drawn from this model): cmb (fixed SED, the thermodynamic-to-RJ conversion),
synchrotron power law a_s (nu/30)^beta_s(p), thermal dust modified blackbody with pivot 545 GHz,
beta_d(p), fixed T_d = 20.93 K (the mean of the d1 temperature map at nside 64; the true map runs
from 14 to 31 K, so the dust model is misspecified by construction: the noise-free fixed-T best-fit
index beta_eff(p) differs from the PySM beta_d by rms 0.028, range -0.13 to +0.08, and is the dust
truth a one-MBB fixed-T model can at best recover; both are stored and both are reported).

Collapsed likelihood: amplitudes integrated out analytically per pixel under the amplitude prior
of E018 (`likelihood: hier`): a_p,cmb ~ N(0, 100^2) (the c1 map rms is 84 uK); a_p,c ~ N(mu_c g_p,c,
(s_c g_p,c)^2) for c in {sync, dust} with the pseudo-templates g defined below, hyperpriors mu_c ~
N(1, 0.5^2), log s_c ~ N(log 0.3, 1^2), guard log s_c >= -8. For a state (beta, mu, s) the density
is the E018 target (M2.4 plan section 1): prod_p N(d_p; M m, N_p + M T M^T) evaluated as 1/2 b'^T
A'^-1 b' - 1/2 log det A' - 1/2 m^T T^-1 m - 1/2 log det T summed over pixels, with A' = M^T N^-1 M
+ T^-1, b' = M^T N^-1 d + T^-1 m. Reviewed and validated in E018 (F19). The amplitude hyperparameters
(mu_c, s_c) are sampled in every option below.

Index priors compared, the same one for both components c in {sync, dust}, sampled jointly with
the amplitude hyperparameters by NUTS (E018 sampler, non-centred field coordinates, initial
diagonal mass from a finite-difference Hessian at the start state, `nuts.init_inv_mass:
fd_hessian`, which changes the adaptation start and not the target). Only the spatial prior
differs between the options; likelihood, amplitude prior and hyperpriors are shared:
1. **GMRF field at nside_beta** in {4, 8, 16}: beta_c = mu_c 1 + sigma_c S(rho_c) z_c, z ~ N(0, I_K),
   K = 12 nside_beta^2 - 1 (graph-Laplacian basis without the constant mode), Matern alpha 2 on the
   HEALPix neighbour graph; mu_c ~ N(mu0_c, 0.3^2), sigma_c ~ half-N(0, 0.3^2), rho_c log-uniform on
   [h, pi] through the logit-transformed u with its Jacobian. Exactly the E016/E018 prior. Its
   marginal per-region prior is N(mu_c, 0.3^2 + sigma_c^2 var_r(rho)) with sigma and rho learned
   from the whole sky (mostly the plane) and applied everywhere (stationarity); that, and the
   coupling of neighbours, is what is being tested.
2. **Exchangeable (iid) hierarchical field at nside_beta** in {4, 8, 16} (`spectral_index_field:
   {type: iid}`): beta_c,r = mu_c + sigma_c z_r, z ~ N(0, I), with the same hyperpriors on mu_c and
   sigma_c as option 1 and no spatial coupling (the rho -> 0 limit of the GMRF; the GMRF's u is
   kept as a decoupled nuisance). The nested analogue: option 1 minus spatial correlation. The
   scalar option 3 is its sigma -> 0 limit.
3. **Independent regions at nside_beta** in {4, 8, 16} (`type: independent`): beta_c,r ~ N(mu0_c,
   1.0^2) independently per region, no index hyperparameters. The fixed-nside ladder of
   fgbuster-style parametric fits with a prior that is flat wherever the data constrain the index
   (Fisher sd per region at nside_beta 16, amplitudes known: median 0.08 for sync with 34 percent
   of regions above 0.1, median 0.013 for dust; at nside_beta 8: 0.040 and 0.006; at 4: 0.019 and
   0.003; full sky 0.0005 and 0.00002). Where a region is unconstrained its posterior is the prior
   (sd 1.0) and the RMSE shows it, which is the honest behaviour of a ladder without shrinkage.
4. **Single index** (nside_beta 0): beta_c ~ N(mu0_c, 0.3^2), as in E017/E018.
mu0_sync = -3.10 and mu0_dust = 1.60: offset from the PySM sky means (-3.000, 1.540) by 1.6 truth
sd, so that shrinkage towards mu0 cannot reproduce the truth by construction; sd 0.3 covers both
truth ranges (-3.28 to -2.74 and 1.26 to 1.79; truth sd 0.063 and 0.037; the s1 index map is the
Miville-Deschenes 2008 model at about 7 degrees, the d1 maps the Planck 2015 Commander fits at 1
degree, see `research/2026-10-27_index_field_recovery_prior_art.md`). Options 1 and 2 learn
their own mean mu_c (hyperprior N(mu0, 0.3^2), weak against the full-sky Fisher sd of 5e-4), so the
offset only matters for options 3 and 4.

Pseudo-templates for the amplitude prior (realistic: built from the data, not from the truth):
g_sync = d(23 GHz) (30/23)^-3.10, g_dust = d(353 GHz) B(545, 20.93 K)/B(353, 20.93 K)
(545/353)^1.60, each smoothed with a 5 degree FWHM Gaussian and floored at its 1st percentile
(floor hits 4.8 percent of the sync and 1.0 percent of the dust pixels, 6.3 and 1.5 percent at
|b| > 20; g/true at |b| > 20, 16/50/84 percentiles: sync 0.41/1.00/1.66, dust 0.85/1.10/1.39).
They carry CMB, noise and the other foreground at the 5 degree scale and a wrong constant index;
mu_c and s_c absorb the scale and scatter (E016/E017 showed hier is calibrated with a 30 percent
large-scale template distortion; this template error is of that order). The same templates serve
every option, so template error cannot separate them.

## Setup
Sky: `pysm3.Sky(nside=64, preset_strings=["d1", "s1", "c1"])` (PySM 3.4.6), emission at the 12
centre frequencies in uK_RJ (PySM evaluates the d1/s1 parameter maps degraded to nside 64; those
degraded maps are the truth we recover against: `pl_index`, `mbb_index`, `mbb_temperature`, and
the c1 CMB map). Noise: white, rms per band as above (0.03 times the E012/E018 nside-16 levels
scaled to nside 64), one noise realization (seed 2027); a second seed for the decisive pair of
runs. Writer: `nextgen/scripts/pysm_intensity_sim.py` (band maps, rms maps, pseudo-templates,
truth.npz with beta_eff, truth amplitudes at the pivots and the per-region Fisher sd, parameter
files, `sim_report.txt`). Mask for the summary statistics: none for the primary numbers; |b| > 20
as a secondary cut; RMSE also per S/N quartile of the pixels.

Runs (Commander4 `nextgen`, 12 MPI ranks, maps-only, one per-pixel hier MCMC group with both
indices, `nuts` max depth 10, target acceptance 0.8, `init_inv_mass: fd_hessian`; 2000 Gibbs
iterations, warm-up 500, burn 600): `gmrf4`, `gmrf8`, `gmrf16`, `iid4`, `iid8`, `iid16`, `indep4`,
`indep8`, `indep16`, `scalar`, in the order gmrf16, iid16, indep16, scalar, gmrf8, iid8, indep8,
then nside_beta 4. Cost: one gradient 30 to 50 ms single process (12 bands, nside 64); the first
timing run without the Hessian start saturated max_depth (1023 gradients, 76 s per iteration) with
the step size pinned at 1e-3 by the hyperparameters; the timing run with `fd_hessian` sets the
budget (see Result).

Metrics (defined before running; `nextgen/scripts/e020_analyze.py`, burn removed):
- **Index error at the data pixels:** per component, beta_hat(p) = posterior mean of the region
  index of pixel p; RMSE over all pixels, over |b| > 20, and per S/N quartile, against the PySM
  truth (both components) and against beta_eff (dust). Uncertainty: the posterior distribution of
  the RMSE of a single draw (mean and sd over draws); two runs tie when their RMSE difference is
  below twice the larger sd. Reference scales: the sd of the truth map (what the single index
  cannot beat), the RMSE of the truth map's own region means, unweighted and S/N-weighted (what a
  region model cannot beat at that nside), and for dust the 0.028 rms of beta_eff - beta_d (what no
  fixed-T model can beat against the PySM map).
- **Index calibration at the region level:** z_r = (posterior mean - truth region mean) / posterior
  sd with the truth region mean unweighted and S/N-weighted; median and MAD over regions (the
  normalized-deviation statistic of Grumitt et al. 2020) and the fractions with |z| < 1 and < 2.
  Reported, not expected to be exactly 0.68 / 0.95: the truth is not drawn from any of the priors.
- **CMB recovery (primary):** the flat-prior GLS CMB map, cmb row of (M^T N^-1 M)^-1 M^T N^-1 d, at
  60 thinned index draws, averaged; minus the input c1 map: RMS over all pixels and |b| > 20, the
  residual power over the input CMB power at l 2 to 30 and 30 to 100, and the regression
  coefficient of the estimate on the input (1 + cross/in) in those bins. References computed once:
  the same GLS at the true indices (the noise floor) and at (beta_s true, beta_eff). Secondary: the
  chain's own posterior-mean CMB alms (shrunk by the N(0, 100^2) prior), same metrics.
- **Fit and cost:** mean reduced chi^2 after burn; divergences; min and median ESS of the index
  regions and of the hyperparameters; gradients per iteration; seconds per iteration (same
  machine, reported, not compared as an algorithmic statement).
- **Hyperparameters:** posterior of (mu, sigma, rho) per GMRF component and (mu, sigma) per iid
  component, against the truth-map mean and sd and the truth map's correlation length.

Decision rule (stated now): the GMRF field is adopted as the Phase 2 index model if, at this noise,
(i) its index RMSE at nside_beta 16 (against the PySM truth for sync, against beta_eff for dust)
is at or below that of every iid and independent run, with differences below two RMSE-draw sd
counted as ties; (ii) its GLS CMB residual RMS and l < 100 residual power are not worse than the
best baseline by more than the difference between the two GLS references; (iii) its RMSE and CMB
residual change between nside_beta 8 and 16 by less than two RMSE-draw sd. If an iid run ties the
GMRF on (i) and (ii), spatial correlation buys nothing at this S/N and the exchangeable
hierarchical prior is the simpler recommendation; if an independent run ties, no hierarchy is
needed (and its best nside, chosen here with knowledge of the truth, is an optimistic baseline a
real analysis cannot select). If the GMRF is worse than the iid or independent runs anywhere, that
is a finding against the stationary smoothness prior on this sky, pointing to the heavy-tailed
variants of Q3. A second noise seed is run for the decisive pair before the conclusion is written.

## Design review
Stats review 2026-10-27 of the first design (independent PySM computation of the sky, Fisher
analysis, code reading): (1) at the first noise level (4 times the E012 nside-16 rms) the sync
Fisher sd per region was 0.65 to 2.7, so every option was prior-dominated and the comparison was
between shrinkage patterns, not priors: noise lowered 33 times (above); (2) the unsmoothed
data-built templates hit the 1e-3 floor in a third of the sky: smoothed to 5 degrees with a
percentile floor; (3) mu0 equal to the truth mean was a truth leak: offset by 1.6 truth sd, and the
iid option added as the nested analogue that learns its mean; (4) the CMB metric on the prior-shrunk
chain amplitudes (shrinkage 0.36 at the old noise) measured the amplitude prior: replaced by the
flat-prior GLS map with the noise-floor reference; (5) truth statistics in the negative control
were wrong (0.11, 0.08): computed by the writer (0.063, 0.037); (6) the fixed T_d puts a 0.028 rms
floor on the dust RMSE against PySM: beta_eff added as the second dust truth; (7) no tolerance in
the decision rule: RMSE-draw sd and a second seed added. Should-fix items taken: floors both ways
and RMSE per S/N bin (8), "no index hyperparameters" wording (10), the gradient budget note (11),
the optimistic-baseline caveat (12); item 9 (a per-region 2-D grid reference for the independent
option, conditional on the amplitude hyperparameters) is a slow L1 test to add in
`tests/test_index_field_l1.py` before the Conclusion is written. Notes 13 to 15 are reflected in
the wording above (noise label, stationarity statement, templates built with the offset indices).

## Reference
Truth: the PySM parameter maps at nside 64, beta_eff, and the c1 CMB map, stored in truth.npz by
the writer. No exact posterior reference exists at L4; the sampler's correctness for every option
is from E018 (GMRF, scalar), from the L1 tests added with the iid and independent options (density
equal to the GMRF likelihood plus the stated prior at matching states to 1e-9, 4th-order finite
differences in every block, parsing and state round trips), and from the per-region grid test of
item 9 above. The `fd_hessian` mass start does not change the target (any positive diagonal mass
is valid); its unit test checks the recovered variances on a Gaussian.

## Command
```bash
python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e020_d1s1 --nside 64 --noise-scale 0.03
for t in gmrf16 iid16 indep16 scalar gmrf8 iid8 indep8 gmrf4 iid4 indep4; do mpirun -n 12 commander4 -p /Users/trygvels/data/e020_d1s1/param_$t.yml; done
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e020_d1s1
```

## Result
(pending)

## Negative control
The single index is the control: its index RMSE must equal the sd of the truth map (0.063 for
beta_s, 0.037 for beta_d at nside 64; the region-mean RMSE floors are 0.037 / 0.024 / 0.013 for
beta_s and 0.028 / 0.024 / 0.019 for beta_d at nside_beta 4 / 8 / 16) up to its own posterior
sd, and its CMB residual must be the largest. A second control: `indep16` with 16 data pixels per region should be the noisiest index
map (posterior sd per region 3 to 10 times the GMRF's where the GMRF shrinks); if it is not, the
GMRF prior is not doing anything on this sky.

## Conclusion
(pending)

## Follow-ups
- Noise scale 3 and a polarization (QU) version once M2.4 stage 2 allows hier for QU streams.
- The remaining Q3 variants: alm GRF with sampled C_l, heavy-tailed increments, clustering.
- d10s5 and d12s7 skies; beams and bandpasses (L4 realism).
