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
from 14 to 31 K, so the dust model is misspecified by construction: the noise-free fixed-T
best-fit index differs from the PySM beta_d by rms 0.028 when fitted to the dust emission alone
(beta_eff) and by rms 0.041 when the cmb and sync amplitudes are profiled as well (beta_eff3, the
index the three-component fit converges to; the two differ by rms 0.012); beta_eff3 is the dust
truth a one-MBB fixed-T model can at best recover. All three are stored and reported).

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
   [h, pi] through the logit-transformed u with its Jacobian. Exactly the E016/E018 prior. Given
   the hyperparameters, a region index has prior N(mu_c, sigma_c^2 var_r(rho_c)) with sigma and
   rho learned from the whole sky (mostly the plane) and applied everywhere (stationarity); that,
   and the coupling of neighbours, is what is being tested.
2. **Exchangeable (iid) hierarchical field at nside_beta** in {4, 8, 16} (`spectral_index_field:
   {type: iid}`): beta_c,r = mu_c + sigma_c z_r, z ~ N(0, I), with the same hyperpriors on mu_c and
   sigma_c as option 1 and no spatial coupling: the rho -> 0 limit of the GMRF up to the constant
   mode (the iid block keeps it, so the field mean is mu_c + sigma_c zbar with prior variance
   0.3^2 + sigma_c^2 / n_beta, at most 0.5 percent more than the GMRF's at nside_beta 4; the
   analyzer reports the field mean next to mu_c). The GMRF's u is kept as a decoupled nuisance
   drawn from its own prior. The nested analogue: option 1 minus spatial correlation. The scalar
   option 4 is its sigma -> 0 limit.
3. **Independent regions at nside_beta** in {4, 8, 16} (`type: independent`): beta_c,r ~ N(mu0_c,
   1.0^2) independently per region, no index hyperparameters. The fixed-nside ladder of
   fgbuster-style parametric fits with a prior that is flat wherever the data constrain the index
   (per-region Fisher sd with the three per-pixel amplitudes profiled out, from the writer's
   `sim_report.txt`: nside_beta 4 median 0.041 sync / 0.024 dust with 22 / 36 percent of regions
   above the truth sd of 0.063 / 0.037; nside_beta 8 median 0.087 / 0.054, 68 / 61 percent above;
   nside_beta 16 median 0.18 / 0.12, 89 / 76 percent above; amplitudes known the medians are 2 to
   9 times smaller, 0.019 / 0.003, 0.040 / 0.006, 0.081 / 0.013; full sky 0.0005 / 0.00002). So
   nside_beta 4 is data-dominated for most regions, 8 is marginal and 16 is prior-dominated for
   most regions: the ladder at 16 returns the prior (sd 1.0) in most regions and its RMSE shows
   it, the honest behaviour of a ladder without shrinkage; the 16 versus 4 comparison is a
   comparison of shrinkage, which is the point of the hierarchical priors.
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
indices, `nuts` max depth 10, target acceptance 0.8, `init_inv_mass: fd_hessian` written into
every parameter file by the writer; 2000 Gibbs
iterations, warm-up 500, burn 600): `gmrf4`, `gmrf8`, `gmrf16`, `iid4`, `iid8`, `iid16`, `indep4`,
`indep8`, `indep16`, `scalar`, in the order gmrf16, iid16, indep16, scalar, gmrf8, iid8, indep8,
then nside_beta 4. Cost: one gradient 30 to 50 ms single process (12 bands, nside 64); the first
timing run without the Hessian start saturated max_depth (1023 gradients, 76 s per iteration) with
the step size pinned at 1e-3 by the hyperparameters; the timing run with `fd_hessian` sets the
budget (see Result).

Metrics (defined before running; `nextgen/scripts/e020_analyze.py`, burn removed):
- **Index error at the data pixels:** per component, beta_hat(p) = posterior mean of the region
  index of pixel p; RMSE over all pixels, over |b| > 20, and per quartile of the per-pixel index
  Fisher weight (a dsed/dbeta / rms summed in quadrature over bands), against the PySM truth (both
  components) and against beta_eff and beta_eff3 (dust). Uncertainty for comparing two runs:
  paired block bootstrap of the RMSE difference over the 192 nside-4 blocks of the sky (the same
  blocks in both runs; blocks are the near-independent units); two runs tie when |difference| is
  below two bootstrap SE. The posterior spread of the RMSE of a single draw is reported as a
  descriptive number only. Reference scales: the sd of the truth map (what the single index
  cannot beat), the RMSE of the truth map's own region means, unweighted and Fisher-weighted
  (what a region model cannot beat at that nside), and for dust the 0.041 rms of beta_eff3 -
  beta_d (what no fixed-T model can beat against the PySM map).
- **Index calibration at the region level:** z_r = (posterior mean - truth region mean) / posterior
  sd with the truth region mean unweighted and Fisher-weighted; median and MAD over regions (the
  normalized-deviation statistic of Grumitt et al. 2020) and the fractions with |z| < 1 and < 2.
  Reported, not expected to be exactly 0.68 / 0.95: the truth is not drawn from any of the priors.
- **CMB recovery (primary):** the flat-prior GLS CMB map, cmb row of (M^T N^-1 M)^-1 M^T N^-1 d, at
  60 thinned index draws, averaged (a plug-in posterior-predictive estimate: index uncertainty
  enters through the average, so a wide index posterior is penalized, which is intended); minus
  the input c1 map: RMS over all pixels and |b| > 20, the (2l+1)-weighted residual power over the
  input CMB power at l 2 to 30 and 30 to 100, the regression coefficient of the estimate on the
  input (1 + cross/in) in those bins, and the spread of the per-draw RMS. References computed once:
  the same GLS with the true dust SED per pixel (`dust_bands / a_dust`; the noise floor no index
  model can beat), and with the fixed-T dust column at beta_d, beta_eff and beta_eff3 (which carry
  the dust misfit). Run comparisons by the same paired block bootstrap on the squared residual.
  Secondary: the chain's own posterior-mean CMB alms (shrunk by the N(0, 100^2) prior).
- **Fit and cost:** mean reduced chi^2 after burn; divergences; min and median ESS of the index
  regions, of the field hyperparameters and of the amplitude hyperparameters; gradients per
  iteration; seconds per iteration (same machine, reported, not compared as an algorithmic
  statement).
- **Hyperparameters:** posterior of (mu, sigma, rho) per GMRF component and (mu, sigma) per iid
  component, against the truth-map mean and sd and the truth map's correlation length.

Decision rule (stated now; "better" and "tie" by the paired block bootstrap, two SE):
(a) **Same-nside comparisons** (the prior at fixed region size): at each nside_beta in {4, 8, 16},
GMRF against iid against independent on sync RMSE (PySM truth), dust RMSE (beta_eff3) and GLS CMB
residual. At nside_beta 4 (data-dominated) all three are expected to tie; a GMRF loss there is a
finding against the stationary smoothness prior on this sky. At 8 and 16 (prior-dominated) the
hierarchical priors are expected to beat the independent ladder; GMRF against iid at 16 is the
test of spatial correlation itself.
(b) **The Phase 2 question** (does the field let us go to the finest nside): the GMRF at
nside_beta 16 is adopted as the Phase 2 index model if its sync RMSE, dust RMSE and GLS CMB
residual are each at or below (tie allowed) the best independent run at any nside (the ladder's
best nside is chosen here with knowledge of the truth, an optimistic baseline a real analysis
cannot select) and at or below the best iid run, and (c) its RMSE and CMB residual at nside_beta 8
and 16 tie (the inferred smoothness makes nside_beta a nuisance). If an iid run ties the GMRF on
every metric, spatial correlation buys nothing at this S/N and the exchangeable prior is the
simpler recommendation; if the best independent run ties, no hierarchy is needed. Expectation to
write the Conclusion against (review 2, note 14): the s1 index map is smooth at about 7 degrees,
so a GMRF win on sync is partly the smoothness of the truth; the d1 map has 1-degree scatter, so
a smaller or no GMRF gain on dust is expected and is not a failure of the method. A second noise
seed (data and pseudo-templates) is run for the decisive pair before the conclusion is written.

## Design review
Second stats review 2026-10-27 of the revised design (independent re-derivation of the four
targets from the code, PySM recomputation of the sky, GLS check on the noiseless sky to 3e-5 uK,
region-map convention check): targets and metrics confirmed; three must-fix items taken before
running: the parameter files omitted `init_inv_mass` (writer fixed); the CMB "floor" used the
fixed-T dust column and so was 6.4 uK against 2.8 uK with the true dust SED, making rule (ii)
unfalsifiable (true-SED floor added, tolerance from the paired bootstrap); the Fisher sd with the
amplitudes profiled is 2 to 9 times the amplitudes-known value, so nside_beta 16 is
prior-dominated (numbers above; decision rule split into same-nside and cross-nside parts).
Should-fix items taken: three-component beta_eff3 (4), paired block bootstrap instead of the
RMSE-draw sd (5), Fisher weights for the weighted region means and the bins (6), iid constant-mode
wording and field-mean report (7), field type written into the chain state so iid runs are not
reported with a rho (8), conditional-prior wording (9); notes 12 and 13 (plug-in averaging over
draws, (2l+1) weights) reflected above; note 15 (ESS of the amplitude hyperparameters) added.

First stats review 2026-10-27 of the first design (independent PySM computation of the sky, Fisher
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
Second seed for the decisive pair: `--dir /Users/trygvels/data/e020_d1s1_seed2 --seed 2028`.

## Result
**Timing run 1, gmrf16, non-centred, fd_hessian (40 Gibbs iterations, 06:53 to 08:44):** not
usable. The finite-difference start costs 12304 gradients (18 min per chain; Commander4 always
runs two chains, so every run pays it twice) and gives an inverse-mass range 1.5e-8 to 2. The
step size then collapses in the first ten transitions to 2e-5 to 1e-3 (chain 2 stays at 3e-5 to
8e-5, chain 1 oscillates up to 3e-3), with 31 to 1023 gradients per transition and 40 to 90 s
per iteration. The analyzer, run on the 35 post-burn iterations of chain 1 (this also exercised
its HDF5 path for the first time; it ran without error), shows the state frozen: both field means
equal to the prior start to four decimals, posterior sd per region 2e-4 to 4e-4, index RMSE 0.118
(sync) and 0.071 (dust), the amplitude hyperparameters unchanged at mu 1, log s -1.2, and the
GLS CMB residual 12.6 uK. Cause: with 0.03 x E012 noise, the profiled Fisher sd of a region index
at nside_beta 16 spans 4e-3 to 0.61 (sync) and 5e-4 to 1.16 (dust), a precision range of about
5e6, and in the non-centred coordinates z these strongly informed regions are tied to mu and
sigma, which a diagonal mass cannot undo (the usual non-centred versus centred trade-off:
non-centred is the right choice for weak data, E016 to E018, and the wrong one here). The same
Fisher sd range is itself a finding: at this noise one third of the nside_beta 16 regions have an
index precision better than the truth's own spread (0.063 / 0.037), the rest are prior-dominated.

**Timing run 2, gmrf16, centred, fd_hessian (30 iterations, 08:47 to 10:26):** the chain moves.
Step size 2e-3 to 3e-2 throughout (three divergent transitions in the first ten, none after),
gradients per transition 1 to 1023 (chain 1 saturates max_depth 10 from iteration 15 on, chain 2
stays at 3 to 255; both still in the first warm-up window, so the windowed mass adaptation has
not acted yet). The analyzer on the 20 post-burn iterations: region z median/MAD -0.30/1.13
(sync) and +1.86/1.51 (dust), posterior sd 0.13/0.03, both field means and all four amplitude
hyperparameters moving, index RMSE 0.34/0.09, GLS CMB residual 7.1 uK. These are warm-up numbers
and prove only that the centred parameterization is usable here; the production runs use
`parameterization: centred` (a sampler coordinate choice, same target, E018 L1 covers both
parameterizations).

**Profile of a Gibbs iteration (from the centred log):** 36 s of every iteration went into
rebuilding the graph-Laplacian eigenbasis (dense eigh of 3072^2, 18 s per component, on every
rank, every iteration); fixed by caching the basis per process (commit cba5c33, same target).
The finite-difference start costs 12304 gradients per chain, 21 min at 2 solver threads, that is
about 98 ms per gradient in the MPI run (two chains per run, so 42 min fixed cost per run).

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
