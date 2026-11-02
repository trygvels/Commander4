# E026: Stage-2 chains on the realistic-noise PySM sky: hier prior with a region-scale template multiplier field against the scalar law, component-separation-quality templates

- **Date started / finished:** 2026-10-30 / 2026-11-02
- **Branch and commit:** nextgen (mu-map prior under stats review before the runs)
- **Open question addressed:** Q16 (D022, D023), Q3 / M2.5 (first chains at realistic noise), M2.4 (dense metric at nside 16 with the reviewed code)
- **Validation level:** L3 (misspecified simulation with known truth; the sky is not from the model class: indices vary inside regions, template errors are log-normal and smooth)

## Question
Do the E025 exact-target results carry over to the posterior under the GMRF index prior: with
component-separation-quality templates (true amplitude times a 10-deg, 30 percent log-normal
error field) is the gmrf16 posterior of the index fields at realistic noise calibrated against
the noise-free region fit when the hier prior carries a template multiplier per nside-16
region (`mu_nside: 16`), and biased when it carries the scalar law; and how close does the
mu-map run get to the attainable reference, the same chain with exact templates? This is the
first index-field chain at realistic noise, with the flat CMB prior (D021) and the dense
Hessian metric (D020) at nside 16 (dim about 12,300 with the mu maps).

## Target and priors
The E021 gmrf16 target on the E023 noise-1 sky (`/Users/trygvels/data/e023_pysm_n1`: PySM
d1s1c1, nside 64, 12 bands 23 to 353 GHz, white noise at the E012 levels) with: CMB amplitude
prior N(0, 1e4^2) (flat in effect, D021); sync and dust hier priors with the templates below,
`mu_prior [1, 0.5]` per mu entry, `log_s_prior [log 0.3, 1]`; GMRF index fields at nside 16
(centred, graph-Laplacian basis, alpha 2, learned mu, sigma, rho; index prior means -3.1 and
1.6 with sd 0.3 on the field mean); dust at the fixed T 20.93 K (E025: the T misfit changes no
shift by more than 0.01 at this depth). Three runs, 1200 Gibbs iterations, warm-up 400, burn
500, dense metric refreshed at 50 and 200:
- gmrf16_cs_mumap: templates `template_{sync,dust}_cs_*.fits` (true amplitude times exp(f), f
  the E025 field, seed 11), `mu_nside: 16` on both hier priors (3072 multipliers each).
- gmrf16_cs_scalar: the same templates, scalar mu (the model E023/E025 showed biased).
- gmrf16_exact_mumap: templates equal to the true amplitudes, `mu_nside: 16`: the attainable
  reference at this depth (its mu map should sit at 1 and its s at the guard).

## Setup
Parameter files `param_gmrf16_{cs_mumap,cs_scalar,exact_mumap}.yml` in the data directory
(from `param_gmrf16.yml` with the template paths, CMB rms 1e4, `init_inv_mass:
fd_hessian_dense`, `mu_nside`), 12 MPI ranks on the Mac mini, one run at a time (about 10 h
each at the E021 gmrf16 rate; the dense Hessian at dim 12,300 costs 24,600 gradients per
refresh plus an eigendecomposition of a 1.2 GB matrix, to be timed). Analysis:
`e020_analyze.py` (index RMSE against the PySM truth and against beta_eff3, GLS CMB residual,
z and coverage against the S/N-weighted truth region mean; the mu map enters its amplitude
trace as its mean) plus a dedicated script for: index RMSE and z against the noise-free
region-fit reference of E023 (`E023/scan_n1_refs.npz`, fixT), per template-S/N quartile; the mu
map against the realized region means of exp(-f); s against the within-region residual; ESS
per gradient of every block (the dense metric at this dimension is itself under test).

## Reference
Truth: the PySM index maps and the E023 region-fit reference (fixT). Pre-registered reading:
(1) calibration: in gmrf16_cs_mumap the z of the region indices against the region-fit
reference has |median| below 0.3 and MAD within 0.5 to 1.5 in every S/N quartile where the
block sd is below 0.05 (dust Q1 to Q3, sync Q1), and no region has |beta_s - reference| above
0.5; (2) bias where the scans predicted it: gmrf16_cs_scalar has a dust z median above +0.5 in
Q2 and Q3 (E025: shifts beyond 0.1 at block sds 0.03 to 0.05) and the sync Q1 z median below
-0.3; (3) attainable: the cs_mumap index RMSE per quartile is within a factor 1.5 of
exact_mumap's; (4) hyperparameters: the cs_mumap mu maps track the realized region means of
exp(-f) (correlation above 0.8 over the regions of the two brightest quartiles) and s sits
near the residual 0.06 (E025); (5) CMB: the GLS CMB residual of cs_mumap within 10 percent of
exact_mumap's. Efficiency (secondary): ESS per gradient of the mu blocks and of the field
means; the dense metric's eigendecomposition time at dim 12,300. If (1) fails while (3) holds,
the posterior widths (the index prior) are the issue, not the amplitude prior; if (1) and (3)
fail, the E025 scans missed something that only the joint posterior shows and the record says
so.

## Command
```bash
for t in gmrf16_cs_mumap gmrf16_cs_scalar gmrf16_exact_mumap; do mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_$t.yml; done
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e023_pysm_n1 --burn 500 --tags gmrf16_cs_mumap gmrf16_cs_scalar gmrf16_exact_mumap
```

## Result
**Centred runs (2026-10-30 04:38 to 2026-10-31 13:17; 7.5, 9.3 and 15.8 h; 1200 iterations each,
700 post-burn; no tracebacks; `E026/e026_analysis.txt`, `e026_e020_analysis.txt`,
`hyper_ess_centred.txt`).** The dense metric at dim 12,298 (two mu maps) swept 6150 coordinates
(12,300 gradients, 15 min) with the 6144 mu coordinates analytic; 2 to 4 negative-curvature
directions at every (re)initialization; 31 / 112 / 230 gradients per iteration post-burn.

The chains are not valid posterior samples of their targets: the GMRF field hyperparameters did
not mix. Sync log sigma: ESS 9 (cs_mumap, sigma 0.0013), 4 (cs_scalar, 0.014), 4 (exact, 0.0020);
dust log sigma ESS 13 / 291 / 12; u ESS 17 / 30 / 4 (cs_mumap, cs_scalar, exact) and 6 / 132 / 5
for dust. The sync field sigma went from its start (0.15, INIT_SIGMA_FRACTION x sigma_scale) to
the funnel neck during warm-up and stayed there (sigma 0.001 to 0.002 against a truth of 0.063),
so the sync field is constant over the sky at the field mean (-2.897 cs_mumap, -2.949 exact,
-3.166 cs_scalar) and every region's posterior sd is the same 0.027 to 0.028; the dust field
has sigma 0.15 / 0.38 / 0.06 with rho 179 / 44 / 93 deg, a near-constant field at 1.60 / 1.74 /
1.58 against a sky mean of 1.54. This is the funnel of E018 (F20: centred NUTS 5 to 240 times
worse on the funnel) at a depth where the per-region data are weak (region Fisher sds 1 to 3.7
for sync, 0.06 to 1.5 for dust), with the dense metric refreshed at transitions 50 and 200 at
states already in the neck. The scans of E023 and E025 did not see this because they fix the
field at the reference and carry no index prior.

What the three chains show, read as diagnostics only (`e026_analyze.py`; z against the region-fit
reference; sync sd median 0.028 in every run):
- gmrf16_exact_mumap (exact templates, mu map: the attainable reference): sync z median +0.75,
  +1.4, +1.5, +1.4 (Q1 to Q4; the constant field 0.05 above the median region fit, 1.9 sd of
  its own sd), dust z +2.3, +1.8, +0.8, +0.4 (RMSE 0.078, 0.083, 0.063, 0.049), sync RMSE 0.12,
  0.078, 0.072, 0.070. The dust mu map sits at 1.028 +- 0.007 (should be 1.000; the fixed-T
  amplitude-index degeneracy, the dust field 0.04 above the sky mean), s_d 0.022, s_s 0.039; the
  mu maps are uncorrelated with the (absent) template error, as they should be (correlation
  0.02 to 0.04).
- gmrf16_cs_mumap: sync z median +2.6 to +3.3 (the constant field 0.09 above the region fits;
  RMSE 0.14, 0.11, 0.11, 0.11), dust z +1.3, +1.3, +1.0, +0.6 (RMSE 0.071, 0.109, 0.114, 0.100);
  mu maps track the realized exp(-f) with correlation 0.72 (sync Q1) and 0.89 (dust Q1), 0.2 to
  0.4 in the fainter quartiles, where the mu posteriors are 0.3 to 0.47 wide; s 0.057 / 0.068.
  Against exact_mumap: paired pixel-level +0.040 +- 0.001 (sync), +0.030 +- 0.004 (dust).
- gmrf16_cs_scalar (negative control): sync z median -4.2 to -4.8 (field at -3.166, amp mu 0.764,
  s 0.21: the template's mean error absorbed by the amplitude mean and the index), dust z +0.9 to
  +1.25 with RMSE 0.27 to 0.38 and sigma_d 0.38 (the field tries to follow the template error);
  29 divergences; paired against exact +0.096 / +0.260. The predicted direction and size of the
  scalar-law bias, on an invalid chain.
- GLS CMB residual 92 uK in all three and 93 uK for the true-SED GLS: at this depth the per-pixel
  flat-prior GLS CMB is noise; the CMB metric is uninformative at realistic noise and the
  pre-registered reading (5) cannot be evaluated this way.

Pre-registered reading: (1) fails and (3) fails, but for a reason the design did not
anticipate (an unmixed chain), so neither "the posterior widths are the issue" nor "the scans
missed something" can be concluded; (2) the negative control shows the predicted bias; (4) the
mu maps track the template error where the data allow (bright quartiles) and s sits at the
E025 residual. Decision (D024): rerun the three chains in the non-centred parameterization
(E018: it reaches and leaves the small-sigma tail), exact templates first, before any reading
of calibration; the centred chains are kept as the record of the failure mode.

**Non-centred runs:** started 2026-10-31 13:28 (`param_*_nc.yml`, output `chains_*_nc/`), exact
first. Early diagnostic at 150 iterations of gmrf16_exact_mumap_nc (15:20): the field
hyperparameters move. Sync log sigma between -2.6 and -1.4 (sigma 0.07 to 0.21 over the last
50 iterations against a truth of 0.063), u between -0.6 and 3.0, field mean -2.80 to -2.90;
dust sigma 0.06 to 0.12, u 0.25 to 2.4, mean 1.56 to 1.61; 15 to 63 leapfrogs per transition
(one at 1023). The dense metric at the non-centred start had 4 negative directions, 2 at the
transition-50 refresh. The chain is out of the neck.

**Non-centred runs complete (13:28 to 2026-11-01 06:48; 7.5, 6.4 and 3.4 h; 1200 iterations,
700 post-burn; `E026/e026nc_analysis.txt`, `e026nc_e020_analysis.txt`,
`hyper_ess_noncentred.txt`).** Hyperparameters mix: log sigma ESS 77 / 24 / 143 (sync) and
60 / 6 / 38 (dust) for exact_mumap, cs_mumap, cs_scalar; u 48 / 13 / 109 and 43 / 5 / 25; field
means 120 / 43 / 570 and 131 / 34 / 169; region index ESS medians 234 / 339 / 1072 (sync) and
136 / 152 / 972 (dust), minima 29 / 27 / 315 and 21 / 14 / 323; 63 / 31 / 31 gradients per
iteration. cs_mumap's dust hyperparameters (ESS 5 to 6) are the weakest; the others are usable.
Fields: exact sigma_s 0.152 +- 0.047 (sky sd 0.063), rho 59 +- 39 deg, mean -2.868 +- 0.042;
sigma_d 0.062 +- 0.009, rho 40 +- 11 deg, mean 1.578 +- 0.012 (sky mean 1.540); cs_mumap
sigma_s 0.35, mean -2.72; sigma_d 0.12, mean 1.585; cs_scalar sigma_s 0.87 (!), mean -3.31,
sigma_d 0.37, mean 1.72. Amplitude hyperparameters: exact mu_s map mean 0.997, mu_d map mean
1.029 +- 0.007 (should be 1.000), s 0.034 / 0.022; cs_mumap mu maps 0.971 / 1.030 with
correlation to the realized exp(-f) 0.71 / 0.89 in the brightest quartiles and 0.2 to 0.4 in
the others (posterior sds 0.25 to 0.48 there), s 0.055 / 0.068 (E025: 0.06); cs_scalar amp mu
0.78 / 1.07, s 0.084 / 0.166.

Against the region-fit reference per template-S/N quartile (Q1 bright to Q4 faint; z median /
|z| < 1 / RMSE; sync posterior sd medians 0.13 to 0.15 for exact, 0.22 to 0.29 for cs_mumap,
0.50 to 0.72 for cs_scalar; dust 0.04 to 0.06, 0.055 to 0.10, 0.08 to 0.27):
- exact_mumap: sync +0.50 / 0.83 / 0.10, +0.82 / 0.65 / 0.13, +0.90 / 0.60 / 0.15, +0.94 / 0.56 /
  0.16 (MAD 0.25 to 0.40: the whole field sits about +0.9 of its sd, +0.13, above the
  reference, with a per-region scatter smaller than its sd); dust +1.44 / 0.30 / 0.072, +1.34 /
  0.28 / 0.080, +0.72 / 0.63 / 0.063, +0.17 / 0.82 / 0.047.
- cs_mumap: sync +0.51 / 0.76 / 0.19, +0.95 / 0.53 / 0.30, +1.07 / 0.45 / 0.37, +1.19 / 0.35 /
  0.41; dust +0.90 / 0.51 / 0.069, +0.88 / 0.54 / 0.102, +0.58 / 0.66 / 0.099, +0.26 / 0.78 /
  0.080. Paired against exact: pixel-level +0.188 +- 0.012 (sync), +0.022 +- 0.003 (dust).
- cs_scalar: sync -1.02 / 0.40 / 0.83, -0.72 / 0.64 / 0.65, -0.42 / 0.84 / 0.50, -0.24 / 0.86 /
  0.50; dust +0.74 / 0.24 / 0.27, +0.99 / 0.40 / 0.29, +1.14 / 0.43 / 0.33, +1.16 / 0.40 / 0.35.
  Paired against exact +0.476 +- 0.025 and +0.244 +- 0.011.
GLS CMB residual 94 to 95 uK in all three against a true-SED floor of 93: uninformative at
this depth, as in the centred runs.

Two features of the exact-template reference need an exact-target check before the reading:
the sync field's aggregate sits 0.13 (3 of its own sd of 0.042) above the sky's mean index
(and above every quartile's reference), and the bright dust quartiles sit +0.06 to +0.07 (1.4
sd per region, systematic) with the dust mu map at 1.029. A whole-sky shift scan of the exact
target on this realization (expected and realization modes, global and per-region mu) is
running (`e026_allsky_scan.txt`). A second finding on the way: the template-S/N quartiles are
built from the data-built 23 GHz template, so on this realization the faint sync quartile is
where the 23 GHz noise is low (mean noise -12 standard errors in that quartile against -0.7
over the sky), which is why the realization-mode faint-sync scans of E023 and E025 preferred a
flatter index by up to +1.4 at 68 nats; the expected-mode scans (noiseless sky) and the chains
(which know no quartiles) are not affected, but every per-quartile statistic on a single
realization carries this selection, on real data too when regions are selected by a template
built from the same data.

**Exact-target checks of the exact-template offsets (2026-11-02; `E026/allsky_scan_exact.txt`,
`state_check_exact_mumap_nc.txt`, `nextgen/scripts/e026_state_check.py`).** (1) Whole-sky
shift scans of the exact target with the index fields held at the region-fit reference (flat
CMB prior; global or per-region mu profiled; expected and realization modes): sync -0.02
(curvature sd 0.029, gain 0.2 to 0.5 nats) and dust 0.000 (sd 0.0037) in every cell except the
realization mode with per-region mu, where the sync optimum is -0.10 (7 nats; the 3072 free
multipliers absorb this realization's noise). The data term at the reference prefers neither
the +0.13 sync nor the +0.04 dust offset. (2) The full target (likelihood, GMRF priors at the
chain's posterior-mean hyperparameters, hier priors with the posterior-mean mu maps) at the
chain's posterior-mean state against the same state with the fields at the reference: the
chain's state is 2696 nats higher (1807 from the dust field, 887 from the sync field); setting
the posterior-mean dust mu map (1.0286) to 1 costs 1191 nats, moving the dust field down by
0.04 costs 375; moving the sync field mean by -0.065 gains 10 nats and by +0.065 costs 23. So
the chain is where its target puts the mass: the offsets are properties of the model at this
depth, not of the sampler, and they are not visible to block-shift scans at a reference field.
Reading: (a) dust: with the exact template at the 545 GHz reference and a fixed T, the SED is
misspecified pixel by pixel; the flat-amplitude reference absorbs that per pixel, the hier
model cannot (s 0.022), so the per-region multipliers (+2.9 percent) and the bright-quartile
index (+0.06) absorb it jointly, 1191 + 375 nats better than the reference; the E025 scans with
one global mu could not take that path (+0.5 percent) and so did not see it. The region-scale
multiplier field makes the fixed-T misfit a bias where a scalar mu kept it a misfit. (b) sync:
every region is prior-dominated (region Fisher sd 1 to 3.7; posterior sd 0.13 to 0.15 set by
the field scale, itself prior-dominated at 0.14 +- 0.05 against a sky sd of 0.063); the data
constrain only the aggregate, and the aggregate posterior mean -2.868 +- 0.042 sits 3 sd above
the data-term optimum -3.00 +- 0.03. A run with a single sync index on the same data (dust
field and mu maps unchanged, `chains_exact_syncscalar_mumap_nc/`, started 13:54) tests whether
the field model's extra freedom produces that offset.

**Where the +0.13 sync aggregate comes from (2026-11-02 16:44 to 17:30; `E026/sync_level_check.txt`,
`nextgen/scripts/e026_sync_level_check.py`; chain-state scan `E026/chainstate_scan_exact.txt`).**
(1) The single-sync-index run (`chains_exact_syncscalar_mumap_nc`, 1200 iterations, 0 post-burn
divergences, dust field and mu maps as in the field run: dust field mean 1.5785, mu_d map 1.028)
gives beta_s = -2.944 +- 0.031 (ESS 130), +0.039 (1.3 sd) above its pseudo-true value, the
Fisher-weighted truth mean -2.983 (the sky's unweighted mean -3.000 is not what a scalar fit
estimates). No significant offset for a single index. (2) The field run's aggregate posterior
is symmetric (skew 0.04 for mu and for the field mean; per-region skew 0.08), so the offset is
not the mean of a skewed marginal. (3) The bright quartile holds 90 percent of the sync Fisher
weight, and there the field agrees with the region-fit reference: the 10 / 30 / 100 / 300
regions with the largest Fisher weight (31 / 50 / 67 / 82 percent of it; region Fisher sds 0.10
to 0.55, posterior sds 0.08 to 0.12) sit -0.02 / +0.03 / +0.03 / +0.03 from the reference with
z medians -0.6 / +0.2 / 0.0 / +0.2 and 77 to 90 percent inside |z| < 1. The Fisher-weighted
level of the posterior-mean field is -2.978 against a Fisher-weighted truth of -2.988 and
reference of -3.009. (4) The unweighted field mean (mu, the constant mode the basis excludes)
is +0.13 because the three faint quartiles (10 percent of the Fisher weight, region Fisher sds
2.5 to 3.7, posterior sds 0.15 set by the field scale) float at +0.12 to +0.14 above their
references, while the bright quartile sits at -0.034 below mu. So at this depth the field's
unweighted mean is not a data-determined quantity; the data determine the Fisher-weighted
level, and that is right to 0.01 (truth) and 0.03 (one sd of the full-sky constraint,
reference). The E026 statement "sync aggregate 3 sd high" compared a prior-dominated quantity
(mu) with the truth; the sd 0.042 is the posterior sd of the floating level, not of a data
constraint, and the comparison that matters is (3). What sets the faint level at +0.13 rather
than at the truth is checked by a whole-sky shift scan of the exact target around the chain's
posterior-mean state (`--ref-npz`): with one global multiplier the data term prefers a shift
of -0.14 (12 nats, curvature sd 0.030) in both modes, which cannot be read directly because the
chain's law has 3072 multipliers; the per-region-multiplier cells are recorded below.

**Chain-state scan, Gibbs-consistency check and the roughness coupling (2026-11-02 16:50 to
17:20; `E026/chainstate_scan_exact.txt`, `conditional_check.txt`, `roughness_check.txt`;
`nextgen/scripts/e026_conditional_check.py`, `e026_roughness_check.py`).** (5) The exact-target
scan around the chain's posterior-mean state prefers a whole-sky sync shift of -0.14 (global
law, 12 nats), -0.12 (per-region multipliers, expected mode, 8 nats) and -0.18 (per-region
multipliers, realization mode, 23 nats; curvature sd 0.028). (6) Gibbs-consistency check of
the chain: for 20 post-burn draws (every 35th) the full Commander4 target was scanned along
the field mean of one component with every other coordinate at the draw's value; the
conditional CDF at the draw is consistent with uniform (sync: mean 0.43 +- 0.065, KS p 0.15;
dust: 0.40, KS p 0.055), and the conditional optimum sits at the draw within +0.005 +- 0.014
(sync; conditional sd 0.018) and +0.0007 +- 0.0014 (dust; sd 0.0014). The chain samples its
target along the mean directions; the sampler is not the cause. (7) Why (5) and (6) differ:
the data-preferred common level depends on the field's roughness. For six draws, scaling the
non-centred coefficients z by 1, 0.5, 0.25 and 0 (the draw's field, sd 0.08 to 0.17 over the
sky, down to a constant field at mu) moves the block-shift optimum from -0.025 to +0.026
(at the draw) to -0.065 to -0.130 (constant field), and at the posterior-mean state (z and
multipliers averaged, a smooth field) from -0.072 to -0.180 (z x 0). The scans of E023, E025
and the whole-sky scans above were all made at smooth fields (the region-fit reference or the
posterior mean) and so measured the level of a smooth field, which is 0.1 steeper than the
level the posterior gives a field of the sampled roughness. The likelihood is nonlinear in
the index over the 0.15 the prior-dominated field fluctuates by, and its optimum for the
common level shifts with that fluctuation; the field scale sigma_s = 0.15 +- 0.05 is itself
prior-dominated (sky 0.063), so the level of the faint sky is set by the prior on sigma
through this coupling, not by the data. The dust field (sigma 0.062 +- 0.009, data-constrained)
shows no such offset: its conditional optima sit at the draws.

**Independent stats review of items 1 to 7 (2026-11-02 19:40; `E026/stats_review_sync_level.md`,
`E026/review_scripts/`).** (8) The Gibbs-consistency check is valid and passes for sync; dust is
borderline (optimum minus draw +0.0007 at 2.3 se, step too coarse); a K = 100 run with a finer
dust step is in progress. (9) The roughness coupling of item 7 is withdrawn: the optimum shift
is linear in the z scale (a first-order effect of the data-pinned bright regions being pulled
off their optima; z x 0 is the scalar model), and with the multipliers integrated analytically
the faint-region data level is -2.79 +- 0.065 for constant, reference, chain, shuffled and iid
patterns alike. The cause of the +0.13 is the marginalization over the 3072 per-region
multipliers with hyperprior N(1, 0.5^2): its volume term -1/2 sum_r log(-Q_r) rises with a
steeper index (slope +130 nats per unit beta over the faint regions; 66 percent of the sync
information sits at 23 GHz), so the exact marginal of the faint level is +0.15 +- 0.05 above the
truth over simulated noise realizations (+0.19 here); the chain's +0.13 is that bias held back
by the field prior. With hyperprior sd 0.2 or 0.1 the bias is +0.06 or zero within noise
(E029 tests this in a chain). (10) The scans with profiled per-region multipliers (E025, E026
items 1 and 5) carry a Neyman-Scott bias toward steeper indices (-0.2 to -0.4), and the
global-multiplier cells are a different model: with the multipliers integrated, the whole-sky
optimum at the posterior-mean field is -0.003 +- 0.025, i.e. the chain sits at its target's
optimum. F28 (scan-based) is flagged for a recheck with the analytic marginal. (11) The scalar
run's pseudo-true value under the model is -2.956 +- 0.003 and the exact 1-d marginal for this
realization -2.951 +- 0.025: the chain (-2.944) agrees with an exact marginal; the earlier
reference -2.983 (amplitude-known Fisher weights) was the wrong yardstick. With the dust at the
truth and mu_d at 1 the scalar pseudo-true moves to -2.83: the sync index is sensitive at the
0.1 level to the fixed-T dust state.

## Negative control
gmrf16_cs_scalar (the scalar law the scans showed biased) must show the predicted dust and sync
biases; if it does not, either the scans or the chains are wrong and the record says which.

## Conclusion
1. **Valid chains at realistic noise need the non-centred parameterization** (D024): the centred
   chains stuck in the funnel neck; the non-centred ones mix (hyperparameter ESS 24 to 143, 2
   divergences in 700) and pass a Gibbs-consistency check along both field means (conditional
   CDF uniform, KS p 0.15 and 0.055; conditional optima at the draws). The check costs 0.1 s
   per evaluation at nside 64 and is adopted as the standard validity check of production
   index-field chains (D025).
2. **Templates.** With the per-region multiplier field, a component-separation-quality template
   (30 percent smooth error) costs +0.19 (sync) and +0.02 (dust) pixel-level RMSE against exact
   templates; the multiplier maps track the realized error where the data allow (correlation
   0.7 to 0.9 in the bright quartiles). For dust the cs template is usable at this depth (bright
   quartile as good as exact, fainter quartiles 1.3 to 1.7 times worse, no extra bias); for sync
   it is not (RMSE 2 to 3 times worse, field scale 0.35). The scalar law is a catastrophe
   (+0.48 and +0.24; sync field scale 0.87): the negative control fails as predicted, in a
   valid chain. F30.
3. **The exact-template chain is right where the data speak and prior-set where they do not.**
   Sync: the data-dominated regions (90 percent of the Fisher weight in the bright quartile)
   agree with the region-fit reference (z medians -0.6 to +0.2, 77 to 90 percent inside |z| <
   1); the Fisher-weighted level is right to 0.01; the single-index run is 1.3 sd from its
   pseudo-true value. The unweighted field mean (+0.13, 3 sd of its posterior sd) is the
   floating level of the 2300 regions with region Fisher sd 2.5 to 3.7, and that level carries
   the Occam bias of the marginalized per-region multipliers under the N(1, 0.5^2) hyperprior
   (+0.15 +- 0.05 over noise realizations by the exact marginal; item 9), partly held back by
   the field prior. F29 (as corrected after the stats review; the roughness mechanism first
   written here is withdrawn). Dust: the bright quartile sits +0.06 with
   the mu_d map at 1.029 because the fixed-T misfit is absorbed jointly by the multipliers and
   the index (1191 + 375 nats better than the reference), a path the global law could not take
   (E025 saw +0.5 percent). Corrected 2026-11-02 (E028 baseline): against the true-T region
   fit the bright quartile is unbiased (z median -0.06); the +0.06 is the difference between
   beta_eff (fixed-T fit) and beta_d (true), and the chain estimates the latter because the
   multipliers absorb the T misfit at the pivot. See F30's correction.
4. **Consequences.** (a) At realistic noise the unweighted mean of the nside-16 sync field is
   not a reportable quantity as the model stands; Paper A reports the Fisher-weighted level and
   the data-dominated regions, with the faint sky labelled prior-dominated. (b) The faint-sky
   bias is a property of the multiplier hyperprior width (sd 0.5), not of the field: E029 tests
   sd 0.1 in a chain (prediction: faint level within noise of the truth); E027 (nside 8 and 4)
   still places the regime boundary. (b2) Exact-target scans must integrate per-region nuisance
   multipliers analytically; profiling them is biased (item 10), and F28 is rechecked so. (c) The per-region
   multiplier field turns the fixed-T misfit into a dust-index bias of 0.06 in the bright sky;
   region-wise or sampled T_d comes before the multiplier field is used for dust on real data
   (M2.6), in the opposite order to D023's deprioritization. (d) Scans at a smooth reference
   field measure a smooth field's level; a prior-dominated field's level needs the chain (or a
   scan at the sampled roughness). (e) The per-pixel GLS CMB metric is uninformative at this
   depth; a posterior-CMB metric is needed (follow-up).
5. **Scope.** One sky (PySM d1s1c1), one noise realization (seed 2027), one template-error
   realization (seed 11), the E012 noise levels at nside 64; the template-S/N quartiles of a
   single realization correlate with the 23 GHz noise (faint sync quartile at -12 SE), which
   touches per-quartile statistics but not the posterior.

## Follow-ups
- E027: the resolution ladder (nside 8, 4) with the Gibbs-consistency check on each rung; does
  the sync level coupling vanish when sigma is data-determined.
- An informative prior on the field scale (or sigma fixed at a measured value) as the
  alternative fix for the faint-sky level; test by the roughness check before a chain.
- Region-wise or sampled T_d before per-region multipliers are used for dust (M2.6).
- The CMB metric at realistic noise: the per-pixel GLS is noise; use the posterior CMB from the
  amplitude step or a smoothed GLS.
- A second noise realization of the PySM noise-1 sky to separate realization effects from
  model properties in the per-quartile numbers.
