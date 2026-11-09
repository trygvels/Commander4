# E032: index fields with the collapsed hier prior at WMAP/LFI depth (noise scale 0.03), two runs of two chains

- **Date started / finished:** 2026-11-09 21:59 (scan), run A from 22:17:49, restarted 22:34:00 writing both in-flight chains / running
- **Branch and commit:** `nextgen` at 52061ad
- **Open question addressed:** Q16, Q18, Q3, and the multi-chain gap of `VALIDATION.md` (referee
  objections 1, 4 and 7 in `../reviews/2026-11-09_referee_report.md`)
- **Validation level:** L3 on a simulation with known truth (PySM d1s1 intensity), first two-chain
  read of the production index-field configuration

## Question

E031 showed that noise scale 1 of the PySM ladder is 24 to 90 times real WMAP/LFI intensity noise,
and that scale 0.03 is 0.7 to 2.7 times the real depth. Every result of E023 to E030 (F26 to F38)
was obtained at scale 1. Does the model of E030 chain 1, run unchanged on the same sky at scale
0.03, (a) recover the synchrotron and dust index fields without the field-mean displacement F38,
(b) give calibrated region uncertainties, (c) beat the trivial baselines (a constant index at the
truth mean, and the truth's own region means), and (d) do two chains with different seeds and
dispersed starts agree?

What changes with the answer: if the displacement is gone and the chains agree, the template-prior
biases of F29 to F38 are weak-data effects of a regime with no real intensity counterpart, and
Paper A can state the model as validated at real intensity depth on one sky. If a displacement
above 0.05 survives, the template prior is biased at any depth with a 30 percent template error and
must not be recommended for intensity.

## Target and priors

Data: d_p = M_p(beta_s,r(p), beta_d,r(p), T_r(p)) a_p + n_p per nside-64 pixel p, twelve bands
(23, 30, 33, 41, 44, 61, 70, 94, 100, 143, 217, 353 GHz), a_p = (CMB, synchrotron at 30 GHz, dust
at 545 GHz) in uK_RJ, n_p white with the known rms per band (8.4, 7.2, 7.2, 6.6, 7.2, 7.2, 8.4,
4.8, 2.88, 1.92, 2.88, 9.6 uK_RJ). Synchrotron SED a power law, dust a modified blackbody with the
oracle region-mean temperature map (nside 16, not sampled).

Priors, all as in E030 chain 1:

- CMB amplitude per pixel N(0, (1e4 uK)^2), flat for practical purposes (D021).
- Synchrotron and dust amplitudes a_p ~ N(mu_r(p) g_p, (s g_p)^2) with g the template, one
  multiplier mu_r per nside-16 region with hyperprior N(1, 0.3^2), one global scatter s per
  component with log s ~ N(log 0.3, 1). Amplitudes and multipliers are integrated out analytically
  (collapsed likelihood, `RegionIndexLikelihood`).
- Templates: the true 30 GHz synchrotron and 545 GHz dust amplitude maps times exp(f), f a smooth
  Gaussian field with rms 0.3 and FWHM 10 deg (seed 11), identical files to E030.
- Index fields: one value per nside-16 region (3072 per component), GMRF prior in the non-centred
  parameterization with hyperparameters (field mean mu, log sigma with sigma_scale 0.3, range
  parameter u with a log-uniform rho prior); field-mean priors N(-3.1, 0.3^2) and N(1.6, 0.3^2).

Sampled: the two non-centred index fields and their hyperparameters, and the two log scatters,
by NUTS with the dense finite-difference Hessian metric (refreshed at iterations 50 and 200), 400
warm-up iterations, 1200 iterations in total, burn 500. The target is the collapsed posterior of
these parameters, identical in form to E030 chain 1; only the noise level differs.

Chain A: seed 7, start at the prior means (-3.1, 1.6). Chain B: seed 8, start at -2.85 and 1.48,
on the other side of the truth means (-3.000, 1.540). The `beta` key sets only the start of the
field mean; the prior is in `spectral_index_prior` and is the same for both.

## Setup

Data dir `/Users/trygvels/data/e020_d1s1` (PySM d1s1c1 sky, seed 2027, the same sky and the same
noise realization as the scale 1 ladder, times 0.03; checked 2026-11-09 22:17: rms maps and
map-minus-sky rms agree with the list above). Parameter files
`param_gmrf16_cs_mumap_nc_Tmap_mu03_s7.yml` and `..._s8.yml` (chain 1's file with the data paths,
output dir, seed and start changed; diff checked before launch). 12 MPI ranks on the Mac mini M4,
chains one after the other (a second job beside a chain slows it threefold). Driver and logs under
`/Users/trygvels/data/e020_d1s1/e032/`.

Exact scan before the chains (`E032/scan_n003_sd03.txt`, 21:59 to 22:11): the noise-averaged
collapsed log posterior as one quartile of regions is shifted by delta from the noise-free
region-fit reference (the fit with the model's own T map, called mapT below), for the exact and
the cs templates, hyperprior sd 0.3, scatters profiled per quartile.

Quartile labels: the scan's (regions ranked by the Fisher sd of the scale 0.03 data-built
template). 1397 of 3072 synchrotron regions and 292 dust regions carry a different label from the
scale 1 labels of E027; tables are given with the scale 0.03 labels, and the comparison with E030
chain 1 repeats them with the E027 labels.

## Reference

Input truth (PySM index maps), and three derived references:

1. mapT region reference: the noise-free flat-prior fit of one index pair per nside-16 region with
   the model's T map (`E032/scan_n003_sd03_refs.npz`, keys `bs_mapT`, `bd_mapT`). It is the
   pseudo-truth of the model class. Against the truth's region means it has rms 0.0156 (sync) and
   0.0114 (dust); means -2.9972 and 1.5358.
2. Baselines for the pixel RMSE against the truth: constant at the truth mean 0.0629 (sync) and
   0.0374 (dust); the truth's own region means 0.0130 and 0.0187 (the floor of any nside-16 field).
3. The exact scan, as the prediction for the quartile-level displacement.

Scan result (read 22:16, before any chain output existed): optimum of the common shift per
quartile, relative to the mapT reference, with the curvature sd in brackets.

| template | s_Q1 | s_Q2 | s_Q3 | s_Q4 | s_all | d_Q1 | d_Q2 | d_Q3 | d_Q4 | d_all |
|---|---|---|---|---|---|---|---|---|---|---|
| exact | 0.000 (0.0009) | +0.004 (0.0040) | +0.007 (0.0054) | +0.011 (0.0063) | 0.000 | 0.000 (0.0001) | +0.002 (0.0017) | +0.008 (0.0042) | above +0.030 | 0.000 |
| cs | 0.000 (0.0009) | -0.007 (0.0040) | -0.009 (0.0054) | -0.016 (0.0063) | -0.001 | 0.000 (0.0001) | +0.005 (0.0017) | +0.012 (0.0042) | -0.018 (0.0093) | 0.000 |

Gains over delta 0 are 0.4 to 4.5 nats, except the exact-template faint dust quartile (8.0 nats at
the scan edge; a parabola through the scan ends puts its optimum at +0.039 with sd 0.0095). The
largest synchrotron displacement of the scale 1 counterpart was -0.30 (pass 7) and +0.3 (pass 3).
At this depth the scan predicts displacements of at most 0.02, one to three curvature sd of the
quartile mean and a fifth of a region's own posterior sd.

**Pre-registered (written 2026-11-09 22:20:01 by the file's mtime; chain A had started 22:17:49 and had
written no chain file by 22:27:18, so no chain output existed when this was fixed):**

- (i) D025 conditional-CDF check passes on both chains by its own rule.
- (ii) Agreement of the chains: split R-hat over four sequences (two chains, post-burn halves of
  350) below 1.05 for the two field means, the two log scatters and 95 percent of the region
  indices; below 1.10 for log sigma and the range parameter. The posterior means of each field mean
  differ between the chains by less than 0.5 posterior sd. The modern criterion R-hat below 1.01
  cannot be met by 700 draws at the effective sample sizes seen so far (30 to 500) and is reported,
  not required.
- (iii) Field level: the sky average of the posterior-mean field minus the sky average of the mapT
  reference within +-0.015 for both components. Against the truth mean: synchrotron within 0.02 of
  -3.000, dust within 0.015 of 1.540.
- (iv) No F38 displacement: per quartile, posterior mean minus mapT reference within +-0.03 for
  both components, and within +-0.01 in the bright quartile. Scan point predictions: synchrotron
  0.000, -0.007, -0.009, -0.016; dust 0.000, +0.005, +0.012, -0.018.
- (v) GMRF hyperparameters: sigma_s in 0.04 to 0.09 and sigma_d in 0.02 to 0.05 (the truth's region
  means have sd 0.0615 and 0.0323).
- (vi) Pixel RMSE against the PySM truth: synchrotron 0.045 to 0.062, dust 0.025 to 0.033, both
  below the constant-at-truth-mean baseline (0.0629, 0.0374). Arithmetic behind the ranges: region
  Fisher sd 0.03 to 0.12 (sync) and 0.002 to 0.064 (dust) shrunk by a field prior of sd 0.06 and
  0.03, in quadrature with the reference offset (0.0156, 0.0114) and the floor (0.0130, 0.0187).
  The gain over the constant is expected to be small for synchrotron (at most 30 percent) and
  larger for dust.
- (vii) Calibration against the mapT reference: |z| below 1 for 55 to 80 percent of the regions in
  every quartile, median z within +-0.3.
- (viii) Multiplier maps against the realized region mean of a/g: correlation above 0.9 in every
  quartile and rms difference below 0.10 over the sky.
- (ix) Amplitude scatters: log s between -3.1 and -2.6 for both components (the within-region
  residual of the template error, 0.06 at scale 1).
- (x) CMB, by the E033 split: bias term of the flat-prior GLS CMB estimate at the posterior-mean
  fields below 7 uK over the sky (floor 5.95 at the truth's region means; constant at the truth
  mean 16.3) and below 0.8 uK at |b| above 20 (floor 0.42, constant 0.95); total residual within 3
  percent of the quadrature sum of the bias and the noise term 2.79 uK.
- (xi) No runaway: no region with |beta_s minus reference| above 0.5 (F21 saw a runaway at this
  depth with a data-built template and one global scatter).
- Cost: 31 to 125 gradients per iteration, 7 to 20 hours per chain.

Alternative outcomes, stated in advance: a field-level displacement above 0.05 in either component
means the template prior is biased at real depth and objection 2 of the referee report stands in
full; a failure of (ii) means the production sampler is not validated at this configuration
(objection 4); a synchrotron RMSE above 0.0629 means the field does not beat a constant index for
this sky at this depth, which would be reported as the answer to objection 7.

Limits known in advance: one sky, one noise realization, one template-error realization, oracle T
map, truth-informed quartile labels only for reporting. Field-mean priors centred 0.10 and 0.06
from the truth means with sd 0.3. None of the hyperpriors was tuned on this data set.

**Amendment, fixed 2026-11-09 22:34:36 by the clock, with no chain file in existence (the first run was stopped
at 22:33:44 with an empty `chains_compsep/`):** the pipeline audit
(`../reviews/2026-11-09_pipeline_assumptions_audit.md`, item G9) says that Commander4 computes two
chains in every run and writes one. Checked by reading the code: `gibbs_schedule` in `cli.py` walks
`chain in (1, 2)` inside every iteration; the random stream is keyed on (seed, chain, iteration,
rank); the NUTS state is kept per (chain, group) in `_STATES` of `compsep/index_field_sampler.py`;
`file_io/chain_writer.py:195` writes only the chains listed in `output.chains.write`, which is
`[1,]` in every parameter file used so far. So every chain of E020 to E030 had a twin with the same
start and an independent random stream that was computed and thrown away. Run A was restarted at
22:34:00 with `output.chains.write: [1, 2]`, and run B has the same setting. The computation is
unchanged, so `chain01` of each run is the chain pre-registered above and (i) to (xi) are read on
the two `chain01` sets as written. The aborted directory is kept in `$D/e032/aborted_s7_2217/`.
Added criteria for the four chains (A1, A2 share start A; B1, B2 share start B):

- (i-b) D025 passes on each of the four chains.
- (ii-b) Split R-hat over eight sequences below 1.05 for the two field means, the two log scatters
  and 95 percent of the region indices; below 1.10 for log sigma and the range parameter.
- (ii-c) The posterior means of each field mean agree between every pair of the four chains to
  within 0.5 posterior sd. With an effective sample size near 100 per chain the expected largest
  pairwise difference is about 0.3 sd.
- (ii-d) The between-start difference (mean of A1, A2 minus mean of B1, B2) of each field mean is
  below 0.35 posterior sd in absolute value. A larger value with (ii-c) satisfied inside each start
  means the chains still remember their start after 500 iterations.

Wall time per written iteration covers both chains. Gradient counts per iteration in the chain
files are per chain.

## Command

```bash
cd /Users/trygvels/code/Commander4 && source .venv/bin/activate
D=/Users/trygvels/data/e020_d1s1
python nextgen/scripts/e023_exact_scan.py --dir $D --out nextgen/experiments/E032/scan_n003_sd03.txt \
  --cells exact_mapT_perregmarg cs_mapT_perregmarg --mu-sd 0.3 --cmb-sd 10000 \
  --fine -0.03 0.03 0.001 --scans s_Q1 s_Q2 s_Q3 s_Q4 s_all d_Q1 d_Q2 d_Q3 d_Q4 d_all
nohup $D/e032/run_e032_chains.sh >> $D/e032/run_e032_chains.log 2>&1 &  # s7 then s8, 12 ranks each, chains 1 and 2 written
python nextgen/scripts/e032_analyze.py --dir $D --tags gmrf16_cs_mumap_nc_Tmap_mu03_s7 \
  gmrf16_cs_mumap_nc_Tmap_mu03_s8 --refs nextgen/experiments/E032/scan_n003_sd03_refs.npz \
  --t-map $D/truth_Td_nside16.fits > nextgen/experiments/E032/chains_n003.txt
```

## Result

Pending (chain A started 2026-11-09 22:17:49).

## Negative control

The constant field at the truth mean and the truth's region means, through the same RMSE and CMB
bias metrics. A scalar-index chain with the same amplitude prior is the non-oracle baseline and is
queued behind the pair if the pair's timing allows.

## Conclusion

Pending.

## Follow-ups

- The same pair with exact templates (the referee's test for F24 on the PySM maps).
- The scalar-index baseline chain.
