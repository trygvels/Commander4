# E032: index fields with the collapsed hier prior at WMAP/LFI depth (noise scale 0.03), two runs of two chains

- **Date started / finished:** 2026-11-09 21:59 (scan), run A from 22:17:49, restarted 22:34:00 writing both in-flight chains, ended 2026-11-10 05:17:06; run B 05:17:06 to 14:12:33 / 2026-11-10 14:17
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

**Run A (seed 7, start at the prior means), interim read written 2026-11-10 05:21; run B started 05:17:06 and
is untouched by this read.** Run A ran 2026-11-09 22:34:00 to 2026-11-10 05:17:06 (6 h 43 min for
both in-flight chains), rc 0, 1200 iterations per chain, 700 post-burn. Outputs:
`E032/runA_two_chains.txt`, `runA_mu_maps.txt`, `runA_cmb_split.txt`, and the two conditional
checks. A1 is `chain01`, the chain the pre-registration names; A2 is its twin with the same start
and an independent random stream. Numbers are for A1, with A2 where it differs.

| Item | Pre-registered | Run A | Verdict |
|---|---|---|---|
| (i) D025 on each chain | passes | A1: KS p 0.66 (sync), 0.81 (dust); A2: 0.09, 0.32; 100 draws each, scan step 0.6 to 0.75 conditional sd | pass, A2 sync marginally |
| (ii) field means, two chains | split R-hat below 1.05, means within 0.5 sd | sync 1.014 (0.23 sd), dust 1.047 (0.20 sd) | pass for A1 and A2; run B pending |
| (ii) amplitude scatters | below 1.05 | 1.000, 0.999 | pass |
| (ii) region indices | 95 percent below 1.05 | sync 100 percent (max 1.018), dust 99.8 percent (max 1.059) | pass |
| (ii) GMRF log sigma and range | below 1.10 | sync 1.006 and 1.009; dust 1.40 and 1.44, chain means 1.4 and 1.5 sd apart, ESS 6 to 13 | fail for dust |
| (iii) field level against the reference | within 0.015 | sync -0.019, dust +0.009 | fail for sync by 0.004 |
| (iii) field level against the truth mean | sync within 0.02, dust within 0.015 | -0.016, +0.005 | pass |
| (iv) quartiles against the reference | within 0.03, bright within 0.01 | sync -0.009, -0.026, -0.024, -0.016; dust +0.003, +0.014, +0.026, -0.006 | pass |
| (v) GMRF sigma | sync 0.04 to 0.09, dust 0.02 to 0.05 | 0.094 +- 0.007; 0.074 +- 0.006 (A2 0.067) | fail, both too large |
| (vi) pixel RMSE | sync 0.045 to 0.062, dust 0.025 to 0.033, below the constant's 0.0629 and 0.0374 | sync 0.0653 (A2 0.0647); dust 0.0618 (A2 0.0590) | fail: neither field beats a constant index |
| (vii) calibration | 55 to 80 percent inside one sd per quartile, median z within 0.3 | sync 69, 66, 67, 70 percent, median z -0.20, -0.36, -0.31, -0.23; dust 68, 60, 49, 36 percent, median z +0.22, +0.49, +0.72, -0.11 | sync: coverage passes, median z misses in Q2 and Q3; dust fails in the faint half |
| (viii) multiplier maps | correlation above 0.9 per quartile, rms below 0.10 | sync 0.997, 0.991, 0.987, 0.988, rms 0.046; dust 0.998, 0.996, 0.983, 0.908, rms 0.087 | pass |
| (ix) amplitude scatters | log s in -3.1 to -2.6 | -2.710, -2.600 | pass, dust at the edge |
| (x) CMB bias | below 7 uK all sky, below 0.8 at abs(b) > 20, closure within 3 percent | 4.55 and 0.78; total 5.35 against 5.34 | pass |
| (xi) no runaway | no region beyond 0.5 | largest 0.33 (sync), 0.26 (dust) | pass |

Cost: 31.0 gradients per iteration after burn-in in both chains (88 and 75 during the 400
warm-up iterations), 38 and 23 divergences, all in warm-up. The conditional sd of a field mean
given the rest of the state is 0.0008 (sync) and 0.0001 (dust), so the scans used half-widths of
0.004 and 0.0006.

Hyperparameters, A1 and A2: sync field mean -3.0160 +- 0.0030 (ESS 50) and -3.0167 +- 0.0031
(101); dust field mean 1.5449 +- 0.0028 (15) and 1.5455 +- 0.0026 (60); dust log sigma -2.601 (13)
and -2.711 (8). Region indices: ESS per chain median 320 (sync) and 155 (dust).

CMB bias split at this noise (noise term 2.79 uK), all sky and abs(b) > 20: chain fields 4.55 and
0.78; constant at the truth means 16.3 and 0.95; constant at the chain's own field means 15.5 and
0.96; truth region means 5.95 and 0.42; region-fit reference 4.43 and 0.45. With the true dust
column the chain's sync field leaves 0.59 and 0.24; with the true sync the chain's dust field
leaves 4.37 and 0.66. So the fields beat a constant index on the CMB bias by 18 percent at high
latitude and by 72 percent over the sky, and what remains is the dust field's.

**Reading of run A.**

1. *The field-mean displacement of F38 is gone at this depth.* The sync field mean is 0.016 from
   the truth where the scale 1 chains were +0.30 and -0.16 away, every quartile is within 0.03 of
   the region-fit reference, and the multiplier maps track the realized template error with
   correlations of 0.91 to 1.00 (0.25 to 0.72 at scale 1). This is the question the experiment
   was run for, and the answer is that F29 to F38 are weak-data effects.
2. *The fields do not beat a constant index.* Pixel RMSE 0.065 against 0.063 for synchrotron and
   0.060 against 0.037 for dust. The stated alternative outcome applies: this is the answer to
   the referee's objection 7 on this sky at this depth. The inferred field amplitudes are too
   large (0.094 against a truth scatter of 0.062; 0.07 against 0.032), so the fields carry noise
   that a tighter field would have shrunk away.
3. *Dust is not calibrated in the faint half* (36 to 49 percent of regions inside one sd), and
   its field amplitude and range do not agree between the two chains (R-hat 1.40 and 1.44). A
   single chain would have shown neither. Every earlier statement about a dust field amplitude
   from one chain (F32 to F38) carries this caveat.
4. A candidate for the dust scatter, not tested: the dust amplitude and template are referenced
   at 545 GHz, above every band, and the faint-quartile multipliers are uncertain by 0.11 to 0.16,
   which a lever arm of 1 / ln(545 / 353) = 2.3 turns into index scatter of that size divided by
   the lever. It is the reference-against-pivot effect of F44 again.
5. Noise scale 0.03 is WMAP and LFI white-noise depth; for the HFI bands it is 2.6 to 33 times
   the real noise (E031), so the dust numbers are not a statement about real dust data.

**Run B (seed 8, start at -2.85 and 1.48) and the four chains, read 2026-11-10 14:17.** Run B ran
2026-11-10 05:17:06 to 14:12:33 (8 h 55 min), rc 0, 1200 iterations per chain. Outputs:
`E032/chains_n003.txt`, `prereg_two_chain_rhat.txt`, `runB_mu_maps.txt`, `cmb_split_four_chains.txt`
and the conditional checks. B1 is `chain01` of run B, B2 its twin.

| Item | Pre-registered | Result | Verdict |
|---|---|---|---|
| (i), (i-b) D025 on every chain | passes | KS p, sync and dust: A1 0.66, 0.81; A2 0.09, 0.32; B1 0.51, 0.17; B2 0.19, 0.62 | pass |
| (ii) field means, A1 against B1 | split R-hat below 1.05, means within 0.5 sd | sync 1.109 and 0.78 sd; dust 1.343 and 1.41 sd | fail |
| (ii) GMRF log sigma and range, A1 against B1 | below 1.10 | sync 1.015, 1.016; dust 1.111, 1.079 | sync pass; dust log sigma fails by 0.011 |
| (ii) amplitude scatters and region indices | below 1.05; 95 percent of regions | 1.000 and 0.999; 100 and 99.2 percent | pass |
| (ii-b) eight half-chains | field means and scatters below 1.05; log sigma and range below 1.10 | sync mean 1.050, dust mean 1.179; dust log sigma 1.418, range 1.410; sync 1.015 and 1.016; regions 100 percent | fail for both field means and for the dust field |
| (ii-c) every pair of chains | field means within 0.5 sd | largest 0.79 sd (sync), 1.41 sd (dust) | fail |
| (ii-d) between starts | below 0.35 sd | +0.42 (sync), -0.83 (dust) | fail |
| (iii) field level, B1 | reference within 0.015; truth within 0.02 and 0.015 | sync -0.021 and -0.018; dust +0.013 and +0.009 | sync against the reference fails, the rest passes |
| (iv) quartiles, B1 | within 0.03, bright within 0.01 | sync -0.009, -0.028, -0.028, -0.020; dust +0.003, +0.014, +0.030, +0.006 | pass |
| (v) GMRF sigma | 0.04 to 0.09; 0.02 to 0.05 | sync 0.093 to 0.096 in the four chains; dust 0.067 to 0.079 | fail |
| (vi) pixel RMSE | below the constant's 0.0629 and 0.0374 | sync 0.0647 to 0.0668; dust 0.0590 to 0.0624 | fail in all four chains |
| (vii) calibration, B1 | 55 to 80 percent, median z within 0.3 | sync 69, 66, 64, 69 percent, median z -0.19, -0.41, -0.34, -0.28; dust 70, 60, 48, 34 percent | as run A |
| (viii) multiplier maps | correlation above 0.9, rms below 0.10 | sync 0.987 to 0.997, rms 0.047; dust 0.907 to 0.998, rms 0.086 | pass |
| (ix) amplitude scatters | -3.1 to -2.6 | -2.710, -2.600 in all four | pass |
| (x) CMB bias | below 7 and 0.8 uK | 4.55 and 0.76 to 0.79 in all four | pass |
| (xi) no runaway | none beyond 0.5 | largest 0.33 | pass |

Chain means of the hyperparameters, A1, A2, B1, B2, with ESS: sync field mean -3.0160 (50),
-3.0167 (101), -3.0184 (55), -3.0169 (129), posterior sd 0.0031; dust field mean 1.5449 (15), 1.5455
(60), 1.5489 (10), 1.5462 (44), sd 0.0028; dust log sigma -2.601, -2.711, -2.544, -2.596 (ESS 8 to
15), sd 0.085. Cost: 31.0 gradients per iteration after burn-in in A1, A2 and B1, and 101 in B2,
which also had 13 divergences after burn-in (none in the others). The two chains of a run adapt
separately, and B2's adaptation ended at a smaller step.

**Reading of the four chains.**

1. *The disagreement is slow mixing, not the start.* Chains with the same start differ as much as
   chains with different starts (A1 against A2: 1.4 sd in the dust log sigma). The chain means of
   the field means scatter by about one posterior sd, where the quoted ESS of 10 to 130 would
   allow 0.1 to 0.3. So the ESS of the hyperparameters is overestimated within a chain, and the
   700 draws hold a few independent samples of them at most.
2. *What disagrees is small in the units that matter.* The four sync field means span 0.0024 and
   the four dust means 0.0040, against offsets from the truth of 0.017 and 0.007. The region
   indices, the multiplier maps, the amplitude scatters and the CMB bias agree across the four
   chains to their last quoted digit. The field means are known to about their posterior sd, not
   better, and the dust field amplitude and range are known to about two posterior sd.
3. *D025 passes on all four chains while R-hat fails.* The conditional check is necessary and
   not sufficient, as D028 states. The four-chain criterion is what found this.

## Conclusion

Pre-registered outcome, in the terms fixed in advance: no field-level displacement above 0.05,
so the template prior is not biased at this depth in the sense of F38; (ii) fails, so the
production sampler is not validated at this configuration for the field hyperparameters; the
synchrotron RMSE is above 0.0629, so the field does not beat a constant index on this sky at this
depth.

1. **The field-mean displacement of F29 to F38 is a weak-data effect.** At WMAP and LFI
   white-noise depth, with the same 30 percent template error, the sync field mean is 0.017 from
   the truth (0.30 and 0.16 at scale 1), every quartile is within 0.03 of its reference, the
   multipliers follow the template error, and the CMB bias is 18 percent below a constant
   index's at high latitude and 72 percent below over the sky.
2. **The index fields are not better than a constant index.** Pixel RMSE 0.065 to 0.067 against
   0.063 (sync) and 0.059 to 0.062 against 0.037 (dust), in four chains. The inferred field
   amplitudes are 1.5 and 2.2 times the truth's scatter. Dust is not calibrated in the faint
   half of the sky.
3. **The sampler is not validated for the field hyperparameters at this configuration.** Four
   chains from two starts fail the split R-hat criterion for both field means and for the dust
   field's amplitude and range, after 700 draws each, at a cost of 7 to 9 hours per run. Region
   indices and amplitude scatters pass. Every single-chain statement about a field mean, a field
   amplitude or a range in F29 to F38 has an unquantified Monte Carlo error of about one
   posterior sd.
4. With this the intensity template-prior sequence is closed (D029, proposed): the model has no
   displacement at real depth and no advantage over a constant index there.

What this does not establish: anything about real dust data (for the HFI bands this noise is 2.6
to 33 times the real one); whether exact templates or a tighter field prior would beat the
constant; the cause of the overestimated field amplitudes (the reference-against-pivot lever of
the dust amplitude at 545 GHz is a candidate, untested). One sky, one noise realization, one
template-error realization, oracle T map.

Changes: F45 added; F38 and F28 get the real-depth result; the referee's objections 4 and 7 are
answered for this configuration (not validated; does not beat a constant).


## Negative control

The constant index at the truth mean and the truth's own region means are the baselines of items
(vi) and (x). No deliberately wrong sampler was run.

## Follow-ups

- Not run under D029 (proposed): the same pair with exact templates, and the scalar-index
  baseline chain. They would say whether the template error or the field prior is what keeps the
  field from beating a constant.
- For Paper B: longer chains or a better update of the field hyperparameters before any field
  amplitude or range is quoted; the four-chain criterion with an ESS floor.
