# E023: Robustness of the hierarchical index target to template and SED misspecification on the PySM sky, by exact-target scans

- **Date started / finished:** 2026-10-29 /
- **Branch and commit:** nextgen, fe50735
- **Open question addressed:** Q16 (robust amplitude-template prior and dust SED at high S/N); feeds M2.5 (PySM skies) and the real-data model (M3)
- **Validation level:** L3 (misspecified simulation with known truth), stage 1 exact target only

## Question
On the PySM d1s1c1 sky at realistic noise (E012 levels, noise scale 1) and at the E020 level
(0.03), how far does the data term of Commander4's collapsed hier index target (the likelihood
with the amplitudes integrated out and the amplitude hyperparameters profiled, no index prior)
move a whole S/N quartile of nside-16 regions away from the noise-free region fit, which of the
two misspecifications (data-built template under one global scatter; fixed T_d) drives it, and
does a per-quartile scatter (or per-quartile scatter and mean) remove it? What it can say: whether
this misspecification type at this depth biases the data term; what it cannot say: the posterior
bias under the GMRF prior (stage 2), the transfer to real data (quartiles here come from
templates, but noise is white, bandpasses are deltas, beams absent), fixed against sampled T
(only fixed against known T), correlated template errors beyond what the scatter candidates
model, and seed dependence (one realization per cell; the expected-log-likelihood scan removes
the noise scatter of the peak, not the template realization). Decides: whether the amplitude
prior needs a change before any stage-2 chain on a misspecified sky, and which candidate goes
first.

## Target and priors
The stage-1 object is the collapsed per-pixel marginal that Commander4's `hier` likelihood
evaluates, as a function of the index fields alone, with the amplitudes integrated out and the
amplitude hyperparameters at their joint MAP (profile) for each index state:

  L(beta) = max over (mu_s, s_s, mu_d, s_d) of
            sum_p [ -1/2 r_p^T C_p^-1 r_p - 1/2 log det C_p ] + log p(mu_s) + log p(log s_s)
            + log p(mu_d) + log p(log s_d),

with, per pixel p (every one of the 49,152 nside-64 pixels, no mask) and the 12 bands, M_p =
[m_cmb, m_sync(beta_s,p), m_dust(beta_d,p)] (uK_RJ per
unit amplitude; CMB in uK_CMB, sync a power law referenced at 30 GHz, dust one MBB referenced at
545 GHz at the fixed T = 20.93 K, or at the PySM T_d,p where stated), prior means m_p = (0,
mu_s g_s,p, mu_d g_d,p), prior variances T_p = diag(100^2, (s_s g_s,p)^2, (s_d g_d,p)^2), C_p =
M_p T_p M_p^T + N_p (white, the rms maps), r_p = d_p - M_p m_p, and the E021 hyperpriors mu_c ~
N(1, 0.5^2), log s_c ~ N(log 0.3, 1^2). This is the E018/E021 target (stats-reviewed there)
minus the index prior and with the hyperparameters profiled instead of marginalized; the
profile is corrected by the Laplace term -1/2 log det H(theta*) (H the Hessian of the profiled
objective, central differences of its analytic gradient) at five deltas per scan, so that the
profile-against-marginal difference is measured where it could move a peak (at noise 1 the
block curvature is weak). The index prior is left out on purpose: the question is whether the
data term is biased; an index prior can shift as well as shrink (F24, iid), which is stage 2.

Two modes. The realization scan evaluates L on the simulated data d; its peak scatters by the
block curvature sd, which at noise 1 is 0.08 to 0.55 for the fainter dust quartiles, larger
than any threshold below, so it is reported only as a secondary line with that sd. The primary
scan is the expected log-likelihood on the noiseless sky s_p from truth.npz with the noise
covariance kept,

  E[l_p] = -1/2 r_p^T C_p^-1 r_p - 1/2 tr(C_p^-1 N_p) - 1/2 log det C_p,  r_p = s_p - M_p m_p,

whose peak is the pseudo-true bias of the target with no noise scatter (the trace term depends
on beta and theta and is kept). The true-T cells use the PySM T_d,p in the dust column and are a
reference outside Commander4's model class (ThermalDust.T is a scalar); the real-data analogue,
sampled T, is not tested here.

Candidate priors (stage 1): (a) the global scatter s_c replaced by one scatter per S/N quartile
of the component, s_c,q, q = 1..4 (ten hyperparameters: two means, eight log scatters, each
with the E021 hyperprior); (b) mean and scatter per quartile, mu_c,q and s_c,q (16
hyperparameters), because F23's mechanism is a trend of E[a | g] / g with g, a mean effect that
a scatter alone cannot absorb. The quartiles come from the Fisher sd of each nside-16 region
computed with the data-built template as the amplitude and the prior-mean indices (computable
on real data), the same quartiles in every cell. The scans alone cannot tell a prior that fixes
a bias from one that merely loosens until the data win (and a profiled s_q that inflates moves
the marginal toward the flat-prior one with its volume term, F3/F4); the record therefore
reports the profiled s_c,q and mu_c,q against the realized sd and mean of a/g per quartile (on
the E021 sky they agree to 0.001), and the stage-2 chains with their z calibration are where
the distinction is made.

Statistic: for each component c and each quartile Q, the block shift delta that maximizes
L(beta_ref + delta 1_Q) over delta in [-0.1, 0.1] in steps of 0.01 with the other component at
its reference, plus a wide sync scan over [-3, +3] in steps of 0.2 in the faintest quartile (the
E020 runaway modes at beta_s below -9 and near -0.6 are at delta -6 and +2.4; the second is in
range, the first is read from the slope at -3); the curvature sd of the block at the peak; the
split into the chi-square, log-determinant, trace and Laplace parts. The block shift is a
projection (the mean bias over 768 regions, blind to zero-mean per-region bias) and a
conditional (the other component held at its reference). The reference beta_ref is the
noise-free region fit defined under Reference, per SED level.

Factor cells (2 x 2 x 2):
- noise scale 1 (new sky, same seed 2027, `--sky pysm --noise-scale 1`) and 0.03 (the E020 sky);
- templates: data-built (the E020 pseudo-templates, g/true 0.41 to 1.66 for sync at |b| > 20 at
  noise 0.03; 0.20 to 1.79 at noise 1) and exact (g = the true amplitude, floored positive, with
  the scatter profiled, which runs it to the guard): the no-template-error reference. The
  model-story law of E021 (a = g r drawn given g) cannot be realized on a fixed sky, and the
  only template law with scatter that keeps the PySM amplitudes is g = a / r, which F23 shows is
  itself a misspecification; so the template factor has the levels "as built from the data" and
  "none".
- dust SED: fixed T = 20.93 K and the PySM T_d,p (the no-SED-misfit reference).

The exact-template, true-T cell is not from the model class either: the PySM indices vary
inside the nside-16 regions, so that cell measures the region model's own misfit (the one
misspecification E021 did not have) and is the baseline the other cells are read against.
At noise scale 1 the profiled Fisher sd of an nside-16 region is 2.7 (sync) and 0.45 (dust) in
the median, so every region is prior-dominated there and the scans measure what the amplitude
prior does with the index; a quartile block of 768 regions still has a Fisher sd of order 0.1.

## Setup
Skies at nside 64, 12 bands 23 to 353 GHz, white noise, from `nextgen/scripts/pysm_intensity_sim.py`:
`/Users/trygvels/data/e020_d1s1` (noise 0.03, E020) and `/Users/trygvels/data/e023_pysm_n1` (noise
scale 1, this experiment). Model-story templates for the PySM sky: built in the scan script from
`truth.npz` (a_sync, a_dust) with the same r law and seed as E021 attempt 2. Scan script:
`nextgen/scripts/e023_exact_scan.py` (v2: numpy 12 x 12 marginal in realization and expected
modes, per-pixel T_d, global / per-quartile-scatter / per-quartile-scatter-and-mean
hyperparameter layouts, L-BFGS-B profile with analytic gradients warm-started along each scan,
log s bounded below by HIER_LOGS_MIN = -8, Laplace term by finite-difference Hessians at five
deltas per scan, noise-free region-fit reference by Nelder-Mead per region). About 3 s per
state plus 5 to 25 s per Laplace point (4 to 16 hyperparameters); seven cells of 199 states in
expected mode and two in realization mode, about 2.5 h per sky on one core of the Mac mini
alongside the E021 seed-2028 queue. Cells: exact_trueT_global (baseline), exact_fixT_global,
data_trueT_global, data_fixT_global (the Commander4 configuration), data_fixT_perq,
data_trueT_perq, data_fixT_perqmu. Smoke test of v1 (`--quick`, 0.03 sky, data_fixT_global,
realization mode): the profile at the region-mean truth gives mu_s 1.144, s_s 0.673, mu_d 0.901,
s_d 0.270 (the E020 chain's s_s 0.36 was at its runaway state, not at the truth).

## Reference
The reference index of a region is the noise-free region fit: for each nside-16 region the
(beta_s, beta_d) minimizing the profile chi-square of the noiseless sky (truth.npz `sky`) over
three free amplitudes per pixel with the rms weights, under the cell's SED (fixed T or the PySM
T_d,p). It is what a region-constant index model with flat amplitude priors would converge to
with infinite data, independent of the templates and of the noise scale; the region means of
the per-pixel truth (beta_s; beta_eff3 for fixed T) are not it on a sky whose indices vary
inside the regions (rms 0.013 sync, 0.020 beta_eff3 within nside-16 regions on the PySM sky), and
are printed next to it. On the E021 sky (region-constant truth) the two coincide to 1e-5.

Code control (done 2026-10-29, `nextgen/scripts/e023_kernel_check.py`, output in
`E023/kernel_check_e021.txt`): the numpy marginal agrees with Commander4's own evaluation
(`marginal_loglike_regions` through `RegionIndexLikelihood.from_arrays` plus the hier trailing
terms) on the E021 sky at four states (the reference field and two hyperparameter settings, a
+0.5 sync shift of the faintest quartile, a +0.05 dust shift of the brightest) to 4e-12 relative
in the state differences and 5e-7 nats absolute on a total of 1.4e9, after setting h/k to
Commander4's astropy value (with the simulator's 7-digit constant the dust-shift difference was
off by 2e-7 relative). Model control: on the E021 sky (model-story templates, region-constant
truth) the expected scan must peak at zero in every quartile with the per-quartile hyperparameters
at the realized a/g statistics (smoke test: argmax 0.000 in all eight scans, s_q and mu_q within
0.001 of the realized values, Laplace variation below 0.01 nats).

Pre-registered reading (stage 1, on the expected scan):
- A cell is "unbiased at this depth" when every quartile's preferred shift is below 0.02 (sync)
  and 0.01 (dust), smaller than the E021 region-level noise error (0.033 / 0.018) at the same
  nside, and the wide sync scan's maximum is within 0.1 of the reference.
- Attribution: compare the data-built / fixed-T cell with the two single-misspecification cells
  (data-built / true-T and exact / fixed-T) and with the exact / true-T baseline; the
  misspecification whose removal brings the cell to "unbiased" is the driver. The exact / true-T
  cell measures the region model's own misfit (indices varying inside the regions) and the
  within-region S/N weighting; it is the baseline the others are read against, not a zero.
- A candidate (per-quartile scatter; per-quartile scatter and mean) is carried to stage 2 only if
  it brings the data-built cell to "unbiased" at noise scale 1 without degrading the exact-template
  cell, and its profiled s_c,q (and mu_c,q) sit within a factor 1.5 (0.1) of the realized sd
  (mean) of a/g in every quartile; otherwise it has only loosened the prior.
- A shift in the exact-template / true-T cell beyond 0.01 in the two brightest quartiles is a
  defect in the scan or in the reference fit (the E021 control bounds the scan code), not a prior
  result; stop and find it.
- The Laplace term's variation over a scan is reported; where it exceeds a tenth of the gain at
  the argmax the profile reading is replaced by the Laplace-corrected argmax.

F21 reproduction (a result, not a control): at noise 0.03 with data-built templates and fixed T
the wide sync scan and the dust Q1 scan should show the E020 preferences (sync away from the
reference in the faint quartile; dust up in the brightest). A one-block scan with s profiled
may not reproduce a joint runaway under the GMRF prior; if it does not, that is recorded as such.

Stage 2 (chains, after the E021 seed-2028 queue): gmrf16 hier with the prior stage 1 selects on
the noise-1 PySM sky, E021 metrics and the Q16 acceptance (index RMSE within two bootstrap SE of
the E021 consistent-sky result per quartile; no region with |beta_s - truth| above 0.5), with the
z calibration against the region-fit reference. Written up as a second stage of this record or
as E024 if the prior changes the code.

## Command
```bash
python nextgen/scripts/e023_kernel_check.py --dir /Users/trygvels/data/e021_consistent_v2 > nextgen/experiments/E023/kernel_check_e021.txt
python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e023_pysm_n1 --nside 64 --noise-scale 1 --seed 2027 --sky pysm --n-iter 1200 --n-warmup 400 --refresh-mass-at 50 200 --no-adapt-mass
# scratchpad run_e023_scans.sh, one sky after the other (started 2026-10-29 11:46):
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e021_consistent_v2 --out nextgen/experiments/E023/scan_e021_control.txt --cells exact_fixT_global data_fixT_global --real-cells data_fixT_global
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e020_d1s1 --out nextgen/experiments/E023/scan_n003.txt
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e023_pysm_n1 --out nextgen/experiments/E023/scan_n1.txt
```

## Result
(pending)

## Negative control
Replaced by the code control (kernel agreement) and the model control (E021 sky) under
Reference, after the stats review; the F21 reproduction at noise 0.03 is reported as a result.

## Design review
stats-reviewer 2026-10-29 (target re-derived from `index_field.py` and `compsep_pymod.cc`,
pseudo-band form shown equal to the record's formula to 2e-15): formula CONFIRMED; design as
first written NOT able to decide its reading at noise 1 (single-realization peaks scatter by the
block curvature sd, 0.08 to 0.55 for faint dust quartiles), reference confounded (region means
of a within-region-varying truth are not the region model's optimum), profiling not shown valid
at noise 1 (Laplace term), candidate counted wrongly (ten hyperparameters), quartiles from the
truth not computable on real data, the per-quartile scatter unable to separate "fixes" from
"loosens", the F21 reproduction not a valid code control, the Question's posterior claim and the
Decides clause too wide. All taken: expected-log-likelihood scan as the primary line with the
trace term, Laplace term at five deltas per scan, noise-free region-fit reference, template-based
quartiles, mu-per-quartile candidate with the realized a/g comparison, kernel-agreement check
(PASS at 4e-12), wide sync range [-3, 3], Question and Decides narrowed. The reviewer's 1e-8
relative agreement is met only with Commander4's h/k constant; with the simulator's 7-digit
constant the dust-shift difference was 2e-7 off, the known constant mismatch of the E021 design
review, not a structural one.

## Conclusion
(pending)

## Follow-ups
