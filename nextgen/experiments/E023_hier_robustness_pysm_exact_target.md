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
**Model control, E021 consistent sky (`E023/scan_e021_control.txt`, 11:46 to 13:00):** the
region fits equal the region means (rms 0.0000; region-constant truth). exact_fixT_global
(templates equal to the true amplitudes; the profiled scatters run to 0.0012 and 0.0003 at the
guard) and data_fixT_global (the model-story templates of E021; profiled mu_s 1.0027, s_s 0.2990,
mu_d 0.9999, s_d 0.2987 against the realized 1.0031 / 0.2989 / 0.9999 / 0.2990) peak at delta
0.000 in all eight fine scans and the wide sync scan, in expected mode; the realization mode of
data_fixT_global peaks at 0.000 in eight scans and at -0.010 in sync Q4 (gain 1.4 nats,
curvature sd 0.0070: a 1.4-sd noise excursion, as expected of a single realization). The
Laplace term varies by at most 0.02 nats over a scan with the model-story templates and by
up to 2.9 nats (dust Q1) with the exact templates at the guard, where it moves no argmax.
Control PASS: the expected-log-likelihood scan, the profile, the reference fit and the
quartile machinery give zero bias where the model is exact.

**PySM sky at noise 0.03 (`E023/scan_n003.txt`, from 13:00; partial, five of seven cells at
16:20):** the region fits differ from the region means (sync rms 0.056 at fixed T, 0.021 at
true T, max 0.45 to 0.49; dust rms 0.012 and 0.009, max 0.12): the region means are not the
region model's optimum on this sky, as the reviewer said. Template-S/N quartiles; the
data-built template's realized a/g per sync quartile (bright to faint) has mean / sd 0.89 /
0.25, 0.83 / 0.31, 1.04 / 0.42 and 77 / 155: in the faintest quartile the data-built template
(23 GHz map scaled and smoothed, CMB included, floored at its first percentile) is not a
template of the synchrotron at all.
- exact_trueT_global (baseline, region misfit only; profiled s_s 0.022, s_d 0.021): argmax
  0.000 in every scan, expected and realization modes, including the wide sync scan. The
  region-fit reference is the target's optimum when the templates and the SED are exact.
- exact_fixT_global (fixed T only; s_s 0.038, s_d 0.021, mu 1.014 / 1.016): sync +0.000,
  +0.010, +0.020, +0.030 (gains 0, 5, 15, 22 nats), dust +0.000, +0.010, +0.010, +0.010 (gains
  0, 112, 50, 3.5); wide sync scan at 0.000. Fixed T with an exact template shifts the hier
  target by up to 0.03 (sync faint quartile) and 0.01 (dust) relative to the flat-prior region
  fit: at the "unbiased" threshold for dust, above it for sync Q4.
- data_trueT_global (data-built template only; mu_s 1.14, s_s 0.68, mu_d 0.89, s_d 0.27): sync
  +0.000, +0.070, +0.050, -0.100 (edge; gains 146, 54, 981), dust +0.000, +0.050, +0.100,
  +0.100 (edges; gains 437, 644, 395); the wide sync scan peaks at -0.6 (beta_s about -3.6 in
  the faintest quartile, gain 3307 nats): the E020 runaway direction, reproduced as a block
  shift.
- data_fixT_global (the Commander4 configuration): the same numbers within 0.01 (sync +0.070,
  +0.060, -0.100; dust +0.050, +0.100, +0.100; wide -0.6). Realization mode identical to the
  expected mode within 0.01. Attribution at this depth: the data-built template under one
  global scatter is the driver; fixed T adds at most 0.03.
- data_fixT_perq (per-quartile scatter; profiled s_s [0.26, 0.32, 0.42, 3.44] against the
  realized sd(a/g) [0.25, 0.31, 0.42, 155], s_d [0.33, 0.22, 0.22, 0.27] against [0.33, 0.22,
  0.22, 0.40]; mu_s 0.91, mu_d 0.91): sync +0.000, +0.040, -0.050, -0.070 (gains 48, 52, 57),
  dust +0.000, +0.030, +0.060, +0.100 (gains 246, 228, 243); the wide sync scan peaks at 0.000
  (runaway removed: with s_4 = 3.4 the prior is inert in the faintest quartile). Not "unbiased":
  the shifts fall by a factor 1.5 to 2 in sync and dust Q2 and stay at 0.06 to 0.10 in the
  fainter dust quartiles, where the global mu (0.91 against realized means 0.89 to 0.99 per
  quartile) and the template's structure, not its scatter, carry the misfit.
- data_trueT_perq: as data_fixT_perq within 0.01 (sync +0.000, +0.040, -0.060, -0.080; dust
  +0.000, +0.030, +0.060, +0.100; wide 0.000): the fixed T plays no role once the template
  drives the fit.
- data_fixT_perqmu (mean and scatter per quartile; profiled mu_s [0.91, 0.83, 1.04, 2.93]
  against the realized means [0.89, 0.83, 1.04, 77], s_s [0.26, 0.31, 0.41, 2.41]; mu_d [0.90,
  0.93, 0.91, 0.84] against [0.89, 0.93, 0.94, 0.99], s_d [0.33, 0.22, 0.22, 0.26] against
  [0.33, 0.22, 0.22, 0.40]): sync +0.000, +0.010, -0.010, -0.060 (gains 1.7, 2.9, 40), dust
  +0.000, +0.050, +0.100, +0.100 (gains 422, 499, 77); wide 0.000. The per-quartile mean brings
  the sync Q2 and Q3 shifts to 0.01 (the mean of a/g per quartile, 0.83 to 1.04, was the sync
  misfit there) but not the faintest quartile (where a/g is not a distribution at all) and not
  the dust at all: the dust template, 353 GHz scaled with a fixed T and smoothed, has an error
  that varies with T_d and hence with the pixel's own index and position, which no scalar per
  quartile absorbs; the dust shifts under perqmu equal those under the global prior.

- data_fixT_perregmu (added 2026-10-29 21:27 after the perqmu result: one template multiplier
  per nside-16 region for each component, 3072 each with the N(1, 0.5^2) hyperprior, plus a
  scatter per quartile; `E023/scan_n003_perregmu.txt`; profiled region means per sync quartile
  0.91, 0.83, 1.03, 2.12 with sd 0.22 to 1.0, dust 0.90, 0.93, 0.91, 0.82; s_s [0.13, 0.08,
  0.16, 1.10], s_d [0.26, 0.18, 0.18, 0.18], below the realized scatters because the region
  means absorb part of them): sync +0.000, +0.000, +0.000, -0.100 (edge; wide -0.2, gain 179)
  and dust +0.000, +0.050, +0.100, +0.100 (gains 389, 624, 112). A free template amplitude per
  region fixes the sync everywhere the template is a template at all and leaves the dust
  shifts exactly where the global, per-quartile and per-quartile-mean laws left them: the dust
  template's error is a within-region, pixel-scale structure (the 353 GHz map scaled with one T
  carries the sky's T_d variation, which tracks the index itself), not an amplitude error at
  any region scale.

Reading of the 0.03 sky (pre-registered rules): the Commander4 configuration (data-built
template, fixed T, global scatter) is biased in every quartile but the brightest, by 0.05 to
above 0.1 in both indices with the faint-sync runaway; the driver is the template under one
global scatter, fixed T adds up to 0.03 (sync faint quartile) and 0.01 (dust) with an exact
template; the per-quartile scatter removes the runaway and nothing else of substance; the
per-quartile mean and scatter also fix the mid-S/N sync but leave the dust at 0.05 to 0.10.
Neither candidate is "unbiased" by the rule, and the dust failure is structural (a
template-error law correlated with the index, F23's mechanism in a second form), so neither
is carried to stage 2 as a fix; the result for Q16 is that a data-built template with any
scalar-hyperparameter law, and with a free multiplier per nside-16 region, is not usable for
dust at this depth: the error lives at the pixel scale and is tied to T_d, so the template must
come from a component separation with its own T_d (or the dust SED must carry a sampled T_d so
that the template error decorrelates from the index), not from one scaled frequency map.

**PySM sky at noise 1, first pass with the E021 CMB prior (`E023/scan_n1.txt`, 17:19 to
20:45):** the "stop and find it" rule fired. The baseline cell exact_trueT_global (exact
templates, true T; profiled s_s 0.032, s_d 0.023) peaks at +0.050 in dust Q2 (gain 7.7 nats,
curvature sd 0.013: 3.8 sd), at the +0.100 edge in dust Q3 and Q4 (9.4 and 2.9 nats) and at
-0.100 in sync Q4 (2.4 nats), all in the chi-square part (log-determinant and trace parts
below 0.2 nats). Three hypotheses were tested at 21:00: (H1) the flat-amplitude region fit is
not the pseudo-true value of a model whose amplitudes are pinned to the template: refuted, a
pinned-amplitude region fit (sync and dust amplitudes at the truth, CMB free) differs from the
flat fit by 0.001 in every dust quartile and 0.002 in sync (scratchpad
`e023_n1_pinned_refs.npz`); (H2) the determinant or trace terms: refuted by the part table;
(H4) the per-pixel CMB prior N(0, 100^2) of the E018/E021 target: confirmed as the cause by
the same cell with the CMB prior sd at 1e4 (`E023/scan_n1_baseline_cmbflat_check.txt`,
`--cmb-sd 1e4`), which peaks at +0.010 in dust Q2, Q3, Q4 (gains 0.02 to 0.04 nats) and at
-0.020 in sync Q4 (0.02 nats), the wide scan at 0.000. The mechanism was identified by the
stats review (independent code, `E023/f26_check.py`, `f26_shuffle.py`, `f26_scaling.py`,
reproduced, `f26_shuffle_output.txt`): not shrinkage but leakage of the fixed, spatially
correlated CMB realization through the finite per-pixel prior into the term of the expected
chi-square that is linear in the realization; the sign flips with the CMB's sign, rolling the
CMB over the sky gives shifts of either sign (dust Q3 from -0.10 to +0.10), shuffling its
pixels removes them (0.002 to 0.010), the leak fraction is 0.25 at 100 uK and the true
variance 84 uK makes it larger (0.058, 11 nats); the shift scales as noise^2 and is about
1e-4 at 0.03 noise, which is why every earlier control passed. Consequence (F26, D021): the
collapsed per-pixel likelihood at realistic noise must carry a flat CMB amplitude prior
(`amp_prior: none`, or sd 1e4), at the cost of the CMB's one degree of freedom per pixel
(faint dust block sd 20 to 25 percent wider); the writer gains `--cmb-prior-sd`, and the
noise-1 cells are rerun with it (`E023/scan_n1_cmbflat.txt`, from 21:04). The first-pass
noise-1 cells are kept as evidence but not read for attribution.

Noted from the first pass regardless of the CMB prior: at noise 1 the block curvature sds are
0.03 to 0.12 (sync quartiles) and 0.003 to 0.05 (dust), so the quartile scans resolve shifts of
0.1 only at 1 to 4 sd for the fainter quartiles; the data-built templates at this noise are not
templates in the faint half of the sky (realized a/g means 2.3 and 141 in sync Q3 and Q4, 5.9
and 48 in dust Q3 and Q4: the 23 and 353 GHz maps at this noise, scaled, smoothed and floored).

**PySM sky at noise 1 with the flat CMB prior (`E023/scan_n1_cmbflat.txt`, 21:04 to 00:09; CMB
prior sd 1e4, seven cells in expected mode, two also in realization mode; gains in nats, 2 nats
about 2 block sd):**
- exact_trueT_global (baseline): every scan within its curvature sd, gains 0.00 to 0.08 (sync
  -0.01 to -0.03 against sds 0.03 to 0.12; dust +0.01 against 0.016 to 0.053; wide 0.000). Clean.
  Realization mode of the same cell: shifts up to the 0.1 edge with gains 0.2 to 9 and the
  wide sync scan at +1.4 (68 nats): one noise realization at this depth does not constrain the
  faint sync quartile at all (its linearized block sd of 0.12 understates a likelihood that is
  flat over a unit of beta_s), which is why the expected scan is the primary line.
- exact_fixT_global: gains 0.00 to 0.40 (sync -0.01 to -0.07 against sds 0.03 to 0.12, dust
  0.00 to +0.02): fixed T is invisible at this depth.
- data_trueT_global and data_fixT_global (the Commander4 configuration; profiled mu_s 0.90 to
  0.91, s_s 0.35, mu_d 0.86 to 0.87, s_d 0.41): sync Q1 -0.090 (4.4 nats, 3 sd at sd 0.030), Q2
  +0.100 edge (8.4 nats), Q3 and Q4 at the -0.100 edge (1.0 and 2.6 nats, not significant), wide
  sync scan at -2.0 to -2.2 (33 nats: beta_s near -5 in the faintest quartile, the E020 runaway
  direction); dust Q1 +0.010 (5 to 6 nats, 3 sd at sd 0.0037), Q2 +0.100 edge (58 to 62 nats), Q3
  +0.100 edge (19 to 21 nats), Q4 -0.100 (1 nat, not significant). Realization mode within 0.02
  of the expected mode where the gain exceeds 2 nats. Fixed against true T: identical within
  0.01. The data-built template biases both indices where the data can still see the index
  (sync Q1 and Q2, dust Q1 to Q3) and runs the faintest sync quartile away.
- data_fixT_perq and data_trueT_perq: identical to the global cells to 0.01 in every scan
  (profiled s_q 0.33 to 0.50 against realized scatters of 0.28 to 212: at this depth the
  profile cannot see that the faint-quartile templates are not templates, and the index
  likelihood there is the prior's).
- data_fixT_perqmu (mu_s [0.93, 0.73, 0.93, 1.58], mu_d [0.90, 0.73, 0.64, 1.02]): sync Q1 -0.040
  (0.7 nats), Q2 to Q4 at the -0.100 edge with 0.6 to 1.0 nats (not significant), wide -0.6 (1.6
  nats, not significant); dust Q1 +0.030 (22 nats, 8 sd), Q2 +0.100 (5.6 nats), Q3 and Q4 not
  significant. The per-quartile mean removes every significant sync shift at this depth and
  creates a bright-quartile dust bias of +0.03 at 8 sd (a quartile-wide multiplier 0.90 with
  scatter 0.41 is the wrong law for a template whose error has pixel-scale structure tied to
  T_d, and the brightest quartile is where the data can tell).

Reading of the noise-1 sky (pre-registered rules): with the flat CMB prior the exact-template
cells are "unbiased at this depth" (every shift inside its block sd); the Commander4
configuration with the data-built template is biased in sync Q1 and Q2 (3 sd and beyond the
scan) and dust Q1 to Q3 (3 sd to beyond the scan) and runs the faintest sync quartile away;
fixed T is not a factor; the per-quartile scatter changes nothing; the per-quartile mean fixes
the sync and breaks the bright dust. No candidate is carried to stage 2.

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
Stage 1 (exact-target scans, no chains), on the PySM d1s1c1 sky at 0.03 x and at 1 x the E012
noise, with the E021 consistent sky as the model control (zero bias in every scan) and the
numpy marginal checked against Commander4's kernel (4e-12).

1. **The data-built template is the driver.** With templates built from one scaled, smoothed
   frequency map (23 GHz for sync, 353 GHz at a fixed T for dust) under the hier prior's law
   a ~ N(mu g, (s g)^2) with global hyperparameters, the collapsed index target is biased at both
   depths: at 0.03 x noise by 0.05 to beyond 0.1 in both indices in every quartile but the
   brightest, with the faintest sync quartile running toward beta_s of -3.6; at realistic noise
   by 0.09 (sync Q1, 3 sd), beyond 0.1 (sync Q2, dust Q2 and Q3), 0.01 (dust Q1, 3 sd), with the
   faintest sync quartile running to beta_s near -5. With exact templates the same target is
   unbiased at realistic noise and biased by at most 0.03 (sync faint quartile, fixed T) at 0.03
   x noise.
2. **Fixed T_d is a second-order effect here:** at most 0.03 (sync) and 0.01 (dust) at 0.03 x
   noise with an exact template, invisible at realistic noise, and identical within 0.01 to the
   true-T cells whenever the template drives the fit.
3. **No scalar-hyperparameter law repairs a single-frequency dust template.** A scatter per
   S/N quartile removes the sync runaway at 0.03 x noise and nothing else; a mean and scatter
   per quartile fixes the mid-S/N sync (to 0.01 at 0.03 x, to insignificance at noise 1) and
   leaves the dust at 0.05 to 0.10 (0.03 x) or moves the bias into the brightest dust quartile
   (+0.03 at 8 sd, noise 1); a free template multiplier per nside-16 region (3072 per
   component) fixes the sync wherever the template is a template and leaves the dust untouched.
   The dust template's error is pixel-scale and tied to T_d (hence to the index), and the
   profiled hyperparameters cannot see it where the data are weak.
4. **The per-pixel CMB prior leaks at realistic noise (F26).** The finite per-pixel CMB prior of
   the E018/E021 target leaks the fixed, correlated CMB realization into the index likelihood
   (dust mid quartiles shifted by 0.05 to 0.1 with a realization-dependent sign, at several of
   the model's own block sd); the flat CMB prior removes it exactly at the cost of one degree of
   freedom per pixel (D021). Inert at 0.03 x noise (about 1e-4).
5. **What this decides for Q16 (D022):** the hier amplitude prior is usable for synchrotron with
   a template multiplier at the region scale (a template-error field), and not for dust with a
   template from one scaled frequency map at any depth tested; the dust template must come from
   a component separation with its own T_d (its error then decorrelates from the pixel's T_d),
   or T_d must be sampled so that the fixed-T scaling is not what the template encodes; both are
   stage-2 candidates, tested first by the same scans with a template equal to the true dust
   amplitude times a smooth error field (the "component-separation-quality" template) and with
   the dust SED carrying the true T_d. At realistic noise the faint half of the sky is
   prior-dominated for both indices at nside 16 (block sds 0.08 to 0.12 for sync, 0.03 to 0.05
   for dust), so what the prior says there is what the posterior says; a wrong template is not
   a small effect at this depth.

Not shown: the posterior bias under the GMRF index prior (stage 2), sampled T (only fixed
against known T), correlated template errors of a component-separation kind, seed dependence
(one realization per cell; the expected-log-likelihood scans remove the noise scatter of the
peak, not the template realization), and anything beyond one sky. F27.

## Follow-ups
- Stage 2 design (E025): the same scans with (i) a component-separation-quality dust template
  (true amplitude times a smooth 30 percent error field at 5 to 10 deg, no T_d coupling), (ii)
  the dust SED with the true T_d per pixel as the stand-in for sampled T, (iii) sync with the
  per-region multiplier; then the gmrf16 chain on the noise-1 sky with the flat CMB prior and
  the dense metric, E021 metrics and the Q16 acceptance.
- Sampled T_d as a second spectral parameter per region in the kernel (M2.6 scope).
- Realistic-noise faint sky: the index posterior there is the prior's; Paper A states this as
  the regime boundary (block sd per quartile from this record).
