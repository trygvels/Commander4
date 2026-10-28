# E022: BeyondPlanck XV polarized beta_s: amplitude-index Gibbs coupling on a low-resolution surrogate

- **Date started / finished:** 2026-10-28 /
- **Branch and commit:** nextgen, 90720f1
- **Open question addressed:** F14 scope (planned E014 follow-up); feeds Paper A Sect. 3. No Q
  number of its own; it bounds what Paper A may say about BeyondPlanck XV.
- **Validation level:** L0 (exact 1-D targets on real data; the Gibbs chain is checked against them)

## Question
On the BeyondPlanck XV configuration for the polarized synchrotron index (four sky regions, LFI
30 to 70, WMAP Ka, Q, V and HFI 353 in Q and U, prior N(-3.3, 0.1), amplitudes held fixed during
the index step), what is the ratio of the collapsed (amplitude-marginalized) posterior sd of
beta_s to the conditional sd at fixed amplitudes, per region, and what Gibbs autocorrelation
length does it imply? BP XV reports a correlation length of 30 to 50 Gibbs samples over about 500
samples and sds of 0.06 (Spur) and 0.07 (Plane). If the implied length is of that order, the BP XV
uncertainties stand and F14 (ratio 180 for a full-sky intensity index at nside 64) does not
transfer; if it is far larger, the BP XV error bars are underestimates and Paper A says so with
numbers. Either way the answer is measured, not assumed.

## Target and priors
All targets are evaluated on one fixed low-resolution surrogate of the BP XV data (Setup), on
HEALPix nside 16, Stokes Q and U, seven bands, internal units uK_RJ. Model per pixel p and Stokes
s: d_{nu,p,s} = A_{nu}(beta_{r(p)}) a_{p,s} + n_{nu,p,s}, with a_{p,s} = (synchrotron, thermal
dust, CMB) amplitudes, A_nu the 7 x 3 mixing matrix: synchrotron (nu / 30 GHz)^{beta_s},
thermal dust a modified blackbody with beta_d = 1.6234 and the BP per-pixel T_d (both from the
BP release dust map, fixed), reference 353 GHz, CMB flat in uK_CMB converted to uK_RJ per band.
Mixing-matrix columns use delta bandpasses at the BP nominal frequencies (28.4, 44.1, 70.1,
33.0, 40.6, 60.8, 353 GHz). Noise is Gaussian, independent between pixels, Stokes and bands,
variance sigma^2_{nu,p,s} (Setup). Regions r are the four BP XV regions (High Latitude, Spur,
Galactic Center, Plane) recovered from the released BP polarized beta_s map; the targets are
separable over regions, so each is a 1-D distribution per region.

- **T1, conditional (the BP XV index step):** P(beta_r | d, a) propto exp(-chi2_r(beta_r) / 2)
  N(beta_r; -3.3, 0.1^2), chi2_r = sum over p in r, s, nu of (d - A(beta_r) a)^2 / sigma^2, with
  a fixed. Amended 2026-10-28 before the main result was read: a is drawn from P(a | beta*, d)
  (beta* = the T2 mean with prior for Spur and Plane, the prior mean for the prior-only
  regions), and the width is averaged over 20 draws, because the Gibbs coupling depends on the
  conditional at amplitudes consistent with the data; T1 at the BP release posterior-mean maps
  is kept as a diagnostic of the surrogate's consistency with those maps only (the first run
  showed it pinned at the grid edge: the release maps were fit with bandpass-integrated SEDs,
  the full-mission 353 map and native resolution, and a conditional of width 0.02 resolves
  per-mille mismatches). Reported with and without the prior factor.
- **T2, collapsed, flat amplitude prior:** P(beta_r | d) propto prod over p in r, s of
  [ det(M_{p,s})^{-1/2} exp(b^T M^{-1} b / 2) ] N(beta_r; -3.3, 0.1^2), with
  M = A^T N^{-1} A (3 x 3) and b = A^T N^{-1} d evaluated at beta_r; the -1/2 log det M term is
  the flat-prior amplitude marginal (Commander4 `marginal` mode). The profile variant without
  the log det term (Commander4 `ridge` in the limit of a wide ridge) is reported next to it.
  Reported with and without the prior.
- **T3, two-block Gibbs on the surrogate:** alternate a_{p,s} ~ N(M^{-1} b, M^{-1}) for every
  pixel and Stokes (independent given beta) and beta_r ~ T1 at the current a (exact inverse-CDF
  draw on a grid, no Metropolis tuning), for all four regions jointly. Its stationary
  distribution for beta is T2 (flat amplitude prior) by construction; its integrated
  autocorrelation time tau_r for beta_r is the measured Gibbs coupling length. The Gaussian
  two-block prediction tau = 2 r^2 - 1 with r = sd_T2 / sd_T1 (lag-1 autocorrelation
  1 - 1/r^2) is reported next to the measured tau.

Prior-only regions: BP XV sampled High Latitude and Galactic Center from the prior alone; those
regions carry no Gibbs coupling and are reported for completeness only.

## Setup
Surrogate built by `nextgen/scripts/prep_bpxv_surrogate.py` into `~/data/bpxv_surrogate/`:
- LFI 30, 44, 70: BP release frequency maps (`BP_0xx_IQU_n0512/1024_v2.fits`, Q_MEAN, U_MEAN,
  uK_CMB), smoothed from the LFI beams (32.29, 27.0, 13.25 arcmin) to 60 arcmin, degraded to
  nside 16; noise sd from the BP low-resolution noise maps
  `BP7_0xx_rms_n0016_60arcmin_uK_true.fits` (BP's own nside-16, 60 arcmin rms used for its
  low-resolution index sampling; they include correlated noise).
- WMAP Ka, Q, V: the nine-year nside-16 band polarization maps and dense QU inverse noise
  covariances from LAMBDA (`wmap_band_qumap_r4_9yr_{Ka,Q,V}_v5.fits`,
  `wmap_band_quninv_r4_9yr_{Ka,Q,V}_v5.fits`), the data BP XV used (Sect. 3); maps as provided.
  For T1 to T3 the noise is the diagonal of N = (N^{-1})^{-1}; the dense N^{-1} is used in a
  conditional-only check (T1 with the dense quadratic form) to show what the diagonal
  approximation does to the conditional width.
- HFI 353: NPIPE A split map (`npipe6v20A_353_map_uK.fits`, uK_CMB, nside 2048) smoothed from
  4.94 to 60 arcmin and degraded; noise sd from 60 white-noise realizations drawn from the
  `mcscaled` rms map put through the same operations (BP used the full-mission NPIPE map; the
  A split has about 1.4 times its noise, which only weakens the dust constraint).
- Components held fixed in T1: BP release maps `BP_synch_IQU_n1024_v2.fits` (Q_MEAN, U_MEAN,
  uK_RJ at 30 GHz, 60 arcmin), `BP_dust_IQU_n1024_v2.fits` (Q, U at 353 GHz, 10 arcmin,
  smoothed to 60; QU_T_MEAN for T_d; QU_BETA_MEAN = 1.6234), `BP_cmb_IQU_n1024_v2.fits` (Q, U,
  uK_CMB, 14 arcmin, smoothed to 60), all degraded to nside 16.
- Regions: from `QU_BETA_MEAN` and `QU_BETA_STDDEV` of the BP synchrotron map at nside 16
  (plateaus: prior-only -3.30 with sd 0.10; Spur about -3.17 with sd 0.06; Plane about -3.00
  with sd 0.05; the 10 degree smoothing BP applied to the region map blurs the boundaries, and
  the nearest-plateau assignment is a surrogate for BP's smoothed-boundary model). Galactic
  Center = prior-only pixels within 25 degrees of the centre; the rest of the prior-only pixels
  are High Latitude.
- No mask beyond non-finite pixels (BP XV's processing mask is not on disk; the regions already
  separate the plane).
- Analysis `nextgen/scripts/bpxv_gibbs_coupling.py`: grids of 2001 points on [-4.5, -1.5] per
  region for T1 and T2; Gibbs chain of 4000 iterations for T3 with tau from the
  initial-positive-sequence estimator and from blocking; machine Mac mini, one process, while
  E021 runs (the surrogate is 7 x 2 x 3072 numbers; seconds per grid, a minute per chain).

## Reference
T1 and T2 are one-dimensional and computed on a grid (exact to the grid spacing of 0.0015, far
below the sds of 0.01 to 0.1 expected). T3's reference is T2: the chain mean and sd of each
beta_r must match the T2 grid within Monte Carlo error, which checks the chain before its tau is
used. tau is compared with the Gaussian two-block formula and with BP XV's reported 30 to 50.
The conditional with dense WMAP N^{-1} is an independent re-evaluation of T1 (no grid
assumption beyond the same beta grid).

## Command
```bash
python -I nextgen/scripts/prep_bpxv_surrogate.py --out ~/data/bpxv_surrogate
python -I nextgen/scripts/bpxv_gibbs_coupling.py --dir ~/data/bpxv_surrogate --iters 4000
```

## Result
**Side measurement for Q7 (2026-10-28, `nextgen/scripts/q7_lfi_rms_check.py`):** BP's nside-16,
60 arcmin LFI noise rms from noise simulations (`BP7_0xx_rms_n0016_60arcmin_uK_true.fits`, the
noise description BP used for its low-resolution index sampling) against white-noise propagation
of the BP release Q_RMS/U_RMS maps through the same smoothing and degrading (40 realizations):

| Band | rms_true median (uK) | white-propagated median (uK) | ratio median | ratio 10 to 90 percent |
|---|---|---|---|---|
| 30 Q / U | 1.649 / 1.648 | 1.023 / 1.019 | 1.59 / 1.60 | 1.37 to 1.87 |
| 44 Q / U | 2.147 / 2.075 | 1.293 / 1.227 | 1.64 / 1.66 | 1.40 to 1.95 |
| 70 Q / U | 1.712 / 1.708 | 1.062 / 1.051 | 1.59 / 1.60 | 1.38 to 1.87 |

Correlated noise raises the LFI noise variance at this resolution by a factor 2.5 to 2.8 over
white propagation, everywhere on the sky (ratio 1.4 to 1.9 between the 10th and 90th
percentiles, maximum 2.9). A per-pixel likelihood at nside 16 built from white-noise propagation
of the native rms maps would therefore understate the LFI noise sd by 1.6 and overstate the
index precision by the same factor. The surrogate uses BP's rms_true for LFI, so this bias is
not in E022; the 353 noise is white-propagated (NPIPE A, dust at S/N above 15 at nside 16, so
it hardly matters for beta_s), and the WMAP noise is the diagonal of the dense covariance.

**Main result (2026-10-28, 10:56; `result.json`, `posteriors_and_chains.png`, chains
`gibbs_chain_{bp,all}.npy`):** surrogate with 7 bands, 3072 pixels, regions High Latitude 2194,
Spur 464, Galactic Center 91, Plane 323 pixels (BP map means there -3.295, -3.183, -3.294,
-3.022; sd map 0.098, 0.064, 0.095, 0.067, matching the paper's prior-only / 0.06 / prior-only
/ 0.07 pattern). Grid spacing 0.0015.

| Region | T1 conditional, sd (mean) | T1 with dense WMAP N^-1, sd | T2 collapsed with prior, mean +- sd | profile with prior, mean +- sd | ratio sd_T2 / sd_T1 | Gaussian tau = 2r^2 - 1 | Gibbs tau, BP scheme | Gibbs tau, all regions |
|---|---|---|---|---|---|---|---|---|
| Spur | 0.0198 (-2.853) | 0.0199 | -2.8725 +- 0.0466 | -3.1431 +- 0.0454 | 2.35 | 10.1 | 14.0 +- 1.8 | 10.2 +- 1.1 |
| Plane | 0.0184 (-3.109) | 0.0187 | -3.1094 +- 0.0389 | -3.2349 +- 0.0387 | 2.11 | 7.9 | 8.1 +- 0.8 | 7.0 +- 0.6 |
| Galactic Center | 0.0130 (-3.134) | 0.0133 | -2.6379 +- 0.0257 | -2.6581 +- 0.0257 | 1.99 | 6.9 | prior-only | 17.1 +- 2.4 |
| High Latitude | 0.0168 (-3.044) | 0.0172 | -1.5013 +- 0.0019 (grid edge) | -2.9523 +- 0.0360 | (edge) | (edge) | prior-only | 3.1 +- 0.2 at the edge |

T1 sds vary by 0.0002 over the 20 amplitude draws (their means by 0.008 to 0.02). The Gibbs
chains (4000 iterations, burn 400) reproduce T2: Spur -2.8759 / sd 0.0465 against
-2.8725 +- 0.0466, Plane -3.1124 / 0.0386 against -3.1094 +- 0.0389 (BP scheme; the all-regions
chain -2.8762 / 0.0455 and -3.1163 / 0.0385). Effective samples 257 to 511 of 3600. The
all-regions chain is a second realization for Spur and Plane and adds the prior-only regions:
High Latitude runs from -3.3 to the grid edge at -1.5 within 50 iterations under the flat
amplitude prior (the collapsed target there is at the edge: the -1/2 log det M term of 4388
pixel-Stokes marginals dominates a region with per-pixel synchrotron S/N below 1), while its
conditional (-3.04 +- 0.017) and profile (-2.95 +- 0.036) are regular. The Galactic Center
collapsed target sits at -2.64 +- 0.026 with the conditional at -3.13 +- 0.013 and tau 17.
Prior-to-posterior width ratios (prior 0.1 over sd_T2 with prior): Spur 2.1, Plane 2.6, GC 3.9.

Surrogate consistency: the conditional at the BP release posterior-mean amplitudes pins at
the grid edge -1.5 in every region (chi2 strongly prefers the flattest synchrotron allowed),
so the release amplitudes do not fit the surrogate at the level a conditional of width 0.02
resolves; expected (bandpass integration, native resolution, full-mission 353, the BP
10-degree-smoothed index map) and the reason T1 is evaluated at draws. The surrogate's own
T2 means (Spur -2.87, Plane -3.11) are not comparable to BP XV's -3.17 +- 0.06 and
-3.03 +- 0.07 for the same reasons and because the flat-prior volume term scales with the
number of amplitude pixels (Follow-ups); the profile values (-3.14, -3.23) happen to lie close
to the paper's, which is noted, not interpreted.

## Negative control
The intensity analogue of F14: the same scripts on the E014 nside-64 intensity configuration are
not rerun here; instead the Gaussian formula is checked on T3 itself (a chain whose tau disagrees
with 2 r^2 - 1 by more than its error flags a non-Gaussian 1-D posterior or a chain bug).

## Conclusion
(draft, stats review of the targets and scripts pending)
- On the BP XV polarization configuration the amplitude-index Gibbs coupling is weak: the
  collapsed posterior sd of beta_s is 2.1 to 2.4 times the conditional sd in the two
  likelihood-sampled regions, and a two-block Gibbs sampler with exact conditional draws has
  integrated autocorrelation times of 7 to 14 iterations, as the Gaussian two-block formula
  predicts (8 to 10). F14's ratio of 180 (full-sky intensity index at nside 64, per-pixel
  amplitude S/N in the hundreds) does not transfer to this configuration, where the polarized
  synchrotron amplitude per nside-16 pixel has S/N of 2 to 3 at 30 GHz (|Q| median 4.4 uK_RJ
  against sigma 1.6).
- BP XV's reported Gibbs correlation length of 30 to 50 is larger than the 7 to 14 the
  amplitude-index coupling alone produces here, so there is no evidence that it understates the
  coupling; the remainder is attributable to the Metropolis step within the conditional and the
  other Gibbs blocks (gain, bandpass, noise), which this cell does not model. With about 500
  samples and tau 30 to 50 the paper had 10 to 17 independent index samples, so its sds of 0.06
  and 0.07 carry a Monte Carlo uncertainty of about 20 percent, not a coupling bias. Paper A
  Sect. 3 states this instead of extrapolating F14.
- The flat-amplitude-prior marginal (the stationary distribution of a Gibbs sampler with flat
  amplitude priors and a conditional index step) is shifted from the profile by +0.27 (Spur),
  +0.13 (Plane), +0.02 (GC) and runs to the flattest allowed index for the High Latitude
  region on this real-data surrogate: the prior-volume effect of Paper A Sect. 4 at a
  publication-relevant configuration, and a direct mechanism for the "prior-dominated" regions
  BP XV sampled from the prior alone. The quantitative shifts belong to this surrogate (nside 16
  amplitudes); at BP's native resolution the volume term is larger by the pixel ratio.
- Limits: diagonal noise for the collapsed target (the dense WMAP covariance changes the
  conditional width by under 2 percent); delta bandpasses; NPIPE A split at 353; regions
  recovered from the smoothed release map; nside 16 instead of multi-resolution; the ratio r is
  approximately resolution-independent for white noise and smooth signals (both Fisher
  quantities scale with the summed inverse variance), the volume term is not.

## Follow-ups
- Resolution scaling of the flat-prior volume term: repeat T2 at nside 32 and 64 (LFI and 353
  only, WMAP excluded) to measure how the collapsed means move with the number of amplitude
  pixels; this is the quantitative link to BP's native-resolution amplitude step and to Paper A
  Sect. 4.
- Dense WMAP covariance in the collapsed target (18432-amplitude Gaussian marginal per beta
  grid point) if the diagonal approximation moves the conditional width by more than 20 percent.
- Bandpass-integrated mixing matrix from the BP bandpass files instead of delta bandpasses.
