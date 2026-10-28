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
  a fixed at the BP release posterior-mean component maps. Reported with and without the prior
  factor (the "without" sd is the conditional likelihood width).
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
(pending)

## Negative control
The intensity analogue of F14: the same scripts on the E014 nside-64 intensity configuration are
not rerun here; instead the Gaussian formula is checked on T3 itself (a chain whose tau disagrees
with 2 r^2 - 1 by more than its error flags a non-Gaussian 1-D posterior or a chain bug).

## Conclusion
(pending)

## Follow-ups
- Dense WMAP covariance in the collapsed target (18432-amplitude Gaussian marginal per beta
  grid point) if the diagonal approximation moves the conditional width by more than 20 percent.
- Bandpass-integrated mixing matrix from the BP bandpass files instead of delta bandpasses.
