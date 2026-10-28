# E022: BeyondPlanck XV polarized beta_s: amplitude-index Gibbs coupling on a low-resolution surrogate

- **Date started / finished:** 2026-10-28 / 2026-10-28
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

Variants (added after the stats review, 2026-10-28): the CMB amplitude prior is either flat
(the released BP configuration, `COMP_CL_TYPE01 = none`) or Gaussian N(0, sigma_cmb^2) per
pixel and Stokes with sigma_cmb = 0.187 uK_CMB from the Planck 2018 EE spectrum at this
resolution (then M gains the prior precision on the CMB diagonal in T2 and T3); the LFI noise
sd is either BP's simulated nside-16 rms or white-noise propagation. All four combinations are
reported; the flat-CMB, BP-rms case is the primary one and the others bound the dependence.

## Setup
Surrogate built by `nextgen/scripts/prep_bpxv_surrogate.py` into `~/data/bpxv_surrogate/`:
- LFI 30, 44, 70: BP release frequency maps (`BP_0xx_IQU_n0512/1024_v2.fits`, Q_MEAN, U_MEAN,
  uK_CMB), smoothed from the LFI beams (32.29, 27.0, 13.25 arcmin) to 60 arcmin, degraded to
  nside 16; noise sd from the BP low-resolution noise maps
  `BP7_0xx_rms_n0016_60arcmin_uK_true.fits` (BP's simulated nside-16, 60 arcmin rms, which
  include correlated noise; correction after the stats review: BP XV's index step evaluated
  Eq. 26 at native resolution with the per-pixel white rms, so these maps are not what BP XV
  weighted with, and a white-noise variant from the release Q_RMS/U_RMS maps propagated through
  the same smoothing is run next to them).
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
  region for T1 and T2; Gibbs chains of 6000 iterations for T3 (two per variant) with tau from
  the initial-positive-sequence estimator and its lag-1 checked against the exact identity
  1 - E[var_cond]/var; machine Mac mini, one process, while E021 runs (the surrogate is
  7 x 2 x 3072 numbers; seconds per grid, a minute per chain).

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
60 arcmin LFI noise rms from noise simulations (`BP7_0xx_rms_n0016_60arcmin_uK_true.fits`,
BP's low-resolution noise description for local index sampling at smoothing scale 3, which the
released polarized beta_s step did not use: it ran the alm sampler at native resolution)
against white-noise propagation
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

**Main result (2026-10-28, 11:40, after the stats review; `result_<variant>.json`,
`posteriors_and_chains_{bp_cmb0,bp_cmb0.187}.png`, chains `gibbs_chain_{bp,all}_<variant>.npy`):**
surrogate with 7 bands, 3072 pixels, regions High Latitude 2194, Spur 464, Galactic Center 91,
Plane 323 pixels (BP map means there -3.295, -3.183, -3.294, -3.022; released-map sd plateaus
0.098, 0.064, 0.095, 0.067; the paper quotes 0.06 for the Spur and 0.07 for the Plane, the
released map's Plane plateau is 0.05 to 0.07). Grid spacing 0.0015. Four variants, requested by
the review: LFI noise from BP's simulated nside-16 rms ("BP sim. rms") or from white-noise
propagation of the release rms maps ("white rms", what a native-resolution white-rms chi2 as in
BP XV Eq. 26 effectively weights with; the two differ by the Q7 factor 1.6), crossed with the
CMB amplitude prior: flat per pixel (the prior of the released BP configuration:
`COMP_CL_TYPE01 = none` in `BP_param_c0001_v2.txt`, read as no power-spectrum prior on the CMB
in the Gibbs chain, consistent with `BP_mean_v2.h5` carrying no `cmb/Dl_*` datasets while the
foregrounds do; the parser semantics of `none` are not verified against the Commander3 source;
the synchrotron and dust amplitudes carry Gaussian smoothing priors at 30 and 10 arcmin, which
do not act on nside-16 amplitudes; the resolution is the surrogate's, not BP's) or Gaussian N(0, 0.187^2 uK_CMB) per pixel and Stokes, the CMB Q/U rms at
nside 16 and 60 arcmin from the Planck 2018 EE spectrum (CAMB, l <= 200, beam and pixel window;
BB negligible). Chains: 6000 iterations, burn 600, exact inverse-CDF index draws, prior-only
regions drawn from the prior in the "BP scheme" chain and from the likelihood in the "all" chain.

Columns: T1 sd (conditional with prior, mean over 20 amplitude draws; the dense WMAP N^-1
changes it by at most 5 percent); T2 and profile means and sds with prior; r = sd_T2 / sd_T1;
Gaussian 2 r^2 - 1; measured tau (BP scheme / all chain); lag-1 measured / exact identity
1 - E[var_cond]/var; chain mean with its Monte Carlo error and the offset from the T2 mean in
units of that error.

| Variant | Region | T1 sd | T2 mean +- sd | profile mean +- sd | r | 2r^2-1 | tau (bp / all) | lag-1 meas. / identity | chain mean (offset) |
|---|---|---|---|---|---|---|---|---|---|
| BP sim. rms, flat CMB prior, nside 16 | Spur | 0.0198 | -2.872 +- 0.0466 | -3.143 +- 0.0454 | 2.35 | 10.1 | 11.3 +- 1.9 / 9.4 +- 1.2 | 0.817 / 0.819 | -2.8709 +- 0.0021 (+0.7) |
| BP sim. rms, flat CMB prior, nside 16 | Plane | 0.0184 | -3.109 +- 0.0389 | -3.235 +- 0.0387 | 2.11 | 7.9 | 8.0 +- 0.9 / 7.3 +- 0.8 | 0.780 / 0.777 | -3.1111 +- 0.0015 (-1.1) |
| BP sim. rms, flat CMB prior, nside 16 | High Latitude (prior-only in BP) | 0.0168 | -1.501 +- 0.0019 (grid edge) | -2.952 +- 0.0360 | 0.11 | -1.0 | - / 2.3 +- 0.2 | 0.372 / 0.405 | -1.5016 +- 0.0000 (-7.8) |
| BP sim. rms, flat CMB prior, nside 16 | Galactic Center (prior-only in BP) | 0.0130 | -2.638 +- 0.0257 | -2.658 +- 0.0257 | 1.99 | 6.9 | - / 16.0 +- 3.7 | 0.843 / 0.845 | -2.6347 +- 0.0014 (+2.4) |
| white rms, flat CMB | Spur | 0.0156 | -2.978 +- 0.0361 | -3.150 +- 0.0362 | 2.31 | 9.7 | 10.9 +- 1.8 / 8.9 +- 1.1 | 0.812 / 0.815 | -2.9769 +- 0.0016 (+0.8) |
| white rms, flat CMB | Plane | 0.0142 | -3.165 +- 0.0305 | -3.245 +- 0.0306 | 2.15 | 8.2 | 8.1 +- 0.9 / 7.5 +- 0.8 | 0.788 / 0.785 | -3.1664 +- 0.0012 (-1.5) |
| white rms, flat CMB | High Latitude (prior-only in BP) | 0.0131 | -2.439 +- 0.0301 | -3.034 +- 0.0281 | 2.30 | 9.6 | - / 19.5 +- 2.8 | 0.916 / 0.918 | -2.4410 +- 0.0018 (-1.0) |
| white rms, flat CMB | Galactic Center (prior-only in BP) | 0.0101 | -2.756 +- 0.0204 | -2.768 +- 0.0203 | 2.01 | 7.1 | - / 14.8 +- 3.5 | 0.831 / 0.833 | -2.7536 +- 0.0010 (+2.4) |
| BP sim. rms, CMB prior 0.19 uK | Spur | 0.0231 | -3.197 +- 0.0236 | -3.193 +- 0.0238 | 1.02 | 1.1 | 1.1 +- 0.1 / 1.1 +- 0.0 | 0.026 / 0.025 | -3.1965 +- 0.0003 (+0.0) |
| BP sim. rms, CMB prior 0.19 uK | Plane | 0.0189 | -3.177 +- 0.0195 | -3.176 +- 0.0196 | 1.03 | 1.1 | 1.1 +- 0.0 / 1.1 +- 0.1 | 0.055 / 0.029 | -3.1770 +- 0.0003 (-0.6) |
| BP sim. rms, CMB prior 0.19 uK | High Latitude (prior-only in BP) | 0.0172 | -3.104 +- 0.0177 | -3.098 +- 0.0181 | 1.03 | 1.1 | - / 1.1 +- 0.0 | 0.040 / 0.079 | -3.1040 +- 0.0003 (-1.0) |
| BP sim. rms, CMB prior 0.19 uK | Galactic Center (prior-only in BP) | 0.0120 | -2.974 +- 0.0123 | -2.973 +- 0.0123 | 1.02 | 1.1 | - / 1.1 +- 0.0 | 0.027 / 0.033 | -2.9730 +- 0.0002 (+3.2) |
| white rms, CMB prior 0.19 uK | Spur | 0.0173 | -3.189 +- 0.0179 | -3.194 +- 0.0180 | 1.03 | 1.1 | 1.2 +- 0.1 / 1.1 +- 0.0 | 0.043 / 0.046 | -3.1888 +- 0.0003 (+0.2) |
| white rms, CMB prior 0.19 uK | Plane | 0.0146 | -3.226 +- 0.0152 | -3.228 +- 0.0153 | 1.04 | 1.2 | 1.2 +- 0.1 / 1.2 +- 0.1 | 0.078 / 0.053 | -3.2258 +- 0.0002 (-1.0) |
| white rms, CMB prior 0.19 uK | High Latitude (prior-only in BP) | 0.0133 | -3.158 +- 0.0137 | -3.176 +- 0.0140 | 1.03 | 1.1 | - / 1.1 +- 0.0 | 0.061 / 0.101 | -3.1586 +- 0.0002 (-0.7) |
| white rms, CMB prior 0.19 uK | Galactic Center (prior-only in BP) | 0.0095 | -3.030 +- 0.0097 | -3.030 +- 0.0097 | 1.01 | 1.1 | - / 1.1 +- 0.0 | 0.027 / 0.032 | -3.0296 +- 0.0001 (+3.3) |

Reading of the table:
- With the flat CMB prior (released configuration) the collapsed-to-conditional ratio is 2.1
  to 2.35 for Spur and Plane and the Gibbs coupling length 7 to 11 iterations, independent of
  the LFI noise model (the noise model changes the widths by 1.3, not the ratio); the Gaussian
  two-block formula predicts 8 to 10 and the measured lag-1 equals the exact identity to 0.005.
  The High Latitude region's collapsed target sits at the grid edge with the BP noise and at
  -2.44 +- 0.03 with white noise, with tau 20 there; the Galactic Center has tau 15 to 16
  against a Gaussian prediction of 7 (its 1-D posterior is skewed; the chain mean sits 2.4 to
  3.3 Monte Carlo errors from the grid mean, the largest offset in the table).
- With a C_l-sized Gaussian prior on the CMB amplitude the coupling disappears: r = 1.02 to
  1.04, tau 1.1 in every region and both noise models, and the High Latitude and Galactic Center
  targets become regular (-3.10 +- 0.02 and -2.98 +- 0.01 with the BP noise). The
  amplitude-index coupling in this configuration is carried by the per-pixel CMB amplitude
  (0.19 uK of true signal against 1.5 to 2 uK of noise per pixel and Stokes), as the reviewer's
  toy found: with a flat prior the CMB column absorbs synchrotron and the marginal over it
  both widens the index posterior and shifts it (Spur -2.87 against -3.20 with the prior; the
  profile, which does not integrate over the CMB, sits at -3.14 in both cases).
- Chain validation: for Spur and Plane the chain means are within 0.7 to 1.5 Monte Carlo errors
  of the grid means in all variants and the sds agree to 1 percent (relative MC error 0.01 to
  0.03); effective samples 476 to 5100 of 5400. In the CMB-prior variants the measured lag-1
  and the identity differ by up to 0.04 for High Latitude (MC error about 0.014), harmless at
  tau 1.1.
- Resolution (stats re-check, 2026-10-28): both the flat-prior stationary distribution and the
  coupling depend on the number of free amplitude elements. The reviewer's two-component toy
  with the same total information and each pixel split into K subpixels (flat CMB prior, prior
  N(-3.3, 0.1), truth -3.1) gives T2 means -3.03 (K = 1), -2.80 (K = 4) and the prior edge for
  K >= 16, with tau 9, 13 and pinned. BP's chain has about 50 times the surrogate's synchrotron
  elements per region (nside 1024 with a 30 arcmin smoothing prior) and about 4000 times its
  CMB elements (lmax 2000, no prior), so the surrogate's r and tau are lower bounds on the
  flat-prior Gibbs coupling at BP's resolution, not estimates of it; the released chain is the
  only direct evidence (next paragraph).

**The released BP chain itself (2026-10-28, 12:20; `nextgen/scripts/bp_chain_beta_s.py`,
`bp_chain_analyze.py`; data `~/data/beyondplanck/chain_beta_s/c000{1,2}_synch_pixreg.npz`):**
every sample of `BP_c0001_v2.h5` (1001 samples, 1.8 TB file read by HTTP range requests,
613 MB fetched) carries `synch/beta_pixreg_val` of shape (2, 4): intensity (a single value,
the other three unused) and polarization (four regions). The parameter file confirms the
polarized step ran the alm sampler on native-resolution bands with per-pixel rms
(`COMP_BETA_POL_LMAX02 = 100`, `ALMSAMP_PIXREG`, `COMP_BETA_POL_FIX_PIXREG02 = 1,3`, pixreg
prior -3.3; WMAP bands 012 to 014 are the r4 maps with the `QUcov` dense N^-1 used here;
353 excluded from the index chi2 by `COMP_BETA_NU_MAX02 = 80`). For a two-block scheme with
exact conditional draws under stationarity var(delta beta)/2 = E[var_cond]; BP's step is a
1000-proposal Metropolis walk of length 6e-4 inside a conditional of sd at least 0.03 (reach
at most 0.019 rms per Gibbs step), which does not decorrelate within the conditional, so the
increment width is a lower bound on the conditional sd, while the other Gibbs blocks (gains,
bandpasses, correlated noise) move the conditional's centre between iterations and push it
the other way; the two cannot be separated from the chain. For any stationary series
var(delta x) = 2 var(x)(1 - rho_1), so "r from the increments" is 1 / sqrt(1 - rho_1), lag-1
restated under the two-block reading, not an independent measurement (the identity was a real
test on the surrogate, where E[var_cond] came from the grid). Both chains, burn 10
(half-chain difference in units of the MC error of a half mean times sqrt 2; trend from a
linear fit with n_eff = n / tau; both chains start from the same state):

| Chain, region | mean | sd (detrended) | increment width sqrt(var(delta)/2) | 1/sqrt(1 - rho_1) | lag-1 rho_1 | tau IPS | AR(1) tau from lag-1 | ESS | half-chain diff (MC err) | trend per 1000 samples | acf at lags 5, 10, 20, 50, 100 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1, 0 (prior-only) | -3.303 | 0.104 (0.104) | 0.103 | 1.01 | 0.03 | 1.1 | 1.1 | 939 | -0.003 (-0.4) | -0.006 +- 0.012 | 0 within 0.04 |
| 1, 1 (Spur) | -3.152 | 0.069 (0.066) | 0.029 | 2.43 | 0.83 | 22.5 +- 11 | 10.8 | 44 | +0.032 (+1.5) | +0.073 +- 0.037 | 0.44, 0.24, 0.17, 0.10, 0.07 |
| 1, 2 (prior-only) | -3.296 | 0.102 (0.102) | 0.103 | 0.98 | -0.03 | 1.0 | 0.9 | 991 | -0.009 (-1.3) | -0.014 +- 0.011 | 0 within 0.04 |
| 1, 3 (Plane) | -2.984 | 0.075 (0.071) | 0.043 | 1.75 | 0.67 | 43 +- 28 | 5.1 | 23 | +0.031 (+1.0) | +0.080 +- 0.054 | 0.46, 0.37, 0.33, 0.17, 0.02 |
| 2, 0 (prior-only) | -3.297 | 0.102 (0.102) | 0.102 | 0.99 | -0.02 | 1.0 | 1.0 | 991 | -0.007 (-1.1) | -0.008 +- 0.011 | 0 within 0.06 |
| 2, 1 (Spur) | -3.167 | 0.070 (0.067) | 0.029 | 2.40 | 0.83 | 16.9 +- 9 | 10.5 | 59 | +0.027 (+1.5) | +0.065 +- 0.032 | 0.40, 0.09, 0.08, 0.02, -0.03 |
| 2, 2 (prior-only) | -3.304 | 0.096 (0.096) | 0.094 | 1.02 | 0.03 | 1.2 | 1.1 | 848 | +0.002 (+0.4) | +0.006 +- 0.012 | 0 within 0.05 |
| 2, 3 (Plane) | -3.022 | 0.074 (0.069) | 0.043 | 1.74 | 0.67 | 40 +- 28 | 5.0 | 25 | +0.027 (+0.9) | +0.092 +- 0.052 | 0.42, 0.34, 0.27, 0.17, 0.04 |
| 3, 1 (Spur) | -3.150 | 0.067 (0.066) | 0.029 | 2.34 | 0.82 | 11.1 +- 3 | 9.9 | 89 | +0.018 (+1.3) | +0.029 +- 0.025 | |
| 3, 3 (Plane) | -2.999 | 0.067 (0.066) | 0.044 | 1.51 | 0.56 | 10.4 +- 4 | 3.6 | 95 | +0.025 (+1.9) | +0.025 +- 0.024 | |
| 4, 1 (Spur) | -3.153 | 0.062 (0.062) | 0.028 | 2.19 | 0.79 | 8.7 +- 2 | 8.5 | 114 | -0.002 (-0.1) | +0.006 +- 0.020 | |
| 4, 3 (Plane) | -2.997 | 0.070 (0.067) | 0.042 | 1.64 | 0.63 | 10.4 +- 4 | 4.4 | 96 | +0.035 (+2.4) | +0.070 +- 0.025 | |

The two prior-only regions are exact prior draws by construction (sd 0.096 to 0.104, lag-1
within 0.03 of zero, trends within 1.3 sigma of zero): they check the reader and the
arithmetic, not the two-block identity. The Spur and Plane reproduce the paper (means
-3.15 / -3.17 and -2.98 / -3.02 in the two chains against the paper's -3.17 and -3.03; sds
0.069 to 0.075 against 0.06 and 0.07; the IPS tau of 17 to 43, with 40 to 65 percent errors
and biased low at this ESS, is consistent with the paper's 30 to 50). Lag-1 is 0.83 (Spur) and
0.67 (Plane) in both chains; under the two-block reading this corresponds to r 2.4 and 1.75,
consistent with the surrogate's lag-1 0.82 for the Spur and below its 0.78 for the Plane, but
the chain's lag-1 also contains the incomplete Metropolis relaxation and the other Gibbs
blocks, so this is consistency, not confirmation. The increment widths 0.029 and 0.043 are a
lower bound on the conditional sd from the Metropolis side and an upper bound from the
other-block side (surrogate 0.020 and 0.018). Beyond lag 10 the autocorrelation decays slower
than the AR(1) implied by lag-1 (0.10 to 0.17 at lag 50 in three of the four likelihood-
sampled chain-regions, about one sigma each), suggesting a slow component, which is why the
IPS tau is 2 to 8 times the AR(1) value. Both chains move towards flatter indices from their
common start, +0.065 to +0.092 per 1000 samples (1.5 to 2.0 sigma per chain-region under an
AR model with a poorly determined tau); the four tests share a start and pair within chains,
so they are not four independent detections, and in chain 1 the movement is a joint excursion
of Spur and Plane in the last quarter (quarter means Spur -3.168, -3.168, -3.159, -3.113;
Plane -3.002, -2.996, -2.996, -2.941) rather than a steady drift. Between-chain means differ
by 1.0 sigma (Spur) and 1.7 sigma (Plane). This is consistent with a chain that has not
settled and with a transient from the common start; it establishes neither. The movement adds
little to the sds (detrended 0.066 to 0.071 against 0.069 to 0.075).

Chains 3 and 4 (added 13:10, rows above; not separately reviewed) start from the same state
as chains 1 and 2, so they do not provide an independent start. They mix better: IPS tau 9 to
11 with 40 percent errors, 89 to 114 effective samples, lag-1 0.82 and 0.79 (Spur) and 0.56
and 0.63 (Plane); the prior-only regions behave as before. The four chain means agree to
0.02 (Spur -3.152, -3.167, -3.150, -3.153, pooled -3.155; Plane -2.984, -3.022, -2.999,
-2.997, pooled -3.000). All eight likelihood-sampled chain-regions move towards flatter
indices (Spur +0.073, +0.065, +0.029, +0.006 per 1000 samples, Plane +0.080, +0.092, +0.025,
+0.070; 0.3 to 2.8 sigma each under the AR model; the prior-only ones +0.012 to -0.014,
mixed signs), on average +0.04 (Spur) and +0.07 (Plane) per 1000 samples, half to one
reported sd. In the quarter means the movement is largest in the last quarters of chains 1, 2
and 4 rather than decaying from the start, which is what a slow approach to flatter values
looks like and what a transient from the common start does not; but with a shared start and
0.3 to 2.8 sigma per chain-region it remains an indication, not a measurement. The paper's
means (-3.17, -3.03) match the first-half means of chains 1 and 2.

**Resolution scaling of the collapsed target (2026-10-28, 12:25; `resolution/result_n*.json`;
first follow-up):** the same surrogate built at nside 16, 32 and 64 with LFI 30 to 70 and 353
only (WMAP exists only at nside 16), white-propagated LFI noise at all three (BP's simulated
rms exists only at nside 16; Q7 factor 1.6 applies to the widths), regions by the same
plateau rule on the finer BP map, CMB prior either flat or Gaussian with the per-nside
Planck-EE rms (0.19, 0.32, 0.44 uK); 2000-iteration chains for Spur and Plane (the JSON
`variant` field says `lfi_noise: bp` for nside 32 and 64 because the builder stores the white
propagation in the primary sigma there; all six rows are white). Cells: T2 mean +- sd with
prior / profile mean with prior / r = sd_T2 / sd_T1 / measured tau.

| nside, CMB prior | Spur pixels | High Latitude | Spur | Galactic Center | Plane |
|---|---|---|---|---|---|
| 16, white LFI noise, flat CMB | 464 | -2.175 +- 0.043 / -3.215 / r 2.83 | -2.880 +- 0.045 / -3.162 / r 2.63 / tau 18 | -3.024 +- 0.028 / -3.046 / r 2.26 | -3.157 +- 0.039 / -3.294 / r 2.42 / tau 9 |
| 32, white, flat CMB | 1863 | -1.500 (edge) / -3.212 | -2.128 +- 0.047 / -3.167 / r 6.00 / tau 77 | -3.034 +- 0.023 / -3.090 / r 2.32 | -2.875 +- 0.031 / -3.238 / r 2.86 / tau 15 |
| 64, white, flat CMB | 7434 | -1.500 (edge) / -3.236 | -1.500 (edge) / -3.234 / pinned | -3.019 +- 0.016 / -3.129 / r 2.35 | -2.651 +- 0.017 / -3.130 / r 3.35 / tau 24 |
| 16, white, CMB prior 0.19 | 464 | -3.229 +- 0.017 / -3.281 / r 1.08 | -3.166 +- 0.021 / -3.182 / r 1.07 / tau 1.2 | -3.143 +- 0.013 / -3.144 / r 1.05 | -3.286 +- 0.019 / -3.295 / r 1.11 / tau 1.5 |
| 32, white, CMB prior 0.32 | 1863 | -3.151 +- 0.013 / -3.272 / r 1.07 | -3.160 +- 0.018 / -3.203 / r 1.07 / tau 1.4 | -3.165 +- 0.010 / -3.167 / r 1.06 | -3.135 +- 0.014 / -3.152 / r 1.10 / tau 1.6 |
| 64, white, CMB prior 0.44 | 7434 | -3.065 +- 0.008 / -3.255 / r 1.06 | -3.137 +- 0.012 / -3.219 / r 1.08 / tau 1.2 | -3.098 +- 0.007 / -3.102 / r 1.04 | -2.875 +- 0.006 / -2.886 / r 1.07 / tau 1.3 |

With flat amplitude priors the collapsed target of the Spur moves from -2.88 (nside 16) to
-2.13 (nside 32) to the uniform-prior boundary (nside 64, 7434 pixels), the High Latitude
region is at the boundary from nside 32, the Plane moves from -3.16 to -2.65, and r and tau
grow with the pixel count (Spur r 2.6 to 6.0, tau 18 to 77 before pinning; Plane 2.4 to 3.4,
9 to 24); the Galactic Center (91 to 1523 pixels at high S/N) is stable at -3.02. The profile
is stable to 0.08 (Spur, High Latitude, Galactic Center) and 0.16 (Plane) across the three
resolutions. With the C_l-sized CMB prior tau is 1.2 to 1.6 everywhere and the collapsed
means move by 0.03 (Spur), 0.07 (Galactic Center), 0.16 (High Latitude; its profile is
stable at -3.28 to -3.26, so this is the synchrotron volume term) and 0.41 (Plane, -3.29 to
-2.88) over nside 16 to 64. The Plane's movement is not a volume term: its profile, which has
none, moves the same way (-3.30, -3.15, -2.89), so it is a resolution-dependent likelihood
effect this cell does not determine (candidates: pixel averaging of the model, nonlinear in
T_d, across the plane's strong gradients; region re-assignment at the plane boundary; a CMB
prior too tight for a flat-spectrum residual there; the flat-CMB profile moves 0.16 too).
This is the reviewer's toy reproduced on the real data (toy: pinned from 16 elements per
nside-16 pixel; data: pinned at 16 per pixel for the Spur). Second follow-up (12:36,
`resolution/result_n{16,32,64}_allpriors.json`): zero-mean Gaussian priors on all three
amplitudes, the CMB at the Planck-EE rms and synchrotron and dust at the sky rms of the BP
release maps at each nside (13 to 15 and 12.6 to 13.7 uK_RJ): collapsed means with prior
Spur -3.182, -3.197, -3.200; High Latitude -3.252, -3.199, -3.140; Galactic Center -3.176,
-3.222, -3.190; Plane -3.302, -3.170, -2.909 over nside 16, 32, 64; tau 1.2 to 1.6
everywhere. Spur and Galactic Center are stable to 0.05, High Latitude to 0.11, and the Plane
still moves 0.39 with every amplitude prior in place, which confirms that its movement is
not a volume term. Two candidate causes excluded (13:13, `result_n64_{reg16,td16}_allpriors.json`):
with the nside-16 region labels upgraded to nside 64 (5168 Plane pixels instead of 5133
re-assigned ones) the Plane's conditional and collapsed means stay at -2.908, and with T_d
upgraded from the nside-16 map they stay at -2.907, against -2.907 with the native nside-64
inputs; the other regions move by at most 0.03. Remaining hypothesis, untested: the Plane's
true index varies within the region (BP's own release map shows -3.0 along the outer plane
and the prior value at the centre), and a region-constant index fitted at finer resolution
weights the bright, flat-index inner plane more heavily than the diluted nside-16 pixels do;
that is a model-misspecification effect of the region-constant model, which the E020/E021
line (index fields at the resolution the data support) is designed to remove. Tested (13:30,
`nextgen/scripts/bpxv_plane_subregions.py`; all-amplitude priors, grid [-4, -2]): the Plane
split into its nside-4 pixels (25 to 39 sub-regions with at least 4 pixels) gives
sub-region indices that are resolution-stable, most within 0.05 over nside 16, 32, 64 (for
example -3.342, -3.339, -3.360; -3.232, -3.225, -3.201; -3.352, -3.345, -3.355), with a
pixel-weighted mean of -3.249, -3.252, -3.245 and a weighted spread of 0.27 to 0.30 (values
from -2.43 to -2.49 in one bright sub-region near the centre to -3.97 to -3.99, the grid
edge, in two others), while the single-index Plane fit moves -3.301, -3.170, -2.909. The
Plane's resolution dependence is therefore the misspecification of one index for a region
whose index spans 1.5, weighted differently by the resolution; the sub-regions at the grid
edge and the one at -2.4 are themselves suspect (bright plane, fixed dust SED) and are not
interpreted. Caveats for nside 32 and 64: the
noise is white and diagonal although the 60 arcmin smoothing correlates neighbours (about
0.3 at nside 64), and the regions are re-assigned per nside, so each resolution is a
different data set and the extrapolation to nside 1024 is a direction, not a value. Inference, labelled as such
(fifth review round): BP's effective number of free synchrotron elements per region is not
determined within three orders of magnitude (7k if the gauss C_l prior with D_l = 200 binds
the synchrotron above l of about 150, 44k for 30 arcmin beams, 4M raw nside-1024 pixels), so
whether the BP XV flat-prior target is shifted as the surrogate's is or neutralized by the
run's component priors is open; the Conclusion states the dichotomy and its discriminator. Whether BP's CMB carried an effective C_l prior at the index step rests on the parameter
file's first token and on the absence of `cmb/Dl_*` in the mean file; the parser semantics of
`COMP_CL_TYPE = none` are not verified, and the chains' lag-1 carries no information on it.

**Conditional width at BP's own resolution from a released sample (2026-10-28, 13:50;
`nextgen/scripts/bp_sample_cond_width.py`, Q17 item 2):** sample 000500 of chain 1 holds the
synchrotron amplitude as E and B alms with lmax 1500 (34 MB read over HTTP for l <= 1500).
The information term of the conditional, sum over bands, Stokes and region pixels of
w (a_s f_nu ln(nu/30))^2, with a_s the sample's Q/U amplitude put through each LFI beam
(32.29, 27.0, 13.25 arcmin) and pixel window at the band's native nside (512, 512, 1024)
against the release Q_RMS/U_RMS maps, plus WMAP Ka, Q, V at nside 16 with the diagonal of
their covariances, and the N(-3.3, 0.1) prior:

| Region | sd_cond, sample amplitudes | sd_cond, release mean map (60 arcmin, no pixel noise) | ratio |
|---|---|---|---|
| High Latitude | 0.0088 | 0.0134 | 1.5 |
| Spur | 0.0160 | 0.0186 | 1.2 |
| Galactic Center | 0.0091 | 0.0109 | 1.2 |
| Plane | 0.0085 | 0.0116 | 1.4 |

The sample's amplitude Q rms in the Spur is 14.7 to 15.7 uK_RJ through the three beams
against 13.2 for the mean map: the sample carries little small-scale amplitude noise, so the
Gaussian C_l prior (`COMP_CL_TYPE02 = gauss`) binds the synchrotron amplitude and BP's
effective element count is at the low end of the 7k to 4M range (the surrogate's curvature
scaling 2150 + 0.4 N gives N of about 1700 Stokes-elements for the Spur's 0.016, less than
the nside-16 surrogate's 928 times two). Consequences, from the numbers alone: (1) the
chain's per-step increments (0.029 Spur, 0.043 Plane) are 1.8 and 5 times the conditional
width, so the amplitude block contributes about 30 percent (Spur) and 4 percent (Plane) of
the increment variance and the other Gibbs blocks the rest; (2) with the chain sds as the
marginal width, r = 0.069 / 0.016 = 4.3 (Spur) and 0.074 / 0.0085 = 8.7 (Plane), Gaussian
two-block tau 36 and 150, against the measured IPS 17 to 43 (biased low at 23 to 114
effective samples); (3) a chain with r of 4 to 9 relaxes towards its target in of order r^2
= 20 to 80 steps, so the 1000-sample chains, which moved by 0.05, are at a target within
about 0.05 of the reported values: the flat-prior shift the surrogate shows (+0.27 at nside
16, the boundary at nside 64) is not present in BP's chain, which means the CMB column does
not act as a free per-pixel amplitude there (an effective constraint from `COMP_CL_TYPE01 =
none` or the CG's treatment of the CMB; the semantics remain unverified). The measurement
uses the information term only (the residual term of the curvature is not included), delta
bandpasses, and the release rms maps; it is an order-of-magnitude discriminator.

Surrogate consistency: the conditional at the BP release posterior-mean amplitudes pins at
the grid edge -1.5 in every region (chi2 prefers the flattest synchrotron allowed), so the
release amplitudes do not fit the surrogate at the level a conditional of width 0.02 resolves;
expected (bandpass integration, native resolution, full-mission 353, the BP 10-degree-smoothed
index map, unsmoothed WMAP r4 maps against 60-arcmin amplitudes) and the reason T1 is
evaluated at draws. The surrogate's T2 means are not comparable to BP XV's -3.17 +- 0.06 and
-3.03 +- 0.07 for the same reasons and because the flat-prior volume term scales with the
number of amplitude pixels (Follow-ups); only widths and ratios are used.

## Negative control
The intensity analogue of F14: the same scripts on the E014 nside-64 intensity configuration are
not rerun here. The chain check is the exact identity lag-1 = 1 - E[var_cond]/var(beta), which
holds for any two-block Gibbs sampler regardless of Gaussianity; the Gaussian formula
2 r^2 - 1 is only a reference. The prior-only regions of the released BP chain are exact
prior draws by construction; they check the reader and the arithmetic (sd 0.10, no
autocorrelation, no trend), not the two-block identity (Result).

## Conclusion
(stats review 2026-10-28, three rounds: targets, code, inputs and the chain numbers
CONFIRMED; the conclusion wording on the chains follows the reviewer's must-fix list)
- On an nside-16 surrogate of the BP XV polarization configuration with the released
  configuration's flat per-pixel CMB prior, the collapsed posterior sd of beta_s is 2.1 to
  2.35 times the conditional sd for the two likelihood-sampled regions and a two-block Gibbs
  sampler with exact conditional draws has integrated autocorrelation times of 7 to 11 (lag-1
  equals the exact identity 1 - E[var_cond]/var), for either LFI noise model; a C_l-sized
  Gaussian CMB prior removes the coupling (ratio 1.02 to 1.04). The coupling is carried by the
  free CMB amplitude, not by the synchrotron amplitude (per-pixel S/N 2 to 3), which is why
  F14's ratio of 180 (intensity, S/N in the hundreds) does not appear. The surrogate's volume
  term and stationary distribution scale with the number of amplitude elements and are lower
  bounds for BP's native resolution.
- The released BP chains (four, 1001 samples each, one shared start) have lag-1 0.79 to
  0.83 (Spur) and 0.56 to 0.67 (Plane), IPS tau 9 to 43 consistent with the paper's 30 to
  50, and 23 to 114 effective samples. The conditional width at BP's own resolution,
  measured from a released sample's synchrotron amplitudes through the LFI beams, is 0.016
  (Spur) and 0.0085 (Plane): the Gaussian C_l prior binds the synchrotron amplitude, so the
  effective element count is at the low end and the amplitude block contributes about 30
  percent (Spur) and 4 percent (Plane) of the chain's per-step increment variance, the rest
  coming from the other Gibbs blocks. With the chain sds as marginal widths the two-block
  r is 4.3 (Spur) and 8.7 (Plane), Gaussian tau 36 and 150; the measured IPS tau of 17 to
  43 is biased low at this effective sample size, so the Plane's coupling may be
  underestimated by the paper's 30 to 50 while the Spur's is about right. A chain with r of
  4 to 9 relaxes in 20 to 80 steps, so the released chains, which moved 0.05 in 1000
  samples, sit within about 0.05 of their target; the flat-prior shift of the surrogate
  (+0.27 at nside 16, the boundary at nside 64) is not present in BP's chain, which means
  the CMB column does not act there as a free per-pixel amplitude (the mechanism in
  Commander3's handling of `COMP_CL_TYPE = none` is unverified). The reported sds of 0.06
  and 0.07 are then posterior widths with 23 to 114 effective samples per chain (10 to 15
  percent Monte Carlo uncertainty), and the slow movement is consistent with other-block
  equilibration. Wording subject to the sixth review round. How to present this in Paper A
  is Trygve's call (STATUS).
- The flat-amplitude-prior marginal is shifted from the profile by +0.27 (Spur) and +0.13
  (Plane) on the surrogate and runs to the flattest allowed index in the High Latitude region;
  a C_l-sized CMB prior removes the runaway and moves the Spur to -3.20. This is the
  prior-volume effect of Paper A Sect. 4 on real polarization data, carried by the free
  per-pixel CMB amplitude at S/N below 0.2, and it grows with resolution (measured: Spur
  -2.88, -2.13, boundary at nside 16, 32, 64). With a C_l-sized CMB prior the coupling is
  gone (tau 1.2 to 1.6) and the collapsed means move by 0.03 (Spur), 0.07 (Galactic Center),
  0.16 (High Latitude; profile stable, the synchrotron volume term) and 0.41 (Plane; the
  profile moves too, a resolution-dependent likelihood effect this cell does not explain)
  over nside 16 to 64. With zero-mean data-scaled priors on all three amplitudes the Spur
  and Galactic Center are stable to 0.05 and High Latitude to 0.11, so priors on every
  component with a data-set scale (the E016 hierarchical form is the proper version) give a
  resolution-stable target where the volume term is the cause; the Plane still moves 0.39
  with all priors in place, which the sub-region test attributes to one index being fitted
  to a region whose index spans 1.5 (sub-region values resolution-stable, their weighted
  mean -3.25 at all three resolutions, the single fit from -3.30 to -2.91): a region-constant
  index is misspecified in the plane, the case for index fields at the resolution the data
  support (M2.x).
- Limits: diagonal noise in the collapsed target (the dense WMAP covariance changes the
  conditional width by under 5 percent); delta bandpasses; NPIPE A split at 353; unsmoothed
  WMAP r4 maps; regions recovered from the smoothed release map; nside 16 instead of
  multi-resolution; the movement in the BP chains is a 0.3 to 2.8 sigma effect per
  chain-region, all eight from one shared start state (the release has no chain started
  elsewhere); the nside 32 and 64 surrogates use white diagonal noise (neighbour correlation
  from the smoothing ignored) and re-assigned regions; the within-conditional Metropolis
  relaxation of BP's polarized alm sampler is not characterized here (the stored
  `beta_pixreg_proplen` and `nprop` of 6e-4 and 1000 belong to the intensity row).

## Follow-ups
- The two Plane sub-regions at the grid edge (-4.0) and the one at -2.4: inspect for
  compact sources and the fixed dust SED before any science use; E022 does not interpret
  them.
- A longer or differently started polarized beta_s chain is the only way to settle the
  movement question; it is not in the release.
- Dense WMAP covariance in the collapsed target (18432-amplitude Gaussian marginal per beta
  grid point) if the diagonal approximation moves the conditional width by more than 20 percent.
- Bandpass-integrated mixing matrix from the BP bandpass files instead of delta bandpasses.
