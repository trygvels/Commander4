# Pipeline assumptions audit: the Commander4 Gibbs sampler as coded

| | |
|---|---|
| Written | finished 2026-11-09 22:30 (the file's modification time; the agent left a placeholder here, filled in by the main session) |
| Repository state | branch `nextgen` at `52061ad` |
| Author | Claude Code subagent, at the request of the main session. Not independently reviewed. |
| Method | Reading and grepping only. Nothing was executed (a production MPI chain holds all cores). |
| Evidence level | Code facts were read at the cited lines. Every magnitude is a hand estimate unless it cites an experiment ID. Every "possible bug" is a hypothesis to be verified. |

Scope read: `src/commander4/cli.py`, `compsep/`, `sky/`, `data_models/`, `file_io/`,
`math_utils/`, `mpi/`, `tod/` (briefly), `src/lib_cpp/compsep_pymod.cc`,
`params/param_explanation.yml` and `params/Planck_mapsonly/param_LFI+HFI+WMAP.yml` (other
parameter files were only grepped).

Cross-checked against `nextgen/FINDINGS.md` (F1 to F38 and sections 1 to 7),
`nextgen/OPEN_QUESTIONS.md` (Q1 to Q18), `nextgen/DECISIONS.md` (D001 to D027),
`nextgen/ROADMAP.md`, `nextgen/DATA.md` and the experiment records E001 to E031.

Path convention: source paths are relative to `src/commander4/` unless they begin with `src/`,
`params/` or `nextgen/`. Line numbers refer to commit `52061ad`.

Sibling document: `nextgen/reviews/2026-11-09_referee_report.md` reviews the findings. This
document audits the code and was written without reading it.

## 0. How to read this

Every entry has the same fields.

- **Assumes**: what the code does, as coded, not as papers describe it.
- **Where**: file and line.
- **Kind**: `exact` (an exact statement of a stated model), `approximation` (a known
  simplification of the stated model), or `silent convention` (a choice that changes answers and
  is not presented to the user as a model choice).
- **Could bias**: CMB map, foreground amplitudes, spectral indices, uncertainties, convergence.
- **Test or relax**: the cheapest check I can think of, and what relaxing it would take.
- **Status**: (a) already tested by a recorded experiment, with the F number; (b) known in
  `nextgen/` but untested; (c) not mentioned anywhere in `nextgen/`.

Definition of (c): a keyword search over `nextgen/**/*.md`, excluding `reviews/` and `research/`,
found no statement of the assumption. Where a passing mention exists (the thing is named but its
consequence is not discussed) I say so and still list it under (c). An assumption discussed under
words I did not search for could have been missed, so (c) means "I found no record", not "there
is no record".

Item labels: G = Gibbs structure, D = data model, S = sky model, A = amplitude solves, T = TOD
side, B = possible bug or inconsistency.

## 1. Headline results

### 1.1 The validated region and the real-data region barely overlap

As far as the experiment records show, every simulation behind F12 to F38 is intensity only, has
white noise, is generated at the analysis nside with no beam and no smoothing, and (E031) sits at
24 to 90 times the white noise of real WMAP and LFI intensity maps. The one real-data analysis
(E014) is intensity only, nside 64, 120 arcmin, per-pixel solver, `ridge`.

It follows that the following code paths carry real-data results or will carry them, and have no
recorded test against a known truth:

1. the map ingest path: degradation, noise-weighted smoothing, the model-error floor (D2, D3, D5,
   D9);
2. everything polarized: QU smoothing, QU noise, the shared I and QU index (D1, D4, D6, S2);
3. the CG amplitude path: truncated solves, the null prior, the preconditioner (A1 to A5);
4. bandpass and unit conversion (D7) and per-band calibration (D10).

This is a statement about coverage, not about correctness. The samplers validated in F16 are
validated for the model they are given. What is unvalidated is whether the maps handed to them
obey that model.

### 1.2 Possible bugs that could move real-data numbers (all to be verified)

| ID | One line | Where |
|---|---|---|
| B1 | Maps are degraded with `hp.ud_grade` before they are smoothed, so sub-pixel power is aliased in a band-dependent way and an unrecorded pixel window is added | `file_io/map_reader.py:156-162` before `:173-175` |
| B14 | The model-error floor is added at native resolution and then averages down like white noise, about 26 times at E014 settings | `file_io/map_reader.py:141-143` before `:156-162`, `:173-175` |
| B3 | `get_sky` evaluates a component with an index field at the mean index; the chi-square diagnostic, the fixed-component residual, the conditional MH likelihood and the TOD-side sky all go through it | `sky/diffuse_components.py:442-450` |
| B5 | "No C(l) prior" is a unit-variance Gaussian prior per alm in (uK_RJ at `nu_ref`)^2, so the answer depends on the unit and on `nu_ref` | `sky/diffuse_components.py:235-236`, `compsep/cg_solver.py:191-193, 262-268` |
| B4 | The per-pixel amplitude draw has no guard for pixels observed by fewer bands than components, while the likelihood kernels skip them | `src/lib_cpp/compsep_pymod.cc:119-121` against `:272`, `:405` |
| B2 | Q and U are smoothed as two independent scalar maps | `data_models/detector_map.py:171-175` |
| B15 | The conditional MH likelihood sums chi-square in single precision | `compsep/mcmc.py:158-160` |

Details and the reasons each might turn out harmless are in section 7.

### 1.3 Steps accepted as samples that are not exact conditional draws

The CG amplitude solve is run for a fixed 100 iterations from a zero start and accepted whatever
the residual (A1; the documentation says this is what every real-data file does). On the TOD
side the white-noise level is a point estimate (T2), the default correlated-noise draw is the
stationary solution on gap-filled data (T3), and the temporal-gain CG has no check after its
200 iteration cap (T1). NUTS warm-up draws are written to the chain like any other (G7). None of
these is recorded in `nextgen/` as a validity question; Q12 treats the CG as a cost question
only.

### 1.4 Where the rest is

The (c) list is section 9. The two rankings are section 10. Ten cheap tests, none of which needs
the cluster, are section 11. The limits of this audit are section 12.

## 2. Gibbs structure and validity

### G1. Update order and two-chain pipelining

- **Assumes**: one Gibbs iteration per chain is: TOD side (glitches, absolute gain, relative
  gain, temporal gain, correlated noise and PSD, sigma0, mapmaking), then CompSep side (all
  amplitude groups in file order, then all MCMC groups). Under `marginal`, `ridge` and `hier` the
  amplitudes are redrawn from their conditional after the index steps. Two chains always run,
  interleaved per iteration (`for iteration ... for chain in (1, 2)`), with TOD and CompSep on
  disjoint MPI ranks, so chain 1's CompSep runs while chain 2's TOD runs. Each chain's TOD step
  conditions on that chain's own latest CompSep sample.
- **Where**: `cli.py:41-49`, `cli.py:66-147`, `compsep/processing.py:822`,
  `tod/processing.py:176-190`.
- **Kind**: exact. A systematic-scan blocked Gibbs sweep is valid if every block is a valid
  conditional draw or a valid MH step. G6 lists the blocks that are not.
- **Could bias**: nothing by itself. Convergence of the joint chain depends on the correlation
  between TOD parameters and the sky (Q10).
- **Test or relax**: an end-to-end run with known truth at the upper levels of
  `nextgen/VALIDATION.md`.
- **Status**: (a) F9, established by code reading only. No end-to-end run is recorded.

### G2. The legacy `resolve_chisq` move and its latent default

- **Assumes**: the legacy move accepts an index proposal on the chi-square after re-solving the
  amplitudes. It is still present as `_run_resolve_chisq` and selectable. It is also the default
  of a dictionary getter: `getattr(compsep_state, "mcmc_likelihoods", {}).get(group.name,
  "resolve_chisq")`.
- **Where**: `compsep/mcmc.py:260`, `compsep/processing.py:878-879`.
- **Kind**: silent convention (a latent fallback to a sampler known to be invalid).
- **Could bias**: spectral indices and their uncertainties, if the fallback were ever reached.
- **Test or relax**: make the getter raise on a missing key. Presumably unreachable today, because
  `_resolve_mcmc_likelihoods` (`compsep/processing.py:434`) fills the dictionary for every group;
  this was not verified by running.
- **Status**: invalidity of the move (a) F1, F2. The latent fallback (c).

### G3. Amplitude prior behind each likelihood mode, and the `auto` rule

- **Assumes**: `auto` resolves to `ridge` where the per-pixel marginal is available and to
  `conditional` otherwise. `marginal` integrates the amplitudes under a flat prior. `ridge`
  accepts on exp(-chi2_min/2), which equals the flat-prior marginal times
  det(M^T N^-1 M)^(1/2) per pixel, so it is the exact marginal for an amplitude measure that
  depends on the index. `conditional` holds the amplitudes fixed.
- **Where**: `compsep/processing.py:434-470`, `compsep/perpix_solver.py:155-237`,
  `src/lib_cpp/compsep_pymod.cc` (`marginal_loglike_perpix`).
- **Kind**: exact for a stated prior. The prior is the assumption.
- **Could bias**: spectral indices through the prior volume term; the CMB through any finite
  amplitude prior.
- **Test or relax**: done at length.
- **Status**: (a) F3, F4, F12, F15, F16, F26; D003.

### G4. The parameter mask enters the index step only

- **Assumes**: `chisq_mask` is applied in the MCMC likelihood and nowhere else. One mask serves Q
  and U, and only field 0 of the mask file is read. The amplitude solve and the noise model use
  all pixels. The documentation states this: "This is a parameter-domain mask only: it does not
  enter the amplitude solve or the noise model."
- **Where**: `compsep/processing.py:233-258`, `compsep/perpix_solver.py:201-202`,
  `compsep/mcmc.py:157`, `params/param_explanation.yml:456-462`.
- **Kind**: silent convention for the statistics. In the per-pixel path without smoothing this is
  a valid "cut" posterior: masked pixels carry free amplitudes and no information on the index.
  In the CG path, and after common-resolution smoothing, pixels are coupled, so the index step
  and the amplitude step are not conditionals of one joint distribution.
- **Could bias**: spectral indices (plane signal leaks into kept pixels through the smoothing and
  through the CG), and uncertainties.
- **Test or relax**: simulation with a bright mis-modelled plane; compare the index posterior
  with the mask applied before and after smoothing. Relax by masking in the data model (zero
  inverse variance) when a group asks for it.
- **Status**: (c) for Commander4. The same structure is noted for Commander3 in
  `nextgen/experiments/E009_c3_nonlin_review/independent_reviews_2026-10-25.md:291`.

### G5. The index proposal

- **Assumes**: a joint Gaussian random walk over all sampled scalars with a fixed, required
  per-component sigma, `numstep` 1 by default, no adaptation.
- **Where**: `compsep/spectral_index.py:71-80, 126-133`, `params/param_explanation.yml:419`.
- **Kind**: exact (a symmetric proposal), with efficiency left to the user.
- **Could bias**: convergence. E014 needed 20 steps per iteration and sigmas of 0.004 and 0.002
  by hand to reach an acceptance of 0.24.
- **Test or relax**: adapt the proposal covariance during a declared warm-up, or use the NUTS path.
- **Status**: (a) F5, F14.

### G6. Blocks that are not draws from an exact conditional

Collected here so the list exists in one place. Each has its own entry.

| Block | What it is instead | Entry |
|---|---|---|
| CG amplitude solve | Fixed iteration count from a zero start, accepted at the cap | A1 |
| CG amplitude solve with `optimize: true` | The Wiener mean, on request | A1 |
| Per-pixel amplitude draw | Exact per pixel, then stored through a band-limited alm round trip | D14 |
| Per-pixel amplitude draw, deficient pixels | Zero or ill-conditioned amplitudes | A7 |
| sigma0 | Deterministic estimator | T2 |
| Correlated noise | Stationary Wiener solution on white-noise-filled data (default) | T3 |
| Noise PSD parameters | Gridded Whittle sampler with 4 internal passes | T4 |
| Temporal gain | CG with a 200 iteration cap and no check after the loop; mean removed after the draw | T1 |
| CG mapmaker | Accepted at the iteration cap | T7 |
| NUTS index field during warm-up | Adaptive, not a fixed kernel | G7 |
| Conditional MH likelihood | Evaluated at the mean index for a component with an index field | S4 |

### G7. Initialization, restart, burn-in

- **Assumes**: the initial sky is `init_from` if given, else zero amplitudes. Both chains start
  from the same state, also after a restart ("A chain the source run never wrote starts from the
  first chain it did write."). The sampler has no notion of burn-in; only a standalone tool has a
  burn-in option. NUTS warm-up runs inside the recorded chain (`adapting = self.iteration <=
  self.n_warmup`, and warm-up draws are returned like any other). The NUTS adaptation state lives
  in a module-level dictionary and is not written to the chain file, so a restart adapts again.
- **Where**: `sky/sky_model.py:7-16`, `file_io/chain_restart.py:136-139`,
  `standalone_tools/constrained_cmb_realizations.py:553`, `math_utils/nuts.py:214-220`,
  `compsep/index_field_sampler.py:427, 544-554`.
- **Kind**: silent convention.
- **Could bias**: uncertainties and any convergence diagnostic. Two chains from one starting
  point cannot show a between-chain disagreement that comes from the start. Warm-up draws in the
  output bias summaries unless the analysis drops them, and after a restart the number to drop is
  not recorded in the file.
- **Test or relax**: write an `is_warmup` flag per sample and the adaptation state to the chain
  file; start chain 2 from a dispersed state.
- **Status**: warm-up and metric handling (b) D012, D015. The shared start and the absence of
  burn-in bookkeeping (c).

### G8. Random number streams

- **Assumes**: one generator per (seed, chain, iteration, world rank):
  `np.random.default_rng([root_seed, chain, iteration, mpi_info.world.rank])`. The default
  `root_seed` is 1995 and the example files set `seed: 1995`.
- **Where**: `cli.py:52-63`, `file_io/chain_restart.py:1-19, 110-134`. Legacy global `np.random`
  calls remain only in dead or standalone code: `compsep/constrained_cmb_loop.py:118`,
  `compsep/constrained_cmb_loop_mpi.py:119`,
  `standalone_tools/constrained_cmb_realizations.py:353`,
  `simulations/inplace_litebird_sim.py:61, 255`.
- **Kind**: exact within one run (the streams are independent). Silent convention across runs.
- **Could bias**: uncertainties when runs are combined or differenced. Two runs with the same
  seed and MPI layout draw from the same streams, and where the sequence of draws is the same
  they share every fluctuation vector. Their difference then understates the sampling scatter
  and their union is not two independent chains. All 19 parameter files under
  `nextgen/params/q1/` set `seed: 1995`. A run started with `gibbs.init_chain` numbers its
  samples from iteration 1 (`file_io/chain_restart.py:9-10, 134`), so with an unchanged seed it
  draws the streams of iterations 1 onward of the source run a second time. A chain extended
  that way is driven by the same random numbers as the stretch it extends, and for a contracting
  sampler (a sweep of Gaussian conditional draws is one) the new stretch converges onto the old
  one, which adds samples without adding information (B16; from reading, not from a run).
  `chain_status: append` continues the numbering and is not affected. Streams also depend on the
  MPI layout, which the documentation states (`params/param_explanation.yml:31-34`).
- **Test or relax**: add the run identifier or the absolute sample index to the seed sequence;
  record the seed in the chain file.
- **Status**: (c).

### G9. Only chain 1 is written

- **Assumes**: `output.chains.write` is `[1,]` in the documentation and examples, so chain 2 is
  computed and discarded.
- **Where**: `file_io/chain_writer.py:195`, `params/param_explanation.yml:92-93`.
- **Kind**: silent convention.
- **Could bias**: convergence assessment. The second chain costs half the compute and is the
  natural input to a between-chain diagnostic that nothing computes.
- **Test or relax**: write both chains and report a split or between-chain R-hat per parameter.
- **Status**: existence and cost (b) E020, E021. No between-chain diagnostic and the shared start
  (c).

## 3. Data model

### D1. Independent Gaussian noise per pixel and per Stokes parameter

- **Assumes**: the noise covariance is diagonal in pixel space and between Q and U. "Pixel-pixel
  and QU covariance matrices are not supported." The TOD side computes the per-pixel IQU
  covariance and passes on only `1/diag(A^-1)`.
- **Where**: `params/param_explanation.yml:501`, `file_io/map_reader.py:144-148`,
  `compsep/mcmc.py:159-160`, `src/lib_cpp/compsep_pymod.cc:103-118`,
  `compsep/cg_solver.py:225-230`, `tod/mapmaking/output.py:37-39`, `tod/mapmaking/binned.py:28`.
- **Kind**: approximation.
- **Could bias**: uncertainties everywhere, and the weighting of large-scale polarization most.
  WMAP and LFI low-multipole polarization noise is dominated by correlated modes that a diagonal
  model cannot represent. Dropping the QU correlation keeps each Stokes variance right but
  miscalibrates any joint statement about polarization amplitude or angle.
- **Test or relax**: at nside 16, compare amplitude draws with a dense noise covariance against
  the diagonal model on simulations with the released WMAP and LFI covariance structure. Relax in
  two steps: a 2 by 2 QU block per pixel in the per-pixel kernel (cheap), then a low-resolution
  dense covariance.
- **Status**: QU covariance (c). The level of correlated noise (b): Q7, E022 (variance 2.5 to 2.8
  times white for BP LFI at nside 16 and 60 arcmin), limits of E031.

### D2. Smoothed noise treated as independent between pixels

- **Assumes**: after common-resolution smoothing, the per-pixel variance is propagated
  analytically and the pixel-pixel covariance created by the smoothing is dropped.
- **Where**: `data_models/detector_map.py:42-69`.
- **Kind**: approximation.
- **Could bias**: uncertainties of every global or region-level parameter, in the overconfident
  direction. The number of pixels per independent noise patch is 4 pi sigma^2 / omega_pix, which
  is 10.8 at nside 64 and 120 arcmin (8.7 for WMAP K, whose extra kernel is 107.8 arcmin). A
  noise-limited global parameter would then have a posterior too narrow by up to sqrt(10.8) = 3.3.
  This factor is my estimate and has not been measured for the E014 setup.
- **Test or relax**: Q7 as written. Cheapest version: white noise at native resolution through
  the ingest path, then the scalar-index posterior width against the scatter over realizations.
  Relax by working at a pixel size matched to the beam, by a harmonic-space likelihood for the
  index step, or by an effective-sample-size correction that is validated rather than assumed.
- **Status**: (a) F8 (E006: 1.9 and 3.7 standard deviations). Q7 is still open.

### D3. Noise-weighted smoothing gives each band its own, position-dependent beam

- **Assumes**: the signal is smoothed as `smooth(d * w) / smooth(w)` with `w` the inverse
  variance. The effective kernel is B_ij w_j / sum_j B_ij w_j, which differs from the Gaussian
  wherever `w` varies across the beam, and differs between bands because `w` does. The per-pixel
  model that follows assumes identical beams, and the documentation says the signal and the rms
  go "through the same kernel".
- **Where**: `data_models/detector_map.py:22-39`, `params/param_explanation.yml:322-326`.
- **Kind**: silent convention.
- **Could bias**: foreground amplitudes and indices near strong weight gradients: scan-pattern
  features, survey edges, and above all any band whose rms tracks the signal. The beam centroid
  moves by about sigma^2 grad(ln w). For a weight proportional to 1/signal^2 the smoothed value
  is low by about twice the fractional variance of the signal inside the beam. My estimates:
  0.1 percent at high latitude, percent level in the plane. Small against D5 and D7, but it is a
  band-dependent multiplicative error, which is what an index fit is sensitive to.
- **Test or relax**: smooth one simulated sky with each band's real weights and difference the
  results. Relax with plain smoothing of the signal and a separate, correct propagation of the
  variance.
- **Status**: (c) for the consequence. Its existence is mentioned in passing in
  `nextgen/experiments/E014_q1_real_maps_lowres.md:145-147`.

### D4. Q and U are smoothed as independent scalar maps

- **Assumes**: the smoothing loop runs over Stokes maps and calls the scalar smoother on each.
- **Where**: `data_models/detector_map.py:171-175`.
- **Kind**: approximation, exact only in the flat-sky limit. Possible bug for polarization (B2).
- **Could bias**: the polarized CMB and polarized foreground amplitudes near the poles of the
  map's coordinate system, where the Q and U basis rotates across the beam. For a locally uniform
  field the polarized amplitude is reduced by about 2 sigma^2 cot^2(theta), with theta the
  colatitude: 0.13 percent at 30 degrees from a pole and 1.4 percent at 10 degrees for a 120
  arcmin kernel; 3.3 percent at 30 degrees from a pole for a 10 degree kernel. There is E to B
  mixing of the same order. For Galactic-coordinate maps the poles are the cleanest sky.
- **Test or relax**: smooth an E-only simulated QU map both ways and look at the B power and at
  the difference against latitude (seconds with healpy). Relax by a spin-2 transform, which the
  code already uses elsewhere.
- **Status**: (c).

### D5. Degradation before smoothing

- **Assumes**: a map whose nside differs from `eval_nside` is first brought to `eval_nside` with
  `hp.ud_grade` (unweighted mean of the signal; inverse variances summed with `power=-2`), and
  only then smoothed to the common resolution.
- **Where**: `file_io/map_reader.py:156-162`, then `:173-175`. Confirmed for the E014 inputs:
  `nextgen/scripts/prep_pr4_maps.py` writes nside 512 and the parameter files set
  `eval_nside: 64`.
- **Kind**: silent convention. Possible bug (B1).
- **Could bias**: everything downstream of the maps, in three ways.
  1. Aliasing. A pixel average is a weak low-pass filter, so power above the nside-64 band limit
     folds into the map before the smoothing can remove it. For the CMB alone I estimate several
     hundred uK^2 of folded variance for arcminute beams and about 100 uK^2 for WMAP K, of which
     the smoothing leaves an rms of roughly 5 to 10 uK and 1 to 3 uK. The estimate is uncertain
     by a factor of a few. The part common to all bands has the CMB spectrum and is absorbed by
     the CMB amplitude. The part that differs between bands, because it scales with the native
     beam, is a few uK, which is the size of the per-pixel rms in the low-frequency bands.
     Compact sources are worse: a narrow-beam band puts a source in one coarse pixel, a
     wide-beam band spreads it according to its true position.
  2. An unrecorded window. The nside-64 pixel average is equivalent to roughly 37 arcmin of extra
     FWHM, so the maps are at about 126 arcmin, not the 120 recorded in `amp_fwhm_rad`. Common to
     all bands, so harmless for indices, but wrong for anything that uses the stated beam.
  3. Variance. The variance of an unweighted mean is sum(s_i^2)/n^2, and the code assigns
     1/sum(1/s_i^2), which is never larger. Equal only for uniform noise inside the coarse pixel.
- **Why I think it may matter**: the E014 high-latitude reduced chi-squares are largest in the
  widest-beam bands (Haslam 6.4, WMAP K 4.5, Planck 30 GHz 3.8, against 1.1 at 70 GHz). That is
  what band-dependent aliasing would do. It is also what low-frequency foreground model error
  would do, and I cannot separate the two by reading.
- **Test or relax**: a CMB-only simulation at nside 512 with two native beams through both
  orders, differenced (a minute with healpy). Relax by smoothing at native resolution and then
  degrading, or by degrading in harmonic space.
- **Status**: (c).

### D6. Beams are symmetric Gaussians, with no pixel window and a spin-0 window for polarization

- **Assumes**: every beam is `hp.gauss_beam` at the stated FWHM. `blm` raises
  `NotImplementedError`. The forward model has no pixel window ("CompSep models the data as
  B*Y*a with this beam and NO pixel window."). `apply_B` multiplies every alm row, E and B
  included, by the spin-0 window, while `_realize_alms_as_map` calls `hp.smoothalm` with its
  default `pol=True`, which gives row 0 the spin-0 window and later rows the spin-2 window. For a
  two-row QU view that means E gets one window and B the other.
- **Where**: `data_models/detector_map.py:127, 179, 213-219, 221-226`,
  `compsep/preconditioners.py:251-255`, `sky/diffuse_components.py:429`,
  `params/param_explanation.yml:510-511`; `hp.pixwin` appears only at
  `data_models/detector_map.py:58`.
- **Kind**: approximation (Gaussian beams, no window) plus a small confirmed inconsistency (B7).
- **Could bias**: the CMB spectrum at high multipoles through beam shape and the missing pixel
  window; foreground amplitudes near bright sources through the non-Gaussian beam wings. The two
  polarization windows differ by exp(2 sigma^2) - 1: 1.1e-4 at 1 degree, 2.7e-3 at 5 degrees,
  1.1e-2 at 10 degrees. Negligible except in a low-resolution polarization analysis.
- **Test or relax**: read the released beam transfer functions; multiply by `hp.pixwin` when the
  input maps are pixel averages; use one polarization window convention everywhere.
- **Status**: (c). Effective beams appear in `nextgen/DATA.md` as data products only.

### D7. No bandpass integration, no colour correction, monochromatic unit conversion

- **Assumes**: every SED is evaluated at one frequency per band, and input units are converted
  at that frequency. `SkyModel.get_sky(band)` and `DetectorMap.bandShape` raise
  `NotImplementedError`. The LFI reader carries "TODO: Re-implement the per-detector bandpass
  shift".
- **Where**: `sky/sky_model.py:25-28`, `units.py:3-5`, `file_io/map_reader.py:104-108`,
  `data_models/detector_map.py:233-236`, `file_io/experiments/planck_lfi.py:12`,
  `params/param_explanation.yml:512-513`.
- **Kind**: approximation, documented upstream.
- **Could bias**: dust amplitude and index most. E014 records a unit conversion at 857 GHz that
  is off by a factor 1.58 (monochromatic 1.44 against band-integrated 2.27 MJy/sr per K_CMB) and
  2.5 percent at 353 GHz. The closing remark of F13 quotes 7 to 14 percent, not quantified
  further.
- **Test or relax**: integrate the SED over the released bandpass per band (a table lookup per
  index value); on simulations, generate with bandpasses and fit with delta functions.
- **Status**: (b) closing remark of F13, E014, section 5 of FINDINGS, `nextgen/DATA.md:369`.
  Bandpass mismatch is on `nextgen/ROADMAP.md:56-58`.

### D8. Missing header information is filled in silently

- **Assumes**: a map without a unit keyword is taken as uK_CMB and a map without an ordering
  keyword as RING, with a warning only.
- **Where**: `file_io/map_reader.py:59-61, 83-85`.
- **Kind**: silent convention.
- **Could bias**: everything, for one band, if a file is in another unit or in NESTED order. A
  NESTED map read as RING is scrambled, which shows up at once. A unit error is a
  calibration-like factor that can pass unnoticed for a foreground-dominated band.
- **Test or relax**: make both required, or abort unless the parameter file states them.
- **Status**: (c).

### D9. The model-error floor `add_signal_fraction_to_rms`

- **Assumes**: the rms map gains, in quadrature, the stated fraction times the full-sky mean
  absolute value of the same map: `map_rms = np.sqrt(map_rms**2 + (fraction *
  np.nanmean(np.abs(map_sky)))**2)`. It is applied at native resolution, before the degradation
  and the smoothing. The documentation calls it "a stand-in for the regularization noise C3 adds
  at the TOD-to-map handoff".
- **Where**: `file_io/map_reader.py:141-143`, `params/param_explanation.yml:504-509`.
- **Kind**: silent convention.
- **Could bias**: the relative weight of bands, hence indices, and all uncertainties.
  1. The scale is a full-sky mean of an absolute value, so it is set by the Galactic plane and by
     whatever zero level the map carries, not by the signal where the fit happens.
  2. It is white. From nside 512 to 64 the inverse variances of 64 pixels are summed (rms down by
     8), and the smoothing takes about another 3.3 (D2), so the floor that reaches the likelihood
     is about 26 times below the nominal one. A nominal 0.05 becomes roughly 0.2 percent of the
     mean absolute map per analysis pixel. Real model error (calibration, bandpass, zero level,
     SED mismatch) is coherent and does not average down.
  3. The E014 high-latitude table lists "rms (with floor)" of 13.4 uK_RJ at 353 GHz and 23 at 857
     GHz, with reduced chi-squares of 0.0. The record says the table comes from a separate
     reconstruction. If those values are the nominal floor, the chain itself saw a floor about 26
     times smaller, and those two bands were far from down-weighted. This is an inference from
     the order of the code lines and must be checked against the inverse variances the chain
     actually used (B14).
- **Test or relax**: one call to `read_data_map_from_file` per band, as
  `nextgen/scripts/q1_grid_check.py:42` already does, printing percentiles of the rms at high
  latitude. Relax by modelling the errors that the floor stands in for: a gain factor per band
  (D10), an offset per band (S3), and a floor defined at the analysis resolution.
- **Status**: its use (b) E014. The three properties above (c).

### D10. No gain or calibration factor per band in the maps-only path

- **Assumes**: a map read from file has `g0 = 0.0` and `gain = 0.0` and is taken as perfectly
  calibrated.
- **Where**: `file_io/map_reader.py:170-171`.
- **Kind**: silent convention.
- **Could bias**: indices directly. A relative calibration error dg between two bands moves a
  two-band index by dg / ln(nu2/nu1): 5 percent at 545 against 353 GHz is 0.12 in beta_d, 5
  percent at 857 against 353 GHz is 0.06, 10 percent at 408 MHz against 23 GHz is 0.025 in
  beta_s. More bands dilute this; the quoted calibration uncertainties of the 408 MHz map and of
  the two highest HFI bands are of this size.
- **Test or relax**: one multiplicative factor per band with a Gaussian prior from the stated
  calibration uncertainty, sampled or marginalized. With a flat-prior or ridge likelihood this is
  one more column in the per-pixel algebra only if linearized; otherwise one more scalar in the
  MH or NUTS block.
- **Status**: (c) as a model element. Q10 covers TOD gain only. Calibration is mentioned in
  passing in E031 (lines 226 to 228) and in `nextgen/DATA.md:157, 188-189`.

### D11. Masks and unobserved pixels

- **Assumes**: there is no analysis mask in the amplitude solve. A pixel is unobserved only if
  its inverse variance is zero. See G4 for the index-step mask.
- **Where**: `file_io/map_reader.py:145-148`.
- **Kind**: silent convention.
- **Could bias**: amplitudes at survey edges. If an input signal map holds zeros rather than the
  HEALPix bad value in unobserved pixels, `hp.ud_grade` averages those zeros into the coarse
  pixel while the summed inverse variance says the pixel is observed (to be verified; it depends
  on the input files). Not relevant for full-sky bands.
- **Test or relax**: a partial-sky band through the reader, inspecting edge pixels.
- **Status**: (c).

### D12. Single precision by default

- **Assumes**: `double_precision` defaults to false, "enough in practice". Maps are read as
  float32. The maps-only reader builds the `DetectorMap` without the switch, so the data and the
  inverse variance are float32 whatever `compsep.double_precision` says. The dot product in the
  CG accumulates in float32 unless double precision is on. The conditional MH likelihood forms
  `np.sum(whitened_residual**2)` on float32 arrays, while the chi-square diagnostic uses
  `dtype=np.float64`. The per-pixel `marginal` and `ridge` paths cast to float64 and are not
  affected.
- **Where**: `compsep/processing.py:603-604`, `params/param_explanation.yml:314-320`,
  `file_io/map_reader.py:47, 168-169`, `data_models/detector_map.py:129-130, 176-177`,
  `math_utils/arithmetic.py:129`, `compsep/mcmc.py:158-160`, `compsep/chisq.py:107`,
  `tod/gain.py:312, 315`.
- **Kind**: approximation. Possible bug for the conditional likelihood (B15).
- **Could bias**: convergence of the `conditional` index step, and CG convergence at high
  dynamic range. A float32 sum S is resolved to about S times 1.2e-7. MH needs the chi-square
  difference to much better than 1. One band's chi-square is about 3e6 at nside 512 with a
  reduced chi-square of 1 (resolution 0.25), 2.4e7 for WMAP K at the E014 full-sky reduced
  chi-square of 481 at nside 64 (resolution 2), and above 1e9 for a badly fitting nside-512 band
  (resolution 64 or worse). In that regime the acceptance decision is partly rounding.
- **Test or relax**: evaluate `local_loglike` on one real band in both precisions. Relax with
  `dtype=np.float64` in that one sum, which costs nothing.
- **Status**: (c).

### D13. Degrees of freedom in the chi-square diagnostic

- **Assumes**: `ndof` is the number of observed band pixels.
- **Where**: `compsep/chisq.py:104-108`.
- **Kind**: silent convention.
- **Could bias**: nothing in the chain. It affects how a reduced chi-square is read. The count is
  right for a sample from the flat-prior conditional (expected chi-square equals N), not for the
  posterior mean or for band-limited draws (F7 found 0.57 = (7 - 3)/7 for a mean).
- **Status**: (c), minor.

### D14. Per-pixel draws are stored as band-limited alms

- **Assumes**: the per-pixel amplitude map is passed through `map2alm_healpix(niter=3)` at the
  component `lmax` and kept as alms.
- **Where**: `compsep/perpix_solver.py:117-128`.
- **Kind**: approximation.
- **Could bias**: the amplitude sample is no longer a draw from its conditional once
  re-synthesized, which matters for the `conditional` likelihood and for the chi-square.
- **Status**: (a) F7 (E007). Q6 is open.

### D15. The TOD side sees per-pixel amplitudes at the wrong beam

- **Assumes**: when the per-pixel solver produced the sky, the TOD side realizes it at
  `max(band FWHM, amp_fwhm_rad)`. The comment says: "we simply have to use that, even though
  it's wrong."
- **Where**: `mpi/transfer.py:73-77`.
- **Kind**: approximation, documented in a comment.
- **Could bias**: gain and correlated-noise estimates in an end-to-end run that uses the
  per-pixel solver, because the sky subtracted from the TOD is smoother than the sky the
  detector saw.
- **Test or relax**: do not combine the per-pixel solver with TOD sampling, or deconvolve to the
  band beam with a regularized inverse.
- **Status**: (c).

## 4. Sky model

Implemented component classes: CMB, ThermalDust, Synchrotron, FreeFree, SpinningDust,
RadioSources (`params/param_explanation.yml:558-560`). All diffuse amplitudes are alms in uK_RJ
at a reference frequency, and `get_sed(nu)` is the Rayleigh-Jeans ratio between `nu` and `nu_ref`
(`params/param_explanation.yml:552-554`).

### S1. Nuisance SED parameters are fixed at parameter-file values

- **Assumes**: thermal dust is a single modified blackbody with a fixed temperature (20 K in the
  reference file) unless the hierarchical path samples a temperature field; free-free has a fixed
  electron temperature (7000 K); spinning dust is one tabulated spectrum shifted to a fixed peak
  frequency (25.3 GHz); synchrotron is a pure power law with no curvature or break; each point
  source has the spectral index of its catalogue entry. Only `beta` of dust and synchrotron can
  be sampled in the scalar path.
- **Where**: `sky/diffuse_components.py:689-690` (dust), `:814` (free-free), `:868-871`
  (spinning dust), `sky/point_sources.py:122-123, 207`; `params/param_explanation.yml:614-615,
  705, 721, 738`.
- **Kind**: approximation.
- **Could bias**: indices absorb what the fixed parameters cannot. A wrong dust temperature moves
  beta_d (F34 to F36 quantify this for the hierarchical path). A wrong electron temperature
  changes the free-free slope between 23 and 70 GHz by a small amount that is degenerate with the
  synchrotron index. A wrong spinning-dust peak is degenerate with both in the 20 to 40 GHz
  range, which E031 identifies as the place where real intensity synchrotron is limited by the
  component model rather than by noise.
- **Test or relax**: simulate with a perturbed fixed parameter (T_e 5000 against 8000 K, peak 20
  against 30 GHz, a curved synchrotron law) and fit with the reference values; report the index
  shift in units of its posterior width. Relax with scalar MH parameters for T_e and the peak,
  which the existing scalar machinery can carry once the components expose them.
- **Status**: dust temperature (a) F34, F35, F36. Curved or broken laws and missing components
  (b) section 5 of FINDINGS. Fixed T_e and fixed spinning-dust peak (c).

### S2. One spectral index for intensity and polarization

- **Assumes**: the I view and the QU view of a component declared `IQU` form one index
  parameter. The grouping key is `(type(comp), comp.comp_name)`, and the likelihood object
  "returns the summed log likelihood of both polarization streams". The hierarchical amplitude
  prior is intensity only: the QU view keeps a flat prior, and the region likelihood raises an
  error for a `hier` prior on a polarized stream.
- **Where**: `compsep/spectral_index.py:22-31, 100-124`, `compsep/perpix_solver.py:164-169`,
  `compsep/index_field_sampler.py:372-377`, `compsep/index_field.py:554-557`.
- **Kind**: silent convention.
- **Could bias**: polarized foreground amplitudes and the polarized CMB. Intensity has far more
  signal to noise, so the shared index is set by intensity, where it is also contaminated by
  free-free and spinning dust that polarization does not see. The polarized synchrotron index is
  then wrong by whatever intensity model error exists, with a posterior width that reflects the
  intensity data. BeyondPlanck and Cosmoglobe fitted the polarized synchrotron index separately
  for this reason.
- **Test or relax**: simulate with different I and QU indices (for example -3.0 and -3.2), fit
  with the shared parameter, and measure the residual in the polarized CMB. Presumably the model
  can already be relaxed by declaring two components, one `I` and one `QU`, with different names;
  I did not check that this works.
- **Status**: (c) as a code assumption. Every experiment in `nextgen/` so far is intensity only
  (F22 is the nearest neighbour).

### S3. No monopole, dipole, zero level or per-band offset, and no positivity

- **Assumes**: the map-level data model is d_b = sum_c SED_c(nu_b) a_c + n_b with nothing else.
  There is no offset per band, no dipole template, no monopole component. In the CG path every
  component carries its own l = 0 and l = 1 alms, which are fitted like any other mode. In the
  per-pixel path whatever zero level a band has is fitted as sky. Amplitudes are unconstrained in
  sign. A search for monopole, dipole, offset or zero level in `compsep/`, `sky/`,
  `file_io/map_reader.py` and `data_models/detector_map.py` finds only a plotting helper
  (`sky/diffuse_components.py:658-670`).
- **Where**: absence; the placeholders are `sky/template_component.py:1-16`.
- **Kind**: silent convention for offsets; a model choice for positivity.
- **Could bias**: indices, strongly. F17 (E014) measured index shifts of 0.03 to 0.4 from
  zero-level choices, hundreds of posterior standard deviations. Without positivity, a noisy
  low-signal pixel can hold a negative free-free or dust amplitude that is compensated by the
  CMB.
- **Test or relax**: decided in D013 (per-band offsets with informative priors), not
  implemented. A per-band offset is linear, so it is one more column in the amplitude solve.
- **Status**: offsets (a) F17, Q15, D013, section 5. Positivity (c).

### S4. `get_sky` evaluates an index-field component at the mean index

- **Assumes**: `DiffuseComponent.get_sky` returns `get_component_map(nside, fwhm) *
  get_sed(nu)`, with the scalar `self.beta` (and `self.T`). For a component that holds a
  `beta_map`, `set_beta_map` sets `self.beta` to the map mean, and its docstring says this is so
  that "chain output, `mixing/<band>` and the scalar code paths see a representative value". The
  projection operators used by the CG and the per-pixel solver go through `_band_sed`, which does
  use the map. `get_sky` does not.
- **Where**: `sky/diffuse_components.py:442-450` (the only diffuse `get_sky`), `:475-491`,
  `:493-508`, `:510-533`. Callers of `get_sky` through `SkyModel.get_sky_at_nu`:
  `compsep/chisq.py:91-92` (the chi-square diagnostic), `compsep/processing.py:541-542` (the
  fixed-component residual), `compsep/mcmc.py:154-156` (the conditional MH likelihood),
  `mpi/transfer.py:77` (the sky handed to the TOD side).
- **Kind**: approximation, stated in a docstring, with consequences the docstring does not list.
  Possible bug (B3).
- **Could bias**: for any run with an index field: (1) the chi-square and reduced chi-square that
  Commander4 reports, which are then computed for a different sky from the one sampled; (2) the
  residual seen by another group when the field component is held fixed; (3) the sky subtracted
  on the TOD side, hence gain and noise; (4) a `conditional` index step on any other parameter.
  The index-field sampling itself goes through the region kernels and is presumably unaffected.
  The error is first order in the index spread: a spread of 0.1 in beta_s moves the predicted
  synchrotron at 408 MHz against 23 GHz by exp(0.1 times 4.0) - 1 = 49 percent.
- **Test or relax**: take one saved sample of an index-field run and compute the chi-square two
  ways, with `get_sky_at_nu` and with a per-pixel SED. If they agree, something I did not see
  handles it. Relax by using `_band_sed` logic inside `get_sky`.
- **Status**: (c).

### S5. The index field is constant inside coarse regions, with a stationary Gaussian prior

- **Assumes**: one index per HEALPix pixel at `spectral_index_nside` (at most 16), copied to
  every data pixel inside it with `hp.ud_grade`. The field prior is a Matern Gaussian Markov
  random field on the 8-neighbour graph, alpha = 2, stationary and isotropic, with the constant
  mode handled separately.
- **Where**: `sky/diffuse_components.py:510-526`, `math_utils/gmrf_field.py:1-14, 145-183`,
  `params/param_explanation.yml:634-652`.
- **Kind**: exact for a stated prior. The prior is the assumption.
- **Could bias**: indices and their widths where the true index varies inside a region or where
  its statistics differ between the plane and high latitude. Sub-region variation acts as an
  unmodelled SED curvature in the region average.
- **Test or relax**: simulate with a truth that varies inside regions and check coverage of the
  region means.
- **Status**: (a) F24, F29, F32, F38 for the prior and the resolution ladder. I found no test
  with truth varying inside a region; D019 withdraws the resolution rule until that exists.

### S6. Point sources are a fixed catalogue

- **Assumes**: amplitudes, positions and indices are read from a file and never sampled ("there
  is no amplitude, position or spectral-index sampling for sources yet"). The CG right-hand side
  carries "FIXME: how will eta_2 work for point sources?". The reference file has the component
  disabled.
- **Where**: `sky/point_sources.py:116-124, 207`, `params/param_explanation.yml:742-746`,
  `compsep/cg_solver.py:221`.
- **Kind**: approximation, documented.
- **Could bias**: every diffuse amplitude near a source, and indices through the band-dependent
  footprint of a source after smoothing (see D5). With the component disabled, sources are
  fitted as diffuse emission at each band's own beam.
- **Test or relax**: mask sources in the index step as a stopgap; sample amplitudes as linear
  parameters.
- **Status**: (c).

### S7. The reference frequency becomes part of the prior as soon as the prior is not flat

- **Assumes**: amplitudes are stored at `nu_ref`, per polarization view (dust `[857, 353]`,
  synchrotron `[0.408, 30]` in the reference file). The CMB class comments that the choice "is
  arbitrary (the sky is invariant to it)".
- **Where**: `params/param_explanation.yml:552-554, 616, 706`,
  `sky/diffuse_components.py:624-626, 638-640`.
- **Kind**: exact as a reparameterization under a flat amplitude prior. Silent convention under
  any other prior: the unit-variance prior of A2, a `gaussian` prior, and the `hier` prior are
  all stated in amplitude units at `nu_ref`, and the Jacobian between two reference frequencies
  depends on the index.
- **Could bias**: indices (through the volume term of F3 and F12, which changes with `nu_ref`)
  and amplitudes (through A2).
- **Test or relax**: the same run at two reference frequencies; any difference beyond Monte
  Carlo error is prior dependence.
- **Status**: `hier` template reference (a) E030 pass 7, and F13 for the DR1 convention. Under
  the CG null prior (c).

### S8. Template components do not exist

- **Assumes**: `TemplateComponent` and `CMBRelQuad` raise `NotImplementedError`; "nothing
  constructs them yet".
- **Where**: `sky/template_component.py:1-16`.
- **Kind**: missing model element.
- **Could bias**: anything a template would absorb: zero levels (S3), zodiacal light, the
  relativistic quadrupole, CO, far sidelobes at map level.
- **Status**: (b) section 5 of FINDINGS for the missing components. Section 5 also lists
  templates among what is covered, which does not match this file (B12).

### S9. Priors on the spectral index

- **Assumes**: optional hard bounds and an optional Gaussian prior per component. The default is
  neither, which is a flat improper prior. The normalization of the Gaussian is dropped, which is
  correct because it is constant.
- **Where**: `compsep/spectral_index.py:22-31, 88-97, 135-162`,
  `params/param_explanation.yml:622-633`.
- **Kind**: exact for a stated prior.
- **Could bias**: indices at low signal to noise, where the flat-prior marginal is improper or
  nearly so and the answer is set by the volume term.
- **Status**: (a) F3, F4, F12, F24, F29.

## 5. Amplitude solves

The CG path solves, per polarization stream,

    (1 + S^(1/2) A^T N^-1 A S^(1/2)) x = S^(1/2) A^T N^-1 d + S^(1/2) A^T N^(-1/2) eta_1
                                         + S^(-1/2) mu + eta_2,        a = S^(1/2) x,

with A = B Y^-1 M Y per band, N diagonal in pixel space, and S the C(l) prior
(`compsep/cg_solver.py:43-53`). The per-pixel path solves (M^T N^-1 M) a = M^T N^-1 d + M^T
N^(-1/2) eta in every pixel (`src/lib_cpp/compsep_pymod.cc:103-121`).

### A1. A CG solve stopped at the iteration cap is accepted as a sample

- **Assumes**: the solve starts from zero every Gibbs iteration and stops when either
  `iter >= max_iter` or `err < err_tol`, where `err = next_rz/rz0` is the squared preconditioned
  residual relative to its starting value. The cap produces a warning only if `err_tol > 0`. The
  documentation recommends `err_tol: 0.0` with 100 iterations, "what every real-data file here
  uses, because the CG residual is not a reliable convergence measure for a constrained
  realization", and the upstream real-data example does exactly that. With `optimize: true` the
  fluctuation terms are dropped and the result is the Wiener mean.
- **Where**: `compsep/cg_solver.py:274, 344-350, 414-415`, `compsep/cg_driver.py:80`,
  `compsep/cg_solver.py:93`, `params/param_explanation.yml:341-351`,
  `params/Planck_mapsonly/param_LFI+HFI+WMAP.yml:69-73`.
- **Kind**: approximation. The returned vector is a draw from the amplitude conditional only if
  the system is solved. After k iterations from zero it is the projection of the exact draw onto
  a Krylov space, so the directions the preconditioner serves worst are biased toward zero and
  under-dispersed. The Gibbs chain then targets a distribution nobody has written down.
- **Could bias**: the CMB map and all foreground amplitudes in the slow directions, which for a
  preconditioner that assumes isotropic noise and constant mixing (A5) are expected to be the
  large scales under strongly varying noise, the masked or unobserved regions, and the
  near-degenerate component combinations. I do not know which of these are unconverged at 100
  iterations; the code prints the residual and stores `CG_residuals`, but nothing records them
  in `nextgen/`. A zero start also forgoes the previous sample as a warm start, which would be
  exact if the solve converged and matters only because it does not.
- **Test or relax**: solve one right-hand side at 100 iterations and to a tight tolerance and
  compare per multipole and per component (same eta, so the difference is pure truncation). At
  nside 8 to 16 compare with a dense solve; `dense_matrix_debug_mode` exists for this but is
  documented as "currently broken". Relax by a convergence rule on a quantity that means
  something (the A-norm error estimate, or the change in the solution on the modes of interest),
  and by reporting the residual in the chain file.
- **Status**: (c). Q12 asks what the CG costs, not whether its output is a sample.

### A2. "No C(l) prior" is a unit-variance Gaussian prior per alm

- **Assumes**: with `Cl_prior_amplitude: null`, `P_Cl_prior` returns ones, so S is the identity:
  a prior N(0, 1) on every alm, in (uK_RJ at `nu_ref`)^2. The identity term is always part of the
  operator ("Adds input vector to output, since (1 + S^{1/2}...)a = a + (S^{1/2}...)a") and the
  prior fluctuation eta_2 is always drawn. The code docstring says this plainly ("a Gaussian
  prior constraining the alms to N(0, C_l)"). The parameter documentation says the null setting
  "switches the prior off entirely (S = identity), which is what most files here do". Fifteen
  parameter files under `params/` use it.
- **Where**: `sky/diffuse_components.py:215-216, 235-236`, `compsep/cg_solver.py:191-193,
  262-268`, `params/param_explanation.yml:580-581, 611`.
- **Kind**: silent convention. Possible bug in the sense of a prior nobody chose (B5).
- **Could bias**: amplitudes of any component, depending on a ratio that the user does not see.
  The data precision per alm is sum_b SED_b^2 b_l^2 / (sigma_b^2 omega_pix), against a prior
  precision of 1. Hand estimates:
  1. In the upstream real-data example (synchrotron at 30 GHz, dust at 545 GHz, CMB at 1 GHz)
     the data precision at low multipoles is hundreds to thousands, so the prior shrinks
     amplitudes by 0.1 percent or less there. It takes over for synchrotron above a multipole of
     roughly 400 to 600, where it acts as a regularizer, and stays at 0.3 percent or less for
     the CMB up to l = 2000.
  2. As a statement about the sky it is a white spectrum: D_l = l(l+1)/(2 pi), which is 17.5,
     1607 and 1.6e5 uK^2 at l = 10, 100 and 1000, and about 100 times the true EE power at
     l = 5. Wherever the data do not constrain a mode (unobserved pixels, multipoles past the
     beam, degenerate component combinations) the sample is this white field, with a pixel rms of
     (lmax+1)/sqrt(4 pi): 289 uK at lmax 1024.
  3. The answer depends on `nu_ref`. The reference parameter file enables the CG group and gives
     synchrotron a null prior with `nu_ref: [0.408, 30]`. At 408 MHz the data precision for the
     intensity amplitude is tiny in these units: about 4e-7 from Haslam if its rms is 0.8 K per
     nside-512 pixel (1e-8 at 5 K), plus about 1e-9 from WMAP K. A prior precision of 1 would
     then outweigh the data by a factor of a million or more and hold the synchrotron intensity
     at essentially zero. This is arithmetic from the code as read, not an observation, and such
     a failure would be obvious in any output, so either nobody runs that combination or I am
     missing something. It must be checked before the reference file is used as a template
     (B5).
  4. The null branch returns before the apodization is applied, so `Cl_prior_l_apod` does
     nothing for a null prior (B6).
- **Test or relax**: a dense solve at nside 8 or 16 with one component at two reference
  frequencies, flat prior against S = identity (pure numpy, seconds). Relax by a true flat-prior
  option (drop the identity term and eta_2 for that component, which needs a different
  preconditioner scaling), or by requiring an explicit prior in physical units.
- **Status**: (c). No file in `nextgen/` discusses what the null prior is.

### A3. The C(l) prior is a fixed power law, the same for E and B, with no TE

- **Assumes**: D_l = amplitude (l/l_pivot)^beta with an optional Gaussian roll-off and a high-l
  taper, given per polarization view, applied identically to both rows of a QU view. It is
  "Fixed, not sampled". The l = 0 entry is set to D_1 (that is, C_0 = C_1/pi).
- **Where**: `sky/diffuse_components.py:211-246, 248-260, 272-284`,
  `params/param_explanation.yml:576-590`.
- **Kind**: approximation, documented.
- **Could bias**: the CMB map in low signal-to-noise modes and, for polarization, B modes: a
  prior with E power applied to B puts E-sized variance into unconstrained B modes. Without TE,
  the temperature data do not inform E. With a fixed prior the CMB samples are not samples from
  the joint posterior of map and spectrum.
- **Test or relax**: C(l) sampling (Q11); separate E and B priors; a TE block.
- **Status**: C(l) sampling (b) section 5, Q11. Same E and B, no TE, and the monopole
  convention (c).

### A4. Component multipoles above every band limit are drawn from the prior alone

- **Assumes**: a component `lmax` above the highest band `lmax` is allowed. The validator logs a
  warning or an error and continues: "Both are logged rather than raised".
- **Where**: `compsep/processing.py:476-520`, `data_models/detector_map.py:122-126`,
  `params/param_explanation.yml:565-569`, `parameters/schema.py:275-295`. The upstream
  real-data example is not affected: its HFI bands are at `eval_nside: 2048`
  (`params/Planck_mapsonly/maps_HFI.yml`), so the default band limit of 3 nside - 1 = 6143 is
  above the component limits of 2048 and 4096.
- **Kind**: silent convention, recognized upstream.
- **Could bias**: amplitude maps and the chi-square at high multipoles; nothing at low
  multipoles.
- **Status**: (c) in `nextgen/`.

### A5. The preconditioner assumes isotropic noise, constant mixing and Gaussian beams

- **Assumes**: the joint preconditioner uses the sky mean of the inverse variance per band, the
  trace part of the polarization weight, scalar SEDs ("Mixing is assumed to be spatially
  constant") and `hp.gauss_beam`. Three other preconditioner names in the documentation "are
  unimplemented and crash mid-solve".
- **Where**: `compsep/preconditioners.py:222-255`, `params/param_explanation.yml:352-359`.
- **Kind**: exact as a preconditioner (it changes speed, not the solution) if the solve
  converges. Under A1 it decides which modes are solved.
- **Could bias**: with A1, the modes where real noise departs most from isotropy: the scan
  pattern on large scales, the mask, and the plane for signal-dependent rms.
- **Test or relax**: as A1. A pseudo-inverse or multigrid preconditioner for the masked
  low-multipole block is the standard remedy.
- **Status**: (c).

### A6. Harmonic analysis by adjoint synthesis

- **Assumes**: `map_to_alm` is adjoint synthesis times 4 pi / npix ("using only the scalar
  normalization factor"), which is not the inverse of synthesis on HEALPix. Initial states read
  from maps use an iterative inverse capped at 5 iterations.
- **Where**: `math_utils/sht.py:74-79`, `sky/comp_io.py:72-74`.
- **Kind**: exact where it is used as the adjoint inside A^T N^-1 A. Approximation where it is
  used as an inverse, which includes the Y^-1 in the mixing operator when the mixing varies over
  the sky.
- **Could bias**: amplitudes at multipoles above about 2 nside in the spatially varying CG path.
- **Status**: (c), minor.

### A7. The per-pixel draw has no guard for pixels with fewer bands than components

- **Assumes**: every pixel with a non-negative region index is solved. The Cholesky routine
  reports failure only for a pivot `<= 0.0`, and failure gives zero amplitudes. The two
  likelihood kernels skip pixels with `nobs < ncomp`; the amplitude kernel does not count
  observed bands at all.
- **Where**: `src/lib_cpp/compsep_pymod.cc:25-37, 97-121` against `:269-272` and `:402-405`;
  conditioning TODO at `compsep/perpix_solver.py:97-98`.
- **Kind**: silent convention. Possible bug (B4).
- **Could bias**: amplitudes in partially observed pixels. With 0 < nobs < ncomp the normal
  matrix is singular and the flat-prior conditional is improper. In exact arithmetic the pivot
  is zero and the pixel is zeroed; in floating point it can be a small positive number, and the
  amplitudes are then enormous. Either way the map goes through `map2alm_healpix` afterwards
  (D14), which spreads the damage. Fully unobserved pixels give an exact zero and are safe. So
  the issue is latent for full-sky bands and live for any partial-sky band (ground-based data,
  masked bands) and for ill-conditioned mixing (the TODO).
- **Test or relax**: a toy call to `cpp_compsep.solve_perpix` with one band's inverse variance
  set to zero in some pixels. Relax with the same `nobs < ncomp` rule as the likelihood kernels
  plus a relative pivot threshold.
- **Status**: (c) for the draw. The likelihood-kernel skip is noted in
  `nextgen/experiments/E009_c3_nonlin_review/independent_reviews_2026-10-25.md:390`.

### A8. How the CMB prior enters each path

- **Assumes**: in the CG path the CMB has whatever S it is given: a fixed power law in the
  reference file, the unit-variance prior of A2 in the upstream real-data example
  (`params/Planck_mapsonly/param_LFI+HFI+WMAP.yml:115`). In the per-pixel path the CMB amplitude
  has a flat prior per pixel, unless a `gaussian` or `hier` amplitude prior is configured, in
  which case the prior enters as a pseudo-band. No path uses a physical C(l).
- **Where**: `compsep/cg_solver.py:43-53`, `compsep/perpix_solver.py:33-36, 133-152`,
  `compsep/index_field.py` (`augment_stream`, `augment_mixing`).
- **Kind**: exact for a stated prior.
- **Could bias**: the CMB map. F26 shows that a finite per-pixel CMB prior leaks the CMB
  realization into the foreground amplitudes.
- **Status**: (a) F26, D021.


## 6. TOD side

Read briefly, as asked: enough to list the modelling choices that reach the maps and the
large-scale polarization, not a line-by-line review. As the comments describe it, the time-domain
model is gain times (sky + orbital dipole + far sidelobes), passed through an optional bolometer
transfer function, plus correlated noise plus white noise, with glitch and jump handling before
that. Every sampling step below is disabled by default and has to be switched on.

### T1. Gain

- **Assumes**: the gain of a detector is the sum of an absolute, a relative and a temporal term.
  All three steps default to `"enabled": False`. The absolute gain is calibrated against the
  orbital dipole by default and is a Gaussian draw with a flat prior; above 400 GHz the code
  only warns that the dipole "makes a poor calibrator". Relative and temporal gain are calibrated
  against the full sky model. The temporal gain is a Wiener-filter draw whose prior spectrum is
  hard-coded (`alpha_gain = -2.5`, `fknee_gain = 1.0`, `sigma0_gain = 1e-4*mean_gain`, with the
  comment "I found the BP prior too weak, and this gives me more sensible results"). It is solved
  by CG from the mean, with a cap of 200 iterations and an exit at `err < 1e-10`, and with no
  warning or fallback when the cap is reached. The mean of the temporal draw is removed
  afterwards. A comment leaves the treatment of masked samples undecided. The relative-gain
  accumulators are float32.
- **Where**: `tod/gain.py:37, 128-145, 269-276, 312-315, 429-434, 482-488, 527-534`.
- **Kind**: approximation (hard-coded prior, unchecked CG) and silent convention (defaults off;
  with all terms off the gain is the value in the data file, with no uncertainty).
- **Could bias**: the absolute calibration of every band, hence indices in the same way as D10;
  and the large-scale polarization, because a gain error turns the dipole and the Galactic
  intensity into spurious large-scale signal that the scan pattern couples into Q and U.
  Calibrating relative and temporal gain against the sky model makes the gain depend on the
  foreground model, which is the degeneracy Q10 asks about. A prior tuned by hand for "more
  sensible results" is a prior fitted to the data it is applied to.
- **Test or relax**: an end-to-end simulation with a known gain drift; write the temporal-gain
  CG residual to the chain; expose the prior spectrum as parameters and state it in the model.
- **Status**: the gain and sky correlation (b) Q10. The specifics above (c).

### T2. The white-noise level is estimated, not drawn

- **Assumes**: sigma0 of each detector and scan is a deterministic estimate: by default from
  first differences of the residual with outlier rejection (`pairwise`), re-estimated after the
  correlated noise is subtracted, or from the bottom of the binned PSD (`binned_psd`). The
  switch is called `sample_sigma0` and the method docstring says "Estimate sigma0 alone".
- **Where**: `tod/noise/sampling.py:63-70, 159-164`, `tod/noise/sample_ncorr.py:302-306,
  347-350`.
- **Kind**: approximation (a point estimate in place of a conditional draw).
- **Could bias**: uncertainties, by a negligible amount through the missing scatter (sigma0 is
  known to about 1/sqrt(2 n) from n samples), and every inverse variance handed to CompSep
  through any bias of the estimator. A first-difference estimator sees the noise at the sampling
  frequency only, so it differs from the level relevant for the map wherever the spectrum is not
  flat there.
- **Test or relax**: compare the two estimators per detector on real data. A draw from the
  inverse-gamma conditional is one line.
- **Status**: (c).

### T3. The default correlated-noise draw ignores the gaps

- **Assumes**: correlated-noise sampling is off by default. When it is on, the default
  `cg.max_iter: 0` takes "the stationary (full-mask) Wiener solution directly": the masked
  samples are first filled in (`fill_all_masked`), and the draw is then made with the mask
  replaced by ones, so the filled values are treated as data. A gap-aware CG exists; it falls
  back to the same stationary solution when it does not converge or when the residual carries
  more power than the data. Each scan is an independent stationary stretch on a mirrored FFT
  grid.
- **Where**: `tod/noise/sampling.py:63-70`, `tod/noise/sample_ncorr.py:317-341`.
- **Kind**: approximation. Exact only for a scan without masked samples.
- **Could bias**: the correlated-noise estimate next to masked stretches (Galactic plane, bright
  sources, glitches), and through it the large-scale maps, which is where residual 1/f noise
  ends up. Treating scans as independent drops noise correlations longer than a scan.
- **Test or relax**: simulate 1/f noise with a realistic mask and compare the stationary draw
  with the converged gap-aware CG. The code already stores `ncorr_converged` and
  `ncorr_high_var` per scan, so the fallback rate can be reported.
- **Status**: (c).

### T4. One noise PSD model, sampled on a grid

- **Assumes**: the only implemented PSD is P(f) = sigma0^2 (1 + (f/fknee)^alpha). Five other
  models are class stubs marked "TODO: Implement". When PSD sampling is on (it is off by
  default), fknee and alpha are drawn by inversion sampling on 150-point grids inside hard
  bounds, alternating for `n_burnin + 1` = 5 passes, from the periodogram of the residual with
  sigma0 held fixed.
- **Where**: `tod/noise/psd.py:256-361, 364-398`, `tod/noise/sampling.py:63-70`.
- **Kind**: approximation: a model restriction, a Whittle likelihood, a grid and a short inner
  chain.
- **Could bias**: the correlated-noise estimate, hence large-scale maps, for detectors whose
  noise is not a single 1/f law. From memory, the BeyondPlanck noise paper needed an extra
  log-normal term for the LFI 30 and 44 GHz radiometers, which is what the `NoisePSDOofGauss`
  stub is for; verify before citing.
- **Test or relax**: a posterior predictive check of the residual periodogram per detector;
  implement the stubs.
- **Status**: (c).

### T5. Known bug with a bolometer transfer function

- **Assumes**: for a band with a time constant (`tf_tau_sec`) the data hold T(sky + orbital
  dipole + sidelobes), "but the models subtracted here and in the other Gibbs steps are not
  passed through T, so the residuals and the CG right-hand side keep (1 - T) times them".
- **Where**: `tod/mapmaking/tod2map.py:120-124`, `tod/mapmaking/cg.py:74-79`.
- **Kind**: known bug, documented in the code (B9).
- **Could bias**: HFI-type bands in an end-to-end run: residuals, gain, correlated noise, and
  the maps along the scan direction.
- **Test or relax**: as the comment proposes, deconvolve the TOD and filter the models with the
  same regularization kernel.
- **Status**: (c) as a code fact. Transfer functions are on `nextgen/ROADMAP.md:56-58`.

### T6. What the maps carry to CompSep

- **Assumes**: samples are weighted by the white-noise level only, `inv_var =
  (gain/view.sigma0)**2`. The per-pixel IQU weight matrix is accumulated, but "Only its inverse
  diagonal survives in `map_inv_var`".
- **Where**: `tod/mapmaking/tod2map.py:144-146`, `tod/mapmaking/output.py:37-39`,
  `tod/mapmaking/binned.py:26-28`.
- **Kind**: exact given the correlated noise when correlated-noise sampling is on (the
  correlated part is then a sampled quantity that has been subtracted, and what is left is
  white). Approximation when it is off, which is the default, and for the correlations between
  I, Q and U in every case.
- **Could bias**: with correlated-noise sampling off the map noise is not white and CompSep is
  told that it is. In every case CompSep loses the Q and U correlation inside a pixel (D1).
  `1/diag(A^-1)` is the marginal variance of each Stokes parameter, so the diagonal that is kept
  is the right one; what is lost is the off-diagonal.
- **Test or relax**: pass the 3 by 3 block per pixel to the per-pixel kernel and to the CG
  operator.
- **Status**: (c).

### T7. The CG mapmaker

- **Assumes**: "The CG's own noise model is white". The solve is accepted at the iteration cap,
  with a warning only when the tolerance is positive. Cross-talk across several detectors is an
  open TODO ("how to deal with multiple detectors here???").
- **Where**: `tod/mapmaking/cg.py:77-82, 206-211, 372-383`.
- **Kind**: approximation.
- **Could bias**: maps of bands that need the transfer function or cross-talk deconvolution,
  when the solve has not converged.
- **Status**: (c).

### T8. The sky that the TOD steps condition on

- **Assumes**: the sky subtracted from the TOD is one map per band, `get_sky_at_nu` at the
  band's single frequency with a symmetric Gaussian beam at `max(band FWHM, amp_fwhm_rad)`. All
  detectors of a band see the same sky, so there is no per-detector bandpass and no model of
  bandpass-mismatch leakage. A component with an index field enters at its mean index (S4).
- **Where**: `mpi/transfer.py:73-77`, with D7, D15 and S4.
- **Kind**: approximation.
- **Could bias**: gain and correlated noise through the residual; the polarized maps through
  intensity leaking into polarization when detectors with different bandpasses are differenced.
  The latter was a leading systematic for LFI, WMAP and HFI polarization in the published
  analyses.
- **Test or relax**: per-detector bandpass shifts (the LFI reader carries the TODO, D7), with
  the leakage templates that follow from them.
- **Status**: (c) as a code fact. Bandpass mismatch is on `nextgen/ROADMAP.md:56-58`.

## 7. Possible bugs and inconsistencies

Every item here is a hypothesis from reading. None was reproduced, because nothing could be run.
For each one I give the reason it may turn out harmless, which I take as seriously as the reason
it may not. Items B1 to B5, B14 and B15 are the ones that could move real-data numbers.

### B1. Maps are degraded before they are smoothed (D5)

- **What**: `read_data_map_from_file` brings a map to `eval_nside` with `hp.ud_grade` and only
  then smooths it to the common resolution. Power below the coarse pixel scale is aliased before
  the smoothing can remove it, by an amount that depends on the band's native beam. An unrecorded
  pixel window is added. The variance assigned to the coarse pixel is the inverse of the summed
  inverse variances, which is not the variance of the unweighted mean that was taken.
- **Where**: `file_io/map_reader.py:156-162` before `:173-175`.
- **Why it may be harmless**: nothing happens when the input is already at `eval_nside`. The
  part of the aliasing that all bands share has the CMB spectrum and goes into the CMB
  amplitude. My size estimate for the band-dependent part at E014 settings (a few uK) is
  uncertain by a factor of a few and may be below the noise at high latitude.
- **Verify**: test 1 in section 11.

### B2. Q and U are smoothed as scalars (D4)

- **What**: the common-resolution smoothing calls the scalar smoother on each Stokes map.
- **Where**: `data_models/detector_map.py:171-175`.
- **Why it may be harmless**: it is exact in the flat-sky limit, and the error for a 120 arcmin
  kernel is 0.13 percent at 30 degrees from a pole. Every experiment in `nextgen/` so far is
  intensity only. It is not harmless for the 5 to 10 degree kernels of low-resolution
  polarization work (3.3 percent at 30 degrees from a pole for 10 degrees).
- **Verify**: test 6.

### B3. `get_sky` uses the mean index for a component with an index field (S4)

- **What**: the chi-square diagnostic, the fixed-component residual, the conditional MH
  likelihood and the sky handed to the TOD side are computed with one scalar index per
  component, also when the sampled sky has an index map.
- **Where**: `sky/diffuse_components.py:442-450`; the mean is set at `:475-491` and `:493-508`.
- **Why it may be harmless**: the docstring of `set_beta_map` says the mean exists for "the
  scalar code paths". The index-field sampler and both amplitude solvers go through `_band_sed`
  and the region kernels, so the sampled posterior is presumably right and only the listed
  consumers are affected. If no index-field run so far has used the Commander4 chi-square, a
  fixed field component or the TOD side, nothing recorded changes.
- **Verify**: test 3.

### B4. No guard for under-determined pixels in the per-pixel draw (A7)

- **What**: the amplitude kernel solves every pixel and relies on a Cholesky pivot `<= 0.0` to
  detect failure, while the two likelihood kernels skip pixels with fewer observed bands than
  components.
- **Where**: `src/lib_cpp/compsep_pymod.cc:25-37, 119-121` against `:269-272` and `:402-405`.
- **Why it may be harmless**: with full-sky bands every pixel has all bands. A pixel with no
  data gives an exact zero. An exactly singular matrix gives a zero pivot and zeroed amplitudes
  when the arithmetic happens to be exact.
- **Verify**: test 7.

### B5. A null C(l) prior is a unit-variance prior, documented as no prior (A2)

- **What**: `Cl_prior_amplitude: null` makes S the identity, which is a Gaussian prior of unit
  variance on every alm in (uK_RJ at `nu_ref`)^2. The parameter documentation says it "switches
  the prior off entirely (S = identity)" and annotates the setting with "No C(l) prior".
- **Where**: `params/param_explanation.yml:580-581, 611` against
  `sky/diffuse_components.py:215-216, 235-236` and `compsep/cg_solver.py:191-193, 262-268`.
- **Why it may be harmless**: wherever the data precision per alm is far above 1 the prior does
  nothing, and by my estimate that covers the low and intermediate multipoles of the upstream
  real-data example (0.1 percent or less). At high multipoles some regularizer is wanted anyway.
- **Why it may not be**: the strength of the prior depends on the unit and on `nu_ref`. For the
  reference file's synchrotron at 0.408 GHz my arithmetic gives a prior that outweighs the data
  by a factor of a million or more (A2, point 3). Unconstrained modes are filled with a white
  field of unit variance per alm. And no CG run has the flat amplitude prior that F3, F12 and
  F15 analyse, so those results do not transfer to the CG path without an argument.
- **Verify**: test 4.

### B6. A null prior bypasses the high-l apodization

- **What**: `P_Cl_prior` returns ones before the line that applies `Cl_prior_apodization`, so
  `Cl_prior_l_apod` does nothing for a null prior. The `lmax` validator nevertheless treats an
  apodization at or below the band limit as protection ("Only Cl_prior_l_apod = ... keeps those
  modes finite") and recommends setting it.
- **Where**: `sky/diffuse_components.py:235-236` before `:246`; `compsep/processing.py:513-520`.
- **Why it may be harmless**: with a null prior the unseen multipoles are drawn with unit
  variance, which is bounded, so the run does not blow up. The advice in the log message is
  simply without effect.
- **Verify**: evaluate `P_Cl_prior` for a component with a null amplitude and an apodization
  set. One line.

### B7. Two beam-window conventions for polarization (D6)

- **What**: `apply_B` multiplies every alm row by the spin-0 Gaussian window.
  `_realize_alms_as_map` calls `hp.smoothalm` with its default `pol=True`, which gives row 0 the
  spin-0 window and later rows the spin-2 window. A two-row QU view therefore has E smoothed
  with one window and B with the other, and neither matches `apply_B` for both rows.
- **Where**: `data_models/detector_map.py:213-219` against `sky/diffuse_components.py:429`.
- **Why it is small**: the windows differ by exp(2 sigma^2) - 1, which is 1.1e-4 at 1 degree and
  1.1e-2 at 10 degrees. This one is confirmed by reading; only its size is in question.
- **Verify**: realize a single E mode and a single B mode through both routines.

### B8. Latent fallback to the invalid legacy move (G2)

- **What**: the likelihood lookup defaults to `"resolve_chisq"` for a group missing from the
  resolved dictionary.
- **Where**: `compsep/processing.py:878-879`.
- **Why it may be harmless**: `_resolve_mcmc_likelihoods` presumably fills the dictionary for
  every group, so the default is never taken.
- **Verify**: make the lookup raise on a missing key and run the test suite.

### B9. Bolometer transfer function not applied to the models (T5)

- **What** and **where**: documented as a known bug at `tod/mapmaking/tod2map.py:120-124`.
- **Why it may be harmless**: it is live only for bands with `tf_tau_sec`, in end-to-end runs.

### B10. The TOD side realizes per-pixel amplitudes at the common beam (D15)

- **What** and **where**: `mpi/transfer.py:73-77`, with the comment "we simply have to use
  that, even though it's wrong".
- **Why it may be harmless**: it is live only when TOD sampling is combined with the per-pixel
  solver and common-resolution smoothing.

### B11. The temporal-gain CG is not checked after its loop (T1)

- **What**: 200 iterations at most, an exit at `err < 1e-10`, and no warning, flag or fallback
  when the loop ends without reaching it.
- **Where**: `tod/gain.py:527-531`.
- **Why it may be harmless**: the system is small and well conditioned for a smooth prior, and
  may always converge in far fewer than 200 iterations. Nothing records whether it does.
- **Verify**: log the final `err` and the iteration count.

### B12. Documentation that does not match the code

1. `params/param_explanation.yml:570-573` says no component returns a map-valued SED "so this is
   False everywhere". On `nextgen` index fields exist (`:634-652`), so the remark is stale.
2. `params/param_explanation.yml:322-326` says the signal and the rms go "through the same
   kernel". The signal goes through a noise-weighted kernel (D3).
3. `params/param_explanation.yml:352-359` offers three preconditioner names that "are
   unimplemented and crash mid-solve" and a debug mode that is "currently broken". The second is
   the tool one would use to check A1.
4. `params/param_explanation.yml:580-581, 611`: B5.
5. `nextgen/FINDINGS.md:141-149` lists templates among the components Commander4 covers.
   `sky/template_component.py:1-16` says both classes are placeholders and "nothing constructs
   them yet". The same paragraph lists "spatially varying spectral parameters" as missing,
   which is no longer true on `nextgen` (index fields under `likelihood: hier`), and per-band
   monopoles and dipoles as missing, which is right.
6. Low confidence, cannot be checked from here: `params/Planck_mapsonly/maps_WMAP.yml:8-11`
   reads the rms (`type: "rms"`) from a file named `..._N_v4.fits` with a first column called
   `II_Stokes`. If that column held a variance or an inverse-noise quantity instead of an rms,
   the WMAP weights would be wrong by a power. The files are on the Oslo cluster and I could not
   look at them. The name alone proves nothing.

### B13. Fixed global seed in a simulation script

- **What**: `np.random.seed(0)` followed by draws from the global generator.
- **Where**: `simulations/inplace_litebird_sim.py:61, 255`.
- **Why it may be harmless**: it is a standalone simulation tool, and a fixed seed is what one
  wants for a regression test. It matters only if several simulations made with it are treated
  as independent noise realizations.

### B14. The model-error floor that reaches the likelihood is far below the nominal one (D9)

- **What**: `add_signal_fraction_to_rms` is added per native pixel, before the degradation and
  the smoothing, so it is averaged down like white noise: about 8 times from nside 512 to 64 and
  about 3.3 times more from the smoothing, 26 in all at E014 settings. Its scale is the full-sky
  mean of the absolute map, which is set by the Galactic plane and the zero level.
- **Where**: `file_io/map_reader.py:141-143` before `:156-162` and `:173-175`. The E014 files
  read nside 512 maps at `eval_nside: 64` with the floor on the HFI bands
  (`nextgen/params/q1/param_q1_n64_ridge_npipetuned.yml:260-288`).
- **Why it may be harmless**: if the "rms (with floor)" column of the E014 high-latitude table
  was produced by the same reader after degradation and smoothing, then it is what the chain
  used and this item is void. The record says the table comes from a separate reconstruction
  and does not say how the column was made, so I cannot tell by reading.
- **Why I think it is not void**: three checks on recorded numbers. The reference values are my
  own, from memory, and good to perhaps 20 percent.
  1. A nominal 5 percent of the full-sky mean absolute map is about what the column lists for
     353 and 857 GHz: 13.4 and 23 uK_RJ correspond to means of about 1 and 10 MJy/sr, which is
     what I expect for those maps. A floor 26 times lower would be 0.5 and 0.9 uK_RJ.
  2. The column lists 3.3 uK_RJ at 30 GHz and 3.1 at 70 GHz. Those are the white-noise levels of
     the Planck maps in one nside-64 pixel (2.5 and 3.5 uK deg in CMB units) with no reduction
     from the smoothing. Haslam fits the same pattern: 1.3e5 uK is the documented 1 percent plus
     0.8 K divided by 8. So the column seems to be the degraded rms plus the nominal floor,
     while the chain propagates the smoothing into the variance
     (`data_models/detector_map.py:42-69`). By that reading the chain had every band about 3
     times below the column and the floored bands 13 to 26 times below it.
  3. Commander4's own total reduced chi-square in that record is 15000, while the table's
     full-sky values average about 100 over the 17 bands. Scaling the table by the factors
     above gives about 7000. The record attributes the gap to the different smoothing in the
     plane, which may also contribute.
- **Why it matters if true**: the E014 reading of which bands carry the fit changes. By the
  table the three HFI bands are slack (reduced chi-squares of 0.0 to 0.1 at high latitude). By
  my reading they were the tightest bands in the chain, at about 0.4, 0.5 and 0.9 uK_RJ per
  pixel against 1 to 2 uK_RJ for WMAP K and the LFI bands, with high-latitude reduced
  chi-squares of order 6 to 13, and the dust index was set by them at full weight. The same
  mechanism applies to the Haslam rms, which is a model-error budget written per nside-512
  pixel and is averaged down about 23 times.
- **Verify**: test 2.

### B15. Single-precision chi-square in the conditional MH likelihood (D12)

- **What**: `SkyModel.get_sky_at_nu` allocates float32, the data and the inverse variance are
  float32, and `local_loglike` returns `-0.5 * float(np.sum(whitened_residual**2))` with the sum
  taken in float32. The sums are combined across bands in double precision afterwards.
- **Where**: `compsep/mcmc.py:153-160`, `sky/sky_model.py:63-68`, `file_io/map_reader.py:47`,
  `data_models/detector_map.py:176-177`.
- **Why it may be harmless**: the rounding is a deterministic function of the state, so the MH
  step keeps detailed balance for a slightly rough version of the target, and for small sums the
  roughness is negligible. With the E014 mask (38165 pixels) a band's sum is 1e5 to 1e6 and the
  float32 spacing is 0.01 to 0.06. The `marginal`, `ridge` and `hier` paths cast to float64 and
  are not affected, and those carry every recorded result.
- **Why it may not be**: the CG path resolves `auto` to `conditional`, and that is the
  high-resolution path. In the upstream beta-sampling example
  (`params/Planck_mapsonly/param_LFI+HFI+WMAP_betasamp.yml:61, 78-89, 98`: single precision, CG
  amplitudes, `likelihood: auto`, HFI at nside 2048, the mask commented out) a band's sum is at
  least 5e7, where the float32 spacing is 4 and the accumulated summation error is of that order
  or larger. Chi-square differences of order 1, which is what a tuned MH step produces, are then
  not resolved. That configuration has a bigger problem first: by F14 the conditional cannot
  mix, and the example's proposal sigma of 0.01 is thousands of conditional standard deviations
  going by the E014 numbers.
- **Verify**: test 8. The fix is `dtype=np.float64` in one sum.

### B16. A run started with `init_chain` reuses the random streams of the run it starts from (G8)

- **What**: the generator is keyed on (seed, chain, iteration, world rank). `init_chain`
  "starts new chains (at iteration 1)" from a saved sample of another run. With the same seed and
  MPI layout the new run therefore draws the same streams as iterations 1 onward of the source
  run. If the new run extends the old one with an unchanged model, the two stretches are driven
  by identical random numbers, and for a contracting sampler the second converges onto the
  first.
- **Where**: `cli.py:52-63`, `file_io/chain_restart.py:9-10, 134`.
- **Why it may be harmless**: the documented way to extend a chain is `chain_status: append`,
  which continues the iteration numbering and is not affected. `init_chain` is meant for starting
  a different analysis from an old state, where shared random numbers correlate the runs but do
  not duplicate samples. No file in `nextgen/` uses `init_chain`.
- **Verify**: test 10. The fix is to put the source iteration or a run identifier into the seed
  sequence, or to refuse an unchanged seed when `init_chain` is set.

## 8. Cross-check against FINDINGS.md

The classes are: (a) a recorded experiment tested the assumption, and the finding is cited; (b)
the assumption is named somewhere in `nextgen/` (a remark in a finding, an open question, a
decision, the roadmap or an experiment record) but no experiment tests it; (c) I found no
mention. The search behind (c) is a keyword grep over the markdown files under `nextgen/`
without `reviews/` and `research/`, so a (c) can be wrong if the topic is recorded under words I
did not think of.

### 8.1 By finding

| Finding | Subject | Audit item | Comment |
|---|---|---|---|
| F1, F2 | The legacy index move is invalid, and the fix | G2 | Tested. The latent default to the legacy move (B8) is new |
| F3, F4 | Flat-prior marginal against ridge | G3, S9 | Tested for the per-pixel path. Does not transfer to the CG path, whose prior is never flat (A2) |
| F5 | Gibbs alternation is the mixing bottleneck | G5 | Tested (E001) |
| F6 | Language is not the bottleneck | none | Efficiency only |
| F7 | Per-pixel draws stored as band-limited alms | D14 | Tested (E007). Q6 open |
| F8 | Smoothed noise treated as white | D2 | Tested on simulations (E006). Not measured for the real-map setup. Q7 open |
| F9 | Two-chain pipelining is valid | G1 | Established by code reading. This audit read the same code and agrees |
| F10, F11 | Commander3 index samplers and bandpass shift | none | Commander3 only. F11 is a neighbour of D7 |
| F12, F15, F16 | Spatially varying index, joint target, calibration of the posterior | G3 | Tested on simulations drawn from the model |
| F13 | Reference frequency of the DR1 dust map | S7 | Shows what a reference-frequency error costs. Its closing remark names bandpass (D7) without testing it |
| F14 | The conditional cannot mix on real maps | G5 | Tested (E014). B15 sits in the same code path |
| F17 | Zero levels set the full-sky indices | S3 | Tested (E014). Offsets are still not in the model (D013) |
| F18, F19, F21, F23, F27, F28, F30, F31, F33, F37 | Hierarchical amplitude prior, template error, multipliers | G3, A8 | Tested on simulations, per-pixel path, intensity |
| F20, F25 | NUTS on the index field | G7 | Efficiency. The warm-up bookkeeping of G7 is not covered |
| F22 | Amplitude and index coupling in the BeyondPlanck polarized setup | S2 | Neighbour only. It is about coupling, not about one index for I and QU |
| F24, F29, F32, F38 | Spatial prior and resolution of the index field | S5 | Tested on simulations |
| F26 | A finite CMB prior leaks the CMB realization | A8 | Tested. Decision D021 |
| F34, F35, F36 | Dust temperature | S1 | Tested for the dust temperature only |
| Section 5 | Coverage against Commander3 | S3, D7, A3, S1, S8 | Names gaps without testing them. Two of its statements are out of date (B12, item 5) |

### 8.2 By audit item

| Item | Class | Reference |
|---|---|---|
| G1 | (a) | F9 |
| G2 | (a) for the invalidity, (c) for the latent default | F1, F2, B8 |
| G3 | (a) | F3, F4, F12, F15, F16, F26, D003 |
| G4 | (c) | The Commander3 analogue is noted in an E009 review |
| G5 | (a) | F5, F14 |
| G6 | collects A1, T2, T3, T4, G7 | none |
| G7 | (b) for warm-up and metric, (c) for the shared start and missing burn-in bookkeeping | D012, D015 |
| G8 | (c) | B16 |
| G9 | (b) for existence and cost, (c) for the missing between-chain diagnostic | E020, E021 |
| D1 | (c) for the QU covariance, (b) for the level of correlated noise | Q7, E022 |
| D2 | (a) | F8. Q7 open |
| D3 | (c) | Passing mention in E014 |
| D4 | (c) | B2 |
| D5 | (c) | B1 |
| D6 | (c) | B7 |
| D7 | (b) | Closing remark of F13, E014, section 5, DATA.md, ROADMAP |
| D8 | (c) | none |
| D9 | (b) for its use, (c) for its properties | E014, B14 |
| D10 | (c) as a model element | Q10 covers TOD gain only |
| D11 | (c) | none |
| D12 | (c) | B15 |
| D13 | (c), minor | none |
| D14 | (a) | F7. Q6 open |
| D15 | (c) | Documented in a code comment. B10 |
| S1 | (a) for dust temperature, (b) for curved laws and missing components, (c) for the fixed electron temperature and spinning-dust peak | F34 to F36, section 5 |
| S2 | (c) | F22 is the nearest neighbour |
| S3 | (a) for offsets, (c) for positivity | F17, Q15, D013, section 5 |
| S4 | (c) | B3 |
| S5 | (a) | F24, F29, F32, F38, D019 |
| S6 | (c) | none |
| S7 | (a) for the hierarchical path, (c) for the CG null prior | E030 pass 7, F13 |
| S8 | (b) | Section 5, inconsistently (B12) |
| S9 | (a) | F3, F4, F12, F24, F29 |
| A1 | (c) | Q12 treats the CG as a cost question only |
| A2 | (c) | B5, B6 |
| A3 | (b) for C(l) sampling, (c) for the shared E and B prior, the missing TE and the monopole convention | Section 5, Q11 |
| A4 | (c) | none |
| A5 | (c) | none |
| A6 | (c), minor | none |
| A7 | (c) | B4 |
| A8 | (a) | F26, D021 |
| T1 | (b) for the gain and sky correlation, (c) for the specifics | Q10, B11 |
| T2 | (c) | none |
| T3 | (c) | none |
| T4 | (c) | none |
| T5 | (c) as a code fact | Documented in a code comment. B9 |
| T6 | (c) | none |
| T7 | (c) | none |
| T8 | (c) as a code fact | ROADMAP names bandpass mismatch |

### 8.3 What the cross-check shows

Of the 48 assumptions (G6 only collects others), 8 are tested, 4 are tested in part, 8 are named
but untested in whole or in part, and 28 are not mentioned.

About 25 of the 38 findings concern the amplitude prior, the prior on the index field or the
dust temperature, all in the per-pixel intensity path (G3, A8, S5, S1). That part of the model is
now well mapped. No finding concerns the ingest steps D3 to D6 and D8 to D11, anything polarized
(D1, D4, D6, S2), the CG path (A1 to A6) or the models on the TOD side (T1 to T8). The record is
deep where the work has been and empty next to it, and the empty parts are the ones every
real-data run passes through first.


## 9. The (c) list

Assumptions that no file in `nextgen/` mentions, subject to the search limits in section 12.
Grouped by where they act, and ordered inside each group by how much I think they matter. Bug
numbers refer to section 7.

**Ingest of maps (every maps-only run on real data passes through this)**

1. D5 (B1). Maps are degraded before they are smoothed: band-dependent aliasing, an unrecorded
   pixel window, and a variance rule that does not match the average taken.
2. D9 (B14). The model-error floor is added per native pixel and averages down like white noise,
   about 26 times at E014 settings; its scale is set by the Galactic plane and the zero level.
   The Haslam error budget goes the same way.
3. D3. Noise-weighted smoothing gives each band its own, position-dependent effective beam,
   while the per-pixel model assumes identical beams.
4. D10. A map read from file has no calibration factor; a relative calibration error dg moves a
   two-band index by dg / ln(nu2/nu1).
5. D4 (B2). Q and U are smoothed as two scalar maps.
6. D11. There is no analysis mask in the amplitude solve; a pixel is excluded only by a zero
   inverse variance.
7. D8. A map without a unit or ordering keyword is taken as uK_CMB and RING, with a warning.
8. D6 (B7). Beams are symmetric Gaussians without a pixel window, and the two polarization beam
   routines use different windows.

**CG amplitude path (the only path at full resolution)**

9. A1. A CG solve stopped at a fixed 100 iterations from a zero start is accepted as a sample,
   with `err_tol: 0.0` as the documented practice.
10. A2 (B5, B6), with S7. "No C(l) prior" is N(0, 1) per alm in (uK_RJ at `nu_ref`)^2, so the
    posterior depends on the unit and on `nu_ref`, the CG path never has the flat prior that F3,
    F12 and F15 analyse, and the apodization is bypassed.
11. A5. The preconditioner assumes isotropic noise, constant mixing and Gaussian beams; the
    alternatives in the documentation crash.
12. A3. One fixed power-law prior for E and B, no TE, and C_0 set from C_1.
13. A4. Multipoles above every band limit are drawn from the prior alone, after a logged error.
14. A6. Harmonic analysis by adjoint synthesis with a scalar normalization.

**Sky model**

15. S4 (B3). `get_sky` evaluates a component with an index field at its mean index, for the
    chi-square diagnostic, fixed-component residuals, the conditional likelihood and the TOD
    side.
16. S2. One index is shared by intensity and polarization, and the hierarchical amplitude prior
    exists for intensity only.
17. S1 in part. The electron temperature and the spinning-dust peak frequency are fixed at
    parameter-file values.
18. S3 in part. No positivity for any amplitude.
19. S6. Point sources are a fixed catalogue, with an open FIXME in the CG right-hand side.

**Per-pixel kernel and precision**

20. A7 (B4). The amplitude draw has no guard for pixels with fewer bands than components.
21. D12 (B15). Data, inverse variance and model sky are single precision by default, and the
    conditional MH likelihood sums chi-square in single precision.
22. D13. The degrees of freedom of the chi-square diagnostic are the observed band pixels.

**Gibbs bookkeeping**

23. G4. The parameter mask enters the index step only; plane signal reaches kept pixels through
    the smoothing and the amplitude solve.
24. G8 (B16). Every run uses seed 1995 unless told otherwise, and a run started with
    `init_chain` reuses the random streams of the run it starts from.
25. G7 in part. Both chains start from the same state, there is no burn-in bookkeeping, and NUTS
    warm-up draws are written like samples.
26. G9 in part. Two chains are computed, one is written, and no between-chain diagnostic exists.
27. G2 in part (B8). The likelihood lookup defaults to the invalid legacy move.

**Polarization noise**

28. D1 in part. The QU covariance within a pixel is dropped.
29. T6. Only the inverse diagonal of the white-noise weight matrix reaches CompSep.

**TOD side**

30. T3. The default correlated-noise draw is the stationary solution on gap-filled data.
31. T2. The white-noise level is a point estimate, although the switch is called
    `sample_sigma0`.
32. T4. One 1/f noise model, sampled on a 150-point grid; five other models are stubs.
33. T1 in part (B11). All gain steps are off by default, the temporal-gain prior is hard-coded,
    and its CG is not checked.
34. T7. The CG mapmaker has a white noise model and is accepted at its iteration cap.
35. T5 (B9), T8, D15 (B10). The transfer-function bug, the single sky per band on the TOD side
    and the beam mismatch with per-pixel amplitudes are documented in code comments only.


## 10. Rankings

These are my own judgement from reading the code and the records. The order weighs how large
the effect is likely to be against how sure I am that it is there. Sizes quoted without an F
number are hand estimates.

### 10.1 Calibrated spectral-index posteriors on real WMAP and Planck maps

1. **S3. Zero levels and per-band offsets are not in the model.** F17 measured index shifts of
   0.03 to 0.4 from zero-level choices, hundreds of posterior standard deviations, and decision
   D013 is not implemented.
2. **D2, D9 (B14), T6. The error model that reaches the likelihood is too tight.** Smoothed
   noise is treated as independent (F8: 1.9 and 3.7 standard deviations; up to 3.3 in width by
   my estimate) and the model-error floor averages down 13 to 26 times, so the widths fail
   before any modelling question is asked.
3. **D7. Delta bandpasses and monochromatic unit conversion.** E014 records a factor 1.58 at
   857 GHz and 2.5 percent at 353 GHz. The dust index takes this directly, and if B14 holds the
   HFI bands carry it at full weight.
4. **D10. No calibration factor per band.** Percent-level calibration errors move a two-band
   index by 0.03 to 0.1, and nothing in the model can absorb them.
5. **D5 (B1). Degradation before smoothing.** Band-dependent aliasing of a few uK by my
   estimate, which is the size of the per-pixel rms in the low-frequency bands. Cheap to test
   and to fix.
6. **D3. Band-dependent effective beams.** Wherever the noise level varies across the smoothing
   kernel the bands stop sharing a beam, which the per-pixel model requires. Size unknown.
7. **G3 and A8, with A2 and S7 in the CG path. The amplitude prior and its volume term.** F3,
   F12 and F26 show that the prior decides the index posterior. The CG path has a third prior,
   undocumented, that depends on `nu_ref`.
8. **S1. Fixed nuisance SED parameters and missing components.** The dust temperature is tested
   (F34 to F36). The electron temperature, the spinning-dust peak and curvature are not, and
   E031 finds that real intensity synchrotron is limited by the component model.
9. **S4 (B3) with B15. Mean-index `get_sky` and single-precision chi-square.** Both act on the
   chi-square diagnostic and on the `conditional` likelihood, which is the one the CG path gets.
10. **G4, G7, G8, G9. Mask coupling and chain bookkeeping.** Plane signal reaches the kept
    pixels, warm-up draws sit in the chain, one chain is written, and there is no convergence
    diagnostic.

### 10.2 The large-scale polarized CMB

1. **D1 with T6. No QU covariance and no correlated-noise modes in the map noise model.**
   Low-multipole polarization noise in WMAP and LFI is dominated by correlated modes, and the
   maps reach CompSep with a white-noise diagonal only.
2. **A1 with A5. A truncated CG accepted as a sample.** The slow directions of a preconditioner
   that assumes isotropic noise are the large scales under uneven coverage, which is where the
   polarized CMB signal is.
3. **A2 with A3. The null prior is a white unit prior on E and B.** It is about 100 times the
   true EE power at l = 5, and there is no physical C(l), no TE and no C(l) sampling.
4. **T3, T2, T4. The noise steps on the TOD side.** The default correlated-noise draw is not
   the masked conditional, sigma0 is a point estimate, and only a 1/f model exists. Residual
   1/f noise ends up in large-scale polarization.
5. **T1. Gain.** Off by default, with a hard-coded temporal prior and an unchecked CG. Gain
   errors turn the dipole and the Galaxy into large-scale polarization.
6. **D7 and T8. No bandpass-mismatch model.** One sky per band on the TOD side, the
   per-detector bandpass shift is an open TODO, and polarization angles and efficiencies are
   fixed inputs, so intensity leaks into polarization with nothing to absorb it.
7. **D4 (B2). Q and U smoothed as scalars at ingest.** Percent-level depolarization and E to B
   mixing for 5 to 10 degree kernels away from the equator.
8. **S2. One index for intensity and polarization.** The polarized synchrotron index is tied to
   the intensity one, which free-free, spinning dust and zero levels pull, and the hierarchical
   prior does not exist for QU.
9. **G9, G7, G8. One written chain, a shared start, no diagnostics.** A low-multipole
   polarization result needs demonstrated convergence more than any other product.
10. **D6 (B7), D12, T5 (B9). Window conventions, single precision, the transfer-function bug.**
    Each is at the percent level or below, which is the level of the signal.

### 10.3 What the two lists have in common

Neither list is led by the sampler. Both are led by the error model of the maps (D1, D2, D9,
T6) and by systematic terms that the data model has no parameter for (S3, D7, D10, T8). The
index samplers validated in F12 to F38 are exact for a model whose inputs these items change.


## 11. Cheap first tests

Ten tests, ordered by what they would settle per hour. None needs the cluster. None was run for
this audit, and all should wait until the production chain has released the cores, or run on a
single thread.

1. **Ingest order (D5, B1).** A CMB-only simulation at nside 512 with two native beams (for
   example 5 and 50 arcmin), brought to nside 64 and 120 arcmin in both orders: degrade then
   smooth (the code), and smooth then degrade. The difference map between the two bands, in
   each order, is the band-dependent aliasing. Minutes with healpy.
2. **The floor and the rms the chain used (D9, B14).** One call to `read_data_map_from_file`
   per band, as `nextgen/scripts/q1_grid_check.py:42` already does, with and without
   `add_signal_fraction_to_rms`. Print percentiles of 1/sqrt(map_inv_var) at high latitude and
   set them beside the E014 column (3.3, 3.1, 5.1, 13.4 and 23 uK_RJ at 30, 70, 100, 353 and
   857 GHz). My reading predicts about 1.0, 0.9, 0.4, 0.5 and 0.9.
3. **Index field in `get_sky` (S4, B3).** For a saved sample of an index-field run, compare the
   chi-square from `get_sky_at_nu` with the one from a per-pixel SED built from `beta_map`.
4. **Null prior and `nu_ref` (A2, S7, B5).** A dense solve at nside 8 or 16 in numpy: one
   component, two values of `nu_ref`, flat prior against S = identity. Include the reference
   file's case of synchrotron at 0.408 GHz. Then one call to `P_Cl_prior` with a null amplitude
   and an apodization set (B6).
5. **CG truncation (A1).** The same right-hand side solved with 100 iterations and to
   convergence, at a resolution where a dense solution is also available. Report the relative
   difference per multipole and per component, and the recorded `CG_residuals` of an existing
   real-data run if one is at hand.
6. **Scalar against spin-2 smoothing (D4, B2).** An E-only QU map smoothed by the ingest
   routine and by `hp.smoothing` with `pol=True`. Report the B power and the difference against
   latitude, for 2, 5 and 10 degree kernels.
7. **Under-determined pixels (A7, B4).** A toy call to `cpp_compsep.solve_perpix` with one
   band's inverse variance set to zero in some pixels, so that fewer bands than components
   remain there. Check what the kernel returns and whether the likelihood kernels skip the same
   pixels.
8. **Single-precision chi-square (B15).** `local_loglike` for one real band in float32 and in
   float64, at nside 64 with the mask and at nside 512 or above without it. Then the difference
   of two such values for an index step of one conditional standard deviation.
9. **Noise correlation after smoothing (D2, Q7).** White noise at native resolution through the
   ingest path, many realizations. Compare the pixel covariance with the diagonal the code
   stores, and the posterior width of a scalar index with its scatter over realizations.
10. **Seed reuse under `init_chain` (G8, B16).** A small maps-only run A of 20 iterations, then
    run B with `init_chain` pointing at A and the same seed. Compare sample k of B with sample
    k of A; I expect them to converge onto each other. Repeat with a changed seed.

Tests 1, 2 and 9 together say whether the maps that reach the likelihood obey the model that F16
validated the sampler for. I would run those three first.


## 12. Limits of this audit

1. Nothing was run. A production chain held every core, so this is a reading audit. No possible
   bug in section 7 has been reproduced, and each stays a hypothesis until one of the tests
   above confirms it.
2. The TOD side was read briefly, as asked: the gain, noise and mapmaking modules and the
   handoff to CompSep. The experiment-specific readers under `file_io/experiments/` were only
   spot checked, and sidelobes, jumps and glitches were not audited.
3. Parameter files: `params/param_explanation.yml` and the upstream `Planck_mapsonly` example
   with its three band files and its index-sampling variant were read. Other parameter files
   were only searched.
4. The (c) class rests on keyword searches over the markdown files under `nextgen/`, without
   `reviews/` and `research/`. E031 and several long experiment records were read in part. An
   item can be mentioned somewhere under words I did not search for.
5. Sizes without an F number are hand estimates with small-angle and white-noise
   approximations. The aliasing estimate in D5 is uncertain by a factor of a few. The Haslam
   rms assumed in A2 spans 0.8 to 5 K. The reference values used in B14 (mean sky brightness at
   353 and 857 GHz, Planck white-noise levels) are from memory.
6. Statements about the published analyses (T4, S2, T8, and the remark on correlated noise in
   D1) are from memory and flagged as such. They must be checked against the papers before
   they are cited anywhere.
7. Line numbers are for commit `52061ad`. Files not tracked in that commit were not used as
   evidence, apart from this report.
8. The second report in this directory, `2026-11-09_referee_report.md`, was written in parallel
   by another agent. I did not read it, so overlaps and disagreements between the two are
   unreconciled.


## 13. Verification status (added by the main session, 2026-11-09 23:01)

The audit ran nothing. Items checked since, by the main session:

| Item | Status | Where |
|---|---|---|
| B14 / D9, the model-error floor averages down | Confirmed on the real E014 maps: 27.6 times below nominal; the audit's point predictions hold to 10 percent | E034, F41 |
| B1 / D5, degrading before smoothing | Confirmed on a simulated CMB: 0.6 to 1.7 uK_CMB between bands, zero in the other order | E034, F41 |
| G9 / G7, two chains computed, one written, shared start | Confirmed by reading `cli.py`, `chain_writer.py` and `index_field_sampler.py`; E032 now writes both | F42 |
| B3 / S4, `get_sky` at the mean index | Confirmed by reading `sky/diffuse_components.py`; no recorded result depends on it | F43 |
| B8 / G2, latent default to the legacy move | Confirmed by reading `compsep/processing.py`; the comment there states the fallback on purpose. Not exercised by any run | none |
| B5 / A2 and B6, the null C(l) prior | Confirmed by reading `sky/diffuse_components.py` (`P_Cl_prior` returns ones before the apodization) and `compsep/cg_solver.py` (the identity term and the prior fluctuation are always present): a null prior is N(0, 1) per alm in (uK_RJ at nu_ref)^2. Its size against the data, including the 408 MHz case, is not yet computed | Q22 |
| B2, B4, B7, B15, B16, A1 | Not yet checked. Tests 5 to 8 and 10 of section 11 are queued behind the E032 chains | none |

The sizes quoted for B14 in section 7 were hand estimates; the measured ones are in E034.

