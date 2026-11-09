# Referee report: the nextgen Commander4 program and the Paper A outline

- **Written:** 2026-11-09 21:33 CET (clock time at the start of writing)
- **Repository state reviewed:** branch `nextgen`. Reading started at commit c475106. The authors
  committed 52061ad at 21:39:51, during the review; it adds the chain 3 raw outputs
  `nextgen/experiments/E030/chain3_*.txt` (written 20:58 to 21:00, not yet read into the E030
  record, FINDINGS.md or the outline) and experiment E031. I read
  `nextgen/experiments/E031_noise_regime_of_the_ladder_vs_real_maps.md` at 21:18 (Result and
  Conclusion empty) and again at 21:39 (rerun of 21:28, conclusion recorded 21:30), with
  `nextgen/experiments/E031/fisher_report.txt` and `fisher_report_run1.txt`. At 21:50 FINDINGS.md
  (last modified 12:21), the outline, OPEN_QUESTIONS.md and STATUS.md (14:06) and DECISIONS.md
  (8 November) were unchanged, so every statement below about "the documents" refers to those
  versions.
- **Role:** external referee (Bayesian CMB component separation, Gibbs sampling, MCMC validation).
  Every conclusion in the documents was treated as a claim and checked against the experiment record
  it cites.
- **Constraints honoured:** read-only, nothing executed (grep, sed and arithmetic only), no git
  operations, no subagents, no outside contact. This file is the only file written.

All paths below are relative to `/Users/trygvels/code/Commander4/` unless they start with `/`.
"Line N" refers to the file as it stood at review time.

Grades: **A** = established against an exact or independent reference, with replication.
**B** = supported, but one sky, one noise realization, one seed or one configuration.
**C** = suggestive only. **D** = could be an artefact of the setup or a bug.

---

## 1. Summary

### 1.1 The ten objections, ranked by severity

**1. The noise regime is mislabelled, and most of Paper A Sect. 5 lives in a regime with no
real-data counterpart.**
`nextgen/scripts/pysm_intensity_sim.py` line 20 and line 159 set the noise to
`RMS16 * (nside / 16) * noise_scale` with RMS16 = 70, 60, 60, 55, 60, 60, 70, 40, 24, 16, 24, 80
uK_RJ. At nside 64 and scale 1 that is 280, 240, 240, 220, 240, 240, 280, 160, 96, 64, 96, 320 uK_RJ
per pixel. The numbers descend from the E006 test harness (`aux/validate_mcmc_likelihood.py`,
`NSIDE = 16`, `RMS = [60, 60, 70, 24, 16, 24, 80]`) and `nextgen/prototypes/indexfield/simulate.py`
lines 23 to 28 says the WMAP-like values "are set here at comparable levels and are part of the
stated setup". They were never derived from an instrument.
E020 (`nextgen/experiments/E020_index_prior_comparison_pysm_d1s1.md` lines 11 and 27 to 28) calls
scale 0.03 "Planck-like white noise" and "within a factor 2 of WMAP/Planck at this resolution".
From D014 on (`nextgen/DECISIONS.md` lines 114, 121, 270 to 272, 339 to 342, 361) scale 1, which is
33 times noisier, is "realistic noise"; F29 to F38, the outline (line 135, "At Planck-like depth"),
E027 (line 212 on) and Q16 (`nextgen/OPEN_QUESTIONS.md` line 340, "realistic noise (noise scale
1)") follow. No record before E031 gives a reason for the change of label.
Cross-check from the program's own real-data record, independent of anything I remember:
`nextgen/experiments/E014_q1_real_maps_lowres.md` line 314 gives a median intensity rms at nside 512
of 11.9 / 4.7 / 6.5 / 21.9 uK_CMB at 100 / 143 / 217 / 353 GHz. Divided by 8 for nside 64 and
converted to RJ (factors 0.777, 0.605, 0.334, 0.077) this is 1.16, 0.36, 0.27, 0.21 uK_RJ. Scale 1
is 83, 180, 355 and 1500 times higher; scale 0.03 is 2.5, 5.3, 10.7 and 46 times higher. Line 76 of
the same record says "noise ~1 uK per pixel" at nside 64.
E031, started by the authors during this review, confirms it from the rms maps on disk
(`nextgen/experiments/E031/fisher_report.txt`, Table 1): scale 1 over real white noise is 24 to 90
for WMAP and LFI, 88 at 100 GHz, 176 and 337 at 143 and 217, 1090 at 353 (NPIPE A split). The median
synchrotron S/N per nside-64 pixel at 23 GHz is 0.31 in the ladder against 16.0 in the real K band.
Its conclusion, recorded at 21:30, withdraws the label for intensity. That is the right call, and it
is to the authors' credit that it was made within the hour. When this report was written the
withdrawal existed in the E031 record and its registry row only: FINDINGS.md (F26 to F38), the
outline (Sect. 5), STATUS.md, DECISIONS.md (D014, D021, D024, D025) and Q16 still carried the label,
and the registry row cited a finding F39 that was not in FINDINGS.md.
Consequences: (a) F26's leak matters only because the per-pixel CMB noise is 58 uK; (b) "every
region is prior-dominated" (`nextgen/experiments/E023_hier_robustness_pysm_exact_target.md` line 98)
and with it F29, F32 and F38 are properties of this noise level: with the same three-component
model at real white noise the profiled region sd of beta_s has a median of 0.078 with 70 percent of
regions below 0.1 (E031 Table 2), against 5.99 and none at scale 1; (c) the ladder is not a clean
polarization proxy either (5.99 against 0.77 for a WMAP K plus LFI polarization set and 2.18 without
K band, E031 Table 3); no record before E031 argues that it was meant as one, and E031's own
criterion for it (item 6, within a factor 5 of the K-band set) failed; (d) at the level that is
close to real WMAP and LFI depth (scale 0.03, ratios 0.7 to 2.7) the hier prior ran away under a
data-built template (F21, E020), and the multiplier remedy was never run as a chain there (no such
chain directory exists under `/Users/trygvels/data/e020_d1s1/`; the evidence at 0.03 is the
profiled-multiplier scans behind F28, which F31 says must be redone); (e) with free-free and AME in
the model the real-noise region sd is 1.36 (0.10 with the 408 MHz map, E031 Table 2), so the
binding constraint on real intensity data is the component model, which a three-component ladder
cannot exercise.

**2. Thesis claim (3) is not supported as worded.**
"Calibrated index posteriors down to index-pixel S/N below 1" (outline thesis, F16) holds for data
drawn from the model's own hyperprior, exact-law templates and one sampled index with the dust SED
fixed at truth (`nextgen/experiments/E016_hierarchical_amplitude_prior.md` line 22). Under a 30
percent template error, which the program treats as the realistic case, four chains on the same
data give synchrotron field means of -2.723 (E026 cs), -2.701 (E030 chain 1), -2.722 (chain 2) and
-3.164 (chain 3) against a truth of -3.000, each with a posterior sd near 0.05, and chain 3 has
region-level coverage |z| < 1 of 0.27 to 0.33
(`nextgen/experiments/E030/chain3_sync_level.txt`). The E030 record concludes (line 599) that the
width "is a bias-variance knob with no calibrated setting under a 30 percent template". The outline
(lines 257 to 258) still calls the hierarchical prior "the calibrated choice for spatially varying
indices".

**3. The remedy already written into the outline failed its own pre-registered test, and the bias
changed sign instead of vanishing.**
Outline lines 200 to 209 carry a "Statement for the paper: pin the amplitude of a template-informed
component at the band that carries its index information", with the chain-level confirmation
"pending (E030 chain 3)". Chain 3 has finished. Against its pre-registration
(`nextgen/experiments/E030_hyperprior_width_vs_template_law.md` line 540 on): field mean -3.164 +-
0.049 (asked: within 0.10 of -3.00); faint quartiles -0.18 against the region reference (asked:
within 0.10); hundred brightest regions -0.126 (asked: within 0.05); multiplier map Q2 correlation
0.45 (asked: above 0.6) and rms 0.27 (asked: below 0.20). Full table in Sect. 7. The displacement
went from +0.30 to -0.16. The pass 7 scan had shown this sign: its faint-quartile optima were -0.18,
-0.24 and -0.30, read as "flat" because the gain was 1.0 to 1.5 nats (E030 lines 497 to 501). A sign
change is also what the program's own F29 mechanism gives when the reference crosses the
information-weighted mean log frequency: the Occam term of the marginalized multipliers has a beta
slope equal to minus the weighted mean of ln(nu_b / nu_ref), which for this band set is +0.14 nats
per region per unit beta at nu_ref = 30 GHz, -0.12 at 23 GHz and zero near 26 GHz (my arithmetic,
single-component weights, Sect. 7.2; not verified with the scan code, and it does not by itself
give the magnitude). Whatever the mechanism, the same data, template error and hyperprior width give
field means of -2.70 and -3.16 depending on where the amplitude is referenced, nine posterior sd
apart. The hier target with free per-region multipliers is therefore not invariant to the reference
frequency, which is the defect that thesis claim (2) holds against flat amplitude priors.

**4. The production sampler is not validated at the configuration whose results are reported
(thesis claim 4).**
The end-to-end check in `nextgen/experiments/E018_c4_index_field_hier_validation.md` is one
replicate at nside_beta 2 (K 47), one index field with the dust index fixed, against the authors'
own prototype (lines 124 to 125 and 146). The L1 field test passes only with "the synchrotron 30
times" brighter than the template; at template brightness the two samplers disagree on log sigma
(-1.83 +- 0.02 against a reference drifting from -1.31 to -1.52), and the discrepancy is assigned to
the reference without a longer run (lines 83 to 98). Production chains (nside_beta 16, two index
fields, 2 x 3072 multipliers, about 12,300 dimensions) are single chains of 700 post-burn draws with
index-field hyperparameter ESS as low as 5 to 8 (`nextgen/experiments/E026/hyper_ess_noncentred.txt`,
cs block: 5.7 and 4.9 for the dust log sigma and u; `nextgen/experiments/E029/hyper_ess.txt`: 8.3
and 7.4; chain 3: 21 to 224), no R-hat anywhere in the program (`nextgen/VALIDATION.md` has no
multi-chain criterion), no independent sampler ("not available", `nextgen/DECISIONS.md` line 378),
and no SBC of the Commander4 implementation. The chain 3 output directory
(`/Users/trygvels/data/e023_pysm_n1/chains_gmrf16_cs_mumap_nc_Tmap_mu03_ref23/chains_compsep/`)
holds 1200 `chain01` files and nothing else.

**5. Single realizations carry quartile-level conclusions.**
Every scale 1 result uses one sky (PySM d1s1c1), one noise realization (seed 2027) and one
template-error field (seed 11; a second seed in scans only):
`nextgen/experiments/E026_stage2_chains_cs_template_mumap.md` line 336,
`nextgen/experiments/E023_hier_robustness_pysm_exact_target.md` line 406, the scope line of F36.
Offsets of 0.01 to 0.06 per quartile are given mechanisms without a between-realization variance,
although E030 pass 4 found the bright-synchrotron shift to be "+0.07, -0.06, +0.01" across error
realizations (F37), and the follow-up "A second noise realization of the PySM noise-1 sky" (E026
line 349) was never run.

**6. Pre-registration is honoured in the recording and eroded in the reading.**
About twenty pre-registered criteria failed or were re-posed between E016 and E030 (table in Sect.
4.4). Misses are recorded as misses, which is to the program's credit. But D027
(`nextgen/DECISIONS.md` lines 410 to 430) adds, after the misses, a rule by which a failed
per-quartile tolerance in a "prior-dominated" quartile is read against a shrinkage prediction
instead. Chain 3's criterion (vi) uses the GLS CMB residual that F34 had declared "not an
index-quality criterion"; chain 3 passes it at 90.9 uK, 2.2 uK below the true-SED "floor" of 93.18,
with a synchrotron field biased by -0.16 everywhere.

**7. Oracle inputs, hyperpriors tuned to the simulation truth, and no trivial baseline.**
The exact templates are the truth; the T map is the truth at nside 16; the multiplier hyperprior is
centred on 1 and its sd was chosen each time knowing the simulated template error (0.5, then 0.1
for an exact template, then 0.3 "matched" to the 30 percent field: F33, F36, F37); the index prior
equals the generator in `nextgen/experiments/E017_scalar_index_coverage.md` line 25. F33 and F36 are
existence proofs under oracle inputs. Sect. 5 reports pixel RMSE without a baseline. A constant
field at the true mean has an RMSE of 0.063 (synchrotron) and 0.037 (dust), the truth sds quoted in
F32, and the authors' own scalar-index chain on the same data (-2.944,
`nextgen/experiments/E030/chain3_sync_level.txt`) corresponds to 0.084 for synchrotron (my
arithmetic). In Table 3 of `nextgen/papers/paperA_tables.md` every chain with the 30 percent
template is worse than that scalar by a factor 2 to 5 (0.335, 0.426, 0.327; chain 3 0.177); the
exact-template fields at sd 0.5 are worse too (0.133 to 0.147); only the oracle-tuned sd 0.1 chains
reach the constant-at-truth level in synchrotron (0.061 to 0.066), and no chain beats it in dust
(best 0.038). Chain 3's 0.177 is what a constant at -3.164 gives (0.175). So "RMSE" in Sect. 5
measures the level offset and prior-induced scatter, and F32's resolution rule (nside 4, 8, 16 at
0.143, 0.135, 0.147) ranks three fields that are all worse than one scalar index.

**8. Thesis claim (1) overreaches beyond Commander4.**
"The acceptance step used in production code" cites "E009/E010 for Commander3", which found the
BeyondPlanck-era samplers free of the F1 structure (F10). The Commander3 code that has it is a 2024
module, established by code reading only and never run. No published result is shown to be
affected. The k^(1/4) "prediction" of 18.5x (F1; E014 lines 191 to 192) divides 1.9e5 amplitude
degrees of freedom by 600 pixels; with amplitude degrees of freedom on both sides (three components
per pixel in the E005 toy, `tests/test_mcmc_likelihood.py`) the prediction is 14.1 to 14.5 against a
measured 16.5 to 17.1. On real maps the legacy chain's mean sits about 3 ridge sd from the exact
posterior (0.001 in beta_s, not significant at ESS 74; 0.00013 in beta_d), against a zero-level
systematic of 0.035 to 0.47.

**9. Thesis claim (2): the effect is known, the paper cites nothing, and its relevance to
production chains is explicitly not established.**
The Jeffreys remedy for prior-volume bias is credited to Eriksen et al. 2008 in the program's own
`nextgen/LITERATURE.md` (line 35) and in F4; the outline contains no citation at all. F15 says the
result "does not transfer to Commander3/4 production chains as they stand".
`nextgen/experiments/E022_bp_xv_polarized_index_gibbs_coupling.md` (lines 401 and 454) finds the
flat-prior shift absent in the BeyondPlanck chains, with the mechanism "not identified". On real
maps the shift is quoted as 54 to 56 sd where the sd comes from a model with a reduced chi-square
of 15000 (E014 line 131), so the unit has no sampling meaning.

**10. F17 and Sect. 7 state more than E014 shows.**
"beta_s from -3.34 to -2.86 depending only on the input zero levels" (outline line 248) compares a
17-band set (NPIPE A-split 353, PR3 100) with a 20-band PR4 set that adds 545 and 857 GHz (E014
lines 305 to 316; `nextgen/OPEN_QUESTIONS.md` line 232 says "between the two HFI sets"). The
like-for-like zero-level effect (pass 3 to pass 4, "same 20 maps, only the zero levels changed",
E014 line 388) is +0.035 in beta_s and -0.016 in beta_d, and the Haslam dipole alone +0.011.

### 1.2 Claims that earn grade A

1. F1, core statement: Commander4's `resolve_chisq` acceptance is not a valid MCMC step for the
   stated posterior, and the posterior widths inflate in controlled tests (E002 with E003 as an
   independent reference, E005, E006, E011 on a second machine, E014 on real maps).
2. F2: the `ridge`, `marginal` and `conditional` modes reproduce their exact grid targets (E005,
   E006, E011, the tuned E014 chains at ESS above 400).
3. F15, with the simulation side of F3, F12 and F18: flat amplitude priors give the collapsed
   marginal with the volume term -1/2 ln det A; on template data it biases a scalar index and
   collapses an index field (derivation, exact grids, 12 to 96 replicates).
4. F31: profiling per-region multipliers is not marginalizing them, and the analytic integral is
   right (brute force to 1e-8, quadrature to 5e-12).
5. F23, mechanism only: the hierarchical amplitude prior is exact only for its own template-error
   law (independent marginal code, replays of both generative stories, a 20 sd effect).
6. F26, algebra only: a finite per-pixel CMB prior leaks the CMB into the index likelihood with
   fraction sigma_c^2 / (tau^2 + sigma_c^2) (sign-flip, roll and shuffle controls).

Nothing in Sect. 5 of the outline beyond F31 reaches grade A, and nothing about real data does.

### 1.3 Paper-level judgement in one paragraph

There is one coherent short paper in this material, and it is not the one outlined. Grade A and B
evidence supports a correctness-and-priors paper: the legacy acceptance step is invalid and the fix
is exact; the amplitude prior (flat, profile, hierarchical) moves spectral indices by more than
their posterior width, with a derivation and exact-grid demonstrations on simulations and on real
low-resolution maps; a collapsed likelihood with a hierarchical template prior is calibrated inside
its model class and fails outside it, through mechanisms that are established in part (the template
law, F23; profiling against marginalizing, F31) and not yet predictive (chain 3). Sect. 5 as it
stands (E023 to E030, 149 of the outline's 275 lines) should leave the paper. It is a
single-realization study at a noise level 24 to 1090 times above the real maps, in which no chain
with a realistic template beats a single scalar index, three scan passes were invalid, the
recommended hyperprior width changed three times, and the last remedy reversed the sign of the bias
instead of removing it. That sequence has been a detour with diminishing returns since E029. The
authors' own E031 conclusion, recorded while this report was being written, reaches the same view
of the regime and proposes to replace further intensity template work with a polarization
configuration. Before more chains are run the program should carry that withdrawal through the
documents, add replication and multi-chain diagnostics, put a scalar baseline in every table, and
decide whether Paper A is a correctness note or a study of template priors. On the present evidence
it cannot be both.

---


## 2. The headline claims of the Paper A thesis

The thesis (`nextgen/papers/paperA_methods_outline.md` lines 12 to 27) makes four numbered claims
and one by-product. I split each into the part the evidence carries and the part it does not. Each
entry gives the grade, the evidence I checked, the strongest objection and the smallest experiment
that would settle it.

### 2.1 Claim (1): the production acceptance step is not exact

- **(1a) Commander4's `resolve_chisq` move is invalid: grade A.**
  Evidence checked: registry rows E002, E003, E005, E006 (`nextgen/experiments/registry.md` lines 10
  to 17); the unit tests `tests/test_mcmc_likelihood.py`; the harness
  `aux/validate_mcmc_likelihood.py`; the second-machine reproduction
  `nextgen/experiments/E011_mac_mini_repro/`; the real-map chains of
  `nextgen/experiments/E014_q1_real_maps_lowres.md` (lines 182 to 196: width ratios 17.1 and 16.5
  against the exact ridge grid, with the valid ridge chain at 1.01 and 1.18 as the control).
  Strongest objection: the invalidity is a derivation and it is right. The demonstrated consequence
  is a posterior that is too wide, for a full-sky scalar index whose real uncertainty is set by zero
  levels (objection 10). No case is shown where the invalid step changes a scientific conclusion.
  E001 to E008 have registry rows but no experiment records, so the simulation evidence for F1 rests
  on the tests and the harness. That is acceptable because they can be rerun, but it is not the
  "question, setup, command, result, conclusion" record that CLAUDE.md principle 5 asks for.
  Smallest settling experiment: none is needed for validity. For relevance: one low S/N or spatially
  varying configuration where the legacy move and an exact move are both run and the index map
  differs by more than its error.
- **(1b) The k^(1/4) width law: grade B.**
  Evidence: three points, 3.2 at 450 amplitude degrees of freedom, 4.4 at 1800 (E005), 16.5 to 17.1
  at 1.9e5 (E014).
  Objection: the E014 check mixes pixels and degrees of freedom (objection 8). With consistent units
  the law under-predicts the real-map inflation by 14 to 21 percent, which is fine for a scaling law
  but is not what the record says ("follows the k^(1/4) scaling ... = 18.5", E014 lines 191 to 192).
  The scaling itself is plausible on general grounds and survives the reduced chi-square of 15000:
  the chi-square is quadratic in the amplitudes, so at freshly drawn amplitudes it fluctuates as a
  chi-square with k degrees of freedom whatever the misfit. The constant depends on the proposal
  tuning.
  Settle: state the law with its derivation and a table in degrees of freedom; add one point between
  1800 and 1.9e5 from the harness.
- **(1c) "used in production code", citing "E009/E010 for Commander3": grade C as worded.**
  Evidence: F10 says the BeyondPlanck-era samplers do not have the F1 structure; the Commander3
  module that has it is from 2024 and was read, not run. F11 (the bandpass-shift acceptance on rank
  0's chi-square) and the sign error in the Commander3 local `marginal` lnL (F10) are also code
  reading.
  Objection: a paper cannot assert errors in another group's code from reading it. Each needs a run
  that shows a wrong stationary distribution, or it should be raised with the developers and left
  out. No published number is shown to be affected.
  Settle: a toy problem driven through the actual Commander3 routine and compared with a grid
  posterior. Contacting the developers needs Trygve's decision (CLAUDE.md, "Ask Trygve before").

### 2.2 Claim (2): the flat-prior volume term

- **(2a) The mathematics and the simulated demonstrations: grade A.**
  Evidence: F15 (derivation, independently reviewed), E012 `study_v1` (108 cells), E017 (96 profile
  and 24 other replicates per noise scale against exact one-dimensional grids).
  Objection: the size and the direction of the effect depend on the amplitude reference frequency.
  F3 reports a rise of about 1800 in log density with dust normalized at 545 GHz, and F13 then
  corrects the DR1 reference to 70 GHz. A flat prior on the amplitude at one reference frequency is
  a beta-dependent prior at any other, so "the flat-prior posterior" is a family and the paper has
  to say which member each production code uses.
  Settle: one figure of the volume-term shift against reference frequency for the E017
  configuration; the grid code exists.
- **(2b) "moves a full-sky index by tens of posterior standard deviations" on real maps: grade B for
  the shift, grade D for the unit.**
  Evidence: E014, +56.05 +- 0.06 ridge sd in beta_s for pass 2 (lines 264 to 270), +16.1 for pass 3
  and +22.7 for pass 4 (lines 386 to 392).
  Objection: the absolute shift is +0.019 in beta_s for pass 2 (F3) and 0.003 for pass 4 (22.7 times
  0.00014, my arithmetic). The "posterior sd" comes from a fit with a reduced chi-square of 15000
  (E014 line 131), so it is a formal curvature and not an uncertainty; the number of sd changes by a
  factor 3.5 between passes. In absolute terms the prior shift is smaller than the like-for-like
  zero-level shift of 0.035.
  Settle: quote absolute shifts, and repeat on a mask where the model fits (reduced chi-square near
  1) so that an sd means something.
- **(2c) Relevance to published chains: grade C to D.**
  F15 itself says the result "does not transfer to Commander3/4 production chains as they stand".
  E022 finds the shift "absent in BP's chain" with the mechanism "not identified"
  (`nextgen/experiments/E022_bp_xv_polarized_index_gibbs_coupling.md` lines 401 and 454).
  Settle: identify the E022 mechanism, that is, which prior the amplitude block effectively has in
  the BeyondPlanck CG step, before the paper says anything about production analyses.

### 2.3 Claim (3): calibrated index posteriors down to index-pixel S/N below 1

- **Inside the model class: grade B.**
  Evidence: `nextgen/experiments/E016_hierarchical_amplitude_prior.md`, SBC on 600 replicates over
  16 quantities and a coverage study on 108 replicates with replicate-clustered errors.
  Objections. (i) SBC draws the data from the model's own prior. It tests the sampler and the code;
  a correct sampler passes for any prior, so SBC says nothing about whether this prior suits a sky.
  (ii) The first 300 replicates failed the pre-registered rule ("Verdict by the pre-registered rule:
  FAIL", lines 259 to 269). The rule that turns a clean replication into "PASS on 600" was written
  after that failure, though before the replication was read (line 283 on). I accept the reading
  (pooled rank mean 0.520, +1.7 se, my arithmetic agrees with F16), but the F16 heading should say
  that the first run failed. (iii) The coverage study has a mean z of +0.35 +- 0.16 at noise 1 and
  one 58-pixel bin at 0.48 / 0.86 (line 166). (iv) NUTS shows "low E-BFMI (0.03 to 0.17) at every
  noise level" (line 177); the program's own E017 analyzer flags values below 0.3
  (`nextgen/experiments/E017_scalar_index_coverage.md` line 130). (v) One index is sampled; the dust
  SED is fixed at truth (E016 line 22). (vi) The "wrong template" variants (E016 lines 88 and 141)
  distort g while the amplitudes scatter about the undistorted mean, which by the program's own
  later F23 is the benign direction of template error.
- **As a general statement: grade D.**
  The same prior fails under the data-story template error (F23, +0.114 at 20 sd), runs away under
  SED or template misspecification at the noise level near real depth (F21, E020), and has "no
  calibrated setting under a 30 percent template"
  (`nextgen/experiments/E030_hyperprior_width_vs_template_law.md` line 599). Objection 2 gives the
  chain numbers.
  Settle: SBC with data generated under the data story (g = a / r at 30 percent) in the existing
  prototype, and a scope clause in the thesis sentence ("for templates that follow the prior's
  law").
- Stale text: the thesis says "E017 for the full-sky case, pending"; F18 reports it done.

### 2.4 Claim (4): implemented in Commander4, reproduces exact targets end to end

- **At nside_beta 2, one replicate: grade B.**
  Evidence: `nextgen/experiments/E018_c4_index_field_hier_validation.md` L3 (lines 100 to 147): 55
  scalars within 1.9 pooled standard errors of the prototype chain, sd ratios 0.86 to 1.09, and a
  good negative control (the same chain is 19 Monte Carlo standard errors from the profile grid and
  1237 grid sd from the flat-marginal grid).
  Objections. (i) The reference is the authors' own prototype of the same derivation, not an exact
  target, so the two share any error in the target. "Reproduces exact targets" holds for the scalar,
  non-hierarchical modes only. (ii) One replicate, with a resolution of about 0.15 posterior sd in
  means and 15 percent in widths (my estimate: 4 pooled standard errors at ESS 1500). (iii) The L1
  field test at template brightness fails on log sigma and the failure is assigned to the reference
  sampler (lines 83 to 98): the tested chain gives -1.83 +- 0.02 while the centred reference moves
  from -1.31 to -1.52 as it is extended. That is unresolved, not explained; the test was then "run
  where the reference mixes", with the synchrotron 30 times brighter. (iv) The comparison was
  iterated after an analyzer bug (same record), which is normal practice but means agreement was the
  stopping rule.
- **At the production configuration: grade D.**
  No test of the sampler exists at nside_beta 16 with two index fields, per-region multipliers and a
  temperature model (objection 4).
  Settle: four dispersed chains with split R-hat at the production configuration; SBC of the
  Commander4 implementation at nside_beta 2 with 100 replicates; one run at template brightness
  against a reference that mixes.
- Stale text: the thesis says "E018, pending L3" (F19 reports it done), and F19 says the L2 gate "is
  still running" (F16 reports it done).

### 2.5 By-product: the zero-level convention sets the full-sky index (F17)

- **Like-for-like effect: grade B.** Pass 3 to pass 4 on the same 20 maps: +0.035 in beta_s and
  -0.016 in beta_d; the Haslam dipole alone +0.011 (E014 lines 386 to 392).
- **"0.03 to 0.4" and "depending only on the input zero levels": grade C.** The 0.47 compares
  different band sets (objection 10).
  Objection: production analyses fit or sample per-band monopoles, so a full-sky index at fixed
  monopoles is not a quantity that the literature reports (my recollection, not verified in this
  review); "hundreds of posterior sd" uses the unit that objection 9 rejects.
  Settle: add free per-band offsets with one anchor to the ridge grid, which is a linear extension,
  and report the index and its width with the offsets marginalized.

---



## 3. Evidence grades for F1 to F38

Findings are in `nextgen/FINDINGS.md`. Related findings are grouped. "Settle" is the smallest
additional experiment. The regime objection (1) and the single-realization objection (5) apply to
every entry of Sect. 3.8 and are not repeated there.

### 3.1 Sampler correctness: F1, F2, F9, F10, F11

- **F1: A for the core statement, B for the k^(1/4) law.** See Sect. 2.1.
- **F2: A.** Evidence: E005, E006 (7 ranks, within 0.03 sd of the exact ridge target, width ratios
  0.94 and 0.97), E011, the tuned E014 chains. Objection: in E006 the exact ridge target that the
  chain reproduces is itself 1.9 and 3.7 sd from the simulation input (F8), so "validated end to
  end" means "the chain samples its target", not "the target recovers the truth". Settle: the F8
  experiment below.
- **F9: C.** "The two-chain pipelining in `cli.py` is valid" has no experiment ID; it is an
  assertion from reading the code. Settle: cite the code lines, or run two chains with a planted
  cross-chain dependence and show the test catches it.
- **F10: B and C.** B for "the BeyondPlanck-era index samplers do not have the F1 structure" (two
  review records, `nextgen/experiments/E009_c3_nonlin_review/` and
  `E010_c3_polarized_and_bandpass_review/`). C for the sign error in the local `marginal` lnL and
  for the 2024 module having F1: both read, neither run.
- **F11: C.** Code reading only (E015). Settle for F10 and F11: Sect. 2.1 (1c).

### 3.2 Amplitude priors: F3, F4, F8, F12, F13, F15, F18

- **F15: A. F12 and F18 on simulations: A.** Derivation plus exact grids plus 12 to 96 replicates
  per cell. See Sect. 2.2.
- **F3: B.** The real-map shift is real and exact for its target (+0.019, tuned chain against grid
  prediction); the unit is objection 9. "Scientifically useless for full-sky indices" is the 545 GHz
  normalization; at the DR1 reference (F13) the shifts are 3 to 23 sd.
- **F4: C.** "Ridge is the principled default" is a recommendation that rests on a literature credit
  (Stompor et al. 2009, taken from `nextgen/LITERATURE.md`, not verified by me) and on the identity
  ridge = flat-prior marginal times det(A)^(1/2). The finding itself ends "This needs a proper
  coverage study (Q2)", and that study shows the profile failing at low S/N: z of +5 to +13 with
  coverage 0.3 to 0.6 below S/N 3 (F12) and coverage 0.33 / 0.66 at S/N 23 (F18). A conditional
  Jeffreys prior on the amplitudes that depends on beta is itself a prior on beta, and the paper has
  to say so in those words. Settle: merge F4 into F12 and F18 with the S/N range over which the
  profile is calibrated.
- **F8: C.** The exact ridge target sits 3.7 sd from the input on the validation simulation, and the
  finding assigns this to correlated noise after smoothing without a test: Q7 is still open
  (`nextgen/OPEN_QUESTIONS.md` line 162). By CLAUDE.md principle 8 it is a bug until shown
  otherwise. Settle: 20 noise seeds through the harness with and without the smoothing step; the z
  distribution of the exact ridge optimum decides it with grid evaluations only.
- **F13: B.** A correction of the program's own earlier assumption about the DR1 reference
  frequency; credit for recording it. It also shows how sensitive F3 is to that convention.

### 3.3 Efficiency: F5, F6, F14, F20, F25

- **F5: C.** One 300-pixel toy with a registry row and no record. The 9 times and 22 times compare
  collapsed NUTS with Gibbs plus MH in amplitude solves per effective sample, which is a fair unit,
  at one problem size and one tuning of the MH proposal.
- **F6: C, and not paper material.** Implementation timings of one kernel.
- **F14: B for the ratio, C for the consequence.** The conditional to marginal sd ratio of 180 for
  beta_s is an exact-grid result on one configuration (nside 64, 38165 pixels, 17 bands). "Of order
  3e4 iterations per independent" sample is the ratio squared; no autocorrelation time was measured.
- **F20: B minus.** "one seed per cell; 2000 samples; our own adaptation schedule"
  (`nextgen/experiments/E019_sampler_ess_per_gradient.md` line 119). Ratios of 0.3 to 1.1 for the
  worst pixel support "no gain seen in these cells", not "there is no case for a second sampler".
- **F25: C.** The factor 27 is a ratio to an ESS of 9 out of 700 draws, which is not known to better
  than a factor of about 2. The third pre-registered clause failed
  (`nextgen/experiments/E024_dense_hessian_metric_index_nuts.md` lines 148 to 152) and the run
  "started (11:56) with the pre-review code" (line 75). Scope by the finding's own words: nside 8,
  one sky, high S/N.
- Settle, for all five: three seeds per cell, and ESS quoted with a between-seed uncertainty.

### 3.4 Implementation: F7

- **F7: B.** One measurement (reduced chi-square 0.57 to about 794 after an alm round trip). A real
  hazard that deserves a sentence, not a finding in the paper.

### 3.5 The hierarchical prior inside its model class: F16, F18 (hier part), F19

- **F16: B.** See Sect. 2.3.
- **F18, hier rows: B minus.** 24 replicates per noise scale; coverage at 68 percent between 0.50
  and 0.83; eight cells collapsed at the log s floor and were rerun after a guard was added to the
  target (E017 lines 127 to 136; the finding does not say whether the collapse was the sampler or
  the target); the index prior equals the generator (E017 line 25). The closing sentence, "with hier
  or ridge the real-data index would agree to a fraction of a sd", is a prediction: hier was not run
  on the real maps.
- **F19: B at K 47.** See Sect. 2.4.

### 3.6 Real data: F17, F22

- **F17: B for the like-for-like effect, C for the headline.** See Sect. 2.5.
- **F22: three statements, three grades.**
  (a) Surrogate ratio of collapsed to conditional sd, 2.35 (Spur) and 2.11 (Plane): B as a
  computation (exact, four review rounds). But the finding's own control shows what it is a property
  of: with a C_l-sized CMB prior the ratio is 1.02 to 1.04. The Gibbs-coupling penalty in this
  configuration comes from the flat per-pixel CMB prior of the surrogate, and whether the
  BeyondPlanck CG step behaves that way is the open point of the record.
  (b) Diagnostics of the released chains (lag-1 autocorrelation 0.79 to 0.83, 23 to 114 effective
  samples per chain, drift): B, read from public files.
  (c) "The BP XV chain's flat-prior stationary target for the Spur lies at the boundary and the
  released chains are far from it": D as a publishable claim. The record calls it "consistency, not
  confirmation" (E022 line 262), "an indication, not a measurement" (line 289) and "an
  order-of-magnitude discriminator" (line 407), and the finding says "inference, labelled; a
  direction, not a value". The registry row (line 38) quotes two chains and the finding four.
  Settle: rerun the BP XV polarized index step for one region from the released parameter files with
  ten times the chain length. That is a Commander3 run and Trygve's call (E022 line 462).

### 3.7 Misspecification at 0.03 times the ladder noise: F21, F23, F24

These three are at the level that E031 finds close to real WMAP and LFI depth, and they are the part
of the misspecification work I would keep.

- **F21: B plus.** One sky and one realization, but the effect is of order 1e5 nats and its sign is
  reproduced by the F23 replays. It is the most important negative result in the program.
- **F23: A for the mechanism, B for the size.** Independent marginal code, a chain that reproduces
  the shift (+0.128 against +0.114), and replays under both generative stories (+0.125 / +0.127
  against -0.001 / -0.009). This is the classical errors-in-variables attenuation and should be
  cited as such.
- **F24: B plus.** Ten runs, a second noise and template seed with the gaps reproduced to 0.0005,
  and the single-index baseline reported (0.0669 / 0.0691 against 0.0333 / 0.0180 for the GMRF).
  This is the one place in the program where an index field is shown to beat a scalar. Objections:
  the truth is "region-constant and spatially smooth at nside 16", that is, drawn inside the GMRF
  model class at the model's own resolution; templates are exact; T is fixed; hyperparameter ESS is
  3 to 14 of 700 and "ESS 2" in the first attempt
  (`nextgen/experiments/E021_index_prior_comparison_consistent_sky.md` line 502), where "the
  pre-registered remedy (extend the chains) was not applied" (lines 484 and 616); the intervals
  over-cover (78 / 98 percent).
  Settle: the same four-way comparison on a truth outside the model class (the PySM index maps at
  0.03 times, first with exact and then with data-story templates), with chains long enough for the
  hyperparameters.

### 3.8 The scale 1 sequence: F26 to F38

- **F26: A for the algebra, artefact for the size.** The leak fraction sigma_c^2 / (tau^2 +
  sigma_c^2) is established by sign-flip, roll and shuffle controls. It reaches 0.25 only because
  the per-pixel CMB noise is 58 uK at scale 1; by the finding's own scaling the index shift is about
  1e-4 at 0.03 times.
- **F27: B for the direction, C for the numbers.** Data-built templates bias the field at both noise
  levels. The per-quartile numbers are "one sky, one realization per cell, the data term without the
  index prior". The clause "a free multiplier per nside-16 region ... fixes the sync" is superseded
  by F29, F37 and F38 and should be struck.
- **F28: D.** Already "qualified" by the authors. It rests on profiled multipliers, which F31 shows
  to be biased by -0.2 to -0.4, and on scans that peaked at the edge of a +-0.16 grid. It should be
  withdrawn like F30, not qualified.
- **F29: B for existence and sign at nu_ref 30 GHz, C for "three quarters of the sky".** +0.15 +-
  0.05 over six simulated noise realizations (`nextgen/experiments/E026/stats_review_sync_level.md`
  line 27) is the only between-realization number in the sequence. The sky fraction is a property of
  the noise level, and after chain 3 the sign is a property of the reference frequency (Sect. 7.2).
- **F30:** withdrawn by the authors; credit.
- **F31: A.** Brute force to 1e-8 and quadrature to 5e-12. General (Sect. 5).
- **F32: D as a rule.** Three RMSE values per component from one realization with no Monte Carlo or
  realization error (0.143, 0.135, 0.147); all three are worse than a scalar index (objection 7);
  the pre-registered expectations were withdrawn or refuted in the record
  (`nextgen/experiments/E027_index_resolution_realistic_noise.md` lines 158, 183, 201); the heading
  says "Planck-like depth".
- **F33: C.** The width 0.1 is tuned to an exact template; two pre-registered tolerances were missed
  (`nextgen/experiments/E029_multiplier_hyperprior_width.md` lines 70 and 81); the GLS CMB residual
  rose by 0.91 +- 0.07 uK.
- **F34: C.** A negative statement from one pair of chains, at odds in emphasis with E029 line 131
  ("Neither hyperprior gives an unbiased dust index with a fixed T") and with F36 (Sect. 6).
- **F35: C.** T = 19.19 +- 0.54 K against a pseudo-true 19.90 +- 0.47 K differ by 1.0 combined sd
  (my arithmetic, treating them as independent); "Pre-registered criteria: both fail"
  (`nextgen/experiments/E028_sampled_dust_temperature_field.md` line 220); the reference value was
  re-read after the fact (line 308).
- **F36: C.** "No likelihood bias in either component at realistic noise" is a null result on one
  realization with two oracle inputs (the true T map, and a width tuned to an exact template). Dust
  offsets of +0.033 and +0.031 remain in Q2 and Q3 (`nextgen/papers/paperA_tables.md`, Table 3), and
  the record words the outcome as "a marginal miss with the mechanism identified, not as a pass"
  (E028 line 444).
- **F37: B minus and C.** B minus for "a width narrower than the template error biases the bright
  regions" (38 nats, 9 curvature sd, with the error-seed scatter +0.07, -0.06, +0.01 at the matched
  width). C for "the matched width is the rule", which chains 1 and 2 undercut (E030 line 599).
- **F38: C.** The mechanism's discriminating prediction failed in both pre-registered forms in chain
  3 (Sect. 7). The finding also labels stage A' numbers as stage A (Sect. 4.1).

---



## 4. Internal inconsistencies

### 4.1 Numbers and labels that disagree between documents

1. **The noise label.** "Planck-like" names scale 0.03 in
   `nextgen/experiments/E020_index_prior_comparison_pysm_d1s1.md` (lines 11 and 27 to 28) and scale
   1 in the outline (line 135) and in F32. By E031 neither holds for HFI (scale 0.03 is 2.6 to 33
   times the real white noise at 100 to 353 GHz, scale 1 is 88 to 1090 times). See objection 1.
2. **F38 compares against the wrong row.** F38 (`nextgen/FINDINGS.md` line 775 on) writes "Pixel
   RMSE 0.426 sync, 0.117 dust against stage A's 0.066, 0.038 with exact templates" and "96.4 uK
   against stage A's 93.1". In Table 3 of `nextgen/papers/paperA_tables.md` stage A (sd 0.5, oracle
   T) has 0.133 / 0.085 and 92.45 uK; 0.066 / 0.038 and 93.08 belong to stage A' (sd 0.1). The slip
   matters because it hides a confound: there is no exact-template chain at sd 0.3 at either
   reference (the chain directories under `/Users/trygvels/data/e023_pysm_n1/` are exact at sd 0.5
   and 0.1 and cs at sd 0.5 and 0.3), so every "template effect" read from chains 1 and 3 is a
   difference in template and in width at once. The like-for-like pair that does exist is chain 2
   against stage A (both sd 0.5, oracle T map): field mean -2.722 against -2.880, RMSE 0.327 against
   0.133.
3. **The CMB "floor" and the CMB column.** 93.18 uK is the "floor" of Table 3 and "the true-SED GLS
   floor" of F34, which already records a chain below it (92.3) with the remark that "a wrong mixing
   matrix trades bias for variance". Chain 3 is at 90.90 and its analysis prints "best cmb run minus
   the true-SED GLS floor: -2.185 +- 0.171 uK"
   (`nextgen/experiments/E030/chain3_e020_analysis.txt`). The column also mixes conventions: for
   chain 1 the outputs hold 96.42 (posterior-mean fields, `chain3_cmb_attribution_Tmap.txt`), 97.66
   (the `cmb rms all` column of `chain3_e020_analysis.txt`) and 98.60 +- 1.72 (per draw). Table 3
   uses 96.42 for chain 1 and an "at draws" 94.45 for chain 2; the outline (lines 194 to 196)
   compares 94.5 with 97.7. Whether the width "lowers the CMB residual" by 2.0 or by 3.2 uK depends
   on which pair is taken, and Sect. 8.7 shows that neither number measures foreground residuals.
4. **F29 contradicts its own sign.** "the volume term ... rises with a steeper index (slope +130
   nats per unit beta ...)" (`nextgen/FINDINGS.md` line 528 on). A slope of +130 per unit beta means
   the term rises toward flatter indices, which is what the stats review says
   (`nextgen/experiments/E026/stats_review_sync_level.md` item 2: "a steeper index gives each mu_r
   more leverage", hence a larger -Q_r and a lower volume term) and what the +0.13 offset requires.
5. **The k^(1/4) prediction** mixes units (objection 8):
   `nextgen/experiments/E014_q1_real_maps_lowres.md` lines 191 to 192.
6. **The resolution rule.** D019 (`nextgen/DECISIONS.md` line 233 on) withdrew "the finest the
   budget allows" because no run had tested it. E027 then measured the synchrotron optimum at nside
   8 (RMSE 0.143, 0.135, 0.147 at nside 4, 8, 16) and recommended "the finest field the budget
   allows under the GMRF prior" all the same
   (`nextgen/experiments/E027_index_resolution_realistic_noise.md` lines 211 to 214). F32 says "the
   optimum at nside 8 for synchrotron"; `nextgen/OPEN_QUESTIONS.md` (ranking of 2026-11-08) lists
   the rule as answered in the "finest" form. For synchrotron the finest rung is the worst of the
   three.
7. **Three width rules in one file.** F33's heading: sd 0.1 "removes the faint-sky synchrotron index
   bias". F37's heading: the width "must be at least the template error's scale", and in its text
   "must be matched and not wider". The E030 conclusion
   (`nextgen/experiments/E030_hyperprior_width_vs_template_law.md` line 599) and the outline (lines
   193 to 194): "a bias-variance knob with no calibrated setting".
   `nextgen/experiments/E029_multiplier_hyperprior_width.md` (lines 136 to 138) still carries
   "Recommendation (needs Trygve, a default change): multiplier hyperprior sd 0.1".
8. **The dust temperature.** `nextgen/experiments/E025_hier_cs_template_and_trueT_scans.md` lines
   102 to 104: "**Fixed T_d is not the problem.** ... Sampled T_d is deprioritized (D023)." Outline
   line 213: "with a fixed T the hier target has no unbiased configuration". F34's heading: "The
   fixed dust temperature is not the cause of the dust multiplier offset". F36: "the T model is what
   moves the dust RMSE from 0.056 to 0.038". The statements can be reconciled (they hold at
   different widths), but D023 rested on profiled scans that F31 later ruled out as evidence, and
   the log has no addendum saying so.
9. **D025's rationale is the refuted one.** The "Why" paragraph of D025 (`nextgen/DECISIONS.md`
   lines 359 to 385) says "the roughness check showed the level is a property of the posterior under
   a prior-dominated field scale". The stats review of the same day found "Roughness coupling: NOT
   supported" (`nextgen/experiments/E026/stats_review_sync_level.md` line 15) and F29 says the first
   version was "refuted by the review". The addendum to D025 (line 381) concerns the negative
   control only.
10. **"Converged" at ESS 4 to 13.** F38 gives "E026's hyperparameter ESS of 4 to 13" and its
    amendment of 12:20 calls E026's cs numbers "a converged chain's", on the ground that chain 2 has
    a similar RMSE.
11. **Ladder labels.** F35: "the sampled-T code path is validated at realistic noise (the L3 rung
    for the stage C T field)". `nextgen/VALIDATION.md` line 29 defines L3 as a run "on simulations
    drawn from the model being sampled, compared to the exact or independently computed target". The
    PySM sky is the L4 construct, and the reference was refined after it was pre-registered (Sect.
    4.4, row 13).
12. **Registry against findings.** `nextgen/experiments/registry.md`: the E012 row (line 21) still
    says "running" while F12 is marked final; the E022 row (line 38) quotes chains 1 and 2 while F22
    speaks of four; the E031 row (line 47) cites a finding F39 that `nextgen/FINDINGS.md` (last
    modified 12:21) does not contain.
13. **Stale headers.** `nextgen/FINDINGS.md` opens with "Established on 2026-10-25 in a single
    exploratory session" above 38 findings dated to 9 November. F19 says the L2 gate "is still
    running" while F16 reports it passed. `nextgen/experiments/E016_hierarchical_amplitude_prior.md`
    line 3 reads "(running: target written, no code yet)" on a record with a 600-replicate SBC.
14. **Two chains or one.** `nextgen/experiments/E021_index_prior_comparison_consistent_sky.md` line
    102: "Commander4 always runs two chains; the analyzer reads chain 1." Every chain directory I
    listed holds `chain01` files only. Either the second chain is written elsewhere and never read,
    or the statement is wrong; in both cases no between-chain diagnostic was ever computed.
15. **Paper scope.** `nextgen/ROADMAP.md` (lines 33 to 35) defines Paper A as "exact
    spectral-parameter sampling ...; the prior-volume problem and the ridge/Jeffreys resolution;
    hierarchical GMRF index fields with inferred smoothness; collapsed sampling efficiency". More
    than half of the outline (Sect. 5, lines 89 to 237) is about template multipliers, hyperprior
    widths and the dust temperature at one noise level, none of which is in that definition. F4
    calls ridge "the principled default"; the outline's thesis calls the hierarchical prior "the
    calibrated choice".

### 4.2 Claims without an experiment ID or a record

- F9 (`nextgen/FINDINGS.md` lines 98 to 99) has no experiment ID.
- E001 to E008 exist as registry rows only (`nextgen/experiments/registry.md` lines 10 to 17). There
  are no record files with question, setup, command, result and conclusion, which is what principle
  5 of `CLAUDE.md` asks for. F1, F2, F5 and F6 rest on them.
- The outline has no literature citation, although `nextgen/LITERATURE.md` exists and credits the
  Jeffreys remedy to earlier work (objection 9).
- F29's "+0.15 +- 0.05 over simulated noise realizations", "+130 nats per unit beta" and the sd 0.2
  and 0.1 values come from a review note (`nextgen/experiments/E026/stats_review_sync_level.md` item
  2) that points to "Scripts in `review_scripts/`" without commands or output files. The note says 6
  realizations; F29's scope line says "six to eight".
- "within a factor 2 of WMAP/Planck at this resolution"
  (`nextgen/experiments/E020_index_prior_comparison_pysm_d1s1.md` lines 27 to 28) has no source. By
  E031 it is about right for WMAP and LFI (0.7 to 2.7) and wrong for HFI.
- The statement that a 2024 Commander3 module carries the F1 structure (F10, F11) rests on code
  reading; the module was never run (objection 8).
- F38's amendment ("E026's cs numbers are a converged chain's") cites no convergence diagnostic.
- `nextgen/experiments/E032/` held scan outputs from 21:59 with no record file and no registry row
  at 22:16, so whatever was predicted for it is not in the repository. E033, opened at 22:09, did
  this in the right order (record and predictions first, output at 22:10).

### 4.3 Conclusions stronger than their evidence

| Statement | Where | What the record supports |
|---|---|---|
| "the calibrated choice for spatially varying indices" | outline lines 257 to 258 | Calibrated on data drawn from its own hyperprior (E016). On the one misspecified sky, four chains with a 30 percent template put the field mean at -2.72, -2.70, -2.72 and -3.16 (truth -3.00). |
| "has no likelihood bias in either component at realistic noise" | F36; outline line 217 | One sky, one noise realization, exact templates, oracle T map, sd 0.1 chosen knowing the template is exact; two of four pre-registered criteria missed and reassigned to the prior. |
| sd 0.1 "removes the faint-sky synchrotron index bias" | F33 | One chain: -2.981 +- 0.028 against -2.868. With a fixed T the same width moves the dust index instead (E029 lines 131 to 133). |
| "the sampled-T code path is validated at realistic noise" | F35 | 19.19 +- 0.54 K against a reference changed from 20.2 +- 0.5 to 19.90 +- 0.47 K after pre-registration; 22 divergences; implementation reviewed after the run. |
| "three valid chains by the D025 check"; "so the chain samples its target" | F32; F38 | The check is necessary, not sufficient (Sect. 8.4); the review that introduced it says so (`E026/stats_review_sync_level.md` item 1). |
| "The remedy is the template reference at the leading band" | outline lines 192 to 193 and 200 to 209; E030 lines 505 to 518 | A scan whose faint-quartile optima were -0.18 to -0.30. The chain then failed (Sect. 7). |
| "makes the hier prior usable with a component-separation-quality template" | E025 line 89 | Profiled scans (ruled out by F31), one error realization; contradicted by F38. |
| "depending only on the input zero levels" | outline line 248 | Two different band sets (objection 10). |
| "+0.019 (54 ridge sd)" | F3 | The sd comes from a fit with reduced chi-square 15000 (objection 9). |
| "the resolution contributing less than 10 percent of the index error" | E027 lines 211 to 214; F32 | Three chains on one realization; the bootstrap errors (0.001) are within-realization; all three rungs are worse than a single index for synchrotron (Sect. 8.3). |
| "what holds at both is the data-dominated regions and the Fisher-weighted level" | F38 amendment; outline lines 178 to 182 | Two chains that differ in width only. Chain 3 moved the hundred brightest regions from +0.09 to -0.126 and the Fisher-weighted level from +0.08 to -0.14 against truth. |
| "samples its stated target at L1 and L3" | F19 | One replicate at nside_beta 2 with one index field; the L1 field test at template brightness disagreed with its reference (objection 4). |

---

**Note added by the main session, 2026-11-09 23:01.** This report is incomplete. The agent was stopped by a
session interruption after Sect. 4.3 (file last written 22:22). Sect. 4.4 and Sect. 5 to 8, which
the text above cites (the chain 3 reading as "Sect. 7", the missing baselines as "Sect. 8.3", the
sufficiency of the D025 check as "Sect. 8.4"), were never written. What exists was read in full
and is acted on in `2026-11-09_program_review.md`, `../FINDINGS.md` (reading guide and grades) and
the E030, E033 and E035 records. Two of its points have since been overtaken: the k^(1/4)
arithmetic of F1 is replaced by a derived law (E035), and E032 has a record with its
pre-registration fixed before any chain output.

