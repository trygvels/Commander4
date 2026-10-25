# E015: BeyondPlanck absolute bandpass shift was accepted on rank 0's chi-square only

- **Date started / finished:** 2026-10-25 / 2026-10-25 (code reading; numerical impact not measured)
- **Branch and commit:** Commander3 BP10 (3aa493e) and HEAD (f760333)
- **Open question addressed:** Q14 (scope: Commander3 production samplers), Q13
- **Validation level:** L0 (code reading by one stats-reviewer agent, confirmed by me line by line)

## Question
Is the TOD-level chi-square used to accept or reject the absolute bandpass shift proposal reduced
over the MPI ranks that hold the band's scans?

## Target and priors
Stated: MH on p(delta_abs | d, sky, gains, n_corr) proportional to exp(-chi2_TOD(delta)/2) over all
scans of the band, flat prior (see E010 addendum).

## Setup
Read `comm_tod_lfi_smod.f90` (process_TOD), `comm_tod_driver_mod.f90` (compute_chisq_abs_bp),
`comm_tod_mapmaking_mod.f90` (finalize_binned_map), `comm_tod_bandpass_mod.f90` (sample_bp) at
both tags.

## Reference
The reduce for the relative-shift statistic inside `finalize_binned_map` (lines 565-570 at BP10)
shows what the code does when the chi-square is meant to be global.

## Result
- `compute_chisq_abs_bp` (driver, no MPI calls) adds each accepted scan-detector chi-square of the
  calling rank into `chisq_S(j, k)`.
- In `process_TOD` the absolute branch calls `finalize_binned_map(self, binmap, handle, rms_out,
  1.d6)` WITHOUT `chisq_S` (BP10 line 849; HEAD line 913), so the reduce at 565-570 does not run.
- `sample_bp` computes `cp = sum(chisq_S(:, k))` on `tod%myid == 0` and broadcasts the decision.
Hence the acceptance uses only the scans owned by rank 0. Same structure at HEAD and in the WMAP,
DIRBE (Cosmoglobe DR2), SPIDER, LiteBIRD and QUIET TOD modules.

Expected effect: the posterior of the 30 GHz absolute shift is the conditional given a fraction
1/N_ranks of the data: wider by roughly sqrt(N_ranks) and centred on rank 0's scans. The published
BP 30 GHz bandpass shift uncertainty is therefore not the full-data conditional. The size of the
effect depends on the number of ranks per band in the production run (not recorded here).

## Negative control
None (code reading). A one-rank versus many-rank run of Commander3 on the same simulated TOD would
show the width difference directly.

## Conclusion
Reproducible code defect, present at HEAD. Reporting it to Cosmoglobe is Trygve's decision (Q13).
Commander4 does not have a bandpass sampler yet (FINDINGS section 5); when one is written, its
test must include a one-rank versus many-rank equivalence check.

## Follow-ups
- Measure the effect on simulated LFI TOD with Commander3 if a cluster becomes available.
