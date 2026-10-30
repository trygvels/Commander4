# E026: Stage-2 chains on the realistic-noise PySM sky: hier prior with a region-scale template multiplier field against the scalar law, component-separation-quality templates

- **Date started / finished:** 2026-10-30 /
- **Branch and commit:** nextgen (mu-map prior under stats review before the runs)
- **Open question addressed:** Q16 (D022, D023), Q3 / M2.5 (first chains at realistic noise), M2.4 (dense metric at nside 16 with the reviewed code)
- **Validation level:** L3 (misspecified simulation with known truth; the sky is not from the model class: indices vary inside regions, template errors are log-normal and smooth)

## Question
Do the E025 exact-target results carry over to the posterior under the GMRF index prior: with
component-separation-quality templates (true amplitude times a 10-deg, 30 percent log-normal
error field) is the gmrf16 posterior of the index fields at realistic noise calibrated against
the noise-free region fit when the hier prior carries a template multiplier per nside-16
region (`mu_nside: 16`), and biased when it carries the scalar law; and how close does the
mu-map run get to the attainable reference, the same chain with exact templates? This is the
first index-field chain at realistic noise, with the flat CMB prior (D021) and the dense
Hessian metric (D020) at nside 16 (dim about 12,300 with the mu maps).

## Target and priors
The E021 gmrf16 target on the E023 noise-1 sky (`/Users/trygvels/data/e023_pysm_n1`: PySM
d1s1c1, nside 64, 12 bands 23 to 353 GHz, white noise at the E012 levels) with: CMB amplitude
prior N(0, 1e4^2) (flat in effect, D021); sync and dust hier priors with the templates below,
`mu_prior [1, 0.5]` per mu entry, `log_s_prior [log 0.3, 1]`; GMRF index fields at nside 16
(centred, graph-Laplacian basis, alpha 2, learned mu, sigma, rho; index prior means -3.1 and
1.6 with sd 0.3 on the field mean); dust at the fixed T 20.93 K (E025: the T misfit changes no
shift by more than 0.01 at this depth). Three runs, 1200 Gibbs iterations, warm-up 400, burn
500, dense metric refreshed at 50 and 200:
- gmrf16_cs_mumap: templates `template_{sync,dust}_cs_*.fits` (true amplitude times exp(f), f
  the E025 field, seed 11), `mu_nside: 16` on both hier priors (3072 multipliers each).
- gmrf16_cs_scalar: the same templates, scalar mu (the model E023/E025 showed biased).
- gmrf16_exact_mumap: templates equal to the true amplitudes, `mu_nside: 16`: the attainable
  reference at this depth (its mu map should sit at 1 and its s at the guard).

## Setup
Parameter files `param_gmrf16_{cs_mumap,cs_scalar,exact_mumap}.yml` in the data directory
(from `param_gmrf16.yml` with the template paths, CMB rms 1e4, `init_inv_mass:
fd_hessian_dense`, `mu_nside`), 12 MPI ranks on the Mac mini, one run at a time (about 10 h
each at the E021 gmrf16 rate; the dense Hessian at dim 12,300 costs 24,600 gradients per
refresh plus an eigendecomposition of a 1.2 GB matrix, to be timed). Analysis:
`e020_analyze.py` (index RMSE against the PySM truth and against beta_eff3, GLS CMB residual,
z and coverage against the S/N-weighted truth region mean; the mu map enters its amplitude
trace as its mean) plus a dedicated script for: index RMSE and z against the noise-free
region-fit reference of E023 (`E023/scan_n1_refs.npz`, fixT), per template-S/N quartile; the mu
map against the realized region means of exp(-f); s against the within-region residual; ESS
per gradient of every block (the dense metric at this dimension is itself under test).

## Reference
Truth: the PySM index maps and the E023 region-fit reference (fixT). Pre-registered reading:
(1) calibration: in gmrf16_cs_mumap the z of the region indices against the region-fit
reference has |median| below 0.3 and MAD within 0.5 to 1.5 in every S/N quartile where the
block sd is below 0.05 (dust Q1 to Q3, sync Q1), and no region has |beta_s - reference| above
0.5; (2) bias where the scans predicted it: gmrf16_cs_scalar has a dust z median above +0.5 in
Q2 and Q3 (E025: shifts beyond 0.1 at block sds 0.03 to 0.05) and the sync Q1 z median below
-0.3; (3) attainable: the cs_mumap index RMSE per quartile is within a factor 1.5 of
exact_mumap's; (4) hyperparameters: the cs_mumap mu maps track the realized region means of
exp(-f) (correlation above 0.8 over the regions of the two brightest quartiles) and s sits
near the residual 0.06 (E025); (5) CMB: the GLS CMB residual of cs_mumap within 10 percent of
exact_mumap's. Efficiency (secondary): ESS per gradient of the mu blocks and of the field
means; the dense metric's eigendecomposition time at dim 12,300. If (1) fails while (3) holds,
the posterior widths (the index prior) are the issue, not the amplitude prior; if (1) and (3)
fail, the E025 scans missed something that only the joint posterior shows and the record says
so.

## Command
```bash
for t in gmrf16_cs_mumap gmrf16_cs_scalar gmrf16_exact_mumap; do mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_$t.yml; done
python nextgen/scripts/e020_analyze.py --dir /Users/trygvels/data/e023_pysm_n1 --burn 500 --tags gmrf16_cs_mumap gmrf16_cs_scalar gmrf16_exact_mumap
```

## Result
(pending)

## Negative control
gmrf16_cs_scalar (the scalar law the scans showed biased) must show the predicted dust and sync
biases; if it does not, either the scans or the chains are wrong and the record says which.

## Conclusion
(pending)

## Follow-ups
