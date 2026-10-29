# Decision log

Append new decisions at the bottom. Never rewrite an old entry; supersede it with a new one that
references it.

## D001 (2026-10-25): build on Commander4, in Python and C++
Commander4 already provides the engineering gains a rewrite would give (Python front end, C++
kernels via nanobind, ducc0, MPI, YAML config, ~860 tests) and is actively developed. Kernel speed is
at parity between Julia, Numba and C++ (F6). A separate Julia rewrite would duplicate an active
official code. Julia stays for fast prototypes only.

## D002 (2026-10-25): spectral-index MCMC must use an exact acceptance
The original resolve-then-chi^2 acceptance is not a valid MH step (F1). New `likelihood` option on
MCMC groups: marginal, ridge, conditional, resolve_chisq (legacy, warns).

## D003 (2026-10-25): `auto` resolves to ridge where available, else conditional
The flat-amplitude-prior marginal can be dominated by its volume term for full-sky indices (F3).
Ridge is the standard spectral likelihood and corresponds to a joint Jeffreys prior (F4). Revisit
after Q2 (coverage study).

## D004 (2026-10-25): spatially varying indices aim for the finest resolution the S/N supports
Trygve's direction. Hard region boundaries are unphysical and alms ring, but the goal is the
statistically optimal model. Leading candidate: hierarchical GMRF index fields with inferred
hyperparameters (Q3). Alms, clustering and fixed ladders are comparison baselines.

## D005 (2026-10-25): research scaffolding lives in `nextgen/`, work happens on branch `nextgen`
`main` keeps tracking upstream. Upstream pull requests are cut from focused branches, and only
after Trygve decides (Q13).

## D006 (2026-10-25): allow macOS development builds on `nextgen`
Upstream `CMakeLists.txt` stops on any non-Linux system. The long-running session will work on a
Mac mini, so `nextgen` allows Darwin with a warning, picks `-march=native` or `-mcpu=native` by
compiler check, and skips `-s` on Apple. Linux behaviour is unchanged (re-verified). Keep this out of
upstream pull requests unless Cosmoglobe wants it.

## D007 (2026-10-25): long-running research happens on the Mac mini, heavy runs on a cluster
The cloud workspace (2 vCPU, 7.8 GB, jobs stopped after each reply) cannot hold research runs. The
Mac handles T0 to T2 (COMPUTE.md); T3/T4 go to ITA or Sigma2 via Slurm.

## D008 (2026-10-25): sampler correctness is tested with a proper Gaussian amplitude prior; prior choices are a coverage study
Neither the flat amplitude prior (marginal, w = 1) nor the profile likelihood (ridge, w = 0) has a
proper generative distribution, so an SBC run with them cannot distinguish a sampler error from a
prior-volume effect (E012 review). The index-field prototype and every later sampler therefore pass
L1/L2 with the exact marginal under a proper Gaussian amplitude prior N(m, T) (a two-line change to
the collapsed kernel), and ridge / marginal / Jeffreys are compared in a separate frequentist
coverage study reported as bias and coverage versus S/N (principle 4). D003 (ridge default) stands.

## D009 (2026-10-25): index-field prior conventions
The hierarchical field drops the constant mode (the field mean is mu, sigma is the sd of
fluctuations about it), uses radians for the correlation length in both the GMRF and the alm family
with the graph-Laplacian eigenvalues rescaled by 1/(3 h^2), and samples the bounded log-uniform
correlation length through a logit transform with its Jacobian. Reasons in E012 v2.

## D010 (2026-10-26): Q1 is answered at nside 64; Q2 (prior choice) moves to the top of the ranking
E014 settled Q1 on real maps: the legacy move is 16 to 17x too wide and shifts beta_d by 3 sd, and the
exact `ridge` and `marginal` chains both reproduce independent grid targets. The same run showed that
the flat-prior volume term moves beta_s by 54 posterior sd with no change in width or chi-square, so
the prior is the dominant systematic for full-sky indices on real data, far above sampler width. The
ranking in OPEN_QUESTIONS therefore becomes Q2 (coverage study, E012 coverage run in progress, then
a scalar-index coverage table at M2.2), Q3/Q4 (index-field model and sampler, M2.x), then the second
half of Q14. Q1's remaining items (IQU, nside 128 scaling, HFI full set) are follow-ups in E014, not
blockers. The `conditional` Gibbs mode is dropped from real-data comparisons (F14: conditional sd
2e-6 against marginal 3.6e-4; of order 3e4 iterations per independent sample); it stays in the code
as a documented option. No default changes.

## D011 (2026-10-27): zero-level convention for real-map inputs; `hier` is an option, not a default; NUTS start state
Real-map inputs for E014-style runs use one zero-level convention for every band: offset_b = mean
over |b| > 30 of (map_b minus the DR1 sky model evaluated at band b), `dr1_model_zero_levels.py`
(pass 4), with the Haslam non-solar dipole left in by default and removed as a documented variant.
This is a working convention for comparing methods on identical inputs, not a claim that it is
right (F17, Q15); a per-band zero-level parameter in the model is the proper fix and waits for
Trygve's decision. The M2.4 stage-1 index-field sampler enters Commander4 as `likelihood: hier`, a
non-default MCMC group option with its own `nuts` block and component keys (`spectral_index_nside`,
`spectral_index_field`, `amp_prior`); `auto` still resolves to ridge (D003) and no existing run
changes. A hier group without a restart state starts the field at its mean (z = 0, mu = the
component's beta, sigma = sigma_scale / 2, rho mid-range) and the amplitude hyperparameters at their
prior means; a restart from a field chain at the same nside resumes from the stored field (NUTS
adaptation restarts). Reason: a deterministic start that warm-up leaves quickly, and no dependence
on the parameter-file proposal widths the MH groups need.

## D012 (2026-10-27): NUTS stays the production sampler for the index field; M2.3 closed for Phase 2
E019 (F20) measured MAMS against our NUTS on the non-centred hierarchical index-field target at
K 11 to 191 with a diagonal preconditioner and found parity in effective samples per gradient
(0.9 to 1.2 at the median pixel, 0.3 to 1.1 at the worst pixel) with a 2 to 20 times dearer
warm-up; the pre-registered threshold for adopting a second sampler was a factor 2. NUTS
(`math_utils/nuts.py`, dual averaging, windowed diagonal mass adaptation) therefore remains the
only index-field sampler in Commander4; the MAMS prototype stays in `nextgen/prototypes/indexfield/
mams.py` with its self-test for the Q5 joint-target experiment, where d is 1e5 and above and the
literature's gains are reported. Non-diagonal preconditioning (Fisher, low-rank) is deferred to
the same experiment. Also fixed by E017/E018: the hier amplitude prior is the recommended option
below full-sky S/N of order 100 (F18), but changing the default from ridge needs Trygve (D003).

## D013 (2026-10-27): per-band zero levels enter the collapsed per-pixel likelihood as offset columns with Gaussian priors
Commander4 gets per-band zero-level (monopole) parameters for maps-only runs on real intensity
data, implemented as one linear offset column per band in the per-pixel mixing matrix of the
collapsed likelihood, with a Gaussian prior per band (M3.3). Not as a separate monopole component
sampled in the CG step (the Commander3 `md` design). Reasons: (1) the offset is linear in the data,
so it is profiled or marginalized by the same kernel as the amplitudes at no sampling cost, and the
index posterior then carries the zero-level uncertainty instead of inheriting an input convention
(F17: the convention moves a full-sky index by hundreds of posterior sd); (2) with spatially
constant SEDs, n_component directions of the n_band offsets are exactly degenerate with the sky
mean of each component amplitude and can only come from priors, while the remaining n_band minus
n_component directions are data-constrained; a spatially varying index weakens even the first
set, which is why zero levels and indices are coupled in the first place, so they must be sampled
jointly; (3) the CG path is not in the collapse yet (Q5, M3.1), and putting the offsets in the
collapsed step keeps Phase 2 and Phase 3 results comparable. Prior centres come from a stated
convention per band (HI cross-correlation at high latitude for the dust-dominated HFI bands, the
DR1 sky model or a masked monopole fit for WMAP and LFI); prior widths from the measured
disagreement between conventions (E014 passes 2 to 4), never zero. Polarization has no monopole,
so Paper B's synchrotron field does not depend on this; every intensity result does. Trygve asked
for the choice to be made here (2026-10-27) rather than deferred to him. Q15 item (2) is closed by
this entry; item (3), propagating the prior widths, follows the implementation.

## D014 (2026-10-28): E020 production stopped; the index-prior comparison moves to a model-consistent sky (E021); the misspecified PySM sky becomes a robustness experiment at realistic noise

**Decision.** Stop the E020 ten-run production (attempt 3 stopped at chain-1 iteration 185 after
6.7 h; 85 h per run projected). Run the same ten options, metrics and analyzer as E021 on a sky
drawn from the model class: region-constant indices at nside 16 equal to the PySM region means,
one modified blackbody at the fixed T, PySM amplitudes at the pivots, the c1 CMB, templates equal
to the true amplitudes with an iid 30 percent scatter, the same bands and noise. Keep the PySM sky
with its misspecifications for a later robustness experiment at realistic noise (Q16), after the
prior comparison is done. The E020 data and chain outputs are kept as the evidence for F21.

**Why.** E020 on the PySM sky was measuring how each prior handles a runaway driven by fixed T_d
and a 60-percent-wrong template (F21), not how well it recovers an index field; every option would
have scored on the same pathology. Principles 3 and 4 say: simulations from the model first, then
misspecification as its own question. The sampler was exonerated by direct computation, so the
code does not change; the design does.

**Sampler settings for E021.** `init_inv_mass: fd_hessian` stays; `adapt_mass` returns to its
default true (the windows were cleared of blame by attempt 3: the saturation was the funnel
climb). The E021 timing run decides between the two settings on gradients per transition after
the first window. Iteration counts are set from that timing run so that the ten runs fit in
about two days on the Mac mini (D-level: a shorter chain with the same burn fraction, logged in
E021, not a change of target).

## D015 (2026-10-28): E021 runs with a warm-up metric refresh, 1200 iterations and two concurrent lanes

**Decision.** The E021 production uses `nuts.refresh_inv_mass_at: [50, 200]` (the diagonal
inverse mass recomputed from the finite-difference Hessian at the current state during warm-up,
step size re-found; c707a34), 1200 Gibbs iterations with warm-up 400 and burn 500 instead of
E020's 2000/500/600, and two runs at a time on the Mac mini.

**Why.** The E021 timing run cost 255 to 1023 gradients per transition on the start-state
Hessian metric; an offline test from the iteration-40 state showed a fresh Hessian brings it to
136 (a fixed hyperparameter block does not help: 187), so the cost was the stale metric, not the
funnel. With 136 gradients at 83 ms each, 2000 iterations are 13 h per nside-16 run and ten runs
over a week; 1200 iterations keep 700 post-burn samples, and the analyzer's MC inflation of the
RMSE says whether that is enough per comparison. The refresh is adaptation only (any positive
metric is a valid kernel; the step size restarts as at a Stan window end), so the target is the
one reviewed for E018 and E020; the E021 Result records it and the stats review of the E021
sampler configuration is run while the production is going, not after.

**Amendment (07:17):** `adapt_mass: false` for the E021 production, so the last refreshed Hessian (transition 200) is the metric for the rest of the run; with the windows on, the one ending at 250 would have replaced it by a 50-sample variance estimate (stats review of c707a34, note 3). Review items on the refresh (terminal buffer, once-per-count firing) applied before the restart.

**Amendment (10:31):** one E021 job at a time. Measured at chain-1 iteration 70: two concurrent runs cost 150 to 220 ms per gradient against 82 to 95 ms for one (the pre-registered 1.5x limit in the E021 record), so two lanes deliver the throughput of one. Lane B is paused with SIGSTOP and the jobs run in start order (scratchpad `e021_serializer.sh`); nothing is restarted, the paused run keeps its state and the kernel is unaffected (a stopped process resumes where it was).

## D016 (2026-10-28): the BP XV polarized index cell run as E022 on a surrogate plus the released chains; Paper A presentation deferred to Trygve

**Decision.** The planned E014 follow-up (F14 scope: does the Gibbs conditional mixing issue
bear on BeyondPlanck XV's reported uncertainties) was run as its own experiment, E022, on an
nside-16 surrogate of the BP XV polarization data (LAMBDA's nside-16 WMAP maps and dense
covariances fetched for it), and then directly on the four released BP chains by reading the
per-region index values over HTTP range requests from the 1.8 TB chain files (613 MB fetched
per chain). The stats-reviewer was run five times on it (targets, variants, chain analysis,
resolution scaling, element count). Paper A keeps F14 restricted to intensity at Planck
resolution; the polarization statement is the dichotomy of E022's conclusion. How the BP XV
chain diagnostics are presented in Paper A is left to Trygve (STATUS, Needs Trygve).

**Why.** The surrogate alone could only bound the coupling from below (the flat-prior volume
term and the stationary target scale with the number of amplitude elements, reviewer rounds
2 and 5), and the released chains alone cannot separate an unsettled flat-prior chain from a
stationary one driven by the other Gibbs blocks (BP's effective element count is
undetermined between 7k and 4M per region). Writing a verdict on a published result of
Trygve's own with that gap open would overstate; the record holds the numbers and the
discriminator (Commander3 prior semantics or one conditional-width measurement on a
released sample). The resolution scaling on the real data (flat-prior Spur target from -2.88
to the prior boundary between nside 16 and 64; stable with data-scaled priors on all
components) goes into Paper A Sect. 4 as the real-data case for the prior-volume effect,
independent of the BP XV question.

## D017 (2026-10-28): E021 attempt 1 stopped after its first complete run; the sky regenerated with the hier prior's own template story and the ten runs restarted one at a time

**Decision.** Stop the E021 queue (gmrf16 complete, iid16 at iteration 200, indep16 just
started) at 17:11, change the writer so that the templates are the PySM amplitude maps and the
true amplitudes are a = g r (`--template-story model`, the new default for the consistent
sky), regenerate the sky in `/Users/trygvels/data/e021_consistent_v2` with the production
sampler settings written directly (1200 iterations, warm-up 400, refresh at 50 and 200,
`adapt_mass: false`), verify the target on the new data (exact marginal at the truth, likelihood
probe, gradient check) and restart the ten runs in one serialized queue (gmrf16 first). Keep
the attempt-1 data and chain as evidence (F23) and as a Q16 data point. Stats review of the
diagnosis and of the new sky runs alongside the first run; if it refutes, the queue stops.

**Why.** The complete gmrf16 chain on the attempt-1 sky put the dust index of the faint three
quarters of the sky 0.13 above the truth at four stated sigma and the amplitude
hyperparameters 2 to 24 sd off their targets, and the exact marginal shows this is the target,
not the sampler: the attempt-1 templates g = a / r give a | g a law different from the prior's.
E021 exists to measure prior behaviour in prior-dominated regions on a sky from the model
class; on that sky it would have measured the template-law mismatch in every option. Ten hours
of compute lost now against forty had the queue run on.

**Alternatives considered.** Keep the attempt-1 sky and report the bias as the result: it is a
useful robustness number (F23, Q16) but not the clean baseline Q3 and M2.5 need. Fit the ten
options on both skies: deferred; Q16 does the misspecified cases at realistic noise.

## D018 (2026-10-29): E021 closed; the hierarchical GMRF is the index prior of record; the index resolution rule is the finest the budget allows; a second realization of the nside-16 triple runs before the paper table is final

**Decision.** Close E021 attempt 2 with its ten runs (F24). Adopt the hierarchical GMRF (graph
Laplacian basis, learned mean, amplitude and correlation length) as the index prior of record
for M2.x and Paper A; drop the exchangeable (iid) prior as a candidate; keep independent
regions and the single index only as references. Set the resolution rule for the index field
as the finest nside the compute budget allows under the spatial prior, not a coarse ladder.
Before the E021 table goes into Paper A, run gmrf16, iid16 and indep16 once more on a second
noise and template realization (seed 2028, about 14.5 h, started 11:13) with the pre-registered
expectation written in the E021 Follow-ups. Treat the slow mixing of the dust field mean in the
centred parameterization as an M2.4 efficiency item.

**Why.** The ordering GMRF, iid, independent holds on both indices at nside 16, 8 and 4, with the
nside-16 gaps at 13 and 9 bootstrap SE against iid and the mechanism of the iid loss verified
on the exact target (plane-weighted global mean). Under the GMRF the region-level error does not
change from nside 4 to 16 while everything the coarse field cannot represent is lost, so there
is no accuracy argument for coarse regions, only cost. One realization is enough to establish
the ordering and the mechanism but not the size of the gaps, hence the second seed for the
pair that matters.

**Alternatives considered.** Run the full ten-run ladder on the second seed (40 h): the nside-8
and 4 runs add nothing to the decision. Run nside 32 on this sky (about 40 h) to test a field
finer than the truth: deferred to Q16's PySM sky, where the truth has structure at every scale
and the question is the real one. Keep iid as a cheap option for low-resolution runs: rejected,
it is biased wherever the sky differs from its global mean, and at nside 4, where it is closest
to the GMRF, it is still the most biased of the three.

## D019 (2026-10-29): D018 amended after the third stats review; the resolution rule is withdrawn until a field finer than the sky's structure is tested

**Decision.** Keep the GMRF as the index prior of record and the seed-2028 replication (D018).
Withdraw the D018 clause "the index resolution rule is the finest the budget allows" as a
result: E021 showed that the nside-16 GMRF keeps most of the per-region averaging of coarse
regions (within 13 / 30 percent of the nside-4 model's per-region error against the
S/N-weighted reference, where independent regions lose factors 3.5 / 4.9) and keeps the
structure coarse fields cannot represent, but no run had a field finer than the truth's
structure, which is the real-data case. The rule returns as a hypothesis for E023 / Q16 (PySM
sky with structure at every scale). Negative controls and cross-run rules in future records
quote the statistic's reference explicitly (weighted or unweighted region mean), and the
MC-inflation rule is applied to every pair with the quadrature arithmetic.

**Why.** The reviewer showed the "flat per-region error" reading rested on the analyzer's
unweighted reference, which measures misfit, not per-region cost; against the weighted
reference the error falls with coarser regions as the data per region grow. The ordering
results stand.

## D020 (2026-10-29): the dense Hessian metric is the production setting for index-field NUTS runs at high S/N; the default stays diagonal

**Decision.** Production index-field runs (E023 stage 2, the seed-2028 follow-ups, Paper A runs)
use `nuts.init_inv_mass: fd_hessian_dense` with `adapt_mass: false` and refreshes at 50 and
200; the code default stays `ones` (principle 7; the simulated parameter writers set the
production value). Adopted on E024 although its third pre-registered clause (no region index
worse than a factor 1.5 per gradient) failed for the sync regions (1.6 to 2.2): the clause was
set against single quantities, and the chain's worst quantity, which bounds its usable length,
improved 14 times; the failure is recorded in the record and here.

**Why.** The dust field mean and the faint dust regions were the bottleneck of every E021 GMRF
run (ESS 5 to 14 of 700, hyperparameter inference unreliable at the run lengths the Mac mini
allows); the dense metric makes them mix at several hundred ESS with the posterior unchanged.

**Alternatives considered.** Non-centred field mean or pixel coordinates: a reparameterization
that changes the geometry for every quantity and would need its own validation; deferred,
since a Gaussian metric already fixes the observed coupling. Longer chains: 30 times the
compute. Per-block step sizes: the next efficiency item for the sync-side loss.

