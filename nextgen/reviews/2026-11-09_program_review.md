# Program review, 2026-11-09: what stands, what was wrong, and the way forward

- **Written by:** the main session, at Trygve's request for a proper review of the work so far and
  a second look at the findings and the path forward.
- **Repository state:** branch `nextgen`. The review started at 52061ad and its corrections are in
  the commits that follow it.
- **Inputs:** four agent reports written without sight of each other
  (`../research/2026-11-09_novelty_check_methods.md`,
  `../research/2026-11-09_science_landscape_and_data.md`, `2026-11-09_referee_report.md`,
  `2026-11-09_pipeline_assumptions_audit.md`), three more commissioned on the conclusions of this
  document (`2026-11-09_way_forward_challenge.md`, `2026-11-09_stats_review_legacy_move.md`,
  `../research/2026-11-09_noisy_mh_and_translation_move_prior_art.md`), and eight experiments opened
  during the review (E031 to E038), which ran into the early hours of 2026-11-10.
- **Evidence levels:** numbers cite an experiment or a finding. Literature statements come from
  agent reports and are unverified against the PDFs unless marked.

Section 0 is the verdict. Sections 5 to 9 answer "are we on our way to publishable and useful
results": what the review found, where effort buys most, the plan, the decisions it needs and the
assumptions it rests on. Sections 2 to 4 are the audit of what exists.

## 0. Verdict

1. **There is a publishable methods paper in this work, and it is narrower and sounder than the
   outline.** Its core: what the legacy index move converges to, the incidental-parameter problem
   of per-pixel amplitudes with a formula for the flat-prior shift and for the profile's width,
   the limits of priors tied to a template, and a polarization section at real depth and on real
   data. Everything in it is A or B grade and computable as exact grids. Six weeks to a draft.
2. **The largest block of effort was in the wrong regime.** Thirteen findings on template
   multipliers at "realistic noise" were obtained at 24 to 90 times real intensity noise, and
   they refine a prior that the sound estimators do not need. They shrink to one paragraph on
   limits. The label had never been compared with a map.
3. **The most interesting results came out of the review itself** (Sect. 5): the legacy move
   converges to the profile posterior at temperature sqrt(2 pi k)/4 + 1/2; a flat amplitude prior
   shifts the index by ln(nu_ref / nu_pivot) times the single-amplitude index variance, which is
   0.06 to 0.9 at real polarization depth with a sign set by the reference frequency and the CMB
   prior, and which the real-data numbers of E022 obey to 0.2 percent; the profile likelihood with a
   sandwich width, or equivalently a zero-mean prior with one free scale per region, is
   calibrated there and needs no template.
4. **The pipeline's assumptions are weaker than its samplers.** Two ingest bugs are confirmed on
   real maps and every real-map width is withdrawn. A second chain was computed and discarded in
   every run. Per-pixel amplitudes are stored in a form that loses a quarter of their scatter and
   shifts any acceptance evaluated on them. The full-resolution amplitude path is untested.
   Neither of the audit's rankings of what limits real-data results is led by a sampler.
5. **Way forward:** Paper A1 as in item 1. A K band against 30 GHz null test with dense
   covariances before any polarized field code, because those two bands carry 91 percent of the
   index information and are known to disagree. Paper B is decided on its outcome. Intensity
   template work, the dust temperature field and sampler efficiency are frozen.
6. **The plan was attacked before it was written down.** An independent reviewer rejected the
   first version ("two half-papers") and proposed the reduced one, which is adopted here with one
   open choice for Trygve.
7. **Not yet known:** the outcome of E032 (running through 2026-11-10), whether a trustworthy LFI
   low-resolution covariance exists publicly, and whether the public data determine a
   high-latitude index at all once K band and 30 GHz are allowed to disagree.

## 1. How the review was done

1. **Independent reads.** A novelty check of the six method claims against the literature. A
   survey of the science landscape and the public data. A hostile referee report that graded
   every finding F1 to F38 against its experiment record (stopped by a session interruption after
   its Sect. 4.3; what exists was read in full). A code audit that listed 48 assumptions of the
   pipeline as coded, 28 of them mentioned nowhere in `nextgen/`.
2. **Experiments to settle what the reads raised.** E031: the noise regime of the simulations
   against real maps. E033: what the CMB metric used since E026 measures. E034: the audit's two
   ingest suspicions, on the real maps. E035: what the legacy index move converges to. E036:
   amplitude priors for signed Q/U amplitudes at real depth. E037: the validation harness over 40
   seeds. E038: first input checks for polarization. E032: the template-prior model at real
   intensity depth with four chains (running).
3. **Corrections to the record.** `FINDINGS.md` now opens with a reading guide and a grade per
   finding. E030 is closed with the chain 3 verdict. The registry, the literature table and the
   data notes carry the flagged corrections. `VALIDATION.md` has eight new rules.
4. **A second look at the conclusions.** Three more independent agents were given the proposed
   way forward, the legacy-move derivation and the prior-art question, with instructions to find
   what is wrong.

## 2. What stands

Grades as in `FINDINGS.md` section 0. "Needs" is what a referee would still ask for.

| Result | Evidence | Grade | Needs |
|---|---|---|---|
| The legacy "re-solve, then accept on chi-square" index move is a pseudo-marginal kernel for the ridge posterior that the Gibbs amplitude draw breaks. With one index step per amplitude draw the chain converges to the ridge posterior tempered by sqrt(2 pi k)/4 + 1/2: width 0.79 k^(1/4) | E035 (derivation, toy to 0.7 percent, exact controls), an independent derivation with an exact-kernel reference, E002, E005, E014 | A, real-map size B | A Commander4 run with map storage to confirm the explanation of the real-map width; citations to Ball et al. 2003 and Andrieu and Roberts 2009 |
| Exact index samplers for stated targets (marginal, ridge, conditional), end to end and on real maps | E005, E006, E014 (0.07 sd, 5 percent in width) | A | Nothing |
| With flat amplitude priors the index posterior carries a volume term. For a spatially varying index the flat-prior marginal is unusable at any S/N, the profile is calibrated only at high S/N | E006, E012, E017 (F3, F12, F15, F18) | A on simulations | Citations to 0709.1058, 0811.4280, 1608.00551. The new part is the coverage curve against S/N |
| A hierarchical amplitude prior passes simulation-based calibration inside its model class | E016, 600 replicates (F16) | B | SBC of the Commander4 implementation; a scope clause: data drawn from its own hyperprior |
| The Gibbs conditional cannot mix where the amplitude S/N is high; collapsing the amplitudes cures it | E001, E014 (ratio 180), E022 (F5, F14, F22 a) | B | Three seeds per efficiency number; the qualitative result is known |
| A GMRF index field with inferred smoothness beats exchangeable and independent-region priors on a model-consistent sky at real depth | E021 (F24) | B plus | The same comparison on a PySM sky at scale 0.03 with a scalar-index baseline and four chains |
| A data-built template makes the collapsed hier posterior run away at real depth; the hier prior is exact only for its own template-error law | E020, E023 (F21, F23) | B plus, A for the mechanism | Nothing; these define the scope |
| At noise scale 1 with a 30 percent template error no hyperprior width and neither reference frequency gives a calibrated field mean | E030 (F37, F38), four chains | B as a description | The cause is not established; E032 decides whether it matters at real depth |
| Zero levels move real-map indices by 0.035 (beta_s) and 0.016 (beta_d) like for like | E014 passes 2 to 4 (F17) | B | A model that fits; per-band offsets (Q15) |
| The Cosmoglobe DR1 dust amplitude map is referenced at 70 GHz, not the 545 GHz of its header | E014 (F13) | B | Confirmation from the authors (Trygve's call) |
| Commander3: the `marginal` likelihood has a sign error; the BP bandpass shift was accepted on rank 0's chi-square | Code reading (F10, F11) | B as code facts, C for impact | Never run |
| Two ingest steps make the real-map error model wrong | E034 (F41) | A, B | A corrected reader and a rerun (Q19) |
| The CMB residual metric of E026 to E030 was the noise term | E033 (F40) | B | Nothing; the metric is replaced |

## 3. What was wrong or overstated, now corrected in the record

1. **The noise regime.** "Realistic noise" for scale 1 of the PySM ladder was wrong by a factor
   of 24 to 90 for WMAP and LFI intensity (F39). Thirteen findings (F26 to F38) and 44 percent of
   the Paper A outline live in a regime whose only real counterpart is polarization, where none
   of them was obtained. The label came from a test harness and was never checked against a map.
2. **The CMB metric.** Models were ranked by a total residual that is the noise term (F40). Every
   such ranking is withdrawn. At high latitude no chain with a component-separation-quality
   template beats a constant index.
3. **Chain 3.** The remedy written into the outline (reference at the leading band) failed, and
   the scan reading it was built on was a misreading: faint-quartile optima of -0.18 to -0.30 were
   called flat (E030).
4. **Real maps.** The E014 chains saw an error model 3 to 28 times tighter per pixel than the
   record's table, with high-latitude reduced chi-squares of 13 to 60 (F41). "54 ridge sd" and
   "hundreds of posterior sd" are dropped. F17's headline of 0.4 compared two band sets.
5. **F1.** "Not a valid MCMC step" and an empirical k^(1/4) fit with a prediction the referee
   could not reproduce are replaced by a derivation (E035).
6. **Single chains.** Every result rests on one chain per configuration, with hyperparameter
   effective sample sizes down to 4 to 13, although each run computed a second chain and threw it
   away (F42).
7. **Rules drawn from single runs.** "The finest field the budget allows" (F32), "sd 0.1 removes
   the bias" (F33), "the sampled-T path is validated" (F35), "no likelihood bias at realistic
   noise" (F36), "close to unbiased for a component-separation-quality template" (F28), "the
   calibrated choice for spatially varying indices" (outline). Each is graded C or withdrawn.
8. **Pre-registration discipline.** Criteria were reassigned after the fact in F35 and F36, and
   one timestamp was written before the time it named (E032, corrected).

The pattern behind items 1, 2 and 7 is the same: a number was read without the trivial baseline
or the real-data anchor that would have shown what it meant. The rules added to `VALIDATION.md`
(regime gate, trivial baselines, four chains, replication, pooled scan readings) are aimed at
exactly this.

## 4. Pipeline and modelling assumptions, challenged

The audit's central finding is that the validated region and the real-data region barely overlap.
Everything behind F12 to F38 is intensity only, white noise, simulated at the analysis
resolution with no beam and no smoothing. The table lists the assumptions that matter most,
with what is now known. "Checked" means by this review, after the audit.

| Assumption in the code | Status | What it can bias | Cost to settle |
|---|---|---|---|
| The model-error floor is a model error (audit B14) | Checked, false: it averages down 27.6 times (E034) | Every real-map width; which bands carry the fit | Reader fix, small |
| Degrading before smoothing is harmless (B1) | Checked, false: 0.6 to 1.7 uK between bands (E034) | Real-map chi-square and indices at high latitude | Reader fix, small |
| Smoothed noise is independent between pixels (D2, F8) | Known, unsized on real maps; up to 3.3 in width by estimate | Every width from smoothed maps | One short test, then a model choice |
| One chain is enough (G9, G7) | Checked: a second chain was always computed and discarded (F42) | Every convergence statement | None; now written |
| `get_sky` of an index-field component (B3) | Checked: evaluated at the mean index (F43) | Chi-square diagnostic, fixed-component residuals, conditional likelihood, TOD-side sky | Small fix |
| Per-pixel amplitudes can be stored as alms (D14, F7) | Checked by the stats review on the E014 inputs: the stored maps cost 1.9 chi-square per amplitude inside the mask and stored samples keep 0.74 of their variance | Every stored amplitude sample of the per-pixel path; any acceptance through `global_loglike` (the legacy move, probably `conditional`) | Store maps; one confirming run |
| The legacy move as latent default (B8) | Checked by reading: fallback exists, not exercised | Nothing recorded | Make the lookup raise |
| Zero levels, calibration and bandpass are preprocessing (S3, D10, D7) | Known (F17), not in the model | Real intensity indices at the 0.03 to 0.1 level | Model components; weeks |
| Noise is diagonal per pixel and Stokes (D1, T6) | Untested | Large-scale polarization, first in the audit's ranking | A dense low-resolution covariance path |
| A CG solve stopped at 100 iterations is a sample (A1) | Untested | CMB and foreground amplitudes in the slow directions | One dense comparison at low nside |
| "No C(l) prior" is no prior (A2, B5) | Read: it is N(0, 1) per alm in units that depend on the reference frequency. Size unchecked | CG-path amplitudes; possibly the reference file's synchrotron at 408 MHz | A numpy check |
| One index for intensity and polarization; hier prior for intensity only (S2) | Untested | Any polarized index result | Design work (Q20) |
| Q and U smoothed as scalars (B2), window conventions (B7), single-precision chi-square (B15), seed reuse under `init_chain` (B16), under-determined pixels (B4) | Untested, each small or latent | Percent-level polarization effects; conditional likelihood at high resolution | Minutes each |
| TOD side: sigma0 a point estimate, stationary correlated-noise draw on gap-filled data, hard-coded gain prior, unchecked CGs (T1 to T7) | Untested, read briefly | Large-scale polarization in end-to-end runs | Cluster-scale; not ours to run now |

Two rankings follow from the audit and I agree with both. For calibrated indices on real
intensity maps the order is zero levels, the error model of the maps, bandpass, calibration,
ingest order. For large-scale polarization it is the noise covariance, the truncated CG, the
C(l) prior, the TOD noise steps, gain. **Neither list is led by a sampler.** The index samplers
validated in F12 to F38 are exact for a model whose inputs these items change.

## 5. Three things the review found that change the science

These came out of checking, not out of the plan, and they are the strongest new material.

**5.1 The legacy move has a law (E035).** The move is the "Monte Carlo within Metropolis" scheme
of the statistics literature: a pseudo-marginal kernel whose auxiliary variable, here the
amplitude draw, is refreshed between accept steps. Alone it targets the profile posterior with
amplitudes of half their posterior covariance. With the Gibbs refresh and one index step per
amplitude draw it converges to the profile posterior at temperature sqrt(2 pi k)/4, k the number
of re-solved amplitudes, so the index width is 0.79 k^(1/4) times too large. Toy chains confirm it
to 0.7 percent, two exact controls pass, and five earlier simulated measurements agree to 5
percent. A stats reviewer derived the same law without seeing the record and took it further:
the exact small-step form T = 1/2 + 1/(4 c_k), a record-process formula for several index steps
per amplitude draw, and the general form T = 1/2 + 0.443 times the sd of the chi-square at a
fresh draw. With that form the one loose end closed the same evening. The real-map width was 22
to 31 percent above the toy, and the chain's own log shows why: the chi-square the move sees has
2.44 times the exact variance, because per-pixel amplitudes are stored as alms and the stored
mean maps are not band-limited (F7). The toy with that factor returns the measured widths and
acceptance, and the gradient of the storage residual predicts the chain's +3 sd shift in the
dust index. The reviewer also found three misreadings in my record, now corrected. The
literature check found the temperature law in print for the sign rule (Ball, Fink and Bowler
2003, T = sqrt(pi/8) sigma; Bardenet, Doucet and Holmes 2017 with a slightly different
constant), so the paper cites them and claims the identification of the move, the chi-square
form and the k^(1/4) consequence.

**5.2 The flat amplitude prior is a different prior at every reference frequency (E036).** With
flat amplitude priors the index posterior shifts by ln(nu_ref / nu_pivot) times the index
variance of a single amplitude, where nu_pivot is the information-weighted mean frequency. The
shift does not average down over the pixels that share an index. At real WMAP/LFI polarization
depth, on a realistic amplitude field, it is -0.06 at a 23 GHz reference and +0.06 at 30 GHz for
30 degree regions with an index error of 0.06, and -0.40 against +0.34 for regions with an error
of 0.16. With a free CMB and a free dust amplitude in every pixel the pivot moves down to the
lowest band, so a flat prior stated at 23 GHz is nearly unbiased and one stated at 30 GHz is
biased flat by +0.17 at an index error of 0.064 and by +0.86 at 0.13. E022 saw the same sign and
size on real data at a 30 GHz reference (+0.27 in the Spur), and its recorded numbers obey the
formula: two noise models give the same pivot to 0.2 percent. With a CMB prior the pivot of the
BeyondPlanck band set is 29.5 to 30.2 GHz, which is the reference BeyondPlanck used. That is a
candidate explanation of why its chains show no flat-prior pull, the open point of F22. Published high-latitude values of
the polarized synchrotron index differ by 0.2 to 0.3 between analyses. This experiment does not
show that the reference frequency explains that spread, and no published analysis was checked
for it. It does show an effect of that size in the estimator that a Gibbs sampler with flat
amplitude priors targets (F15).

**5.3 The profile likelihood is sound, and a proper prior is its Bayesian form (E036, challenge
Sect. E).** For linear amplitudes with known noise the profile score is unbiased at any S/N. Its
width needs a sandwich factor sqrt(1 + 1/u), with u the mean amplitude S/N squared, which is 2
percent at the real-depth median. A zero-mean Gaussian amplitude prior with one free scale per
index region gives the same answers, needs no template and does not depend on the reference. In
E036 the two agree to 0.02 in the index in every class at real depth, with coverage 0.64 to 0.72
against a nominal 0.68, and they stay calibrated with free CMB and dust amplitudes. The profile
holds down to a per-element S/N of 1.1, four times the real noise, except in the faintest 30
degree regions, where its intervals are 1.4 times too narrow. That is the size of the sandwich
factor for a mean amplitude S/N near 1; the factor was not computed region by region here.
The free-scale prior shows a bias of 0.3 posterior sd in one class there, which a scale map
removes. The same analysis says what is unsound: a scale shared across index
regions and tied to a template. A scale map 30 percent high in the faint sky moves the faint
index by +0.12 or -0.12 depending on the reference in the challenge's two-zone toy. E036 then
tested it on the realistic field: the sign structure and the null at the pivot hold, and the
size is 0.03 for the regions the data inform at real depth. That is the pattern of F38 and of
chain 3, the displaced faint sky with a sign that flips with the reference, and it is a
candidate mechanism for them, not an established one.

The consequence for the program: the thirteen findings on template multipliers (F26 to F38)
were refinements of a prior whose mean is tied to a template. The sound estimators need no
template at all.

## 6. Impact weighting: where the next unit of effort buys most

Ordered by expected value per effort. "Evidence" is what says it will work.

| # | Work | What it buys | Evidence | Cost | Main risk | Verdict |
|---|---|---|---|---|---|---|
| 1 | Paper A1: the legacy move, the incidental-parameter problem (flat, profile with sandwich, proper priors), limits under template error, a Q/U section | A submittable methods paper from A and B grade evidence only | E035, E036, F2, F12, F15, F16, F18, F21, F23 | Six weeks, exact grids on one thread | A referee wants a positive real-data recommendation; the Q/U section supplies it | Do now |
| 2 | K band against 30 GHz null test with dense covariances: region indices from {K, Ka}, {K, 30}, {30, 44, Ka} | Decides whether Paper B is about an index field or about systematics; a result by itself | The discrepancy is published (Weiland et al. 2022, unverified extract); E022 infrastructure; the K-band map and covariance are fetched and pass their first checks (E038) | Two weeks, grids only | LFI low-resolution covariance quality; a few WMAP modes at multipoles 2 to 6 carry 20 to 90 times the median noise variance | Do now, before any field code |
| 3 | Report the two ingest bugs, the alm storage of per-pixel amplitudes, the unwritten chain, `get_sky` and the null-prior behaviour upstream; open the pull request for the index-sampler fix | Useful to every Commander4 user today; an outside check of our reading | E034 (measured), F42, F43 (read), the fix branch (tested) | Hours | None technical | Needs Trygve's yes |
| 4 | E032 (running): the template-prior model at WMAP/LFI white-noise depth, four chains | Closes the intensity template question, pass or fail | Scan predicts shifts below 0.02 | Already running | A fail leaves no at-depth evidence for the template prior, which the plan no longer needs | Finish; no follow-ups |
| 5 | The many-seed test of F8 (ridge target 1.9 and 3.7 sd from the input on the validation simulation) | Removes a possible bug under the harness that carries F2 | Referee and challenge both ask for it | An hour of grids; pre-registered as E037 and queued | It might be a bug | Running |
| 6 | Polarized index field (GMRF, NUTS) at nside 8 with an nside-2 field, on the gates of the challenge report Sect. F | Paper B | F24 (B plus), E036; feasibility: one factorization per second at nside 8 | Six to ten weeks after item 2 | The data may not determine a field: per-region error 0.09 at nside 2 against a sky scatter near 0.06 | After item 2, if its gate passes |
| 7 | tau from WMAP plus LFI as a paired difference between index models on the same data | The first number for how much the index model moves tau there | Dense likelihood exists once item 2 does | Two weeks after item 2 | Low power for the absolute value, which is not the claim | By-product of items 2 and 6 |
| 8 | Reader fix behind a switch and the E014 grid rerun (Q19) | A real-map demonstration that fits | E034 | Three days | D2, D3, D7, D10 remain | Only if the real-intensity section is kept |
| 9 | Per-band offsets, calibration factors and bandpass integration as model components | Needed for real intensity, and the same linear columns are what a K against 30 GHz model is made of | D013 decided, not implemented | Weeks in Commander4; days in the low-resolution code | None | In the low-resolution code as part of item 2 |
| 10 | Template-law refinements, the dust temperature field, sampler efficiency, the CG-path collapse, the translation move | Nothing that gates a paper | F39, F40 | Open-ended | Sunk-cost continuation | Freeze |

Item 10 is where two weeks went. The time box the challenge proposes is adopted: two follow-ups
on a failed pre-registration, then a design review before a third.

## 7. Way forward

**Paper A1 (six weeks).** Working title: spectral indices and incidental amplitudes in
pixel-based Bayesian component separation. Sections: (1) the model and its collapsed
likelihoods; (2) the legacy move as Monte Carlo within Metropolis, its temperature law, exact
alternatives; (3) the flat-prior volume term and its reference-frequency form, the profile with
its sandwich width, proper priors, with coverage against S/N from E012, E016, E017 and E036 and
the sandwich factor overlaid; (4) limits: a prior mean tied to a template (F21, F23, and F38 in
one paragraph), a scale shared across regions; (5) polarization at real depth: E036 with CMB and
dust, then real WMAP and LFI data in regions with dense covariance, flat against profile against
proper priors; (6) implementation and validation in Commander4. Not in it: the GMRF field as a
main claim, the real-intensity section, the template-multiplier sequence beyond one paragraph.

Remaining work for it, all exact grids unless stated:
1. E036 is complete (F44): the noise ladder to four times real noise, free CMB and dust, a CMB
   prior, off-centre truths, and the shared-scale arm at three references. Left for the paper:
   the sandwich factor per region, a truth with a latitude step, a data-built scale map.
2. The sandwich factor overlaid on the E012 and E017 coverage tables, and E017 at three
   reference frequencies.
3. The dense-covariance upgrade of E022 with K band on the four BP XV regions and on nside-2
   regions (shared with the null test below).
4. The many-seed F8 test (E037, queued); the Commander4 legacy mode with 1 and 20 steps per
   draw, amplitudes as maps and as alms (Q6, Q23; needs cores).
5. An SBC of the Commander4 `ridge` and `marginal` modes at small size, 300 replicates.

**The K band against 30 GHz null test (two weeks, overlapping).** Gate 0 of the challenge report
first (E038): the K-band res 4 map and inverse covariance are fetched and checked; still to do
are E022's conditional with dense matrices, an LFI covariance that passes a null test, and one
window for WMAP and LFI. Then region indices from
the three band pairs with the profile likelihood and its sandwich error, on WMAP9 with PR3, PR4
and DR1 30 GHz maps, and the residual projected on a few named modes (a relative gain, the WMAP
loss-imbalance modes, an LFI leakage template). Outcome A: the pairs agree within 2 combined sd,
and Paper B is an index-field paper on the gates of Sect. F. Outcome B: they do not and a
five-mode model absorbs it, and Paper B carries that model. Outcome C: neither, and the paper is
the three values and the statement that the public data do not determine the high-latitude
index. Each outcome is publishable.

**Paper B and tau.** Not planned in detail until the null test is read. Position against Dunkley
et al. 2009, which ran this architecture on WMAP5 with 48 index pixels. The polarized analysis
runs in a standalone low-resolution module with dense algebra, not in Commander4, which has no
Q/U covariance. That is a recorded decision for Trygve (D029).

**Commander4.** Queued fixes after E032: reader order and floor behind a switch, `get_sky` with
an index field, the legacy lookup raising on a missing key. Audit tests 4 to 8 and 10 when cores
are free.

**Process.** D028, plus from the challenge: the regime gate also names the components a
simulation lacks, its noise correlation, its window and how its template was built; an effective
sample size floor of 400 beside the R-hat bound; coverage statements with a replicate-clustered
error below 0.03; no truth at the centre of a grid; literature numbers checked in the PDF before
they gate anything.

## 8. Decisions for Trygve

1. **Adopt the plan of Sect. 7** (D029, proposed): Paper A1 as scoped, the null test next, Paper
   B decided on its outcome, intensity template work frozen after E032.
2. **Where the GMRF index field goes.** In Paper A1 as an outlook section with F24 and E032, or
   wholly in Paper B. I recommend Paper B: it keeps A1 on exact-grid evidence.
3. **Upstream.** Open the pull request for the index-sampler fix, and report the ingest bugs, the
   alm storage of per-pixel amplitudes, the null-prior behaviour, the unwritten chain and
   `get_sky` to the Commander4 developers. Also
   whether to tell the Cosmoglobe team about the DR1 dust reference frequency (F13) and the
   Commander3 code facts (F10, F11). Nothing has been sent.
4. **The polarized analysis outside Commander4**, as a standalone low-resolution module in this
   repository.
5. **The real-intensity section of the paper.** Cut it, or keep it and spend three days on the
   reader fix and a rerun.
6. **Data.** The WMAP K-band res 4 products (150 MB) were fetched on 2026-11-09. Planck
   low-resolution covariances (a few GB) will be fetched without asking, being far below the
   200 GB limit.

## 9. The assumptions this plan rests on, and how each is tested

| Assumption | Status | Test | Result or date |
|---|---|---|---|
| The profile likelihood and the per-region zero-mean prior are calibrated for Q/U at real depth | Tested for synchrotron only, white noise | E036 | Pass at real depth; harder cells below |
| They stay calibrated with free CMB and dust amplitudes | Tested | E036 second stage | Pass (E036 second stage): mean z within 0.08, sd of z 0.91 to 1.08, coverage 0.62 to 0.74; the region index error grows by 1.63 |
| Their calibration does not come from the truth sitting at the grid centre | Tested | E036 with truths -2.8 and -3.4 | Pass: nothing moves below a region index error of 0.25; above it the posterior mean reverts 13 percent toward the grid centre |
| A scale shared across regions on a template is biased with a reference-dependent sign | Shown in a two-zone toy by the challenge; tested | E036 with a shared-scale arm at 23 GHz, the pivot and 30 GHz | Sign and null confirmed, size small: a scale map 30 percent high at high latitude moves the high-latitude index by -0.026, 0.000 and +0.029 |
| The CMB prior is the prior that matters in polarization | Shown on the E022 surrogate with diagonal noise | E036 with and without a CMB prior; the dense upgrade of E022 | Pass at real depth (E036): with a 0.2 uK CMB prior the region index error equals the synchrotron-only value (0.059 against 0.060); with a flat CMB prior it is 1.63 times larger and the pivot of the flat-prior shift moves to the lowest band. The dense upgrade of E022 is still to do |
| The WMAP K, Ka, Q, V inverse covariances are usable as released | K and Ka checked (E038): symmetric, positive definite, pixel rms as estimated; a few modes at multipoles 2 to 6 carry 20 to 90 times the median variance | Q and V through the same check; E022's conditional with dense matrices | Week 1 |
| A trustworthy LFI 30 GHz low-resolution covariance exists publicly | Unverified | Gate 0 half-split null test | Week 1 |
| K band and 30 GHz can be told apart from the index | Doubtful: they carry 91 percent of the index information and only Ka separates a relative gain from the index | The null test; Fisher with a gain and large-scale modes added | Weeks 1 to 3 |
| The data determine an index field, not only a mean and a latitude trend | Marginal. Hand estimate from the E036 region errors with white noise: the variance of the field is detectable at 2 to 3.5 sigma for a sky scatter of 0.06 and at 3.5 to 6 sigma for 0.10 to 0.15, with 48 or 192 regions. A correlation length is hardly constrained. Real large-scale noise is worse (E038) | The exact posterior of the field amplitude and range for a field seen with those errors and the real covariances | With the null test |
| Dense NUTS is affordable | Supported at nside 8 (1 s per factorization), not at nside 16 (a minute) | Gate 5 | Timed by the challenge |
| The legacy-move law describes Commander4 | Supported on simulations and, with the storage noise measured from the chain's log, on real maps to 3 percent | A Commander4 run with amplitudes kept as maps (Q6, Q23); predicted widths 13.1 and 13.6, no dust shift | After E032 |
| The validation harness has no bug behind F8 | Untested | E037: 40 seeds, exact grids, a white-noise control | Queued behind E036 |
| The literature numbers that gate decisions are right | Unverified extracts | Five or six PDFs to read: Weiland 2018 and 2022, Dunkley 2009, QUIJOTE XII, BP XV, Ball et al. 2003 | Before they are cited |

## 10. Limits of this review

1. One session, one machine, with a production chain holding the cores, so every new experiment
   here is a single-thread grid or a toy.
2. E032 is unread. Its outcome changes one paragraph of Paper A1 and nothing in the plan.
3. The referee report is incomplete. The challenge report and this document cover its missing
   sections in part.
4. No literature statement in this document was checked against a PDF by the main session.
5. The reference-frequency shift is confirmed by one experiment (E036). The temperature law is
   derived twice and confirmed by two toys and an exact-kernel grid. The sandwich factor and the
   two-zone bias of a shared scale come from the challenge report's algebra and toy, which no
   experiment record holds yet. The explanation of the real-map legacy width re-implements the
   storage round trip and has not been confirmed by a Commander4 run.
6. E036 uses white noise, one index for the whole sky, and a truth amplitude field that is a
   posterior mean of real data and therefore carries noise of its own.

