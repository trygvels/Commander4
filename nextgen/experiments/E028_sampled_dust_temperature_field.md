# E028: A sampled dust temperature field (log T_d per region) in the collapsed hier index-field likelihood

- **Date started / finished:** 2026-11-02 /
- **Branch and commit:** nextgen, 4282505 (design), code from 71c3be7 (SED) and fbcfb6c (layout, likelihood, sampler, component state)
- **Open question addressed:** Q16 (F30: the per-region multiplier field turns the fixed-T misfit into a bright-sky dust-index bias), M2.6 (sampled T_d per region), Q3 (joint fields for beta_d and T_d)
- **Validation level:** L0 to L1 for the code (finite-difference and exact-target checks), then L4 on the PySM noise-1 sky against the fixed-T chain of E026

## Question
Does sampling the dust temperature as a per-region field remove the bright-quartile dust-index
offset (+0.06) and the 3 percent multiplier offset of the fixed-T chain at realistic noise (E026,
F30) without degrading the synchrotron side, and is the joint (beta_d, T_d) field posterior
well-behaved (the beta-T degeneracy per region) at nside 16 with Planck plus WMAP intensity?
If yes, the hier prior with per-region multipliers is admissible for dust on real data; if the
degeneracy makes T prior-dominated everywhere, region-wise T at a coarser nside or a fixed
measured T map is the fallback.

## Target and priors
Everything as in E026 gmrf16_exact_mumap_nc (its Target and priors: 12 bands 23 to 353 GHz at
nside 64, per-pixel amplitudes (a_cmb, a_s, a_d) marginalized analytically with the flat CMB
prior N(0, 1e4^2) and the hier priors a_c ~ N(mu_r g_c, (s_c g_c)^2) with exact templates and
one multiplier per nside-16 region, hyperpriors mu_r ~ N(1, 0.5^2), log s_c ~ N(log 0.3, 1);
GMRF index fields for beta_s and beta_d at nside 16, non-centred, graph-Laplacian basis,
alpha 2, sigma_scale 0.3, rho loguniform, field-mean hyperpriors N(-3.1, 0.3^2) and N(1.6,
0.3^2)), with one change: the dust SED in pixel p of region r is

  sed_d(nu; beta_r, T_r) = (nu / nu_ref)^(beta_r + 1) expm1(h nu_ref / k T_r) / expm1(h nu / k T_r),

with T_r = exp(tau_r) and tau a third GMRF field on the same nside-16 regions (first version:
the same nside as beta_d; a coarser T nside is a later extension), non-centred,
tau = mu_tau + sigma_tau V z_tau, with the same basis and alpha, hyperpriors mu_tau ~ N(log 20,
0.2^2), sigma_tau half-normal with scale 0.15 (the HyperPrior convention used for beta), rho
loguniform on the same range. The sampled coordinates of the index-field group become
(z_s, mu_s, log sigma_s, u_s; z_d, mu_d, log sigma_d, u_d; z_tau, mu_tau, log sigma_tau, u_tau;
mu maps and log s of the two hier priors). The collapsed likelihood is unchanged in form; its
gradient with respect to tau_r is the kernel's per-column gradient with the derivative column
d sed_d / d tau = T d sed_d / dT = sed_d [x e^x / expm1(x) - x0 e^x0 / expm1(x0)], x = h nu / k T,
x0 = h nu_ref / k T (stable form x / (1 - e^-x)). The CG amplitude step and the band projections
use the region T map through the same mixing matrices (the fixed scalar T stays the parameter
file's value when T is not sampled; default behaviour unchanged).

## Stats review of the target (2026-11-02, independent stats-reviewer, before any chain)
Verdicts: target and derivative sound with three specification gaps (now fixed above and
here); validation plan insufficient as first written; pre-registered reading not falsifiable
in part and, decisively, the experiment as designed cannot answer its question on this data set.
- Log T versus T: the prior is Gaussian in tau, so the induced prior on T_r is log-normal
  (median exp(mu_tau), one-sd range 16.4 to 24.4 K at the field mean); no Jacobian is needed
  in tau coordinates. The field is tau = mu_tau + sigma_tau V diag(s_i(rho)) z, the Matern
  weights of the index field (`gmrf_field.py`), not an iid field.
- Degeneracy: with s_d at the guard the dust model per region is mu_r g_p sed(nu; beta_r, T_r)
  with three free numbers, and mu_r is the amplitude at 545 GHz, outside the 23 to 353 GHz
  bands. Reviewer's Fisher analysis on the E023 sky (true amplitudes, CMB profiled, sync
  pinned; its known-amplitude beta sds reproduce `fisher_sd_dust_16`): corr(log mu_r, tau_r)
  +0.91, corr(beta_r, tau_r) -0.97; Fisher sd of T_r with (mu_r, beta_r) free 52 / 310 / 852 /
  3141 K per quartile (medians), minimum 2.75 K over 3072 regions, against a prior sd of
  about 3 K from sigma_tau 0.15: **at nside 16 the temperature field is prior-dominated in
  every region of this sky.** A sky-wide scalar tau with per-region (mu, beta) free has sd(T)
  0.7 K; at nside 4, 5 of 192 regions are below 2 K; adding 545 and 857 GHz (rms 1/2 and 1/4
  of the 353 GHz rms in uK_RJ) gives 157 of 3072 regions below 2 K at nside 16 (corr still
  -0.97). A 10 percent T error moves sed(353)/sed(545) by 3 percent: the F30 multiplier
  offset of 1.029 is this mechanism.
- Validation: (i) a T-independent target (a key mismatch or the cached scalar-T SED) would pass
  a finite-difference test with 0 = 0; a sensitivity check and the independent covariance-form
  target at several T states are required (added: `test_region_index_likelihood_with_temperature_field_matches_covariance_form`).
  (ii) The amplitude step (`mixing_regions`, `get_sed_map`, `needs_region_solver`) and the
  index step build the mixing matrix separately; if T reached only one, the Gibbs chain would
  have no joint target (added: the state round-trip test checks `mixing_regions` against the
  component's SED at the T map, and the covariance-form test checks a fixed `T_map` against the
  sampled field at the same values). (iii) The toy grid must be the numpy marginal times the
  priors, not the code's own density. (iv) The most informative additional check is an
  exact-model recovery run at the production nside (sky drawn from the model, full Gibbs
  chain, coverage of (T_r, beta_r, mu_r) and the D025 check along tau).
- Reading: criterion (c) was vacuous (no region with Fisher sd of T below 1 K); the beta-only
  criterion (a) against `bd_trueT` measures the tau prior centre when T is prior-dominated
  (corr -0.97); both chains must be read against the same reference.

**Design after the review (D026).** The question "does a sampled T field remove the F30 dust
offsets" cannot be answered by a prior-dominated field with 23 to 353 GHz. Three stages:
- **Stage A (oracle fixed T map; runs first, after E027):** `T_map` = the true T_d region means
  at nside 16 (`/Users/trygvels/data/e023_pysm_n1/truth.npz` T_d, ud_grade to 16), everything
  else as E026 gmrf16_exact_mumap_nc. Tests whether the fixed-T misfit is the cause of the
  +0.06 bright-quartile dust offset and the 1.029 multipliers. Pre-registered: both the
  fixed-T chain and this chain read against `bd_trueT` (region fit at the true per-pixel T) in
  the bright dust quartile; the oracle chain's z median must be within 0.5 and the mu_d map
  mean within 0.01 of 1; if the offsets remain, the fixed T is not their cause and F30's
  mechanism statement is withdrawn.
- **Stage A' (added 2026-11-03 01:52 after E029, F33):** stage A's oracle T map with the
  E029 multiplier hyperprior N(1, 0.1^2) on both components
  (`param_gmrf16_exact_mumap_nc_Tmap_mu01.yml`, queued after stage B by
  `run_e028_stageAprime.sh`). E029 showed the tight hyperprior fixes the sync side and
  exposes the fixed-T dust misfit; this chain tests the combination. Pre-registered (against
  E029 as the control): sync field mean and faint quartiles unchanged from E029 (within
  0.03); dust bright quartile z median within 0.5 of `bd_trueT` and every quartile within
  0.03 of `bd_trueT` (E029: -0.048, -0.017, -0.022, -0.045); mu_d map 1.00 +- 0.01; the
  fixed-T GLS CMB residual at or below the E026 control's 93.6 uK (E029: 94.5). If the faint
  dust quartiles stay 0.04 below `bd_trueT` with the true T, the shrinkage is a field-prior
  effect and not the T misfit.
- **Stage B (sampled scalar log T):** `sample_T: true`, `T_nside: 0`, prior N(log 20, 0.2^2):
  identifiable (sd 0.7 K) and the production exercise of the code path with the D025 check
  along tau. Pre-registered: posterior T within 2 sd of the Fisher-weighted truth T (to be
  computed from the truth file with the dust amplitude weights); the Gibbs-consistency check
  passes; a scalar T cannot remove a bias from spatial T variation, so (a) and (b) are not
  criteria here.
- **Stage C (sampled T field):** needs 545 and 857 GHz in the simulation (a new E023-type sky
  `e023_pysm_n1_hfi` with 14 bands) and T at nside 4 to 8; the exact-model recovery run is the
  L3 rung before it. Deferred until E027 and stage A are read. For real data this is the
  configuration that matters (Planck 545 and 857 GHz exist in intensity).

## Setup
Code: `ThermalDust.get_sed(nu, beta, T)` and `get_dsed_dT`; a per-(component, parameter)
`IndexSpec` with prefixed parameter keys (`sample_T`, `T_nside`, `T_prior`, `T_field`,
documented in `params/param_explanation.yml`); `StateLayout` keyed by (component, parameter);
one derivative column per parameter in `_stream_terms` and in `mixing_regions` (kernel called
once per parameter, option (b) of the code map
`nextgen/research/2026-11-02_sampled_Td_field_code_map.md`, unless profiling shows the second
call matters); chain output nested by parameter (`index_field/beta_map` kept for beta,
`T_field/` for T). Data: `/Users/trygvels/data/e023_pysm_n1`. Run: `param_gmrf16_exact_mumap_nc_T.yml`
(beta_d and log T_d fields at nside 16), 1200 iterations, warm-up 400, burn 500, dense metric.

## Reference
1. L0: finite-difference tests of `get_dsed_dT` and of the likelihood gradient with respect to
   the tau coordinates (extension of `test_region_index_likelihood_gradient_matches_finite_differences`).
2. L1: the numpy per-pixel 12 x 12 marginal of `e023_exact_scan.py` extended with a per-region
   T must agree with `RegionIndexLikelihood` at random (beta, tau, hyperparameter) states to
   1e-9 relative, as in the E023 kernel check; and on a 12-pixel, 2-region toy sky the NUTS
   posterior of (beta_r, tau_r) with fixed hyperparameters must match a grid posterior
   (means within 0.05 sd, sds within 5 percent; the E012 L1 pattern).
3. L4: on the PySM noise-1 sky against E026 gmrf16_exact_mumap_nc (fixed T): (a) the dust
   index in the bright quartile against the true-T region-fit reference (`bd_trueT` in
   `E027/refs_nside16.npz`) with z median within 0.5 (fixed T: +1.4 against the fixed-T
   reference); (b) the dust multiplier map mean within 0.01 of 1 (fixed T: 1.029 +- 0.007);
   (c) the T field against the region mean of the PySM T_d (sd 2.15 K over the sky, 1.7 at
   |b| > 20) per dust S/N quartile: where the region Fisher sd of T (to be computed, amplitudes
   profiled) is below 1 K the posterior must cover the region mean (|z| < 2 for 90 percent),
   elsewhere it may be prior-dominated and must say so (posterior sd at the field scale);
   (d) the sync field unchanged within its posterior sd per quartile; (e) the Gibbs-consistency
   check (D025) along all three field means, and the roughness check for the tau field.
Pre-registered risk: the beta-T degeneracy per region is strong below 353 GHz; the posterior
may put T on the prior in the faint half and let beta carry the fit, in which case (a) and (b)
decide whether the bright sky, where the data constrain both, is repaired.

## Command
Stage A and B (queued behind E027 by `run_e028_queue.sh`, 2026-11-02 17:49; parameter files
differ from `param_gmrf16_exact_mumap_nc.yml` only in the output directory and the ThermalDust
lines shown):
```bash
# stage A: T_map: "/Users/trygvels/data/e023_pysm_n1/truth_Td_nside16.fits" (true T_d region means
#          at nside 16: mean 20.93 K, sd 2.05, range 15.5 to 26.5)
mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_gmrf16_exact_mumap_nc_Tmap.yml
# stage B: sample_T: true, T_nside: 0, T_prior gaussian mean log 20 rms 0.2
mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_gmrf16_exact_mumap_nc_Tscalar.yml
```
Stage B's pseudo-true scalar T is not the sky mean: the unweighted truth T is 20.93 K, the
a_dust^2-weighted 22.64 K, the |a_dust|-weighted 20.10 K; the value a scalar T estimates under
the pinned-amplitude model is computed before the chain is read (noise-free fit of one T with
per-region (mu_r, beta_r) free, the reviewer's construction), and the 2-sd criterion refers to
that value.
Grid scan (`E028/Tgrid_scan_global_and_profiled.txt`, 18:17, expected mode, region fit redone
at every T): under the global law (one multiplier per component) the optimum is at 22.1 K with
a curvature sd of 0.8 K (the a_dust^2-weighted truth is 22.6 K; the dust index fit drops from
1.600 at 18 K to 1.478 at 24 K, the beta-T degeneracy). The per-region-profiled rows of that
file are invalid (F31: the profiled multipliers prefer the low end monotonically); the
pseudo-true value under the chain's own law (multipliers integrated, `perregmarg`) is being
computed (`Tgrid_scan_marginal.txt`) and is the stage B reference. Result (18:36,
`E028/Tgrid_scan_marginal_sd05.txt`): under the chain's law the expected target rises
monotonically over 19 to 24 K (7.3 nats from 19 to 20 K, 1.8 nats from 23 to 24 K, still rising
at 24 K): the scalar T is not pinned near the sky's 20.9 K by this target; the marginalized
multipliers and the re-fitted index field trade against T (the beta-T-multiplier degeneracy of
the review, now in the target's optimum). The grid is being extended to 34 K and repeated with
the sd-0.1 hyperprior. Stage B prediction, written before the chain: under sd 0.5 the scalar T
posterior sits well above 22 K (where exactly, the extended grid says), with the dust index
field lower than the truth in step (the region fit gives 1.478 at 24 K against 1.531 at 20.9 K);
this is the identifiability failure of the review in its scalar form, and stage B is read as a
measurement of it, not as a recovery test. Extended grids (19:22; `Tgrid_scan_marginal_sd05_extended.txt`,
`Tgrid_scan_marginal_sd01.txt`): under the sd-0.5 multiplier hyperprior the optimum is at 27 K
(curvature sd about 1.7 K; totals -83.97, -82.38, -82.50, -83.65, -87.54 at 24, 26, 28, 30, 34 K),
a +6 K Occam push of the marginalized multipliers on the temperature, with the region-fit dust
index at 1.44 against a truth of 1.54; under the sd-0.1 hyperprior the optimum is at 20.2 K
(sd 0.5 K; -598.1, -597.0, -604.1 at 19, 20.5, 22 K), at the sky's 20.9 K within 1.5 sd. Stage
B's pseudo-true scalar T is therefore 27 +- 2 K under the queued parameter file (sd 0.5) and
the chain is predicted there; a stage B run with `mu_prior [1, 0.1]` (the E029 setting) is
predicted at 20.2 +- 0.5 K and is the configuration that can recover T. Decision: stage B runs
with the sd-0.1 hyperprior (parameter file updated before the queue reaches it), the sd-0.5
version is dropped (its prediction is recorded and does not need a chain).
```bash
python -m pytest -q tests/test_components.py tests/test_index_field.py tests/test_index_field_sampler.py tests/test_hier_group_wiring.py
mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_gmrf16_exact_mumap_nc_T.yml
```

## Result
Baseline for stage A (fixed-T chain against both references, E027 quartile labels):
`E028/fixedT_baseline_vs_trueT.txt`, computed before the oracle chain ran. The fixed-T chain's
bright dust quartile is already at the true-T reference (z median -0.06, mean +0.006 against
`bd_trueT`; +1.46, +0.060 against `bd_fixT`): the +0.06 of F30 is beta_eff minus beta_d, and
the chain with free multipliers estimates beta_d. Stage A's prediction is therefore sharpened
before it runs: the oracle T map should leave the bright-quartile index where it is (z within
0.5 of bd_trueT, as now) and move the mu_d map from 1.029 to 1.00 +- 0.01; the fainter
quartiles (+0.06 against bd_trueT, z +1.0 to +1.2) are expected to stay, being a field-prior
effect (E027). If the mu_d map stays at 1.03 with the true T, the multipliers absorb something
other than the T misfit and F30's mechanism is wrong.

## Negative control
The fixed-T chain (E026 gmrf16_exact_mumap_nc) is the control: its +0.06 / 1.029 must be
reproduced by the new code with `sample_T: false` (regression) before the sampled-T run is read.

## Conclusion
**Stage A (oracle T map, sd 0.5), read 2026-11-03 08:10 to 08:19.** Chain
`gmrf16_exact_mumap_nc_Tmap`, 1200 iterations in 6.4 h (01:43 to 08:05), 700 post-burn draws,
63 gradients per iteration, 2 divergences, no tracebacks. Files `E028/stageA_*.txt`.

*Validity.* The chain used the T map and is a valid draw from its target: the D025
conditional-CDF check against the T-map target (K = 100, `--t-map` added to
`e026_conditional_check.py`) passes, sync KS p 0.12 (mean CDF 0.450, optimum +0.20
conditional sd), dust p 0.39 (mean 0.512, optimum -0.01 sd); against the fixed-T target
(control, K = 40) the same draws fail, dust KS D 0.98 p 0.000 (optimum -5.2 conditional sd),
sync p 0.000, and the log target at the first draw is 1269 nats lower (1274181 against
1275450). The region-mean T map is not a coarse stand-in: the region fit at the map T differs
from the fit at the true per-pixel T by -0.003 in the dust index (rms 0.008; `stageA_scan_*`
reference lines), so `bd_trueT` is the reference for this chain.

*Pre-registered criteria: both fail.* The mu_d map mean is 1.0254 +- 0.0073 (required 1.00 +-
0.01; fixed-T chain 1.0286 +- 0.0070) and the bright dust quartile sits at z median +1.31
against `bd_trueT` (required within 0.5; fixed-T chain -0.06). By the pre-registration the
fixed T is not the cause of the 1.029 multipliers, and F30's mechanism statement ("the
multipliers absorb the fixed-T misfit at the 545 GHz pivot") is withdrawn. The oracle T
makes the dust index worse, not better:

| dust, against `bd_trueT` (mean offset / z median) | Q1 bright | Q2 | Q3 | Q4 faint | field mean | pixel RMSE vs truth |
|---|---|---|---|---|---|---|
| fixed T, sd 0.5 (E026) | +0.006 / -0.06 | +0.061 / +1.17 | +0.059 / +0.97 | +0.030 / +0.46 | 1.578 +- 0.012 | 0.060 |
| oracle T map, sd 0.5 (stage A) | +0.047 / +1.31 | +0.095 / +2.17 | +0.098 / +2.09 | +0.066 / +1.37 | 1.616 +- 0.010 | 0.085 |
| fixed T, sd 0.1 (E029) | -0.048 / -1.55 | -0.017 / -0.29 | -0.022 / -0.30 | -0.045 / -0.67 | 1.506 +- 0.011 | 0.056 |
| truth (sky mean 1.540) | | | | | | |

Against each model's own reference (`bd_fixT` for the fixed-T chain, `bd_trueT` for the
oracle chain) the bright quartile is +0.060 and +0.047: the sd-0.5 chain sits about +0.05
above its own pseudo-true index whatever the T model, and the fixed-T chain's agreement with
`bd_trueT` (the F30 correction of 2026-11-02) was a coincidence of two offsets of opposite
sign (beta_eff3 minus beta_d is -0.054 in Q1: `bd_fixT` 1.515, `bd_trueT` 1.570). The multiplier offset lives in the bright half
(`stageA_mu_map_quartiles.txt`: dust mu map per quartile 1.038, 1.041, 1.017, 1.005 against
1.062, 1.038, 1.012, 1.002 fixed-T; posterior sds 0.10, 0.34, 0.45, 0.49): the true T halves
the bright-quartile multiplier offset (1.062 to 1.038) and raises the index there by 0.04.
In the sd-0.1 chain the map is 1.012, 1.000, 1.000, 1.000.

*Sync unchanged (criterion (d) met).* Field mean -2.880 +- 0.042 against -2.868 +- 0.042,
quartile offsets against the reference +0.06, +0.11, +0.12, +0.12 (fixed T: +0.05 to +0.14),
sigma_s 0.122 against 0.152, pixel RMSE 0.133 against 0.147; sync mu map 0.994.

*CMB.* The flat-prior GLS CMB residual at the posterior-mean fields is 92.45 uK with the
fixed-T GLS and 92.31 with the T-map GLS (fixed-T chain 93.60, sd-0.1 chain 94.52, true-SED
floor 93.18): 1.12 +- 0.07 uK below the fixed-T chain and 0.69 +- 0.13 uK below the floor.
A biased index field beats the floor because the GLS with a wrong mixing matrix trades bias
for variance; replacing the chain's dust field by `bd_trueT` everywhere raises the residual
to 94.38 (`stageA_cmb_attribution_Tmap.txt`). The GLS CMB residual rms therefore rewards a
steeper dust index and is not an index-quality criterion. **Amendment to the stage A'
pre-registration (08:19, before A' runs):** the criterion "GLS CMB residual at or below the
E026 control's 93.6 uK" is dropped as a quality criterion; the residual is reported, and the
index criteria against `bd_trueT` decide.

*Mixing.* Dust region ESS min 93, median 765 (fixed T: 21, 136); dust hyperparameters ESS 39
to 120; sync region ESS min 179.

*Exact-target scans (`stageA_scan_mapT_vs_fixT_sd05.txt`, `_sd01.txt`; the F29 machinery with
a `mapT` model added to `e023_exact_scan.py`, multipliers integrated, scatters profiled at
the reference, dust quartile shifts against `bd_trueT`; 08:19 to 08:27):* the expected-target
optimum per dust quartile (Q1 to Q4; curvature sd in Q1 0.003 to 0.004, in Q4 0.05 to 0.06)
and the whole-sky sync shift:

| model | d_Q1 | d_Q2 | d_Q3 | d_Q4 | d_all | s_all |
|---|---|---|---|---|---|---|
| sd 0.5, fixed T | +0.02 | +0.10 (edge) | +0.09 | +0.08 | +0.03 | +0.10 (edge) |
| sd 0.5, region T map | +0.01 | +0.10 (edge) | +0.09 | +0.07 | +0.01 | +0.08 |
| sd 0.1, fixed T | -0.01 | -0.03 | -0.01 | +0.01 | -0.01 | +0.03 |
| sd 0.1, region T map | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | +0.01 |
| chain, sd 0.5, fixed T (E026) | +0.006 | +0.061 | +0.059 | +0.030 | | |
| chain, sd 0.5, region T map (stage A) | +0.047 | +0.095 | +0.098 | +0.066 | | |
| chain, sd 0.1, fixed T (E029) | -0.048 | -0.017 | -0.022 | -0.045 | | |

Reading. (i) The F29 Occam push of the marginalized multipliers acts on the dust index under
sd 0.5 with either T model: +0.07 to +0.10 in the three fainter quartiles (the logdet term
rises with a steeper index; in Q2 the gain is 12 nats at the +0.10 edge), and this is what
the stage A chain shows in Q2 to Q4 (+0.07 to +0.10). The fixed-T chain is held back to
+0.03 to +0.06 by the field prior, as the sync field was in F29. (ii) The bright quartile's
pseudo-true offset is +0.01 to +0.02 under sd 0.5; the chain's +0.047 with the true T is the
GMRF field prior transmitting the faint-sky push (field mean 1.616 against a Q1 level of
1.57; region posterior sd 0.039, field sigma 0.046), while with the fixed T the field mean
1.578 sits at the Q1 level and there is no pull. (iii) Under sd 0.1 with the true T map the
exact marginal has its optimum at zero in every dust quartile and +0.01 for the sync level:
**the stage A' prediction is the dust index at `bd_trueT` in every quartile and the sync field
as in E029.** Under sd 0.1 with the fixed T the optimum is -0.01 to -0.03 (the fixed-T misfit,
read against beta_d), consistent with E029's -0.02 to -0.05 after the field prior. Both
F34's mechanism and the A' prediction are therefore written before A' runs (A' is queued
behind stage B).

## Follow-ups
