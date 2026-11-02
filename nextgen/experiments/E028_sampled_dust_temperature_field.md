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
Stage A and B (queued behind E027 by `run_e028_queue.sh`, 2026-11-02 19:05; parameter files
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
(pending)

## Follow-ups
