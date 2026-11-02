# E028: A sampled dust temperature field (log T_d per region) in the collapsed hier index-field likelihood

- **Date started / finished:** 2026-11-02 /
- **Branch and commit:** nextgen, 4282505 (design; no code yet)
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
```bash
python -m pytest -q tests/test_components.py tests/test_index_field.py tests/test_index_field_sampler.py tests/test_hier_group_wiring.py
mpirun -n 12 commander4 -p /Users/trygvels/data/e023_pysm_n1/param_gmrf16_exact_mumap_nc_T.yml
```

## Result
(design stage)

## Negative control
The fixed-T chain (E026 gmrf16_exact_mumap_nc) is the control: its +0.06 / 1.029 must be
reproduced by the new code with `sample_T: false` (regression) before the sampled-T run is read.

## Conclusion
(pending)

## Follow-ups
