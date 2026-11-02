# Code map: sampling a per-region dust temperature field alongside beta_d (agent report, 2026-11-02)

Explore-agent reading of the nextgen code at commit 4282505; verified by the author only at the
points marked (checked). Nothing was modified or run.

**Summary.** Every layer of the hier/NUTS path assumes one sampled parameter per component,
named `beta`, and one region nside per sampling group. The C++ kernel is the smallest part of the
change (its maths handles any parameter that touches a single mixing column; its output has one
gradient slot per component). Nothing supports a per-pixel or per-region T today, and
`ThermalDust` has no dT derivative (checked: `src/commander4/sky/diffuse_components.py:645-667`).

## 1. Discovery and layout
- `discover_index_specs` (`compsep/index_field_sampler.py:254-285`): gates on `hasattr(comp,
  "beta")` and `sample_spectral_index`; reads `spectral_index_nside`, `spectral_index_prior`;
  one `IndexSpec(comp.comp_name, ...)` per component; `views` keyed by component name.
  Change: loop over a per-class list of sampleable parameters (e.g. `("beta", "T")`), prefixed
  parameter keys, one spec per (component, parameter).
- `IndexSpec` (`compsep/index_field.py:340-372`): no parameter name. Add `param` and a key
  `(comp_name, param)`.
- `StateLayout` (`index_field.py:375-441`): raises on a repeated component (382-383); slices
  keyed by component (388, 392, 408-418, 433-436); one nside for the whole group (398-403),
  scalars broadcast to it. Change: key by (comp, param); a fine common nside with per-spec
  parent maps if T is to live at a coarser nside than beta.
- Hier hyperparameter slices (`hier_component_dims`, `hier_slices`, `hypers()`) are per
  component and need no change.
- `RegionIndexLikelihood`: per-spec `hasattr(comp, "beta")` (509-512); `region_map` at
  `layout.nside_beta` (532, 744); `g_beta` dicts per component (573, 648); `_stream_terms`
  (597-605) calls `get_sed(freqs, beta=...)` and `get_dsed_dbeta` once per component with one
  dM column; gradient accumulation per component (617-619); the allreduce packs
  `layout.nreg` per spec (656, 667-669). Change: one dM column per parameter, per-spec gradient
  lengths.
- `initial_state` (405-445) reads `comp.beta` / `comp.beta_map`; `apply_state` (448-461) writes
  `comp.beta` or `set_beta_map(..., layout.nside_beta)`; `unpacked()` (619-645) writes
  `index_field/{beta_map, beta_nside, mu, log_sigma, u, z}` per component (a second parameter
  would overwrite). Component chain I/O (`sky/diffuse_components.py:338-340, 365-376`) writes
  and reads only `index_field/beta_map`, `beta_nside`.
- Validators: `processing.py:364-383` (`_validate_index_fields`), 386-395, 339-361;
  `standalone_tools/validate_params.py:121-128` whitelists `spectral_index_*` for ThermalDust;
  the MH path `spectral_index.py:114-118` gates on `beta`.

## 2. SED, derivatives, kernel
- `ThermalDust.get_sed` (645-662) uses scalar `self.T`; `_broadcast_nu_beta` (441-448)
  broadcasts only beta. `get_dsed_dbeta` (664-667) is `sed ln(nu/nu_ref)`. No `get_dsed_dT`
  exists. With x = h nu / k T and x0 at nu_ref: d sed / dT = sed / T [x e^x / expm1(x) - x0
  e^x0 / expm1(x0)], stable form x / (1 - e^-x).
- `mixing_regions` (`index_field.py:311-333`) reads only `beta_map` and gives one dM column per
  component; `augment_mixing` (262-272) assumes dM has M's shape; `field_nside` (283-290)
  reads `beta_nside` only; `needs_region_solver` (275-280) and `_band_sed` / `get_sed_map`
  (`diffuse_components.py:468-485`) trigger only on `beta_map` (a T map would be silently
  ignored by the CG/SHT projections 504-567 unless added).
- C++ `Py_marginal_loglike_regions` (`src/lib_cpp/compsep_pymod.cc:319-496`): `dM_reg (nreg,
  nband, ncomp)` asserted equal to `M_reg`'s shape (366-367); doc at 324-325 states one
  parameter per component; gradient `(nreg, ncomp)` (371, 451-477); the per-parameter formula
  only needs the column c the derivative belongs to. Options: (a) accept `dM_reg (nreg, nband,
  npar)` plus `par_comp (npar,)`, return `(nreg, npar)`; (b) call the kernel once per parameter
  (same likelihood value; doubles the cost).

## 3. Field priors
- `field_prior_from_params(comp_params, nside_beta, centred)` (`index_field_sampler.py:203-251`)
  reads `spectral_index_prior` and `spectral_index_field` from the component block; builds
  `IndependentIndexPrior`, iid `IndexFieldPrior`, or `MaternField(graph_laplacian_basis(nside),
  alpha)` with `HyperPrior` (sigma_scale default 0.3, rho loguniform, rho_max pi); capped at
  `NSIDE_BETA_MAX = 16` (`math_utils/gmrf_field.py:48`). The prior classes are
  parameter-agnostic; refactor the signature to take the prior and field blocks. T needs its
  own scale defaults and positivity (sample log T; chain rule dT = T dlogT).
- Parameter documentation: `params/param_explanation.yml:615` (T), 621-652 (index keys),
  408-418 (hier group text says "sampled indices").

## 4. Where T lives
- Scalar `self.T` (`diffuse_components.py:642`); `sed_param_names = ("beta", "T", "nu_ref")`
  (629); `Component.write_chain` writes `sed/T` (`sky/component.py:196-199`); `read_chain`
  restores only beta (221-228): T is never read back. No `T_map` anywhere in `src/` or `tests/`;
  the prototype `nextgen/prototypes/indexfield/model.py:45-63` has a fixed scalar T_d.

## 5. Tests pinning the layout
- `tests/test_index_field.py`: 347-370 (repeated component and differing nside must raise),
  173-187, 160-171, 206-232 (T-map cases), 109-145 and 652-694 (kernel call shapes), 401-419
  (natural home for a beta+T finite-difference test), 537-583 (covariance-form check with its
  own `dust_sed(nu, b, T)` helper).
- `tests/test_index_field_sampler.py`: 70-92, 97-124, 129-172 (asserts the unpacked keys),
  180-240.
- `tests/test_hier_group_wiring.py`: 72-88, 109-145 (chain round trip), 147-188.
- `tests/test_index_field_mpi.py`: both tests hand-pack a single-spec layout.
- `tests/test_components.py:856-866`: add a `get_dsed_dT` finite-difference test.

## Minimal change set, in order
1. SED: `get_sed(nu, beta=None, T=None)` and `get_dsed_dT`. 2. Component state: `T_map`,
`T_nside`, chain write/read. 3. Specs per (component, parameter) with prefixed keys;
`field_prior_from_params` taking blocks. 4. Layout keyed by (comp, param). 5. `mixing_regions`,
`augment_mixing`, kernel dM with a parameter axis (or two calls). 6. Per-spec gradient folding
and packing. 7. Chain layout nested by parameter. 8. Validators and documentation.
