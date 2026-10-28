"""The index-field sampling group: NUTS over spectral-index fields and amplitude hyperparameters
(M2.4 stage 1, `likelihood: hier`).

One MCMC group with `likelihood: hier` samples the joint target of `RegionIndexLikelihood`
(compsep/index_field.py, plan section 1) with the NUTS sampler of math_utils/nuts.py. The root
rank of the CompSep communicator drives the transitions; the other ranks serve likelihood
evaluations in lock-step, since `RegionIndexLikelihood.logp_and_grad` is collective (every rank
holds part of the bands). After the transitions every rank applies the final state to its
component views (`beta_map` or `beta`, and the hier hyperparameters), so the amplitude Gibbs
draw that follows (`solve_compsep_perpix`, prior- and region-aware) conditions on the new indices.

The sampler state (NUTS adaptation and the current packed vector) persists across Gibbs
iterations per (chain, group) in `IndexSamplerState`, so the warm-up spans the first `n_warmup`
transitions of a run regardless of `numstep` (transitions per Gibbs iteration).
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from mpi4py import MPI
from numpy.typing import NDArray
from pixell.bunch import Bunch

from commander4.compsep.index_field import (AmpPrior, IndexSpec, RegionIndexLikelihood,
                                            StateLayout, component_prior)
from commander4.math_utils import gmrf_field
from commander4.math_utils.nuts import NUTSSampler
from commander4.sky.comp_list import CompList
from commander4.sky.diffuse_components import DiffuseComponent

logger = logging.getLogger(__name__)

NSIDE_BETA_MAX = gmrf_field.NSIDE_BETA_MAX
# Starting point of a field block when a run does not restart from a field chain: z = 0 (the
# field at its mean), mu = the component's current beta, sigma = sigma_scale / 2, u = 0 (rho at
# the middle of its log range). Warm-up moves away from it; recorded in DECISIONS (M2.4).
INIT_SIGMA_FRACTION = 0.5


@dataclass
class NUTSConfig:
    """`nuts:` block of an MCMC group with `likelihood: hier`."""
    n_warmup: int = 500
    target_accept: float = 0.8
    max_depth: int = 10
    parameterization: str = "noncentred"
    init_inv_mass: str = "ones"
    adapt_mass: bool = True
    refresh_inv_mass_at: tuple = ()

    @classmethod
    def from_block(cls, block: Bunch | dict | None) -> "NUTSConfig":
        values = dict(block) if block else {}
        unknown = set(values) - {"n_warmup", "target_accept", "max_depth", "parameterization",
                                 "init_inv_mass", "adapt_mass", "refresh_inv_mass_at"}
        if unknown:
            raise ValueError(f"Unknown nuts settings {sorted(unknown)}.")
        if "refresh_inv_mass_at" in values:
            at = values["refresh_inv_mass_at"]
            at = [at] if isinstance(at, (int, np.integer)) and not isinstance(at, bool) else at
            if not isinstance(at, (list, tuple)) or not all(
                    isinstance(v, (int, np.integer)) and not isinstance(v, bool) and v > 0
                    for v in at):
                raise ValueError("nuts.refresh_inv_mass_at must be a list of positive integers "
                                 "(transition counts during warm-up).")
            values["refresh_inv_mass_at"] = tuple(int(v) for v in at)
        cfg = cls(**values)
        buffer = 50 if cfg.n_warmup >= 100 else cfg.n_warmup // 2
        if any(v > cfg.n_warmup - buffer for v in cfg.refresh_inv_mass_at):
            raise ValueError(f"nuts.refresh_inv_mass_at entries must be below n_warmup by at "
                             f"least {buffer} transitions (the step size re-adapts after a "
                             f"refresh; a refresh at the end of warm-up leaves it at the "
                             f"heuristic value).")
        if cfg.parameterization not in ("noncentred", "centred"):
            raise ValueError(f"nuts.parameterization must be 'noncentred' or 'centred', got "
                             f"{cfg.parameterization!r}.")
        if cfg.init_inv_mass not in ("ones", "fd_hessian"):
            raise ValueError(f"nuts.init_inv_mass must be 'ones' or 'fd_hessian', got "
                             f"{cfg.init_inv_mass!r}.")
        if cfg.n_warmup < 0 or not 0.0 < cfg.target_accept < 1.0 or cfg.max_depth < 1:
            raise ValueError("nuts: n_warmup >= 0, 0 < target_accept < 1, max_depth >= 1.")
        if not isinstance(cfg.adapt_mass, bool):
            raise ValueError(f"nuts.adapt_mass must be true or false, got {cfg.adapt_mass!r}.")
        return cfg


def _read_param(container, key: str, default=None):
    return container[key] if key in container else default


FD_MASS_STEP = 1e-3
FD_MASS_MAX = 10.0


def fd_inverse_mass(logp_and_grad, x: NDArray, h: float = FD_MASS_STEP,
                    cap: float = FD_MASS_MAX) -> NDArray:
    """Initial diagonal inverse mass from a central finite-difference estimate of the diagonal
    of the negative Hessian at x: inv_mass_j = 1 / H_jj where H_jj = -(g_j(x + h e_j) - g_j(x -
    h e_j)) / 2h > 0, else 1 (the prior scale of the non-centred blocks), capped at `cap`.
    This only sets the starting point of the mass adaptation; the target is unchanged. Costs
    2 dim gradient evaluations."""
    x = np.asarray(x, dtype=float)
    diag = np.empty(x.size)
    for j in range(x.size):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        gp = logp_and_grad(xp)[1][j]
        gm = logp_and_grad(xm)[1][j]
        diag[j] = -(gp - gm) / (2.0 * h)
    inv = np.ones(x.size)
    ok = np.isfinite(diag) & (diag > 0.0)
    inv[ok] = np.minimum(1.0 / diag[ok], cap)
    return inv


def field_prior_from_params(comp_params: Bunch, nside_beta: int, centred: bool
                            ) -> gmrf_field.IndexFieldPrior | gmrf_field.IndependentIndexPrior:
    """The spatial prior of a component's index field from its `spectral_index_field` block and
    the required `spectral_index_prior` (mu0, sd0). `type: gmrf` (default): `IndexFieldPrior`,
    basis graph, alpha 2, sigma_scale 0.3, rho prior log-uniform on [h, pi], (mu0, sd0) the
    hyperprior of the field mean. `type: independent`: `IndependentIndexPrior`, every region
    index iid N(mu0, sd0^2), no index hyperparameters (the fixed-nside ladder). `type: iid`:
    `IndexFieldPrior` over an `IidField`, beta = mu 1 + sigma z with the GMRF hyperpriors on mu
    and sigma and no spatial coupling (the exchangeable analogue of the GMRF)."""
    prior = _read_param(comp_params, "spectral_index_prior")
    if prior is None or _read_param(prior, "mean") is None or _read_param(prior, "rms") is None:
        raise ValueError("A spectral-index field needs 'spectral_index_prior' with 'mean' and "
                         "'rms' (the hyperprior of the field mean, or the per-region prior of "
                         "an independent field).")
    settings = dict(_read_param(comp_params, "spectral_index_field", {}) or {})
    kind = str(settings.pop("type", "gmrf"))
    if kind == "independent":
        if settings:
            raise ValueError(f"spectral_index_field type independent takes no other settings; "
                             f"got {sorted(settings)}.")
        return gmrf_field.IndependentIndexPrior(nside_beta, float(prior["mean"]),
                                                float(prior["rms"]))
    if kind == "iid":
        unknown = set(settings) - {"sigma_scale"}
        if unknown:
            raise ValueError(f"spectral_index_field type iid takes only sigma_scale; got "
                             f"{sorted(unknown)}.")
        hyper = gmrf_field.HyperPrior(mu0=float(prior["mean"]), s_mu=float(prior["rms"]),
                                      sigma_scale=float(settings.get("sigma_scale", 0.3)))
        return gmrf_field.IndexFieldPrior(gmrf_field.IidField(gmrf_field.IidBasis(nside_beta)),
                                          hyper, centred=centred)
    if kind != "gmrf":
        raise ValueError(f"spectral_index_field.type must be 'gmrf', 'iid' or 'independent', "
                         f"got {kind!r}.")
    unknown = set(settings) - {"basis", "alpha", "sigma_scale", "rho_prior", "rho_max"}
    if unknown:
        raise ValueError(f"Unknown spectral_index_field settings {sorted(unknown)}.")
    if settings.get("basis", "graph") != "graph":
        raise ValueError("spectral_index_field.basis: only 'graph' is implemented.")
    if nside_beta > NSIDE_BETA_MAX:
        raise ValueError(f"spectral_index_nside {nside_beta} exceeds {NSIDE_BETA_MAX} (dense "
                         "eigendecomposition).")
    basis = gmrf_field.graph_laplacian_basis(nside_beta)
    fld = gmrf_field.MaternField(basis, alpha=float(settings.get("alpha", 2.0)))
    hyper = gmrf_field.HyperPrior(mu0=float(prior["mean"]), s_mu=float(prior["rms"]),
                                  sigma_scale=float(settings.get("sigma_scale", 0.3)),
                                  rho_prior=settings.get("rho_prior", "loguniform"),
                                  rho_max=float(settings.get("rho_max", np.pi)))
    return gmrf_field.IndexFieldPrior(fld, hyper, centred=centred)


def discover_index_specs(comp_list: CompList, selected_comps: list[str] | None, centred: bool
                         ) -> tuple[list[IndexSpec], dict[str, list[DiffuseComponent]]]:
    """The sampled indices of a hier group: one `IndexSpec` per logical component (its I and QU
    views share the index) with `sample_spectral_index: true` (and in `selected_comps` if
    given); `spectral_index_nside` > 0 gives a field, 0 or absent a single index. Returns the
    specs and the views of each sampled component by name."""
    selected = None if selected_comps is None else set(selected_comps)
    views: dict[str, list[DiffuseComponent]] = {}
    specs: list[IndexSpec] = []
    for comp in comp_list:
        if not hasattr(comp, "beta"):
            continue
        if selected is not None and comp.comp_name not in selected:
            continue
        if not bool(_read_param(comp.comp_params, "sample_spectral_index", False)):
            continue
        if comp.comp_name in views:
            views[comp.comp_name].append(comp)
            continue
        views[comp.comp_name] = [comp]
        nside_beta = int(_read_param(comp.comp_params, "spectral_index_nside", 0) or 0)
        prior = _read_param(comp.comp_params, "spectral_index_prior")
        if prior is None:
            raise ValueError(f"Component {comp.comp_name!r}: 'likelihood: hier' needs a "
                             "'spectral_index_prior' (gaussian, mean and rms).")
        mean, rms = float(prior["mean"]), float(prior["rms"])
        if nside_beta > 0:
            fld = field_prior_from_params(comp.comp_params, nside_beta, centred)
        else:
            fld = None
        specs.append(IndexSpec(comp.comp_name, fld, mean, rms))
    return specs, views


def hier_component_names(comp_list: CompList, amp_comps: list[str] | None) -> list[str]:
    """Names of the amplitude group's components carrying a hier `amp_prior`, in list order."""
    names = []
    for comp in comp_list:
        if amp_comps is not None and comp.comp_name not in amp_comps:
            continue
        if component_prior(comp).is_hier and comp.comp_name not in names:
            names.append(comp.comp_name)
    return names


def amp_prior_from_params(comp: DiffuseComponent) -> AmpPrior | None:
    """The component's `amp_prior` block as an `AmpPrior`, or None when absent.

    `{type: gaussian, mean, rms}` (scalars, uK_RJ at the component's nu_ref) or `{type: hier,
    template: <fits>, template_floor, mu_prior: [mean, sd], log_s_prior: [mean, sd]}`. The hier
    template is a HEALPix map (field 0) in uK_RJ at the component's nu_ref; it is floored at
    `template_floor` (default: 1e-3 of its mean) so that it is strictly positive, and ud_graded
    to the data nside when the prior is used.
    """
    block = _read_param(comp.comp_params, "amp_prior")
    if not block:
        return None
    kind = _read_param(block, "type")
    if kind == "gaussian":
        mean, rms = _read_param(block, "mean"), _read_param(block, "rms")
        if mean is None or rms is None:
            raise ValueError(f"Component {comp.comp_name!r}: gaussian amp_prior needs 'mean' and "
                             "'rms'.")
        unknown = set(block) - {"type", "mean", "rms"}
        if unknown:
            raise ValueError(f"Component {comp.comp_name!r}: unknown amp_prior keys "
                             f"{sorted(unknown)}.")
        return AmpPrior("gaussian", mean=float(mean), rms=float(rms))
    if kind == "hier":
        import healpy as hp
        unknown = set(block) - {"type", "template", "template_floor", "mu_prior", "log_s_prior"}
        if unknown:
            raise ValueError(f"Component {comp.comp_name!r}: unknown amp_prior keys "
                             f"{sorted(unknown)}.")
        path = _read_param(block, "template")
        if not path:
            raise ValueError(f"Component {comp.comp_name!r}: hier amp_prior needs a 'template' "
                             "map.")
        if comp.eval_pol != "I":
            # Intensity-only prior (Q and U are signed, the template is positive): the QU view
            # keeps the flat prior; the component's hyperparameters are those of its I view.
            logger.info(f"Component {comp.comp_name!r} (QU): hier amp_prior applies to intensity "
                        "only; the QU amplitudes keep a flat prior.")
            return None
        template = np.asarray(hp.read_map(path, field=0), dtype=np.float64)
        floor = _read_param(block, "template_floor")
        floor = 1e-3 * float(np.mean(np.abs(template))) if floor is None else float(floor)
        if floor <= 0.0:
            raise ValueError(f"Component {comp.comp_name!r}: template_floor must be positive.")
        template = np.maximum(template, floor)
        mu_prior = tuple(float(v) for v in _read_param(block, "mu_prior", (1.0, 0.5)))
        log_s_prior = tuple(float(v) for v in _read_param(block, "log_s_prior",
                                                          (float(np.log(0.3)), 1.0)))
        if len(mu_prior) != 2 or len(log_s_prior) != 2 or mu_prior[1] <= 0 or log_s_prior[1] <= 0:
            raise ValueError(f"Component {comp.comp_name!r}: mu_prior and log_s_prior are "
                             "[mean, sd] with sd > 0.")
        return AmpPrior("hier", template=template, mu=mu_prior[0], log_s=log_s_prior[0],
                        mu_prior=mu_prior, log_s_prior=log_s_prior)
    raise ValueError(f"Component {comp.comp_name!r}: amp_prior type must be 'gaussian' or "
                     f"'hier', got {kind!r}.")


def load_amp_priors(comp_list: CompList, amplitude_method: str | None) -> None:
    """Attach `amp_prior` objects to every component view that configures one. Only the per-pixel
    amplitude solver applies them; with CG amplitudes they are an error rather than silently
    ignored."""
    for comp in comp_list:
        if not isinstance(comp, DiffuseComponent):
            continue
        prior = amp_prior_from_params(comp)
        if prior is None:
            continue
        if amplitude_method != "per_pixel":
            raise ValueError(f"Component {comp.comp_name!r} sets amp_prior, which only the "
                             "per-pixel amplitude solver applies (per_pixel_sampling_groups).")
        comp.amp_prior = prior
        logger.info(f"Component {comp.comp_name!r} ({comp.eval_pol}): amp_prior {prior.type}.")


@dataclass
class IndexSamplerState:
    """NUTS sampler and current packed vector of one (chain, group), kept across Gibbs
    iterations."""
    sampler: NUTSSampler
    x: NDArray
    n_transitions: int = 0
    mass_initialized: bool = False
    refreshed: set = field(default_factory=set)   # refresh_inv_mass_at entries already applied


_STATES: dict[tuple[int, str], IndexSamplerState] = {}


def reset_sampler_states() -> None:
    _STATES.clear()


def initial_state(layout: StateLayout, views: dict[str, list[DiffuseComponent]],
                  hier_priors: dict[str, AmpPrior]) -> NDArray:
    """The packed vector from the components' current state: a field block from `beta_map` when
    the component has one at the layout's nside (a restart), else the constant start described
    at the top; a scalar block from `beta`; (mu, log s) from the hier priors."""
    blocks = {}
    for spec in layout.specs:
        comp = views[spec.comp_name][0]
        if spec.field is None:
            blocks[spec.comp_name] = np.array([float(comp.beta)])
            continue
        fld = spec.field
        beta_map = getattr(comp, "beta_map", None)
        if isinstance(fld, gmrf_field.IndependentIndexPrior):
            # Independent regions: the state is beta itself; restart from the stored map, else
            # the constant start at the component's beta.
            if beta_map is not None and beta_map.shape == (fld.n_beta,):
                blocks[spec.comp_name] = np.asarray(beta_map, dtype=float).copy()
            else:
                blocks[spec.comp_name] = np.full(fld.n_beta, float(comp.beta))
            continue
        hyper = fld.hyper
        if beta_map is not None and beta_map.shape == (fld.n_beta,):
            mu = float(beta_map.mean())
            resid = beta_map - mu
            # Centred coordinates y = V^T (beta - mu); sigma from the field's realized rms at
            # unit variance, rho at the middle of its range.
            u = 0.0
            rho = float(fld.transform.rho(u))
            y = fld.basis.V.T @ resid
            var_unit = float(np.mean(fld.field.pixel_variances(rho)))
            sigma = max(float(np.sqrt(max(np.mean(resid**2), 1e-12) / var_unit)), 1e-3)
            x_c = fld.pack(y, mu, np.log(sigma), u)
            blocks[spec.comp_name] = x_c if fld.centred else fld.to_noncentred(x_c)
        else:
            sigma = INIT_SIGMA_FRACTION * hyper.sigma_scale
            blocks[spec.comp_name] = fld.pack(np.zeros(fld.K), float(comp.beta), np.log(sigma),
                                              0.0)
    hypers = {name: (hier_priors[name].mu, hier_priors[name].log_s)
              for name in layout.hier_names}
    return layout.pack(blocks, hypers)


def apply_state(x: NDArray, layout: StateLayout, views: dict[str, list[DiffuseComponent]],
                hier_views: dict[str, list[DiffuseComponent]]) -> None:
    """Write a packed state into the component views on this rank: `beta_map` (field) or
    `beta` (scalar), and the hier (mu, log s)."""
    betas = layout.betas(x)
    for spec in layout.specs:
        for comp in views.get(spec.comp_name, []):
            if spec.field is None:
                comp.beta = float(x[layout.index_slices[spec.comp_name]][0])
            else:
                comp.set_beta_map(betas[spec.comp_name], layout.nside_beta)
    for name, (mu, log_s) in layout.hypers(x).items():
        for comp in hier_views.get(name, []):
            comp.amp_prior.mu, comp.amp_prior.log_s = mu, log_s


_TOKEN_EVAL, _TOKEN_STOP = 1, 0


class IndexFieldSamplingGroup:
    """Drives `numstep` NUTS transitions of the joint index-field target on the root rank, with
    every CompSep rank serving the collective likelihood, then applies the state everywhere.

    Args:
        compsep_comm: the CompSep communicator (all band views).
        likelihood: the `RegionIndexLikelihood` built for this group (collective object).
        layout: its `StateLayout`.
        views: sampled component views on this rank by name (`discover_index_specs`).
        hier_views: hier component views on this rank by name.
        hier_priors: one `AmpPrior` per hier component name (hyperprior settings), on every rank.
        nuts: the group's `NUTSConfig`.
        root: rank that runs the sampler.
    """

    def __init__(self, compsep_comm: MPI.Comm, likelihood: RegionIndexLikelihood,
                 layout: StateLayout, views: dict[str, list[DiffuseComponent]],
                 hier_views: dict[str, list[DiffuseComponent]],
                 hier_priors: dict[str, AmpPrior], nuts: NUTSConfig, root: int = 0):
        self.comm = compsep_comm
        self.likelihood = likelihood
        self.layout = layout
        self.views = views
        self.hier_views = hier_views
        self.hier_priors = hier_priors
        self.nuts = nuts
        self.root = root
        self.is_root = compsep_comm.Get_rank() == root
        likelihood.set_hier_defaults(hier_priors)

    def state(self, chain: int, group_name: str, rng: np.random.Generator) -> IndexSamplerState:
        key = (chain, group_name)
        if key not in _STATES:
            x0 = initial_state(self.layout, self.views, self.hier_priors) if self.views else None
            x0 = self.comm.bcast(x0, root=self.root)
            sampler = NUTSSampler(self.layout.dim, rng, n_warmup=self.nuts.n_warmup,
                                  target_accept=self.nuts.target_accept,
                                  max_depth=self.nuts.max_depth,
                                  adapt_mass=self.nuts.adapt_mass)
            _STATES[key] = IndexSamplerState(sampler, x0)
        return _STATES[key]

    def _serve(self) -> None:
        """Non-root ranks: evaluate the collective likelihood until the root sends the stop
        token."""
        while True:
            token, x = self.comm.bcast(None, root=self.root)
            if token == _TOKEN_STOP:
                return
            self.likelihood.logp_and_grad(x)

    def _logp_and_grad_root(self, x: NDArray) -> tuple[float, NDArray]:
        self.comm.bcast((_TOKEN_EVAL, x), root=self.root)
        return self.likelihood.logp_and_grad(x)

    def run(self, *, numstep: int, chain: int, group_name: str, rng: np.random.Generator) -> dict:
        """Take `numstep` transitions and apply the final state. Returns the diagnostics of the
        iteration (root) or {} (other ranks)."""
        st = self.state(chain, group_name, rng)
        st.sampler.rng = rng
        stats: dict = {}
        if self.is_root:
            x = st.x
            infos = []
            try:
                if not st.mass_initialized and self.nuts.init_inv_mass == "fd_hessian":
                    # Serving ranks are in _serve() by now; every evaluation below is collective.
                    inv = fd_inverse_mass(self._logp_and_grad_root, x)
                    st.sampler.inv_mass = inv
                    logger.info(f"Index-field group {group_name!r}: initial inverse mass from "
                                f"the finite-difference Hessian, {2 * x.size} gradients, range "
                                f"{inv.min():.3g} to {inv.max():.3g}, "
                                f"{int(np.sum(inv == 1.0))} coordinates left at 1")
                st.mass_initialized = True
                due = [v for v in self.nuts.refresh_inv_mass_at
                       if v <= st.n_transitions and v not in st.refreshed]
                if due and st.sampler.iteration < st.sampler.n_warmup:
                    # Warm-up only, once per configured count, at the first Gibbs iteration whose
                    # transition count has reached it (numstep need not divide the count): a
                    # fresh finite-difference Hessian at the current state replaces the metric
                    # (the start-state Hessian is stale once the chain has left the transient)
                    # and the step size is re-found; same target.
                    st.refreshed.update(due)
                    inv = fd_inverse_mass(self._logp_and_grad_root, x)
                    st.sampler.reset_metric(inv, x, self._logp_and_grad_root)
                    logger.info(f"Index-field group {group_name!r}: inverse mass refreshed from "
                                f"the finite-difference Hessian at transition "
                                f"{st.n_transitions} (configured {sorted(due)}), "
                                f"{2 * x.size} gradients, range {inv.min():.3g} to "
                                f"{inv.max():.3g}, step size {st.sampler.step_size:.3g}")
                for _ in range(numstep):
                    x, info = st.sampler.transition(x, self._logp_and_grad_root)
                    infos.append(info)
            finally:
                self.comm.bcast((_TOKEN_STOP, None), root=self.root)
            st.x = x
            st.n_transitions += numstep
            stats = {
                "step_size": float(st.sampler.step_size),
                "accept_stat": float(np.mean([i.accept_stat for i in infos])),
                "n_leapfrog": int(sum(i.n_leapfrog for i in infos)),
                "tree_depth": int(max(i.tree_depth for i in infos)),
                "divergent": int(sum(i.divergent for i in infos)),
                "adapting": bool(infos[-1].adapting),
                "n_transitions": int(st.n_transitions),
                "n_rejected": int(self.likelihood.n_rejected),
                "x": st.x.copy(),
            }
            logger.info(f"Index-field group {group_name!r}: {numstep} NUTS transition(s), step "
                        f"size {stats['step_size']:.3g}, mean accept {stats['accept_stat']:.2f}, "
                        f"{stats['n_leapfrog']} gradients, {stats['divergent']} divergent"
                        + (" (warm-up)" if stats["adapting"] else ""))
        else:
            self._serve()
        x = self.comm.bcast(st.x if self.is_root else None, root=self.root)
        st.x = x
        apply_state(x, self.layout, self.views, self.hier_views)
        return stats

    def unpacked(self, x: NDArray) -> dict:
        """Per-component arrays of a packed state for the chain file: `index_field/{beta_map,
        beta_nside, mu, log_sigma, u, z}` per field spec, `beta` per scalar spec, and
        `amp_prior/{mu, log_s}` per hier component."""
        out: dict = {}
        betas = self.layout.betas(x)
        for spec in self.layout.specs:
            block = x[self.layout.index_slices[spec.comp_name]]
            if spec.field is None:
                out[spec.comp_name] = {"beta": float(block[0])}
            elif isinstance(spec.field, gmrf_field.IndependentIndexPrior):
                out[spec.comp_name] = {"index_field": {
                    "beta_map": betas[spec.comp_name], "beta_nside": self.layout.nside_beta,
                    "parameterization": "independent"}}
            else:
                fld = spec.field
                x_nc = fld.to_noncentred(block) if fld.centred else block
                z, mu, ls, u = fld.split(x_nc)
                out[spec.comp_name] = {"index_field": {
                    "beta_map": betas[spec.comp_name], "beta_nside": self.layout.nside_beta,
                    "mu": mu, "log_sigma": ls, "u": u, "z": z.copy(),
                    "parameterization": "centred" if fld.centred else "noncentred",
                    "field_type": "iid" if isinstance(fld.field, gmrf_field.IidField)
                    else "gmrf"}}
        for name, (mu, log_s) in self.layout.hypers(x).items():
            out.setdefault(name, {})["amp_prior"] = {"mu": mu, "log_s": log_s}
        return out
