"""Spatially varying spectral indices with a hierarchical amplitude prior (M2.4 stage 1).

This module holds the pieces that make the per-pixel amplitude model prior- and region-aware:

- `AmpPrior`: the per-pixel amplitude prior of one component, `none` (flat, today's model),
  `gaussian` (a_p ~ N(m_p, t_p^2), fixed) or `hier` (a_p ~ N(mu g_p, (s g_p)^2) with a fixed
  strictly positive template g and sampled hyperparameters mu and log s).
- Pseudo-band augmentation: every component's prior enters the per-pixel kernels as one extra
  band with map m_p, inverse variance 1/t_p^2 (0 for a flat component) and mixing row e_c. The
  kernels' (A, b) then are exactly (A', b') = (A + T^-1, b + T^-1 m), so
  `marginal_loglike_regions` returns the Gaussian-prior collapsed likelihood
  1/2 b'^T A'^-1 b' - 1/2 log det A' and `solve_perpix_regions` draws N(A'^-1 b', A'^-1), the
  amplitude conditional under the prior (its fluctuation on a pseudo-band is T^-1/2 eta).
- Region maps and per-region mixing matrices from the components' `beta_map` (one index value
  per pixel at `beta_nside`, each data pixel belonging to its parent index pixel).

The target distribution is stated in nextgen/plans/M2.4_commander4_index_field.md, section 1.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace

import healpy as hp
import numpy as np
from mpi4py import MPI
from numpy.typing import NDArray

from commander4.backend import compsep as cpp_compsep
from commander4.data_models.detector_map import DetectorMap
from commander4.math_utils.gmrf_field import IndexFieldPrior
from commander4.sky.diffuse_components import DiffuseComponent

logger = logging.getLogger(__name__)

# Hyperprior defaults of the hier prior (E016): mu ~ N(1, 0.5^2), log s ~ N(log 0.3, 1^2), and
# the guard below which log s has zero density (E017: the prior precision 1/(s g)^2 swamps
# M^T N^-1 M in double precision; hyperprior mass below the guard is 5e-12).
HIER_MU_PRIOR = (1.0, 0.5)
HIER_LOG_S_PRIOR = (float(np.log(0.3)), 1.0)
HIER_LOGS_MIN = -8.0
AMP_PRIOR_TYPES = ("none", "gaussian", "hier")


@dataclass
class AmpPrior:
    """The per-pixel amplitude prior of one component, in the component's amplitude unit (uK_RJ
    at its reference frequency) on the data pixelization.

    Attributes:
        type: "none", "gaussian" or "hier".
        mean: gaussian: prior mean m_p, scalar or (npix,).
        rms: gaussian: prior sd t_p > 0, scalar or (npix,).
        template: hier: strictly positive shape template g_p (npix,).
        mu, log_s: hier: current hyperparameters (sampled by the index-field group).
        mu_prior, log_s_prior: hier: (mean, sd) of the Gaussian hyperpriors on mu and log s.
    """
    type: str = "none"
    mean: float | NDArray | None = None
    rms: float | NDArray | None = None
    template: NDArray | None = None
    mu: float = 1.0
    log_s: float = HIER_LOG_S_PRIOR[0]
    mu_prior: tuple[float, float] = HIER_MU_PRIOR
    log_s_prior: tuple[float, float] = HIER_LOG_S_PRIOR
    # Maps regridded to the stream's pixel count, by (name, npix); the configured maps may be at
    # another nside than the data (ud_grade keeps a positive template positive).
    _regrid: dict = field(default_factory=dict, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.type not in AMP_PRIOR_TYPES:
            raise ValueError(f"amp_prior type {self.type!r} is not one of {AMP_PRIOR_TYPES}.")
        if self.type == "gaussian":
            if self.mean is None or self.rms is None:
                raise ValueError("A gaussian amp_prior needs 'mean' and 'rms'.")
            if np.any(np.asarray(self.rms, dtype=np.float64) <= 0.0):
                raise ValueError("A gaussian amp_prior needs rms > 0.")
        if self.type == "hier":
            if self.template is None:
                raise ValueError("A hier amp_prior needs a strictly positive 'template'.")
            self.template = np.asarray(self.template, dtype=np.float64)
            if self.template.ndim != 1 or np.any(self.template <= 0.0):
                raise ValueError("The hier amp_prior template must be a strictly positive 1-D "
                                 "map at the data nside.")

    def _map_at(self, name: str, npix: int) -> NDArray:
        """The map attribute `name` (scalar or 1-D map) as an (npix,) array, ud_graded when its
        pixel count differs from npix."""
        arr = np.asarray(getattr(self, name), dtype=np.float64)
        if arr.ndim == 0:
            return np.broadcast_to(arr, (npix,))
        if arr.shape == (npix,):
            return arr
        key = (name, npix)
        if key not in self._regrid:
            try:
                nside_out = hp.npix2nside(npix)
                hp.npix2nside(arr.size)
            except (ValueError, TypeError) as error:
                raise ValueError(f"amp_prior {name} has {arr.size} pixels, the stream {npix}; "
                                 "both must be HEALPix maps.") from error
            self._regrid[key] = hp.ud_grade(arr, nside_out)
        return self._regrid[key]

    @property
    def is_flat(self) -> bool:
        return self.type == "none"

    @property
    def is_hier(self) -> bool:
        return self.type == "hier"

    def pseudo_band(self, npix: int) -> tuple[NDArray, NDArray]:
        """(m, inv_T) of this prior as a pseudo-band: the prior mean map and inverse variance map
        (npix,), both zero for a flat prior."""
        if self.type == "none":
            return np.zeros(npix), np.zeros(npix)
        if self.type == "gaussian":
            m = self._map_at("mean", npix).copy()
            t = self._map_at("rms", npix)
            return m, 1.0 / t**2
        g = self._map_at("template", npix)
        return self.mu * g, 1.0 / (np.exp(self.log_s) * g)**2

    def trailing_terms(self, used: NDArray) -> float:
        """The prior's beta-independent terms -1/2 sum_p m^T T^-1 m - 1/2 sum_p log det T over
        the pixels `used` (bool, npix). Zero for `none`; for `gaussian` a constant that the
        samplers may drop, returned for completeness; for `hier` it depends on (mu, log s) and
        must be kept, summed over exactly the pixel set the kernel used."""
        if self.type == "none":
            return 0.0
        m, inv_T = self.pseudo_band(used.size)
        return float(-0.5 * np.sum(m[used]**2 * inv_T[used]) + 0.5 * np.sum(np.log(inv_T[used])))

    def hyper_grad(self, x: NDArray, diag_Ainv: NDArray, used: NDArray) -> tuple[float, float]:
        """(dF/dmu, dF/dlog s) of the hier target for this component, from the kernel extras of
        the augmented stream: x the conditional means A'^-1 b' and diag_Ainv the diagonal of A'^-1
        of this component's column (npix,), on the used pixels. With P = T^-1 (E016):
            dF/dmu    = sum_p P g (x - m),
            dF/dlog s = sum_p [P (x - m)^2 + P (A'^-1)_cc] - n_used.
        """
        if self.type != "hier":
            raise ValueError("hyper_grad is defined for the hier prior only.")
        m, P = self.pseudo_band(used.size)
        dx = (x - m)[used]
        P, g = P[used], self._map_at("template", used.size)[used]
        return (float(np.sum(P * g * dx)),
                float(np.sum(P * dx**2 + P * diag_Ainv[used]) - np.count_nonzero(used)))

    def log_hyperprior(self) -> float:
        """log p(mu) + log p(log s) of the hier hyperpriors (-inf below the log s guard), 0 for
        the other types. Normalizing constants dropped."""
        if self.type != "hier":
            return 0.0
        if self.log_s < HIER_LOGS_MIN:
            return -np.inf
        return -0.5 * ((self.mu - self.mu_prior[0]) / self.mu_prior[1])**2 \
            - 0.5 * ((self.log_s - self.log_s_prior[0]) / self.log_s_prior[1])**2

    def log_hyperprior_grad(self) -> tuple[float, float]:
        """(d/dmu, d/dlog s) of `log_hyperprior`."""
        if self.type != "hier":
            return 0.0, 0.0
        return (-(self.mu - self.mu_prior[0]) / self.mu_prior[1]**2,
                -(self.log_s - self.log_s_prior[0]) / self.log_s_prior[1]**2)


def component_prior(comp: DiffuseComponent) -> AmpPrior:
    """The component's `amp_prior`, or the flat prior when it has none."""
    prior = getattr(comp, "amp_prior", None)
    return prior if prior is not None else AmpPrior("none")


def augment_stream(maps_sky: NDArray, maps_inv_var: NDArray, comps: list[DiffuseComponent]
                   ) -> tuple[NDArray, NDArray]:
    """Append one pseudo-band per component to a (nband, npix) stream: the prior mean map and
    inverse variance map of `component_prior(comp)` (zeros for a flat prior). The returned arrays
    have shape (nband + ncomp, npix); the pseudo-band of component c is row nband + c."""
    npix = maps_sky.shape[1]
    rows_m, rows_w = [], []
    for comp in comps:
        m, inv_T = component_prior(comp).pseudo_band(npix)
        rows_m.append(m)
        rows_w.append(inv_T)
    return (np.concatenate([maps_sky, np.array(rows_m)], axis=0),
            np.concatenate([maps_inv_var, np.array(rows_w)], axis=0))


def augment_mixing(M_reg: NDArray, dM_reg: NDArray | None = None
                   ) -> NDArray | tuple[NDArray, NDArray]:
    """Append the identity rows of the pseudo-bands to per-region mixing matrices
    (nreg, nband, ncomp) -> (nreg, nband + ncomp, ncomp); the derivative rows are zero."""
    nreg, _, ncomp = M_reg.shape
    rows = np.broadcast_to(np.eye(ncomp), (nreg, ncomp, ncomp))
    M_aug = np.ascontiguousarray(np.concatenate([M_reg, rows], axis=1))
    if dM_reg is None:
        return M_aug
    dM_aug = np.ascontiguousarray(np.concatenate([dM_reg, np.zeros_like(rows)], axis=1))
    return M_aug, dM_aug


def needs_region_solver(comps: list[DiffuseComponent]) -> bool:
    """Whether the amplitude draw must go through the region kernel: some component has a
    `beta_map` or a non-flat amplitude prior. Otherwise the plain `solve_perpix` path is used, so
    runs without fields or priors reproduce today's draws for the same random state."""
    return any(getattr(c, "beta_map", None) is not None or not component_prior(c).is_flat
               for c in comps)


def field_nside(comps: list[DiffuseComponent]) -> int:
    """The common `beta_nside` of the components carrying a `beta_map` (0 if none). Stage 1
    requires one index resolution per amplitude group."""
    nsides = {int(c.beta_nside) for c in comps if getattr(c, "beta_map", None) is not None}
    if len(nsides) > 1:
        raise ValueError(f"The components of one amplitude group must share one "
                         f"spectral_index_nside; got {sorted(nsides)}.")
    return nsides.pop() if nsides else 0


def region_map(nside_beta: int, nside_data: int, keep: NDArray | None = None) -> NDArray:
    """(npix_data,) int64 region index of every data pixel: the ring index of its parent pixel at
    `nside_beta` (0 everywhere for nside_beta 0), -1 where `keep` (bool, npix_data) is False.
    Requires nside_data >= nside_beta > 0, or nside_beta 0."""
    npix = 12 * nside_data * nside_data
    if nside_beta == 0:
        region = np.zeros(npix, dtype=np.int64)
    else:
        if nside_data < nside_beta:
            raise ValueError(f"spectral_index_nside {nside_beta} exceeds the data nside "
                             f"{nside_data}.")
        labels = np.arange(12 * nside_beta * nside_beta, dtype=np.float64)
        region = np.rint(hp.ud_grade(labels, nside_data)).astype(np.int64)
    if keep is not None:
        region = np.where(keep, region, -1)
    return region


def mixing_regions(freqs: NDArray, comps: list[DiffuseComponent], nreg: int,
                   with_derivative: bool = False) -> NDArray | tuple[NDArray, NDArray]:
    """Per-region mixing matrices M_reg (nreg, nband, ncomp) from the components' current state:
    column c is `comp.get_sed(freqs, beta=comp.beta_map)` for a component with a field (its
    `beta_map` indexed by region) and `comp.get_sed(freqs)` broadcast over regions otherwise.
    With `with_derivative`, also dM_reg (zero columns for components without a field)."""
    nband = len(freqs)
    M_reg = np.empty((nreg, nband, len(comps)))
    dM_reg = np.zeros_like(M_reg)
    for c, comp in enumerate(comps):
        beta_map = getattr(comp, "beta_map", None)
        if beta_map is None:
            M_reg[:, :, c] = comp.get_sed(freqs)[None, :]
        else:
            if beta_map.shape != (nreg,):
                raise ValueError(f"Component {comp.comp_name!r}: beta_map has {beta_map.size} "
                                 f"pixels, the amplitude group has {nreg} regions.")
            M_reg[:, :, c] = comp.get_sed(freqs, beta=beta_map)
            if with_derivative:
                dM_reg[:, :, c] = comp.get_dsed_dbeta(freqs, beta=beta_map)
    if with_derivative:
        return M_reg, dM_reg
    return M_reg


# ----------------------------------------------------------------------------------------------
# The joint target over index fields and amplitude hyperparameters (plan section 1)
# ----------------------------------------------------------------------------------------------

@dataclass
class IndexSpec:
    """One sampled spectral index: the component (by `comp_name`, shared by its I and QU views),
    its field prior (`IndexFieldPrior` at spectral_index_nside > 0) or None for a single index,
    and the Gaussian index prior (mu0, sd0) of `spectral_index_prior`, which for a field is the
    hyperprior on the field mean (already inside `field.hyper`) and for a scalar is the prior on
    the index itself."""
    comp_name: str
    field: IndexFieldPrior | None
    prior_mean: float
    prior_sd: float

    def __post_init__(self) -> None:
        # For a field the (mu0, sd0) live in field.hyper; the spec's copy must agree, so that one
        # source of configuration cannot drift from the other.
        if self.field is not None and (self.prior_mean != self.field.hyper.mu0
                                       or self.prior_sd != self.field.hyper.s_mu):
            raise ValueError(f"IndexSpec {self.comp_name!r}: prior ({self.prior_mean}, "
                             f"{self.prior_sd}) differs from the field hyperprior "
                             f"({self.field.hyper.mu0}, {self.field.hyper.s_mu}).")

    @property
    def dim(self) -> int:
        return self.field.dim if self.field is not None else 1

    @property
    def n_beta(self) -> int:
        return self.field.n_beta if self.field is not None else 1

    @property
    def nside_beta(self) -> int:
        return self.field.basis.nside_beta if self.field is not None else 0


class StateLayout:
    """Packing of the sampled vector x = (index block of every IndexSpec in order, then (mu, log
    s) of every hier component in `hier_names` order). The index block of a field spec is the
    IndexFieldPrior state (field block, mu, log sigma, u); of a scalar spec the index itself."""

    def __init__(self, specs: list[IndexSpec], hier_names: list[str]):
        if len({s.comp_name for s in specs}) != len(specs):
            raise ValueError("Each component may carry one sampled index.")
        self.specs = list(specs)
        self.hier_names = list(hier_names)
        self.index_slices: dict[str, slice] = {}
        self.hier_slices: dict[str, slice] = {}
        pos = 0
        for spec in self.specs:
            self.index_slices[spec.comp_name] = slice(pos, pos + spec.dim)
            pos += spec.dim
        for name in self.hier_names:
            self.hier_slices[name] = slice(pos, pos + 2)
            pos += 2
        self.dim = pos
        nsides = {s.nside_beta for s in self.specs if s.nside_beta > 0}
        if len(nsides) > 1:
            raise ValueError(f"All index fields of one group must share spectral_index_nside; "
                             f"got {sorted(nsides)}.")
        self.nside_beta = nsides.pop() if nsides else 0
        self.nreg = 12 * self.nside_beta**2 if self.nside_beta > 0 else 1

    def spec(self, comp_name: str) -> IndexSpec:
        return next(s for s in self.specs if s.comp_name == comp_name)

    def betas(self, x: NDArray) -> dict[str, NDArray]:
        """The index per region (nreg,) of every sampled component at state x (a scalar index is
        broadcast over the regions)."""
        out = {}
        for spec in self.specs:
            block = x[self.index_slices[spec.comp_name]]
            if spec.field is None:
                out[spec.comp_name] = np.full(self.nreg, float(block[0]))
            else:
                out[spec.comp_name] = spec.field.beta(block)
        return out

    def hypers(self, x: NDArray) -> dict[str, tuple[float, float]]:
        """(mu, log s) of every hier component at state x."""
        return {name: (float(x[sl][0]), float(x[sl][1])) for name, sl in self.hier_slices.items()}

    def pack(self, index_blocks: dict[str, NDArray], hypers: dict[str, tuple[float, float]]
             ) -> NDArray:
        x = np.empty(self.dim)
        for spec in self.specs:
            x[self.index_slices[spec.comp_name]] = np.atleast_1d(index_blocks[spec.comp_name])
        for name in self.hier_names:
            x[self.hier_slices[name]] = hypers[name]
        return x


class RegionIndexLikelihood:
    """Log density and gradient of the joint target over index fields and amplitude
    hyperparameters, with the per-pixel amplitudes integrated out (plan section 1):

        log p(x | d) = sum_streams sum_p [1/2 b'^T A'^-1 b' - 1/2 log det A']
                       + sum_{hier c} [-1/2 sum_p m^T T^-1 m - 1/2 sum_p log det T + log p(mu_c)
                                       + log p(log s_c)]
                       + sum_{field specs} [log N(z; 0, I) or the centred equivalent
                                            + hyperpriors and Jacobians]
                       + sum_{scalar specs} log N(beta; mu0, sd0^2),

    where the pixel sums run over the pixels the kernel used (region >= 0 under the chisq mask,
    at least ncomp observed bands counting the proper pseudo-bands) and A' = A + T^-1,
    b' = b + T^-1 m with the amplitude priors as pseudo-bands (module docstring). That pixel set
    is state independent; a state at which A' fails to be positive definite in some used pixel
    has zero density (below), it is never dropped from the sum. With proper priors, pixels
    observed in fewer than ncomp data bands enter the sum (their proper marginal exists), unlike
    in `PerPixelMarginalLikelihood`, whose flat-prior marginal excludes them. The log det term
    is always kept. The hier prior is intensity-only: a QU view of a hier component carries the
    flat prior (index_field_sampler.load_amp_priors), the hyperparameters are shared through the
    component name and their terms come from the intensity stream alone. The
    gradient with respect to the index of region r and component c is the kernel's analytic
    gradient; the chain rule to the field coordinates and the field prior are those of
    `IndexFieldPrior.logp_and_grad`; the (mu, log s) gradients are `AmpPrior.hyper_grad`.

    Construction is collective over ``compsep_comm`` and mirrors `PerPixelMarginalLikelihood`: the
    residual band maps are gathered once to the root of ``solver_comm``, which then holds this
    rank's polarization streams. `logp_and_grad` is also collective: every rank evaluates its
    streams, the likelihood value and its gradients are summed over ``compsep_comm``, and every
    rank then adds the prior terms, so all ranks return the same (logp, grad). A state at which
    any kernel skipped a pixel (A' not positive definite) or produced a non-finite value, or with
    log s below `HIER_LOGS_MIN`, has zero density: (-inf, zeros), the NUTS contract.

    Args:
        compsep_comm, solver_comm, residual_data, active_comps, keep_mask, nthreads: as for
            `PerPixelMarginalLikelihood`.
        layout: the `StateLayout` of the sampled vector; sampled components are matched to the
            active components by `comp_name` (an I and a QU view share one index).
    """

    def __init__(self, compsep_comm: MPI.Comm, solver_comm: MPI.Comm | None,
                 residual_data: DetectorMap | None, active_comps: list[DiffuseComponent],
                 layout: StateLayout, keep_mask: NDArray | None, nthreads: int):
        self.compsep_comm = compsep_comm
        self.layout = layout
        self.active_comps = list(active_comps)
        self.nthreads = nthreads
        self.is_solver_root = solver_comm is not None and solver_comm.Get_rank() == 0
        self.streams: list[tuple[NDArray, NDArray, NDArray]] = []
        self.region = None
        self.n_rejected = 0
        self.nskipped_last = 0
        self._hier_defaults: dict[str, AmpPrior] = {}
        for comp in self.active_comps:
            if not isinstance(comp, DiffuseComponent):
                raise ValueError(f"The index-field likelihood does not support "
                                 f"{type(comp).__name__} components.")
        names = [c.comp_name for c in self.active_comps]
        for spec in layout.specs:
            if spec.comp_name in names and not hasattr(
                    self.active_comps[names.index(spec.comp_name)], "beta"):
                raise ValueError(f"Component {spec.comp_name!r} has no spectral index to sample.")
        for name in layout.hier_names:
            if name in names:
                view = self.active_comps[names.index(name)]
                if view.eval_pol == "I" and not component_prior(view).is_hier:
                    raise ValueError(f"Component {name!r} is in the hier block but its intensity "
                                     "view has no hier amp_prior.")
        if solver_comm is None:
            return
        inv_var = np.array(residual_data.map_inv_var, dtype=np.float64, copy=True)
        all_freq = solver_comm.gather(residual_data.nu, root=0)
        all_map_sky = solver_comm.gather(residual_data.map_sky, root=0)
        all_inv_var = solver_comm.gather(inv_var, root=0)
        if not self.is_solver_root:
            return
        nside = residual_data.nside
        if residual_data.npol > 1 and any(component_prior(c).is_hier for c in self.active_comps):
            raise ValueError("The hier amp_prior is intensity-only in this version (its template "
                             "must be strictly positive; Q and U are signed). Give the QU view "
                             "amp_prior none or gaussian.")
        self.region = region_map(layout.nside_beta, nside, keep_mask)
        for ipol in range(residual_data.npol):
            freqs, maps_sky, maps_inv_var = [], [], []
            for iband in range(len(all_freq)):
                if all_map_sky[iband][ipol] is not None:
                    freqs.append(all_freq[iband])
                    maps_sky.append(all_map_sky[iband][ipol])
                    maps_inv_var.append(all_inv_var[iband][ipol])
            self.streams.append((np.array(freqs, dtype=np.float64),
                                 np.array(maps_sky, dtype=np.float64),
                                 np.array(maps_inv_var, dtype=np.float64)))

    # -- likelihood part on this rank's streams ----------------------------------------------

    def _hier_views(self, hypers: dict[str, tuple[float, float]]) -> list[int]:
        """Indices of this rank's active views that carry a hier prior with sampled hypers."""
        return [c for c, comp in enumerate(self.active_comps)
                if comp.comp_name in hypers and component_prior(comp).is_hier]

    def _set_hypers(self, hypers: dict[str, tuple[float, float]]) -> dict[int, tuple]:
        """Write the evaluated (mu, log s) into the hier priors; returns the previous values so the
        caller can restore them (the components must not keep a rejected trajectory point)."""
        previous = {}
        for c in self._hier_views(hypers):
            prior = self.active_comps[c].amp_prior
            previous[c] = (prior.mu, prior.log_s)
            prior.mu, prior.log_s = hypers[self.active_comps[c].comp_name]
        return previous

    def _restore_hypers(self, previous: dict[int, tuple]) -> None:
        for c, (mu, log_s) in previous.items():
            self.active_comps[c].amp_prior.mu, self.active_comps[c].amp_prior.log_s = mu, log_s

    def _local_terms(self, betas: dict[str, NDArray], hypers: dict[str, tuple[float, float]]
                     ) -> tuple[float, dict[str, NDArray], dict[str, NDArray], bool]:
        """(log L, dlog L/dbeta per sampled component (nreg,), (dF/dmu, dF/dlog s) per hier
        component, rejected) summed over this rank's streams."""
        layout = self.layout
        logL = 0.0
        g_beta = {name: np.zeros(layout.nreg) for name in layout.index_slices}
        g_hyper = {name: np.zeros(2) for name in layout.hier_names}
        rejected = False
        self.nskipped_last = 0
        if not self.is_solver_root:
            return logL, g_beta, g_hyper, rejected
        previous = self._set_hypers(hypers)
        try:
            return self._stream_terms(betas, hypers, logL, g_beta, g_hyper)
        finally:
            self._restore_hypers(previous)

    def _stream_terms(self, betas, hypers, logL, g_beta, g_hyper):
        layout = self.layout
        comps = self.active_comps
        ncomp = len(comps)
        hier_cols = self._hier_views(hypers)
        any_hier = bool(hier_cols)
        for freqs, maps_sky, maps_inv_var in self.streams:
            nband = freqs.size
            M_reg = np.empty((layout.nreg, nband, ncomp))
            dM_reg = np.zeros_like(M_reg)
            with np.errstate(over="ignore", invalid="ignore"):
                for c, comp in enumerate(comps):
                    if comp.comp_name in betas:
                        M_reg[:, :, c] = comp.get_sed(freqs, beta=betas[comp.comp_name])
                        dM_reg[:, :, c] = comp.get_dsed_dbeta(freqs, beta=betas[comp.comp_name])
                    else:
                        M_reg[:, :, c] = comp.get_sed(freqs)[None, :]
            if not (np.all(np.isfinite(M_reg)) and np.all(np.isfinite(dM_reg))):
                return logL, g_beta, g_hyper, True   # zero density at this state
            sky_aug, w_aug = augment_stream(maps_sky, maps_inv_var, comps)
            M_aug, dM_aug = augment_mixing(M_reg, dM_reg)
            out = cpp_compsep.marginal_loglike_regions(
                sky_aug, w_aug, self.region, M_aug, dM_aug, self.nthreads, True, True, any_hier)
            ll_reg, grad_reg, nskipped = out[0], out[1], out[2]
            self.nskipped_last += int(nskipped)
            if nskipped or not np.all(np.isfinite(ll_reg)) or not np.all(np.isfinite(grad_reg)):
                return logL, g_beta, g_hyper, True
            logL += float(ll_reg.sum())
            for c, comp in enumerate(comps):
                if comp.comp_name in betas:
                    g_beta[comp.comp_name] += grad_reg[:, c]
            if any_hier:
                x_cond, diag_Ainv, used = out[3], out[4], out[5]
                for c in hier_cols:
                    prior = comps[c].amp_prior
                    logL += prior.trailing_terms(used)
                    g_hyper[comps[c].comp_name] += prior.hyper_grad(x_cond[c], diag_Ainv[c], used)
        return logL, g_beta, g_hyper, False

    # -- the collective evaluation -------------------------------------------------------------

    def logp_and_grad(self, x: NDArray) -> tuple[float, NDArray]:
        layout = self.layout
        x = np.asarray(x, dtype=np.float64)
        if x.shape != (layout.dim,):
            raise ValueError(f"x must have shape ({layout.dim},), got {x.shape}")
        hypers = layout.hypers(x)
        guard = any(log_s < HIER_LOGS_MIN for _, log_s in hypers.values())
        rejected = guard or not np.all(np.isfinite(x))
        if rejected:
            betas = {}
        else:
            with np.errstate(over="ignore", invalid="ignore"):
                betas = layout.betas(x)
            rejected = not all(np.all(np.isfinite(b)) for b in betas.values())
        if rejected:
            logL, g_beta, g_hyper = 0.0, {n: np.zeros(layout.nreg) for n in layout.index_slices}, \
                {n: np.zeros(2) for n in layout.hier_names}
        else:
            logL, g_beta, g_hyper, rejected = self._local_terms(betas, hypers)
        # One allreduce: [rejected flag, log L, g_beta blocks..., g_hyper blocks...].
        local = np.concatenate([[1.0 if rejected else 0.0, logL]]
                               + [g_beta[s.comp_name] for s in layout.specs]
                               + [g_hyper[n] for n in layout.hier_names])
        total = np.empty_like(local)
        self.compsep_comm.Allreduce(local, total, op=MPI.SUM)
        if total[0] > 0.0 or not np.isfinite(total[1]):
            self.n_rejected += 1
            return -np.inf, np.zeros(layout.dim)
        logp = float(total[1])
        grad = np.zeros(layout.dim)
        pos = 2
        for spec in layout.specs:
            g = total[pos:pos + layout.nreg]
            pos += layout.nreg
            sl = layout.index_slices[spec.comp_name]
            if spec.field is None:
                r = (float(x[sl][0]) - spec.prior_mean) / spec.prior_sd
                logp += -0.5 * r * r
                grad[sl] = g.sum() - r / spec.prior_sd
            else:
                # Field prior, hyperpriors, Jacobians and the chain rule from dlog L/dbeta to the
                # sampled coordinates, with the likelihood value already counted above.
                lp_f, g_f = spec.field.logp_and_grad(x[sl], loglike=lambda _b, g=g: (0.0, g))
                if not np.isfinite(lp_f):
                    self.n_rejected += 1
                    return -np.inf, np.zeros(layout.dim)
                logp += lp_f
                grad[sl] = g_f
        for name in layout.hier_names:
            sl = layout.hier_slices[name]
            # A copy carrying the evaluated hyperparameters: the component's own prior keeps
            # its state (a trajectory point is not a sample).
            prior = replace(self._hier_prior(name), mu=hypers[name][0], log_s=hypers[name][1])
            lp_h = prior.log_hyperprior()
            if not np.isfinite(lp_h):
                self.n_rejected += 1
                return -np.inf, np.zeros(layout.dim)
            logp += lp_h
            grad[sl] = total[pos:pos + 2] + np.array(prior.log_hyperprior_grad())
            pos += 2
        if not np.isfinite(logp):
            self.n_rejected += 1
            return -np.inf, np.zeros(layout.dim)
        return logp, grad

    def _hier_prior(self, name: str) -> AmpPrior:
        """The hier prior object of a component by name, with its hyperparameters at the state
        being evaluated. Ranks without that component's view evaluate the hyperprior through a
        stand-in carrying the same hyperprior settings, which every rank must be given."""
        for comp in self.active_comps:
            if comp.comp_name == name and component_prior(comp).is_hier:
                return comp.amp_prior
        return self._hier_defaults[name]

    def set_hier_defaults(self, defaults: dict[str, AmpPrior]) -> None:
        """Hier prior objects (hyperprior settings only) for components this rank does not hold,
        so the hyperprior terms can be evaluated on every rank; their (mu, log s) are set from the
        state before use."""
        self._hier_defaults = dict(defaults)

    @classmethod
    def from_arrays(cls, freqs: NDArray, maps_sky: NDArray, maps_inv_var: NDArray,
                    comps: list[DiffuseComponent], layout: StateLayout, nside: int,
                    keep_mask: NDArray | None = None, nthreads: int = 1,
                    comm: MPI.Comm = MPI.COMM_SELF) -> "RegionIndexLikelihood":
        """In-process construction from one intensity stream (nband, npix) already gathered, for
        tests and the L1 harness; the evaluation is the same as in the MPI run."""
        obj = cls(comm, None, None, comps, layout, keep_mask, nthreads)
        obj.is_solver_root = True
        obj.region = region_map(layout.nside_beta, nside, keep_mask)
        obj.streams = [(np.asarray(freqs, dtype=np.float64),
                        np.asarray(maps_sky, dtype=np.float64),
                        np.asarray(maps_inv_var, dtype=np.float64))]
        return obj
