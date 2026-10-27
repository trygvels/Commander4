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

from dataclasses import dataclass, field

import healpy as hp
import numpy as np
from numpy.typing import NDArray

from commander4.sky.diffuse_components import DiffuseComponent

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
            m = np.broadcast_to(np.asarray(self.mean, dtype=np.float64), (npix,)).copy()
            t = np.broadcast_to(np.asarray(self.rms, dtype=np.float64), (npix,))
            return m, 1.0 / t**2
        g = self.template
        if g.shape != (npix,):
            raise ValueError(f"hier template has shape {g.shape}, expected ({npix},).")
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
        P, g = P[used], self.template[used]
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
