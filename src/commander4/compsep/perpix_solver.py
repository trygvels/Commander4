"""The per-pixel amplitude solver: an independent least-squares fit in each pixel.

Alternative to the global CG solve for the case where every band has been brought to a common
resolution. The mixing matrix is then the same at every multipole, so the amplitudes follow from a
small (nband x ncomp) normal-equation solve per pixel, done in C++. Unlike the CG solver, the
amplitudes it produces carry the common beam rather than being deconvolved.
"""
import time
import logging
import numpy as np
from mpi4py import MPI
from pixell import curvedsky

from commander4.sky.component import Component
from commander4.sky.diffuse_components import DiffuseComponent
from commander4.backend import compsep as cpp_compsep
from commander4.data_models.detector_map import DetectorMap
from commander4.diagnostics.performance import benchmark


def solve_compsep_perpix(proc_comm: MPI.Comm, detector_data: DetectorMap,
                         comp_list: list[Component],
                         double_precision: bool, nthreads: int,
                         rng: np.random.Generator) -> list[Component]:
    """ A pixel-by-pixel solver for the component separation problem. Requires uniform nside and a
        common beam across all bands (smoothing to a common resolution is done at the data sources,
        controlled by ``compsep.common_res_fwhm``). This solver ignores beams entirely, so bands
        arriving at differing resolutions are silently mixed; a warning is logged below.
        The per-pixel solves run on the master rank only, with ``nthreads`` threads, and draw
        their fluctuation terms from the master's ``rng``.
    """
    for component in comp_list:
        if not isinstance(component, DiffuseComponent):
            raise ValueError(
                f"Per-pixel component separation does not support "
                f"{type(component).__name__} components."
            )
    logger = logging.getLogger(__name__)
    if proc_comm.Get_rank() == 0:
        logger.verbose("Starting pixel-by-pixel component separation.")
    if double_precision:
        complex_dtype = np.complex128
        real_dtype = np.float64
    else:
        complex_dtype = np.complex64
        real_dtype = np.float32

    npol = detector_data.npol
    pol = detector_data.pol
    spin = 2 if pol else 0
    map_sky = detector_data.map_sky
    band_freq = detector_data.nu
    map_inv_var = detector_data.map_inv_var  # Zero where unobserved.

    # This solver has no beam model, so all bands should already share one resolution (smoothed at
    # ingest via compsep.common_res_fwhm); differing FWHMs silently mix resolutions, so warn.
    all_fwhm = proc_comm.allgather(detector_data.fwhm)
    if not np.allclose(all_fwhm, all_fwhm[0]):
        logger.warning(f"Per-pixel solver received bands at differing resolutions {all_fwhm}; the "
                       "result mixes resolutions. Set compsep.common_res_fwhm to a common beam.")

    ncomp = len(comp_list)
    with benchmark("perpix-gather"):
        all_freq = proc_comm.gather(band_freq, root=0)
        all_map_sky = proc_comm.gather(map_sky, root=0)
        all_map_inv_var = proc_comm.gather(map_inv_var, root=0)

    nside = detector_data.nside
    npix = 12*nside**2
    comp_maps = [None, None] if pol else [None]
    if proc_comm.Get_rank() == 0:  # Unfortunately, only master rank does the work.
        map_shapes = np.array([_map.shape for _map in all_map_sky])
        if not np.all(map_shapes == map_shapes[0]):
            raise ValueError("Per-pixel solver requires all maps to have the same shape, but "
                             f"received {map_shapes}.")
        for ipol in range(npol):
            t0 = time.time()
            freqs = []
            maps_sky = []
            maps_inv_var = []
            for iband in range(len(all_freq)):
                if all_map_sky[iband][ipol] is not None:
                    freqs.append(all_freq[iband])
                    maps_sky.append(all_map_sky[iband][ipol])
                    maps_inv_var.append(all_map_inv_var[iband][ipol])
            freqs = np.array(freqs)
            maps_sky = np.array(maps_sky)
            maps_inv_var = np.array(maps_inv_var)
            nband = len(freqs)
            M = np.empty((nband, ncomp))
            for icomp in range(ncomp):
                M[:, icomp] = comp_list[icomp].get_sed(freqs)
            rand = rng.standard_normal((npix, nband))
            # TODO: Should scale M to make solution more well-conditioned, and then adjust
            # solution with the scaling factor used.
            with benchmark("perpix-solve"):
                comp_maps[ipol] = cpp_compsep.solve_perpix(
                    maps_sky.astype(np.float64, copy=False),
                    maps_inv_var.astype(np.float64, copy=False), M, rand, nthreads)
            logger.info(f"Finished pixel-by-pixel component separation in {time.time()-t0:.2f}s "\
                        f"for polarization {ipol+1} of {npol}.")

    comp_maps = proc_comm.bcast(comp_maps, root=0)
    with benchmark("perpix-alm"):
        for icomp in range(ncomp):
            if pol:
                input_map = np.array([comp_maps[0][icomp], comp_maps[1][icomp]], dtype=real_dtype)
            else:
                input_map = np.array([comp_maps[0][icomp]], dtype=real_dtype)
            comp_alms = curvedsky.map2alm_healpix(input_map, niter=3, spin=spin,
                                                  lmax=comp_list[icomp].lmax)
            comp_list[icomp].alms = comp_alms.astype(complex_dtype, copy=False)
            # These amplitudes are component maps at the data resolution; record that beam so the
            # forward model applies no extra smoothing when predicting a band at this resolution.
            comp_list[icomp].amp_fwhm_rad = detector_data.fwhm_rad

    return comp_list


class PerPixelMarginalLikelihood:
    """Log likelihood of the spectral parameters with the per-pixel amplitudes integrated out.

    For the per-pixel data model (common resolution, diagonal noise, flat amplitude prior; the
    model `solve_compsep_perpix` samples from) the amplitudes of the active components can be
    marginalized analytically in every pixel, see the `marginal_loglike_perpix` kernel. An MCMC
    group accepting on this quantity, and redrawing the amplitudes from their conditional after its
    steps, samples the joint posterior exactly without an amplitude solve per proposal.

    Construction is collective over ``compsep_comm``: the residual band maps (data minus the
    components the amplitude group holds fixed) are gathered once to the root of ``solver_comm``.
    This is valid for as long as the fixed components do not change, i.e. for one MCMC group run
    whose sampled components are all active in the amplitude group (`init_compsep_processing`
    enforces this before selecting the marginal likelihood). Calling the object is also collective
    over ``compsep_comm``; it returns the summed log likelihood of both polarization streams.

    Args:
        compsep_comm: Communicator spanning every CompSep rank (both polarization streams).
        solver_comm: This rank's per-polarization communicator of bands taking part in the amplitude
            group, or None if this rank's band does not take part.
        residual_data: This band's residual map (fixed components subtracted). Ignored when
            ``solver_comm`` is None.
        active_comps: The amplitude group's active components for this rank's polarization stream.
        keep_mask: Optional boolean pixel mask (True = keep) applied to every band.
        nthreads: Threads for the kernel on the solver root.
        use_logdet: Keep the -1/2 ln det term (the marginal likelihood, flat amplitude prior) or
            drop it (the ridge / profile likelihood, see the kernel's docstring).
    """

    def __init__(self, compsep_comm: MPI.Comm, solver_comm: MPI.Comm | None,
                 residual_data: DetectorMap | None, active_comps: list[Component],
                 keep_mask: np.ndarray | None, nthreads: int, use_logdet: bool = True):
        self.compsep_comm = compsep_comm
        self.use_logdet = use_logdet
        self.active_comps = list(active_comps)
        self.nthreads = nthreads
        self.is_solver_root = solver_comm is not None and solver_comm.Get_rank() == 0
        self.nskipped_warned = False
        self.streams = []  # One (freqs, maps, inv_var) entry per polarization, on the root only.
        for component in self.active_comps:
            if not isinstance(component, DiffuseComponent):
                raise ValueError(f"The per-pixel marginal likelihood does not support "
                                 f"{type(component).__name__} components.")
        if solver_comm is None:
            return
        inv_var = np.array(residual_data.map_inv_var, dtype=np.float64, copy=True)
        if keep_mask is not None:
            inv_var[:, ~keep_mask] = 0.0
        all_freq = solver_comm.gather(residual_data.nu, root=0)
        all_map_sky = solver_comm.gather(residual_data.map_sky, root=0)
        all_inv_var = solver_comm.gather(inv_var, root=0)
        if not self.is_solver_root:
            return
        for ipol in range(residual_data.npol):
            freqs, maps_sky, maps_inv_var = [], [], []
            for iband in range(len(all_freq)):
                if all_map_sky[iband][ipol] is not None:
                    freqs.append(all_freq[iband])
                    maps_sky.append(all_map_sky[iband][ipol])
                    maps_inv_var.append(all_inv_var[iband][ipol])
            self.streams.append((np.array(freqs),
                                 np.array(maps_sky, dtype=np.float64),
                                 np.array(maps_inv_var, dtype=np.float64)))

    def mixing_matrix(self, freqs: np.ndarray) -> np.ndarray:
        """(nband, ncomp) mixing matrix of the active components at their current parameters."""
        return np.stack([comp.get_sed(freqs) for comp in self.active_comps], axis=1)

    def __call__(self) -> float:
        local = 0.0
        if self.is_solver_root:
            for freqs, maps_sky, maps_inv_var in self.streams:
                M = np.ascontiguousarray(self.mixing_matrix(freqs), dtype=np.float64)
                loglike, _, nskipped = cpp_compsep.marginal_loglike_perpix(
                    maps_sky, maps_inv_var, M, self.nthreads, self.use_logdet)
                if nskipped and not self.nskipped_warned:
                    logging.getLogger(__name__).warning(
                        f"Per-pixel marginal likelihood skipped {nskipped} observed pixel(s) whose "
                        "normal matrix is not positive definite at the current parameters; the "
                        "MCMC acceptance ratio is then not exact. Check the mixing matrix scaling.")
                    self.nskipped_warned = True
                local += loglike
        return float(self.compsep_comm.allreduce(local, op=MPI.SUM))
