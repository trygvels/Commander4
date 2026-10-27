"""The collective evaluation of RegionIndexLikelihood under MPI (from the 2026-10-27 stats
review). Meaningful under

    mpirun -n 4 python -m pytest tests/test_index_field_mpi.py

(three one-band ranks plus a rank without a stream); it also runs, trivially, on one rank."""
import numpy as np
from mpi4py import MPI
from pixell.bunch import Bunch

from commander4.compsep import index_field as ifd
from commander4.data_models.detector_map import DetectorMap
from commander4.math_utils import gmrf_field
from commander4.sky.diffuse_components import CMB, Synchrotron, ThermalDust

NSIDE, NSIDE_BETA = 4, 1
NPIX = 12 * NSIDE**2
FREQS = np.array([30.0, 70.0, 353.0])


def _comp(cls, name, **sed):
    p = Bunch(polarization="I", shortname=name, spatially_varying_MM=False, Cl_prior_amplitude=None,
              lmax=2, units="uK_RJ", **sed)
    object.__setattr__(p, "_name", name)
    return cls(p, Bunch(nside=NSIDE, double_precision=True), allocate_empty_alms=True,
               eval_pol="I", comp_name=name)


def test_distributed_evaluation_equals_single_process() -> None:
    comm = MPI.COMM_WORLD
    rank, size = comm.Get_rank(), comm.Get_size()
    rng = np.random.default_rng(5)                    # identical on every rank
    g_s = 10**(1 + 0.5 * rng.random(NPIX))
    g_d = 10**(1.5 + 0.5 * rng.random(NPIX))
    nband = min(FREQS.size, size)
    freqs = FREQS[:nband]
    d = rng.normal(size=(FREQS.size, NPIX)) * 30 + 50
    w = 1.0 / (5.0 * (0.5 + rng.random((FREQS.size, NPIX))))**2
    w[:, 0] = 0.0
    w[1:, 1] = 0.0
    keep = rng.random(NPIX) > 0.3

    def make_comps():
        cs = [_comp(CMB, "cmb", nu_ref=1e-3), _comp(Synchrotron, "sync", beta=-3.0, nu_ref=30.0),
              _comp(ThermalDust, "dust", beta=1.6, T=20.0, nu_ref=545.0)]
        cs[0].amp_prior = ifd.AmpPrior("gaussian", mean=0.0, rms=70.0)
        cs[1].amp_prior = ifd.AmpPrior("hier", template=g_s)
        cs[2].amp_prior = ifd.AmpPrior("hier", template=g_d)
        return cs

    field = gmrf_field.IndexFieldPrior(
        gmrf_field.MaternField(gmrf_field.GraphLaplacianBasis(NSIDE_BETA)),
        gmrf_field.HyperPrior(mu0=-3.0, s_mu=0.2))
    layout = ifd.StateLayout([ifd.IndexSpec("sync", field, -3.0, 0.2)], ["sync", "dust"])
    ref = ifd.RegionIndexLikelihood.from_arrays(freqs, d[:nband], w[:nband], make_comps(), layout,
                                                NSIDE, keep_mask=keep)
    holds = rank < nband
    solver_comm = comm.Split(0 if holds else MPI.UNDEFINED, rank)
    if holds:
        data = DetectorMap(d[rank], w[rank], nu=FREQS[rank], fwhm=60.0, nside=NSIDE,
                           double_precision=True)
        like = ifd.RegionIndexLikelihood(comm, solver_comm, data, make_comps(), layout, keep, 1)
        solver_comm.Free()
    else:
        like = ifd.RegionIndexLikelihood(comm, None, None, [], layout, keep, 1)
    like.set_hier_defaults({"sync": ifd.AmpPrior("hier", template=np.ones(1)),
                            "dust": ifd.AmpPrior("hier", template=np.ones(1))})
    rs = np.random.default_rng(9)
    for _ in range(3):
        x = np.concatenate([rs.normal(size=field.K),
                            [-3.0 + 0.1 * rs.normal(), np.log(0.2), rs.normal()],
                            np.array([1.0, np.log(0.3), 1.0, np.log(0.3)]) + 0.2 * rs.normal(size=4)])
        lp, g = like.logp_and_grad(x)
        lp_r, g_r = ref.logp_and_grad(x)
        assert abs(lp - lp_r) <= 1e-9 * abs(lp_r)
        np.testing.assert_allclose(g, g_r, rtol=1e-9, atol=1e-9)
    bad = x.copy()
    bad[-1] = -9.0
    assert like.logp_and_grad(bad)[0] == -np.inf
