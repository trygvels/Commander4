"""The MCMC acceptance likelihoods: the per-pixel marginal kernel, the step loop of the
'marginal' and 'conditional' modes, how each group's mode is resolved, and a statistical check
that both modes sample the exact posterior of a spectral index.
"""
import numpy as np
import pytest
from mpi4py import MPI
from pixell.bunch import Bunch

from commander4.backend import compsep as cpp_compsep
from commander4.compsep.mcmc import MCMCSamplingGroup
from commander4.compsep.processing import (
    CGSamplingGroupConfig,
    MCMCSamplingGroupConfig,
    PerPixelSamplingGroupConfig,
    _resolve_mcmc_likelihoods,
)
from commander4.data_models.detector_map import DetectorMap
from commander4.sky.comp_list import CompList


# ---------------------------------------------------------------- kernel

def _reference_marginal(d, w, M, use_logdet=True):
    """Independent NumPy reference: sum over well-determined pixels of b^T A^-1 b/2 - ln det A/2."""
    total, nused = 0.0, 0
    ncomp = M.shape[1]
    for p in range(d.shape[1]):
        if np.count_nonzero(w[:, p] > 0) < ncomp:
            continue
        A = (M.T * w[:, p]) @ M
        b = M.T @ (w[:, p] * d[:, p])
        total += 0.5 * b @ np.linalg.solve(A, b)
        if use_logdet:
            total -= 0.5 * np.linalg.slogdet(A)[1]
        nused += 1
    return total, nused


@pytest.mark.parametrize("nthreads", [1, 3])
def test_marginal_kernel_matches_numpy_reference(nthreads):
    rng = np.random.default_rng(1)
    nband, npix, ncomp = 7, 300, 4
    M = rng.uniform(0.1, 3.0, (nband, ncomp))
    d = rng.normal(size=(nband, npix)) * 10
    w = rng.uniform(0.2, 4.0, (nband, npix))
    w[:, :7] = 0.0            # unobserved pixels
    w[:4, 7:20] = 0.0         # observed by fewer bands than components: never in the sum
    loglike, nused, nskipped = cpp_compsep.marginal_loglike_perpix(d, w, M, nthreads)
    ref, nref = _reference_marginal(d, w, M)
    assert nused == nref == npix - 20
    assert nskipped == 0
    assert loglike == pytest.approx(ref, rel=1e-12, abs=1e-9)


def test_ridge_kernel_drops_only_the_log_determinant():
    rng = np.random.default_rng(5)
    M = rng.uniform(0.1, 3.0, (6, 3))
    d = rng.normal(size=(6, 80)) * 5
    w = rng.uniform(0.5, 2.0, (6, 80))
    ridge = cpp_compsep.marginal_loglike_perpix(d, w, M, 2, False)[0]
    assert ridge == pytest.approx(_reference_marginal(d, w, M, use_logdet=False)[0], rel=1e-12)


def test_marginal_kernel_ignores_masked_bands_in_a_pixel():
    rng = np.random.default_rng(2)
    M = rng.uniform(0.5, 2.0, (5, 2))
    d = rng.normal(size=(5, 10))
    w = np.ones((5, 10))
    w_masked = w.copy()
    w_masked[0, :] = 0.0
    d_garbage = d.copy()
    d_garbage[0, :] = 1e30  # A band with zero weight must not enter at all.
    a = cpp_compsep.marginal_loglike_perpix(d_garbage, w_masked, M, 1)[0]
    b = cpp_compsep.marginal_loglike_perpix(d[1:], w[1:], M[1:], 1)[0]
    assert a == pytest.approx(b)


# ---------------------------------------------------------------- toy model shared by the tests

class PerPixelToy:
    """Per-pixel data d = M(beta) a + n over a few bands, three components, one scalar index.

    Component SEDs: flat (CMB-like), a fixed power law, and a power law with the sampled index.
    """

    def __init__(self, npix=150, beta_true=1.6, seed=3, noise=1.0):
        rng = np.random.default_rng(seed)
        self.freqs = np.array([30.0, 44.0, 70.0, 100.0, 143.0, 217.0, 353.0])
        self.w = np.full((self.freqs.size, npix), 1.0 / noise**2)
        a_true = np.stack([rng.normal(0, 20, npix),
                           10**rng.uniform(0, 2, npix),
                           10**rng.uniform(0.5, 2, npix)])
        self.d = self.mixing(beta_true) @ a_true + rng.normal(size=(self.freqs.size, npix)) * noise
        self.beta = 1.5
        self.prior = (1.5, 0.3)
        self.amps = np.zeros((3, npix))
        self.rng = np.random.default_rng(seed + 100)

    def mixing(self, beta):
        nu = self.freqs
        return np.stack([np.ones_like(nu), (nu / 30.0)**-3.0, (nu / 353.0)**beta], axis=1)

    def marginal(self, use_logdet=True):
        return cpp_compsep.marginal_loglike_perpix(self.d, self.w, self.mixing(self.beta), 1,
                                                   use_logdet)[0]

    def ridge(self):
        return self.marginal(use_logdet=False)

    def chisq_loglike(self):
        r = self.d - self.mixing(self.beta) @ self.amps
        return -0.5 * float(np.sum(self.w * r**2))

    def resolve(self):
        """Exact draw of every pixel's amplitudes from p(a | beta, d), written in place."""
        rand = self.rng.standard_normal((self.d.shape[1], self.d.shape[0]))
        np.copyto(self.amps, cpp_compsep.solve_perpix(
            self.d, self.w, self.mixing(self.beta), rand, 1))

    def _logpost(self, betas, use_logdet=True):
        """Unnormalized log posterior of beta, vectorized NumPy (independent of the C++ kernel)."""
        out = []
        for beta in betas:
            M = self.mixing(beta)
            A = np.einsum("bp,bi,bj->pij", self.w, M, M)
            b = np.einsum("bp,bi->pi", self.w * self.d, M)
            quad = np.einsum("pi,pi->p", b, np.linalg.solve(A, b[..., None])[..., 0])
            logdet = np.linalg.slogdet(A)[1] if use_logdet else 0.0
            out.append(np.sum(0.5 * quad - 0.5 * logdet)
                       - 0.5 * ((beta - self.prior[0]) / self.prior[1])**2)
        return np.array(out)

    def exact_posterior(self, use_logdet=True):
        """Mean and sd of beta by quadrature on a grid around the mode."""
        grid = np.linspace(self.prior[0] - 1.0, self.prior[0] + 1.0, 201)
        peak = grid[np.argmax(self._logpost(grid, use_logdet))]
        for width in (0.05, None):  # zoom in once, then integrate over +-10 sd
            if width is None:
                width = 10 * sd
            grid = np.linspace(peak - width, peak + width, 801)
            logp = self._logpost(grid, use_logdet)
            p = np.exp(logp - logp.max())
            p /= p.sum()
            mean = np.sum(p * grid)
            sd = np.sqrt(np.sum(p * (grid - mean)**2))
            peak = mean
        return mean, sd


class ScalarIndexGroup(MCMCSamplingGroup):
    """An MCMC group over the toy's scalar index, on COMM_SELF."""

    class _AmpHolder:  # What the original 'resolve_chisq' path buffers and restores.
        def __init__(self, data):
            self._data = data

    def __init__(self, toy, sigma):
        dummy = DetectorMap(np.zeros((1, 12)), np.ones((1, 12)), nu=100.0, fwhm=0.0, nside=1)
        super().__init__(MPI.COMM_SELF, dummy, [self._AmpHolder(toy.amps)], target_pol="I",
                         chisq_active=True)
        self.toy, self.sigma = toy, sigma

    def has_parameters(self):
        return True

    def capture_state(self):
        return float(self.toy.beta)

    def propose(self, current_state, rng):
        return float(rng.normal(current_state, self.sigma)), True

    def apply_state(self, state):
        self.toy.beta = state

    def log_prior(self, state):
        return -0.5 * ((state - self.toy.prior[0]) / self.toy.prior[1])**2

    def global_loglike(self):
        return self.toy.chisq_loglike()


# ---------------------------------------------------------------- step loop

class CountingCalls:
    def __init__(self):
        self.n = 0

    def __call__(self):
        self.n += 1


def test_marginal_mode_never_resolves_inside_the_loop_and_once_after_acceptance():
    toy = PerPixelToy(npix=40)
    toy.resolve()
    group = ScalarIndexGroup(toy, sigma=0.01)
    resolve = CountingCalls()
    calls = CountingCalls()

    def marginal():
        calls()
        return toy.marginal()

    stats = group.run(numstep=25, resolve_amplitudes=resolve, rng=np.random.default_rng(0),
                      likelihood="marginal", marginal_loglike=marginal)
    assert stats["n_accept"] > 0
    assert resolve.n == 1          # one redraw for the final accepted state
    assert calls.n == 26           # once for the starting state, once per proposal


def test_conditional_mode_keeps_amplitudes_fixed_while_stepping():
    toy = PerPixelToy(npix=40)
    toy.resolve()
    amps_before = toy.amps.copy()
    group = ScalarIndexGroup(toy, sigma=1e-4)

    def resolve():
        np.testing.assert_array_equal(toy.amps, amps_before)  # untouched until the steps end
        toy.resolve()

    stats = group.run(numstep=10, resolve_amplitudes=resolve, rng=np.random.default_rng(1),
                      likelihood="conditional")
    assert stats["n_accept"] > 0


def test_no_acceptance_means_no_amplitude_redraw():
    toy = PerPixelToy(npix=40)
    toy.resolve()
    group = ScalarIndexGroup(toy, sigma=50.0)  # every proposal lands far out and is rejected
    resolve = CountingCalls()
    beta0 = toy.beta
    stats = group.run(numstep=10, resolve_amplitudes=resolve, rng=np.random.default_rng(2),
                      likelihood="marginal", marginal_loglike=toy.marginal)
    assert stats["n_accept"] == 0
    assert resolve.n == 0
    assert toy.beta == beta0


def test_marginal_mode_requires_the_callable():
    group = ScalarIndexGroup(PerPixelToy(npix=10), sigma=0.01)
    with pytest.raises(ValueError, match="marginal_loglike"):
        group.run(numstep=1, resolve_amplitudes=lambda: None, rng=np.random.default_rng(0),
                  likelihood="marginal")


# ---------------------------------------------------------------- statistical correctness

def _run_chain(toy, likelihood, n_iter, numstep, sigma, seed):
    group = ScalarIndexGroup(toy, sigma=sigma)
    rng = np.random.default_rng(seed)
    samples = np.empty(n_iter)
    for it in range(n_iter):
        toy.resolve()  # the amplitude group, run first in every Gibbs iteration
        marginal = {"marginal": toy.marginal, "ridge": toy.ridge}.get(likelihood)
        group.run(numstep=numstep, resolve_amplitudes=toy.resolve, rng=rng,
                  likelihood=likelihood, marginal_loglike=marginal)
        samples[it] = toy.beta
    return samples


@pytest.mark.parametrize("likelihood,n_iter,numstep", [("marginal", 1500, 4),
                                                        ("ridge", 1500, 4),
                                                        ("conditional", 2500, 4)])
def test_valid_modes_sample_the_exact_posterior(likelihood, n_iter, numstep):
    toy = PerPixelToy(npix=150)
    mean, sd = toy.exact_posterior(use_logdet=likelihood != "ridge")
    samples = _run_chain(toy, likelihood, n_iter, numstep, sigma=1.5 * sd, seed=7)
    kept = samples[n_iter // 5:]
    assert abs(kept.mean() - mean) < 0.3 * sd
    assert 0.8 < kept.std() / sd < 1.2


# ---------------------------------------------------------------- mode resolution

def _spectral_comp_list() -> CompList:
    params = Bunch(compsep=Bunch(nside=2, double_precision=False))
    sync = Bunch(enabled=True, component_class="Synchrotron", params=Bunch(
        lmax=2, polarization="I", shortname="sync", spatially_varying_MM=False,
        Cl_prior_amplitude=None, beta=-3.1, nu_ref=30.0, sample_spectral_index=True,
        spectral_index_proposal_sigma=0.02))
    object.__setattr__(sync, "_name", "Synchrotron")
    cmb = Bunch(enabled=True, component_class="CMB", params=Bunch(
        lmax=2, polarization="I", shortname="cmb", spatially_varying_MM=False,
        Cl_prior_amplitude=None))
    object.__setattr__(cmb, "_name", "CMB")
    return CompList.init_from_params(Bunch({"CMB": cmb, "Synchrotron": sync}), params)


def _resolve(mcmc_cfg, amp_groups, method):
    return _resolve_mcmc_likelihoods({"m": mcmc_cfg}, amp_groups, method, _spectral_comp_list(),
                                     report=False)["m"]


def test_auto_picks_marginal_for_one_per_pixel_group_covering_the_sampled_components():
    amp = {"pp": PerPixelSamplingGroupConfig(name="pp")}
    cfg = MCMCSamplingGroupConfig(name="m", update_amplitude_groups=("pp",))
    assert _resolve(cfg, amp, "per_pixel") == "marginal"


def test_auto_falls_back_to_conditional():
    cfg = MCMCSamplingGroupConfig(name="m", update_amplitude_groups=("cg",))
    assert _resolve(cfg, {"cg": CGSamplingGroupConfig(name="cg")}, "cg") == "conditional"
    assert _resolve(MCMCSamplingGroupConfig(name="m"), {}, "per_pixel") == "conditional"
    amp = {"pp": PerPixelSamplingGroupConfig(name="pp", comps=("CMB",))}  # sync not solved
    cfg = MCMCSamplingGroupConfig(name="m", update_amplitude_groups=("pp",))
    assert _resolve(cfg, amp, "per_pixel") == "conditional"


def test_explicit_marginal_that_cannot_be_honoured_is_an_error():
    cfg = MCMCSamplingGroupConfig(name="m", update_amplitude_groups=("cg",), likelihood="marginal")
    with pytest.raises(ValueError, match="unavailable"):
        _resolve(cfg, {"cg": CGSamplingGroupConfig(name="cg")}, "cg")


def test_explicit_ridge_needs_the_per_pixel_marginal_too():
    cfg = MCMCSamplingGroupConfig(name="m", update_amplitude_groups=("cg",), likelihood="ridge")
    with pytest.raises(ValueError, match="unavailable"):
        _resolve(cfg, {"cg": CGSamplingGroupConfig(name="cg")}, "cg")


def test_explicit_settings_are_kept_and_unknown_values_rejected():
    amp = {"pp": PerPixelSamplingGroupConfig(name="pp")}
    for value in ("ridge", "conditional", "resolve_chisq"):
        cfg = MCMCSamplingGroupConfig(name="m", update_amplitude_groups=("pp",), likelihood=value)
        assert _resolve(cfg, amp, "per_pixel") == value
    with pytest.raises(ValueError, match="likelihood"):
        MCMCSamplingGroupConfig(name="m", likelihood="profile")


def test_original_resolve_chisq_mode_does_not_sample_the_posterior():
    """Documents why 'resolve_chisq' is no longer the default: accepting on the chi-squared at
    freshly drawn amplitudes widens and shifts the posterior (here by more than 2x in width)."""
    toy = PerPixelToy(npix=150)
    mean, sd = toy.exact_posterior()
    samples = _run_chain(toy, "resolve_chisq", 2500, 4, sigma=1.5 * sd, seed=7)[500:]
    assert samples.std() / sd > 2.0


# ---------------------------------------------------------- end to end through process_compsep

def test_process_compsep_runs_a_marginal_mcmc_group(monkeypatch):
    """One band, one component, COMM_SELF: the per-pixel amplitude group, then an MCMC group whose
    'auto' likelihood resolves to 'marginal', exercised through the real process_compsep wiring."""
    from commander4.compsep import processing
    from commander4.compsep.processing import CompSepState, process_compsep

    params = Bunch(compsep=Bunch(nside=2, double_precision=True))
    sync = Bunch(enabled=True, component_class="Synchrotron", params=Bunch(
        lmax=5, polarization="I", shortname="sync", spatially_varying_MM=False,
        Cl_prior_amplitude=None, beta=-3.0, nu_ref=30.0, sample_spectral_index=True,
        spectral_index_proposal_sigma=0.05))
    object.__setattr__(sync, "_name", "Synchrotron")
    comp_list = CompList.init_from_params(Bunch({"Synchrotron": sync}), params)

    nside, npix = 2, 48
    rng = np.random.default_rng(4)
    detector_data = DetectorMap(map_sky=rng.uniform(5, 50, (1, npix)),
                                map_inv_var=np.ones((1, npix)), nu=30.0, fwhm=0.0, nside=nside)
    amp = {"pp": PerPixelSamplingGroupConfig(name="pp")}
    mcmc = {"beta": MCMCSamplingGroupConfig(name="beta", update_amplitude_groups=("pp",),
                                            numstep=20)}
    likelihoods = _resolve_mcmc_likelihoods(mcmc, amp, "per_pixel", comp_list, report=False)
    assert likelihoods == {"beta": "marginal"}
    state = CompSepState(
        band_name="Band30", band_identifier="Band30_I", target_pol="I",
        amplitude_method="per_pixel", amplitude_groups=amp, mcmc_groups=mcmc,
        double_precision=True, num_threads=1, chisq_masks={}, nside_chisq=1,
        include_chisq_map=False, include_residual_maps=False, mcmc_likelihoods=likelihoods)
    mpi_info = Bunch(compsep=Bunch(comm=MPI.COMM_SELF, subcomm=MPI.COMM_SELF, rank=0, master=0,
                                   size=1))
    built = []
    real_make = processing._make_marginal_loglike

    def spy(*args, **kwargs):
        built.append(real_make(*args, **kwargs))
        return built[-1]

    monkeypatch.setattr(processing, "_make_marginal_loglike", spy)
    monkeypatch.setattr(processing, "write_compsep_chain_to_file", lambda *args: None)
    process_compsep(mpi_info, state, detector_data, iter=1, chain=1, params=Bunch(),
                    comp_list=comp_list, rng=np.random.default_rng(0))
    assert len(built) == 1
    assert np.isfinite(built[0]())
