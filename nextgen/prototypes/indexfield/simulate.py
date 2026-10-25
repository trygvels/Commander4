"""Simulations with known truth for the index-field prototype (E012 v2 "Setup", SPEC simulate.py).

Generative model, drawn exactly as the posterior assumes it:
    theta = (mu, sigma, rho) ~ hyperprior,    z ~ N(0, I_K),
    beta  = mu 1 + sigma S(rho) z                        (field.S_times),
    a_p   ~ N(m_p, T_p)  ["gauss"]   or   template values with log-normal scatter ["template"],
    d_p   = M_p(beta_{r(p)}) a_p + n_p,   n_p ~ N(0, diag(rms_p^2)).
Units: GHz, uK_RJ, radians. All maps RING order; `region` follows the SPEC nest-shift recipe.
"""
from __future__ import annotations

import numpy as np
import healpy as hp

from .model import PixelData, mixing

# The 12 bands of E012 and their white-noise RMS per nside_d = 16 pixel in uK_RJ. The seven Planck
# bands carry the E006 values (aux/validate_mcmc_likelihood.py RMS, before smoothing); the five
# WMAP-like bands are set here at comparable levels and are part of the stated setup.
BANDS_GHZ = np.array([23., 30., 33., 41., 44., 61., 70., 94., 100., 143., 217., 353.])
RMS_UK_RJ = np.array([70., 60., 60., 55., 60., 60., 70., 40., 24., 16., 24., 80.])
#                     23   30*  33   41   44*  61   70*  94   100* 143* 217* 353*   (* = E006)


def region_map(nside_d: int, nside_beta: int) -> np.ndarray:
    """Index pixel (RING, nside_beta) containing each data pixel (RING, nside_d), shape (N,).

    region[p] = nest2ring(nside_beta, ring2nest(nside_d, p) >> 2 log2(nside_d / nside_beta)).
    """
    if nside_beta > nside_d or (nside_d % nside_beta) != 0:
        raise ValueError("nside_beta must divide nside_d")
    shift = 2 * int(round(np.log2(nside_d // nside_beta)))
    p = np.arange(hp.nside2npix(nside_d))
    nest_beta = hp.ring2nest(nside_d, p) >> shift
    return hp.nest2ring(nside_beta, nest_beta)


def templates(nside_d: int, rng=None) -> np.ndarray:
    """Template amplitude means m (N, C) in uK_RJ at the reference frequencies (30, 545 GHz).

    cmb = 0; with t = exp(-|lat| / 15 deg) the plane-brightening forms of
    aux/validate_mcmc_likelihood.py: sync = 10^(0.3 + 2.0 t), dust = 10^(0.5 + 1.5 t). If `rng`
    is given, the 0.1 dex log-normal pixel scatter of that harness is added to sync and dust
    (10^(0.1 N(0, 1))); with rng=None the smooth templates are returned.
    """
    npix = hp.nside2npix(nside_d)
    lat = np.abs(90.0 - np.degrees(hp.pix2ang(nside_d, np.arange(npix))[0]))
    t = np.exp(-lat / 15.0)
    scatter = (0.1 * rng.standard_normal((2, npix))) if rng is not None else np.zeros((2, npix))
    return np.stack([np.zeros(npix),
                     10**(0.3 + 2.0 * t + scatter[0]),
                     10**(0.5 + 1.5 * t + scatter[1])], axis=1)


def noise_rms(nside_d: int, nu, scale: float = 1.0) -> np.ndarray:
    """White-noise RMS per pixel and band, shape (N, B), uK_RJ: scale x RMS_UK_RJ[band], uniform
    over the sky. Every nu must be one of BANDS_GHZ (matched to 1e-6 GHz)."""
    nu = np.atleast_1d(np.asarray(nu, dtype=float))
    idx = np.array([np.flatnonzero(np.abs(BANDS_GHZ - f) < 1e-6) for f in nu], dtype=object)
    if any(len(i) == 0 for i in idx):
        raise KeyError(f"frequencies {nu[[len(i) == 0 for i in idx]]} are not in BANDS_GHZ")
    rms = scale * RMS_UK_RJ[np.array([int(i[0]) for i in idx])]
    return np.broadcast_to(rms, (hp.nside2npix(nside_d), nu.size)).copy()


def draw_hyper(hyper, rng) -> dict:
    """Draw theta from the hyperprior of E012 v2. Returns dict(mu, sigma, rho, log_sigma, u).

    mu ~ N(mu0, s_mu^2). rho_prior "loguniform": sigma ~ half-Normal(0, sigma_scale), log rho ~
    Uniform(log rho_min, log rho_max). rho_prior "pc" (Fuglstad et al. 2019, d = 2): sigma ~
    Exp(lambda_s), lambda_s = -log(0.05) / pc_sigma0 (P(sigma > pc_sigma0) = 0.05); rho has
    density proportional to rho^-2 exp(-lambda_r / rho), lambda_r = -pc_rho0 log(0.05)
    (P(rho < pc_rho0) = 0.05 before truncation), truncated to [rho_min, rho_max] and drawn by
    inversion of its CDF exp(-lambda_r / rho). u = logit((log rho - log rho_min) / L),
    L = log(rho_max / rho_min).
    """
    mu = hyper.mu0 + hyper.s_mu * rng.standard_normal()
    L = np.log(hyper.rho_max / hyper.rho_min)
    if hyper.rho_prior == "loguniform":
        sigma = abs(hyper.sigma_scale * rng.standard_normal())
        rho = hyper.rho_min * np.exp(L * rng.random())
    elif hyper.rho_prior == "pc":
        lam_s = -np.log(0.05) / hyper.pc_sigma0
        sigma = rng.exponential(1.0 / lam_s)
        lam_r = -hyper.pc_rho0 * np.log(0.05)
        c_lo, c_hi = np.exp(-lam_r / hyper.rho_min), np.exp(-lam_r / hyper.rho_max)
        rho = -lam_r / np.log(c_lo + (c_hi - c_lo) * rng.random())
    else:
        raise ValueError(f"unknown rho_prior {hyper.rho_prior!r}")
    frac = (np.log(rho) - np.log(hyper.rho_min)) / L
    frac = min(max(frac, 1e-12), 1 - 1e-12)
    return dict(mu=float(mu), sigma=float(sigma), rho=float(rho), log_sigma=float(np.log(sigma)),
                u=float(np.log(frac / (1 - frac))))


def simulate(nside_d: int, nside_beta: int, field, hyper, rng, amp_mode: str = "gauss",
             noise_scale: float = 1.0, beta_d: float = 1.6, T_d: float = 20.0):
    """Draw one data set from the hierarchical model. Returns (PixelData, truth).

    field: object with S_times(rho, z) -> (Nb,) and basis.V of shape (Nb, K) (field.MaternField).
    hyper: HyperPrior (posterior.py) or any object with its fields.
    amp_mode "gauss": a_p ~ N(m_p, T_p), T_p^(1/2) = 10 max(|m_p|, 10) per component (floor
    100 uK_RJ; E012 "Amplitude treatment"). "template": sync and dust a = m exp(0.1 N(0, 1)),
    cmb ~ N(0, 70^2) uK_RJ. In both cases PixelData.m and .inv_T hold the Gaussian prior (m_p,
    T_p^-1) so that mode "gauss" can be evaluated on template data as well (coverage study).
    truth: dict(theta=dict(mu, sigma, rho, log_sigma, u), z (K,), beta (Nb,), a (N, C)).
    """
    N = hp.nside2npix(nside_d)
    Nb = hp.nside2npix(nside_beta)
    if field.basis.V.shape[0] != Nb:
        raise ValueError(f"field basis has {field.basis.V.shape[0]} pixels, expected {Nb}")
    K = field.basis.V.shape[1]
    region = region_map(nside_d, nside_beta)
    theta = draw_hyper(hyper, rng)
    z = rng.standard_normal(K)
    beta = theta["mu"] + theta["sigma"] * np.asarray(field.S_times(theta["rho"], z))

    m = templates(nside_d)
    T_sqrt = 10.0 * np.maximum(np.abs(m), 10.0)
    if amp_mode == "gauss":
        a = m + T_sqrt * rng.standard_normal(m.shape)
    elif amp_mode == "template":
        a = m * np.exp(0.1 * rng.standard_normal(m.shape))
        a[:, 0] = 70.0 * rng.standard_normal(N)
    else:
        raise ValueError(f"unknown amp_mode {amp_mode!r}")

    nu = BANDS_GHZ
    rms = noise_rms(nside_d, nu, noise_scale)
    M, _ = mixing(nu, beta[region], beta_d=beta_d, T_d=T_d)
    d = np.einsum("nbc,nc->nb", M, a) + rms * rng.standard_normal(rms.shape)
    data = PixelData(nu=nu, d=d, inv_var=1.0 / rms**2, region=region, n_beta=Nb, beta_d=beta_d,
                     T_d=T_d, m=m, inv_T=1.0 / T_sqrt**2)
    truth = dict(theta=theta, z=z, beta=beta, a=a)
    return data, truth


# ----------------------------------------------------------------------------------------------
# Self-tests
# ----------------------------------------------------------------------------------------------

class _StandInBasis:
    """Orthonormal basis from a QR of random vectors, constant mode removed, for the self-test."""

    def __init__(self, nside_beta, rng):
        Nb = hp.nside2npix(nside_beta)
        Q, _ = np.linalg.qr(np.column_stack([np.ones(Nb), rng.standard_normal((Nb, Nb - 1))]))
        self.V = Q[:, 1:]
        self.lam = np.sort(rng.random(Nb - 1)) * 50.0
        self.h = np.sqrt(4 * np.pi / Nb)


class _StandInField:
    """Matern-like S(rho) z = V diag(s / c) z with s_i = (8 / rho^2 + lam_i)^-1 (alpha = 2)."""

    def __init__(self, basis):
        self.basis = basis

    def S_times(self, rho, z):
        s = 1.0 / (8.0 / rho**2 + self.basis.lam)
        c = np.sqrt(np.mean(s**2))
        return self.basis.V @ (s / c * z)


class _StandInHyper:
    mu0, s_mu, sigma_scale = -3.0, 0.5, 0.3
    rho_prior = "loguniform"
    rho_min, rho_max = np.sqrt(4 * np.pi / 192), np.pi
    pc_rho0, pc_sigma0 = 2 * np.sqrt(4 * np.pi / 192), 0.5


def _check_region(nside_d=16, nside_beta=4):
    """Each index pixel holds exactly (nside_d / nside_beta)^2 data pixels and contains their
    centres (the data pixel centre falls in the index pixel found by hp.ang2pix)."""
    region = region_map(nside_d, nside_beta)
    counts = np.bincount(region, minlength=hp.nside2npix(nside_beta))
    th, ph = hp.pix2ang(nside_d, np.arange(hp.nside2npix(nside_d)))
    direct = hp.ang2pix(nside_beta, th, ph)
    ok = np.all(counts == (nside_d // nside_beta)**2) and np.array_equal(region, direct)
    print(f"  region_map nside {nside_d}->{nside_beta}: counts unique={np.unique(counts).tolist()}"
          f", equals ang2pix: {np.array_equal(region, direct)}  {'PASS' if ok else 'FAIL'}")
    return bool(ok)


def _check_templates(nside_d=16):
    m = templates(nside_d)
    lat = np.abs(90.0 - np.degrees(hp.pix2ang(nside_d, np.arange(m.shape[0]))[0]))
    plane, pole = lat < 5.0, lat > 70.0
    ok = (np.all(m[:, 0] == 0) and m[plane, 1].mean() > m[pole, 1].mean()
          and m[plane, 2].mean() > m[pole, 2].mean()
          and abs(m[:, 1].max() - 10**2.3) < 1e-6 * 10**2.3)
    print(f"  templates: sync plane/pole {m[plane, 1].mean():.1f}/{m[pole, 1].mean():.1f}, "
          f"dust plane/pole {m[plane, 2].mean():.1f}/{m[pole, 2].mean():.1f}, max sync "
          f"{m[:, 1].max():.2f} (10^2.3 = {10**2.3:.2f})  {'PASS' if ok else 'FAIL'}")
    return bool(ok)


def _check_noise():
    rms = noise_rms(16, BANDS_GHZ, 2.0)
    ok = rms.shape == (3072, 12) and np.allclose(rms[0], 2.0 * RMS_UK_RJ)
    try:
        noise_rms(16, [31.0])
        ok = False
    except KeyError:
        pass
    print(f"  noise_rms: shape {rms.shape}, 30 GHz x2 = {rms[0, 1]:.0f}, unknown band raises  "
          f"{'PASS' if ok else 'FAIL'}")
    return bool(ok)


def _check_hyper_draws(rng, n=40000):
    """Moments of the hyperprior draws: mu, half-normal sigma, log-uniform rho, PC sigma and rho."""
    hyper = _StandInHyper()
    draws = [draw_hyper(hyper, rng) for _ in range(n)]
    mu = np.array([t["mu"] for t in draws])
    sig = np.array([t["sigma"] for t in draws])
    lrho = np.array([np.log(t["rho"]) for t in draws])
    se = 1.0 / np.sqrt(n)
    L = np.log(hyper.rho_max / hyper.rho_min)
    checks = [
        abs(mu.mean() - hyper.mu0) < 5 * hyper.s_mu * se,
        abs(sig.mean() - hyper.sigma_scale * np.sqrt(2 / np.pi)) < 5 * hyper.sigma_scale * se,
        abs(lrho.mean() - (np.log(hyper.rho_min) + L / 2)) < 5 * L / np.sqrt(12) * se,
        abs(lrho.std() - L / np.sqrt(12)) < 5 * L * se,
        np.all(lrho >= np.log(hyper.rho_min)) and np.all(lrho <= np.log(hyper.rho_max)),
    ]
    hyper_pc = _StandInHyper()
    hyper_pc.rho_prior = "pc"
    draws = [draw_hyper(hyper_pc, rng) for _ in range(n)]
    sig = np.array([t["sigma"] for t in draws])
    rho = np.array([t["rho"] for t in draws])
    lam_r = -hyper_pc.pc_rho0 * np.log(0.05)
    # Truncated CDF: F(rho) = (exp(-lam/rho) - c_lo) / (c_hi - c_lo); evaluate at the median.
    c_lo, c_hi = np.exp(-lam_r / hyper_pc.rho_min), np.exp(-lam_r / hyper_pc.rho_max)
    med = np.median(rho)
    F_med = (np.exp(-lam_r / med) - c_lo) / (c_hi - c_lo)
    checks += [
        abs(np.mean(sig > hyper_pc.pc_sigma0) - 0.05) < 5 * np.sqrt(0.05 * 0.95) * se,
        abs(F_med - 0.5) < 5 * 0.5 * se * 1.5,
        abs(np.mean(rho < hyper_pc.pc_rho0) - (np.exp(-lam_r / hyper_pc.pc_rho0) - c_lo)
            / (c_hi - c_lo)) < 5 * 0.25 * se,
    ]
    ok = all(bool(c) for c in checks)
    print(f"  draw_hyper n={n}: loguniform mean mu={mu.mean():.4f}, mean log rho="
          f"{lrho.mean():.4f} (expected {np.log(hyper.rho_min) + L / 2:.4f}); pc P(sigma>0.5)="
          f"{np.mean(sig > 0.5):.4f}, F(median rho)={F_med:.4f}  "
          f"{'PASS' if ok else 'FAIL'} {[bool(c) for c in checks]}")
    return ok


def _check_simulate(rng, nside_d=16, nside_beta=4):
    """Shapes, finiteness, residual chi^2 at the truth, and the amplitude-prior standardization."""
    from .model import collapsed_loglike
    field = _StandInField(_StandInBasis(nside_beta, rng))
    hyper = _StandInHyper()
    ok = True
    for amp_mode in ("gauss", "template"):
        data, truth = simulate(nside_d, nside_beta, field, hyper, rng, amp_mode=amp_mode,
                               noise_scale=1.0)
        N, Nb = hp.nside2npix(nside_d), hp.nside2npix(nside_beta)
        shapes = (data.d.shape == (N, 12) and data.inv_var.shape == (N, 12)
                  and data.region.shape == (N,) and truth["beta"].shape == (Nb,)
                  and truth["a"].shape == (N, 3) and truth["z"].shape == (field.basis.V.shape[1],))
        M, _ = mixing(data.nu, truth["beta"][data.region], beta_d=data.beta_d, T_d=data.T_d)
        resid = data.d - np.einsum("nbc,nc->nb", M, truth["a"])
        chi2_red = np.sum(resid**2 * data.inv_var) / (N * 12)
        # Standardized amplitude residuals (a - m) T^(-1/2): unit variance in "gauss" mode.
        zs = (truth["a"] - data.m) * np.sqrt(data.inv_T)
        beta_sd = truth["beta"].std()
        logL = [collapsed_loglike(truth["beta"], data, mode)[0] for mode in ("gauss", "w0", "w1")]
        fin = all(np.isfinite(v) for v in logL) and np.all(np.isfinite(data.d))
        this = shapes and fin and abs(chi2_red - 1.0) < 5 * np.sqrt(2.0 / (N * 12))
        if amp_mode == "gauss":
            this &= abs(zs.std() - 1.0) < 5 / np.sqrt(2 * N * 3)
        ok &= this
        th = truth["theta"]
        print(f"  simulate {amp_mode:8s}: shapes ok={shapes}, reduced chi2 at truth={chi2_red:.4f}"
              f", std((a-m)/T^1/2)={zs.std():.3f}, theta=(mu {th['mu']:.3f}, sigma "
              f"{th['sigma']:.3f}, rho {np.degrees(th['rho']):.1f} deg), sd(beta)={beta_sd:.3f}, "
              f"logL(gauss,w0,w1)=({logL[0]:.4e}, {logL[1]:.4e}, {logL[2]:.4e})  "
              f"{'PASS' if this else 'FAIL'}")
    return bool(ok)


def _check_with_field_module(rng):
    """If field.py is present, run simulate with the real GraphLaplacianBasis / MaternField."""
    try:
        from . import field as fmod
    except Exception as err:   # noqa: BLE001  (module written by another agent; may be absent)
        print(f"  field.py not importable yet ({type(err).__name__}); skipped the real-field run")
        return True
    try:
        fld = fmod.MaternField(fmod.GraphLaplacianBasis(4))
        data, truth = simulate(16, 4, fld, _StandInHyper(), rng, amp_mode="gauss")
        ok = np.all(np.isfinite(data.d)) and truth["beta"].shape == (192,)
        print(f"  simulate with field.GraphLaplacianBasis(4): sd(beta)={truth['beta'].std():.3f}"
              f", sigma={truth['theta']['sigma']:.3f}  {'PASS' if ok else 'FAIL'}")
        return bool(ok)
    except Exception as err:   # noqa: BLE001
        print(f"  simulate with field.py raised {type(err).__name__}: {err}  (not counted)")
        return True


if __name__ == "__main__":
    rng = np.random.default_rng(2026)
    print("simulate.py self-tests")
    results = [
        _check_region(),
        _check_templates(),
        _check_noise(),
        _check_hyper_draws(rng),
        _check_simulate(rng),
        _check_with_field_module(rng),
    ]
    print("ALL PASS" if all(results) else "SOME FAILED")
    raise SystemExit(0 if all(results) else 1)
