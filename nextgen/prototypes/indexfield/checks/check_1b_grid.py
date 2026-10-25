"""E012 v2 check 1b: exact grids with the real SED against NUTS, N_beta = 1 and N_beta = 2.

Target with the hyperparameters fixed: for N_beta index pixels with a fixed Gaussian prior,

    log post(beta) = log L_mode(beta) + log N(beta; mu_vec, C),      mode in {gauss, w0, w1},

with log L the collapsed likelihood of model.collapsed_loglike (dropped constant -1/2 d^T N^-1 d
in every mode, identical in reference.py). The reference is the exact grid of reference.grid_1d /
grid_2d (plain-loop per-pixel likelihood, no shared code with model.py), fine enough that the
quadrature error on the moments is negligible (801 points over +-9 posterior sd per axis). NUTS
runs on a one- or two-dimensional target whose log density is model.collapsed_loglike plus the
same Gaussian prior (the field classes need a HEALPix basis and do not apply to N_beta = 1 or 2;
checks 1c and 1d cover them).

Statistics (per mode and N_beta): (i) pointwise agreement of the prototype log density with the
grid log density at 25 grid points, spread (max - min) of the difference relative to |log L|;
(ii) NUTS mean within 4 MC standard errors (sd / sqrt(ESS), Geyer ESS) of the grid mean;
(iii) sd ratio NUTS / grid within [0.95, 1.05]; (iv) zero likelihood errors in the sampling phase.
The margin of the check is reported as the grid-mean difference w1 - w0 in units of the posterior
sd (the log det term that 1b is meant to catch).

Setting: nside_d 16 (N = 3072), 12 bands, amplitudes a ~ N(m, T) as in the SBC generator, white
noise at 1 x (the SBC setting) and 5 x the E006 levels (where the log det term dominates the w1
target), true beta -3.0 (N_beta 1) and (-3.1, -2.8) (N_beta 2, northern/southern hemisphere by
RING index).
N_beta = 2 uses a correlated prior N((-3, -3), 0.3^2 [[1, 0.5], [0.5, 1]]) so the 2-D grid is
not a product. Writes checks/out/check_1b_grid.json.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from indexfield import model, reference, sampler, simulate      # noqa: E402

OUT = os.path.join(HERE, "out", "check_1b_grid.json")
NSIDE_D = 16
NOISE_SCALES = (1.0, 5.0)
N_WARMUP, N_SAMPLES = 500, 2500
Z_MAX, SD_LO, SD_HI = 4.0, 0.95, 1.05
POINTWISE_TOL = 1e-8          # relative spread of (prototype - grid) log density


def to_py(o):
    if isinstance(o, dict):
        return {k: to_py(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_py(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def make_data(n_beta, beta_true, rng, amp_mode="gauss", noise_scale=1.0):
    """PixelData at nside_d with n_beta regions (1: all pixels; 2: hemispheres by RING index) and
    the SBC amplitude generator: a ~ N(m, T), T^(1/2) = 10 max(|m|, 10) (or template scatter)."""
    N = 12 * NSIDE_D ** 2
    region = np.zeros(N, dtype=np.int64) if n_beta == 1 else (np.arange(N) >= N // 2).astype(int)
    m = simulate.templates(NSIDE_D)
    T_sqrt = 10.0 * np.maximum(np.abs(m), 10.0)
    if amp_mode == "gauss":
        a = m + T_sqrt * rng.standard_normal(m.shape)
    else:
        a = m * np.exp(0.1 * rng.standard_normal(m.shape))
        a[:, 0] = 70.0 * rng.standard_normal(N)
    nu = simulate.BANDS_GHZ
    rms = simulate.noise_rms(NSIDE_D, nu, noise_scale)
    beta_true = np.asarray(beta_true, dtype=float)
    M, _ = model.mixing(nu, beta_true[region])
    d = np.einsum("nbc,nc->nb", M, a) + rms * rng.standard_normal(rms.shape)
    return model.PixelData(nu=nu, d=d, inv_var=1.0 / rms ** 2, region=region, n_beta=n_beta,
                           beta_d=1.6, T_d=20.0, m=m, inv_T=1.0 / T_sqrt ** 2)


class FixedPriorPosterior:
    """x = beta (N_beta,); log pi = collapsed_loglike(beta) + log N(beta; mu_vec, C)."""

    def __init__(self, data, mode, mu_vec, cov):
        self.data, self.mode = data, mode
        self.mu_vec = np.asarray(mu_vec, dtype=float)
        self.cov = np.atleast_2d(np.asarray(cov, dtype=float))
        self.prec = np.linalg.inv(self.cov)
        self.dim = self.mu_vec.size
        self.n_loglike_errors = 0
        self.log_norm = -0.5 * np.linalg.slogdet(self.cov)[1] - 0.5 * self.dim * np.log(2 * np.pi)

    def logp_and_grad(self, x):
        try:
            logL, g = model.collapsed_loglike(x, self.data, self.mode)
        except (ValueError, np.linalg.LinAlgError):
            self.n_loglike_errors += 1
            return -np.inf, np.zeros(self.dim)
        if not np.isfinite(logL):
            self.n_loglike_errors += 1
            return -np.inf, np.zeros(self.dim)
        r = x - self.mu_vec
        pr = self.prec @ r
        return float(logL - 0.5 * r @ pr + self.log_norm), g - pr

    def random_init(self, rng):
        return self.mu_vec + 0.2 * rng.standard_normal(self.dim)


def grid_moments(logpost, axes):
    """Mean and sd per axis and the grid cell weights of exp(logpost) on a product grid."""
    w = np.exp(logpost - logpost.max())
    w /= w.sum()
    means, sds = [], []
    for k, ax in enumerate(axes):
        wk = w.sum(axis=tuple(i for i in range(len(axes)) if i != k))
        mk = float(np.sum(wk * ax))
        means.append(mk)
        sds.append(float(np.sqrt(np.sum(wk * (ax - mk) ** 2))))
    return np.array(means), np.array(sds), w


def locate(data, mode, theta_fixed, n_beta):
    """Coarse-to-fine localization: returns per-axis fine grids over mean +- 9 sd (801 points).

    The coarse grid spans [-8, 4] (spacing 0.01) because the w1 target at low S/N can sit far
    from the truth (the -1/2 N log det A term rewards near-degenerate mixing matrices); the
    posterior sd is below the coarse spacing, so a second pass over +- 0.05 (1001 points) around
    the coarse mean fixes the location and width before the final grid."""
    coarse = np.linspace(-8.0, 4.0, 1201)
    if n_beta == 1:
        lp = reference.grid_1d(data, mode, coarse, theta_fixed=theta_fixed)
        m, s, _ = grid_moments(lp, [coarse])
    else:
        lp = reference.grid_2d(data, mode, coarse, coarse, theta_fixed=theta_fixed)
        m, s, _ = grid_moments(lp, [coarse, coarse])
    axes = [np.linspace(m[k] - 0.05, m[k] + 0.05, 1001) for k in range(n_beta)]
    if n_beta == 1:
        lp = reference.grid_1d(data, mode, axes[0], theta_fixed=theta_fixed)
        m, s, _ = grid_moments(lp, axes)
    else:
        lp = reference.grid_2d(data, mode, axes[0], axes[1], theta_fixed=theta_fixed)
        m, s, _ = grid_moments(lp, axes)
    return [np.linspace(m[k] - 9 * s[k], m[k] + 9 * s[k], 801) for k in range(n_beta)]


def run_case(n_beta, data, mode, mu_vec, cov, rng, noise_scale):
    theta_fixed = (float(mu_vec[0]), float(np.sqrt(cov[0, 0]))) if n_beta == 1 \
        else (np.asarray(mu_vec), np.asarray(cov))
    t0 = time.time()
    axes = locate(data, mode, theta_fixed, n_beta)
    if n_beta == 1:
        lp_grid = reference.grid_1d(data, mode, axes[0], theta_fixed=theta_fixed)
    else:
        lp_grid = reference.grid_2d(data, mode, axes[0], axes[1], theta_fixed=theta_fixed)
    g_mean, g_sd, w = grid_moments(lp_grid, axes)
    edge = float(max(w[[0, -1]].sum(), w[:, [0, -1]].sum() if n_beta == 2 else 0.0))
    t_grid = time.time() - t0

    post = FixedPriorPosterior(data, mode, mu_vec, cov)
    # (i) pointwise: prototype log density minus grid log density at 25 points.
    idx = np.linspace(0, 800, 25).astype(int)
    diffs, scale = [], 0.0
    for i in idx:
        if n_beta == 1:
            x = np.array([axes[0][i]])
            lg = lp_grid[i]
        else:
            x = np.array([axes[0][i], axes[1][800 - i]])
            lg = lp_grid[i, 800 - i]
        lp = post.logp_and_grad(x)[0]
        diffs.append(lp - lg)
        scale = max(scale, abs(lp))
    spread = float(np.max(diffs) - np.min(diffs)) / scale

    t0 = time.time()
    res = sampler.run_nuts(post, rng, N_WARMUP, N_SAMPLES)
    t_nuts = time.time() - t0
    s = res.samples
    n_mean = s.mean(axis=0)
    n_sd = s.std(axis=0, ddof=1)
    ess = res.ess
    z = (n_mean - g_mean) / (n_sd / np.sqrt(ess))
    sd_ratio = n_sd / g_sd
    ok_point = spread < POINTWISE_TOL
    ok_mean = bool(np.all(np.abs(z) < Z_MAX))
    ok_sd = bool(np.all((sd_ratio > SD_LO) & (sd_ratio < SD_HI)))
    ok_err = res.n_loglike_errors_sampling == 0
    ok = ok_point and ok_mean and ok_sd and ok_err
    gm = np.array2string(g_mean, precision=5)
    print(f"  noise x{noise_scale:g} N_beta={n_beta} {mode:5s}: grid mean {gm} sd "
          f"{np.array2string(g_sd, precision=2)} edge mass {edge:.1e} [{t_grid:.0f} s]; "
          f"pointwise spread {spread:.1e} (tol {POINTWISE_TOL:g}); NUTS n={N_SAMPLES} ESS "
          f"{np.array2string(ess, precision=0)} mean z {np.array2string(z, precision=2)} "
          f"sd ratio {np.array2string(sd_ratio, precision=3)} div={res.n_divergent} "
          f"errors={res.n_loglike_errors_sampling} [{t_nuts:.0f} s]  {'PASS' if ok else 'FAIL'}"
          + ("" if ok else f"  (point {ok_point} mean {ok_mean} sd {ok_sd} err {ok_err})"))
    return dict(n_beta=n_beta, mode=mode, noise_scale=noise_scale, grid_mean=g_mean,
                grid_sd=g_sd, grid_edge_mass=edge,
                grid_points=801, pointwise_spread_rel=spread, pointwise_n=len(idx),
                nuts_mean=n_mean, nuts_sd=n_sd, ess=ess, z_mean=z, sd_ratio=sd_ratio,
                n_divergent=int(res.n_divergent), step_size=res.step_size, ebfmi=res.ebfmi,
                n_loglike_errors_sampling=int(res.n_loglike_errors_sampling),
                n_warmup=N_WARMUP, n_samples=N_SAMPLES, seconds_grid=t_grid,
                seconds_nuts=t_nuts, passed=ok)


def main():
    t_all = time.time()
    rng = np.random.default_rng(1202)
    print(f"check 1b: nside_d={NSIDE_D}, 12 bands, gauss amplitudes, noise scales {NOISE_SCALES}; "
          f"NUTS {N_WARMUP}+{N_SAMPLES}; thresholds |z| < {Z_MAX}, sd ratio in "
          f"[{SD_LO}, {SD_HI}], pointwise spread < {POINTWISE_TOL:g}")
    rows = []
    cases = [(1, np.array([-3.0]), np.array([-3.0]), np.array([[0.25]])),
             (2, np.array([-3.1, -2.8]), np.array([-3.0, -3.0]),
              0.09 * np.array([[1.0, 0.5], [0.5, 1.0]]))]
    for noise_scale in NOISE_SCALES:
        for n_beta, beta_true, mu_vec, cov in cases:
            data = make_data(n_beta, beta_true, rng, noise_scale=noise_scale)
            print(f"noise x{noise_scale:g} N_beta={n_beta}: truth {beta_true}, prior mean "
                  f"{mu_vec}, "
                  f"prior sd {np.sqrt(np.diag(cov))}, pixels per region {np.bincount(data.region)}")
            means = {}
            for mode in ("gauss", "w0", "w1"):
                row = run_case(n_beta, data, mode, mu_vec, cov, rng, noise_scale)
                row["beta_true"] = beta_true
                rows.append(row)
                means[mode] = (row["grid_mean"], row["grid_sd"])
            margin = (means["w1"][0] - means["w0"][0]) / means["w0"][1]
            margin_g = (means["gauss"][0] - means["w0"][0]) / means["w0"][1]
            print(f"  margin: grid mean (w1 - w0) / sd_w0 = "
                  f"{np.array2string(margin, precision=2)}; (gauss - w0) / sd_w0 = "
                  f"{np.array2string(margin_g, precision=2)}")
            for r in rows:
                if r["n_beta"] == n_beta and r["noise_scale"] == noise_scale:
                    r["margin_w1_minus_w0_in_sd_w0"] = margin
    passed = all(r["passed"] for r in rows)
    worst_z = max(float(np.max(np.abs(r["z_mean"]))) for r in rows)
    sd_range = (min(float(np.min(r["sd_ratio"])) for r in rows),
                max(float(np.max(r["sd_ratio"])) for r in rows))
    min_ess = min(float(np.min(r["ess"])) for r in rows)
    summary = dict(check="1b_grid", setting=dict(nside_d=NSIDE_D, noise_scales=NOISE_SCALES,
                                                 amp_mode="gauss", n_warmup=N_WARMUP,
                                                 n_samples=N_SAMPLES, modes=["gauss", "w0", "w1"],
                                                 n_beta=[1, 2]),
                   statistic=dict(worst_abs_z=worst_z, sd_ratio_range=sd_range, min_ess=min_ess,
                                  worst_pointwise_spread=max(r["pointwise_spread_rel"]
                                                             for r in rows)),
                   threshold=f"|z| < {Z_MAX}, sd ratio in [{SD_LO}, {SD_HI}], pointwise < "
                             f"{POINTWISE_TOL:g}, zero sampling-phase likelihood errors",
                   passed=passed, rows=rows, seconds=time.time() - t_all)
    with open(OUT, "w") as f:
        json.dump(to_py(summary), f, indent=1)
    print(f"check 1b grid: worst |z| = {worst_z:.2f} (threshold {Z_MAX}), sd ratio range "
          f"[{sd_range[0]:.3f}, {sd_range[1]:.3f}] (threshold [{SD_LO}, {SD_HI}]), min ESS "
          f"{min_ess:.0f}  {'PASS' if passed else 'FAIL'}  [{time.time() - t_all:.0f} s, {OUT}]")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
