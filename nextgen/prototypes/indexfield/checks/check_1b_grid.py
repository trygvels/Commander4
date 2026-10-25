"""E012 v2 check 1b: exact grids with the real SED against NUTS, N_beta = 1 and N_beta = 2.

Target with the hyperparameters fixed: for N_beta index pixels with a fixed Gaussian prior,

    log post(beta) = log L_mode(beta) + log N(beta; mu_vec, C),      mode in {gauss, w0, w1},

with log L the collapsed likelihood of model.collapsed_loglike (dropped constant -1/2 d^T N^-1 d
in every mode, identical in reference.py). The reference is the exact grid of reference.grid_1d /
grid_2d (plain-loop per-pixel likelihood, no shared code with model.py). NUTS runs on a one- or
two-dimensional target whose log density is model.collapsed_loglike plus the same Gaussian prior
(the field classes need a HEALPix basis and do not apply to N_beta = 1 or 2; checks 1c and 1d
cover them).

Grid construction (coarse to fine, basin by basin). A coarse grid over [-8, 4] with spacing 0.01
locates every strict local maximum of the log posterior within BASIN_DROP log units of the top
(the w1 target at low S/N is multimodal: the -1/2 sum_p log det A_p term rewards near-degenerate
CMB/synchrotron columns at beta ~ -0.25 and keeps a secondary basin near the truth). The
posterior sd (0.002 to 0.06 here) is below 5 coarse spacings in every case, so the coarse moments
are quantized and each basin is refined locally with spacing 1e-4 over +-0.05 around its maximum;
from the refined moments two further passes place an 801-point window over mean +- 9 sd per axis
(the final grid, with quadrature error on the moments far below the MC error). Each basin's mass
is its log-sum-exp times the cell volume; the reference mean and sd are those of the mixture of
basins weighted by mass (equal to the dominant basin's when the others carry negligible mass),
and the mass of every basin is reported.

NUTS: N_CHAINS = 4 chains from dispersed starts, the 5, 35, 65 and 95 percent quantiles of the
grid's 1-D marginals (mass-weighted over all basins), run in parallel processes. Each chain's
mean is assigned to the nearest basin maximum. When all chains sit in one basin and agree
pairwise (|mean_i - mean_j| < 5 combined SE on every axis) they are pooled (ESS = sum of the
Geyer ESS of the chains) and compared with the grid. When the chains do not agree the case is
reported as MULTIMODAL (chains split between basins: the fraction of draws per basin is compared
with the grid mass per basin, and the chains in the dominant basin are compared with that
basin's grid moments) or FAIL (chains in one basin that disagree beyond 5 SE: a mixing problem),
with the numbers, instead of a bare FAIL.

Statistics (per mode and N_beta): (i) pointwise agreement of the prototype log density with the
grid log density at 25 points of the dominant basin's grid, spread (max - min) of the difference
relative to |log L|; (ii) pooled NUTS mean within 4 MC standard errors (sd / sqrt(ESS)) of the
grid mean; (iii) sd ratio NUTS / grid within [0.95, 1.05]; (iv) zero likelihood errors in the
sampling phase. The margin of the check is reported as the grid-mean difference w1 - w0 in units
of the posterior sd (the log det term that 1b is meant to catch).

Setting: nside_d 16 (N = 3072), 12 bands, amplitudes a ~ N(m, T) as in the SBC generator
("wide" width, the first-build setting), white noise at 1 x (the SBC setting) and 5 x the E006
levels (where the log det term dominates the w1 target), true beta -3.0 (N_beta 1) and
(-3.1, -2.8) (N_beta 2, northern/southern hemisphere by RING index).
N_beta = 2 uses a correlated prior N((-3, -3), 0.3^2 [[1, 0.5], [0.5, 1]]) so the 2-D grid is
not a product. Writes checks/out/check_1b_grid.json. Exit status 0 PASS, 1 FAIL, 2 MULTIMODAL.
"""
from __future__ import annotations

import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import multiprocessing as mp

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np                                               # noqa: E402
from scipy.ndimage import maximum_filter                         # noqa: E402
from scipy.special import logsumexp                              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from indexfield import model, reference, sampler, simulate      # noqa: E402

OUT = os.path.join(HERE, "out", "check_1b_grid.json")
NSIDE_D = 16
NOISE_SCALES = (1.0, 5.0)
N_WARMUP, N_SAMPLES = 500, 2500
N_CHAINS = 4
START_QUANTILES = (0.05, 0.35, 0.65, 0.95)
Z_MAX, SD_LO, SD_HI = 4.0, 0.95, 1.05
Z_CHAIN = 5.0                 # pairwise chain agreement (combined SE) required for pooling
POINTWISE_TOL = 1e-8          # relative spread of (prototype - grid) log density
COARSE = np.linspace(-8.0, 4.0, 1201)      # spacing 0.01
FINE_SPACING, FINE_HALF = 1e-4, 0.05       # local refinement around every basin maximum
BASIN_DROP = 1e4              # local maxima more than this below the top are ignored
MAX_BASINS = 6
WINDOW_SD, WINDOW_POINTS = 9.0, 801


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


# ------------------------------------------------------------------------------------- grids
def grid_logpost(data, mode, theta_fixed, axes):
    if len(axes) == 1:
        return reference.grid_1d(data, mode, axes[0], theta_fixed=theta_fixed)
    return reference.grid_2d(data, mode, axes[0], axes[1], theta_fixed=theta_fixed)


def grid_moments(logpost, axes):
    """Mean and sd per axis and the normalized cell weights of exp(logpost) on a product grid."""
    w = np.exp(logpost - logpost.max())
    w /= w.sum()
    means, sds = [], []
    for k, ax in enumerate(axes):
        wk = w.sum(axis=tuple(i for i in range(len(axes)) if i != k))
        mk = float(np.sum(wk * ax))
        means.append(mk)
        sds.append(float(np.sqrt(np.sum(wk * (ax - mk) ** 2))))
    return np.array(means), np.array(sds), w


def log_mass(logpost, axes):
    """log of the integral of exp(logpost) over the grid window (cell volume x log-sum-exp)."""
    vol = float(np.prod([ax[1] - ax[0] for ax in axes]))
    return float(logsumexp(logpost) + np.log(vol))


def edge_mass(w):
    if w.ndim == 1:
        return float(w[[0, -1]].sum())
    return float(max(w[[0, -1], :].sum(), w[:, [0, -1]].sum()))


def find_maxima(lp, axes):
    """Coordinates of the strict local maxima of lp on the coarse grid within BASIN_DROP of the
    top (interior points only), sorted by decreasing lp, at most MAX_BASINS."""
    if lp.ndim == 1:
        inner = (lp[1:-1] > lp[:-2]) & (lp[1:-1] > lp[2:])
        idx = [(i + 1,) for i in np.flatnonzero(inner)]
    else:
        foot = np.ones((3, 3), dtype=bool)
        foot[1, 1] = False
        nb_max = maximum_filter(lp, footprint=foot, mode="constant", cval=-np.inf)
        mask = lp > nb_max
        mask[[0, -1], :] = False
        mask[:, [0, -1]] = False
        idx = [tuple(ij) for ij in np.argwhere(mask)]
    idx = [ij for ij in idx if lp[ij] > lp.max() - BASIN_DROP]
    idx.sort(key=lambda ij: -lp[ij])
    return [np.array([axes[k][ij[k]] for k in range(lp.ndim)]) for ij in idx[:MAX_BASINS]]


def refine_basin(data, mode, theta_fixed, centre, n_beta):
    """Local refinement at spacing FINE_SPACING over +-FINE_HALF around a basin maximum, then
    two windows of WINDOW_POINTS points over mean +- WINDOW_SD sd per axis. Returns the basin
    dict (axes, logpost, weights, mean, sd, log_mass, edge_mass, lp_max)."""
    n_fine = int(round(2 * FINE_HALF / FINE_SPACING)) + 1
    axes = [np.linspace(c - FINE_HALF, c + FINE_HALF, n_fine) for c in centre]
    lp = grid_logpost(data, mode, theta_fixed, axes)
    m, s, w = grid_moments(lp, axes)
    for _ in range(2):
        axes = [np.linspace(m[k] - WINDOW_SD * s[k], m[k] + WINDOW_SD * s[k], WINDOW_POINTS)
                for k in range(n_beta)]
        lp = grid_logpost(data, mode, theta_fixed, axes)
        m, s, w = grid_moments(lp, axes)
    return dict(centre_coarse=centre, axes=axes, logpost=lp, weights=w, mean=m, sd=s,
                log_mass=log_mass(lp, axes), edge_mass=edge_mass(w), lp_max=float(lp.max()))


def build_reference(data, mode, theta_fixed, n_beta):
    """Coarse grid, basin detection, per-basin refinement, mixture moments and marginal
    quantiles for the chain starts."""
    axes_c = [COARSE] * n_beta
    lp_c = grid_logpost(data, mode, theta_fixed, axes_c)
    centres = find_maxima(lp_c, axes_c)
    if not centres:                                   # maximum on the coarse boundary
        ij = np.unravel_index(int(np.argmax(lp_c)), lp_c.shape)
        centres = [np.array([axes_c[k][ij[k]] for k in range(n_beta)])]
    basins = []
    for c in centres:
        # skip maxima that fall inside an already refined (higher) basin's window
        if any(np.all((c >= b["axes_lo"]) & (c <= b["axes_hi"])) for b in basins):
            continue
        b = refine_basin(data, mode, theta_fixed, c, n_beta)
        b["axes_lo"] = np.array([ax[0] for ax in b["axes"]])
        b["axes_hi"] = np.array([ax[-1] for ax in b["axes"]])
        basins.append(b)
    lm = np.array([b["log_mass"] for b in basins])
    p = np.exp(lm - lm.max())
    p /= p.sum()
    order = np.argsort(-p)
    basins = [basins[i] for i in order]
    p = p[order]
    for b, pb in zip(basins, p):
        b["mass_frac"] = float(pb)
    means = np.array([b["mean"] for b in basins])
    sds = np.array([b["sd"] for b in basins])
    mean = (p[:, None] * means).sum(axis=0)
    sd = np.sqrt((p[:, None] * (sds ** 2 + (means - mean) ** 2)).sum(axis=0))
    # Mass-weighted 1-D marginals over all basins, for the chain starts.
    quant = np.empty((len(START_QUANTILES), n_beta))
    for k in range(n_beta):
        pts, wts = [], []
        for b, pb in zip(basins, p):
            wk = b["weights"].sum(axis=tuple(i for i in range(n_beta) if i != k))
            pts.append(b["axes"][k])
            wts.append(pb * wk)
        pts, wts = np.concatenate(pts), np.concatenate(wts)
        o = np.argsort(pts)
        cdf = np.cumsum(wts[o])
        cdf /= cdf[-1]
        quant[:, k] = np.interp(START_QUANTILES, cdf, pts[o])
    return dict(basins=basins, mean=mean, sd=sd, starts=quant, n_coarse=lp_c.size,
                coarse_max=float(lp_c.max()))


# ------------------------------------------------------------------------------------- chains
def _chain_worker(args):
    post, x0, seed, n_warmup, n_samples = args
    rng = np.random.default_rng(seed)
    res = sampler.run_nuts(post, rng, n_warmup, n_samples, x0=np.asarray(x0, dtype=float))
    return dict(samples=res.samples, ess=res.ess, n_divergent=int(res.n_divergent),
                n_loglike_errors_sampling=int(res.n_loglike_errors_sampling),
                step_size=res.step_size, ebfmi=res.ebfmi, accept_stat=res.accept_stat)


def run_chains(post, starts, seeds, n_warmup, n_samples, workers):
    jobs = [(post, x0, int(seed), n_warmup, n_samples) for x0, seed in zip(starts, seeds)]
    if workers <= 1:
        return [_chain_worker(j) for j in jobs]
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as ex:
        return list(ex.map(_chain_worker, jobs))


def nearest_basin(points, basins):
    centres = np.array([b["mean"] for b in basins])            # (n_basins, n_beta)
    d = np.linalg.norm(points[:, None, :] - centres[None, :, :], axis=2)
    return np.argmin(d, axis=1)


def chains_agree(chains):
    """Pairwise |mean_i - mean_j| < Z_CHAIN combined SE on every axis; returns (bool, worst z)."""
    worst = 0.0
    for i in range(len(chains)):
        for j in range(i + 1, len(chains)):
            ci, cj = chains[i], chains[j]
            se = np.sqrt(ci["sd"] ** 2 / ci["ess"] + cj["sd"] ** 2 / cj["ess"])
            worst = max(worst, float(np.max(np.abs(ci["mean"] - cj["mean"]) / se)))
    return worst < Z_CHAIN, worst


def pooled_stats(chains, ref_mean, ref_sd):
    s = np.concatenate([c["samples"] for c in chains])
    ess = np.sum([c["ess"] for c in chains], axis=0)
    mean, sd = s.mean(axis=0), s.std(axis=0, ddof=1)
    z = (mean - ref_mean) / (sd / np.sqrt(ess))
    return dict(mean=mean, sd=sd, ess=ess, z_mean=z, sd_ratio=sd / ref_sd, n=s.shape[0],
                n_chains=len(chains))


def fmt(a, prec=4):
    return np.array2string(np.asarray(a), precision=prec, suppress_small=False)


def run_case(n_beta, data, mode, mu_vec, cov, rng, noise_scale, workers):
    theta_fixed = (float(mu_vec[0]), float(np.sqrt(cov[0, 0]))) if n_beta == 1 \
        else (np.asarray(mu_vec), np.asarray(cov))
    t0 = time.time()
    ref = build_reference(data, mode, theta_fixed, n_beta)
    basins = ref["basins"]
    g_mean, g_sd = ref["mean"], ref["sd"]
    t_grid = time.time() - t0
    top = basins[0]
    post = FixedPriorPosterior(data, mode, mu_vec, cov)
    # (i) pointwise: prototype log density minus grid log density at 25 points of the top basin.
    axes, lp_grid = top["axes"], top["logpost"]
    last = WINDOW_POINTS - 1
    idx = np.linspace(0, last, 25).astype(int)
    diffs, scale = [], 0.0
    for i in idx:
        if n_beta == 1:
            x, lg = np.array([axes[0][i]]), lp_grid[i]
        else:
            x, lg = np.array([axes[0][i], axes[1][last - i]]), lp_grid[i, last - i]
        lp = post.logp_and_grad(x)[0]
        diffs.append(lp - lg)
        scale = max(scale, abs(lp))
    spread = float(np.max(diffs) - np.min(diffs)) / scale

    # NUTS: N_CHAINS dispersed starts in parallel.
    t0 = time.time()
    seeds = rng.integers(0, 2 ** 31 - 1, size=N_CHAINS)
    raw = run_chains(post, ref["starts"], seeds, N_WARMUP, N_SAMPLES, workers)
    t_nuts = time.time() - t0
    chains = []
    for k, c in enumerate(raw):
        c = dict(c, start=ref["starts"][k], seed=int(seeds[k]))
        c["mean"], c["sd"] = c["samples"].mean(axis=0), c["samples"].std(axis=0, ddof=1)
        chains.append(c)
    basin_of_chain = nearest_basin(np.array([c["mean"] for c in chains]), basins)
    for c, b in zip(chains, basin_of_chain):
        c["basin"] = int(b)
    all_draws = np.concatenate([c["samples"] for c in chains])
    draw_basin = nearest_basin(all_draws, basins)
    chain_frac = np.bincount(draw_basin, minlength=len(basins)) / all_draws.shape[0]
    grid_frac = np.array([b["mass_frac"] for b in basins])
    agree, worst_pair_z = chains_agree(chains)
    one_basin = len(set(basin_of_chain.tolist())) == 1
    n_div = sum(c["n_divergent"] for c in chains)
    n_err = sum(c["n_loglike_errors_sampling"] for c in chains)

    pooled = pooled_stats(chains, g_mean, g_sd)
    ok_point = spread < POINTWISE_TOL
    ok_err = n_err == 0
    if one_basin and agree:
        ok_mean = bool(np.all(np.abs(pooled["z_mean"]) < Z_MAX))
        ok_sd = bool(np.all((pooled["sd_ratio"] > SD_LO) & (pooled["sd_ratio"] < SD_HI)))
        status = "PASS" if ok_point and ok_mean and ok_sd and ok_err else "FAIL"
        within = None
    else:
        ok_mean = ok_sd = False
        status = "MULTIMODAL" if not one_basin else "FAIL"
        # chains in the dominant basin (by draw count) against that basin's grid moments
        dom = int(np.argmax(chain_frac))
        in_dom = [c for c in chains if c["basin"] == dom]
        within = dict(basin=dom, n_chains=len(in_dom))
        if in_dom:
            within.update(pooled_stats(in_dom, basins[dom]["mean"], basins[dom]["sd"]))
            within["agree"], within["worst_pair_z"] = chains_agree(in_dom)
            within["ok"] = bool(within["agree"]
                                and np.all(np.abs(within["z_mean"]) < Z_MAX)
                                and np.all((within["sd_ratio"] > SD_LO)
                                           & (within["sd_ratio"] < SD_HI)))

    print(f"  noise x{noise_scale:g} N_beta={n_beta} {mode:5s}: grid mean {fmt(g_mean, 5)} sd "
          f"{fmt(g_sd, 5)} [{len(basins)} basin(s), {t_grid:.0f} s]; pointwise spread "
          f"{spread:.1e} (tol {POINTWISE_TOL:g}); NUTS {N_CHAINS}x{N_SAMPLES} pooled ESS "
          f"{fmt(pooled['ess'], 0)} mean z {fmt(pooled['z_mean'], 2)} sd ratio "
          f"{fmt(pooled['sd_ratio'], 3)} div={n_div} errors={n_err} chains agree "
          f"{'yes' if agree else 'NO'} (worst pair z {worst_pair_z:.1f}) [{t_nuts:.0f} s]  {status}"
          + ("" if status == "PASS" else
             f"  (point {ok_point} mean {ok_mean} sd {ok_sd} err {ok_err})"))
    for j, b in enumerate(basins):
        print(f"      basin {j}: max {fmt(b['centre_coarse'], 2)} mean {fmt(b['mean'], 5)} sd "
              f"{fmt(b['sd'], 5)} log mass rel {b['log_mass'] - basins[0]['log_mass']:+.1f} "
              f"(grid frac {b['mass_frac']:.3g}, chain draws frac {chain_frac[j]:.3g}), lp_max "
              f"{b['lp_max'] - basins[0]['lp_max']:+.1f}, edge mass {b['edge_mass']:.1e}")
    for k, c in enumerate(chains):
        print(f"      chain {k}: start {fmt(c['start'], 4)} -> mean {fmt(c['mean'], 5)} sd "
              f"{fmt(c['sd'], 5)} ESS {fmt(c['ess'], 0)} basin {c['basin']} div "
              f"{c['n_divergent']} step {c['step_size']:.3g}")
    if within is not None and within["n_chains"]:
        print(f"      within dominant basin {within['basin']} ({within['n_chains']} chains): mean z"
              f" {fmt(within['z_mean'], 2)} sd ratio {fmt(within['sd_ratio'], 3)} ESS "
              f"{fmt(within['ess'], 0)}  {'ok' if within['ok'] else 'NOT ok'}")
    return dict(n_beta=n_beta, mode=mode, noise_scale=noise_scale, grid_mean=g_mean,
                grid_sd=g_sd, n_basins=len(basins),
                basins=[dict(centre_coarse=b["centre_coarse"], mean=b["mean"], sd=b["sd"],
                             log_mass=b["log_mass"], mass_frac=b["mass_frac"],
                             lp_max=b["lp_max"], edge_mass=b["edge_mass"],
                             chain_draw_frac=float(chain_frac[j]))
                        for j, b in enumerate(basins)],
                grid_edge_mass=top["edge_mass"], grid_points=WINDOW_POINTS,
                fine_spacing=FINE_SPACING, pointwise_spread_rel=spread, pointwise_n=len(idx),
                nuts_mean=pooled["mean"], nuts_sd=pooled["sd"], ess=pooled["ess"],
                z_mean=pooled["z_mean"], sd_ratio=pooled["sd_ratio"],
                chains=[dict(start=c["start"], seed=c["seed"], mean=c["mean"], sd=c["sd"],
                             ess=c["ess"], basin=c["basin"], n_divergent=c["n_divergent"],
                             n_loglike_errors_sampling=c["n_loglike_errors_sampling"],
                             step_size=c["step_size"], ebfmi=c["ebfmi"],
                             accept_stat=c["accept_stat"]) for c in chains],
                chains_agree=agree, worst_pair_z=worst_pair_z, one_basin=one_basin,
                chain_basin_frac=chain_frac, grid_basin_frac=grid_frac,
                within_dominant_basin=within, n_divergent=int(n_div),
                n_loglike_errors_sampling=int(n_err), n_warmup=N_WARMUP, n_samples=N_SAMPLES,
                n_chains=N_CHAINS, seconds_grid=t_grid, seconds_nuts=t_nuts, status=status,
                passed=status == "PASS")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=N_CHAINS, help="parallel chain processes")
    ap.add_argument("--seed", type=int, default=1202)
    args = ap.parse_args()
    t_all = time.time()
    rng = np.random.default_rng(args.seed)
    print(f"check 1b: nside_d={NSIDE_D}, 12 bands, gauss amplitudes, noise scales {NOISE_SCALES}; "
          f"NUTS {N_CHAINS} chains x ({N_WARMUP}+{N_SAMPLES}) from grid quantiles "
          f"{START_QUANTILES}; thresholds |z| < {Z_MAX}, sd ratio in [{SD_LO}, {SD_HI}], "
          f"pointwise spread < {POINTWISE_TOL:g}; coarse spacing {COARSE[1] - COARSE[0]:.3g}, "
          f"basin refinement {FINE_SPACING:g}, final {WINDOW_POINTS} points over "
          f"+-{WINDOW_SD:g} sd")
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
                row = run_case(n_beta, data, mode, mu_vec, cov, rng, noise_scale, args.workers)
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
    statuses = [r["status"] for r in rows]
    overall = "FAIL" if "FAIL" in statuses else ("MULTIMODAL" if "MULTIMODAL" in statuses
                                                 else "PASS")
    unimodal = [r for r in rows if r["status"] != "MULTIMODAL"]
    worst_z = max(float(np.max(np.abs(r["z_mean"]))) for r in unimodal)
    sd_range = (min(float(np.min(r["sd_ratio"])) for r in unimodal),
                max(float(np.max(r["sd_ratio"])) for r in unimodal))
    min_ess = min(float(np.min(r["ess"])) for r in rows)
    summary = dict(check="1b_grid", setting=dict(nside_d=NSIDE_D, noise_scales=NOISE_SCALES,
                                                 amp_mode="gauss", n_warmup=N_WARMUP,
                                                 n_samples=N_SAMPLES, n_chains=N_CHAINS,
                                                 start_quantiles=START_QUANTILES,
                                                 modes=["gauss", "w0", "w1"], n_beta=[1, 2],
                                                 seed=args.seed),
                   statistic=dict(worst_abs_z=worst_z, sd_ratio_range=sd_range, min_ess=min_ess,
                                  worst_pointwise_spread=max(r["pointwise_spread_rel"]
                                                             for r in rows),
                                  statuses=statuses,
                                  multimodal_cases=[f"noise x{r['noise_scale']:g} "
                                                    f"N_beta {r['n_beta']} {r['mode']}"
                                                    for r in rows
                                                    if r["status"] == "MULTIMODAL"]),
                   threshold=f"|z| < {Z_MAX}, sd ratio in [{SD_LO}, {SD_HI}], pointwise < "
                             f"{POINTWISE_TOL:g}, zero sampling-phase likelihood errors, chains "
                             f"in one basin agreeing within {Z_CHAIN:g} combined SE",
                   status=overall, passed=overall == "PASS", rows=rows,
                   seconds=time.time() - t_all)
    with open(OUT, "w") as f:
        json.dump(to_py(summary), f, indent=1)
    print(f"check 1b grid: statuses {statuses}; over the pooled cases worst |z| = {worst_z:.2f} "
          f"(threshold {Z_MAX}), sd ratio range [{sd_range[0]:.3f}, {sd_range[1]:.3f}] "
          f"(threshold [{SD_LO}, {SD_HI}]), min pooled ESS {min_ess:.0f}  {overall}  "
          f"[{time.time() - t_all:.0f} s, {OUT}]")
    return {"PASS": 0, "FAIL": 1, "MULTIMODAL": 2}[overall]


if __name__ == "__main__":
    sys.exit(main())
