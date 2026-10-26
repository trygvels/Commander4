"""E017: coverage of one full-sky synchrotron index under four amplitude treatments (M2.2).

One index for the whole sky (region = 0 for every data pixel, n_beta = 1) on the template skies
of coverage.simulate_template (amp_mode "template": a_p = m_p 10^(0.1 N(0,1)) for sync and dust,
cmb ~ N(0, 70^2)), beta_true ~ N(-3, 0.2^2), noise scale per cell. Posteriors:
    w0, w1, gauss : one-dimensional, computed exactly on an adaptive grid (no sampler), with the
                    prior beta ~ N(BETA_MU, BETA_SD^2) (the generator's);
    hier, hier_distort : the E016 target with a scalar index, five parameters
                    (beta, mu_sync, mu_dust, log s_sync, log s_dust), sampled with NUTS.
Per cell a JSON with truth, posterior mean, sd, central quantiles and (for hier) ESS and the
sampler diagnostics, in the layout of coverage.py so that the same summary statistics apply.
Replicates are independent (own truth, own noise), so the plain SE over replicates is correct.

Self-test (--selftest): with the amplitude hyperpriors pinned at the generator moments the hier
posterior of beta must match the gauss-mode grid with m = mu* g, T = (s* g)^2 on the same data
(E017 "Reference").

Run from nextgen/prototypes as
    python -m indexfield.scalar_coverage --n-rep 24 --noise 1 3 10 30 --jobs 2 --tag study_v1
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np
import healpy as hp

from indexfield import coverage, model, posterior, sampler, simulate
from indexfield.model import PixelData, mixing

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.join(HERE, "checks", "out", "scalar")
MODES = ("w0", "w1", "gauss", "hier", "hier_distort")
GRID_MODES = ("w0", "w1", "gauss")
BETA_MU, BETA_SD = -3.0, 0.2
QUANTILES = coverage.QUANTILES
# Generator moments of the log-normal scatter 10^(0.1 N): used only by the pinned self-test.
GEN_MU, GEN_S = 1.0269, 0.2336


def simulate_scalar(nside_d: int, rng, noise_scale: float, amp_prior: str = "default",
                    beta_d: float = 1.6, T_d: float = 20.0):
    """One replicate with a single index: (PixelData with region 0 everywhere, truth dict)."""
    N = hp.nside2npix(nside_d)
    beta = BETA_MU + BETA_SD * rng.standard_normal()
    m = simulate.templates(nside_d)
    sqrt_T = coverage.amp_prior_sqrt_T(m, amp_prior)
    a = m * 10.0 ** (0.1 * rng.standard_normal(m.shape))
    a[:, 0] = 70.0 * rng.standard_normal(N)
    nu = simulate.BANDS_GHZ
    rms = simulate.noise_rms(nside_d, nu, noise_scale)
    M, _ = mixing(nu, np.full(N, beta), beta_d=beta_d, T_d=T_d)
    d = np.einsum("nbc,nc->nb", M, a) + rms * rng.standard_normal(rms.shape)
    data = PixelData(nu=nu, d=d, inv_var=1.0 / rms ** 2, region=np.zeros(N, dtype=int),
                     n_beta=1, beta_d=beta_d, T_d=T_d, m=m, inv_T=1.0 / sqrt_T ** 2)
    return data, dict(beta=float(beta), a=a)


def snr_full_sky(data: PixelData, beta_true: float) -> float:
    """Full-sky analogue of coverage.snr_proxy's snr_index: sqrt of the largest per-band sum over
    pixels of (template sync amplitude at that band / noise rms)^2."""
    f = (data.nu[None, :] / 30.0) ** beta_true
    s2 = (np.abs(data.m[:, 1])[:, None] * f) ** 2 * data.inv_var
    return float(np.sqrt(s2.sum(axis=0).max()))


def log_prior_beta(beta):
    r = (beta - BETA_MU) / BETA_SD
    return -0.5 * r * r, -r / BETA_SD


def _logpost_grid(data: PixelData, mode: str, grid: np.ndarray) -> np.ndarray:
    out = np.empty(grid.size)
    with np.errstate(over="ignore", invalid="ignore"):
        for i, b in enumerate(grid):
            try:
                out[i] = model.collapsed_loglike(np.array([b]), data, mode)[0] \
                    + log_prior_beta(b)[0]
            except ValueError:
                out[i] = -np.inf
    return out


def grid_summary(grid: np.ndarray, logp: np.ndarray) -> dict:
    w = np.exp(logp - logp.max())
    w /= np.trapezoid(w, grid)
    mean = float(np.trapezoid(w * grid, grid))
    sd = float(np.sqrt(np.trapezoid(w * (grid - mean) ** 2, grid)))
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(grid))])
    cdf /= cdf[-1]
    q = np.interp(np.array(QUANTILES) / 100.0, cdf, grid)
    return dict(mean=mean, sd=sd, quantiles=q)


def grid_posterior(data: PixelData, mode: str, lo=-4.5, hi=-1.5, n_coarse=3001, n_fine=1601,
                   half_sd=8.0):
    """Exact 1D posterior of beta on a coarse-to-fine grid; the fine window is widened until it
    holds `half_sd` posterior sd on each side of the mean and the spacing is below sd / 50."""
    coarse = np.linspace(lo, hi, n_coarse)
    lp = _logpost_grid(data, mode, coarse)
    i = int(np.nanargmax(lp))
    if i in (0, coarse.size - 1):
        # The posterior piles up at the boundary of the support (the flat-prior w1 collapse of
        # F12/F15 for a scalar index): report the posterior truncated to [lo, hi] and flag it.
        s = grid_summary(coarse, lp)
        s.update(grid_lo=float(lo), grid_hi=float(hi), n_grid=int(n_coarse), n_evals=int(n_coarse),
                 sd_coarse=float("nan"), at_edge=True)
        return s
    # curvature estimate from the three points around the maximum
    h = coarse[1] - coarse[0]
    curv = (lp[i - 1] - 2 * lp[i] + lp[i + 1]) / h ** 2
    sd0 = 1.0 / np.sqrt(-curv) if curv < 0 else 5 * h
    centre, half = coarse[i], max(half_sd * sd0, 3 * h)
    for _ in range(8):
        fine = np.linspace(centre - half, centre + half, n_fine)
        lpf = _logpost_grid(data, mode, fine)
        s = grid_summary(fine, lpf)
        ok_width = (s["mean"] - half_sd * s["sd"] > fine[0]) and \
            (s["mean"] + half_sd * s["sd"] < fine[-1])
        ok_spacing = (fine[1] - fine[0]) < s["sd"] / 50.0
        if ok_width and ok_spacing:
            s.update(grid_lo=float(fine[0]), grid_hi=float(fine[-1]), n_grid=int(n_fine),
                     n_evals=int(n_coarse + n_fine), sd_coarse=float(sd0), at_edge=False)
            return s
        centre = s["mean"]
        half = (half_sd + 1.0) * s["sd"] if ok_spacing else half
        n_fine = n_fine if ok_spacing else 2 * n_fine - 1
    raise RuntimeError(f"{mode}: fine grid did not converge (sd {s['sd']:.3g}, half {half:.3g})")


class ScalarHierPosterior:
    """E016 hier target with one index: x = (beta, mu_sync, mu_dust, log s_sync, log s_dust).
    log p = log L(beta, mu, s) + log N(beta; BETA_MU, BETA_SD^2) + log p(mu) + log p(log s), with
    the amplitude hyperpriors of posterior._HierMixin (mu_c ~ N(1, 0.5^2), log s_c ~
    N(log 0.3, 1^2)) unless `pin` = (mu*, s*) replaces them by N(mu*, 1e-4^2), N(log s*, 1e-4^2)
    (self-test only)."""

    dim = 5
    K = 1
    n_beta = 1

    def __init__(self, data: PixelData, g: np.ndarray, pin=None):
        self.data = data
        self.g = np.asarray(g, dtype=float)
        H = posterior._HierMixin
        self.mu0, self.mu_sd = H.HIER_MU0, H.HIER_MU_SD
        self.ls0, self.ls_sd = np.log(H.HIER_S0), H.HIER_LOGS_SD
        if pin is not None:
            self.mu0, self.mu_sd = float(pin[0]), 1e-4
            self.ls0, self.ls_sd = float(np.log(pin[1])), 1e-4
        self.n_errors = 0

    def logp_and_grad(self, x):
        x = np.asarray(x, dtype=float)
        if x.shape != (5,):
            raise ValueError("x must have shape (5,)")
        beta, mu, ls = x[0], x[1:3], x[3:5]
        if np.any(ls < posterior.HIER_LOGS_MIN):
            self.n_errors += 1
            return -np.inf, np.zeros(5)
        try:
            with np.errstate(over="ignore", invalid="ignore"):
                logL, g_beta, g_mu, g_ls = model.hier_loglike(np.array([beta]), self.data,
                                                              self.g, mu, ls)
        except ValueError:
            self.n_errors += 1
            return -np.inf, np.zeros(5)
        if not np.isfinite(logL):
            self.n_errors += 1
            return -np.inf, np.zeros(5)
        lpb, gb = log_prior_beta(beta)
        r_mu, r_ls = (mu - self.mu0) / self.mu_sd, (ls - self.ls0) / self.ls_sd
        lp = logL + lpb - 0.5 * float(r_mu @ r_mu + r_ls @ r_ls)
        grad = np.concatenate([[g_beta[0] + gb], g_mu - r_mu / self.mu_sd,
                               g_ls - r_ls / self.ls_sd])
        return float(lp), grad

    def random_init(self, rng):
        return np.array([rng.normal(BETA_MU, BETA_SD), rng.normal(self.mu0, self.mu_sd),
                         rng.normal(self.mu0, self.mu_sd), rng.normal(self.ls0, self.ls_sd),
                         rng.normal(self.ls0, self.ls_sd)])


HIER_NAMES = ["beta", "mu_sync", "mu_dust", "log_s_sync", "log_s_dust"]


def sample_hier(data, g, rng, n_warmup, n_samples, target_accept=0.8, max_depth=10, pin=None):
    post = ScalarHierPosterior(data, g, pin=pin)
    for n_init in range(1, 21):
        x0 = post.random_init(rng)
        if np.isfinite(post.logp_and_grad(x0)[0]):
            break
    else:
        raise RuntimeError("no finite random initial state in 20 draws")
    # Pinned hyperpriors (self-test) have sd 1e-4 against a beta scale of 0.1: start the diagonal
    # mass matrix at those scales so the first warm-up window is not at tree depth 10.
    inv_mass = None if pin is None else np.array([BETA_SD ** 2] + [1e-8] * 4)
    with np.errstate(over="ignore"):
        res = sampler.run_nuts(post, rng, n_warmup, n_samples, x0=x0, inv_mass=inv_mass,
                               target_accept=target_accept, max_depth=max_depth)
    s = res.samples
    summ = dict(mean=s.mean(axis=0), sd=s.std(axis=0, ddof=1),
                quantiles=np.percentile(s, QUANTILES, axis=0), ess=sampler.ess_geyer(s),
                n_divergent=res.n_divergent, n_divergent_warmup=res.n_divergent_warmup,
                n_loglike_errors=res.n_loglike_errors,
                n_loglike_errors_sampling=res.n_loglike_errors_sampling,
                step_size=res.step_size, ebfmi=res.ebfmi, accept_stat=res.accept_stat,
                n_leapfrog=res.n_leapfrog, n_leapfrog_warmup=res.n_leapfrog_warmup,
                tree_depth_mean=res.tree_depth_mean, tree_depth_max=res.tree_depth_max,
                n_init_draws=n_init)
    return summ


def run_cell(job):
    cfg, rep, noise_scale, mode = job
    out = os.path.join(cfg["out_dir"], f"rep{rep:04d}_noise{noise_scale:g}_{mode}.json")
    if os.path.exists(out) and not cfg["overwrite"]:
        return dict(rep=rep, noise=noise_scale, mode=mode, skipped=True)
    t0 = time.time()
    rng_data = np.random.default_rng([cfg["seed"], rep, int(round(1000 * noise_scale))])
    data, truth = simulate_scalar(cfg["nside_d"], rng_data, noise_scale, cfg["amp_prior"])
    snr = snr_full_sky(data, truth["beta"])
    record = dict(rep=rep, noise_scale=noise_scale, mode=mode,
                  setting={k: v for k, v in cfg.items() if k not in ("out_dir",)},
                  n_beta=1, snr_full_sky=snr, quantile_levels=list(QUANTILES))
    if mode in GRID_MODES:
        s = grid_posterior(data, mode)
        record.update(names=["beta"], truth=[truth["beta"]], mean=[s["mean"]], sd=[s["sd"]],
                      quantiles=np.asarray(s["quantiles"])[:, None], at_edge=s["at_edge"],
                      grid=dict(lo=s["grid_lo"], hi=s["grid_hi"], n=s["n_grid"],
                                n_evals=s["n_evals"]))
        status = dict(z=(s["mean"] - truth["beta"]) / s["sd"], sd=s["sd"],
                      extra=" AT EDGE of the support" if s["at_edge"] else "")
    else:
        rng_g = np.random.default_rng([cfg["seed"], rep, int(round(1000 * noise_scale)), 11])
        g = coverage.shape_template(mode, cfg["nside_d"], rng_g)
        rng = np.random.default_rng([cfg["seed"], rep, int(round(1000 * noise_scale)),
                                     MODES.index(mode), 7])
        s = sample_hier(data, g, rng, cfg["n_warmup"], cfg["n_samples"], cfg["target_accept"],
                        cfg["max_depth"])
        record.update(names=HIER_NAMES, truth=[truth["beta"]] + [None] * 4,
                      mean=s["mean"], sd=s["sd"], quantiles=s["quantiles"], ess=s["ess"],
                      **{k: s[k] for k in s if k not in ("mean", "sd", "quantiles", "ess")})
        status = dict(z=(s["mean"][0] - truth["beta"]) / s["sd"][0], sd=s["sd"][0],
                      extra=f", ESS min {s['ess'].min():.0f}, div {s['n_divergent']}, depth "
                            f"{s['tree_depth_mean']:.1f}, E-BFMI {s['ebfmi']:.2f}, amp "
                            f"({s['mean'][1]:.3f}, {s['mean'][2]:.3f}, {s['mean'][3]:.2f}, "
                            f"{s['mean'][4]:.2f})")
    record["seconds_total"] = time.time() - t0
    with open(out, "w") as f:
        json.dump(coverage.to_py(record), f, indent=1)
    return dict(rep=rep, noise=noise_scale, mode=mode, skipped=False,
                seconds=record["seconds_total"], snr=snr, **status)


def _selftest(nside_d=16, noise_scale=3.0, seed=20261026, n_samples=1000):
    """Pinned hier posterior against the gauss grid with the generator's (mu*, s*)."""
    rng = np.random.default_rng([seed, 0, int(1000 * noise_scale)])
    data, truth = simulate_scalar(nside_d, rng, noise_scale)
    print(f"full-sky S/N at noise x{noise_scale:g}: {snr_full_sky(data, truth['beta']):.1f}",
          flush=True)
    g = simulate.templates(nside_d)
    m_star, inv_T_star = model.hier_prior(g, [GEN_MU, GEN_MU], [np.log(GEN_S)] * 2)
    data_star = PixelData(nu=data.nu, d=data.d, inv_var=data.inv_var, region=data.region,
                          n_beta=1, beta_d=data.beta_d, T_d=data.T_d, m=m_star, inv_T=inv_T_star)
    t0 = time.time()
    grid = grid_posterior(data_star, "gauss")
    t_grid = time.time() - t0
    t0 = time.time()
    s = sample_hier(data, g, np.random.default_rng(seed + 1), 500, n_samples,
                    pin=(GEN_MU, GEN_S))
    t_nuts = time.time() - t0
    dmean = (s["mean"][0] - grid["mean"]) / grid["sd"]
    ratio = s["sd"][0] / grid["sd"]
    ess = s["ess"][0]
    ok = abs(dmean) < max(0.1, 3.0 / np.sqrt(ess)) and abs(ratio - 1) < 0.1 and ess > 300
    print(f"pinned hier vs gauss grid (noise x{noise_scale:g}, truth {truth['beta']:.4f}): grid "
          f"{grid['mean']:.5f} +- {grid['sd']:.2e} ({grid['n_evals']} evals, {t_grid:.0f} s); "
          f"hier {s['mean'][0]:.5f} +- {s['sd'][0]:.2e}, ESS {ess:.0f}, div {s['n_divergent']}, "
          f"{t_nuts:.0f} s; mean diff {dmean:+.3f} grid sd, width ratio {ratio:.3f} "
          f"{'PASS' if ok else 'FAIL'}")
    # the three grid modes on the same data, for orientation
    for mode in GRID_MODES:
        t0 = time.time()
        gm = grid_posterior(data, mode)
        print(f"  {mode:5s}: {gm['mean']:.5f} +- {gm['sd']:.2e}, z {((gm['mean'] - truth['beta'])
              / gm['sd']):+.2f}, {gm['n_evals']} evals, {time.time() - t0:.0f} s"
              f"{' AT EDGE' if gm['at_edge'] else ''}", flush=True)
    print("ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--selftest-samples", type=int, default=1000)
    ap.add_argument("--n-rep", type=int, default=4)
    ap.add_argument("--rep-start", type=int, default=0)
    ap.add_argument("--noise", type=float, nargs="+", default=[1.0, 3.0, 10.0, 30.0])
    ap.add_argument("--modes", nargs="+", default=list(MODES), choices=MODES)
    ap.add_argument("--nside-d", type=int, default=16)
    ap.add_argument("--amp-prior", default="default", choices=["default", "wide"])
    ap.add_argument("--n-warmup", type=int, default=500)
    ap.add_argument("--n-samples", type=int, default=1000)
    ap.add_argument("--target-accept", type=float, default=0.8)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20261026)
    ap.add_argument("--tag", default="default")
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)
    if args.selftest:
        return _selftest(n_samples=args.selftest_samples)
    out_dir = os.path.join(args.out_dir, args.tag)
    os.makedirs(out_dir, exist_ok=True)
    cfg = dict(nside_d=args.nside_d, amp_prior=args.amp_prior, n_warmup=args.n_warmup,
               n_samples=args.n_samples, target_accept=args.target_accept,
               max_depth=args.max_depth, seed=args.seed, beta_prior=[BETA_MU, BETA_SD],
               out_dir=out_dir, overwrite=args.overwrite)
    with open(os.path.join(out_dir, "config.json"), "w") as f:
        json.dump(coverage.to_py(dict(cfg, n_rep=args.n_rep, rep_start=args.rep_start,
                                      noise=args.noise, modes=args.modes)), f, indent=1)
    jobs = [(cfg, rep, ns, mode) for rep in range(args.rep_start, args.rep_start + args.n_rep)
            for ns in args.noise for mode in args.modes]
    print(f"scalar coverage: {len(jobs)} cells ({args.n_rep} replicates x noise {args.noise} x "
          f"modes {args.modes}), nside_d {args.nside_d}, NUTS {args.n_warmup}+{args.n_samples}, "
          f"{args.jobs} workers -> {out_dir}", flush=True)
    t0 = time.time()
    ctx = mp.get_context("spawn")
    n_done = 0
    with ctx.Pool(min(args.jobs, len(jobs))) as pool:
        for r in pool.imap_unordered(run_cell, jobs):
            if r.get("skipped"):
                print(f"  rep {r['rep']} noise {r['noise']:g} {r['mode']:12s}: exists, skipped")
                continue
            n_done += 1
            print(f"  rep {r['rep']} noise {r['noise']:g} {r['mode']:12s}: {r['seconds']:5.0f} s, "
                  f"S/N {r['snr']:.0f}, sd {r['sd']:.2e}, z {r['z']:+.2f}{r['extra']}",
                  flush=True)
    print(f"done: {n_done} cells in {(time.time() - t0) / 60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
