"""E012 v2 Reference 3: frequentist coverage study of the index-field posterior on template data.

Generative model of the replicates (NOT the sampler's amplitude prior; that is the point):
    theta = (mu, sigma, rho) ~ hyperprior (simulate.draw_hyper), or fixed by --theta-fixed,
    z ~ N(0, I_K),   beta = mu 1 + sigma S(rho) z                       (field.MaternField),
    a_p = (cmb ~ N(0, 70^2),  sync = m_sync,p 10^(0.1 N(0,1)),  dust = m_dust,p 10^(0.1 N(0,1))),
    d_p = M_p(beta_{r(p)}) a_p + n_p,  n_p ~ N(0, diag((noise_scale rms_b)^2)),
with m the smooth plane-brightening templates of simulate.templates (0.1 dex log-normal scatter,
the harness convention of aux/validate_mcmc_likelihood.py).

For every replicate and every likelihood mode in --modes the posterior
    p(beta, mu, log sigma, u | d)  proportional to  L_mode(beta) p(beta | theta) p(theta)
is sampled with NUTS (sampler.run_nuts) on posterior.CentredPosterior (default; 1c found it the
better geometry at strong data) or NonCentredPosterior. L_mode is model.collapsed_loglike:
    gauss: amplitude prior a_p ~ N(m_p, T_p), T_p^(1/2) = max(0.5 |m_p|, 20) uK_RJ (sync, dust),
           70 uK_RJ (cmb) by default (SPEC "Decisions"); --amp-prior wide gives 10 max(|m|, 10);
    w1:    flat-prior marginal, -1/2 chi2min - 1/2 log det A;
    w0:    profile (ridge), -1/2 chi2min.

Recorded per (replicate, noise scale, mode), one JSON each in --out-dir: posterior mean, sd and
the 2.5, 16, 84, 97.5 percent quantiles of mu, log sigma, u and every beta pixel; the truth; the
per-index-pixel synchrotron S/N proxy; Geyer ESS per quantity; divergences; zero-density
evaluations in the sampling phase; step size; E-BFMI; leapfrog count; wall time.

S/N proxy (AME-free, the sync signal only), per index pixel r at the true index beta_r:
    snr_r = max_b sqrt( sum_{p in r} ( m_sync,p (nu_b / nu_0)^beta_r / sigma_{b,p} )^2 ),
the signal-to-noise of the synchrotron amplitude of the index pixel in its best band (template
mean, no scatter). The per-data-pixel maximum max_b |m_sync,p f_b| / sigma_{b,p} averaged over the
index pixel is stored as well (snr_pix_mean).

Nothing here is a pass/fail test (E012 principle 4). coverage_analyze.py turns the JSONs into bias
and coverage curves versus S/N.

Run from nextgen/prototypes as
    python -m indexfield.coverage --n-rep 4 --noise 1 --modes gauss w0 w1 --n-warmup 200 \
        --n-samples 400 --tag smoke
or as a script from the package directory.
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

HERE = os.path.dirname(os.path.abspath(__file__))
if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(HERE))

from indexfield import field, posterior, sampler, simulate      # noqa: E402
from indexfield.model import PixelData, mixing                   # noqa: E402

DEFAULT_OUT = os.path.join(HERE, "checks", "out", "coverage")
MODES = ("gauss", "w0", "w1")
QUANTILES = (2.5, 16.0, 84.0, 97.5)


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


def amp_prior_sqrt_T(m: np.ndarray, kind: str) -> np.ndarray:
    """Diagonal sd T^(1/2) (N, C) of the Gaussian amplitude prior used by mode gauss.

    "default" (SPEC Decisions): max(0.5 |m|, 20) uK_RJ for sync and dust, 70 uK_RJ for the cmb.
    "wide" (first build): 10 max(|m|, 10) uK_RJ for every component (cmb 100 uK_RJ).
    """
    if kind == "default":
        T = np.maximum(0.5 * np.abs(m), 20.0)
        T[:, 0] = 70.0
        return T
    if kind == "wide":
        return 10.0 * np.maximum(np.abs(m), 10.0)
    raise ValueError(f"unknown amp prior {kind!r}")


def theta_from_values(hyper, mu: float, sigma: float, rho: float) -> dict:
    """theta dict in the convention of simulate.draw_hyper from (mu, sigma, rho [rad])."""
    L = np.log(hyper.rho_max / hyper.rho_min)
    frac = (np.log(rho) - np.log(hyper.rho_min)) / L
    if not 0.0 < frac < 1.0:
        raise ValueError(f"rho = {np.degrees(rho):.2f} deg is outside [rho_min, rho_max]")
    return dict(mu=float(mu), sigma=float(sigma), rho=float(rho), log_sigma=float(np.log(sigma)),
                u=float(np.log(frac / (1.0 - frac))))


def simulate_template(nside_d: int, nside_beta: int, fld, hyper, rng, noise_scale: float,
                      theta: dict | None, amp_prior: str, beta_d: float = 1.6, T_d: float = 20.0):
    """One replicate with template amplitudes (module docstring). Returns (PixelData, truth).

    theta None: drawn from the hyperprior (simulate.draw_hyper); otherwise the given dict. The
    PixelData carries (m, inv_T) of the gauss-mode amplitude prior, which the generator does not
    use.
    """
    N, Nb = hp.nside2npix(nside_d), hp.nside2npix(nside_beta)
    K = fld.basis.V.shape[1]
    region = simulate.region_map(nside_d, nside_beta)
    th = simulate.draw_hyper(hyper, rng) if theta is None else dict(theta)
    z = rng.standard_normal(K)
    beta = th["mu"] + th["sigma"] * np.asarray(fld.S_times(th["rho"], z))
    m = simulate.templates(nside_d)
    a = m * 10.0 ** (0.1 * rng.standard_normal(m.shape))
    a[:, 0] = 70.0 * rng.standard_normal(N)
    nu = simulate.BANDS_GHZ
    rms = simulate.noise_rms(nside_d, nu, noise_scale)
    M, _ = mixing(nu, beta[region], beta_d=beta_d, T_d=T_d)
    d = np.einsum("nbc,nc->nb", M, a) + rms * rng.standard_normal(rms.shape)
    sqrt_T = amp_prior_sqrt_T(m, amp_prior)
    data = PixelData(nu=nu, d=d, inv_var=1.0 / rms ** 2, region=region, n_beta=Nb, beta_d=beta_d,
                     T_d=T_d, m=m, inv_T=1.0 / sqrt_T ** 2)
    return data, dict(theta=th, z=z, beta=beta, a=a)


def snr_proxy(data: PixelData, beta_true: np.ndarray):
    """(snr_index (Nb,), snr_pix_mean (Nb,)) of the module docstring, template sync amplitude."""
    f = (data.nu[None, :] / 30.0) ** beta_true[data.region][:, None]          # (N, B)
    s = np.abs(data.m[:, 1])[:, None] * f * np.sqrt(data.inv_var)              # (N, B)
    per_band = np.zeros((data.n_beta, data.nu.size))
    np.add.at(per_band, data.region, s ** 2)
    snr_index = np.sqrt(per_band.max(axis=1))
    counts = np.bincount(data.region, minlength=data.n_beta)
    snr_pix_mean = np.bincount(data.region, weights=s.max(axis=1), minlength=data.n_beta) / counts
    return snr_index, snr_pix_mean


def posterior_summary(post, res):
    """Mean, sd, quantiles and Geyer ESS of (mu, log sigma, u, beta_0..beta_{Nb-1})."""
    s = res.samples
    K = post.K
    beta = np.array([post.beta(x) for x in s])
    q = np.column_stack([s[:, K], s[:, K + 1], s[:, K + 2], beta])
    return dict(mean=q.mean(axis=0), sd=q.std(axis=0, ddof=1),
                quantiles=np.percentile(q, QUANTILES, axis=0), ess=sampler.ess_geyer(q),
                beta_mean_field=float(beta.mean(axis=1).mean()),
                beta_sd_field=float(beta.std(axis=1).mean()))


def run_cell(job):
    """One (replicate, noise scale, mode): simulate (seeded by replicate and noise scale only, so
    the three modes see the same data), sample, write the JSON. Returns a short status dict."""
    cfg, rep, noise_scale, mode = job
    out = os.path.join(cfg["out_dir"], f"rep{rep:04d}_noise{noise_scale:g}_{mode}.json")
    if os.path.exists(out) and not cfg["overwrite"]:
        return dict(rep=rep, noise=noise_scale, mode=mode, skipped=True)
    t0 = time.time()
    basis = field.GraphLaplacianBasis(cfg["nside_beta"]) if cfg["family"] == "graph_laplacian" \
        else field.HarmonicBasis(cfg["nside_beta"])
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior(rho_prior=cfg["rho_prior"]).resolved(basis.h)
    theta = None
    if cfg["theta_fixed"] is not None:
        mu, sigma, rho_deg = cfg["theta_fixed"]
        theta = theta_from_values(hyper, mu, sigma, np.radians(rho_deg))
    rng_data = np.random.default_rng([cfg["seed"], rep, int(round(1000 * noise_scale))])
    data, truth = simulate_template(cfg["nside_d"], cfg["nside_beta"], fld, hyper, rng_data,
                                    noise_scale, theta, cfg["amp_prior"])
    snr_index, snr_pix = snr_proxy(data, truth["beta"])
    rng = np.random.default_rng([cfg["seed"], rep, int(round(1000 * noise_scale)),
                                 MODES.index(mode), 7])
    cls = posterior.CentredPosterior if cfg["param"] == "centred" else \
        posterior.NonCentredPosterior
    post = cls(fld, hyper, data=data, mode=mode)
    th = truth["theta"]
    if cfg["init"] == "truth":
        x_nc = post.pack(truth["z"], th["mu"], th["log_sigma"], th["u"])
        x0 = post.from_noncentred(x_nc) if cfg["param"] == "centred" else x_nc
        n_init = 1
    else:
        x0, n_init = None, 0
        for n_init in range(1, 21):
            x0 = post.random_init(rng)
            if np.isfinite(post.logp_and_grad(x0)[0]):
                break
        else:
            raise RuntimeError("no finite random initial state in 20 draws")
    t1 = time.time()
    # Warm-up step-size search visits absurd beta where the sync SED overflows; the posterior
    # reports those states as zero density (counted), so the overflow warning is only noise.
    with np.errstate(over="ignore"):
        res = sampler.run_nuts(post, rng, cfg["n_warmup"], cfg["n_samples"], x0=x0,
                               target_accept=cfg["target_accept"], max_depth=cfg["max_depth"])
    t_nuts = time.time() - t1
    summ = posterior_summary(post, res)
    names = ["mu", "log_sigma", "u"] + [f"beta_{r}" for r in range(post.n_beta)]
    truth_vec = np.concatenate([[th["mu"], th["log_sigma"], th["u"]], truth["beta"]])
    record = dict(
        rep=rep, noise_scale=noise_scale, mode=mode, setting={k: v for k, v in cfg.items()
                                                                if k not in ("out_dir",)},
        K=post.K, n_beta=post.n_beta, names=names, truth=truth_vec,
        theta_truth=th, snr_index=snr_index, snr_pix_mean=snr_pix,
        mean=summ["mean"], sd=summ["sd"], quantile_levels=list(QUANTILES),
        quantiles=summ["quantiles"], ess=summ["ess"],
        beta_mean_field=summ["beta_mean_field"], beta_sd_field=summ["beta_sd_field"],
        n_divergent=res.n_divergent, n_divergent_warmup=res.n_divergent_warmup,
        n_loglike_errors=res.n_loglike_errors,
        n_loglike_errors_sampling=res.n_loglike_errors_sampling,
        step_size=res.step_size, ebfmi=res.ebfmi, accept_stat=res.accept_stat,
        n_leapfrog=res.n_leapfrog, n_leapfrog_warmup=res.n_leapfrog_warmup,
        tree_depth_mean=res.tree_depth_mean, tree_depth_max=res.tree_depth_max,
        n_init_draws=n_init, seconds_nuts=t_nuts, seconds_total=time.time() - t0)
    with open(out, "w") as f:
        json.dump(to_py(record), f, indent=1)
    i_ls = 1
    return dict(rep=rep, noise=noise_scale, mode=mode, skipped=False,
                seconds=record["seconds_total"],
                ess_min=float(np.min(summ["ess"])), ess_log_sigma=float(summ["ess"][i_ls]),
                n_div=res.n_divergent, n_err=res.n_loglike_errors_sampling, step=res.step_size,
                ebfmi=res.ebfmi, depth=res.tree_depth_mean,
                z_mu=float((summ["mean"][0] - th["mu"]) / summ["sd"][0]),
                z_ls=float((summ["mean"][1] - th["log_sigma"]) / summ["sd"][1]),
                z_u=float((summ["mean"][2] - th["u"]) / summ["sd"][2]),
                snr_range=(float(snr_index.min()), float(snr_index.max())), out=out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--n-rep", type=int, default=4, help="replicates per cell")
    ap.add_argument("--rep-start", type=int, default=0, help="first replicate index")
    ap.add_argument("--noise", type=float, nargs="+", default=[1.0, 3.0, 10.0],
                    help="noise scales relative to simulate.RMS_UK_RJ")
    ap.add_argument("--modes", nargs="+", default=list(MODES), choices=MODES)
    ap.add_argument("--nside-d", type=int, default=16)
    ap.add_argument("--nside-beta", type=int, default=2)
    ap.add_argument("--family", default="graph_laplacian", choices=["graph_laplacian", "harmonic"])
    ap.add_argument("--rho-prior", default="loguniform", choices=["loguniform", "pc"])
    ap.add_argument("--param", default="centred", choices=["centred", "noncentred"])
    ap.add_argument("--amp-prior", default="default", choices=["default", "wide"],
                    help="Gaussian amplitude prior of mode gauss")
    ap.add_argument("--theta-fixed", action="store_true",
                    help="fix the truth theta at --mu --sigma --rho-deg instead of drawing it")
    ap.add_argument("--mu", type=float, default=-3.0)
    ap.add_argument("--sigma", type=float, default=0.2)
    ap.add_argument("--rho-deg", type=float, default=30.0)
    ap.add_argument("--init", default="random", choices=["random", "truth"])
    ap.add_argument("--n-warmup", type=int, default=500)
    ap.add_argument("--n-samples", type=int, default=1000)
    ap.add_argument("--target-accept", type=float, default=0.8)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--jobs", type=int, default=4, help="parallel worker processes")
    ap.add_argument("--seed", type=int, default=20261025)
    ap.add_argument("--tag", default="default", help="subdirectory of --out-dir")
    ap.add_argument("--out-dir", default=DEFAULT_OUT)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args(argv)

    out_dir = os.path.join(args.out_dir, args.tag)
    os.makedirs(out_dir, exist_ok=True)
    cfg = dict(nside_d=args.nside_d, nside_beta=args.nside_beta, family=args.family,
               rho_prior=args.rho_prior, param=args.param, amp_prior=args.amp_prior,
               theta_fixed=([args.mu, args.sigma, args.rho_deg] if args.theta_fixed else None),
               init=args.init, n_warmup=args.n_warmup, n_samples=args.n_samples,
               target_accept=args.target_accept, max_depth=args.max_depth, seed=args.seed,
               amp_mode="template", out_dir=out_dir, overwrite=args.overwrite)
    with open(os.path.join(out_dir, "config.json"), "w") as f:
        json.dump(to_py(dict(cfg, n_rep=args.n_rep, rep_start=args.rep_start, noise=args.noise,
                             modes=args.modes)), f, indent=1)
    jobs = [(cfg, rep, ns, mode) for rep in range(args.rep_start, args.rep_start + args.n_rep)
            for ns in args.noise for mode in args.modes]
    print(f"coverage: {len(jobs)} cells ({args.n_rep} replicates x noise {args.noise} x modes "
          f"{args.modes}), nside_d {args.nside_d}, nside_beta {args.nside_beta}, {args.family}, "
          f"{args.param}, amp prior {args.amp_prior}, theta "
          f"{'fixed ' + str(cfg['theta_fixed']) if args.theta_fixed else 'drawn'}, init "
          f"{args.init}, NUTS {args.n_warmup}+{args.n_samples}, {args.jobs} workers -> {out_dir}")
    t0 = time.time()
    results = []
    ctx = mp.get_context("spawn")
    with ctx.Pool(min(args.jobs, len(jobs))) as pool:
        for r in pool.imap_unordered(run_cell, jobs):
            results.append(r)
            if r.get("skipped"):
                print(f"  rep {r['rep']} noise {r['noise']:g} {r['mode']:5s}: exists, skipped")
                continue
            print(f"  rep {r['rep']} noise {r['noise']:g} {r['mode']:5s}: {r['seconds']:6.0f} s, "
                  f"ESS min {r['ess_min']:5.0f} (log sigma {r['ess_log_sigma']:5.0f}), div "
                  f"{r['n_div']}, zero-density {r['n_err']}, step {r['step']:.3g}, depth "
                  f"{r['depth']:.1f}, E-BFMI {r['ebfmi']:.2f}, z(mu, log sigma, u) = "
                  f"({r['z_mu']:+.2f}, {r['z_ls']:+.2f}, {r['z_u']:+.2f}), S/N index "
                  f"[{r['snr_range'][0]:.1f}, {r['snr_range'][1]:.1f}]", flush=True)
    done = [r for r in results if not r.get("skipped")]
    wall = time.time() - t0
    if done:
        sec = np.array([r["seconds"] for r in done])
        per_rep = sec.sum() / max(1, len({(r["rep"], r["noise"]) for r in done}))
        print(f"coverage: {len(done)} cells in {wall:.0f} s wall with {args.jobs} workers; per "
              f"cell {sec.mean():.0f} s (min {sec.min():.0f}, max {sec.max():.0f}); per replicate "
              f"and noise scale (all modes) {per_rep:.0f} s CPU; total divergences "
              f"{sum(r['n_div'] for r in done)}, zero-density evaluations "
              f"{sum(r['n_err'] for r in done)}")
        for mode in args.modes:
            rr = [r for r in done if r["mode"] == mode]
            if rr:
                print(f"  {mode:5s}: {np.mean([r['seconds'] for r in rr]):.0f} s per cell, ESS min "
                      f"median {np.median([r['ess_min'] for r in rr]):.0f}, divergences "
                      f"{sum(r['n_div'] for r in rr)}")
    print(f"next: python -m indexfield.coverage_analyze --dir {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
