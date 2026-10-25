"""E012 v2 check 1c: hyperparameters, NUTS against the independent centred Gibbs reference.

Setting: nside_beta 1 (12 index pixels, graph Laplacian family, K = 11), nside_d 16 (N = 3072,
256 data pixels per index pixel), 12 bands, mode "gauss" (the sampler-correctness target), data
drawn from the hierarchical model by simulate.simulate with the log-uniform rho prior.

Target (both paths): p(beta, mu, log sigma, u | d) with beta = mu 1 + V y, prior
N(y; 0, sigma^2 diag((s/c)^2)), the E012 hyperpriors and Jacobians. The reference chain
(reference.centred_gibbs_reference) is written independently of model.py, field.py and
posterior.py: plain-loop likelihood, own basis, random-walk MH on (y, mu) and an exact 2-D grid
draw of (log sigma, u) | y. The prototype chains are NUTS on posterior.CentredPosterior and
posterior.NonCentredPosterior (sampler.run_nuts), started at the truth like the reference.

Sub-checks:
  (i)  pointwise, zero MC error: the prototype centred log density minus the reference joint log
       density is constant over random states (relative spread < 1e-10), and the non-centred
       density differs from the centred one by exactly -K log sigma - 1/2 logdet K(rho) with the
       reference's logdet_K;
  (ii) chain comparison on mu, log sigma, u, log rho, the 12 beta pixels, the field mean and
       spatial sd of beta: |mean_nuts - mean_ref| / sqrt(SE_nuts^2 + SE_ref^2) < 4 with
       SE = sd / sqrt(ESS) (Geyer ESS, each path with its own estimator), plus the sd ratio and
       its z-score sqrt(sd_n^2 / (2 ESS_n) + sd_r^2 / (2 ESS_r)) reported;
  (iii) zero likelihood errors in the sampling phase of every NUTS chain.

Precondition (E012 1c): ESS > 800 on log sigma in BOTH chains of a comparison; when the budget
does not reach it the comparison is reported as PROVISIONAL (exit code 2), never as PASS, and the
threshold is not lowered. --scale N multiplies the chain lengths (default 2.5: Gibbs 75000
iterations after 8000 burn-in, NUTS centred 500+7500, non-centred 400+3000). The three chains run
in parallel processes (default budget about 10 minutes on a Mac mini).

Writes checks/out/check_1c_hyper.json.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from indexfield import field, posterior, reference, sampler, simulate      # noqa: E402

OUT = os.path.join(HERE, "out", "check_1c_hyper.json")
NSIDE_D, NSIDE_BETA, MODE = 16, 1, "gauss"
Z_MAX, ESS_MIN = 4.0, 800.0
POINTWISE_TOL = 1e-10
NAMES = ["mu", "log_sigma", "u", "log_rho"] + [f"beta_{r}" for r in range(12)] + \
    ["beta_mean", "beta_sd"]


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


def quantities(beta, mu, ls, u, log_rho):
    """(n, 18) matrix of the compared quantities."""
    return np.column_stack([mu, ls, u, log_rho, beta, beta.mean(axis=1), beta.std(axis=1)])


def worker_gibbs(data, truth, n_iter, n_burn, seed):
    fspec = reference.field_spec(NSIDE_BETA, "graph_laplacian")
    hyp = reference.default_hyper(fspec)
    rng = np.random.default_rng(seed)
    th = truth["theta"]
    res = reference.centred_gibbs_reference(data, fspec, hyp, MODE, rng, n_iter, n_burn=n_burn,
                                            n_mh_sub=5, init=dict(beta=truth["beta"],
                                                                  log_sigma=th["log_sigma"],
                                                                  u=th["u"]))
    q = quantities(res["beta"], res["mu"], res["log_sigma"], res["u"], res["log_rho"])
    ess = np.array([reference.ess_geyer(q[:, j]) for j in range(q.shape[1])])
    return dict(chain="gibbs_reference", n=n_iter, n_burn=n_burn, mean=q.mean(axis=0),
                sd=q.std(axis=0, ddof=1), ess=ess, accept_rate=res["accept_rate"],
                seconds=res["seconds"], K=res["K"])


def worker_nuts(data, truth, centred, n_warmup, n_samples, seed):
    basis = field.GraphLaplacianBasis(NSIDE_BETA)
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior()
    rng = np.random.default_rng(seed)
    post_nc = posterior.NonCentredPosterior(fld, hyper, data=data, mode=MODE)
    th = truth["theta"]
    x0 = post_nc.pack(truth["z"], th["mu"], th["log_sigma"], th["u"])
    if centred:
        post = posterior.CentredPosterior(fld, hyper, data=data, mode=MODE)
        x0 = post.from_noncentred(x0)
    else:
        post = post_nc
    t0 = time.time()
    res = sampler.run_nuts(post, rng, n_warmup, n_samples, x0=x0)
    beta = np.array([post.beta(x) for x in res.samples])
    s = res.samples
    log_rho = np.log(post.transform.rho(s[:, -1]))
    q = quantities(beta, s[:, -3], s[:, -2], s[:, -1], log_rho)
    ess = sampler.ess_geyer(q)
    return dict(chain="nuts_centred" if centred else "nuts_noncentred", n=n_samples,
                n_warmup=n_warmup, mean=q.mean(axis=0), sd=q.std(axis=0, ddof=1), ess=ess,
                n_divergent=res.n_divergent, n_divergent_warmup=res.n_divergent_warmup,
                step_size=res.step_size, ebfmi=res.ebfmi, accept_stat=res.accept_stat,
                n_leapfrog=res.n_leapfrog, tree_depth_mean=res.tree_depth_mean,
                n_loglike_errors=res.n_loglike_errors,
                n_loglike_errors_sampling=res.n_loglike_errors_sampling, seconds=time.time() - t0,
                K=post.K)


def pointwise(data, rng, n_points=8):
    """Prototype densities against the reference joint density at random states."""
    basis = field.GraphLaplacianBasis(NSIDE_BETA)
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior()
    post_c = posterior.CentredPosterior(fld, hyper, data=data, mode=MODE)
    post_nc = posterior.NonCentredPosterior(fld, hyper, data=data, mode=MODE)
    fspec = reference.field_spec(NSIDE_BETA, "graph_laplacian")
    hyp = reference.default_hyper(fspec)
    diffs, jac_err, scale = [], [], 0.0
    for _ in range(n_points):
        x_nc = post_nc.random_init(rng)
        x_c = post_c.from_noncentred(x_nc)
        up = post_c.unpack(x_c)
        beta, mu, ls, u = up["beta"], up["mu"], up["log_sigma"], up["u"]
        y_ref = fspec.V.T @ (beta - mu)
        lp_ref = (reference.loglike_pixel_loops(beta[data.region], data, MODE)
                  + float(reference._lp_mu(mu, hyp))
                  + float(reference.hyper_conditional_logdensity(np.array([ls]), np.array([u]),
                                                                 y_ref, fspec, hyp)[0, 0]))
        lp_c = post_c.logp_and_grad(x_c)[0]
        lp_nc = post_nc.logp_and_grad(x_nc)[0]
        diffs.append(lp_c - lp_ref)
        rho = reference.rho_from_u(u, hyp.rho_min, hyp.rho_max)
        expected = -fspec.K * ls - 0.5 * float(reference.logdet_K(fspec, rho))
        jac_err.append(abs((lp_c - lp_nc) - expected))
        scale = max(scale, abs(lp_c))
    spread = (max(diffs) - min(diffs)) / scale
    ok = spread < POINTWISE_TOL and max(jac_err) < 1e-8 * max(1.0, abs(expected))
    print(f"  (i) pointwise: centred - reference = {np.mean(diffs):.6f} constant, relative "
          f"spread {spread:.1e} over {n_points} states (|log pi| ~ {scale:.2e}, tol "
          f"{POINTWISE_TOL:g}); centred - noncentred vs -K log sigma - 1/2 logdet K: max "
          f"|err| {max(jac_err):.1e}  {'PASS' if ok else 'FAIL'}")
    return dict(constant=float(np.mean(diffs)), spread_rel=spread, scale=scale,
                jacobian_max_err=float(max(jac_err)), n_points=n_points, passed=ok)


def compare(nuts, ref):
    se2 = nuts["sd"] ** 2 / nuts["ess"] + ref["sd"] ** 2 / ref["ess"]
    z = (nuts["mean"] - ref["mean"]) / np.sqrt(se2)
    sd_ratio = nuts["sd"] / ref["sd"]
    z_sd = (nuts["sd"] - ref["sd"]) / np.sqrt(nuts["sd"] ** 2 / (2 * nuts["ess"])
                                              + ref["sd"] ** 2 / (2 * ref["ess"]))
    i_ls = NAMES.index("log_sigma")
    precondition = bool(nuts["ess"][i_ls] > ESS_MIN and ref["ess"][i_ls] > ESS_MIN)
    agree = bool(np.all(np.abs(z) < Z_MAX)) and bool(np.all(np.abs(z_sd) < Z_MAX))
    errors_ok = nuts["n_loglike_errors_sampling"] == 0
    status = ("PASS" if agree and errors_ok else "FAIL") if precondition else "PROVISIONAL"
    print(f"  (ii) {nuts['chain']} vs {ref['chain']}: ESS(log sigma) nuts {nuts['ess'][i_ls]:.0f} "
          f"ref {ref['ess'][i_ls]:.0f} (precondition > {ESS_MIN:.0f}: "
          f"{'met' if precondition else 'NOT met'}); max |z| mean {np.abs(z).max():.2f} "
          f"({NAMES[int(np.abs(z).argmax())]}), max |z| sd {np.abs(z_sd).max():.2f}, sd ratio "
          f"[{sd_ratio.min():.3f}, {sd_ratio.max():.3f}], sampling errors "
          f"{nuts['n_loglike_errors_sampling']}  {status}")
    for j, nm in enumerate(NAMES):
        print(f"       {nm:10s} nuts {nuts['mean'][j]:+.5f} +- {nuts['sd'][j]:.5f} (ESS "
              f"{nuts['ess'][j]:6.0f})  ref {ref['mean'][j]:+.5f} +- {ref['sd'][j]:.5f} (ESS "
              f"{ref['ess'][j]:6.0f})  z {z[j]:+.2f}  sd ratio {sd_ratio[j]:.3f} "
              f"(z {z_sd[j]:+.2f})")
    return dict(nuts=nuts["chain"], z_mean=z, z_sd=z_sd, sd_ratio=sd_ratio,
                ess_log_sigma=(float(nuts["ess"][i_ls]), float(ref["ess"][i_ls])),
                precondition_met=precondition, max_abs_z=float(np.abs(z).max()),
                max_abs_z_sd=float(np.abs(z_sd).max()), status=status)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", type=float, default=2.5, help="multiply all chain lengths")
    ap.add_argument("--seed", type=int, default=1203)
    ap.add_argument("--tag", default="", help="suffix of the output JSON name")
    args = ap.parse_args()
    out = OUT.replace(".json", f"{args.tag}.json")
    t_all = time.time()
    rng = np.random.default_rng(args.seed)
    basis = field.GraphLaplacianBasis(NSIDE_BETA)
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior().resolved(basis.h)
    data, truth = simulate.simulate(NSIDE_D, NSIDE_BETA, fld, hyper, rng, amp_mode="gauss",
                                    noise_scale=1.0)
    th = truth["theta"]
    print(f"check 1c: nside_d={NSIDE_D} nside_beta={NSIDE_BETA} K={basis.K} mode={MODE}; truth "
          f"mu={th['mu']:.4f} log_sigma={th['log_sigma']:.4f} u={th['u']:.4f} "
          f"(sigma {th['sigma']:.3f}, rho {np.degrees(th['rho']):.1f} deg); sd(beta) "
          f"{truth['beta'].std():.4f}; scale {args.scale}")
    point = pointwise(data, rng)

    budgets = dict(gibbs=(int(30000 * args.scale), 8000),
                   centred=(500, int(3000 * args.scale)),
                   noncentred=(400, int(1200 * args.scale)))
    ctx = mp.get_context("spawn")
    with ctx.Pool(3) as pool:
        jobs = dict(
            gibbs=pool.apply_async(worker_gibbs, (data, truth, *budgets["gibbs"], args.seed + 1)),
            centred=pool.apply_async(worker_nuts, (data, truth, True, *budgets["centred"],
                                                   args.seed + 2)),
            noncentred=pool.apply_async(worker_nuts, (data, truth, False, *budgets["noncentred"],
                                                      args.seed + 3)))
        chains = {k: j.get() for k, j in jobs.items()}
    ref = chains["gibbs"]
    print(f"  reference: {ref['n']} iterations after {ref['n_burn']} burn-in, accept "
          f"{ref['accept_rate']:.3f}, {ref['seconds']:.0f} s, min ESS {ref['ess'].min():.0f}")
    for k in ("centred", "noncentred"):
        c = chains[k]
        print(f"  {c['chain']}: {c['n_warmup']}+{c['n']} transitions, step {c['step_size']:.3g}, "
              f"accept {c['accept_stat']:.3f}, div {c['n_divergent']} (warmup "
              f"{c['n_divergent_warmup']}), E-BFMI {c['ebfmi']:.3f}, leapfrog {c['n_leapfrog']}, "
              f"depth {c['tree_depth_mean']:.1f}, loglike errors {c['n_loglike_errors']} "
              f"(sampling {c['n_loglike_errors_sampling']}), {c['seconds']:.0f} s, min ESS "
              f"{c['ess'].min():.0f}")
    comps = [compare(chains[k], ref) for k in ("centred", "noncentred")]
    statuses = [c["status"] for c in comps]
    overall = "FAIL" if "FAIL" in statuses or not point["passed"] else \
        ("PASS" if "PASS" in statuses else "PROVISIONAL")
    summary = dict(check="1c_hyper", setting=dict(nside_d=NSIDE_D, nside_beta=NSIDE_BETA,
                                                  mode=MODE, family="graph_laplacian",
                                                  rho_prior="loguniform", amp_mode="gauss",
                                                  noise_scale=1.0, budgets=budgets,
                                                  init="truth", seed=args.seed),
                   truth=dict(theta=th, beta=truth["beta"]), names=NAMES, pointwise=point,
                   chains=chains, comparisons=comps,
                   statistic=dict(max_abs_z=[c["max_abs_z"] for c in comps],
                                  ess_log_sigma=[c["ess_log_sigma"] for c in comps]),
                   threshold=f"ESS(log sigma) > {ESS_MIN:.0f} in both chains, then |z| < {Z_MAX} "
                             f"on means and sds; pointwise spread < {POINTWISE_TOL:g}",
                   status=overall, passed=overall == "PASS", seconds=time.time() - t_all)
    with open(out, "w") as f:
        json.dump(to_py(summary), f, indent=1)
    print(f"check 1c hyper: centred {statuses[0]} (max |z| {comps[0]['max_abs_z']:.2f}, ESS log "
          f"sigma {comps[0]['ess_log_sigma'][0]:.0f}/{comps[0]['ess_log_sigma'][1]:.0f}), "
          f"noncentred {statuses[1]} (max |z| {comps[1]['max_abs_z']:.2f}, ESS "
          f"{comps[1]['ess_log_sigma'][0]:.0f}/{comps[1]['ess_log_sigma'][1]:.0f}); threshold "
          f"|z| < {Z_MAX}, ESS > {ESS_MIN:.0f}  {overall}  [{time.time() - t_all:.0f} s, {out}]")
    return {"PASS": 0, "FAIL": 1, "PROVISIONAL": 2}[overall]


if __name__ == "__main__":
    sys.exit(main())
