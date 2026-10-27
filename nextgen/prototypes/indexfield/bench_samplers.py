"""M2.3 / Q4: effective samples per gradient evaluation of gradient samplers on the E016 hier
target (index field with the hierarchical amplitude prior), on template skies.

    cd nextgen/prototypes && python -m indexfield.bench_samplers --sampler nuts --nside-beta 2 \
        --noise 1 --n-warmup 500 --n-samples 2000 --tag bench_v1

Samplers are registered in SAMPLERS by name; each takes (post, rng, n_warmup, n_samples) and
returns an object with `samples` (n, dim), `n_leapfrog` (gradient evaluations in the sampling
phase), `n_leapfrog_warmup`, `n_divergent` and `step_size`. ESS is Geyer's on each coordinate
and on beta per index pixel; the figure of merit is min and median ESS per gradient over the
field block and over the hyperparameters, as principle 2 of CLAUDE.md asks (cost = gradient
evaluations, never wall time across implementations). Output: JSON per run under
checks/out/bench/<tag>/.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

from indexfield import coverage, posterior, sampler, simulate
from indexfield.field import GraphLaplacianBasis, MaternField
from indexfield.sampler import ess_geyer

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "checks", "out", "bench")


def make_target(nside_beta: int, noise_scale: float, seed: int, parameterization: str,
                nside_d: int = 16):
    """The E016 hier target on one template replicate: (post, truth)."""
    rng = np.random.default_rng(seed)
    fld = MaternField(GraphLaplacianBasis(nside_beta))
    hyper = posterior.HyperPrior(mu0=-3.0, s_mu=0.5).resolved(fld.basis.h)
    data, truth = coverage.simulate_template(nside_d, nside_beta, fld, hyper, rng, noise_scale,
                                             theta=None, amp_prior="default", amp_mode="hier")
    g = simulate.templates(nside_d)
    post = posterior.hier_posterior_class(parameterization)(fld, hyper, data, g)
    return post, truth


def run_nuts(post, rng, n_warmup, n_samples, **kw):
    return sampler.run_nuts(post, rng, n_warmup, n_samples, max_depth=kw.get("max_depth", 10),
                            target_accept=kw.get("target_accept", 0.8))


def run_mams(post, rng, n_warmup, n_samples, **kw):
    from indexfield import mams
    return mams.run_mams(post, rng, n_warmup, n_samples,
                         target_accept=kw.get("target_accept_mams", mams.TARGET_ACCEPT))


SAMPLERS = {"nuts": run_nuts, "mams": run_mams}


def summarize(post, res, truth, seconds: float) -> dict:
    s = res.samples
    K = post.K
    betas = np.array([post.beta(x[:K + 3]) for x in s])
    ess_beta = ess_geyer(betas)
    ess_x = ess_geyer(s)
    n_grad = int(res.n_leapfrog)
    hyper_idx = np.arange(K, s.shape[1])
    z = (betas.mean(0) - truth["beta"]) / betas.std(0)
    return dict(
        n_samples=int(s.shape[0]), n_grad=n_grad, n_grad_warmup=int(res.n_leapfrog_warmup),
        n_divergent=int(res.n_divergent), step_size=float(res.step_size), seconds=seconds,
        ess_beta_min=float(ess_beta.min()), ess_beta_median=float(np.median(ess_beta)),
        ess_field_block_min=float(ess_x[:K].min()), ess_hyper_min=float(ess_x[hyper_idx].min()),
        ess_hyper_median=float(np.median(ess_x[hyper_idx])),
        ess_per_grad_beta_min=float(ess_beta.min() / n_grad),
        ess_per_grad_beta_median=float(np.median(ess_beta) / n_grad),
        ess_per_grad_hyper_min=float(ess_x[hyper_idx].min() / n_grad),
        z_beta_rms=float(np.sqrt(np.mean(z**2))), z_beta_max=float(np.abs(z).max()),
        hyper_mean=s[:, hyper_idx].mean(0).tolist(), hyper_sd=s[:, hyper_idx].std(0).tolist(),
        truth_theta=truth["theta"])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--sampler", default="nuts", choices=sorted(SAMPLERS))
    ap.add_argument("--target-accept-mams", type=float, default=0.9)
    ap.add_argument("--nside-beta", type=int, default=2)
    ap.add_argument("--nside-d", type=int, default=16)
    ap.add_argument("--noise", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--parameterization", default="noncentred", choices=["noncentred", "centred"])
    ap.add_argument("--n-warmup", type=int, default=500)
    ap.add_argument("--n-samples", type=int, default=2000)
    ap.add_argument("--max-depth", type=int, default=10)
    ap.add_argument("--target-accept", type=float, default=0.8)
    ap.add_argument("--tag", default="bench")
    args = ap.parse_args(argv)
    post, truth = make_target(args.nside_beta, args.noise, args.seed, args.parameterization,
                              args.nside_d)
    rng = np.random.default_rng(1000 + args.seed)
    t0 = time.time()
    res = SAMPLERS[args.sampler](post, rng, args.n_warmup, args.n_samples,
                                 max_depth=args.max_depth, target_accept=args.target_accept,
                                 target_accept_mams=args.target_accept_mams)
    out = summarize(post, res, truth, time.time() - t0)
    out.update(sampler=args.sampler, nside_beta=args.nside_beta, nside_d=args.nside_d,
               noise=args.noise, seed=args.seed, parameterization=args.parameterization,
               n_warmup=args.n_warmup, K=int(post.K), dim=int(post.dim))
    os.makedirs(os.path.join(OUT, args.tag), exist_ok=True)
    fn = os.path.join(OUT, args.tag, f"{args.sampler}_nb{args.nside_beta}_noise{args.noise:g}_"
                      f"{args.parameterization}_seed{args.seed}.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=1, default=coverage.to_py)
    print(f"{args.sampler} nside_beta {args.nside_beta} (K {post.K}, dim {post.dim}) noise "
          f"{args.noise:g} {args.parameterization}: {out['n_grad']} gradients for "
          f"{out['n_samples']} samples ({out['seconds']:.0f} s); ESS beta min/median "
          f"{out['ess_beta_min']:.0f}/{out['ess_beta_median']:.0f}, hyper min "
          f"{out['ess_hyper_min']:.0f}; ESS per gradient beta min {out['ess_per_grad_beta_min']:.4f} "
          f"median {out['ess_per_grad_beta_median']:.4f} hyper min {out['ess_per_grad_hyper_min']:.4f}; "
          f"divergences {out['n_divergent']}, z(beta) rms {out['z_beta_rms']:.2f}; wrote {fn}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
