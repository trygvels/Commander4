"""E012 v2 check 1d: linear-Gaussian control with known fixed amplitudes, closed form vs NUTS.

Model: d_pb = offset_pb + F_pb beta_r(p) + n_pb, n ~ N(0, N) (reference.simulate_linear, nside_d
16, B = 4, F ~ 50 U(0.5, 1.5), noise sd 30: each index pixel is pinned to about 0.02 by its 256
data pixels while sigma is 0.3, the E012 1c geometry). Prior on beta and theta as in E012 v2.

Reference (zero MC error): reference.linear_gaussian_closed_form on a 3-D grid of theta =
(mu, log sigma, u), giving theta | d (Woodbury marginal times the hyperprior) and the beta
marginal as the grid mixture of the Gaussian conditionals beta | theta, d. The grid is located in
two passes (coarse, then mean +- 8 sd per axis) and its edge mass is reported. The reference
basis is reference.field_spec, written independently.

Prototype: posterior.NonCentredPosterior and posterior.CentredPosterior with the linear
likelihood supplied as a callable (log L = -1/2 sum inv_var (d - offset - F beta_region)^2),
field.MaternField on the prototype bases, sampler.run_nuts.

Cases:
  A. graph Laplacian, nside_beta 1 (K = 11): prototype and reference bases built independently
     (both are eigenbases of the same G, so the targets are identical);
  B. harmonics, nside_beta 1, lmax 2 (K = 8), literal SPEC construction
     (HarmonicBasis(remove_constant=False)) against the reference's own harmonic basis;
  C. harmonics, nside_beta 1, lmax 2, the prototype default remove_constant=True; the reference
     has no independent construction of this basis, so the closed form is evaluated with the
     prototype's V and lam (labelled "shared basis": tests posterior/field/sampler, not the basis).

Sub-checks per case: (i) pointwise, zero MC error: prototype centred log density minus
[log p(theta | d) + log N(y; E[y | theta, d], P^-1)] from the closed form is constant over random
states (relative spread < 1e-10); (ii) NUTS means of mu, log sigma, u, log rho and the beta pixels
within 4 SE (sd / sqrt(ESS)) of the closed form; (iii) sd z-scores
(sd_nuts - sd_cf) / (sd_nuts / sqrt(2 ESS)) within 4; (iv) zero sampling-phase likelihood
errors; (v) closed-form grid edge mass < 1e-4.
Reports ESS; the comparison counts only when ESS(log sigma) > 800 (otherwise PROVISIONAL).

Writes checks/out/check_1d_linear.json.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from indexfield import field, posterior, reference, sampler      # noqa: E402

OUT = os.path.join(HERE, "out", "check_1d_linear.json")
NSIDE_D, NSIDE_BETA, B = 16, 1, 4
TRUTH = dict(mu=-3.0, sigma=0.3, rho=1.5)     # rho_min = h = 1.023 rad at nside_beta 1
Z_MAX, ESS_MIN, POINTWISE_TOL = 4.0, 800.0, 1e-10


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


def linear_loglike(data):
    """Callable beta -> (log L, dlog L/dbeta) for d = offset + F beta_region + n."""
    d, F, off, iv, region = data.d, data.F, data.offset, data.inv_var, np.asarray(data.region)
    nb = int(data.n_beta)

    def loglike(beta):
        r = d - off - F * beta[region][:, None]
        g = np.bincount(region, weights=np.sum(F * iv * r, axis=1), minlength=nb)
        return -0.5 * float(np.sum(r * r * iv)), g
    return loglike


def closed_form(data, fspec, hyp):
    """Two-pass closed form. Pass 1: mu on mu_hat +- 30 sd_hat (likelihood-only estimate of the
    field mean and its sd, 241 points), log sigma on [-6, 2], u on [-12, 12] (81 each); pass 2:
    mean +- 8 sd per axis from pass 1 with (81, 121, 121) points. Returns (cf, edge_mass)."""
    info = np.bincount(np.repeat(data.region, data.d.shape[1]),
                       weights=(data.F ** 2 * data.inv_var).ravel(), minlength=fspec.n_beta)
    bhat = np.bincount(np.repeat(data.region, data.d.shape[1]),
                       weights=(data.F * (data.d - data.offset) * data.inv_var).ravel(),
                       minlength=fspec.n_beta) / info
    mu_hat, sd_hat = float(bhat.mean()), float(1.0 / np.sqrt(info.sum()))
    cf0 = reference.linear_gaussian_closed_form(
        data.d, data.F, data.offset, data.inv_var, data.region, fspec, hyp,
        np.linspace(mu_hat - 30 * sd_hat, mu_hat + 30 * sd_hat, 241), np.linspace(-6.0, 2.0, 81),
        np.linspace(-12.0, 12.0, 81))
    m, s = cf0.theta_mean[:3], cf0.theta_sd[:3]
    grids = [np.linspace(m[k] - 8 * s[k], m[k] + 8 * s[k], n) for k, n in enumerate((81, 121, 121))]
    cf = reference.linear_gaussian_closed_form(data.d, data.F, data.offset, data.inv_var,
                                               data.region, fspec, hyp, *grids)
    w = cf.weights
    edge = float(max(w[[0, -1]].sum(), w[:, [0, -1]].sum(), w[:, :, [0, -1]].sum()))
    return cf, edge


def pointwise(post_c, data, fspec, hyp, rng, n_points=8):
    """Prototype centred log density against log p(theta | d) + log N(y; mean, P^-1)."""
    diffs, scale = [], 0.0
    for _ in range(n_points):
        x = post_c.random_init(rng)
        up = post_c.unpack(x)
        mu, ls, u = up["mu"], up["log_sigma"], up["u"]
        cf1 = reference.linear_gaussian_closed_form(
            data.d, data.F, data.offset, data.inv_var, data.region, fspec, hyp,
            np.array([mu]), np.array([ls]), np.array([u]))
        bmean, bcov = cf1.beta_given_theta(mu, ls, u)
        y_ref = fspec.V.T @ (up["beta"] - mu)
        y_mean = fspec.V.T @ (bmean - mu)
        cov_y = fspec.V.T @ bcov @ fspec.V
        r = y_ref - y_mean
        log_cond = (-0.5 * float(r @ np.linalg.solve(cov_y, r))
                    - 0.5 * np.linalg.slogdet(cov_y)[1] - 0.5 * fspec.K * math.log(2 * math.pi))
        lp = post_c.logp_and_grad(x)[0]
        diffs.append(lp - (float(cf1.logpost[0, 0, 0]) + log_cond))
        scale = max(scale, abs(lp))
    spread = (max(diffs) - min(diffs)) / scale
    return dict(constant=float(np.mean(diffs)), spread_rel=spread, scale=scale,
                n_points=n_points, passed=spread < POINTWISE_TOL)


def run_nuts_case(post, x0, rng, n_warmup, n_samples, cf, names):
    t0 = time.time()
    res = sampler.run_nuts(post, rng, n_warmup, n_samples, x0=x0)
    s = res.samples
    beta = np.array([post.beta(x) for x in s])
    log_rho = np.log(post.transform.rho(s[:, -1]))
    q = np.column_stack([s[:, -3], s[:, -2], s[:, -1], log_rho, beta])
    ref_mean = np.concatenate([cf.theta_mean, cf.beta_mean])
    ref_sd = np.concatenate([cf.theta_sd, cf.beta_sd])
    mean, sd, ess = q.mean(axis=0), q.std(axis=0, ddof=1), sampler.ess_geyer(q)
    z = (mean - ref_mean) / (sd / np.sqrt(ess))
    z_sd = (sd - ref_sd) / (sd / np.sqrt(2 * ess))
    precondition = bool(ess[1] > ESS_MIN)
    agree = bool(np.all(np.abs(z) < Z_MAX) and np.all(np.abs(z_sd) < Z_MAX))
    status = ("PASS" if agree and res.n_loglike_errors_sampling == 0 else "FAIL") \
        if precondition else "PROVISIONAL"
    name = type(post).__name__.replace("Posterior", "").lower()
    print(f"    nuts {name:10s}: {n_warmup}+{n_samples}, step {res.step_size:.3g}, accept "
          f"{res.accept_stat:.3f}, div {res.n_divergent} (warmup {res.n_divergent_warmup}), E-BFMI "
          f"{res.ebfmi:.3f}, depth {res.tree_depth_mean:.1f}, errors "
          f"{res.n_loglike_errors_sampling}, {time.time() - t0:.0f} s; ESS(log sigma) {ess[1]:.0f} "
          f"min ESS {ess.min():.0f}; max |z| mean {np.abs(z).max():.2f} "
          f"({names[int(np.abs(z).argmax())]}), max |z| sd {np.abs(z_sd).max():.2f}, sd ratio "
          f"[{(sd / ref_sd).min():.3f}, {(sd / ref_sd).max():.3f}]  {status}")
    for j in range(4):
        print(f"       {names[j]:9s} nuts {mean[j]:+.5f} +- {sd[j]:.5f} (ESS {ess[j]:6.0f})  cf "
              f"{ref_mean[j]:+.5f} +- {ref_sd[j]:.5f}  z {z[j]:+.2f}  sd z {z_sd[j]:+.2f}")
    print(f"       beta pixels: |z| max {np.abs(z[4:]).max():.2f}, sd ratio "
          f"[{(sd[4:] / ref_sd[4:]).min():.3f}, {(sd[4:] / ref_sd[4:]).max():.3f}], ESS min "
          f"{ess[4:].min():.0f}")
    return dict(chain=f"nuts_{name}", n_warmup=n_warmup, n_samples=n_samples, mean=mean, sd=sd,
                ess=ess, z_mean=z, z_sd=z_sd, sd_ratio=sd / ref_sd, step_size=res.step_size,
                accept_stat=res.accept_stat, n_divergent=int(res.n_divergent),
                n_divergent_warmup=int(res.n_divergent_warmup), ebfmi=res.ebfmi,
                n_leapfrog=int(res.n_leapfrog), tree_depth_mean=res.tree_depth_mean,
                n_loglike_errors_sampling=int(res.n_loglike_errors_sampling),
                precondition_met=precondition, max_abs_z=float(np.abs(z).max()),
                max_abs_z_sd=float(np.abs(z_sd).max()), status=status, seconds=time.time() - t0)


def run_case(label, basis, fspec, rng, n_warmup, n_samples):
    print(f"case {label}: K={basis.K} (reference K={fspec.K}), h={basis.h:.4f}")
    hyp = reference.default_hyper(fspec)
    data, truth = reference.simulate_linear(fspec, hyp, NSIDE_D, B, rng, TRUTH["mu"],
                                            TRUTH["sigma"], TRUTH["rho"], f_scale=50.0,
                                            noise_sd=30.0)
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior()
    loglike = linear_loglike(data)
    post_nc = posterior.NonCentredPosterior(fld, hyper, loglike=loglike)
    post_c = posterior.CentredPosterior(fld, hyper, loglike=loglike)
    if abs(post_c.hyper.rho_min - hyp.rho_min) > 1e-12:
        raise RuntimeError("hyperprior bounds differ between prototype and reference")
    t0 = time.time()
    cf, edge = closed_form(data, fspec, hyp)
    info = np.bincount(np.repeat(data.region, B), weights=(data.F ** 2 * data.inv_var).ravel(),
                       minlength=fspec.n_beta)
    print(f"    data: N={data.d.shape[0]} B={B}, likelihood sd per index pixel "
          f"{1 / np.sqrt(info).mean():.4f}; truth sd(beta) {truth['beta'].std():.4f}; closed form "
          f"{time.time() - t0:.0f} s, edge mass {edge:.1e}, theta mean (mu, log sigma, u, log rho) "
          f"{np.array2string(cf.theta_mean, precision=4)} sd "
          f"{np.array2string(cf.theta_sd, precision=4)}")
    point = pointwise(post_c, data, fspec, hyp, rng)
    print(f"    (i) pointwise: centred - closed form = {point['constant']:.6f} constant, relative "
          f"spread {point['spread_rel']:.1e} (|log pi| ~ {point['scale']:.2e}, tol "
          f"{POINTWISE_TOL:g})  {'PASS' if point['passed'] else 'FAIL'}")
    names = ["mu", "log_sigma", "u", "log_rho"] + [f"beta_{r}" for r in range(fspec.n_beta)]
    # start at the truth, mapped into the prototype's coordinates
    u_true = float(post_nc.transform.u(TRUTH["rho"]))
    s_, c_ = reference.matern_spectrum(fspec, TRUTH["rho"])
    y_true = basis.V.T @ (truth["beta"] - TRUTH["mu"])
    w_true = field.MaternField(basis).weights(TRUTH["rho"])[0]
    x0_nc = post_nc.pack(y_true / (TRUTH["sigma"] * w_true), TRUTH["mu"], math.log(TRUTH["sigma"]),
                         u_true)
    x0_c = post_c.pack(y_true, TRUTH["mu"], math.log(TRUTH["sigma"]), u_true)
    chains = [run_nuts_case(post_c, x0_c, rng, n_warmup, n_samples, cf, names),
              run_nuts_case(post_nc, x0_nc, rng, n_warmup, n_samples, cf, names)]
    statuses = [c["status"] for c in chains]
    status = "FAIL" if "FAIL" in statuses or not point["passed"] or edge > 1e-4 else \
        ("PASS" if "PASS" in statuses else "PROVISIONAL")
    return dict(case=label, family=basis.kind, K=basis.K, names=names, truth_theta=TRUTH,
                truth_beta=truth["beta"], closed_form=dict(theta_mean=cf.theta_mean,
                                                           theta_sd=cf.theta_sd,
                                                           beta_mean=cf.beta_mean,
                                                           beta_sd=cf.beta_sd, edge_mass=edge,
                                                           grid=(81, 121, 121)),
                pointwise=point, chains=chains, status=status)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=1204)
    ap.add_argument("--tag", default="", help="suffix of the output JSON name")
    args = ap.parse_args()
    out = OUT.replace(".json", f"{args.tag}.json")
    t_all = time.time()
    rng = np.random.default_rng(args.seed)
    n_warmup, n_samples = 1000, int(8000 * args.scale)
    print(f"check 1d: nside_d={NSIDE_D} nside_beta={NSIDE_BETA} B={B}, truth {TRUTH}, NUTS "
          f"{n_warmup}+{n_samples}; thresholds |z| < {Z_MAX}, ESS(log sigma) > {ESS_MIN:.0f}, "
          f"pointwise spread < {POINTWISE_TOL:g}")
    cases = []
    cases.append(run_case("A graph_laplacian (independent bases)",
                          field.GraphLaplacianBasis(NSIDE_BETA),
                          reference.field_spec(NSIDE_BETA, "graph_laplacian"), rng, n_warmup,
                          n_samples))
    cases.append(run_case("B harmonic lmax 2, literal SPEC (remove_constant=False, independent)",
                          field.HarmonicBasis(NSIDE_BETA, lmax=2, remove_constant=False),
                          reference.field_spec(NSIDE_BETA, "harmonic", lmax=2), rng, n_warmup,
                          n_samples))
    hb = field.HarmonicBasis(NSIDE_BETA, lmax=2, remove_constant=True)
    fspec_shared = reference.FieldSpec(hb.V.copy(), hb.lam.copy(), hb.h, hb.K, hb.n_beta,
                                       NSIDE_BETA, "harmonic", hb.ells)
    cases.append(run_case("C harmonic lmax 2, remove_constant=True (SHARED basis)", hb,
                          fspec_shared, rng, n_warmup, n_samples))
    statuses = [c["status"] for c in cases]
    overall = "FAIL" if "FAIL" in statuses else ("PASS" if all(s == "PASS" for s in statuses)
                                                 else "PROVISIONAL")
    worst_z = max(ch["max_abs_z"] for c in cases for ch in c["chains"])
    summary = dict(check="1d_linear", setting=dict(nside_d=NSIDE_D, nside_beta=NSIDE_BETA, B=B,
                                                   f_scale=50.0, noise_sd=30.0, truth=TRUTH,
                                                   n_warmup=n_warmup, n_samples=n_samples,
                                                   seed=args.seed),
                   statistic=dict(worst_abs_z=worst_z,
                                  ess_log_sigma=[[ch["ess"][1] for ch in c["chains"]]
                                                 for c in cases],
                                  worst_pointwise_spread=max(c["pointwise"]["spread_rel"]
                                                             for c in cases)),
                   threshold=f"|z| < {Z_MAX} on means and sds with ESS(log sigma) > {ESS_MIN:.0f}; "
                             f"pointwise spread < {POINTWISE_TOL:g}",
                   cases=cases, status=overall, passed=overall == "PASS",
                   seconds=time.time() - t_all)
    with open(out, "w") as f:
        json.dump(to_py(summary), f, indent=1)
    print(f"check 1d linear: cases {statuses}, worst |z| {worst_z:.2f} (threshold {Z_MAX})  "
          f"{overall}  [{time.time() - t_all:.0f} s, {out}]")
    return {"PASS": 0, "FAIL": 1, "PROVISIONAL": 2}[overall]


if __name__ == "__main__":
    sys.exit(main())
