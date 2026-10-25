"""E012 v2 check 1e: negative controls. Three WRONG targets that must fail check 1d visibly.

Reference: the linear-Gaussian closed form of check 1d, case A (graph Laplacian, nside_beta 1,
K = 11, nside_d 16, B = 4; data from reference.simulate_linear with the check 1d truth and seed).
It has zero Monte Carlo error, so a |z| above 4 on a mean or sd with the stated ESS is a failure
of the target, not of the comparison. The data construction, closed form and loglike are imported
from check_1d_linear so that the control runs on exactly the setting that passes with the correct
target; the correct centred target is rerun here as the positive control.

Wrong targets (thin subclasses of posterior.CentredPosterior / NonCentredPosterior that add or
remove one term in logp_and_grad; posterior.py is not modified):
  (i)   centred parameterization with -K log sigma - 1/2 log det K(rho) omitted: the classic
        error of writing the prior N(beta; mu 1, sigma^2 K) without its normalization,
        log pi_wrong = log pi + K log sigma + 1/2 log det K(rho);
  (ii)  the log rho (logit) Jacobian dropped, log pi_wrong = log pi - log[L s(u)(1 - s(u))], so
        the prior is flat in u instead of flat in log rho. In the u coordinate this target is
        IMPROPER: the likelihood depends on u only through rho(u), which saturates at rho_min and
        rho_max, so log pi_wrong is constant in both tails of u (the Jacobian's -|u| decay is what
        makes the correct target proper). NUTS therefore drifts to |u| of hundreds with saturated
        tree depth; this control runs on a capped budget (500 + 2000, max tree depth 6) and
        reports the drift (max |u|, fraction of draws with |u| > 8, where rho is within 0.03
        percent of its range from a bound) besides the z-scores. Run in both parameterizations;
  (iii) the log sigma Jacobian dropped, log pi_wrong = log pi - log sigma, so the half-normal
        prior is placed on log sigma instead of on sigma. E012 v2 estimates a shift of about
        0.2 sd of log sigma at K = 12, so this control carries the precondition ESS(log sigma)
        > 800 in the NUTS chain (the closed form has no MC error); the margin is reported even
        when the precondition is not met. Run in both parameterizations.

Each wrong target's gradient is first checked against central finite differences at random
states, so that a failure cannot be attributed to a broken gradient. Then NUTS (check 1d budget:
1000 warm-up, 8000 draws, except the capped (ii) runs) is compared with the closed form:
z = (mean_nuts - mean_cf) / (sd_nuts / sqrt(ESS)) on mu, log sigma, u, log rho and the 12 beta
pixels; z_sd = (sd_nuts - sd_cf) / (sd_nuts / sqrt(2 ESS)); sd ratios. Verdict per chain:
FAIL-AS-EXPECTED when max |z| or max |z_sd| exceeds 4 (the check 1d threshold) with the
precondition met; UNEXPECTED-PASS when everything is within 4 with the precondition met;
INCONCLUSIVE when the precondition is not met (margin still printed). The positive control must
PASS. Exit status 0 only if every wrong target fails as expected and the positive control passes.

Writes checks/out/check_1e_negative.json.
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
sys.path.insert(0, HERE)

from indexfield import field, posterior, reference, sampler      # noqa: E402
import check_1d_linear as c1d                                      # noqa: E402

OUT = os.path.join(HERE, "out", "check_1e_negative.json")
Z_MAX, ESS_MIN = c1d.Z_MAX, c1d.ESS_MIN
FD_TOL = 1e-6      # relative finite-difference tolerance on the wrong-target gradients (+1e-4 abs)


# ----------------------------------------------------------------------------- wrong targets
class _NoNormalization:
    """Mixin for the centred posterior: remove -K log sigma - 1/2 log det K(rho).

    d/dlog sigma gains +K; d/du gains +1/2 dlog det K/dlog rho dlog rho/du = sum_i (dw_i/w_i)
    dlog rho/du, with w = s/c and log det K = 2 sum_i log w_i.
    """

    wrong = "(i) centred, -K log sigma - 1/2 logdet K omitted"

    def logp_and_grad(self, x):
        lp, g = super().logp_and_grad(x)
        if not np.isfinite(lp):
            return lp, g
        _, _, ls, u = self._split(x)
        rho = float(self.transform.rho(u))
        w, dw = self.field.weights(rho)
        lp += self.K * ls + float(np.sum(np.log(w)))
        g = g.copy()
        g[self.K + 1] += self.K
        g[self.K + 2] += float(np.sum(dw / w)) * float(self.transform.dlogrho_du(u))
        return float(lp), g


class _NoRhoJacobian:
    """Mixin: drop log |dlog rho / du| = log[L s(u)(1 - s(u))] and its derivative 1 - 2 s(u)."""

    wrong = "(ii) log rho (logit) Jacobian dropped"

    def logp_and_grad(self, x):
        lp, g = super().logp_and_grad(x)
        if not np.isfinite(lp):
            return lp, g
        u = float(x[self.K + 2])
        g = g.copy()
        g[self.K + 2] -= float(self.transform.dlog_jacobian_du(u))
        return float(lp - self.transform.log_jacobian(u)), g


class _NoSigmaJacobian:
    """Mixin: drop the + log sigma Jacobian of sampling log sigma (derivative 1)."""

    wrong = "(iii) log sigma Jacobian dropped"

    def logp_and_grad(self, x):
        lp, g = super().logp_and_grad(x)
        if not np.isfinite(lp):
            return lp, g
        g = g.copy()
        g[self.K + 1] -= 1.0
        return float(lp - x[self.K + 1]), g


class CentredNoNormalization(_NoNormalization, posterior.CentredPosterior):
    pass


class CentredNoRhoJacobian(_NoRhoJacobian, posterior.CentredPosterior):
    pass


class NonCentredNoRhoJacobian(_NoRhoJacobian, posterior.NonCentredPosterior):
    pass


class CentredNoSigmaJacobian(_NoSigmaJacobian, posterior.CentredPosterior):
    pass


class NonCentredNoSigmaJacobian(_NoSigmaJacobian, posterior.NonCentredPosterior):
    pass


def term_check(post, base, rng, n_points=4):
    """The wrong target minus the correct one equals the stated removed term at random states."""
    worst = 0.0
    for _ in range(n_points):
        x = post.random_init(rng)
        lp_w = post.logp_and_grad(x)[0]
        lp_c = base.logp_and_grad(x)[0]
        _, _, ls, u = post._split(x)
        if isinstance(post, _NoNormalization):
            expected = post.K * ls + 0.5 * post.field.logdet_K(float(post.transform.rho(u)))
        elif isinstance(post, _NoRhoJacobian):
            expected = -float(post.transform.log_jacobian(u))
        else:
            expected = -ls
        worst = max(worst, abs((lp_w - lp_c) - expected) / max(1.0, abs(expected)))
    return worst


def gradient_check(post, rng, n_points=3):
    """max over random states of |g_fd - g| / (1e-6 |g| + 1e-4) (check 1a tolerance)."""
    worst = 0.0
    for _ in range(n_points):
        x = post.random_init(rng)
        _, g = post.logp_and_grad(x)
        g_fd = posterior._fd4_grad(lambda v: post.logp_and_grad(v)[0], x, 1e-4)
        worst = max(worst, float(np.max(np.abs(g_fd - g) / (FD_TOL * np.abs(g) + 1e-4))))
    return worst


# ---------------------------------------------------------------------------------- the check
def compare(post, x0, rng, n_warmup, n_samples, cf, names, label, expect_fail, ess_precond,
            max_depth=10):
    """Run NUTS on `post`, compare with the closed form, print and return the verdict."""
    t0 = time.time()
    with np.errstate(over="ignore"):          # saturated sigmoid at |u| > 700 in the improper (ii)
        res = sampler.run_nuts(post, rng, n_warmup, n_samples, x0=x0, max_depth=max_depth)
        s = res.samples
        beta = np.array([post.beta(x) for x in s])
        log_rho = np.log(post.transform.rho(s[:, -1]))
    u_abs_max = float(np.abs(s[:, -1]).max())
    frac_saturated = float(np.mean(np.abs(s[:, -1]) > 8.0))
    q = np.column_stack([s[:, -3], s[:, -2], s[:, -1], log_rho, beta])
    ref_mean = np.concatenate([cf.theta_mean, cf.beta_mean])
    ref_sd = np.concatenate([cf.theta_sd, cf.beta_sd])
    mean, sd, ess = q.mean(axis=0), q.std(axis=0, ddof=1), sampler.ess_geyer(q)
    z = (mean - ref_mean) / (sd / np.sqrt(ess))
    z_sd = (sd - ref_sd) / (sd / np.sqrt(2 * ess))
    ratio = sd / ref_sd
    shift_sd = (mean - ref_mean) / ref_sd          # bias in units of the true posterior sd
    agree = bool(np.all(np.abs(z) < Z_MAX) and np.all(np.abs(z_sd) < Z_MAX))
    precondition = bool(ess[1] > ESS_MIN) if ess_precond else True
    clean = res.n_loglike_errors_sampling == 0
    if expect_fail:
        if not precondition:
            verdict = "INCONCLUSIVE"
        else:
            verdict = "UNEXPECTED-PASS" if agree else "FAIL-AS-EXPECTED"
    else:
        verdict = ("PASS" if agree and clean else "FAIL") if precondition else "PROVISIONAL"
    worst = int(np.abs(z).argmax())
    print(f"  {label}")
    print(f"    nuts {n_warmup}+{n_samples}: step {res.step_size:.3g}, accept {res.accept_stat:.3f}"
          f", div {res.n_divergent} (warmup {res.n_divergent_warmup}), E-BFMI {res.ebfmi:.3f}, "
          f"depth {res.tree_depth_mean:.1f} (max {res.tree_depth_max}, cap {max_depth}), errors "
          f"{res.n_loglike_errors_sampling}, {time.time() - t0:.0f} s; ESS(log sigma) {ess[1]:.0f}"
          f" (precondition > {ESS_MIN:.0f}: {'met' if ess[1] > ESS_MIN else 'NOT met'}), min ESS "
          f"{ess.min():.0f}; max |u| {u_abs_max:.1f}, fraction |u| > 8: {frac_saturated:.3f}")
    for j in range(4):
        print(f"    {names[j]:9s} nuts {mean[j]:+.5f} +- {sd[j]:.5f} (ESS {ess[j]:6.0f})  cf "
              f"{ref_mean[j]:+.5f} +- {ref_sd[j]:.5f}  z {z[j]:+8.2f}  shift/sd {shift_sd[j]:+7.3f}"
              f"  sd ratio {ratio[j]:.3f}  sd z {z_sd[j]:+7.2f}")
    print(f"    beta pixels: |z| max {np.abs(z[4:]).max():.2f}, |shift/sd| max "
          f"{np.abs(shift_sd[4:]).max():.3f}, sd ratio [{ratio[4:].min():.3f}, "
          f"{ratio[4:].max():.3f}], ESS min {ess[4:].min():.0f}")
    print(f"    max |z| mean {np.abs(z).max():.2f} ({names[worst]}), max |z| sd "
          f"{np.abs(z_sd).max():.2f}, threshold {Z_MAX}  ==> {verdict}")
    return dict(label=label, chain=type(post).__name__, n_warmup=n_warmup, n_samples=n_samples,
                names=names, mean=mean, sd=sd, ess=ess, z_mean=z, z_sd=z_sd, sd_ratio=ratio,
                shift_in_sd=shift_sd, step_size=res.step_size, accept_stat=res.accept_stat,
                n_divergent=int(res.n_divergent), n_divergent_warmup=int(res.n_divergent_warmup),
                ebfmi=res.ebfmi, tree_depth_mean=res.tree_depth_mean,
                tree_depth_max=int(res.tree_depth_max), max_depth=max_depth,
                u_abs_max=u_abs_max, frac_u_saturated=frac_saturated,
                n_loglike_errors_sampling=int(res.n_loglike_errors_sampling),
                ess_log_sigma=float(ess[1]), precondition_required=ess_precond,
                precondition_met=precondition, max_abs_z=float(np.abs(z).max()),
                max_abs_z_sd=float(np.abs(z_sd).max()), agree_within_threshold=agree,
                verdict=verdict, seconds=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1204, help="check 1d seed (same data, case A)")
    ap.add_argument("--scale", type=float, default=1.0, help="multiplies the sample budgets")
    ap.add_argument("--nc-scale", type=float, default=1.0,
                    help="extra factor on the non-centred (iii) budget, whose ESS(log sigma) is "
                         "about 100 per 8000 draws here (10 gives the precondition, about 4 min)")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    out = OUT.replace(".json", f"{args.tag}.json")
    t_all = time.time()
    rng = np.random.default_rng(args.seed)
    n_warmup = 1000
    n_c, n_nc = int(8000 * args.scale), int(8000 * args.scale)
    capped = dict(n_warmup=500, n_samples=int(2000 * args.scale), max_depth=6)   # improper (ii)
    print(f"check 1e (negative controls): check 1d case A, nside_d={c1d.NSIDE_D} "
          f"nside_beta={c1d.NSIDE_BETA} B={c1d.B}, truth {c1d.TRUTH}, seed {args.seed}; "
          f"threshold |z| < {Z_MAX}, ESS(log sigma) > {ESS_MIN:.0f} where required; budgets "
          f"centred {n_warmup}+{n_c}, non-centred (iii) {n_warmup}+{int(n_nc * args.nc_scale)}, "
          f"capped (ii) {capped['n_warmup']}+{capped['n_samples']} depth {capped['max_depth']}")

    # Data and closed form exactly as check 1d case A (same rng sequence: field_spec, simulate).
    basis = field.GraphLaplacianBasis(c1d.NSIDE_BETA)
    fspec = reference.field_spec(c1d.NSIDE_BETA, "graph_laplacian")
    hyp = reference.default_hyper(fspec)
    data, truth = reference.simulate_linear(fspec, hyp, c1d.NSIDE_D, c1d.B, rng, c1d.TRUTH["mu"],
                                            c1d.TRUTH["sigma"], c1d.TRUTH["rho"], f_scale=50.0,
                                            noise_sd=30.0)
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior()
    loglike = c1d.linear_loglike(data)
    t0 = time.time()
    cf, edge = c1d.closed_form(data, fspec, hyp)
    print(f"  closed form: {time.time() - t0:.0f} s, edge mass {edge:.1e}, theta mean (mu, log "
          f"sigma, u, log rho) {np.array2string(cf.theta_mean, precision=4)} sd "
          f"{np.array2string(cf.theta_sd, precision=4)}; K={basis.K}, h={basis.h:.4f}")
    names = ["mu", "log_sigma", "u", "log_rho"] + [f"beta_{r}" for r in range(fspec.n_beta)]

    base_c = posterior.CentredPosterior(fld, hyper, loglike=loglike)
    base_nc = posterior.NonCentredPosterior(fld, hyper, loglike=loglike)
    u_true = float(base_c.transform.u(c1d.TRUTH["rho"]))
    y_true = basis.V.T @ (truth["beta"] - c1d.TRUTH["mu"])
    w_true = fld.weights(c1d.TRUTH["rho"])[0]
    x0_c = base_c.pack(y_true, c1d.TRUTH["mu"], math.log(c1d.TRUTH["sigma"]), u_true)
    x0_nc = base_nc.pack(y_true / (c1d.TRUTH["sigma"] * w_true), c1d.TRUTH["mu"],
                         math.log(c1d.TRUTH["sigma"]), u_true)

    full_c = dict(n_warmup=n_warmup, n_samples=n_c, max_depth=10)
    full_nc = dict(n_warmup=n_warmup, n_samples=int(n_nc * args.nc_scale), max_depth=10)
    # (name, wrong class, base posterior, x0, expect_fail, ESS precondition required, budget)
    targets = [
        ("positive control: correct centred target", posterior.CentredPosterior, base_c, x0_c,
         False, True, full_c),
        (_NoNormalization.wrong, CentredNoNormalization, base_c, x0_c, True, False, full_c),
        (_NoRhoJacobian.wrong + ", centred (improper in u: capped budget)", CentredNoRhoJacobian,
         base_c, x0_c, True, False, capped),
        (_NoRhoJacobian.wrong + ", non-centred (improper in u: capped budget)",
         NonCentredNoRhoJacobian, base_nc, x0_nc, True, False, capped),
        (_NoSigmaJacobian.wrong + ", centred (precondition ESS(log sigma) > 800)",
         CentredNoSigmaJacobian, base_c, x0_c, True, True, full_c),
        (_NoSigmaJacobian.wrong + ", non-centred (precondition ESS(log sigma) > 800)",
         NonCentredNoSigmaJacobian, base_nc, x0_nc, True, True, full_nc),
    ]
    results, ok = [], True
    for label, cls, base, x0, expect_fail, precond, budget in targets:
        post = cls(fld, hyper, loglike=loglike)
        pre = dict(term_spread=None, gradient_err_over_tol=None)
        if expect_fail:
            pre["term_spread"] = term_check(post, base, rng)
            pre["gradient_err_over_tol"] = gradient_check(post, rng)
            pre_ok = pre["term_spread"] < 1e-10 and pre["gradient_err_over_tol"] < 1.0
            print(f"  {label}: removed term reproduced to {pre['term_spread']:.1e} (rel), "
                  f"gradient FD err/tol {pre['gradient_err_over_tol']:.2e}  "
                  f"{'ok' if pre_ok else 'BAD'}")
            ok &= pre_ok
        r = compare(post, x0, rng, budget["n_warmup"], budget["n_samples"], cf, names, label,
                    expect_fail, precond, max_depth=budget["max_depth"])
        r.update(pre)
        results.append(r)
        ok &= r["verdict"] in ("FAIL-AS-EXPECTED", "PASS", "INCONCLUSIVE")

    verdicts = [r["verdict"] for r in results]
    n_inconclusive = verdicts.count("INCONCLUSIVE")
    overall = "PASS" if ok and n_inconclusive == 0 else ("PROVISIONAL" if ok else "FAIL")
    summary = dict(check="1e_negative",
                   setting=dict(reference="check 1d case A closed form (graph Laplacian, "
                                "nside_beta 1, K 11)", nside_d=c1d.NSIDE_D,
                                nside_beta=c1d.NSIDE_BETA, B=c1d.B, truth=c1d.TRUTH,
                                seed=args.seed, budget_full_centred=full_c,
                                budget_full_noncentred=full_nc, budget_capped_ii=capped,
                                closed_form_edge_mass=edge),
                   closed_form=dict(theta_mean=cf.theta_mean, theta_sd=cf.theta_sd,
                                    beta_mean=cf.beta_mean, beta_sd=cf.beta_sd),
                   statistic={r["label"]: dict(max_abs_z=r["max_abs_z"],
                                               max_abs_z_sd=r["max_abs_z_sd"],
                                               ess_log_sigma=r["ess_log_sigma"],
                                               verdict=r["verdict"]) for r in results},
                   threshold=f"wrong targets must give |z| > {Z_MAX} on a mean or sd; (iii) "
                             f"requires ESS(log sigma) > {ESS_MIN:.0f}",
                   results=results, verdicts=verdicts, status=overall,
                   passed=overall == "PASS", seconds=time.time() - t_all)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as f:
        json.dump(c1d.to_py(summary), f, indent=1)
    print(f"check 1e negative controls: {verdicts}  {overall}  [{time.time() - t_all:.0f} s, "
          f"{out}]")
    return {"PASS": 0, "FAIL": 1, "PROVISIONAL": 2}[overall]


if __name__ == "__main__":
    sys.exit(main())
