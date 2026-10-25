"""E012 v2 check 1a: gradients of log pi against central finite differences.

Tested:
  (1) model.collapsed_loglike gradient per index pixel, modes gauss, w0, w1, at nside_d 16 and
      nside_beta 2 (48 index pixels), at the truth and at a perturbed field;
  (2) the per-region w0/w1 log likelihood and gradient against the C++ kernel
      commander4._cmdr4_backend.compsep.marginal_loglike_regions, when importable (same
      dropped constant: -1/2 d^T N^-1 d);
  (3) both posteriors (non-centred, centred), both basis families (graph Laplacian, harmonics),
      three likelihood modes, at random states drawn by random_init(rng): |g_fd - g| <
      1e-6 |g| + 1e-4 per coordinate.

Finite differences: fourth-order central stencil, (8 (f(x+h) - f(x-h)) - (f(x+2h) - f(x-2h))) /
(12 h), h = 1e-3 (1e-4 for the hyperparameter coordinates mu, log sigma, u is not needed: the log
density is smooth there). The roundoff floor of the difference is about eps |log pi| / h, which
with |log pi| ~ 1e5 to 1e6 (template amplitudes, noise scale 1) is 1e-8 to 1e-7, well below the
1e-4 absolute term of the tolerance. Gauss-amplitude data (10 x template) give |log pi| ~ 1e8 and
an FD floor at the tolerance itself; that is why template amplitudes are used here (the gradient
check does not depend on how the data were generated).

Writes checks/out/check_1a_gradient.json. Exit code 0 on PASS, 1 on FAIL.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from indexfield import field, model, posterior, simulate      # noqa: E402

OUT = os.path.join(HERE, "out", "check_1a_gradient.json")
TOL_REL, TOL_ABS = 1e-6, 1e-4


def fd4(f, x, h):
    """Fourth-order central-difference gradient of scalar f at x."""
    g = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h
        g[i] = (8.0 * (f(x + e) - f(x - e)) - (f(x + 2 * e) - f(x - 2 * e))) / (12.0 * h)
    return g


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


def check_model_gradient(data, beta_points, h=1e-3):
    """collapsed_loglike gradient against FD per index pixel at each point in beta_points."""
    rows = []
    for mode in ("gauss", "w0", "w1"):
        for label, beta in beta_points.items():
            logL, grad = model.collapsed_loglike(beta, data, mode)
            g_fd = fd4(lambda b: model.collapsed_loglike(b, data, mode)[0], beta, h)
            err = np.abs(g_fd - grad)
            tol = TOL_REL * np.abs(grad) + TOL_ABS
            ok = bool(np.all(err < tol))
            rows.append(dict(check="model.collapsed_loglike", mode=mode, state=label,
                             logL=float(logL), max_abs_grad=float(np.abs(grad).max()),
                             max_err=float(err.max()), max_err_over_tol=float((err / tol).max()),
                             max_rel_err=float(np.max(err / np.maximum(np.abs(grad), 1e-300))),
                             n_coord=int(beta.size), passed=ok))
            print(f"  model {mode:5s} {label:9s}: logL={logL:.4e} max|g|={np.abs(grad).max():.3e} "
                  f"max|g_fd-g|={err.max():.2e} max err/tol={(err / tol).max():.2e}  "
                  f"{'PASS' if ok else 'FAIL'}")
    return rows


def check_cpp_kernel(data, beta):
    """w0/w1 per-region log likelihood and gradient against marginal_loglike_regions."""
    try:
        from commander4._cmdr4_backend import compsep
    except ImportError as err:
        print(f"  C++ kernel not importable ({err}); skipped")
        return [dict(check="cpp_kernel", skipped=True, reason=str(err))]
    rows = []
    M_reg, dM_reg = model.mixing(data.nu, beta, beta_d=data.beta_d, T_d=data.T_d)   # (Nb, B, C)
    dM_full = np.zeros_like(M_reg)
    dM_full[:, :, model.SYNC_COL] = dM_reg
    map_sky = np.ascontiguousarray(data.d.T)
    map_inv_var = np.ascontiguousarray(data.inv_var.T)
    region = np.ascontiguousarray(data.region, dtype=np.int64)
    for mode, use_logdet in (("w0", False), ("w1", True)):
        ll, gr, nskip = compsep.marginal_loglike_regions(map_sky, map_inv_var, region,
                                                         np.ascontiguousarray(M_reg),
                                                         np.ascontiguousarray(dM_full),
                                                         1, use_logdet, True)
        logL, grad = model.collapsed_loglike(beta, data, mode)
        d_ll = abs(float(np.sum(ll)) - logL)
        d_gr = float(np.max(np.abs(gr[:, model.SYNC_COL] - grad)))
        tol_ll = 1e-9 * abs(logL)
        tol_gr = 1e-9 * np.abs(grad).max() + 1e-9
        ok = nskip == 0 and d_ll < tol_ll and d_gr < tol_gr
        rows.append(dict(check="cpp_kernel", mode=mode, logL_numpy=float(logL),
                         logL_cpp=float(np.sum(ll)), abs_diff_logL=d_ll, max_abs_diff_grad=d_gr,
                         max_abs_grad=float(np.abs(grad).max()), nskipped=int(nskip), passed=ok))
        print(f"  C++ kernel {mode}: |dlogL|={d_ll:.2e} (|logL|={abs(logL):.3e}) "
              f"max|dgrad|={d_gr:.2e} (max|g|={np.abs(grad).max():.3e}) nskipped={nskip}  "
              f"{'PASS' if ok else 'FAIL'}")
    return rows


def check_posterior_gradients(data, rng, n_states=2, h=1e-3):
    rows = []
    for basis in (field.GraphLaplacianBasis(2), field.HarmonicBasis(2)):
        fld = field.MaternField(basis)
        for mode in ("gauss", "w0", "w1"):
            for cls in (posterior.NonCentredPosterior, posterior.CentredPosterior):
                post = cls(fld, posterior.HyperPrior(), data=data, mode=mode)
                worst, worst_rel, lps = 0.0, 0.0, []
                for _ in range(n_states):
                    x = post.random_init(rng)
                    lp, grad = post.logp_and_grad(x)
                    if not np.isfinite(lp):
                        continue
                    g_fd = fd4(lambda v: post.logp_and_grad(v)[0], x, h)
                    err = np.abs(g_fd - grad)
                    tol = TOL_REL * np.abs(grad) + TOL_ABS
                    worst = max(worst, float((err / tol).max()))
                    worst_rel = max(worst_rel, float(np.max(err / np.maximum(np.abs(grad),
                                                                             1e-300))))
                    lps.append(float(lp))
                ok = worst < 1.0 and len(lps) == n_states and post.n_loglike_errors == 0
                name = cls.__name__.replace("Posterior", "").lower()
                rows.append(dict(check="posterior", family=basis.kind, mode=mode,
                                 parameterization=name, K=int(post.K), dim=int(post.dim),
                                 n_states=len(lps), logp=lps, max_err_over_tol=worst,
                                 max_rel_err=worst_rel, n_loglike_errors=post.n_loglike_errors,
                                 passed=ok))
                print(f"  posterior {basis.kind:8s} {mode:5s} {name:10s} dim={post.dim:3d}: "
                      f"|logp|~{np.mean(np.abs(lps)):.2e} max err/tol={worst:.2e} "
                      f"max rel={worst_rel:.1e}  {'PASS' if ok else 'FAIL'}")
    return rows


def main():
    t0 = time.time()
    rng = np.random.default_rng(1201)
    nside_d, nside_beta = 16, 2
    basis = field.GraphLaplacianBasis(nside_beta)
    fld = field.MaternField(basis)
    hyper = posterior.HyperPrior().resolved(basis.h)
    data, truth = simulate.simulate(nside_d, nside_beta, fld, hyper, rng, amp_mode="template",
                                    noise_scale=1.0)
    beta_true = truth["beta"]
    print(f"check 1a: nside_d={nside_d} nside_beta={nside_beta} N={data.d.shape[0]} "
          f"B={data.nu.size} Nb={data.n_beta}; template amplitudes, noise scale 1; truth "
          f"mu={truth['theta']['mu']:.3f} sigma={truth['theta']['sigma']:.3f} "
          f"rho={np.degrees(truth['theta']['rho']):.1f} deg")
    print(f"tolerance |g_fd - g| < {TOL_REL:g} |g| + {TOL_ABS:g}, 4th-order central FD h=1e-3")
    rows = []
    print("(1) model.collapsed_loglike gradient per index pixel")
    points = {"truth": beta_true, "perturbed": beta_true + 0.1 * rng.standard_normal(data.n_beta)}
    rows += check_model_gradient(data, points)
    print("(2) C++ kernel marginal_loglike_regions (w0, w1)")
    rows += check_cpp_kernel(data, points["perturbed"])
    print("(3) posterior gradients at random prior states")
    rows += check_posterior_gradients(data, rng)
    tested = [r for r in rows if not r.get("skipped")]
    passed = all(r["passed"] for r in tested)
    worst = max(r.get("max_err_over_tol", 0.0) for r in tested)
    summary = dict(check="1a_gradient", setting=dict(nside_d=nside_d, nside_beta=nside_beta,
                                                     amp_mode="template", noise_scale=1.0,
                                                     fd="4th-order central, h=1e-3",
                                                     tolerance=f"{TOL_REL} |g| + {TOL_ABS}"),
                   statistic=dict(max_err_over_tol=worst, n_subchecks=len(tested)),
                   threshold="max err/tol < 1 in every sub-check", passed=passed, rows=rows,
                   seconds=time.time() - t0)
    with open(OUT, "w") as f:
        json.dump(to_py(summary), f, indent=1)
    print(f"check 1a gradient: {len(tested)} sub-checks, worst err/tol = {worst:.2e} "
          f"(threshold 1)  {'PASS' if passed else 'FAIL'}  [{time.time() - t0:.0f} s, {OUT}]")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
