"""Fitting helpers: optimum of an objective in the region indices, and the sandwich covariance."""
import numpy as np
from scipy.optimize import minimize


def fit_profile(like, beta0, bounds=(-4.6, -1.6), tol=1e-10):
    """Minimize the profile objective with analytic gradients. Returns beta_hat and the optimizer result."""
    res = minimize(lambda b: like.profile(b, grad=True), np.asarray(beta0, float), jac=True, method="L-BFGS-B",
                   bounds=[bounds] * like.nreg, options={"ftol": tol, "gtol": 1e-9, "maxiter": 200})
    return res.x, res


def fit_marginal(like, beta0, log_tau0=None, bounds=(-4.6, -1.6), tau_bounds=(-5.0, 9.0), tol=1e-10):
    """Minimize the marginal objective: flat-prior marginal, or zero-mean-prior marginal jointly in
    (beta, log_tau) when the likelihood has a synchrotron prior."""
    nr = like.nreg
    if not like.sync_prior:
        res = minimize(lambda b: like.marginal(b, grad=True), np.asarray(beta0, float), jac=True, method="L-BFGS-B",
                       bounds=[bounds] * nr, options={"ftol": tol, "gtol": 1e-9, "maxiter": 200})
        return res.x, None, res

    def fun(z):
        f, gb, gt = like.marginal(z[:nr], z[nr:], grad=True)
        return f, np.concatenate([gb, gt])
    z0 = np.concatenate([np.asarray(beta0, float), np.asarray(log_tau0, float)])
    res = minimize(fun, z0, jac=True, method="L-BFGS-B", bounds=[bounds] * nr + [tau_bounds] * nr,
                   options={"ftol": tol, "gtol": 1e-9, "maxiter": 300})
    return res.x[:nr], res.x[nr:], res


def sandwich(like, beta):
    """Naive covariance J^-1 and sandwich covariance J^-1 (J + T) J^-1 of the profile estimator."""
    J, T, sol = like.information(beta)
    Ji = np.linalg.inv(J)
    return Ji, Ji @ (J + T) @ Ji, J, T, sol
