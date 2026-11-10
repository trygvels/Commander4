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
    """Naive covariance J^-1 and sandwich covariance J^-1 (J + T) J^-1 of the profile estimator, with J
    evaluated at the amplitudes that `like` fits. Exact when those amplitudes are the truth (noiseless
    data). On noisy data the fitted amplitudes carry noise and J(a_hat) has expectation J + T: use
    `sandwich_plugin`."""
    J, T, sol = like.information(beta)
    Ji = np.linalg.inv(J)
    return Ji, Ji @ (J + T) @ Ji, J, T, sol


def sandwich_plugin(like, beta):
    """Covariance of the profile estimator from one data set.

    The information evaluated at the fitted amplitudes, J(a_hat) = (M_r a_hat)^T P (M_r' a_hat), has
    expectation J + T, because a_hat = a + e with Cov(e) = A^-1 and tr(A^-1 M_r^T P M_r') = T. So
    J_c = J(a_hat) - T estimates J without bias, J(a_hat) itself estimates the score covariance K = J + T,
    and the estimator covariance J^-1 K J^-1 is estimated by J_c^-1 J(a_hat) J_c^-1.

    Returns (curvature covariance J_c^-1, sandwich covariance, J(a_hat), T, sol, ok). `ok` is False for a
    region whose J_c is not positive (the fitted amplitudes are consistent with noise: no measurement);
    the covariances are then computed on the remaining regions and filled with inf for the others.
    """
    Jp, T, sol = like.information(beta)
    Jc = Jp - T
    ok = np.diag(Jc) > 0
    # drop regions without information, and check that what remains is positive definite
    idx = np.where(ok)[0]
    while idx.size:
        w = np.linalg.eigvalsh(Jc[np.ix_(idx, idx)])
        if w.min() > 0:
            break
        idx = idx[np.argsort(np.diag(Jc)[idx])[1:]]
        idx.sort()
    ok = np.zeros(Jp.shape[0], dtype=bool)
    ok[idx] = True
    nr = Jp.shape[0]
    cov_c, cov_s = np.full((nr, nr), np.inf), np.full((nr, nr), np.inf)
    if idx.size:
        Ji = np.linalg.inv(Jc[np.ix_(idx, idx)])
        cov_c[np.ix_(idx, idx)] = Ji
        cov_s[np.ix_(idx, idx)] = Ji @ Jp[np.ix_(idx, idx)] @ Ji
    return cov_c, cov_s, Jp, T, sol, ok



def marginal_hessian(like, beta, h=2e-3):
    """Hessian in beta of the flat-prior marginal objective, by central differences of its analytic
    gradient. To first order the optimum of the marginal stated at a reference nu' sits at the optimum for
    nu plus H^-1 (n ln(nu'/nu)): the response uses this Hessian, J + Hess (1/2) ln det A, and not the
    profile information J (2026-11-10 stats review; the two agree when T / J is small)."""
    beta = np.asarray(beta, float)
    H = np.zeros((beta.size, beta.size))
    for i in range(beta.size):
        e = np.zeros(beta.size)
        e[i] = h
        _, gp = like.marginal(beta + e, grad=True)
        _, gm = like.marginal(beta - e, grad=True)
        H[:, i] = (gp - gm) / (2 * h)
    return 0.5 * (H + H.T)
