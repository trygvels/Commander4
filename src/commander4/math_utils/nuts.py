"""No-U-Turn Sampler (NUTS) for gradient-based sampling of non-linear parameters.

Implements Algorithm 6 of Hoffman & Gelman (2014, JMLR 15, 1593) -- NUTS with slice sampling and
dual-averaging step-size adaptation -- with a diagonal mass matrix, and Stan-style windowed
adaptation of that mass matrix during warmup. It is written to be driven one transition at a time,
so it can sit inside a Gibbs loop where the target may change slowly between iterations (e.g. the
band maps change as TOD processing updates them).

Usage::

    sampler = NUTSSampler(dim, rng, n_warmup=500)
    for it in range(n_iter):
        theta, info = sampler.transition(theta, logp_and_grad)   # adapts while it < n_warmup

where ``logp_and_grad(theta) -> (float, ndarray)`` returns the log target density (up to a
constant) and its gradient. A non-finite log density is treated as zero density.
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

LogpGrad = Callable[[NDArray], tuple[float, NDArray]]

_DELTA_MAX = 1000.0  # Energy error beyond which a trajectory is declared divergent.


@dataclass
class TransitionInfo:
    """Diagnostics of one NUTS transition."""

    accept_stat: float      # Mean Metropolis acceptance over the trajectory (adaptation target).
    n_leapfrog: int         # Gradient evaluations used.
    tree_depth: int
    divergent: bool
    step_size: float
    adapting: bool


class _DualAveraging:
    """Step-size adaptation of Hoffman & Gelman (2014), section 3.2."""

    def __init__(self, step_size: float, target: float, gamma: float = 0.05, t0: float = 10.0,
                 kappa: float = 0.75):
        self.target, self.gamma, self.t0, self.kappa = target, gamma, t0, kappa
        self.restart(step_size)

    def restart(self, step_size: float) -> None:
        self.mu = np.log(10.0 * step_size)
        self.h_bar = 0.0
        self.log_eps_bar = 0.0
        self.m = 0
        self.log_eps = np.log(step_size)

    def update(self, accept_stat: float) -> float:
        self.m += 1
        w = 1.0 / (self.m + self.t0)
        self.h_bar = (1.0 - w) * self.h_bar + w * (self.target - accept_stat)
        self.log_eps = self.mu - np.sqrt(self.m) / self.gamma * self.h_bar
        eta = self.m ** (-self.kappa)
        self.log_eps_bar = eta * self.log_eps + (1.0 - eta) * self.log_eps_bar
        return float(np.exp(self.log_eps))

    @property
    def final_step_size(self) -> float:
        return float(np.exp(self.log_eps_bar))


def _warmup_windows(n_warmup: int) -> list[int]:
    """Iterations (1-based) at which a slow mass-matrix window ends, following Stan's schedule:
    an initial fast buffer, doubling slow windows, and a terminal fast buffer."""
    if n_warmup < 20:
        return []
    if n_warmup >= 150:
        init, term, base = 75, 50, 25
    else:
        init, term = int(0.15 * n_warmup), int(0.1 * n_warmup)
        base = n_warmup - init - term
    ends, start, size = [], init, base
    last = n_warmup - term
    while start + size <= last:
        next_size = 2 * size
        if start + size + next_size > last:  # Absorb the remainder into this window.
            size = last - start
        ends.append(start + size)
        start += size
        size = next_size
    return ends


class NUTSSampler:
    """NUTS with a diagonal mass matrix, adapting step size and mass matrix during warmup.

    Args:
        dim: Number of parameters.
        rng: Random generator; every draw of the sampler comes from it.
        n_warmup: Number of initial transitions that adapt step size and mass matrix.
        step_size: Initial leapfrog step size; found heuristically on the first call if None.
        inv_mass: Initial diagonal inverse mass matrix (the parameter variances); ones if None.
        target_accept: Dual-averaging target of the mean acceptance statistic.
        max_depth: Maximum tree depth (at most 2**max_depth - 1 leapfrog steps per transition).
    """

    def __init__(self, dim: int, rng: np.random.Generator, n_warmup: int = 500,
                 step_size: float | None = None, inv_mass: NDArray | None = None,
                 target_accept: float = 0.8, max_depth: int = 10):
        self.dim = dim
        self.rng = rng
        self.n_warmup = n_warmup
        self.step_size = step_size
        self.inv_mass = np.ones(dim) if inv_mass is None else np.asarray(inv_mass, dtype=float)
        self.target_accept = target_accept
        self.max_depth = max_depth
        self.iteration = 0
        self._dual = None
        self._window_ends = _warmup_windows(n_warmup)
        self._window_samples: list[NDArray] = []

    # ---------------------------------------------------------------- public interface
    def transition(self, theta: NDArray, logp_and_grad: LogpGrad) -> tuple[NDArray, TransitionInfo]:
        """One NUTS transition from ``theta``; adapts while within the warmup."""
        theta = np.asarray(theta, dtype=float)
        if self.step_size is None:
            self.step_size = self._find_reasonable_step_size(theta, logp_and_grad)
        if self._dual is None:
            self._dual = _DualAveraging(self.step_size, self.target_accept)
        self.iteration += 1
        adapting = self.iteration <= self.n_warmup
        new_theta, info = self._nuts_step(theta, logp_and_grad, self.step_size)
        info.adapting = adapting
        if adapting:
            self._adapt(new_theta, info.accept_stat)
        return new_theta, info

    # ---------------------------------------------------------------- adaptation
    def _adapt(self, theta: NDArray, accept_stat: float) -> None:
        self.step_size = self._dual.update(accept_stat)
        n_init = 75 if self.n_warmup >= 150 else int(0.15 * self.n_warmup)
        if self._window_ends and n_init < self.iteration <= self._window_ends[-1]:
            self._window_samples.append(theta.copy())
        if self.iteration in self._window_ends:
            samples = np.array(self._window_samples)
            n = len(samples)
            if n > 2:
                var = samples.var(axis=0, ddof=1)
                # Stan's regularization towards a small isotropic metric.
                self.inv_mass = (n / (n + 5.0)) * var + 1e-3 * (5.0 / (n + 5.0))
            self._window_samples = []
            self._dual.restart(self.step_size)
        if self.iteration == self.n_warmup:
            self.step_size = self._dual.final_step_size

    def _find_reasonable_step_size(self, theta: NDArray, logp_and_grad: LogpGrad) -> float:
        """Heuristic initial step size (Hoffman & Gelman 2014, Algorithm 4)."""
        eps = 1.0
        logp, grad = logp_and_grad(theta)
        r = self.rng.normal(size=self.dim) / np.sqrt(self.inv_mass)
        joint0 = logp - self._kinetic(r)

        def log_ratio(eps):
            _, r1, logp1, _ = self._leapfrog(theta, r, grad, eps, logp_and_grad)
            joint1 = logp1 - self._kinetic(r1)
            return joint1 - joint0 if np.isfinite(joint1) else -np.inf

        direction = 1.0 if log_ratio(eps) > np.log(0.5) else -1.0
        for _ in range(100):
            lr = log_ratio(eps)
            if direction * lr <= direction * np.log(0.5):
                break
            eps *= 2.0 ** direction
        return float(eps)

    # ---------------------------------------------------------------- dynamics
    def _kinetic(self, r: NDArray) -> float:
        return 0.5 * float(np.dot(r, self.inv_mass * r))

    def _leapfrog(self, theta, r, grad, eps, logp_and_grad):
        r = r + 0.5 * eps * grad
        theta = theta + eps * self.inv_mass * r
        logp, grad = logp_and_grad(theta)
        if not np.isfinite(logp):
            return theta, r, -np.inf, np.zeros_like(theta)
        r = r + 0.5 * eps * grad
        return theta, r, float(logp), np.asarray(grad, dtype=float)

    def _no_uturn(self, theta_minus, theta_plus, r_minus, r_plus) -> bool:
        dtheta = theta_plus - theta_minus
        return (np.dot(dtheta, self.inv_mass * r_minus) >= 0.0
                and np.dot(dtheta, self.inv_mass * r_plus) >= 0.0)

    def _nuts_step(self, theta0, logp_and_grad, eps):
        logp0, grad0 = logp_and_grad(theta0)
        grad0 = np.asarray(grad0, dtype=float)
        r0 = self.rng.normal(size=self.dim) / np.sqrt(self.inv_mass)
        joint0 = logp0 - self._kinetic(r0)
        log_u = joint0 + np.log(self.rng.uniform())
        tm, tp, rm, rp, gm, gp = theta0, theta0, r0, r0, grad0, grad0
        theta_new, n, keep_going, depth = theta0, 1, True, 0
        sum_alpha, n_alpha, n_leapfrog, divergent = 0.0, 0, 0, False
        while keep_going and depth < self.max_depth:
            v = 1 if self.rng.uniform() < 0.5 else -1
            if v == -1:
                tm, rm, gm, _, _, _, t_prop, n_sub, s_sub, a, na, nl, div = self._build_tree(
                    tm, rm, gm, log_u, v, depth, eps, joint0, logp_and_grad)
            else:
                _, _, _, tp, rp, gp, t_prop, n_sub, s_sub, a, na, nl, div = self._build_tree(
                    tp, rp, gp, log_u, v, depth, eps, joint0, logp_and_grad)
            n_leapfrog += nl
            divergent |= div
            sum_alpha += a
            n_alpha += na
            if s_sub and self.rng.uniform() < n_sub / n:
                theta_new = t_prop
            n += n_sub
            keep_going = s_sub and self._no_uturn(tm, tp, rm, rp)
            depth += 1
        accept = sum_alpha / n_alpha if n_alpha > 0 else 0.0
        return theta_new, TransitionInfo(accept_stat=float(accept), n_leapfrog=n_leapfrog,
                                         tree_depth=depth, divergent=divergent,
                                         step_size=float(eps), adapting=False)

    def _build_tree(self, theta, r, grad, log_u, v, depth, eps, joint0, logp_and_grad):
        """Returns (theta-, r-, grad-, theta+, r+, grad+, proposal, n, s, sum_alpha, n_alpha,
        n_leapfrog, divergent)."""
        if depth == 0:
            t1, r1, logp1, g1 = self._leapfrog(theta, r, grad, v * eps, logp_and_grad)
            joint = logp1 - self._kinetic(r1) if np.isfinite(logp1) else -np.inf
            n1 = 1 if log_u <= joint else 0
            s1 = log_u < joint + _DELTA_MAX
            alpha = float(min(1.0, np.exp(joint - joint0))) if np.isfinite(joint) else 0.0
            return t1, r1, g1, t1, r1, g1, t1, n1, s1, alpha, 1, 1, not s1
        (tm, rm, gm, tp, rp, gp, t_prop, n1, s1, a1, na1, nl1, div1) = self._build_tree(
            theta, r, grad, log_u, v, depth - 1, eps, joint0, logp_and_grad)
        if not s1:
            return tm, rm, gm, tp, rp, gp, t_prop, n1, s1, a1, na1, nl1, div1
        if v == -1:
            tm, rm, gm, _, _, _, t2, n2, s2, a2, na2, nl2, div2 = self._build_tree(
                tm, rm, gm, log_u, v, depth - 1, eps, joint0, logp_and_grad)
        else:
            _, _, _, tp, rp, gp, t2, n2, s2, a2, na2, nl2, div2 = self._build_tree(
                tp, rp, gp, log_u, v, depth - 1, eps, joint0, logp_and_grad)
        if n1 + n2 > 0 and self.rng.uniform() < n2 / (n1 + n2):
            t_prop = t2
        s = s2 and self._no_uturn(tm, tp, rm, rp)
        return (tm, rm, gm, tp, rp, gp, t_prop, n1 + n2, s, a1 + a2, na1 + na2, nl1 + nl2,
                div1 or div2)
