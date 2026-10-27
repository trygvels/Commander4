# MAMS (Metropolis-Adjusted Microcanonical Sampler): algorithm extraction for a NumPy port

Date: 2026-10-27. Literature extraction only; no experiment records were read.

Sources and read level:

| Source | What | How read |
|---|---|---|
| arXiv 2503.01707v2 (Robnik, Cohn-Gordon, Seljak, "Metropolis Adjusted Microcanonical Hamiltonian Monte Carlo", v2 2025-05-17) | MAMS | Abstract page plus the full HTML text (`arxiv.org/html/2503.01707v2`) read through a text-extraction tool in several passes. Equations were returned in LaTeX; prose was paraphrased by the tool. Footnote 5 could not be read verbatim (see Section 3). |
| arXiv 2212.08549v3 (Robnik, De Luca, Silverstein, Seljak, "Microcanonical Hamiltonian Monte Carlo", v3 2026-05-27) | MCHMC / MCLMC dynamics, integrators, tuning | Abstract page plus full HTML text, same method. |
| arXiv 2303.18221v3 (Robnik, Seljak, "Fluctuation without dissipation: Microcanonical Langevin Monte Carlo", PMLR 253, v3 2025-05-04) | MCLMC SDE, ergodicity, phi^4 benchmark | Abstract page plus full HTML text, same method. |
| arXiv 2412.08876 (Robnik, Cohn-Gordon, Seljak, "Practical and Scalable HMC Without the Metropolis Test") | Cited by 2503.01707 for the 10 percent tuning fractions | Abstract only. |
| BlackJAX `main` (github.com/blackjax-devs/blackjax): `blackjax/mcmc/integrators.py`, `adjusted_mclmc.py`, `adjusted_mclmc_dynamic.py`, `proposal.py`, `dynamic_hmc.py`, `blackjax/adaptation/adjusted_mclmc_adaptation.py`, `mclmc_adaptation.py`, `laps.py` | Reference implementation | Raw source files read through the extraction tool; key lines quoted by the tool. |
| github.com/reubenharry/blackjax, branch `working_branch`: `blackjax/adaptation/adjusted_alba.py`, `unadjusted_alba.py` | The paper's tuning scheme (not on BlackJAX main) | Raw source, same method. |
| github.com/reubenharry/sampler-benchmarks: `sampler-comparison/sampler_comparison/samplers/microcanonicalmontecarlo/adjusted.py`, `mchmc.py`, `general.py` | The paper's benchmark driver | Raw source, same method. |

Notation below: target p(x) = exp(-L(x)) on R^d, L = -log p, g(x) = grad L(x). "(inference)" marks my own
derivation or synthesis. Where the extraction tool paraphrased rather than quoted, I say "as extracted".

## 1. Dynamics

**Continuous dynamics (2503.01707 Eq. 3; 2212.08549 Eq. 12; 2303.18221 Eq. 1).** Position x in R^d and
unit velocity u on S^(d-1):

    dx/dt = u
    du/dt = -(I - u u^T) g(x) / (d - 1)

(I - u u^T) is the projection onto the tangent plane of the sphere, so |u| = 1 is preserved. 2503.01707
uses the factor 1/(d-1); 2212.08549 Eq. 12 writes the same equations with 1/d (footnote 1 of 2303.18221,
as extracted, says d-1 is used instead of d). For the MH-adjusted algorithm the d-1 convention is the
one consistent with the energy formula in Section 3, so use d-1 throughout.

**Stationary distribution (2503.01707 Sec. B.4, Eq. 29).** The exact flow preserves p(x) U_{S^(d-1)}(u),
so the x-marginal is the target. In 2212.08549 Eq. 8 the underlying Hamiltonian is
H = (d/2) log(Pi^2/d) + L(x) (the q = 0 "variable mass" case); 2503.01707 Sec. B.2 writes the kinetic
term as T(Pi) = (d-1) log|Pi|. The weight w = |Pi|/d of 2212.08549 Eq. 13 does not appear in the MAMS
state; its logarithmic change becomes the per-step kinetic energy change of Section 3.

**Exact velocity map under a frozen force (2503.01707 Eq. 28; 2212.08549 Eq. 16).** With

    e = -g(x)/|g(x)|  (unit vector along grad log p),   delta = eps |g(x)| / (d - 1),

the solution of du/dt = -(I - u u^T) g/(d-1) over a time eps at fixed x is

    B_eps(u) = [ u + ( sinh(delta) + (e.u)(cosh(delta) - 1) ) e ] / [ cosh(delta) + (e.u) sinh(delta) ].

2503.01707 Eqs. 26 and 27 derive it by the substitution u = y'/(g_bar . y) with g_bar = -g/(d-1), which
linearises the equation to y'' = (g_bar g_bar^T) y.

Numerically stable form used by BlackJAX (`esh_dynamics_momentum_update_one_step`, verified algebraically
equal to Eq. 28 by me, (inference)): with zeta = exp(-delta) and c = e.u,

    u_raw = e (1 - zeta) (1 + zeta + c (1 - zeta)) + 2 zeta u,     u_new = u_raw / |u_raw|.

**Partial momentum refreshment (2503.01707 Eq. 11; 2212.08549 Eqs. 28 to 30).** 2503.01707 Eq. 11:

    O_eps(u) = (c1 u + c2 Z/sqrt(d)) / |c1 u + c2 Z/sqrt(d)|,   c1 = exp(-eps/L_partial),  c2 = sqrt(1 - c1^2),

Z ~ N(0, I_d). 2212.08549 Eq. 28 and 29 write the same step as u <- (u + nu z)/|u + nu z| with
nu = sqrt((exp(2 eps/L) - 1)/d); the two forms are the same map up to the normalisation (inference:
c2/(c1 sqrt(d)) = nu). 2212.08549 Eq. 30 gives the meaning of L: <u_n . u_0> = exp(-n eps / L), so L
is the decoherence length of the velocity in units of integration time (which equals path length, since
|u| = 1). In MAMS the parameter L_partial controls this refreshment, and 2503.01707 Appendix D (Eq. 51)
sets L_partial / L = 1.25, where L is the proposal trajectory length. The default MAMS of 2503.01707
has no partial refreshment (full refreshment of u once per MH step); the "MAMS-Langevin" variant adds
O steps and is Algorithm 1 of the paper.

## 2. Integrator

**Leapfrog, the one used in 2503.01707 (Eq. 5):**

    Phi_eps = B_{eps/2} o A_eps o B_{eps/2},     A_eps(x, u) = (x + eps u, u).

Read right to left: half velocity kick, full position drift, half velocity kick. The paper calls it
"leapfrog (or velocity Verlet)" (Sec. 4 as extracted). One gradient per step once the gradient at the
current x is cached (the trailing B_{eps/2} of step i and the leading B_{eps/2} of step i+1 use the same
gradient) (inference). 2212.08549 Eq. 14 is the same scheme, and notes it is not symplectic in these
variables.

**Minimal norm (McLachlan / Omelyan two-stage), 2212.08549 Eq. 17:**

    Phi_eps = B_{lambda eps} o A_{eps/2} o B_{(1-2 lambda) eps} o A_{eps/2} o B_{lambda eps},
    lambda = 0.19318...  (BlackJAX: b1 = 0.1931833275037836, a1 = 0.5, b2 = 1 - 2 b1).

Two gradients per step. 2212.08549 Table 1 shows it as slightly better than leapfrog for the unadjusted
sampler; 2303.18221 Sec. 7 uses it for the phi^4 runs. 2503.01707 does not mention the minimal norm
integrator in the text I could read; its benchmarks use leapfrog, and the authors' benchmark driver
(`adjusted.py`, `adjusted_mclmc(integrator_type="velocity_verlet", ...)`) confirms this. BlackJAX's
`adjusted_mclmc` and `adjusted_mclmc_dynamic` kernels nevertheless default to `isokinetic_mclachlan`.

**Code reading of the maps (per sub-step with coefficient c, step size eps):**

    # velocity kick, returns updated u and the kinetic energy change dK (Section 3)
    gnorm = |g|; e = -g/gnorm; delta = c*eps*gnorm/(d-1); zeta = exp(-delta); cu = e.u
    u = e*(1-zeta)*(1+zeta+cu*(1-zeta)) + 2*zeta*u;  u /= |u|
    dK = (d-1)*(delta - log(2) + log((1+zeta**2) + cu*(1-zeta**2)))
    # position drift
    x = x + c*eps*u

With a diagonal preconditioner sigma (Section 4) the drift becomes x += c eps sigma*u and the kick uses
g_y = sigma*g(x) in place of g (this is what BlackJAX does via `forward_L` and `adjoint_L`).

## 3. The Metropolis adjustment

**Proposal.** One MH step integrates n steps of Phi_eps from z_0 = (x_0, u_0) with u_0 freshly drawn
from U_{S^(d-1)} (2503.01707 Sec. 4, Eq. 4: phi = T o Phi_{t/n}^n with t = n eps, T(x,u) = (x,-u)). The
time reversal T makes phi an involution, phi o phi = id, which is what makes the proposal symmetric in
the MH sense (deterministic involutive proposal). The paper states T "does not affect the dynamics in
any way, because a full velocity refreshment is performed after the Metropolis step" (Sec. 4 as
extracted), so in code u is simply discarded after the accept/reject.

**Accept-reject quantity.** Because the flow is not volume preserving in (x, u), the MH ratio needs a
Jacobian. Lemma 4.1 (2503.01707):

    W(z', z) = -log[ p(z')/p(z) ] - int_0^T div F(z(s)) ds,

with F the drift field, and the acceptance probability is min(1, exp(-W)). Lemma 4.2 makes W additive
over sub-steps:

    position drift A_eps:    W_A = V(x') - V(x) = -log p(x') + log p(x) = L(x') - L(x)         (Eq. 6)
    velocity kick B_eps:     W_B = K(B_eps z) - K(z) = (d-1) log[ cosh(delta) + (e.u) sinh(delta) ]  (Eq. 8)

where u in Eq. 8 is the velocity before the kick, and delta, e as in Section 1. Eq. 8 is the exact
log-Jacobian of the sphere map B_eps plus the (zero) density change of the uniform velocity marginal;
Appendix B.6 (Eqs. 31 to 37) derives it directly in spherical coordinates, and Appendix B.7 (Eqs. 38 to
42) via the divergence div F_B = -|g| (u.e) (the A step has zero divergence). Appendix B.5 shows W equals
the energy change Delta E = Delta V + Delta K of a rescaled Hamiltonian, which is why the paper calls it
an energy error; for HMC the analogous term is Eq. 7, (1/2)|u'|^2 - (1/2)|u|^2.

For a trajectory of n leapfrog steps, the A terms telescope, so (inference)

    W = L(x_n) - L(x_0) + sum over all 2n half-kicks of (d-1) log[cosh(delta_i) + (e_i . u_i) sinh(delta_i)],

accept with probability min(1, exp(-W)). This is exactly what BlackJAX computes:
`delta_energy = -logp(x_0) + logp(x_n) - sum(dK)` = -W, accept with `min(1, exp(delta_energy))`.

**Random trajectory length (2503.01707 Sec. 5).** The k-th proposal uses n_k = ceil(2 h_k L/eps) steps,
h_k ~ U(0,1) or the k-th Halton element. Footnote 5 corrects for the ceiling so that E[n_k] = L/eps.
I could not read footnote 5 verbatim; the extracted content was: Y = floor(2L/eps - 1), E[n_k] =
(Y+1)(y - Y/2)/y, and y = Y(Y+1)/(Y+1 - L/eps), with n_k = ceil(2 y h_k L/eps). That last pair is not
self-consistent with the stated expectation. My own derivation (inference): for n = ceil(h s), h ~ U(0,1),
Y = floor(s), E[n] = (Y+1)(s - Y/2)/s; solving E[n] = mu = L/eps gives s = Y(Y+1)/(2(Y+1-mu)), and
n_k = ceil(h_k s). This is precisely BlackJAX `rescale(mu)`: `k = floor(2 mu - 1); x = k (mu - 0.5(k+1))
/(k+1-mu); return k + x` (which equals k(k+1)/(2(k+1-mu))), with `num_steps = ceil(U * rescale(mu))`
clipped to at least 1. Use the code form.

**After rejection.** The chain keeps x_0 (and its cached logp and gradient); u is not part of the chain
state and a new u ~ U_{S^(d-1)} is drawn at the start of the next transition (Algorithm 1 as extracted:
"On rejection x^{I+1} stays at the start z^0[0]"; Theorem 5.2). There is therefore no momentum flip and no
persistence of velocity across MH steps in the default MAMS. The Langevin variant (Algorithm 1, scheme
T o O B A B O) also redraws u at the start of each MH step (as extracted); the O steps act only within
the trajectory and do not enter W (Theorem 5.1: acceptance min(1, exp(-Delta)) with Delta the energy
error of the deterministic part only). Whether the paper's O within OBABO carries parameter eps or eps/2
was not resolvable from the extraction; BlackJAX `with_isokinetic_maruyama` applies the refresh with
`step_size * 0.5` before and after each deterministic step, which is the natural Strang splitting.

## 4. Tuning

All from 2503.01707 Sec. 6 unless noted.

**Order and budget.** Three stages: step size, then preconditioning matrix, then trajectory length, each
given 10 percent of the total sampling budget (citing 2412.08876 for the 10 percent rule). Sec. 7: "We take
the adaptation steps as our burn-in" and sampling starts from the final adapted state.

**Step size eps.** Dual averaging (Nesterov 2009, as in Hoffman and Gelman 2014) on the MH acceptance
rate. Appendix C derives the optimum: with W ~ N(mu, sigma^2) and the Jarzynski identity E[exp(-W)] = 1
forcing sigma^2 = 2 mu (Eq. 43), P(accept) = 2 Phi(-sqrt(mu/2)) (Eq. 44); ESS per gradient scales as
mu^(1/4) Phi(-sqrt(mu/2)) d^(-1/4) (Eq. 48), maximised at mu = 0.41, i.e. 65 percent acceptance, with
eps proportional to d^(-1/4) ... (Eq. 47: sigma^2/d proportional to eps^4/d^2). Despite this the paper
sets the target to 90 percent ("Empirically ... 90 percent works well" for the second-order integrator;
80 percent would be optimal for a fourth-order integrator, mu = 0.13). The benchmark code passes
`target_acc_rate=0.9`. For the unadjusted MCLMC the step size is instead set from the energy error
variance per dimension (2212.08549 Sec. 2.6.1: target Var[E]/d around 0.0005, update
eps = eps_0 (0.0005 d / Var[E])^(1/4); the text also says 0.0003 and Fig. 5 says 0.001 is near optimal;
BlackJAX `mclmc_find_L_and_step_size(desired_energy_var=5e-4)`).

**Diagonal preconditioning.** "A simple choice of diagonal preconditioning matrix is obtained by
estimating variance along each parameter" (Sec. 6 as extracted); the paper obtains the estimates from a
sampler run (a NUTS run in the grid-search check of Sec. 7). The paper does not write out how the matrix
enters the dynamics. BlackJAX (`generate_isokinetic_integrator`, `inverse_mass_matrix` diagonal): the
drift is x += eps c sqrt(invM) * u and the kick uses sqrt(invM) * grad logp(x). With sqrt(invM) = sigma
(per-parameter standard deviation) this is the MAMS dynamics applied to y = x / sigma (inference). In the
authors' ALBA code the diagonal is the per-dimension sample variance accumulated in the unadjusted
warm-up windows (`mass_matrix_adaptation(is_mass_matrix_diagonal=True)`).

**Trajectory length L.** Initial guess L = sqrt(d) (Sec. 6: after preconditioning the per-dimension
variances are 1, and L is set to sqrt of the sum of eigenvalues of the covariance; BlackJAX
`unadjusted_alba`: `L = sqrt(d)/v`, and `sqrt(sum(inverse_mass_matrix))` when not preconditioning).
Refined rule (Eqs. 12 to 14, "ALBA"): L_ALBA proportional to the expected time between effective
samples, L tau_int, with tau_int[f] = 1 + 2 sum_t rho_t[f] the integrated autocorrelation time of
f(x_i) = x_i, harmonically averaged over i; constant 0.3 for MAMS without Langevin noise and 0.23 with
it, calibrated on standard Gaussians. In the authors' code the same quantity is computed from the
unadjusted MCLMC warm-up: `L = alba_factor * step_size * mean(num_alba_steps / ess)` with
`ess = blackjax.diagnostics.effective_sample_size(samples)`, then clipped to at most 500 step sizes; the
benchmark driver passes `alba_factor=0.3` for MAMS (0.4 is the code default and is used for the grid
search and the unadjusted sampler; 2212.08549 Sec. 2.6.2 uses L = 0.4 l with l from Eq. 33). Note the
code's arithmetic mean of N/ess_i equals N over the harmonic mean of ESS, while the paper's harmonic mean
of tau_i equals N over the arithmetic mean of ESS; these differ (inference). L is then held fixed while
the adjusted step size is dual averaged (`da_adaptation`), and the per-proposal step count is
`ceil(U * rescale(L/eps))` recomputed from the current eps at every iteration.

**L_partial (Langevin variant).** Appendix D, Fig. 4 and Eqs. 49 to 51: optimum for ill-conditioned
Gaussians at damping beta = 0.567 and T = 1.413, so (1/beta)/T = 1.25 and L_partial = 1.25 L. In
BlackJAX this is `L_proposal_factor=1.25` (the kernel multiplies it by n eps); the benchmark driver's
default `L_proposal_factor=jnp.inf` is the no-Langevin MAMS.

## 5. Reported efficiency

Metric (2503.01707 Sec. 7): b^2(f) = (E_sampler[f] - E[f])^2 / Var[f]; b_max^2 = max over i of
b^2(x_i^2) (average over i for product targets; entropy for Cauchy). Report the number of gradient calls
n to reach median-over-128-chains b_max^2 < 0.01, "an accuracy equivalent of 100 effective samples", so
n = 100 / (ESS per gradient). NUTS: dual averaging to 80 percent acceptance. Tuning steps are the
burn-in; whether their gradients are included in n is not stated in the text I read.

Table 1 of 2503.01707 (gradient calls to b_max^2 < 0.01, lower is better; dimensions from Appendix E.1):

| Problem (d) | NUTS | MAMS | MAMS-Langevin | MAMS grid search |
|---|---|---|---|---|
| Gaussian, kappa = 100 (100) | 19652 | 3249 | 3172 | 3121 |
| Banana (2) | 95519 | 14078 | 14818 | 15288 |
| Bimodal (50) | 210758 | 139418 | 136770 | 123295 |
| Rosenbrock (36) | 161359 | 94184 | 103545 | 93782 |
| Cauchy (100) | 171373 | 110404 | 155963 | 87900 |
| Brownian (32) | 29816 | 13528 | 15232 | 14015 |
| German credit (51) | 88975 | 55748 | 49979 | 52265 |
| Item response (501) | 76043 | 45371 | 56902 | 45640 |
| Stochastic volatility (2429) | 843768 | 430088 | 510190 | 431957 |
| Neal's funnel (20) | > 1e8 | 2346899 | 1765311 | 1013048 |

In ESS-per-gradient terms (100/n, inference): Gaussian 0.031 vs 0.0051; stochastic volatility 0.00023
vs 0.00012; item response 0.0022 vs 0.0013. The text says MAMS beats NUTS "by a factor of two at worst and
seven at best"; from the table the worst ratios are about 1.5 (Bimodal) and 1.6 (Cauchy, German credit).
No problem in Table 1 has MAMS worse than NUTS. MAMS-Langevin is slower than plain MAMS on most problems
(Cauchy notably: 155963 vs 110404). Appendix E.2 Table 2 gives bootstrap uncertainties, mostly below
2 percent (Bimodal 6.3 percent, Cauchy 5.4 percent for MAMS). Fig. 2: both samplers scale as
kappa^(-1/2) on 100-d Gaussians, MAMS about 4x faster. Fig. 3: both scale as d^(1/4) in gradients per
effective sample on Gaussian and Rosenbrock, MAMS with the better constant. Unadjusted MCLMC numbers are
not in 2503.01707 beyond Fig. 1 (bias illustration). For reference, 2212.08549 Table 1 (ESS = 200/n at
b_2 = 0.1, tuning excluded) gives unadjusted MCLMC (MN integrator) vs NUTS: ill-conditioned Gaussian
0.110 vs 0.012, bimodal 0.064 vs 0.008, Rosenbrock 0.0033 vs 0.0015, funnel 0.021 vs 0.006, German credit
0.0099 vs 0.0014, stochastic volatility 0.023 vs 0.006; Table 2 (tuning included, leapfrog) 0.075 vs
0.006 on the Gaussian down to 0.011 vs 0.001 on stochastic volatility.

## 6. Minimal pseudocode of one MAMS transition

Consistent with Sections 1 to 4; leapfrog; diagonal preconditioner sigma (set sigma = 1 to disable);
set L_partial = None for the default MAMS or 1.25 L for MAMS-Langevin.

    def mams_transition(x, lp, g, eps, L, sigma, rng, L_partial=None):
        d = x.size
        u = rng.standard_normal(d); u /= norm(u)                      # full refresh, Sec. 4
        n = max(1, ceil(rng.uniform() * rescale(L / eps)))            # Sec. 5, footnote 5 via BlackJAX rescale
        x0, lp0, g0 = x, lp, g; dK = 0.0
        for i in range(n):
            if L_partial: u = refresh(u, eps/2, L_partial, rng)       # O_{eps/2}, Eq. 11 (not in W)
            u, dk = kick(u, sigma*g, eps/2, d); dK += dk              # B_{eps/2}, Eq. 28 and Eq. 8
            x = x + eps * sigma * u                                   # A_eps
            lp, g = logp_and_grad(x); g = -g  # g = grad L = -grad logp
            u, dk = kick(u, sigma*g, eps/2, d); dK += dk              # B_{eps/2}
            if L_partial: u = refresh(u, eps/2, L_partial, rng)       # O_{eps/2}
        W = (-lp) - (-lp0) + dK                                       # Lemma 4.2; A terms telescope
        if log(rng.uniform()) < -W: return x, lp, g, True             # accept, min(1, e^{-W})
        return x0, lp0, g0, False                                     # reject: keep x0, u discarded

    def kick(u, g, h, d):           # g = grad L in preconditioned coordinates, h = sub-step size
        gn = norm(g); e = -g/gn; delta = h*gn/(d-1); z = exp(-delta); c = e @ u
        u = e*(1-z)*(1+z+c*(1-z)) + 2*z*u; u /= norm(u)
        return u, (d-1)*(delta - log(2) + log((1+z*z) + c*(1-z*z)))

    def refresh(u, h, Lp, rng):     # 2212.08549 Eq. 28 and 29 form of 2503.01707 Eq. 11
        nu = sqrt((exp(2*h/Lp) - 1)/u.size); u = u + nu*rng.standard_normal(u.size); return u/norm(u)

Cost: exactly n gradient evaluations per transition with leapfrog (the gradient at x0 is cached, and the
accepted state carries the gradient at x_n) (inference). For the minimal norm integrator replace the
loop body by the five-map sequence of Section 2 (2 gradients per step) and keep summing dk.

Ambiguities flagged: (a) footnote 5 formula (resolved by derivation and the BlackJAX code, Section 3);
(b) whether O carries eps or eps/2 in the paper's OBABO (BlackJAX uses eps/2); (c) the sign convention
of "delta" in Algorithm 1 (the extraction says the accept test is U < exp(-delta) with delta the summed
energy error, consistent with W above); (d) whether tuning gradients count toward n in Table 1; (e) the
paper's harmonic averaging of tau_int versus the code's arithmetic mean of N/ESS (Section 4).

## BlackJAX implementation versus the paper

Kernels on BlackJAX `main`:

- `blackjax.mcmc.adjusted_mclmc_dynamic.build_kernel(integrator=isokinetic_mclachlan, divergence_threshold=1000, next_random_arg_fn=..., integration_steps_fn=lambda key: randint(key, (), 1, 10))`; inner kernel `(rng_key, state, logdensity_fn, step_size, L_proposal_factor=jnp.inf, inverse_mass_matrix=1.0, integration_steps_params=())`. Module docstring: "Public API for the Metropolis Hastings Microcanonical Hamiltonian Monte Carlo (MHMCHMC) Kernel"; the kernel is recommended over the static one. Helpers `make_random_trajectory_length_fn(bool)` (ceil(U * rescale(avg)) or ceil(avg), min 1) and `trajectory_length(t, mu) = rint(0.5 + halton_sequence(t) * rescale(mu))`. `rescale` and `halton_sequence` as in Section 3.
- `blackjax.mcmc.adjusted_mclmc` (static): `build_kernel(integrator=isokinetic_mclachlan, divergence_threshold=1000)`, `as_top_level_api(logdensity_fn, step_size, L_proposal_factor=jnp.inf, inverse_mass_matrix=1.0, *, divergence_threshold=1000, integrator=isokinetic_mclachlan, num_integration_steps=None, integration_steps_params=None)`. `adjusted_mclmc_proposal` runs `fori_loop` over the steps, accumulates the per-step `kinetic_energy_change`, sets `delta_energy = -logp(x0) + logp(x_end) - kinetic_energy` (NaN to -inf; divergent if `-delta_energy > divergence_threshold`) and calls `static_binomial_sampling`, which accepts with `p_accept = min(1, exp(delta_energy))` and returns the previous state on rejection. A fresh unit-vector momentum is drawn every call (`generate_unit_vector`); the state holds no momentum. So the accept test and the rejection behaviour match 2503.01707 (delta_energy = -W).
- `blackjax.mcmc.integrators`: `esh_dynamics_momentum_update_one_step` implements the Section 2 kick with `delta = step_size*coef*gradient_norm/(dims-1)` and `kinetic_energy_change = (delta - log(2) + log_a)*(dims-1)`, `log_a = log(s_plus(1-zeta^2) + 2 zeta^2)` where `s_plus = 0.5|u + e|^2 = 1 + e.u`; I checked this equals Eq. 8 (inference). Coefficient lists: `velocity_verlet [0.5, 1.0, 0.5]`, `mclachlan [b1, 0.5, 1-2b1, 0.5, b1]` with `b1 = 0.1931833275037836`, plus `yoshida` and `omelyan`; even entries are kicks, odd entries drifts, so the first and last maps are kicks as in Eq. 5 and 2212.08549 Eq. 17. Preconditioning via `inverse_mass_matrix` (diagonal or `LowRankInverseMassMatrix`): gradient through `adjoint_L` (multiply by sqrt(invM)), drift by `forward_L(u)`. `with_isokinetic_maruyama` wraps a step as `partially_refresh_momentum(step_size*0.5, L)`, deterministic step, refresh again, with `nu = sqrt((exp(2 h/L) - 1)/dim)` and no refresh when L is infinite; the kernels pass `L = L_proposal_factor * num_integration_steps * step_size`, so the paper's L_partial = 1.25 L is `L_proposal_factor = 1.25` (used in `blackjax/adaptation/laps.py`).
- Differences from the paper's experiments: default integrator mclachlan (paper: leapfrog); default `integration_steps_fn` is uniform on 1 to 9 (paper: ceil(2 h L/eps) with the footnote-5 correction; the paper's driver passes `make_random_trajectory_length_fn(True)` with `L/step_size`).

Adaptation:

- `blackjax.adaptation.adjusted_mclmc_adaptation.adjusted_mclmc_find_L_and_step_size(..., target, frac_tune1=0.1, frac_tune2=0.1, frac_tune3=0.0, diagonal_preconditioning=True, max="avg", tuning_factor=1.3, target_num_integration_steps=2.0)` on `main` is not the paper's scheme: it pins L to `target_num_integration_steps * step_size`, rescales by sqrt(sum of variances) * 1.3, and its optional stage 3 uses `L = 0.5 L / mean(ESS per sample)`. It has no 0.3 or 0.23 constant and no Halton sequence.
- The paper's scheme is `adjusted_alba` in `reubenharry/blackjax`, branch `working_branch`, `blackjax/adaptation/adjusted_alba.py` (imported by the paper's driver as `blackjax.adaptation.adjusted_alba.adjusted_alba`; it is not on BlackJAX `main` nor on the fork's `main`). Signature: `adjusted_alba(unadjusted_mcmc_kernel, unadjusted_init, logdensity_fn, target_eevpd, v, adjusted_mcmc_kernel, adjusted_init, integrator, target_acceptance_rate=0.80, num_alba_steps=500, alba_factor=0.4, preconditioning=True, L_proposal_factor=jnp.inf, **extra)`, `run(rng_key, position, num_steps=1000, num_unadjusted_steps=20000)`. It runs `unadjusted_alba` (unadjusted MCLMC warm-up: initial step size sqrt(d)/5, initial L = sqrt(d)/v, window schedule from `window_adaptation.build_schedule`, energy-variance step size via `robnik_step_size_tuning(desired_energy_var=target_eevpd)`, diagonal variances via `mass_matrix_adaptation`, L reset to sqrt(d)/v after each slow window when preconditioning, then `num_alba_steps` steps to compute `L = alba_factor*step_size*mean(num_alba_steps/ess)` clipped to 500 step sizes), then `da_adaptation` (dual averaging of the adjusted step size at `target_acceptance_rate`, L held fixed, step count ceil(U * rescale(L/eps))). The paper's driver `adjusted.py`: `adjusted_mclmc(integrator_type="velocity_verlet", diagonal_preconditioning=True, L_proposal_factor=jnp.inf, target_acc_rate=0.9, random_trajectory_length=True, num_adjusted_tuning_steps=20000, num_unadjusted_tuning_steps=20000, alba_factor=0.3, num_alba_steps=None (-> num_unadjusted_tuning_steps//3), target_eevpd=5e-4)`. The `v` argument is undocumented in the code I read, and the Langevin constant 0.23 does not appear in the code (only 0.3 and 0.4). The 10 percent per stage of Sec. 6 is not visible as such in this code path (fixed 20000-step budgets are used instead).

## Checks against nextgen/LITERATURE.md and nextgen/research/

- `nextgen/LITERATURE.md` line 55 ("MCLMC; MAMS (adjusted), up to 7x over NUTS") is consistent with Table 1 of 2503.01707.
- `nextgen/research/04_adjacent_field_methods.md` line 9 attributes "on hierarchical models it needed about 3x fewer gradients than NUTS" to 2303.18221. That paper only benchmarks the phi^4 lattice (12x and 32x speedups, abstract). The hierarchical-model comparisons (funnel 0.021 vs 0.006, stochastic volatility 0.023 vs 0.006, German credit 0.0099 vs 0.0014) are in 2212.08549 Table 1, and the ratios are 3.5x to 7x, not 3x. Worth correcting.

## Gaps and suggested follow-ups

Looked for and not found:
- A verbatim reading of footnote 5 and of Algorithm 1 of 2503.01707 (only tool paraphrases; the PDF was not fetched). Follow-up: read `arxiv.org/pdf/2503.01707v2` Sec. 5 and footnote 5 directly.
- An explicit statement in 2503.01707 of how the diagonal preconditioner enters the dynamics; taken from the BlackJAX code instead.
- Whether tuning gradients are counted in Table 1 of 2503.01707.
- The 0.23 Langevin ALBA constant in any code.
- The meaning of `v` in `unadjusted_alba` (searched the fork's `working_branch` file only).
- An `adjusted_alba` module on BlackJAX `main` (directory listing checked: absent), on `reubenharry/blackjax` `main`, `current`, `benchmarks`, `laps` (404); present on `working_branch` only.
- Searches run: web search "blackjax adaptation alba adjusted mclmc Halton L_proposal_factor github"; "github reubenharry sampler-benchmarks OR sampler-comparison alba adjusted mclmc tuning"; "alba adaptation adjusted_mclmc mclmc autocorrelation 0.3 trajectory length blackjax"; GitHub API listings of `blackjax-devs/blackjax/blackjax/adaptation`, `reubenharry/blackjax` branches, and `reubenharry/sampler-benchmarks` tree (first 200k characters of a 970k-character listing).

Suggested next searches: INSPIRE or arXiv citation list of 2503.01707 for follow-up work on MAMS tuning (for example 2412.08876 v3 of 2026-06 and any 2026 ALBA paper); the PDF of 2503.01707 for footnote 5 and Algorithm 1; `blackjax/adaptation/step_size.py` for `dual_averaging_adaptation` defaults to mirror in NumPy.
