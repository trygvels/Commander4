# indexfield prototype: interface specification (M2.1, E012 v2)

Pure NumPy/SciPy/healpy package under `nextgen/prototypes/indexfield/`. The only Commander4 imports
allowed are `commander4.math_utils.nuts.NUTSSampler` (sampler) and, in one optional check, the C++
kernel `commander4._cmdr4_backend...marginal_loglike_regions` for a cross-check. Everything else is
self-contained so the target is auditable. Run with `/Users/trygvels/code/Commander4/.venv/bin/python`.
Style: 100-character lines, docstrings state the mathematics, no em-dashes.

The mathematics is fixed by `nextgen/experiments/E012_gmrf_index_field_target.md` (v2). This file
fixes names, shapes and units so that modules written by different people fit together. Where this
file and E012 disagree, E012 wins and this file is corrected.

Units: frequencies in GHz, amplitudes and data in uK_RJ, angles in radians.

## Shapes and conventions
- `B` bands, `C = 3` components in the order (cmb, sync, dust), `N` data pixels (nside_d, RING
  order throughout), `Nb = N_beta` index pixels (nside_beta, RING), `K` basis modes.
- `region[p]` in `[0, Nb)`: index pixel of data pixel p. Computed as
  `hp.ring2nest(nside_d, p) >> (2 * log2(nside_d / nside_beta))` then `hp.nest2ring(nside_beta, .)`.
- A "field" is `beta` in `R^Nb`; `beta_data = beta[region]` in `R^N`.

## model.py
```python
def cmb_to_rj(nu_ghz) -> float | ndarray            # uK_CMB -> uK_RJ factor x^2 e^x/(e^x-1)^2
def mixing(nu_ghz, beta_s, beta_d=1.6, T_d=20.0, nu0_s=30.0, nu0_d=545.0) -> (M, dM_dbeta_s)
    # M: (B, C) for scalar beta_s, dM_dbeta_s: (B,) derivative of the sync column
    # columns: cmb = cmb_to_rj(nu), sync = (nu/nu0_s)**beta_s, dust = MBB in RJ normalized at nu0_d:
    #   (nu/nu0_d)**(beta_d + 1) * (exp(h nu0_d / k T_d) - 1) / (exp(h nu / k T_d) - 1)

@dataclass
class PixelData:
    nu: ndarray (B,)
    d: ndarray (N, B)                # data, uK_RJ
    inv_var: ndarray (N, B)          # 1/sigma^2, uK_RJ^-2
    region: ndarray (N,) int
    n_beta: int
    beta_d: float; T_d: float
    m: ndarray (N, C) | None         # Gaussian amplitude prior mean (mode "gauss")
    inv_T: ndarray (N, C) | None     # Gaussian amplitude prior precision, diagonal (mode "gauss")

def collapsed_loglike(beta, data, mode) -> (logL: float, grad: ndarray (Nb,))
    # mode in {"gauss", "w0", "w1"}; beta: (Nb,)
    # Per pixel with A = M^T N^-1 M (C x C), b = M^T N^-1 d:
    #  "w1":    -1/2 (d^T N^-1 d - b^T A^-1 b) - 1/2 log det A
    #  "w0":    -1/2 (d^T N^-1 d - b^T A^-1 b)
    #  "gauss": with A' = A + inv_T, b' = b + inv_T m:  1/2 b'^T A'^-1 b' - 1/2 log det A'
    #           (constant terms -1/2 d^T N^-1 d dropped in all modes is acceptable but must be
    #            the same choice for logL used in checks against grids: state it in the docstring)
    # grad per index pixel r = sum over p in r of the pixel gradient (E012 "Gradient"), vectorized:
    #   group pixels by region with np.add.at or bincount.
    # Must raise if any A (or A') is not positive definite (np.linalg.cholesky on the stack).
```
Vectorization: stack `M` as `(N, B, C)` from `mixing` evaluated at `beta[region]` (the sync column
depends on the pixel; cmb and dust columns are shared), use `np.linalg.solve` / `np.linalg.cholesky`
on `(N, C, C)` stacks. Target: one evaluation at N = 3072, B = 12 under 20 ms.

## field.py
```python
class Basis:                   # abstract: V (Nb, K) orthonormal columns, lam (K,) >= 0, h (float)
class GraphLaplacianBasis(Basis):
    def __init__(self, nside_beta):
        # 8-neighbour HEALPix graph (hp.get_all_neighbours, drop -1 entries, symmetrize),
        # G = D - W, eigh(G); DROP the eigenvector with eigenvalue 0 (constant mode);
        # h = sqrt(4 pi / Nb); lam = eigenvalues / (3 h^2)   (D009: continuum Laplacian scale)
class HarmonicBasis(Basis):
    def __init__(self, nside_beta, lmax=None):   # default lmax = 2 nside_beta
        # real Y_lm (1 <= l <= lmax) evaluated at pixel centres, orthonormalized by QR (columns
        # ordered by l, then m); lam = l(l+1) of the column's l; K = (lmax+1)^2 - 1

class MaternField:
    def __init__(self, basis, alpha=2.0): ...
    def spectrum(self, rho) -> (s: ndarray (K,), c: float, ds_dlogrho: ndarray (K,), dc_dlogrho: float)
        # kappa = sqrt(8 (alpha - 1)) / rho  [for alpha = 2: sqrt(8)/rho], rho in radians
        # s_i = (kappa^2 + lam_i)^(-alpha/2) ... NOTE E012 v2 writes s_i = (kappa^2 + lam_i)^(-1)
        #   for alpha = 2; use the general form (kappa^2 + lam_i)^(-alpha/2) and default alpha = 2.
        # c^2 = (1/Nb) sum_i s_i^2
    def S_times(self, rho, z) -> ndarray (Nb,)           # V diag(s/c) z
    def ST_times(self, rho, g) -> ndarray (K,)           # diag(s/c) V^T g
    def dS_dlogrho_times(self, rho, z) -> ndarray (Nb,)  # V diag(d(s/c)/dlogrho) z
    def logdet_K(self, rho) -> float                     # sum_i 2 log(s_i / c)   (centred param.)
    def Kinv_quadform(self, rho, v) -> float             # v^T V diag((c/s)^2) V^T v (centred param.)

class RhoTransform:
    def __init__(self, rho_min, rho_max): ...
    def rho(self, u) -> float; def u(self, rho) -> float
    def log_jacobian(self, u) -> float     # log |d log rho / d u| = log[(L) s(u)(1-s(u))], L = log(rho_max/rho_min)
    def dlogrho_du(self, u) -> float
```

## posterior.py
```python
@dataclass
class HyperPrior:
    mu0: float = -3.0; s_mu: float = 0.5
    sigma_scale: float = 0.3                 # half-normal scale on sigma
    rho_prior: str = "loguniform"            # or "pc"
    rho_min: float; rho_max: float           # radians (h and pi by default)
    pc_rho0: float | None; pc_sigma0: float | None   # PC prior settings (E012)

class NonCentredPosterior:
    """x = (z[K], mu, log_sigma, u). logp_and_grad(x) -> (float, ndarray (K+3,)) of E012 'Full
    posterior sampled' with the stated Jacobians. Also beta(x) -> (Nb,), theta(x) -> (mu, sigma, rho)."""
class CentredPosterior:
    """x = (beta[Nb], mu, log_sigma, u) with prior N(beta; mu 1, sigma^2 K(rho)) restricted to the
    K-dimensional range of V: log prior = -1/2 Kinv_quadform(beta - mu 1)/sigma^2 - K log sigma
    - 1/2 logdet_K(rho) (+ the mu, sigma, u terms as above). The component of beta - mu 1 outside
    range(V) has zero prior variance: enforce by parameterizing beta = mu 1 + V y, y in R^K, i.e.
    x = (y[K], mu, log_sigma, u). State this in the docstring."""
```
Both expose `dim`, `logp_and_grad(x)`, `unpack(x)` and a `random_init(rng)`.

## sampler.py
```python
def run_nuts(post, rng, n_warmup, n_samples, x0=None, target_accept=0.8) -> Result
# Result: samples (n_samples, dim), logp (n_samples,), n_divergent, step_size, ebfmi, ess (dim,)
# ESS by initial positive sequence / Geyer on each coordinate (own implementation, no arviz).
```
Uses `commander4.math_utils.nuts.NUTSSampler` (read `src/commander4/math_utils/nuts.py` for its
API: `transition(theta, logp_and_grad) -> (theta, TransitionInfo)`, adaptation handled inside
during the first `n_warmup` transitions).

## simulate.py
```python
def templates(nside_d, rng) -> m: ndarray (N, C)      # cmb 0; sync, dust plane-brightening as
                                                        # aux/validate_mcmc_likelihood.py (10**(0.3+2.0 t), 10**(0.5+1.5 t))
def noise_rms(nside_d, nu, scale=1.0) -> (N, B)        # E006 RMS list extended to the 12 bands (state values)
def simulate(nside_d, nside_beta, field, hyper, rng, amp_mode, noise_scale, beta_d, T_d)
    -> (data: PixelData, truth: dict(theta, z, beta, a))
# amp_mode "gauss": a ~ N(m, T) with T^(1/2) = 10 * max(|m|, 10) per component (floor 100 uK_RJ)
# amp_mode "template": a = m * exp(0.1 * N(0,1)) for sync, dust; cmb ~ N(0, 70^2)
```

## reference.py  (INDEPENDENT code path: do not import model.py, field.py or posterior.py)
```python
def loglike_pixel_loops(beta_field_data, data, mode) -> float
    # plain loops / small matrices, no shared code with model.py; same definition of dropped constants
def grid_1d(data, mode, beta_grid, theta_fixed=None) -> (logpost on grid)   # N_beta = 1
def grid_2d(data, mode, b1_grid, b2_grid) -> (logpost on grid)              # N_beta = 2
def centred_gibbs_reference(data, field_spec, hyper, mode, rng, n_iter) -> samples of (beta, mu, log sigma, log rho)
    # beta | theta, d : random-walk MH (adaptive during a discarded burn-in only)
    # mu | sigma, rho, beta : exact Gaussian conditional
    # (log sigma, log rho) | beta, mu : exact draw from a fine 2-D grid of the conditional density
    #   (includes -K log sigma - 1/2 logdet K(rho) and the hyperpriors and Jacobians)
    # field_spec: its OWN construction of V and lam from the same recipe (graph Laplacian or
    #   harmonics), written independently
def linear_gaussian_closed_form(F, m, N_inv, field_spec, hyper, mu_grid, logsig_grid, u_grid) -> ...
    # d = m + F beta_region + n; returns beta | theta, d (mean, cov) and theta | d on the 3-D grid
```

## checks/ (scripts; each prints a PASS/FAIL line with the numbers and writes a JSON to checks/out/)
- `check_1a_gradient.py`: both posteriors, both families, three modes; tolerance |g_fd - g| <
  1e-6 |g| + 1e-4; also model.collapsed_loglike gradient against finite differences per index
  pixel, and against the C++ kernel `marginal_loglike_regions` for w0/w1 when importable.
- `check_1b_grid.py`: N_beta = 1 and 2 with real SED, grid vs NUTS (means within 4 SE by ESS, sd
  ratio 0.95 to 1.05), for gauss, w0, w1.
- `check_1c_hyper.py`: nside_beta 1, nside_d 16, NUTS vs centred_gibbs_reference, combined SE,
  ESS > 800 on log sigma required.
- `check_1d_linear.py`: closed form vs NUTS on the fixed-amplitude linear model.
- `check_1e_negative.py`: the three wrong targets of E012 v2 1e; must FAIL 1c or 1d visibly;
  report the margin.
- `sbc.py --config default --n 300`: writes ranks and the rank-histogram tests of E012 v2 L2.
- `coverage.py`: w0 vs w1 vs gauss on template amplitudes, bias and coverage versus S/N.

## Test quantities and reporting
Every check prints: the setting, ESS, the statistic, the threshold, PASS or FAIL. Nothing is
tuned to pass; a FAIL is logged in E012 as a result.
