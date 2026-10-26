"""E012 v2 L2: simulation-based calibration (SBC) of the index-field sampler (Talts et al. 2018).

For replicate i (seed s_i) the joint model is drawn exactly as the posterior assumes it,

    theta ~ hyperprior (simulate.draw_hyper),  z ~ N(0, I_K),  beta = mu 1 + sigma S(rho) z,
    a_p ~ N(m_p, T_p),  d_p = M_p(beta_r(p)) a_p + n_p,  n_p ~ N(0, N_p),

and NUTS (sampler.run_nuts) is run on posterior.NonCentredPosterior or CentredPosterior with the
collapsed likelihood of the same mode ("gauss" is the only mode for which generator and target
agree; w0 and w1 are accepted and flagged generator_matches_target = False). If the posterior is
sampled exactly, the rank of the true value of any scalar test quantity q among L independent
posterior draws is uniform on {0, ..., L} (L + 1 values). Test quantities (E012 v2 L2):

    mu, log sigma, u, four beta pixels (the two brightest and two faintest index pixels by the
    template synchrotron amplitude), the field mean of beta, z^T z,
    (beta - mu 1)^T K(rho)^-1 (beta - mu 1) / sigma^2  [identical to z^T z for beta in
    mu 1 + range(V); kept as a consistency check between unpack and field.Kinv_quadform],
    log p(beta | theta) = -K log sigma - 1/2 log det K(rho) - 1/2 z^T z  (constant dropped),
    the collapsed log likelihood log L(beta) (Modrak et al. 2022), and TARP (Lemos et al. 2023):
    for a reference point beta_ref the count of draws with |beta_i - beta_ref| < |beta* - beta_ref|
    among the L draws, which is again uniform on {0, ..., L}. Two reference schemes: "prior"
    (beta_ref an independent draw of the field prior) and "local" (beta_ref = posterior mean +
    5 posterior sd x Uniform(-1, 1) per pixel, with the mean and sd taken from a SEPARATE short
    pilot chain of pilot_samples draws, never from the chain whose draws are ranked, so the
    reference point is independent of the ranked draws given the data).

Thinning and the escalation loop: the chain of n_samples states is reduced to L = 99 evenly
spaced draws. An attempt is clean when the Geyer ESS of every test quantity exceeds ess_min
(400), the sampling phase had no zero-density evaluations and no divergences. When an attempt is
not clean, the replicate is rerun (up to max_tries attempts, default 5) with the sample budget
doubled and target_accept raised along 0.8 -> 0.9 -> 0.95 -> 0.99 (then held) at every step,
whatever the reason for the previous failure. The ranks of the LAST attempt are always recorded;
a replicate whose last attempt is still not clean gets status "unclean" (sbc_analyze.py includes
it in the primary analysis and reports it alongside). Replicates are never dropped.

Amplitude prior widths (SPEC "Decisions", simulate.amplitude_prior_sd): "default" T^(1/2) =
max(0.5 |m|, 20) uK_RJ for synchrotron and dust and 70 uK_RJ for the CMB; "wide" T^(1/2) =
10 max(|m|, 10) uK_RJ (the first-build setting).

Output: one JSON line per replicate appended to <out>/<config>.jsonl (truth, ranks, ESS, TARP
counts, divergences, n_loglike_errors_sampling, seconds, per attempt); --resume skips indices
already present (--redo-unclean redoes those whose latest record is not "ok", with a 2x budget).
Run from nextgen/prototypes:  python -m indexfield.sbc --config default --n 300 --workers 8
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
import traceback

# One BLAS thread per process: the batched 3x3 solves of the likelihood gain nothing from
# threading, and parallel replicates (--workers) must not oversubscribe the machine.
for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

import numpy as np                 # noqa: E402
import healpy as hp                # noqa: E402

if __package__ in (None, ""):                         # run as a script from indexfield/
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from indexfield import field, posterior, sampler, simulate            # noqa: E402
from indexfield.model import collapsed_loglike, hier_loglike, hier_prior   # noqa: E402

QUANTITIES = ["mu", "log_sigma", "u", "beta_bright0", "beta_bright1", "beta_faint0",
              "beta_faint1", "beta_mean", "ztz", "quad_Kinv", "logprior_beta", "loglike",
              "tarp_prior_dist0", "tarp_local_dist0"]
HIER_QUANTITIES = ["mu_sync", "mu_dust", "log_s_sync", "log_s_dust"]   # mode "hier" (E016)
TARP_SCHEMES = ("prior", "local")


def quantities_for(mode: str) -> list:
    return QUANTITIES + (HIER_QUANTITIES if mode == "hier" else [])
DEFAULTS = dict(config="default", param="noncentred", family="graph", mode="gauss",
                rho_prior="loguniform", nside_beta=4, nside_d=16, noise_scale=1.0,
                amp_prior="default", warmup=500, samples=2000, ess_min=400.0, L=99,
                n_tarp_ref=4, max_depth=10, target_accept=0.8, max_tries=5, pilot_samples=500)
ACCEPT_LADDER = (0.8, 0.9, 0.95, 0.99)     # target_accept of attempts 1, 2, 3, 4 and later
UNCLEAN_STATUSES = ("unclean", "rerun")    # "rerun" is the name used by records before this fix


# ------------------------------------------------------------------------------------ utilities
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


def config_name(cfg: dict) -> str:
    """<config> plus a suffix for every option that departs from the default configuration."""
    parts = [cfg["config"]]
    for key, tag in (("param", None), ("family", None), ("mode", None), ("rho_prior", None),
                     ("nside_beta", "nb"), ("nside_d", "nd"), ("noise_scale", "noise"),
                     ("amp_prior", None)):
        if cfg[key] != DEFAULTS[key]:
            parts.append(f"{tag}{cfg[key]:g}" if tag else str(cfg[key]))
    return "_".join(parts)


def make_field(cfg: dict):
    basis = (field.GraphLaplacianBasis(cfg["nside_beta"]) if cfg["family"] == "graph"
             else field.HarmonicBasis(cfg["nside_beta"]))
    return field.MaternField(basis)


def draw_hier_replicate(cfg: dict, fld, hyper, rng):
    """E016 generator: (theta, z, beta) as simulate.simulate, then (mu_c, log s_c) from the
    HierPosterior hyperpriors and a_p,c ~ N(mu_c g_p,c, (s_c g_p,c)^2) with g the templates,
    cmb ~ N(0, 70^2); so the "hier" target is the exact marginal of the generator."""
    N = hp.nside2npix(cfg["nside_d"])
    Nb = hp.nside2npix(cfg["nside_beta"])
    region = simulate.region_map(cfg["nside_d"], cfg["nside_beta"])
    theta = simulate.draw_hyper(hyper, rng)
    z = rng.standard_normal(fld.K)
    beta = theta["mu"] + theta["sigma"] * np.asarray(fld.S_times(theta["rho"], z))
    g = simulate.templates(cfg["nside_d"])
    HP = posterior.HierPosterior
    amp = np.array([rng.normal(HP.HIER_MU0, HP.HIER_MU_SD), rng.normal(HP.HIER_MU0, HP.HIER_MU_SD),
                    rng.normal(np.log(HP.HIER_S0), HP.HIER_LOGS_SD),
                    rng.normal(np.log(HP.HIER_S0), HP.HIER_LOGS_SD)])
    m, inv_T = hier_prior(g, amp[:2], amp[2:])
    a = m + rng.standard_normal(m.shape) / np.sqrt(inv_T)
    nu = simulate.BANDS_GHZ
    rms = simulate.noise_rms(cfg["nside_d"], nu, cfg["noise_scale"])
    M, _ = simulate.mixing(nu, beta[region], beta_d=1.6, T_d=20.0)
    d = np.einsum("nbc,nc->nb", M, a) + rms * rng.standard_normal(rms.shape)
    data = simulate.PixelData(nu=nu, d=d, inv_var=1.0 / rms ** 2, region=region, n_beta=Nb,
                              beta_d=1.6, T_d=20.0)
    return data, dict(theta=theta, z=z, beta=beta, a=a, amp=amp, g=g)


def draw_replicate(cfg: dict, fld, hyper, rng):
    """One joint draw (theta, z, beta, a, d) of the hierarchical model; returns (PixelData, truth).

    simulate.simulate with amp_mode "gauss" and prior_width = cfg["amp_prior"]: PixelData.m and
    .inv_T carry exactly the (m, T) the amplitudes were drawn with, so the "gauss" target
    matches the generator (the only configuration in which the SBC tests the sampler).
    """
    if cfg["mode"] == "hier":
        return draw_hier_replicate(cfg, fld, hyper, rng)
    return simulate.simulate(cfg["nside_d"], cfg["nside_beta"], fld, hyper, rng,
                             amp_mode="gauss", noise_scale=cfg["noise_scale"],
                             prior_width=cfg["amp_prior"])


def test_pixels(cfg: dict) -> dict:
    """Index pixels of the two brightest and two faintest template synchrotron amplitudes
    (template sync at 30 GHz averaged over the data pixels of each index pixel)."""
    region = simulate.region_map(cfg["nside_d"], cfg["nside_beta"])
    sync = simulate.templates(cfg["nside_d"])[:, 1]
    nb = hp.nside2npix(cfg["nside_beta"])
    mean_sync = np.bincount(region, weights=sync, minlength=nb) / np.bincount(region, minlength=nb)
    order = np.argsort(mean_sync)
    return dict(bright=[int(order[-1]), int(order[-2])], faint=[int(order[0]), int(order[1])],
                mean_sync_bright=[float(mean_sync[order[-1]]), float(mean_sync[order[-2]])],
                mean_sync_faint=[float(mean_sync[order[0]]), float(mean_sync[order[1]])])


def prior_field_draw(fld, hyper, rng) -> np.ndarray:
    th = simulate.draw_hyper(hyper, rng)
    return th["mu"] + th["sigma"] * fld.S_times(th["rho"], rng.standard_normal(fld.K))


# ------------------------------------------------------------------------- chain -> quantities
def chain_quantities(post, fld, samples, logp, pix, refs):
    """Test-quantity matrix (n, len(QUANTITIES)) of a chain, plus the beta matrix (n, Nb).

    refs: dict(prior=(n_ref, Nb), local=(n_ref, Nb)) reference points; only the first of each
    scheme enters QUANTITIES (as a distance) for the ESS precondition. The collapsed log
    likelihood is recovered from the cached log posterior by subtracting the prior and
    hyperprior terms (exact, see posterior.py); it is verified against direct evaluation on the
    thinned draws by the caller.
    """
    n = samples.shape[0]
    K = post.K
    beta = np.empty((n, post.n_beta))
    mu = samples[:, K]
    ls = samples[:, K + 1]
    u = samples[:, K + 2]
    ztz = np.empty(n)
    quad = np.empty(n)
    logdet = np.empty(n)
    lp_h = np.empty(n)
    for i in range(n):
        up = post.unpack(samples[i])
        rho = up["rho"]
        w, _ = fld.weights(rho)
        beta[i] = up["beta"]
        z = up["z"] if "z" in up else up["y"] / (up["sigma"] * w)
        ztz[i] = z @ z
        quad[i] = fld.Kinv_quadform(rho, beta[i] - up["mu"]) / up["sigma"] ** 2
        logdet[i] = 2.0 * np.sum(np.log(w))
        lp_h[i] = post._hyper_terms(up["mu"], up["log_sigma"], up["u"])[0]
        if "mu_amp" in up:
            lp_h[i] += post.amp_hyper_terms(samples[i, K + 3:K + 7])[0]
    logprior = -K * ls - 0.5 * logdet - 0.5 * ztz
    if post.centred:
        loglike = logp + 0.5 * ztz + K * ls + 0.5 * logdet - lp_h
    else:
        loglike = logp + 0.5 * ztz - lp_h
    d_prior = np.linalg.norm(beta - refs["prior"][0], axis=1)
    d_local = np.linalg.norm(beta - refs["local"][0], axis=1)
    Q = np.column_stack([mu, ls, u, beta[:, pix["bright"][0]], beta[:, pix["bright"][1]],
                         beta[:, pix["faint"][0]], beta[:, pix["faint"][1]], beta.mean(axis=1),
                         ztz, quad, logprior, loglike, d_prior, d_local])
    if post.dim == K + 7:                                   # mode "hier": the amplitude block
        Q = np.column_stack([Q, samples[:, K + 3:K + 7]])
    return Q, beta


def truth_quantities(fld, truth, data, mode, pix, refs):
    th, z, beta = truth["theta"], truth["z"], truth["beta"]
    K = fld.K
    w, _ = fld.weights(th["rho"])
    ztz = float(z @ z)
    quad = fld.Kinv_quadform(th["rho"], beta - th["mu"]) / th["sigma"] ** 2
    logprior = -K * th["log_sigma"] - np.sum(np.log(w)) - 0.5 * ztz
    if mode == "hier":
        loglike = hier_loglike(beta, data, truth["g"], truth["amp"][:2], truth["amp"][2:])[0]
    else:
        loglike = collapsed_loglike(beta, data, mode)[0]
    q = [th["mu"], th["log_sigma"], th["u"], beta[pix["bright"][0]],
         beta[pix["bright"][1]], beta[pix["faint"][0]], beta[pix["faint"][1]],
         beta.mean(), ztz, quad, logprior, loglike,
         np.linalg.norm(beta - refs["prior"][0]),
         np.linalg.norm(beta - refs["local"][0])]
    if mode == "hier":
        q += list(truth["amp"])
    return np.array(q)


def tarp_counts(beta_thin, beta_true, refs):
    """For each reference point, the number of thinned draws closer to it than the truth."""
    out = {}
    for scheme, R in refs.items():
        d_true = np.linalg.norm(beta_true[None, :] - R, axis=1)                 # (n_ref,)
        d = np.linalg.norm(beta_thin[:, None, :] - R[None, :, :], axis=2)       # (L, n_ref)
        out[scheme] = (d < d_true[None, :]).sum(axis=0).astype(int)
    return out


# ------------------------------------------------------------------------------- one replicate
def run_one_chain(post, rng, n_warmup, n_samples, target_accept, max_depth):
    """NUTS from a random initial state (redrawn if it has zero density)."""
    for _ in range(50):
        x0 = post.random_init(rng)
        if np.isfinite(post.logp_and_grad(x0)[0]):
            break
    else:
        raise RuntimeError("no finite random initial state in 50 tries")
    return sampler.run_nuts(post, rng, n_warmup, n_samples, x0=x0, target_accept=target_accept,
                            max_depth=max_depth)


def escalate(cfg: dict, attempt: int):
    """(n_samples, target_accept) of attempt number `attempt` (1-based): the budget doubles and
    target_accept climbs ACCEPT_LADDER at every retry, independent of why the previous attempt
    was not clean; the user's --target-accept is a floor."""
    n_samples = cfg["samples"] * 2 ** (attempt - 1)
    ladder = ACCEPT_LADDER[min(attempt, len(ACCEPT_LADDER)) - 1]
    return n_samples, max(cfg["target_accept"], ladder)


def run_replicate(cfg: dict, index: int) -> dict:
    """Draw replicate `index`, sample it (escalating retries, see the module docstring) and
    return its record. The ranks are those of the last attempt, clean or not."""
    t_start = time.time()
    rng = np.random.default_rng([cfg["seed"], index])
    fld = make_field(cfg)
    hyper = posterior.HyperPrior(rho_prior=cfg["rho_prior"]).resolved(fld.basis.h)
    data, truth = draw_replicate(cfg, fld, hyper, rng)
    pix = test_pixels(cfg)
    names = quantities_for(cfg["mode"])
    if cfg["mode"] == "hier":
        post = posterior.hier_posterior_class(cfg["param"])(fld, hyper, data, truth["g"])
    else:
        Post = posterior.CentredPosterior if cfg["param"] == "centred" else \
            posterior.NonCentredPosterior
        post = Post(fld, hyper, data=data, mode=cfg["mode"])
    rng_ref = np.random.default_rng([cfg["seed"], index, 7])
    refs = dict(prior=np.array([prior_field_draw(fld, hyper, rng_ref)
                                for _ in range(cfg["n_tarp_ref"])]))
    # Local TARP references from a separate pilot chain (posterior scale of beta only); the
    # pilot's draws are not ranked and the ranked chains below start from a fresh random state.
    t0 = time.time()
    pilot = run_one_chain(post, rng, cfg["warmup"], cfg["pilot_samples"], cfg["target_accept"],
                          cfg["max_depth"])
    beta_pilot = np.array([post.beta(x) for x in pilot.samples])
    refs["local"] = (beta_pilot.mean(axis=0)
                     + 5.0 * beta_pilot.std(axis=0) * rng_ref.uniform(-1, 1, (cfg["n_tarp_ref"],
                                                                           post.n_beta)))
    pilot_info = dict(n_warmup=cfg["warmup"], n_samples=cfg["pilot_samples"],
                      n_divergent=pilot.n_divergent,
                      n_loglike_errors_sampling=pilot.n_loglike_errors_sampling,
                      step_size=pilot.step_size, seconds=time.time() - t0)
    q_true = truth_quantities(fld, truth, data, cfg["mode"], pix, refs)
    L = cfg["L"]
    attempts = []
    n_warmup = cfg["warmup"]
    for attempt in range(1, cfg["max_tries"] + 1):
        n_samples, target_accept = escalate(cfg, attempt)
        t0 = time.time()
        res = run_one_chain(post, rng, n_warmup, n_samples, target_accept, cfg["max_depth"])
        Q, beta_all = chain_quantities(post, fld, res.samples, res.logp, pix, refs)
        ess = sampler.ess_geyer(Q)
        i_min = int(np.nanargmin(ess))
        reasons = []
        if not ess[i_min] > cfg["ess_min"]:
            reasons.append("ess")
        if res.n_divergent > 0:
            reasons.append("divergent")
        if res.n_loglike_errors_sampling > 0:
            reasons.append("loglike_errors")
        attempts.append(dict(attempt=attempt, n_warmup=n_warmup, n_samples=n_samples,
                             target_accept=target_accept, min_ess=float(ess[i_min]),
                             min_ess_quantity=names[i_min], n_divergent=res.n_divergent,
                             n_divergent_warmup=res.n_divergent_warmup,
                             n_loglike_errors_sampling=res.n_loglike_errors_sampling,
                             n_loglike_errors=res.n_loglike_errors, step_size=res.step_size,
                             ebfmi=res.ebfmi, accept_stat=res.accept_stat,
                             n_leapfrog=res.n_leapfrog, n_leapfrog_warmup=res.n_leapfrog_warmup,
                             tree_depth_mean=res.tree_depth_mean,
                             tree_depth_max=res.tree_depth_max, reasons=reasons,
                             seconds=time.time() - t0))
        if not reasons:
            break
    status = "ok" if not reasons else "unclean"
    # Thin the LAST attempt to L evenly spaced draws and rank the truth among them.
    idx = np.round(np.linspace(0, Q.shape[0] - 1, L)).astype(int)
    Q_thin, beta_thin = Q[idx], beta_all[idx]
    ranks = (Q_thin < q_true[None, :]).sum(axis=0).astype(int)
    tarp = tarp_counts(beta_thin, truth["beta"], refs)
    # Consistency of the recovered log likelihood with direct evaluation on the thinned draws.
    if cfg["mode"] == "hier":
        x_thin = res.samples[idx]
        direct = np.array([hier_loglike(b, data, truth["g"], x[post.K + 3:post.K + 5],
                                        x[post.K + 5:post.K + 7])[0]
                           for b, x in zip(beta_thin[::7], x_thin[::7])])
    else:
        direct = np.array([collapsed_loglike(b, data, cfg["mode"])[0] for b in beta_thin[::7]])
    ll_err = float(np.max(np.abs(direct - Q_thin[::7, names.index("loglike")])))
    ll_scale = float(np.max(np.abs(direct)))
    quad_vs_ztz = float(np.max(np.abs(Q[:, names.index("quad_Kinv")]
                                      - Q[:, names.index("ztz")])))
    th = truth["theta"]
    return dict(index=index, seed=cfg["seed"], config_name=config_name(cfg), config=cfg,
                status=status, reasons=reasons, attempts=attempts, n_attempts=len(attempts),
                pilot=pilot_info,
                generator_matches_target=cfg["mode"] in ("gauss", "hier"), K=fld.K, dim=post.dim,
                pixels=pix, quantities=names, L=L, n_samples_final=int(Q.shape[0]),
                truth=dict(theta=th, beta=truth["beta"], q=dict(zip(names, q_true)),
                           amp=truth.get("amp")),
                ranks=dict(zip(names, ranks)),
                tarp_counts={k: v for k, v in tarp.items()},
                ess=dict(zip(names, ess)), min_ess=float(ess[i_min]),
                post_mean=dict(zip(names, Q.mean(axis=0))),
                post_sd=dict(zip(names, Q.std(axis=0, ddof=1))),
                n_divergent=res.n_divergent, n_divergent_warmup=res.n_divergent_warmup,
                n_loglike_errors_sampling=res.n_loglike_errors_sampling,
                n_loglike_errors=res.n_loglike_errors, ebfmi=res.ebfmi, step_size=res.step_size,
                accept_stat=res.accept_stat, n_leapfrog=res.n_leapfrog,
                tree_depth_mean=res.tree_depth_mean,
                checks=dict(loglike_recovery_maxabs=ll_err, loglike_scale=ll_scale,
                            quad_minus_ztz_maxabs=quad_vs_ztz),
                seconds=time.time() - t_start)


def _worker(args):
    cfg, index = args
    try:
        return run_replicate(cfg, index)
    except Exception as err:                      # noqa: BLE001  (recorded, never dropped)
        return dict(index=index, seed=cfg["seed"], config_name=config_name(cfg), config=cfg,
                    status="error", error=f"{type(err).__name__}: {err}",
                    traceback=traceback.format_exc(), seconds=float("nan"))


# -------------------------------------------------------------------------------------- driver
def load_records(path: str) -> dict:
    """Latest record per replicate index from a JSON-lines file (empty dict if absent)."""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                rec = json.loads(line)
                out[rec["index"]] = rec
    return out


def report_line(rec: dict) -> str:
    if rec["status"] == "error":
        return f"  rep {rec['index']:4d}  ERROR {rec['error'][:80]}"
    a = rec["attempts"][-1]
    return (f"  rep {rec['index']:4d}  {rec['status']:7s} attempts {rec['n_attempts']} "
            f"n {a['n_samples']} acc {a['target_accept']:.2f} minESS {rec['min_ess']:6.0f} "
            f"({a['min_ess_quantity']}) "
            f"div {rec['n_divergent']} err {rec['n_loglike_errors_sampling']} "
            f"ebfmi {rec['ebfmi']:.2f} depth {rec['tree_depth_mean']:.1f} "
            f"leapfrog {rec['n_leapfrog']} step {rec['step_size']:.3g} "
            f"{rec['seconds']:6.0f} s")


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--config", default=DEFAULTS["config"], help="base configuration name")
    ap.add_argument("--n", type=int, default=300, help="number of replicates")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "sbc_out"))
    ap.add_argument("--param", choices=["centred", "noncentred"], default=DEFAULTS["param"])
    ap.add_argument("--family", choices=["graph", "harmonic"], default=DEFAULTS["family"])
    ap.add_argument("--mode", choices=["gauss", "w0", "w1", "hier"], default=DEFAULTS["mode"])
    ap.add_argument("--rho-prior", choices=["loguniform", "pc"], default=DEFAULTS["rho_prior"])
    ap.add_argument("--nside-beta", type=int, default=DEFAULTS["nside_beta"])
    ap.add_argument("--nside-d", type=int, default=DEFAULTS["nside_d"])
    ap.add_argument("--noise-scale", type=float, default=DEFAULTS["noise_scale"])
    ap.add_argument("--amp-prior", choices=["default", "wide"], default=DEFAULTS["amp_prior"])
    ap.add_argument("--warmup", type=int, default=DEFAULTS["warmup"])
    ap.add_argument("--samples", type=int, default=DEFAULTS["samples"])
    ap.add_argument("--ess-min", type=float, default=DEFAULTS["ess_min"])
    ap.add_argument("--L", type=int, default=DEFAULTS["L"], help="thinned draws per replicate")
    ap.add_argument("--n-tarp-ref", type=int, default=DEFAULTS["n_tarp_ref"])
    ap.add_argument("--max-depth", type=int, default=DEFAULTS["max_depth"])
    ap.add_argument("--target-accept", type=float, default=DEFAULTS["target_accept"])
    ap.add_argument("--max-tries", type=int, default=DEFAULTS["max_tries"],
                    help="attempts per replicate (budget x2 and target_accept up at each retry)")
    ap.add_argument("--pilot-samples", type=int, default=DEFAULTS["pilot_samples"],
                    help="draws of the separate pilot chain that sets the local TARP references")
    ap.add_argument("--resume", action="store_true", help="skip replicates already in the file")
    ap.add_argument("--redo-unclean", action="store_true",
                    help="with --resume: redo replicates whose latest status is not ok, 2x budget")
    ap.add_argument("--workers", type=int, default=1, help="parallel replicate processes")
    return ap.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    cfg = dict(config=a.config, param=a.param, family=a.family, mode=a.mode,
               rho_prior=a.rho_prior, nside_beta=a.nside_beta, nside_d=a.nside_d,
               noise_scale=a.noise_scale, amp_prior=a.amp_prior, warmup=a.warmup,
               samples=a.samples, ess_min=a.ess_min, L=a.L, n_tarp_ref=a.n_tarp_ref,
               max_depth=a.max_depth, target_accept=a.target_accept,
               max_tries=max(1, a.max_tries), pilot_samples=a.pilot_samples, seed=a.seed)
    name = config_name(cfg)
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, f"{name}.jsonl")
    done = load_records(path) if a.resume else {}
    jobs = []
    for i in range(a.n):
        rec = done.get(i)
        if rec is None:
            jobs.append((cfg, i))
        elif a.redo_unclean and rec["status"] != "ok":
            jobs.append((dict(cfg, samples=2 * cfg["samples"]), i))
    print(f"sbc: config {name}, {a.n} replicates ({len(jobs)} to run, {len(done)} in file), "
          f"param {a.param}, family {a.family}, mode {a.mode}, nside_beta {a.nside_beta}, "
          f"nside_d {a.nside_d}, warmup {a.warmup}, samples {a.samples}, ess_min {a.ess_min:g}, "
          f"L {a.L}, max_tries {cfg['max_tries']}, pilot {a.pilot_samples}, workers "
          f"{a.workers}, seed {a.seed} -> {path}")
    if cfg["mode"] not in ("gauss", "hier"):
        print("  WARNING: generator (gauss amplitudes) does not match the target mode; this run "
              "is not a sampler test (E012 v2)")
    t0 = time.time()
    n_ok = n_unclean = n_err = 0
    secs = []

    def record(rec):
        nonlocal n_ok, n_unclean, n_err
        with open(path, "a") as f:
            f.write(json.dumps(to_py(rec)) + "\n")
            f.flush()
        print(report_line(rec), flush=True)
        if rec["status"] == "ok":
            n_ok += 1
        elif rec["status"] in UNCLEAN_STATUSES:
            n_unclean += 1
        else:
            n_err += 1
        if np.isfinite(rec["seconds"]):
            secs.append(rec["seconds"])

    if a.workers <= 1:
        for job in jobs:
            record(_worker(job))
    else:
        ctx = mp.get_context("spawn")
        with ctx.Pool(a.workers) as pool:
            for rec in pool.imap_unordered(_worker, jobs):
                record(rec)
    wall = time.time() - t0
    if secs:
        print(f"sbc done: {len(jobs)} replicates in {wall:.0f} s wall ({a.workers} workers); per "
              f"replicate mean {np.mean(secs):.0f} s, median {np.median(secs):.0f} s, max "
              f"{np.max(secs):.0f} s; ok {n_ok}, unclean {n_unclean}, error {n_err}; "
              f"CPU-seconds total {np.sum(secs):.0f}")
    return 0 if n_err == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
