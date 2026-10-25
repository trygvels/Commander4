"""Analysis of the SBC ranks written by sbc.py (E012 v2 L2 pass criteria).

Input: <out>/<config>.jsonl, the latest record per replicate index; replicates with status
"ok" enter the tests, "rerun" and "error" records are counted and listed, never dropped silently
(--include-rerun adds the rerun ranks to the tests, flagged in the summary).

For each test quantity the n ranks r_i are in {0, ..., L} (L + 1 = 100 values for L = 99) and are
uniform under a correct sampler. Tests:

1. Rank histogram with 20 equal bins (bin = r // 5), chi-square statistic X^2 = sum (O - E)^2 / E
   with E = n / 20, p = P(chi^2_19 > X^2). Bonferroni over the q distinct quantities of the
   configuration: fail if min p < alpha / q (alpha = 0.01, the E012 v2 power statement).
   Quantities with identical rank vectors (z^T z and the K^-1 quadratic form, which coincide
   exactly) count once.
2. Uniformity of the set of per-quantity p-values: Kolmogorov-Smirnov against Uniform(0, 1)
   (low power with q ~ 13; reported, and it enters the pass only through its own alpha).
3. ECDF simultaneous confidence band (Sailynoja, Burkner and Vehtari 2022, simulation method):
   with C_k = #{r_i <= k}, k = 0..L-1, C_k ~ Binomial(n, p_k), p_k = (k + 1) / (L + 1) under
   uniformity. For M simulated uniform rank sets the smallest pointwise two-sided binomial tail
   probability over k is recorded; gamma is its alpha-quantile, and the band is
   [Binom.ppf(gamma / 2, n, p_k), Binom.ppf(1 - gamma / 2, n, p_k)] for every k, which the ECDF
   of a uniform sample stays inside with probability about 1 - alpha simultaneously over k (the
   realized simultaneous coverage is re-estimated on fresh simulations and printed). Plotted as
   the ECDF difference C_k / n - p_k.
4. TARP (Lemos et al. 2023): the counts c_i of draws closer to the reference than the truth are
   uniform on {0, ..., L}; the expected coverage probability ECP(alpha) = #{c_i / L <= alpha} / n
   is plotted against alpha with the same simultaneous band, for both reference schemes; the
   first reference of each scheme is also a test quantity in 1 and 3.

Per quantity the summary also prints the mean-rank z-score (bias: z = (mean r - L / 2) /
sqrt(((L + 1)^2 - 1) / (12 n)), positive when the posterior sits below the truth) and the rank
variance ratio var(r) / (((L + 1)^2 - 1) / 12) (below 1: posterior too wide; above 1: too
narrow or biased), which are the interpretable directions of a failure.

Outputs in <out>: <config>_summary.txt, <config>_analysis.json, <config>_rank_hist.png,
<config>_ecdf.png, <config>_tarp.png, <config>_pvalues.png.
Run from nextgen/prototypes:  python -m indexfield.sbc_analyze --config default
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy import stats

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from indexfield.sbc import QUANTITIES, load_records, to_py          # noqa: E402

N_BINS = 20


# --------------------------------------------------------------------------------------- tests
def rank_tests(ranks: np.ndarray, L: int) -> dict:
    """Chi-square on N_BINS bins, mean-rank bias z and variance ratio of one rank vector."""
    n = ranks.size
    bins = ranks // ((L + 1) // N_BINS)
    counts = np.bincount(bins, minlength=N_BINS)
    expected = n / N_BINS
    chi2 = float(np.sum((counts - expected) ** 2 / expected))
    p = float(stats.chi2.sf(chi2, N_BINS - 1))
    var_u = ((L + 1) ** 2 - 1) / 12.0
    z_mean = float((ranks.mean() - L / 2.0) / np.sqrt(var_u / n))
    var_ratio = float(ranks.var(ddof=1) / var_u) if n > 1 else float("nan")
    return dict(n=n, counts=counts, chi2=chi2, p=p, z_mean=z_mean, var_ratio=var_ratio)


def ecdf_band(n: int, L: int, alpha: float, n_sims: int, rng) -> dict:
    """Simultaneous 1 - alpha band for C_k = #{r_i <= k}, k = 0..L-1 (Sailynoja et al. 2022)."""
    k = np.arange(L)
    p_k = (k + 1) / (L + 1.0)

    def min_pointwise_p(counts):           # counts: (sims, L)
        lo = stats.binom.cdf(counts, n, p_k[None, :])
        hi = stats.binom.sf(counts - 1, n, p_k[None, :])
        return np.min(np.minimum(1.0, 2.0 * np.minimum(lo, hi)), axis=1)

    def simulate_counts(m):
        r = rng.integers(0, L + 1, size=(m, n))
        return np.stack([(r <= kk).sum(axis=1) for kk in k], axis=1)

    gamma = float(np.quantile(min_pointwise_p(simulate_counts(n_sims)), alpha))
    lower = stats.binom.ppf(gamma / 2.0, n, p_k)
    upper = stats.binom.ppf(1.0 - gamma / 2.0, n, p_k)
    fresh = simulate_counts(n_sims)
    inside = np.all((fresh >= lower[None, :]) & (fresh <= upper[None, :]), axis=1)
    return dict(p_k=p_k, lower=lower, upper=upper, gamma=gamma, n_sims=n_sims,
                realized_coverage=float(inside.mean()))


def ecdf_counts(ranks: np.ndarray, L: int) -> np.ndarray:
    return np.array([(ranks <= kk).sum() for kk in range(L)])


def ecdf_inside(ranks: np.ndarray, band: dict) -> tuple[bool, float]:
    """Whether the ECDF stays inside the band, and max |C_k / n - p_k|."""
    C = ecdf_counts(ranks, band["p_k"].size)
    ok = bool(np.all((C >= band["lower"]) & (C <= band["upper"])))
    return ok, float(np.max(np.abs(C / ranks.size - band["p_k"])))


# --------------------------------------------------------------------------------------- plots
def plot_all(prefix, names, rank_sets, results, band, tarp, L, n):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    q = len(names)
    ncol = 4
    nrow = int(np.ceil(q / ncol))
    # 1. Rank histograms with the 99 percent pointwise binomial band on the bin counts.
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.4 * nrow), squeeze=False)
    e = n / N_BINS
    half = stats.norm.ppf(0.995) * np.sqrt(n * (1 / N_BINS) * (1 - 1 / N_BINS))
    for ax, nm in zip(axes.ravel(), names):
        r = results[nm]
        ax.bar(np.arange(N_BINS), r["counts"], width=1.0, color="#4C72B0", edgecolor="white")
        ax.axhspan(e - half, e + half, color="grey", alpha=0.25, lw=0)
        ax.axhline(e, color="k", lw=0.8)
        ax.set_title(f"{nm}  p={r['p']:.3f}", fontsize=9)
        ax.set_xticks([])
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[q:]:
        ax.axis("off")
    fig.suptitle(f"SBC rank histograms, n = {n}, L = {L}, {N_BINS} bins (grey: 99% pointwise)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(prefix + "_rank_hist.png", dpi=130)
    plt.close(fig)
    # 2. ECDF differences with the simultaneous band.
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.4 * nrow), squeeze=False)
    p_k = band["p_k"]
    for ax, nm in zip(axes.ravel(), names):
        C = ecdf_counts(rank_sets[nm], L)
        ax.fill_between(p_k, band["lower"] / n - p_k, band["upper"] / n - p_k, color="grey",
                        alpha=0.3, lw=0)
        ax.plot(p_k, C / n - p_k, color="#C44E52", lw=1.2)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"{nm}  {'inside' if results[nm]['ecdf_inside'] else 'OUTSIDE'}",
                     fontsize=9)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[q:]:
        ax.axis("off")
    fig.suptitle(f"ECDF difference with the simultaneous {100 * (1 - band['alpha']):.0f}% band "
                 f"(gamma = {band['gamma']:.4f})", fontsize=10)
    fig.tight_layout()
    fig.savefig(prefix + "_ecdf.png", dpi=130)
    plt.close(fig)
    # 3. TARP expected coverage against credibility level.
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
    for ax, scheme in zip(axes, ("prior", "local")):
        t = tarp[scheme]
        ax.fill_between(p_k, band["lower"] / n, band["upper"] / n, color="grey", alpha=0.3, lw=0)
        ax.plot([0, 1], [0, 1], "k", lw=0.6)
        for j, ecp in enumerate(t["ecp"]):
            ax.plot(p_k, ecp, lw=1.0, alpha=0.9 if j == 0 else 0.5,
                    label=f"ref {j}" + (" (tested)" if j == 0 else ""))
        ax.set_title(f"TARP, {scheme} references: {t['n_inside']}/{t['n_ref']} inside band",
                     fontsize=9)
        ax.set_xlabel("credibility level alpha")
        ax.set_ylabel("expected coverage")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(prefix + "_tarp.png", dpi=130)
    plt.close(fig)
    # 4. p-value set.
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    pv = np.sort([results[nm]["p"] for nm in names])
    ax.step(np.concatenate([[0], pv]), np.linspace(0, 1, pv.size + 1), where="post")
    ax.plot([0, 1], [0, 1], "k", lw=0.6)
    ax.set_xlabel("chi-square p-value")
    ax.set_ylabel("ECDF over quantities")
    ax.set_title("p-value set (KS against uniform)", fontsize=9)
    fig.tight_layout()
    fig.savefig(prefix + "_pvalues.png", dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------------------------------- main
def analyze(path: str, alpha: float, n_sims: int, include_rerun: bool, seed: int,
            make_plots: bool = True) -> dict:
    recs = load_records(path)
    if not recs:
        raise SystemExit(f"no records in {path}")
    by_status = {}
    for r in recs.values():
        by_status.setdefault(r["status"], []).append(r["index"])
    used = [r for r in recs.values() if r["status"] == "ok"
            or (include_rerun and r["status"] == "rerun")]
    used.sort(key=lambda r: r["index"])
    if not used:
        raise SystemExit("no usable replicates")
    L = used[0]["L"]
    names = list(used[0]["quantities"])
    n = len(used)
    rank_sets = {nm: np.array([r["ranks"][nm] for r in used]) for nm in names}
    n_ref = len(used[0]["tarp_counts"]["prior"])
    tarp_counts = {s: np.array([r["tarp_counts"][s] for r in used]) for s in ("prior", "local")}
    for s in ("prior", "local"):
        rank_sets[f"tarp_{s}"] = tarp_counts[s][:, 0]
        names.append(f"tarp_{s}")
    # Distinct rank vectors (z^T z and the K^-1 form coincide exactly; count once).
    distinct, duplicates = [], {}
    for nm in names:
        twin = next((d for d in distinct if np.array_equal(rank_sets[d], rank_sets[nm])), None)
        if twin is None:
            distinct.append(nm)
        else:
            duplicates[nm] = twin
    rng = np.random.default_rng(seed)
    band = ecdf_band(n, L, alpha, n_sims, rng)
    band["alpha"] = alpha
    results = {}
    for nm in names:
        res = rank_tests(rank_sets[nm], L)
        res["ecdf_inside"], res["ecdf_max_abs_diff"] = ecdf_inside(rank_sets[nm], band)
        res["duplicate_of"] = duplicates.get(nm)
        results[nm] = res
    q = len(distinct)
    pvals = np.array([results[nm]["p"] for nm in distinct])
    bonferroni_level = alpha / q
    bonferroni_ok = bool(pvals.min() >= bonferroni_level)
    ks = stats.kstest(pvals, "uniform")
    ks_ok = bool(ks.pvalue >= alpha)
    ecdf_ok = all(results[nm]["ecdf_inside"] for nm in distinct)
    tarp = {}
    for s in ("prior", "local"):
        ecp = [ecdf_counts(tarp_counts[s][:, j], L) / n for j in range(n_ref)]
        inside = [ecdf_inside(tarp_counts[s][:, j], band)[0] for j in range(n_ref)]
        tarp[s] = dict(ecp=ecp, inside=inside, n_inside=int(sum(inside)), n_ref=n_ref)
    verdict = "PASS" if bonferroni_ok and ecdf_ok and ks_ok else "FAIL"
    # Sampler diagnostics over the used replicates.
    diag = dict(
        seconds_mean=float(np.mean([r["seconds"] for r in used])),
        seconds_median=float(np.median([r["seconds"] for r in used])),
        seconds_max=float(np.max([r["seconds"] for r in used])),
        min_ess_median=float(np.median([r["min_ess"] for r in used])),
        min_ess_min=float(np.min([r["min_ess"] for r in used])),
        ess_median_by_quantity={nm: float(np.median([r["ess"][nm] for r in used]))
                                for nm in used[0]["quantities"]},
        n_divergent_total=int(sum(r["n_divergent"] for r in used)),
        n_loglike_errors_sampling_total=int(sum(r["n_loglike_errors_sampling"] for r in used)),
        ebfmi_median=float(np.median([r["ebfmi"] for r in used])),
        ebfmi_min=float(np.min([r["ebfmi"] for r in used])),
        tree_depth_mean=float(np.mean([r["tree_depth_mean"] for r in used])),
        n_leapfrog_mean=float(np.mean([r["n_leapfrog"] for r in used])),
        n_attempts_mean=float(np.mean([r["n_attempts"] for r in used])),
        loglike_recovery_maxabs=float(np.max([r["checks"]["loglike_recovery_maxabs"]
                                              for r in used])),
        quad_minus_ztz_maxabs=float(np.max([r["checks"]["quad_minus_ztz_maxabs"] for r in used])),
        config=used[0]["config"], n_samples_final_median=float(np.median(
            [r["n_samples_final"] for r in used])))
    rerun_reasons = {}
    for r in recs.values():
        if r["status"] == "rerun":
            for reason in r["reasons"]:
                rerun_reasons[reason] = rerun_reasons.get(reason, 0) + 1
    out = dict(path=path, n_used=n, L=L, alpha=alpha, include_rerun=include_rerun,
               status_counts={k: len(v) for k, v in by_status.items()},
               rerun_indices=sorted(by_status.get("rerun", [])),
               error_indices=sorted(by_status.get("error", [])), rerun_reasons=rerun_reasons,
               quantities=names, distinct=distinct, duplicates=duplicates, results=results,
               bonferroni=dict(level=bonferroni_level, min_p=float(pvals.min()),
                               min_p_quantity=distinct[int(pvals.argmin())], ok=bonferroni_ok),
               ks_pvalues=dict(statistic=float(ks.statistic), p=float(ks.pvalue), ok=ks_ok),
               ecdf_band=dict(gamma=band["gamma"], realized_coverage=band["realized_coverage"],
                              n_sims=n_sims, all_inside=ecdf_ok),
               tarp={s: dict(inside=tarp[s]["inside"], n_inside=tarp[s]["n_inside"],
                             n_ref=n_ref) for s in tarp},
               diagnostics=diag, verdict=verdict)
    prefix = os.path.splitext(path)[0]
    if make_plots:
        try:
            plot_all(prefix, names, rank_sets, results, band, tarp, L, n)
            out["plots"] = [prefix + s for s in ("_rank_hist.png", "_ecdf.png", "_tarp.png",
                                                 "_pvalues.png")]
        except Exception as err:                 # noqa: BLE001
            out["plot_error"] = f"{type(err).__name__}: {err}"
    summary = write_summary(out, prefix + "_summary.txt")
    with open(prefix + "_analysis.json", "w") as f:
        json.dump(to_py(out), f, indent=1)
    print(summary)
    return out


def write_summary(out: dict, path: str) -> str:
    cfg = out["diagnostics"]["config"]
    d = out["diagnostics"]
    lines = [f"SBC summary: {os.path.basename(out['path'])}",
             f"  configuration: param {cfg['param']}, family {cfg['family']}, mode {cfg['mode']}, "
             f"rho_prior {cfg['rho_prior']}, nside_beta {cfg['nside_beta']}, nside_d "
             f"{cfg['nside_d']}, noise x{cfg['noise_scale']:g}, amp_prior {cfg['amp_prior']}, "
             f"warmup {cfg['warmup']}, samples {cfg['samples']} (x2 on retry), ess_min "
             f"{cfg['ess_min']:g}, L {out['L']}",
             f"  replicates: {out['status_counts']} ; used in tests: {out['n_used']}"
             f"{' (rerun included)' if out['include_rerun'] else ''}; rerun reasons "
             f"{out['rerun_reasons']}; rerun indices {out['rerun_indices']}; error indices "
             f"{out['error_indices']}",
             f"  sampler: seconds per replicate mean {d['seconds_mean']:.0f} / median "
             f"{d['seconds_median']:.0f} / max {d['seconds_max']:.0f}; attempts mean "
             f"{d['n_attempts_mean']:.2f}; final n median {d['n_samples_final_median']:.0f}; "
             f"min ESS median {d['min_ess_median']:.0f} (min {d['min_ess_min']:.0f}); "
             f"divergences {d['n_divergent_total']}; sampling loglike errors "
             f"{d['n_loglike_errors_sampling_total']}; E-BFMI median {d['ebfmi_median']:.2f} "
             f"(min {d['ebfmi_min']:.2f}); tree depth {d['tree_depth_mean']:.1f}; leapfrogs per "
             f"chain {d['n_leapfrog_mean']:.0f}",
             f"  checks: loglike recovery max |diff| {d['loglike_recovery_maxabs']:.2e}; "
             f"K^-1 form minus z^T z max |diff| {d['quad_minus_ztz_maxabs']:.2e}",
             f"  median ESS by quantity: " + ", ".join(
                 f"{k} {v:.0f}" for k, v in d["ess_median_by_quantity"].items()),
             f"  {'quantity':18s} {'chi2':>7s} {'p':>7s} {'z_mean':>7s} {'var_ratio':>9s} "
             f"{'ecdf_maxdiff':>12s}  ecdf"]
    for nm in out["quantities"]:
        r = out["results"][nm]
        dup = f"  (= {r['duplicate_of']}, not counted)" if r["duplicate_of"] else ""
        lines.append(f"  {nm:18s} {r['chi2']:7.2f} {r['p']:7.4f} {r['z_mean']:+7.2f} "
                     f"{r['var_ratio']:9.3f} {r['ecdf_max_abs_diff']:12.4f}  "
                     f"{'inside' if r['ecdf_inside'] else 'OUTSIDE'}{dup}")
    b, k, e = out["bonferroni"], out["ks_pvalues"], out["ecdf_band"]
    lines += [f"  Bonferroni over {len(out['distinct'])} distinct quantities at alpha "
              f"{out['alpha']:g}: level {b['level']:.5f}, min p {b['min_p']:.4f} "
              f"({b['min_p_quantity']})  {'PASS' if b['ok'] else 'FAIL'}",
              f"  KS of the p-value set against uniform: D {k['statistic']:.3f}, p {k['p']:.3f}  "
              f"{'PASS' if k['ok'] else 'FAIL'}",
              f"  ECDF simultaneous {100 * (1 - out['alpha']):.0f}% band (gamma {e['gamma']:.5f}, "
              f"realized coverage {e['realized_coverage']:.3f} on {e['n_sims']} sims): "
              f"{'all inside  PASS' if e['all_inside'] else 'some OUTSIDE  FAIL'}",
              f"  TARP: prior references {out['tarp']['prior']['n_inside']}/"
              f"{out['tarp']['prior']['n_ref']} ECP curves inside the band, local references "
              f"{out['tarp']['local']['n_inside']}/{out['tarp']['local']['n_ref']}",
              f"  VERDICT {out['verdict']} (n = {out['n_used']}; E012 v2 power at n = 300: sd "
              f"off by 0.8 detected with prob 0.84, bias of 0.5 sd with 1.0)"]
    if out.get("plot_error"):
        lines.append(f"  plots failed: {out['plot_error']}")
    text = "\n".join(lines)
    with open(path, "w") as f:
        f.write(text + "\n")
    return text


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("path", nargs="?", help="JSON-lines file (overrides --out/--config)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "sbc_out"))
    ap.add_argument("--config", default="default", help="config name, <out>/<config>.jsonl")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--band-sims", type=int, default=4000)
    ap.add_argument("--include-rerun", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)
    path = a.path or os.path.join(a.out, f"{a.config}.jsonl")
    out = analyze(path, a.alpha, a.band_sims, a.include_rerun, a.seed, not a.no_plots)
    return 0 if out["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
