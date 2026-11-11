"""Paper A1: tables and numbers of the real-data section from the E044 result files (named regions, mask).

Writes into nextgen/papers/paperA1/: gen_tab_same.tex, gen_tab_cmbprior.tex, gen_tab_moves.tex (complete tabular
environments) and gen_numbers.tex, which defines every number the prose of the section uses as \\N{key}. A key
whose result file does not exist yet is defined as a red "pending". No number of the section is typed by hand.

    python -I nextgen/scripts/paperA1_make_e044.py
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E44 = os.path.join(HERE, "..", "experiments", "E044")
OUT = os.path.join(HERE, "..", "papers", "paperA1")
NAMES = ["High Latitude", "Spur", "Intermediate", "Plane"]
SHORT = ["HL", "Spur", "Int", "Plane"]
HEAD = "\\begin{tabular}{lrrrr}\n\\toprule\n & High Latitude & Spur & Intermediate & Plane \\\\\n\\midrule\n"
FOOT = "\\bottomrule\n\\end{tabular}\n"
LAB = {"lcdm": r"$\Lambda$CDM covariance", "white": r"white, 0.2 $\mu$K", "flat": "flat"}
nums, missing = {}, []


def load(tag):
    p = os.path.join(E44, tag + ".npz")
    if not os.path.exists(p):
        missing.append(tag)
        return None
    return np.load(p, allow_pickle=True)


def fm(val, fmt):
    """Format a number; a value that rounds to zero is printed without a sign."""
    s = fmt.format(val)
    try:
        if float(s) == 0.0:
            s = s.replace("-", "").replace("+", "")
    except ValueError:
        pass
    return s.replace("-", "$-$")


def put(key, val, fmt="{:.3f}"):
    nums[key] = fm(val, fmt) if not isinstance(val, str) else val


def put4(key, vals, fmt="{:.3f}"):
    for s, v in zip(SHORT, vals):
        put(f"{key}.{s}", float(v), fmt)
    lo, hi = float(np.min(np.abs(vals))), float(np.max(np.abs(vals)))
    a, b = fmt.format(lo).replace("+", ""), fmt.format(hi).replace("+", "")
    nums[f"{key}.range"] = a if a == b else f"{a} to {b}"
    nums[f"{key}.list"] = ", ".join(fm(float(v), fmt) for v in vals[:-1]) + " and " + fm(float(vals[-1]), fmt)


def row(label, vals, f):
    return f"{label} & " + " & ".join(f(v) for v in vals) + " \\\\\n"


def val(x):
    return f"${x:.3f}$"


def sgn(x):
    return "$0.000$" if abs(x) < 0.0005 else f"${x:+.3f}$"


def tab_same(z):
    b, s = z["beta"], z["sd_sand"]
    t = HEAD + row("elements", z["nper"], lambda x: f"{int(x)}")
    t += row("profile likelihood", b, val)
    t += row(r"\quad sandwich error", s, lambda x: f"{x:.3f}")
    if "blk_CR3_own" in z.files:
        own, ne = z["blk_CR3_own"], z["blk_neff"]
        t += "\\quad block error, own scores & " + " & ".join((f"{x:.3f}" if n_ >= 4.0 and np.isfinite(x) else "") for x, n_ in zip(own, ne)) + " \\\\\n"
        t += row(r"\quad effective number of patches", ne, lambda x: f"{x:.1f}")
    if "flat23" in z.files:
        t += row("flat prior at 23 GHz", z["flat23"], val)
        t += row("flat prior at 30 GHz", z["flat30"], val)
    if "zm" in z.files:
        t += row("zero-mean prior, free scale", z["zm"], val)
    if "flat23" in z.files:
        t += "\\midrule\n" + row("flat 30 GHz minus flat 23 GHz", z["flat30"] - z["flat23"], sgn)
        t += row(r"\quad first order, curvature of the marginal", z["pred_marg"], sgn)
        t += row(r"\quad first order, $\J_{\rm c}$ (no flat-prior fit)", z["pred_jc"], sgn)
        t += row(r"\quad in sandwich errors", (z["flat30"] - z["flat23"]) / s, lambda x: f"{x:.1f}")
    t += "\\midrule\n" + row(r"$\chi^2$ per element and band", z["chi2_region"], lambda x: f"{x:.3f}")
    tj = np.array([z["T"][k, k] / (z["J"][k, k] - z["T"][k, k]) for k in range(4)])
    t += row(r"$\Tm/\J_{\rm c}$", tj, lambda x: f"{x:.3f}" if x >= 0.001 else f"{x:.4f}")
    return t + FOOT


UPPER = -1.6                                                 # upper bound of the index in lowres.fit


def at_bound(x):
    return abs(float(x) - UPPER) < 1e-6


def tab_cmbprior(R):
    """Absolute indices under three CMB priors; an optimum at the upper bound of the allowed range is marked."""
    mark = lambda x: f"${x:.3f}$" + (r"$^\dagger$" if at_bound(x) else "")
    t = HEAD + "\\multicolumn{5}{l}{\\emph{profile likelihood, with the sandwich error}} \\\\\n"
    for k, z in R.items():
        t += LAB[k] + " & " + " & ".join(f"${b:.3f} \\pm {e:.3f}$" for b, e in zip(z["beta"], z["sd_sand"])) + " \\\\\n"
    t += "\\multicolumn{5}{l}{\\emph{flat amplitude priors, synchrotron amplitude defined at 23 GHz}} \\\\\n"
    for k, z in R.items():
        t += row(LAB[k], z["flat23"], mark)
    t += "\\multicolumn{5}{l}{\\emph{flat amplitude priors, synchrotron amplitude defined at 30 GHz}} \\\\\n"
    for k, z in R.items():
        t += row(LAB[k], z["flat30"], mark)
    t += "\\multicolumn{5}{l}{\\emph{$\\chi^2$ per element and band (expectation 0.60, 0.60 and 0.40)}} \\\\\n"
    for k, z in R.items():
        t += row(LAB[k], z["chi2_region"], lambda x: f"{x:.3f}")
    return t + FOOT


def main():
    z = load("lcdm_all")
    if z is None:
        print("lcdm_all.npz missing: nothing written")
        return
    b0, s0 = z["beta"], z["sd_sand"]
    with open(os.path.join(OUT, "gen_tab_same.tex"), "w") as f:
        f.write("% generated by nextgen/scripts/paperA1_make_e044.py; do not edit\n" + tab_same(z))
    put4("prof", b0); put4("prof2", b0, "{:.2f}"); put4("sd", s0); put4("nper", z["nper"], "{:.0f}")
    put("nkept", int(z["keep"].sum()), "{:d}"); put("ntot", int(z["keep"].size), "{:d}")
    put4("chi", z["chi2_region"])
    # the total against its expectation: chi2_min (prior penalty included) has bands * n - 2 n - regions degrees of freedom
    chi2, dof, nb = float(z["chi2"]), int(z["dof"]), int(z["resid"].shape[0])
    nper = np.asarray(z["nper"], dtype=float)
    pen = chi2 - float(np.sum(z["chi2_region"] * nper * nb))
    put("chi.total", chi2, "{:.0f}"); put("chi.dof", dof, "{:d}"); put("chi.ratio", chi2 / dof); put("chi.nsig", abs(chi2 - dof) / np.sqrt(2.0 * dof), "{:.1f}")
    put("chi.noisepercent", 100 * (1 - np.sqrt(chi2 / dof)), "{:.1f}"); put("chi.penalty", pen, "{:.0f}")
    cexp = (dof - pen) / (nper.sum() * nb)
    cse = np.sqrt(2.0 * (nb - 2) * nper) / (nper * nb)
    put("chi.exp", cexp); put4("chi.se", cse); put4("chi.z", (z["chi2_region"] - cexp) / cse, "{:+.1f}")
    zf_, zk_ = load("flat"), load("KKa")
    if zf_ is not None:
        cl, cf = z["chi2_region"] * nper * nb, zf_["chi2_region"] * nper * nb          # residual chi-square by region, LCDM-prior and flat-prior fit
        dexp = nper - pen * nper / nper.sum()                                          # the prior penalty of the LCDM fit shared out by elements
        dr, dse = (cl - cf) / dexp, np.sqrt(2.0 / nper)
        fr, fse = cf / (2.0 * nper), np.sqrt(1.0 / nper)
        put4("chi.dir", dr); put4("chi.dir.se", dse); put4("chi.dir.z", (dr - 1.0) / dse, "{:+.1f}"); put4("chi.flat2n", fr); put4("chi.flat2n.se", fse)
        put4("chi.flat2n.z", (fr - 1.0) / fse, "{:+.1f}")
        put("chi.dir.Plane.percent", 100 * (float(dr[3]) - 1.0), "{:.0f}"); put("chi.dir.Plane.absz", abs(float((dr[3] - 1.0) / dse[3])), "{:.1f}")
        put("chi.dir.Plane.obs", float(cl[3] - cf[3]), "{:.0f}"); put("chi.dir.Plane.exp", float(dexp[3]), "{:.0f}")
        put("chi.dir.diffsig", abs((float(z["chi2"]) - float(zf_["chi2"])) / nper.sum() - float(zf_["chi2"]) / int(zf_["dof"])) / np.sqrt(2.0 / nper.sum() + 2.0 / int(zf_["dof"])), "{:.1f}")
        tot = lambda x, e: f"${x:.3f} \\pm {e:.3f}$"
        tc = "\\begin{tabular}{lrrrrr}\n\\toprule\n & High Latitude & Spur & Intermediate & Plane & all \\\\\n\\midrule\n"
        tc += "$\\Lambda$CDM prior: $\\chi^2$ per element and band & " + " & ".join(f"{x:.3f}" for x in z["chi2_region"]) + f" & {float(cl.sum()) / (nper.sum() * nb):.3f} \\\\\n"
        tc += "\\quad expectation & \\multicolumn{4}{c}{" + f"{cexp:.3f}, with sampling errors " + ", ".join(f"{x:.3f}" for x in cse) + "} & \\\\\n"
        tc += "flat CMB prior: $\\chi^2 / 2n$ & " + " & ".join(tot(x, e) for x, e in zip(fr, fse)) + " & " + tot(float(zf_["chi2"]) / int(zf_["dof"]), float(np.sqrt(2.0 / int(zf_["dof"])))) + " \\\\\n"
        tc += "along the CMB spectrum: difference of the two fits over its expectation & " + " & ".join(tot(x, e) for x, e in zip(dr, dse)) + " & " + tot((float(z["chi2"]) - float(zf_["chi2"])) / nper.sum(), float(np.sqrt(2.0 / nper.sum()))) + " \\\\\n"
        if zk_ is not None:
            ck = zk_["chi2_region"] * nper * int(zk_["resid"].shape[0])
            pk = float(zk_["chi2"]) - float(ck.sum())
            kr = ck / (nper - pk * nper / nper.sum())
            put4("chi.kka", kr); put("chi.kka.total", float(zk_["chi2"]) / int(zk_["dof"])); put("chi.kka.se", float(np.sqrt(2.0 / int(zk_["dof"])))); put("chi.kka.nsig", abs(float(zk_["chi2"]) - int(zk_["dof"])) / np.sqrt(2.0 * int(zk_["dof"])), "{:.1f}")
            tc += "K, Ka and 353 GHz only, $\\Lambda$CDM prior: $\\chi^2$ over its expectation & " + " & ".join(tot(x, e) for x, e in zip(kr, dse)) + " & " + tot(float(zk_["chi2"]) / int(zk_["dof"]), float(np.sqrt(2.0 / int(zk_["dof"])))) + " \\\\\n"
        with open(os.path.join(OUT, "gen_tab_chi.tex"), "w") as f:
            f.write("% generated by nextgen/scripts/paperA1_make_e044.py; do not edit\n" + tc + FOOT)
    zfl, zff = load("ffmask_lcdm"), load("ffmask_flat")
    if zfl is not None and zff is not None:                  # H16: the brightest thermal sources masked as well
        npf = np.asarray(zfl["nper"], dtype=float)
        clf, cff = zfl["chi2_region"] * npf * nb, zff["chi2_region"] * npf * nb
        penf = float(zfl["chi2"]) - float(clf.sum())
        exf = npf - penf * npf / npf.sum()
        drf = (clf - cff) / exf
        put4("ff.chi.dir", drf); put4("ff.chi.dir.se", np.sqrt(2.0 / npf)); put4("ff.chi.dir.z", (drf - 1.0) / np.sqrt(2.0 / npf), "{:+.1f}")
        put("ff.chi.dir.Plane.obs", float(clf[3] - cff[3]), "{:.0f}"); put("ff.chi.dir.Plane.exp", float(exf[3]), "{:.0f}")
        put("ff.npix.new", int(round((nper[3] - npf[3]) / 2)), "{:d}"); put("ff.npix.Plane", int(round(npf[3] / 2)), "{:d}")
        if zf_ is not None:
            gone_obs, gone_exp = float((cl - cf)[3] - (clf - cff)[3]), float(dexp[3] - exf[3])
            put("ff.gone.obs", gone_obs, "{:.0f}"); put("ff.gone.exp", gone_exp, "{:.0f}"); put("ff.gone.nsig", (gone_obs - gone_exp) / np.sqrt(2.0 * gone_exp), "{:.0f}")
        dff = zfl["beta"] - b0
        put4("ff.move", dff, "{:+.3f}"); put4("ff.prof", zfl["beta"]); put4("ff.sd", zfl["sd_sand"]); put("ff.move.Plane.abs", abs(float(dff[3])))
        put("ff.nest.Plane.abs", abs(float(dff[3])) / float(np.sqrt(zfl["sd_sand"][3] ** 2 - s0[3] ** 2)), "{:.1f}"); put("ff.movesd.Plane.abs", abs(float(dff[3])) / float(s0[3]), "{:.1f}")
        put("ff.chi.ratio", float(zfl["chi2"]) / int(zfl["dof"])); put("ff.flat.chi.ratio", float(zff["chi2"]) / int(zff["dof"]))
        put("ff.chi.dir.all", (float(zfl["chi2"]) - float(zff["chi2"])) / npf.sum()); put4("ff.cr3ownratio", zfl["blk_CR3_own"] / zfl["sd_sand"], "{:.2f}")
        if os.path.exists(os.path.join(OUT, "gen_tab_chi.tex")) and zf_ is not None:
            tcf = open(os.path.join(OUT, "gen_tab_chi.tex")).read()
            rowf = ("\\midrule\nalong the CMB spectrum, thermal sources masked as well & " + " & ".join(f"${x:.3f} \\pm {e:.3f}$" for x, e in zip(drf, np.sqrt(2.0 / npf)))
                    + f" & ${(float(zfl['chi2']) - float(zff['chi2'])) / npf.sum():.3f} \\pm {np.sqrt(2.0 / npf.sum()):.3f}$ \\\\\n")
            with open(os.path.join(OUT, "gen_tab_chi.tex"), "w") as f:
                f.write(tcf.replace(FOOT, rowf + FOOT))
    zsp = load("sim_n8_chi2split")
    if zsp is not None and zsp["cl"].shape[0] >= 10:
        npr = zsp["nper"]
        pn = (zsp["tl"] - zsp["cl"].sum(1))[:, None]
        Dn = (zsp["cl"] - zsp["cf"]) / (npr - pn * npr / npr.sum())
        Fn = zsp["cf"] / (2.0 * npr)
        put("chisplit.n", int(Dn.shape[0]), "{:d}"); put4("chisplit.dir.mean", Dn.mean(0)); put4("chisplit.dir.sd", Dn.std(0, ddof=1)); put4("chisplit.dir.sdover", Dn.std(0, ddof=1) / np.sqrt(2.0 / npr), "{:.2f}")
        put4("chisplit.flat.mean", Fn.mean(0)); put4("chisplit.flat.sdover", Fn.std(0, ddof=1) / np.sqrt(1.0 / npr), "{:.2f}")
        put("chisplit.dir.mean.se", float(np.sqrt(2.0 / npr[3]) / np.sqrt(Dn.shape[0])))
        if zf_ is not None:                                  # the Plane's excess with the sampling error measured under the model, scaled to the elements of the map fit
            put("chi.dir.Plane.absz.null", abs(float((dr[3] - 1.0) / (dse[3] * max(1.0, Dn.std(0, ddof=1)[3] / np.sqrt(2.0 / npr[3]))))), "{:.1f}")
    zc2 = load("sim_n8_chi2null")
    if zc2 is not None:
        put("chinull.n", int(zc2["ratio"].size), "{:d}"); put("chinull.mean", float(zc2["ratio"].mean())); put("chinull.sd", float(zc2["ratio"].std(ddof=1)))
        put("chinull.se", float(zc2["ratio"].std(ddof=1) / np.sqrt(zc2["ratio"].size))); put("chinull.theory", float(np.sqrt(2.0 / int(zc2["dof"]))))
        put4("chinull.creg", zc2["creg"].mean(0)); put4("chinull.cregsd", zc2["creg"].std(0, ddof=1))
    tj = np.array([z["T"][k, k] / (z["J"][k, k] - z["T"][k, k]) for k in range(4)])
    put("tj.max", float(tj.max()))
    put4("tj", tj)
    put4("snr2", 1.0 / tj, "{:.0f}")
    put4("nsig2", z["nper"] * s0 ** 2, "{:.1f}")
    if "flat23" in z.files:
        d23, d30, dd = z["flat23"] - b0, z["flat30"] - b0, z["flat30"] - z["flat23"]
        put4("f23", z["flat23"]); put4("f30", z["flat30"]); put4("f23m", d23, "{:+.3f}"); put4("f30m", d30, "{:+.3f}"); put4("d3023", dd, "{:+.3f}")
        put4("f23sd", d23 / s0, "{:+.1f}"); put4("f30sd", d30 / s0, "{:+.1f}"); put4("d3023sd", dd / s0, "{:.1f}")
        put4("predH", z["pred_marg"], "{:+.3f}"); put4("predJ", z["pred_jc"], "{:+.3f}")
        put4("ratioH", dd / z["pred_marg"], "{:.2f}"); put4("ratioJ", dd / z["pred_jc"], "{:.2f}")
        put("ratioH.percentHL", float(100 * (z["pred_marg"][0] / dd[0] - 1)), "{:.0f}")
        put4("f23two", z["flat23"], "{:.2f}"); put4("f30two", z["flat30"], "{:.2f}"); put4("d3023two", dd, "{:.2f}")
        put4("predH.over", 100 * (z["pred_marg"] / dd - 1), "{:.0f}"); put4("predJ.over", 100 * (z["pred_jc"] / dd - 1), "{:.0f}")
        hn = z["pred_marg"] / np.log(30 / 23)
        put4("pivot", 23.0 / np.exp(d23 / hn), "{:.1f}")
        put("flatconv", "yes" if bool(np.all(z["flat_converged"])) else "NO", "{}")
    if "zm" in z.files:
        put4("zm", z["zm"]); put4("zmsd", (z["zm"] - b0) / s0, "{:+.1f}"); put4("zmm", z["zm"] - b0, "{:+.3f}")
        put4("zmsd.abs", np.abs(z["zm"] - b0) / s0, "{:.1f}"); put4("zmm.abs", np.abs(z["zm"] - b0), "{:.3f}")
    if "blk_CR3" in z.files:
        put4("cr3", z["blk_CR3"]); put4("cr3own", z["blk_CR3_own"]); put4("neff", z["blk_neff"], "{:.1f}"); put4("levmax", z["blk_lev_max"], "{:.2f}")
        put4("cr3ratio", z["blk_CR3"] / s0, "{:.1f}"); put4("cr3ownratio", z["blk_CR3_own"] / s0, "{:.1f}")
        put("cross", float(z["blk_cross"]), "{:.2f}")
        with np.errstate(invalid="ignore", divide="ignore"):
            G, Dm = z["blk_G"], z["blk_D"]
            for j, sname in enumerate(SHORT):
                pos = Dm[:, j] > 0
                put(f"nblk.{sname}", int(pos.sum()), "{:d}")
    R = {"lcdm": z}
    for k in ("white", "flat"):
        zz = load(k)
        if zz is not None and "flat23" in zz.files:
            R[k] = zz
            put4(f"{k}.prof", zz["beta"]); put4(f"{k}.sd", zz["sd_sand"]); put4(f"{k}.move", zz["beta"] - b0, "{:+.3f}")
            put4(f"{k}.prof2", zz["beta"], "{:.2f}"); put4(f"{k}.sd2", zz["sd_sand"], "{:.2f}")
            put4(f"{k}.movesd", (zz["beta"] - b0) / s0, "{:+.1f}"); put4(f"{k}.sdratio", zz["sd_sand"] / s0, "{:.2f}")
            put4(f"{k}.d3023", zz["flat30"] - zz["flat23"], "{:+.3f}"); put4(f"{k}.f23m", zz["flat23"] - zz["beta"], "{:+.3f}")
            put4(f"{k}.d3023two", zz["flat30"] - zz["flat23"], "{:.2f}")
            put4(f"{k}.f30", zz["flat30"]); put4(f"{k}.f23", zz["flat23"]); put4(f"{k}.chi", zz["chi2_region"])
            nest = np.sqrt(np.abs(zz["sd_sand"] ** 2 - s0 ** 2))
            put4(f"{k}.movenest", (zz["beta"] - b0) / nest, "{:+.2f}"); put4(f"{k}.movenest.abs", np.abs((zz["beta"] - b0) / nest), "{:.2f}")
            inner = np.array([not (at_bound(p_) or at_bound(q_)) for p_, q_ in zip(zz["flat23"], zz["flat30"])])
            if inner.all():
                # first-order quantities (the Hessian is evaluated at the flat-prior optimum for 23 GHz): only where no optimum is on a bound
                put4(f"{k}.ratioH", (zz["flat30"] - zz["flat23"]) / zz["pred_marg"], "{:.2f}"); put4(f"{k}.ratioJ", (zz["flat30"] - zz["flat23"]) / zz["pred_jc"], "{:.2f}")
                hn = zz["pred_marg"] / np.log(30 / 23)
                put4(f"{k}.pivot", 23.0 / np.exp((zz["flat23"] - zz["beta"]) / hn), "{:.1f}")
            else:
                d23i, d30i = np.abs(zz["flat23"] - zz["beta"])[inner], np.abs(zz["flat30"] - zz["beta"])[inner]
                nums[f"{k}.f23m.inner.range"] = f"{d23i.min():.2f} to {d23i.max():.2f}"; nums[f"{k}.f30m.inner.range"] = f"{d30i.min():.2f} to {d30i.max():.2f}"
            put4(f"{k}.nsig2", zz["nper"] * zz["sd_sand"] ** 2, "{:.1f}")
            put(f"{k}.nbound23", int(sum(at_bound(x) for x in zz["flat23"])), "{:d}"); put(f"{k}.nbound30", int(sum(at_bound(x) for x in zz["flat30"])), "{:d}")
            put4(f"{k}.Jc", np.diag(zz["J"] - zz["T"]), "{:.0f}")
            c2, df = float(zz["chi2"]), int(zz["dof"])
            put(f"{k}.chi.total", c2, "{:.0f}"); put(f"{k}.chi.dof", df, "{:d}"); put(f"{k}.chi.ratio", c2 / df); put(f"{k}.chi.nsig", abs(c2 - df) / np.sqrt(2.0 * df), "{:.1f}")
            put(f"{k}.chi.noisepercent", 100 * (1 - np.sqrt(c2 / df)), "{:.1f}")
            if k == "flat":
                # the degrees of freedom that the LCDM-prior fit has beyond the flat-prior fit: one per element, chi-square_n under the model
                nel = int(z["keep"].sum())
                put("chi.cmbdir.ratio", (float(z["chi2"]) - c2) / nel); put("chi.cmbdir.se", float(np.sqrt(2.0 / nel))); put("flat.chi.se", float(np.sqrt(2.0 / df)))
        if zz is not None and "mm_bias" in zz.files:
            sdt = np.sqrt(np.diag(zz["mm_cov_true"]))
            put4("mm.bias", zz["mm_bias"], "{:+.4f}"); put4("mm.biassd", zz["mm_bias"] / zz["sd_sand"], "{:+.2f}"); put4("mm.sdtrue", sdt)
            put4("mm.ratio", sdt / zz["sd_sand"], "{:.2f}"); put4("mm.percent", 100 * (sdt / zz["sd_sand"] - 1), "{:.0f}")   # true error above the quoted one
            put4("mm.small.percent", 100 * (1 - zz["sd_sand"] / sdt), "{:.0f}")                                               # quoted error below the true one
            put4("mm.movepair.abs", np.abs(zz["beta"] - b0 - zz["mm_bias"]) / np.sqrt(np.clip(sdt ** 2 - s0 ** 2, 1e-30, None)), "{:.2f}")
            pair = np.sqrt(np.clip(sdt ** 2 - s0 ** 2, 1e-30, None))
            put4("mm.pairsd", pair); put4("mm.movepair", (zz["beta"] - b0 - zz["mm_bias"]) / pair, "{:+.2f}")
            put4("mm.overlcdm.percent", 100 * (sdt / s0 - 1), "{:.0f}"); put4("mm.quotedbelow.percent", 100 * (1 - zz["sd_sand"] / s0), "{:.0f}")
            put("mm.corr.HLInt", float(zz["mm_cov_true"][0, 2] / (sdt[0] * sdt[2])), "{:.1f}")
    if len(R) == 3:
        with open(os.path.join(OUT, "gen_tab_cmbprior.tex"), "w") as f:
            f.write("% generated by nextgen/scripts/paperA1_make_e044.py; do not edit\n" + tab_cmbprior(R))
    # arms
    t = HEAD + row(r"baseline: all scales, $\Lambda$CDM CMB prior", b0, val) + row("statistical error (sandwich)", s0, lambda x: f"{x:.3f}") + "\\midrule\n"
    if "flat30" in z.files:
        t += row("flat amplitude prior at 30 GHz", z["flat30"] - b0, sgn)
    arms = [("white", r"CMB prior: white, 0.2 $\mu$K"), ("flat", "CMB prior: flat"), ("lp30", r"index from $\ell \le 30$"), ("lp20", r"index from $\ell \le 20$"),
            ("lp12", r"index from $\ell \le 12$"), ("KKa", "K and Ka only"), ("bd148", "dust index 1.48 in place of 1.53"), ("bd158", "dust index 1.58"),
            ("nodrift", "no bandpass drift correction"), ("bp353full", "353 GHz: all bolometers"), ("delta", "delta bandpasses, nominal frequencies"),
            ("floor011", r"$\Lambda$CDM prior with a floor of 0.11 $\mu$K")]
    def with_units(label, d, u):
        """A row of moves, each followed in brackets by the move in units of its expected standard deviation."""
        z1 = lambda y: "0.0" if abs(y) < 0.05 else f"{y:+.1f}"
        return label + " & " + " & ".join(f"{sgn(x)} (${z1(y)}$)" for x, y in zip(d, u)) + " \\\\\n"

    def k_ka_lowpass():
        """K and Ka only against the four bands, both read from multipoles up to 20: one row and its numbers."""
        a, c = load("lp20_KKa"), load("lp20")
        if a is None or c is None:
            return ""
        d = a["beta"] - c["beta"]
        nest = np.sqrt(np.clip(a["sd_sand"] ** 2 - c["sd_sand"] ** 2, 1e-30, None))
        put4("lp20KKa.move", d, "{:+.3f}"); put4("lp20KKa.nest", d / nest, "{:+.1f}"); put("lp20KKa.nest.max", float(np.abs(d / nest).max()), "{:.1f}")
        put("lp20KKa.conv", "yes" if bool(a["converged"]) else "NO", "{}")
        worst.append(float(np.abs(d / nest).max()))
        return with_units(r"K and Ka only, both fits from $\ell \le 20$", d, d / nest)

    worst = []
    for tag, label in arms:
        zz = load(tag)
        if zz is None:
            continue
        d = zz["beta"] - b0
        put4(f"{tag}.move", d, "{:+.3f}"); put4(f"{tag}.prof", zz["beta"]); put4(f"{tag}.sd", zz["sd_sand"]); put4(f"{tag}.movesd", d / s0, "{:+.1f}")
        put4(f"{tag}.move.abs", np.abs(d), "{:.3f}")
        put(f"{tag}.move.abs.max", float(np.abs(d).max())); put(f"{tag}.movesd.abs.max", float(np.abs(d / s0).max()), "{:.1f}")
        put(f"{tag}.move.abs.max2", float(np.abs(d).max()), "{:.2f}")
        put4(f"{tag}.movesd.abs", np.abs(d / s0), "{:.1f}")
        if "chi2" in zz.files and "dof" in zz.files and not tag.startswith("lp"):
            put(f"{tag}.chi.ratio", float(zz["chi2"]) / int(zz["dof"])); put(f"{tag}.chi.nsig", abs(float(zz["chi2"]) - int(zz["dof"])) / np.sqrt(2.0 * int(zz["dof"])), "{:.1f}")
        put(f"{tag}.conv", "yes" if bool(zz["converged"]) else "NO", "{}")
        if tag.startswith("lp") or tag in ("KKa", "flat"):
            # the arm uses less information on the same data: the difference has variance sd_arm^2 - sd_base^2 under the model
            nest = np.sqrt(np.clip(zz["sd_sand"] ** 2 - s0 ** 2, 1e-30, None))
            put4(f"{tag}.nest", d / nest, "{:+.1f}"); put4(f"{tag}.nest2", d / nest, "{:+.2f}"); put4(f"{tag}.nestsd", nest)
            put(f"{tag}.nest.max", float(np.abs(d / nest).max()), "{:.1f}")
            put4(f"{tag}.nest.abs", np.abs(d / nest), "{:.1f}"); put4(f"{tag}.twonest", 2.0 * nest)
            worst.append(float(np.abs(d / nest).max()))
            t += with_units(label, d, d / nest)
            if tag == "KKa":
                t += k_ka_lowpass()
        elif tag == "white" and "mm_bias" in zz.files:
            sdt_ = np.sqrt(np.diag(zz["mm_cov_true"]))
            pair_ = np.sqrt(np.clip(sdt_ ** 2 - s0 ** 2, 1e-30, None))
            t += with_units(label, d, (d - zz["mm_bias"]) / pair_)
        else:
            t += row(label, d, sgn)
    if worst:
        put("arms.nest.max", max(worst), "{:.1f}")
    bd = [load("bd148"), load("bd158")]
    if all(x is not None for x in bd):
        dmove = 0.5 * (np.abs(bd[0]["beta"] - b0) + np.abs(bd[1]["beta"] - b0))
        slope = (bd[1]["beta"] - bd[0]["beta"]) / 0.10
        put4("dust.move", dmove); put4("dust.slope", slope, "{:.2f}"); put4("dust.movesd", dmove / s0, "{:.1f}")
        # the dust index that the fit itself prefers: a parabola through chi2_min at 1.48, 1.53 and 1.58
        c = np.array([float(bd[0]["chi2"]), float(z["chi2"]), float(bd[1]["chi2"])])
        a2 = (c[0] + c[2] - 2.0 * c[1]) / (2.0 * 0.05 ** 2)
        b1 = (c[2] - c[0]) / (2.0 * 0.05)
        best, bsd = 1.53 - b1 / (2.0 * a2), 1.0 / np.sqrt(a2)
        put("dust.best", float(best), "{:.2f}"); put("dust.bestsd", float(bsd), "{:.2f}"); put("dust.best3", float(best)); put("dust.bestsd3", float(bsd))
        put("dust.dchi148", float(c[0] - c[1]), "{:.1f}"); put("dust.dchi158", float(c[2] - c[1]), "{:+.1f}"); put("dust.dchi158.abs", abs(float(c[2] - c[1])), "{:.1f}")
        put4("dust.at002", np.abs(slope) * 0.02); put4("dust.at002sd", np.abs(slope) * 0.02 / s0, "{:.1f}")
    # the same low-pass arm on the earlier regions without a mask (E042 and E040), for the before and after
    po, pl = os.path.join(HERE, "..", "experiments", "E040", "s4_bpxv_n16_wm_lcdm.npz"), os.path.join(HERE, "..", "experiments", "E042", "lowpass20_lcdm.npz")
    if os.path.exists(po) and os.path.exists(pl):
        zo, zl = np.load(po), np.load(pl)
        bl, sl_ = zl["K_Ka_Q_V_353"][:4], zl["K_Ka_Q_V_353"][4:]
        dl = bl - zo["beta"]
        nl = np.sqrt(np.clip(sl_ ** 2 - zo["sd_sand"] ** 2, 1e-30, None))
        put4("old.lp20.move", dl, "{:+.3f}"); put4("old.lp20.move.abs", np.abs(dl)); put4("old.lp20.nest.abs", np.abs(dl / nl), "{:.1f}")
        put("old.lp20.move.absPlane", abs(float(dl[3]))); put4("old.lp20.nestsd", nl)
    with open(os.path.join(OUT, "gen_tab_moves.tex"), "w") as f:
        f.write("% generated by nextgen/scripts/paperA1_make_e044.py; do not edit\n" + t + FOOT)
    # the old regions with only Tau A masked (H6), against the unmasked fit of E040
    zt = load("bpxv_tauA")
    old = os.path.join(HERE, "..", "experiments", "E040", "s4_bpxv_n16_wm_lcdm.npz")
    if zt is not None and os.path.exists(old):
        zo = np.load(old)
        for j, sname in enumerate(("HL", "Spur", "GC", "Plane")):
            put(f"tauA.move.{sname}", float(zt["beta"][j] - zo["beta"][j]), "{:+.3f}"); put(f"tauA.prof.{sname}", float(zt["beta"][j])); put(f"old.prof.{sname}", float(zo["beta"][j]))
            put(f"tauA.absmove.{sname}", abs(float(zt["beta"][j] - zo["beta"][j])), "{:.2f}")
            put(f"tauA.chi.{sname}", float(zt["chi2_region"][j])); put(f"tauA.cr3own.{sname}", float(zt["blk_CR3_own"][j])); put(f"tauA.sd.{sname}", float(zt["sd_sand"][j]))
    # null simulations of the block-error ratio at the configuration of the fit
    s16 = [q for q in (load("sim_n16_blocknull"), load("sim_n16_blocknull_b")) if q is not None]
    if s16 and "blk_CR3_own" in z.files:
        q = np.concatenate([x["CR3_own"] / x["s"] for x in s16])
        zz_ = np.concatenate([(x["b"] - x["truth"]) / x["s"] for x in s16])
        m16 = q.shape[0]
        put("null16.n", int(m16), "{:d}")
        put4("null16.median", np.nanmedian(q, axis=0), "{:.2f}"); put4("null16.max", np.nanmax(q, axis=0), "{:.2f}"); put4("null16.p975", np.nanpercentile(q, 97.5, axis=0), "{:.2f}")
        mapr = z["blk_CR3_own"] / s0
        put4("null16.above", [(q[:, j] > mapr[j]).sum() for j in range(4)], "{:.0f}")
        put4("null16.sdz", zz_.std(0, ddof=1), "{:.2f}"); put4("null16.meanz", zz_.mean(0), "{:+.2f}")
        put("null16.se", 1.0 / np.sqrt(m16), "{:.2f}"); put("null16.sdz.se", 1.0 / np.sqrt(2.0 * (m16 - 1)), "{:.2f}")
        put4("null16.cover", np.mean(np.abs(zz_) < 1, axis=0), "{:.2f}")
        put("null16.pooled.meanz", float(zz_.mean()), "{:+.2f}"); put("null16.pooled.sdz", float(zz_.std(ddof=1)), "{:.2f}"); put("null16.pooled.se", 1.0 / np.sqrt(zz_.size), "{:.2f}")
    # the low-pass estimator and the zero-mean prior at nside 16 (few realizations), against the all-scale profile fits of the same realizations
    zl16, zb16 = load("sim_n16_lowpass"), load("sim_n16_blocknull")
    if zl16 is not None and zb16 is not None:
        m = min(zl16["b"].shape[0], zb16["b"].shape[0])
        nest = np.sqrt(np.clip(zl16["s"][:m] ** 2 - zb16["s"][:m] ** 2, 1e-30, None))
        pz = (zl16["b"][:m] - zb16["b"][:m]) / nest
        zt = (zl16["b"][:m] - zl16["truth"]) / zl16["s"][:m]
        put("lp16.n", int(m), "{:d}"); put("lp16.nvals", int(pz.size), "{:d}")
        put("lp16.paired.mean", float(pz.mean()), "{:+.2f}"); put("lp16.paired.rms", float(np.sqrt((pz ** 2).mean())), "{:.2f}"); put("lp16.paired.max", float(np.abs(pz).max()), "{:.1f}")
        put("lp16.z.mean", float(zt.mean()), "{:+.2f}"); put("lp16.z.rms", float(np.sqrt((zt ** 2).mean())), "{:.2f}")
    zz16 = load("sim_n16_zm")
    if zz16 is not None:
        pr = (zz16["b"] - zz16["prof"]) / zz16["s"]
        put("zm16.n", int(pr.shape[0]), "{:d}")
        put4("zm16.paired", pr.mean(0), "{:+.2f}"); put4("zm16.paired.abs", np.abs(pr.mean(0)), "{:.2f}"); put4("zm16.pairedindex.abs", np.abs((zz16["b"] - zz16["prof"]).mean(0)), "{:.3f}")
        put4("zm16.spread", pr.max(0) - pr.min(0), "{:.2f}")
        if "zm" in z.files:
            put4("zm16.overmap", pr.mean(0) / ((z["zm"] - b0) / s0), "{:.2f}")
    zc16 = load("sim_n16_chi2null")
    if zc16 is not None:
        put("chinull16.n", int(zc16["ratio"].size), "{:d}"); put("chinull16.mean", float(zc16["ratio"].mean())); put("chinull16.sd", float(zc16["ratio"].std(ddof=1)))
        put("chinull16.se", float(zc16["ratio"].std(ddof=1) / np.sqrt(zc16["ratio"].size))); put4("chinull16.creg", zc16["creg"].mean(0))
        if "cband" in zc16.files:
            nums["chinull16.cband.list"] = ", ".join(f"{x:.3f}" for x in zc16["cband"].mean(0)[:4])
    for tag, key in (("sim_n8_zm", "zmsim"),):
        zs = load(tag)
        if zs is None:
            continue
        dz = (zs["b"] - zs["truth"]) / zs["s"]
        dp = (zs["prof"] - zs["truth"]) / zs["s"]
        m = dz.shape[0]
        put(f"{key}.n", int(m), "{:d}")
        put4(f"{key}.meanz", dz.mean(0), "{:+.2f}"); put4(f"{key}.profmeanz", dp.mean(0), "{:+.2f}")
        put4(f"{key}.paired", ((zs["b"] - zs["prof"]) / zs["s"]).mean(0), "{:+.2f}")
        put4(f"{key}.paired.abs", np.abs(((zs["b"] - zs["prof"]) / zs["s"]).mean(0)), "{:.2f}")
        put4(f"{key}.pairedse", ((zs["b"] - zs["prof"]) / zs["s"]).std(0, ddof=1) / np.sqrt(m), "{:.2f}")
        put4(f"{key}.pairedscatter", ((zs["b"] - zs["prof"]) / zs["s"]).std(0, ddof=1), "{:.2f}")
        put4(f"{key}.pairedindex", (zs["b"] - zs["prof"]).mean(0), "{:+.3f}")
        put(f"{key}.se", 1.0 / np.sqrt(m), "{:.2f}")
    zb = load("sim_n8_blocknull_b")
    za = load("sim_n8_blocknull")
    zc = load("sim_n8_blocknull_c")
    if za is not None:
        both = [za] + [q for q in (zb, zc) if q is not None]
        zall = np.concatenate([(q["b"] - q["truth"]) / q["s"] for q in both])
        put("nullall.n", int(zall.shape[0]), "{:d}")
        put4("nullall.meanz", zall.mean(0), "{:+.2f}"); put4("nullall.sdz", zall.std(0, ddof=1), "{:.2f}")
        put("nullall.se", 1.0 / np.sqrt(zall.shape[0]), "{:.2f}")
        put4("nullall.cover2", np.mean(np.abs(zall) < 2, axis=0), "{:.2f}")
        put("nullall.sdz.se", 1.0 / np.sqrt(2.0 * (zall.shape[0] - 1)), "{:.2f}")
        put4("nullall.cover", np.mean(np.abs(zall) < 1, axis=0), "{:.2f}")
        put("nullall.cover.se", float(np.sqrt(0.683 * 0.317 / zall.shape[0])), "{:.2f}")
        qall = np.concatenate([q["CR3_own"] / q["s"] for q in both])
        put4("nullall.median", np.nanmedian(qall, axis=0), "{:.2f}"); put4("nullall.p975", np.nanpercentile(qall, 97.5, axis=0), "{:.2f}")
        put4("nullall.p84", np.nanpercentile(qall, 84, axis=0), "{:.2f}"); put4("nullall.max", np.nanmax(qall, axis=0), "{:.2f}")
        if "blk_CR3_own" in z.files:
            mapr = z["blk_CR3_own"] / s0
            put4("nullall.above", [(qall[:, j] > mapr[j]).sum() for j in range(4)], "{:.0f}")
    # scans of the flat-prior marginal along the High Latitude index (H11); slice 0 has the others at the profile optimum
    for k in ("flat", "lcdm"):
        zs = load("scan_" + k)
        if zs is None:
            continue
        for sl, sname in ((0, ""), (1, ".b")):
            rw = zs["rows"][zs["rows"][:, 0] == sl]
            if rw.shape[0] < 3:
                continue
            x, lp_, l23, l30 = rw[:, 1], rw[:, 2], rw[:, 3], rw[:, 4]
            key = f"scan.{k}{sname}"
            put(key + ".npts", int(x.size), "{:d}"); put(key + ".lo", float(x.min()), "{:.1f}"); put(key + ".hi", float(x.max()), "{:.1f}")
            put(key + ".mono23", "yes" if bool(np.all(np.diff(l23) > 0)) else "no", "{}"); put(key + ".mono30", "yes" if bool(np.all(np.diff(l30) > 0)) else "no", "{}")
            put(key + ".rise23", float(l23[-1]), "{:.0f}"); put(key + ".rise30", float(l30[-1]), "{:.0f}"); put(key + ".proffall", float(-lp_[-1]), "{:.0f}")
            put(key + ".low23", float(-l23[0]), "{:.0f}"); put(key + ".low30", float(-l30[0]), "{:.0f}")
            b_opt, n_r = float(zs["profile_opt"][int(zs["region"])]), float(zs["nper"][int(zs["region"])])
            j = int(np.searchsorted(x, b_opt))               # grid points j - 1 and j bracket the profile optimum
            if 2 <= j < x.size - 1:
                h = x[j] - x[j - 1]
                vol = l23 - lp_                              # the volume term, -(1/2) ln det A, for a reference of 23 GHz
                slope = float((vol[j] - vol[j - 1]) / h)
                put(key + ".volslope", slope, "{:.0f}"); put(key + ".volslope.per", slope / n_r, "{:.3f}"); put(key + ".pivot", 23.0 * np.exp(-slope / n_r), "{:.1f}")
                c_lo, c_hi = -(l23[j - 2] - 2 * l23[j - 1] + l23[j]) / h ** 2, -(l23[j - 1] - 2 * l23[j] + l23[j + 1]) / h ** 2
                put(key + ".curv23", 0.5 * float(c_lo + c_hi), "{:.0f}"); put(key + ".curvprof", -0.5 * float((lp_[j - 2] - 2 * lp_[j - 1] + lp_[j]) + (lp_[j - 1] - 2 * lp_[j] + lp_[j + 1])) / h ** 2, "{:.0f}")
                put(key + ".volgain", float(vol[-1] - np.interp(b_opt, x, vol)), "{:.0f}")
            for nm, y in (("23", l23), ("30", l30)):
                j = int(np.argmax(y))
                put(key + f".argmax{nm}", float(x[j]), "{:.2f}"); put(key + f".max{nm}", float(y[j]), "{:.1f}")
                put(key + f".interior{nm}", "yes" if 0 < j < x.size - 1 else "no", "{}")
                if 0 < j < x.size - 1:                       # vertex of the parabola through the three grid points around the maximum
                    y0, y1, y2 = y[j - 1], y[j], y[j + 1]
                    h = x[j + 1] - x[j]
                    put(key + f".vertex{nm}", float(x[j] + 0.5 * h * (y0 - y2) / (y0 - 2 * y1 + y2)), "{:.2f}")
    with open(os.path.join(OUT, "gen_numbers.tex"), "w") as f:
        f.write("% generated by nextgen/scripts/paperA1_make_e044.py; do not edit\n\\makeatletter\n")
        f.write("\\renewcommand{\\N}[1]{\\@ifundefined{num@#1}{\\textcolor{red}{[#1 pending]}}{\\csname num@#1\\endcsname}}\n")
        for k in sorted(nums):
            f.write(f"\\expandafter\\def\\csname num@{k}\\endcsname{{{nums[k]}}}\n")
        f.write("\\makeatother\n")
    print(f"wrote gen_tab_same.tex, gen_tab_moves.tex, gen_numbers.tex ({len(nums)} numbers)" + ("; gen_tab_cmbprior.tex" if len(R) == 3 else ""))
    if missing:
        print("missing result files:", sorted(set(missing)))


if __name__ == "__main__":
    main()
