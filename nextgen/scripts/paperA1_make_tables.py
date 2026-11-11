"""Paper A1: write the tabular bodies of the data-section tables from the result files of E040 and E042, so
that no number in those tables is typed by hand. Files written: nextgen/papers/paperA1/gen_tab_e040.tex,
gen_tab_cmbprior.tex, gen_tab_moves.tex. Rows whose result file does not exist yet are left out and listed.

    python -I nextgen/scripts/paperA1_make_tables.py
"""
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E40 = os.path.join(HERE, "..", "experiments", "E040")
E42 = os.path.join(HERE, "..", "experiments", "E042")
OUT = os.path.join(HERE, "..", "papers", "paperA1")
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
HEAD = "\\begin{tabular}{lrrrr}\n\\toprule\n & High Latitude & Spur & Galactic Center & Plane \\\\\n\\midrule\n"
FOOT = "\\bottomrule\n\\end{tabular}\n"
FILES = {"lcdm": "s4_bpxv_n16_wm_lcdm", "white": "bpxv_n16_wm", "flat": "s4_bpxv_n16_wm_cmbflat"}
LAB = {"lcdm": r"$\Lambda$CDM covariance", "white": r"white, 0.2 $\mu$K", "flat": "flat"}
missing = []


def full(tag):
    """All estimators of one run: dict of arrays (4,), plus chi2 per region from the text file."""
    npz, txt = os.path.join(E40, tag + ".npz"), os.path.join(E40, tag + ".txt")
    if not os.path.exists(npz):
        missing.append(tag)
        return None
    z = np.load(npz, allow_pickle=True)
    d = {k: np.asarray(z[k], dtype=float) for k in ("beta", "sd_sand", "flat23", "flat30", "zm", "pred_marg") if k in z.files}
    J, T = z["J"], z["T"]
    d["tj"] = np.array([T[k, k] / (J[k, k] - T[k, k]) for k in range(4)])
    for ln in open(txt):
        if ln.startswith("# chi2 by region"):
            d["chi"] = np.array([float(x) for x in re.findall(r"(?:High Latitude|Spur|Galactic Center|Plane) ([0-9.]+)", ln)])
    return d


def prof(path, key=None):
    if not os.path.exists(path):
        missing.append(os.path.basename(path))
        return None
    z = np.load(path, allow_pickle=True)
    if key is not None:
        v = np.asarray(z[key], dtype=float)
        return v[:4], v[4:]
    return np.asarray(z["beta"], dtype=float), np.asarray(z["sd_sand"], dtype=float)


def row(label, vals, f):
    return f"{label} & " + " & ".join(f(v) for v in vals) + " \\\\\n"


def val(x):
    return f"${x:.3f}$"


def sgn(x):
    return f"${x:+.3f}$"


def tab_e040():
    d = full(FILES["lcdm"])
    s = HEAD + "elements & 4388 & 928 & 182 & 646 \\\\\n"
    s += row("profile likelihood", d["beta"], val)
    s += "\\quad sandwich error & " + " & ".join((f"{x:.3f}" if k != 2 else f"({x:.3f})") for k, x in enumerate(d["sd_sand"])) + " \\\\\n"
    blk = os.path.join(E40, "block_v2_maps_lcdm.npz")
    if os.path.exists(blk):
        z = np.load(blk)
        cr3, neff = z["CR3"], z["neff"]
        s += "\\quad block error & " + " & ".join((f"{x:.3f}" if k != 2 else "") for k, x in enumerate(cr3)) + " \\\\\n"
        s += "\\quad effective number of patches & " + " & ".join(f"{x:.1f}" for x in neff) + " \\\\\n"
    else:
        missing.append("block_v2_maps_lcdm.npz")
    s += row("flat prior at 23 GHz", d["flat23"], val)
    s += row("flat prior at 30 GHz", d["flat30"], val)
    s += row(r"\quad minus profile, in units of the sandwich error", (d["flat30"] - d["beta"]) / d["sd_sand"], lambda x: f"${x:+.1f}$")
    s += row("zero-mean prior, free scale", d["zm"], val)
    s += "\\midrule\n"
    s += row("flat 30 GHz minus flat 23 GHz", d["flat30"] - d["flat23"], sgn)
    s += row(r"\quad Eq.~(\ref{eq:coupled})", d["pred_marg"], sgn)
    s += "\\midrule\n"
    s += row(r"$\chi^2$ per element and band", d["chi"], lambda x: f"{x:.2f}" if abs(x * 1000 % 10 - 5) > 1e-6 else f"{x:.3f}")
    s += row(r"$\Tm/\J_{\rm c}$", d["tj"], lambda x: f"{x:.3f}" if x >= 0.001 else f"{x:.4f}")
    return s + FOOT


def tab_cmbprior():
    R = {k: full(FILES[k]) for k in ("lcdm", "white", "flat")}
    if not all(v is not None for v in R.values()):
        return None
    s = HEAD
    s += "\\multicolumn{5}{l}{\\emph{profile likelihood, with the sandwich error}} \\\\\n"
    for k in R:
        s += LAB[k] + " & " + " & ".join(f"${b:.3f} \\pm {e:.3f}$" for b, e in zip(R[k]["beta"], R[k]["sd_sand"])) + " \\\\\n"
    s += "\\multicolumn{5}{l}{\\emph{flat amplitude prior at 23 GHz, minus profile}} \\\\\n"
    for k in R:
        s += row(LAB[k], R[k]["flat23"] - R[k]["beta"], sgn)
    s += "\\multicolumn{5}{l}{\\emph{flat amplitude prior at 30 GHz minus at 23 GHz: measured (Eq.~\\ref{eq:coupled})}} \\\\\n"
    for k in R:
        s += LAB[k] + " & " + " & ".join(f"${a:+.3f}$ (${b:+.3f}$)" for a, b in zip(R[k]["flat30"] - R[k]["flat23"], R[k]["pred_marg"])) + " \\\\\n"
    s += "\\multicolumn{5}{l}{\\emph{$\\chi^2$ per element and band (expectation 0.6, 0.6 and 0.4)}} \\\\\n"
    for k in R:
        s += row(LAB[k], R[k]["chi"], lambda x: f"{x:.3f}")
    return s + FOOT


def tab_moves():
    d = full(FILES["lcdm"])
    b0, s0 = d["beta"], d["sd_sand"]
    s = HEAD
    s += row(r"baseline: all scales, $\Lambda$CDM CMB prior", b0, val)
    s += row("statistical error (sandwich)", s0, lambda x: f"{x:.3f}")
    s += "\\midrule\n"
    s += row(r"flat amplitude prior at 30 GHz (Table~\ref{tab:e040})", d["flat30"] - b0, sgn)
    rows = [
        (r"CMB prior: white, 0.2 $\mu$K", os.path.join(E40, "bpxv_n16_wm.npz"), None),
        ("CMB prior: flat", os.path.join(E40, "s4_bpxv_n16_wm_cmbflat.npz"), None),
        (r"index from $\ell \le 30$", os.path.join(E42, "lowpass30_lcdm.npz"), "K_Ka_Q_V_353"),
        (r"index from $\ell \le 20$", os.path.join(E42, "lowpass20_lcdm.npz"), "K_Ka_Q_V_353"),
        (r"index from $\ell \le 12$", os.path.join(E42, "lowpass12_lcdm.npz"), "K_Ka_Q_V_353"),
        ("K and Ka only", os.path.join(E40, "l_bpxv_n16_KKa_wm.npz"), None),
        ("dust index 1.48 in place of 1.53", os.path.join(E40, "l_bpxv_n16_bd148_wm.npz"), None),
        ("dust index 1.58", os.path.join(E40, "l_bpxv_n16_bd158_wm.npz"), None),
        ("no bandpass drift correction", os.path.join(E40, "l_bpxv_n16_nodrift_wm.npz"), None),
        ("353 GHz: all bolometers", os.path.join(E40, "l_bpxv_n16_bp353full_wm.npz"), None),
        ("delta bandpasses at the nominal frequencies", os.path.join(E40, "l_bpxv_n16_delta_wm.npz"), None),
    ]
    notes = []
    # window matching: the LCDM prior is built on the window operator, so this arm exists with the white prior only
    w0, w1 = prof(os.path.join(E40, "bpxv_n16.npz")), prof(os.path.join(E40, "bpxv_n16_wm.npz"))
    if w0 is not None and w1 is not None:
        s += row("no window matching (white CMB prior, against its own baseline)", w0[0] - w1[0], sgn)
    for label, path, key in rows:
        r = prof(path, key)
        if r is None:
            continue
        s += row(label, r[0] - b0, sgn)
        notes.append(f"% {label}: value " + " ".join(f"{x:.4f}" for x in r[0]) + "; sd " + " ".join(f"{x:.4f}" for x in r[1])
                     + "; move in baseline sd " + " ".join(f"{x / e:+.1f}" for x, e in zip(r[0] - b0, s0)))
    lp = os.path.join(E42, "lowpass20_lcdm.npz")
    if os.path.exists(lp):
        a, b = prof(lp, "K_Ka_353"), prof(lp, "K_Ka_Q_V_353")
        dlt = a[0] - b[0]
        s += row(r"K and Ka only, both from $\ell \le 20$ (against the row above for $\ell \le 20$)", dlt, sgn)
        notes.append("% K, Ka against four bands at l <= 20, in units of sqrt(sd_sub^2 - sd_full^2): "
                     + " ".join(f"{x / np.sqrt(max(u * u - v * v, 1e-30)):+.1f}" for x, u, v in zip(dlt, a[1], b[1])))
    st = prof(os.path.join(E40, "l_bpxv_n16_start27_wm.npz"))
    if st is not None:
        notes.append("% start at -2.7: change " + " ".join(f"{x:+.5f}" for x in st[0] - b0))
    return s + FOOT + "\n".join(notes) + "\n"


def main():
    for name, fn in (("gen_tab_e040.tex", tab_e040), ("gen_tab_cmbprior.tex", tab_cmbprior), ("gen_tab_moves.tex", tab_moves)):
        s = fn()
        if s is None:
            print(name, "not written (inputs missing)")
            continue
        with open(os.path.join(OUT, name), "w") as f:
            f.write("% generated by nextgen/scripts/paperA1_make_tables.py; do not edit\n" + s)
        print("wrote", name)


def moves_text():
    """The prose of Sect. `sec:moves` with every number computed from the result files. Returns None until the
    low-pass arms with the LCDM prior exist."""
    base = full(FILES["lcdm"])
    b0, s0 = base["beta"], base["sd_sand"]
    lp = {}
    for L in (30, 20, 12):
        path = os.path.join(E42, f"lowpass{L}_lcdm.npz")
        if not os.path.exists(path):
            return None
        z = np.load(path, allow_pickle=True)
        lp[L] = {k: (np.asarray(z[k], float)[:4], np.asarray(z[k], float)[4:]) for k in z.files}
    four = {L: lp[L]["K_Ka_Q_V_353"] for L in lp}
    pair = {L: lp[L]["K_Ka_353"] for L in lp}
    d = {L: four[L][0] - b0 for L in lp}
    nest = {L: np.sqrt(np.clip(four[L][1] ** 2 - s0 ** 2, 1e-30, None)) for L in lp}
    sig20 = d[20] / nest[20]
    # stability of the three cuts: largest pairwise difference in units of the larger error of the pair
    stab = []
    for r in range(4):
        m = 0.0
        for a_, b_ in ((30, 20), (20, 12), (30, 12)):
            m = max(m, abs(four[a_][0][r] - four[b_][0][r]) / np.sqrt(abs(four[b_][1][r] ** 2 - four[a_][1][r] ** 2) + 1e-30))
        stab.append(m)
    kka = prof(os.path.join(E40, "l_bpxv_n16_KKa_wm.npz"))
    dk = kka[0] - b0
    dk_sig = dk / np.sqrt(kka[1] ** 2 - s0 ** 2)
    dk20 = pair[20][0] - four[20][0]
    dk20_sig = dk20 / np.sqrt(pair[20][1] ** 2 - four[20][1] ** 2)
    bd = [prof(os.path.join(E40, f"l_bpxv_n16_bd{x}_wm.npz"))[0] - b0 for x in ("148", "158")]
    slope = (bd[1] - bd[0]) / 0.10
    drift = prof(os.path.join(E40, "l_bpxv_n16_nodrift_wm.npz"))[0] - b0
    b353 = prof(os.path.join(E40, "l_bpxv_n16_bp353full_wm.npz"))[0] - b0
    delta = prof(os.path.join(E40, "l_bpxv_n16_delta_wm.npz"))[0] - b0
    white = prof(os.path.join(E40, "bpxv_n16_wm.npz"))[0] - b0
    flat = prof(os.path.join(E40, "s4_bpxv_n16_wm_cmbflat.npz"))[0] - b0

    def rng(v, nd=3, signed=False):
        lo, hi = f"{float(np.min(np.abs(v))):.{nd}f}", f"{float(np.max(np.abs(v))):.{nd}f}"
        return lo if lo == hi else f"{lo} to {hi}"

    def lst(v, nd=3, signed=True):
        f = (lambda x: f"${x:+.{nd}f}$") if signed else (lambda x: f"{x:.{nd}f}")
        return ", ".join(f(x) for x in v[:-1]) + " and " + f(v[-1])

    def lst1(v):
        return ", ".join(f"{abs(x):.1f}" for x in v[:-1]) + f" and {abs(v[-1]):.1f}"

    s = rf"""Table~\ref{{tab:moves}} lists what else moves the answer in the central configuration. The rows are
changes of the index; the statistical error of the baseline is given for scale.

\paragraph{{Angular scale.}} To ask at which scales two bands disagree we first marginalized the modes of
one map at multipoles 2 to $L$ and refitted, with the same cut in the K, Ka fit as a control fixed in
advance, where nothing should move. The control failed: the K, Ka index moved toward flat in every
region, by $+0.08$, $+0.16$, $+0.11$ and $+0.37$ at $L = 20$ (white CMB prior; \exref{{E042}}). Read from
the other side, with the pixel scale excluded and the index informed by multipoles up to 30, 20 or 12
only, the four-band index steepens (Table~\ref{{tab:moves}}): by {rng(np.array([d[30][0], d[20][0], d[12][0]]))} at High
Latitude, {rng(np.array([d[30][1], d[20][1], d[12][1]]))} in the Spur, {rng(np.array([d[30][2], d[20][2], d[12][2]]))} in the Galactic Center
and {rng(np.array([d[30][3], d[20][3], d[12][3]]))} in the Plane. The estimator is sound: on simulations with known indices it
is unbiased and calibrated (standardized errors of standard deviation 0.95 to 1.03 over 100
realizations at $N_{{\rm side}}=8$), an independent review re-derived it and confirmed it on its own
simulations, and the difference between the large-scale and the all-scale estimate of one simulated
sky has zero mean and a spread of 0.9 to 1.1 times $(\sigma_{{\rm low}}^2 - \sigma_{{\rm all}}^2)^{{1/2}}$
(\exref{{E039}}). On the maps the differences at $L = 20$ are {lst1(sig20)} times that number. The
shifts are the same within 0.025 with the white CMB prior. So the dependence on scale is a property of
these maps under this model, in every region, and not of the estimator or of the CMB prior. An index
read from the three cuts is not stable either: the largest difference between two cuts is
{lst1(np.array(stab))} times the standard deviation expected for nested fits.

We do not know the cause. One candidate is instrumental: the maps are pixel averages of skies that are
not band-limited, structure below the pixel scale folds into it through each band's own beam, K band
has the widest, and the window operator works to multipole 47 with Gaussian beams. Against it, the one
simulation we have of sub-pixel structure leaves a residual of 0.002, though that sky had no compact
sources. The other candidate is the sky: compact structure in the plane may be flatter than the diffuse
emission, and an index that varies from patch to patch (previous subsection) makes an index from
large scales and an index from all scales weight the sky differently. \citet{{fuskeland2014}}, with maps
smoothed to a common $1^\circ$ beam before any pixel averaging, find $-2.98 \pm 0.01$ in the plane,
close to our all-scale K, Ka value there (${kka[0][3]:.2f}$) and not to the large-scale one (${pair[20][0][3]:.2f}$ at
$L = 20$); that favours the sky. The rows of Table~\ref{{tab:moves}} are sizes of an effect and not
corrections.

\paragraph{{Band subset.}} K and Ka alone give indices that differ from the four bands together by
{lst(dk)}. The two fits share K and Ka, so the standard deviation of their difference
under the model is $(\sigma_{{\rm sub}}^2 - \sigma_{{\rm full}}^2)^{{1/2}}$, and the differences are
{", ".join(f"${x:+.1f}$" for x in dk_sig[:-1])} and ${dk_sig[-1]:+.1f}$ of it; with both fits read from multipoles up to 20 they are
{", ".join(f"${x:+.1f}$" for x in dk20_sig[:-1])} and ${dk20_sig[-1]:+.1f}$. With the white CMB prior the shifts were larger and all of one
sign ($+0.032$, $+0.023$, $+0.041$ and $+0.015$). What remains with the $\Lambda$CDM prior is a flatter K
to Ka index at High Latitude and in the Galactic Center, where the model does not fit in any case. A
spectrum that steepens with frequency would do this, and so would a dust model that predicts too much
polarized dust at 41 and 61 GHz from 353 GHz.

\paragraph{{Dust.}} The dust spectrum is fixed. Changing its index by 0.05 either way moves the
synchrotron index by {lst(0.5 * (np.abs(bd[0]) + np.abs(bd[1])), 3, False)}, linearly ($\partial\beta_{{\rm s}}/\partial\beta_{{\rm d}} =
{slope[0]:.2f}$, ${slope[1]:.2f}$, ${slope[2]:.2f}$ and ${slope[3]:.2f}$): one statistical standard deviation at High Latitude
and two in the Plane, more than we had written down in advance for three of the four regions. With the
uncertainty of the Planck value, 0.02,
it is 0.005 at High Latitude and 0.010 in the Plane.

\paragraph{{Bandpasses.}} Leaving out the correction for the drift of the WMAP band centres changes the
index by {rng(drift)}, and taking the 353 GHz band of all bolometers in place of the
polarization-sensitive ones by {rng(b353)}. Replacing the measured bandpasses by delta functions at
the nominal WMAP frequencies steepens it by {rng(delta)}, {rng(delta / s0, 1)} statistical standard
deviations: the effective frequency of K band for a spectrum this steep is 22.4 GHz, not 22.8.

\paragraph{{Poorly measured modes.}} The released K-band covariance has 20 eigenmodes with more than
ten times the median variance. Projecting all of them out moves the K, Ka index by 0.002 to 0.015
(white CMB prior). This tests the modes the matrix knows about, not the claim that it understates
them.

Taken together: for these maps and regions the profile index of the simple model is $-3.08$ at High
Latitude, in the Spur and in the Galactic Center and $-2.96$ in the Plane. Its error is 0.01 to 0.02
if each region has one index and 0.04 to 0.10 from the scatter between patches. It moves by
{rng(np.concatenate([white, flat]), 2)} with the prior on the CMB, by {rng(np.concatenate([d[30], d[20], d[12]]), 2)} with the angular scales that
inform it, by up to {np.max(np.abs(dk)):.2f} with the band subset, by {rng(0.5 * (np.abs(bd[0]) + np.abs(bd[1])))} for 0.05 in the dust index, and by
{rng(delta, 2)} if the bands are taken at their nominal frequencies. A number with an error bar that
means what it says needs a model in which a region is more than one index and the scale dependence
understood. That is a separate piece of work, and we do not quote an index here.
"""
    return s


def write_moves_text():
    s = moves_text()
    if s is None:
        print("gen_moves_text.tex not written (low-pass arms with the LCDM prior missing)")
        return
    with open(os.path.join(OUT, "gen_moves_text.tex"), "w") as f:
        f.write("% generated by nextgen/scripts/paperA1_make_tables.py (prose template there; numbers from the result files)\n" + s)
    print("wrote gen_moves_text.tex")


if __name__ == "__main__":
    main()
    write_moves_text()
    if missing:
        print("missing inputs (rows left out):", sorted(set(missing)))
