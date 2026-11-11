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
    s += row(r"\quad sandwich error", d["sd_sand"], lambda x: f"{x:.3f}")
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
        ("no window matching", os.path.join(E40, "l_bpxv_n16_nowm.npz"), None),
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
    if missing:
        print("missing inputs (rows left out):", sorted(set(missing)))


if __name__ == "__main__":
    main()
