"""Paper A1: the rows of the central table (Table tab:e040) and of the CMB-prior comparison, printed from the
E040 output files so that no number is typed by hand.

    python -I nextgen/scripts/paperA1_tab_e040.py
"""
import os
import re

import numpy as np

E = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "experiments", "E040")
NAMES = ["High Latitude", "Spur", "Galactic Center", "Plane"]
FILES = {"lcdm": "s4_bpxv_n16_wm_lcdm.txt", "white": "bpxv_n16_wm.txt", "flat": "s4_bpxv_n16_wm_cmbflat.txt"}
COLS = ["n_el", "rms", "profile", "sd_naive", "sd", "f23", "f30", "d3023", "predJ", "predH", "zm", "tau"]


def read(path):
    """Numbers from the .npz written beside the text output (full precision); the chi-square per region from the
    text file, which prints it to three decimals."""
    npz = path[:-4] + ".npz"
    if not (os.path.exists(path) and os.path.exists(npz)):
        return None
    z = np.load(npz, allow_pickle=True)
    out = {}
    for i, nm in enumerate(NAMES):
        d = {"profile": float(z["beta"][i]), "sd": float(z["sd_sand"][i])}
        if "flat23" in z.files:
            d.update(f23=float(z["flat23"][i] - z["beta"][i]), f30=float(z["flat30"][i] - z["beta"][i]),
                     d3023=float(z["flat30"][i] - z["flat23"][i]), predH=float(z["pred_marg"][i]))
        if "zm" in z.files:
            d["zm"] = float(z["zm"][i] - z["beta"][i])
        out[nm] = d
    J, T = z["J"], z["T"]
    tj = [float(T[k, k] / (J[k, k] - T[k, k])) for k in range(4)]
    chi = None
    for ln in open(path):
        if ln.startswith("# chi2 by region"):
            chi = [float(x) for x in re.findall(r"(?:High Latitude|Spur|Galactic Center|Plane) ([0-9.]+)", ln)]
    return out, chi, tj


def row(label, vals, fmt):
    return f"{label} & " + " & ".join(fmt(v) for v in vals) + r" \\"


def main():
    for key, title in (("lcdm", r"$\Lambda$CDM covariance"), ("white", r"white, 0.2 $\mu$K per element"), ("flat", "flat")):
        r = read(os.path.join(E, FILES[key]))
        print(f"% ---- CMB prior: {title}  ({FILES[key]})")
        if r is None:
            print("% not available yet")
            continue
        t, chi, tj = r
        g = lambda c: [t[n][c] for n in NAMES]  # noqa: E731
        sg = lambda x: f"${x:+.3f}$"  # noqa: E731
        print(row("profile likelihood", g("profile"), lambda x: f"${x:.3f}$"))
        print(row(r"\quad sandwich error", g("sd"), lambda x: f"{x:.3f}"))
        if "f23" in t[NAMES[0]]:
            print(row("flat prior at 23 GHz", [a + b for a, b in zip(g("profile"), g("f23"))], lambda x: f"${x:.3f}$"))
            print(row("flat prior at 30 GHz", [a + b for a, b in zip(g("profile"), g("f30"))], lambda x: f"${x:.3f}$"))
            print(row(r"\quad minus profile, in units of the error", [b / s for b, s in zip(g("f30"), g("sd"))], lambda x: f"${x:+.1f}$"))
            if "zm" in t[NAMES[0]]:
                print(row("zero-mean prior, free scale", [a + b for a, b in zip(g("profile"), g("zm"))], lambda x: f"${x:.3f}$"))
                print("%   zm - profile in sd: " + " ".join(f"{b / s:+.2f}" for b, s in zip(g("zm"), g("sd"))))
            print(row("flat 30 GHz minus flat 23 GHz", g("d3023"), sg))
            print(row(r"\quad Eq.~(\ref{eq:coupled})", g("predH"), sg))
            print("%   measured / predicted: " + " ".join(f"{a / b:.3f}" for a, b in zip(g("d3023"), g("predH"))))
            print("%   flat23 - profile in sd: " + " ".join(f"{b / s:+.2f}" for b, s in zip(g("f23"), g("sd"))))
        if chi:
            print(row(r"$\chi^2$ per element and band", chi, lambda x: f"{x:.3f}"))
        if tj:
            print(row(r"$\Tm/\J_{\rm c}$", tj, lambda x: f"{x:.3f}" if x >= 0.001 else f"{x:.4f}"))
    a, b = read(os.path.join(E, FILES["lcdm"])), read(os.path.join(E, FILES["white"]))
    if a and b:
        d = [a[0][n]["profile"] - b[0][n]["profile"] for n in NAMES]
        print("% LCDM minus white, profile: " + " ".join(f"{x:+.4f}" for x in d) + "; in LCDM sd: "
              + " ".join(f"{x / a[0][n]['sd']:+.2f}" for x, n in zip(d, NAMES)))
    c = read(os.path.join(E, FILES["flat"]))
    if a and c:
        d = [c[0][n]["profile"] - a[0][n]["profile"] for n in NAMES]
        print("% flat minus LCDM, profile: " + " ".join(f"{x:+.4f}" for x in d) + "; in flat sd: "
              + " ".join(f"{x / c[0][n]['sd']:+.2f}" for x, n in zip(d, NAMES)) + "; error ratio flat / LCDM: "
              + " ".join(f"{c[0][n]['sd'] / a[0][n]['sd']:.2f}" for n in NAMES))
    print("% ==== Table tab:cmbprior")
    R = {k: read(os.path.join(E, FILES[k])) for k in ("lcdm", "white", "flat")}
    lab = {"lcdm": r"$\Lambda$CDM covariance", "white": r"white, 0.2 $\mu$K", "flat": "flat"}
    if all(R.values()):
        print(r"\multicolumn{5}{l}{\emph{profile likelihood, with the sandwich error}} \\")
        for k in ("lcdm", "white", "flat"):
            print(lab[k] + " & " + " & ".join(f"${R[k][0][n]['profile']:.3f} \\pm {R[k][0][n]['sd']:.3f}$" for n in NAMES) + r" \\")
        print(r"\multicolumn{5}{l}{\emph{flat amplitude prior at 23 GHz, minus profile}} \\")
        for k in ("lcdm", "white", "flat"):
            print(lab[k] + " & " + " & ".join(f"${R[k][0][n]['f23']:+.3f}$" for n in NAMES) + r" \\")
        print(r"\multicolumn{5}{l}{\emph{flat amplitude prior at 30 GHz minus at 23 GHz: measured (Eq.~\ref{eq:coupled})}} \\")
        for k in ("lcdm", "white", "flat"):
            print(lab[k] + " & " + " & ".join(f"${R[k][0][n]['d3023']:+.3f}$ (${R[k][0][n]['predH']:+.3f}$)" for n in NAMES) + r" \\")
        print(r"\multicolumn{5}{l}{\emph{$\chi^2$ per element and band (expectation 0.6, 0.6 and 0.4)}} \\")
        for k in ("lcdm", "white", "flat"):
            print(lab[k] + " & " + " & ".join(f"{x:.3f}" for x in R[k][1]) + r" \\")


if __name__ == "__main__":
    main()
