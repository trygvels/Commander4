"""Paper A1: the rows of Table tab:moves (what moves the profile index besides the estimator), printed from
the result files of E040 and E042 so that no number is typed by hand. Central configuration: the LCDM CMB prior.

    python -I nextgen/scripts/paperA1_tab_moves.py
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E40 = os.path.join(HERE, "..", "experiments", "E040")
E42 = os.path.join(HERE, "..", "experiments", "E042")


def load(path, key=None):
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if key is not None:
        if key not in z.files:
            return None
        v = np.asarray(z[key], dtype=float)
        return v[:4], v[4:]
    return np.asarray(z["beta"], dtype=float), np.asarray(z["sd_sand"], dtype=float)


def fmt(v, signed=True):
    return " & ".join((f"${x:+.3f}$" if signed else f"${x:.3f}$") for x in v)


def main():
    base = load(os.path.join(E40, "s4_bpxv_n16_wm_lcdm.npz"))
    b0, s0 = base
    print(r"baseline, all scales, $\Lambda$CDM CMB prior & " + fmt(b0, False) + r" \\")
    print(r"statistical error & " + " & ".join(f"{x:.3f}" for x in s0) + r" \\")
    z = np.load(os.path.join(E40, "s4_bpxv_n16_wm_lcdm.npz"), allow_pickle=True)
    print(r"flat amplitude prior at 30 GHz (Table~\ref{tab:e040}) & " + fmt(z["flat30"] - z["beta"]) + r" \\")
    rows = [
        (r"CMB prior: white, 0.2 $\mu$K", os.path.join(E40, "bpxv_n16_wm.npz"), None),
        (r"CMB prior: flat", os.path.join(E40, "s4_bpxv_n16_wm_cmbflat.npz"), None),
        ("no window matching", os.path.join(E40, "l_bpxv_n16_nowm.npz"), None),
        (r"index from $\ell \le 30$", os.path.join(E42, "lowpass30_lcdm.npz"), "K_Ka_Q_V_353"),
        (r"index from $\ell \le 20$", os.path.join(E42, "lowpass20_lcdm.npz"), "K_Ka_Q_V_353"),
        (r"index from $\ell \le 12$", os.path.join(E42, "lowpass12_lcdm.npz"), "K_Ka_Q_V_353"),
        ("K and Ka only", os.path.join(E40, "l_bpxv_n16_KKa_wm.npz"), None),
        (r"K and Ka only, both from $\ell \le 20$", (os.path.join(E42, "lowpass20_lcdm.npz"), "K_Ka_353", "K_Ka_Q_V_353"), "pair"),
        ("dust index 1.48", os.path.join(E40, "l_bpxv_n16_bd148_wm.npz"), None),
        ("dust index 1.58", os.path.join(E40, "l_bpxv_n16_bd158_wm.npz"), None),
        ("no bandpass drift correction", os.path.join(E40, "l_bpxv_n16_nodrift_wm.npz"), None),
        ("353 GHz: all bolometers", os.path.join(E40, "l_bpxv_n16_bp353full_wm.npz"), None),
        ("delta bandpasses at nominal frequencies", os.path.join(E40, "l_bpxv_n16_delta_wm.npz"), None),
        ("optimizer started at $-2.7$", os.path.join(E40, "l_bpxv_n16_start27_wm.npz"), None),
    ]
    for label, path, key in rows:
        if key == "pair":
            a, b = load(path[0], path[1]), load(path[0], path[2])
            if a is None or b is None:
                print(f"% {label}: not available yet")
                continue
            d = a[0] - b[0]
            sd = np.sqrt(np.clip(a[1] ** 2 - b[1] ** 2, 0, None))
            print(f"{label} & " + fmt(d) + r" \\" + "   % in units of sqrt(sd_sub^2 - sd_full^2): " + " ".join(f"{x / s:+.1f}" for x, s in zip(d, sd)))
            continue
        r = load(path, key)
        if r is None:
            print(f"% {label}: not available yet ({os.path.basename(path)})")
            continue
        d = r[0] - b0
        extra = "   % value " + " ".join(f"{x:.4f}" for x in r[0]) + "; sd " + " ".join(f"{x:.4f}" for x in r[1]) + "; move in baseline sd " + " ".join(f"{x / s:+.1f}" for x, s in zip(d, s0))
        print(f"{label} & " + fmt(d) + r" \\" + extra)
    for tag in ("lcdm", "white"):
        p = os.path.join(E40, f"block_v2_maps_{tag}.npz")
        if os.path.exists(p):
            zz = np.load(p)
            print(f"% block errors, CMB prior {tag}: n_eff " + " ".join(f"{x:.1f}" for x in zz["neff"]) + "; sandwich " + " ".join(f"{x:.4f}" for x in zz["sd_model"])
                  + "; CR0 " + " ".join(f"{x:.4f}" for x in zz["CR0"]) + "; CR2 " + " ".join(f"{x:.4f}" for x in zz["CR2"]) + "; CR3 " + " ".join(f"{x:.4f}" for x in zz["CR3"]))
        else:
            print(f"% block errors, CMB prior {tag}: not available yet")


if __name__ == "__main__":
    main()
