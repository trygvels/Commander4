"""E041: coverage ladder for one region index with error estimates taken from each data set.

    python nextgen/scripts/e041_coverage_ladder.py --out FILE.npz [--n 32 128 512 2048] [--s 0.5 1 2 4 8 16]
        [--nreal 400] [--common-noise]

One power-law component, five delta bands, diagonal noise. See the E041 record for targets and predictions.
"""
import argparse
import time

import numpy as np
from scipy.optimize import minimize, minimize_scalar

NU = np.array([22.8, 28.4, 33.0, 40.6, 44.1])
SIG = np.array([1.93, 1.15, 1.95, 1.77, 1.38])
LO, HI = -4.6, -1.6
TRUTH = -2.9


class Cell:
    def __init__(self, a, w, d, nu_ref):
        self.a, self.w, self.d, self.ln = a, w, d, np.log(NU / nu_ref)

    def stats(self, beta):
        f = np.exp(beta * self.ln)                           # (nb,)
        wf = self.w * f[None, :]                             # (n, nb)
        A = (wf * f[None, :]).sum(1)
        b = (wf * self.d).sum(1)
        return f, wf, A, b

    def prof(self, beta):
        _, _, A, b = self.stats(beta)
        return -0.5 * np.sum(b * b / A)

    def flat(self, beta):
        _, _, A, b = self.stats(beta)
        return -0.5 * np.sum(b * b / A) + 0.5 * np.sum(np.log(A))

    def zm(self, x):
        _, _, A, b = self.stats(x[0])
        t2 = np.exp(2 * x[1])
        return -0.5 * np.sum(b * b * t2 / (1 + t2 * A)) + 0.5 * np.sum(np.log1p(t2 * A))

    def info(self, beta, amp=None):
        """J (with u from `amp`, or from the fitted amplitudes), T, and the gradient of the volume term at
        this reference: sum_p ln(nu_ref / nu_p) = -sum_p <ln(nu / nu_ref)>_p."""
        f, wf, A, b = self.stats(beta)
        m1 = (wf * f[None, :] * self.ln[None, :]).sum(1) / A
        V = (wf * f[None, :] * self.ln[None, :] ** 2).sum(1) / A - m1 ** 2
        u = (b * b / A) if amp is None else amp ** 2 * A
        return np.sum(u * V), np.sum(V), -np.sum(m1)


def d2(fun, x, h):
    return (fun(x + h) - 2 * fun(x) + fun(x - h)) / h ** 2


def opt1(fun, x0, width):
    """Bounded scalar optimum: bracket around x0, falling back to the whole range."""
    lo, hi = max(LO, x0 - 8 * width), min(HI, x0 + 8 * width)
    r = minimize_scalar(fun, bounds=(lo, hi), method="bounded", options={"xatol": 1e-7})
    if r.x - lo < 1e-5 or hi - r.x < 1e-5:
        r = minimize_scalar(fun, bounds=(LO, HI), method="bounded", options={"xatol": 1e-7})
    return r.x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=[32, 128, 512, 2048])
    ap.add_argument("--s", type=float, nargs="+", default=[0.5, 1, 2, 4, 8, 16])
    ap.add_argument("--nreal", type=int, default=400)
    ap.add_argument("--common-noise", action="store_true")
    ap.add_argument("--seed", type=int, default=41)
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    t0 = time.time()
    cols = ["n", "s", "snr_med", "sd_sand_truth", "t_truth", "pivot",
            "p_bias", "p_edge", "p_mz", "p_sdz_curv", "p_sdz_plain", "p_sdz_deb", "p_cov_curv", "p_cov_plain", "p_cov_deb",
            "p_cov95_deb", "p_nodeb",
            "f23_shift", "f23_pred", "f30_shift", "f30_pred", "fpiv_shift", "f23_mz", "f30_mz", "f23_cov", "f30_cov",
            "z_shift", "z_mz", "z_sdz", "z_cov", "z_cov95", "z_minus_p_rms"]
    rows = []
    print("# " + " ".join(cols))
    for n in a.n:
        rng0 = np.random.default_rng([a.seed, n])
        amp = 8.4 * np.exp(1.2 * rng0.standard_normal(n)) * rng0.choice([-1.0, 1.0], n)
        hfac = np.ones((n, NU.size)) if a.common_noise else np.exp(0.2 * rng0.standard_normal((n, NU.size)))
        for s in a.s:
            sig = SIG[None, :] * s * hfac
            w = 1.0 / sig ** 2
            mean = amp[:, None] * (NU[None, :] / 23.0) ** TRUTH
            c0 = Cell(amp, w, mean, 23.0)
            J0, T0, g23 = c0.info(TRUTH, amp=amp)
            g30 = g23 + n * np.log(30.0 / 23.0)
            piv = 23.0 * np.exp(-g23 / n)                    # reference at which the volume gradient vanishes
            sd_truth = np.sqrt((J0 + T0)) / J0
            width = max(sd_truth, 1e-4)
            rng = np.random.default_rng([a.seed, n, int(round(100 * s))])
            R = {k: [] for k in ("p", "sc", "sp", "sd", "f23", "s23", "f30", "s30", "fpiv", "z", "sz")}
            for _ in range(a.nreal):
                d = mean + sig * rng.standard_normal((n, NU.size))
                c = Cell(amp, w, d, 23.0)
                bp = opt1(c.prof, TRUTH, width)
                h = min(0.02, max(2e-4, 0.3 * width))
                cur = d2(c.prof, bp, h)
                Jp, T, _ = c.info(bp)
                Jc = Jp - T
                R["p"].append(bp)
                R["sc"].append(1 / np.sqrt(cur) if cur > 0 else np.inf)
                R["sp"].append(np.sqrt(Jp + T) / Jp)
                R["sd"].append(np.sqrt(Jp) / Jc if Jc > 0 else np.inf)
                for key, ref in (("23", 23.0), ("30", 30.0)):
                    cf = Cell(amp, w, d, ref)
                    bf = opt1(cf.flat, bp, 6 * width + abs(n * np.log(ref / piv) / J0))
                    cu = d2(cf.flat, bf, h)
                    R["f" + key].append(bf)
                    R["s" + key].append(1 / np.sqrt(cu) if cu > 0 else np.inf)
                cf = Cell(amp, w, d, piv)
                R["fpiv"].append(opt1(cf.flat, bp, 6 * width))
                # zero-mean prior, tau maximized jointly
                _, _, A, b = c.stats(bp)
                lt0 = 0.5 * np.log(max(np.mean((b / A) ** 2 - 1 / A), 1e-6))
                rz = minimize(c.zm, np.array([bp, lt0]), method="L-BFGS-B", bounds=[(LO, HI), (-8.0, 10.0)],
                              options={"ftol": 1e-13, "gtol": 1e-9})
                x = rz.x
                hb, ht = h, 0.02
                H = np.zeros((2, 2))
                e0, e1 = np.array([hb, 0.0]), np.array([0.0, ht])
                f0 = c.zm(x)
                H[0, 0] = (c.zm(x + e0) - 2 * f0 + c.zm(x - e0)) / hb ** 2
                H[1, 1] = (c.zm(x + e1) - 2 * f0 + c.zm(x - e1)) / ht ** 2
                H[0, 1] = H[1, 0] = (c.zm(x + e0 + e1) - c.zm(x + e0 - e1) - c.zm(x - e0 + e1) + c.zm(x - e0 - e1)) / (4 * hb * ht)
                det = H[0, 0] * H[1, 1] - H[0, 1] ** 2
                R["z"].append(x[0])
                R["sz"].append(np.sqrt(H[1, 1] / det) if (det > 0 and H[1, 1] > 0) else np.inf)
            R = {k: np.array(v) for k, v in R.items()}
            edge = np.mean((R["p"] - LO < 1e-4) | (HI - R["p"] < 1e-4))

            def zs(est, sd):
                z = (est - TRUTH) / sd
                ok = np.isfinite(z)
                return z, ok

            zc, okc = zs(R["p"], R["sc"])
            zp, okp = zs(R["p"], R["sp"])
            zd, okd = zs(R["p"], R["sd"])
            z23, ok23 = zs(R["f23"], R["s23"])
            z30, ok30 = zs(R["f30"], R["s30"])
            zz, okz = zs(R["z"], R["sz"])
            cov = lambda z, ok, k=1.0: float(np.mean(np.abs(z[ok]) < k)) if ok.any() else np.nan   # noqa: E731
            sdz = lambda z, ok: float(np.std(z[ok], ddof=1)) if ok.sum() > 2 else np.nan           # noqa: E731
            row = [n, s, np.median(np.abs(amp)) / (SIG[0] * s), sd_truth, T0 / J0, piv,
                   np.mean(R["p"]) - TRUTH, edge, np.mean(zd[okd]) if okd.any() else np.nan,
                   sdz(zc, okc), sdz(zp, okp), sdz(zd, okd), cov(zc, okc), cov(zp, okp), cov(zd, okd), cov(zd, okd, 2.0),
                   1 - np.mean(okd),
                   np.mean(R["f23"] - R["p"]), g23 / J0, np.mean(R["f30"] - R["p"]), g30 / J0, np.mean(R["fpiv"] - R["p"]),
                   np.mean(z23[ok23]), np.mean(z30[ok30]), cov(z23, ok23), cov(z30, ok30),
                   np.mean(R["z"]) - TRUTH, np.mean(zz[okz]) if okz.any() else np.nan, sdz(zz, okz), cov(zz, okz),
                   cov(zz, okz, 2.0), np.sqrt(np.mean((R["z"] - R["p"]) ** 2)) / sd_truth]
            rows.append(row)
            print(" ".join(f"{v:.4g}" for v in row), f"[{time.time() - t0:.0f} s]", flush=True)
    if a.out:
        np.savez(a.out, cols=np.array(cols), rows=np.array(rows))


if __name__ == "__main__":
    main()
