#!/usr/bin/env python
"""E023 stage 1: exact-target block-shift scans of the collapsed hier index likelihood.

For one simulated sky (band maps, rms maps, templates and truth.npz from pysm_intensity_sim.py)
this evaluates, as a function of the index fields alone, the per-pixel marginal that Commander4's
`likelihood: hier` samples, with the three amplitudes (CMB, sync, dust) integrated out under
their Gaussian priors and the amplitude hyperparameters (mu_c, log s_c) at their joint MAP:

    L(beta) = max_theta sum_p [-1/2 r_p^T C_p^-1 r_p - 1/2 log det C_p] + log p(theta),
    C_p = M_p T_p M_p^T + N_p,  r_p = d_p - M_p m_p,
    m_p = (0, mu_s g_s,p, mu_d g_d,p),  T_p = diag(cmb_sd^2, (s_s,p g_s,p)^2, (s_d,p g_d,p)^2),

with hyperpriors mu_c ~ N(1, 0.5^2) and log s_c ~ N(log 0.3, 1^2) (index_field.HIER_*), log s
bounded below by HIER_LOGS_MIN. The index prior is left out on purpose (E023 record). The
scatter s_c,p is one global value per component ("global") or one per S/N quartile of the
nside-16 regions ("perq"). The dust SED is one MBB referenced at 545 GHz at the fixed T or at the
PySM T_d per pixel. Templates are the data-built pseudo-templates on disk or the true amplitudes.

Scans: the index fields are set to the nside-16 region means of the truth (beta_s; beta_eff3 for
dust at fixed T, beta_d at the true T), one S/N quartile of regions is shifted by delta, and the
profiled log-density is tabulated against delta with its chi-square and log-determinant parts;
plus a wide sync scan in the faintest quartile for the E020 runaway. The hyperparameters are
re-profiled at every delta (L-BFGS-B with analytic gradients, warm-started).
"""
import argparse, json, sys, time
import numpy as np
import healpy as hp
from scipy.optimize import minimize

H_K = 4.799243e-11      # h/k in K s
T_CMB = 2.7255
FREQS = np.array([23., 30., 33., 41., 44., 61., 70., 94., 100., 143., 217., 353.])
MU_PRIOR, LOGS_PRIOR, LOGS_MIN = (1.0, 0.5), (float(np.log(0.3)), 1.0), -8.0
CMB_SD = 100.0


def cmb_mixing(nu):
    x = H_K * nu * 1e9 / T_CMB
    return x**2 * np.exp(x) / np.expm1(x)**2          # uK_RJ per uK_CMB


def sync_sed(bs):
    return (FREQS[None, :] / 30.0) ** bs[:, None]


def dust_sed(bd, T):
    x = H_K * FREQS[None, :] * 1e9 / T[:, None]
    xr = H_K * 545e9 / T[:, None]
    return (FREQS[None, :] / 545.0) ** (bd[:, None] + 1.0) * np.expm1(xr) / np.expm1(x)


class Sky:
    def __init__(self, d, rms, gs, gd, T, quart_s, quart_d):
        self.d, self.rms2, self.gs, self.gd, self.T = d, rms**2, gs, gd, T
        self.qs, self.qd = quart_s, quart_d
        self.n = d.shape[0]
        self.cm = np.broadcast_to(cmb_mixing(FREQS), (self.n, 12))

    def unpack(self, theta, perq):
        if perq:
            mus, mud = theta[0], theta[5]
            ls_s, ls_d = theta[1:5][self.qs], theta[6:10][self.qd]
        else:
            mus, ls_s, mud, ls_d = theta[0], np.full(self.n, theta[1]), theta[2], np.full(self.n, theta[3])
        return mus, ls_s, mud, ls_d

    def hyperprior(self, theta, perq):
        mus = [theta[0], theta[5]] if perq else [theta[0], theta[2]]
        lss = list(theta[1:5]) + list(theta[6:10]) if perq else [theta[1], theta[3]]
        lp = sum(-0.5 * ((m - MU_PRIOR[0]) / MU_PRIOR[1])**2 for m in mus)
        lp += sum(-0.5 * ((l - LOGS_PRIOR[0]) / LOGS_PRIOR[1])**2 for l in lss)
        g = np.zeros_like(theta)
        if perq:
            g[0], g[5] = [-(theta[i] - MU_PRIOR[0]) / MU_PRIOR[1]**2 for i in (0, 5)]
            g[1:5] = -(theta[1:5] - LOGS_PRIOR[0]) / LOGS_PRIOR[1]**2
            g[6:10] = -(theta[6:10] - LOGS_PRIOR[0]) / LOGS_PRIOR[1]**2
        else:
            g[0], g[2] = [-(theta[i] - MU_PRIOR[0]) / MU_PRIOR[1]**2 for i in (0, 2)]
            g[1], g[3] = [-(theta[i] - LOGS_PRIOR[0]) / LOGS_PRIOR[1]**2 for i in (1, 3)]
        return lp, g

    def terms(self, Ms, Md, theta, perq, grad=False):
        """(chi2 part, logdet part) of sum_p log N(d_p; M_p m_p, C_p), and the gradient wrt theta."""
        mus, ls_s, mud, ls_d = self.unpack(theta, perq)
        M = np.stack([self.cm, Ms, Md], -1)                                  # (n, 12, 3)
        m = np.stack([np.zeros(self.n), mus * self.gs, mud * self.gd], -1)
        Tv = np.stack([np.full(self.n, CMB_SD**2), (np.exp(ls_s) * self.gs)**2,
                       (np.exp(ls_d) * self.gd)**2], -1)
        C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + self.rms2[:, None, :] * np.eye(12)[None]
        r = self.d - np.einsum("nij,nj->ni", M, m)
        L = np.linalg.cholesky(C)
        sol = np.linalg.solve(C, r[..., None])[..., 0]
        chi = -0.5 * np.einsum("ni,ni->n", r, sol)
        ld = -np.log(np.einsum("nii->ni", L)).sum(1)
        if not grad:
            return chi, ld, None
        # d/dmu_c: r depends on mu_c through m; d/dlog s_c: C depends through Tv_c (dTv = 2 Tv).
        u = np.einsum("nij,ni->nj", M, sol)                                   # M^T C^-1 r, (n, 3)
        CiM = np.linalg.solve(C, M)                                            # (n, 12, 3)
        w = np.einsum("nij,nij->nj", M, CiM)                                   # diag of M^T C^-1 M
        dmu_s, dmu_d = u[:, 1] * self.gs, u[:, 2] * self.gd                   # per pixel
        dls_s = Tv[:, 1] * (u[:, 1]**2 - w[:, 1])                              # chi + logdet parts
        dls_d = Tv[:, 2] * (u[:, 2]**2 - w[:, 2])
        g = np.zeros_like(theta)
        if perq:
            g[0], g[5] = dmu_s.sum(), dmu_d.sum()
            g[1:5] = np.bincount(self.qs, dls_s, 4)
            g[6:10] = np.bincount(self.qd, dls_d, 4)
        else:
            g[:] = [dmu_s.sum(), dls_s.sum(), dmu_d.sum(), dls_d.sum()]
        return chi, ld, g

    def profile(self, Ms, Md, theta0, perq):
        def f(th):
            chi, ld, g = self.terms(Ms, Md, th, perq, grad=True)
            lp, gp = self.hyperprior(th, perq)
            return -(chi.sum() + ld.sum() + lp), -(g + gp)
        nls = 8 if perq else 2
        bounds = [(None, None)] + [(LOGS_MIN, 3.0)] * (4 if perq else 1) + [(None, None)] + \
                 [(LOGS_MIN, 3.0)] * (4 if perq else 1)
        res = minimize(f, theta0, jac=True, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 200, "ftol": 1e-12, "gtol": 1e-4})
        chi, ld, _ = self.terms(Ms, Md, res.x, perq)
        lp, _ = self.hyperprior(res.x, perq)
        return res.x, float(chi.sum()), float(ld.sum()), float(lp), res.nfev


def region_means(x64, nside_r=16):
    return hp.ud_grade(x64, nside_r, order_in="RING", order_out="RING", power=0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells", nargs="+", default=["data_fixT_global", "data_trueT_global",
                                                    "exact_fixT_global", "exact_trueT_global",
                                                    "data_fixT_perq", "data_trueT_perq"])
    ap.add_argument("--fine", type=float, nargs=3, default=[-0.10, 0.10, 0.01])
    ap.add_argument("--wide", type=float, nargs=3, default=[-3.0, 1.0, 0.1])
    ap.add_argument("--quick", action="store_true", help="one cell, three deltas (timing)")
    a = ap.parse_args()
    D = a.dir
    t = np.load(f"{D}/truth.npz")
    nside = int(t["nside"]); npix = hp.nside2npix(nside)
    d = np.array([hp.read_map(f"{D}/band{int(f):03d}_map.fits", field=0) for f in FREQS]).T
    rms = np.array([hp.read_map(f"{D}/band{int(f):03d}_rms.fits", field=0) for f in FREQS]).T
    g_data = {"sync": np.maximum(hp.read_map(f"{D}/template_sync_30GHz.fits", field=0), 1e-3),
              "dust": np.maximum(hp.read_map(f"{D}/template_dust_545GHz.fits", field=0), 1e-3)}
    a_true = {"sync": np.asarray(t["a_sync"], float), "dust": np.asarray(t["a_dust"], float)}
    g_exact = {k: np.maximum(v, 1e-3 * np.mean(np.abs(v))) for k, v in a_true.items()}
    T_fixed = float(t["T_fixed"]); T_true = np.asarray(t["T_d"], float)
    # nside-16 regions, S/N quartiles from the profiled Fisher sds (as in the E021 analyzer)
    reg = hp.ud_grade(np.arange(hp.nside2npix(16), dtype=float), nside).astype(int)
    fs16, fd16 = np.asarray(t["fisher_sd_sync_16_profiled"]), np.asarray(t["fisher_sd_dust_16_profiled"])
    def quartiles(f16):
        q = np.empty(f16.size, int); order = np.argsort(f16); n4 = f16.size // 4
        for k in range(4): q[order[k * n4:(k + 1) * n4 if k < 3 else None]] = k
        return q                                                               # 0 = highest S/N
    qs16, qd16 = quartiles(fs16), quartiles(fd16)
    qs64, qd64 = qs16[reg], qd16[reg]
    truth16 = {"s": region_means(np.asarray(t["beta_s"], float)),
               "d_fixT": region_means(np.asarray(t["beta_eff3"], float)),
               "d_trueT": region_means(np.asarray(t["beta_d"], float))}
    sky = Sky(d, rms, None, None, None, qs64, qd64)
    out = open(a.out, "w")
    def say(s):
        print(s); out.write(s + "\n"); out.flush()
    say(f"E023 exact-target scans on {D}: nside {nside}, noise scale {float(t['noise_scale'])}, "
        f"T_fixed {T_fixed:.3f} K; truth sd beta_s {np.asarray(t['beta_s']).std():.4f}, "
        f"beta_eff3 {np.asarray(t['beta_eff3']).std():.4f}, beta_d {np.asarray(t['beta_d']).std():.4f}; "
        f"region-mean truth sd (nside 16) s {truth16['s'].std():.4f}, d_fixT {truth16['d_fixT'].std():.4f}, "
        f"d_trueT {truth16['d_trueT'].std():.4f}")
    for name in ("sync", "dust"):
        hi = np.abs(hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]) > 20
        r = g_data[name][hi] / a_true[name][hi]
        say(f"  data template {name}: g/true at |b|>20 16/50/84 pct {np.percentile(r,16):.3f}/"
            f"{np.percentile(r,50):.3f}/{np.percentile(r,84):.3f}")
    summary = {}
    fine = np.arange(a.fine[0], a.fine[1] + 1e-9, a.fine[2])
    wide = np.arange(a.wide[0], a.wide[1] + 1e-9, a.wide[2])
    if a.quick:
        a.cells, fine, wide = a.cells[:1], np.array([-0.01, 0.0, 0.01]), np.array([])
    for cell in a.cells:
        tmpl, sed, prior = cell.split("_")
        sky.gs = (g_data if tmpl == "data" else g_exact)["sync"]
        sky.gd = (g_data if tmpl == "data" else g_exact)["dust"]
        T = np.full(npix, T_fixed) if sed == "fixT" else T_true
        perq = prior == "perq"
        bs0 = truth16["s"][reg]; bd0 = truth16["d_fixT" if sed == "fixT" else "d_trueT"][reg]
        theta0 = np.array([1.0, LOGS_PRIOR[0], 1.0, LOGS_PRIOR[0]]) if not perq else \
            np.array([1.0] + [LOGS_PRIOR[0]] * 4 + [1.0] + [LOGS_PRIOR[0]] * 4)
        t0 = time.time()
        th0, chi0, ld0, lp0, nfev = sky.profile(sync_sed(bs0), dust_sed(bd0, T), theta0, perq)
        tot0 = chi0 + ld0 + lp0
        def fmt(th):
            if perq:
                return (f"mu_s {th[0]:.4f} s_s {np.exp(th[1:5]).round(3).tolist()} mu_d {th[5]:.4f} "
                        f"s_d {np.exp(th[6:10]).round(3).tolist()}")
            return f"mu_s {th[0]:.4f} s_s {np.exp(th[1]):.4f} mu_d {th[2]:.4f} s_d {np.exp(th[3]):.4f}"
        say(f"\n== cell {cell}: profiled at the region-mean truth ({nfev} evals, {time.time()-t0:.1f} s): "
            f"{fmt(th0)}; chi2 part {chi0:.1f}, logdet part {ld0:.1f}, hyperprior {lp0:.2f}")
        summary[cell] = {"theta0": th0.tolist(), "scans": {}}
        scans = [(f"{c}_Q{q+1}", c, q, fine) for c in ("s", "d") for q in range(4)]
        if wide.size:
            scans.append(("s_Q4_wide", "s", 3, wide))
        for label, comp, q, deltas in scans:
            sel = (qs64 if comp == "s" else qd64) == q
            th = th0.copy(); rows = []
            t1 = time.time()
            for dl in deltas:
                bs = bs0 + dl * sel if comp == "s" else bs0
                bd = bd0 + dl * sel if comp == "d" else bd0
                th, chi, ld, lp, nfev = sky.profile(sync_sed(bs), dust_sed(bd, T), th, perq)
                rows.append((float(dl), chi - chi0, ld - ld0, lp - lp0, chi + ld + lp - tot0, th.copy()))
            tot = np.array([r[4] for r in rows]); i = int(np.argmax(tot))
            sd = float("nan")
            if 0 < i < len(rows) - 1:
                h = deltas[1] - deltas[0]; curv = (tot[i + 1] - 2 * tot[i] + tot[i - 1]) / h**2
                sd = float(np.sqrt(-1.0 / curv)) if curv < 0 else float("nan")
            say(f"-- {cell} {label}: {sel.sum()} pixels; argmax delta {deltas[i]:+.3f} (curvature sd "
                f"{sd:.4f}), gain over delta 0 {tot[i] - tot[np.argmin(abs(deltas))]:.1f}; "
                f"hyperparameters at argmax {fmt(rows[i][5])}; {time.time()-t1:.0f} s")
            for dl, c_, l_, p_, tt, _ in rows:
                say(f"   {dl:+.3f} chi2 {c_:12.1f} logdet {l_:12.1f} hyper {p_:8.2f} total {tt:12.1f}")
            summary[cell]["scans"][label] = {"deltas": [r[0] for r in rows], "total": tot.tolist(),
                                             "chi2": [r[1] for r in rows], "logdet": [r[2] for r in rows],
                                             "argmax": float(deltas[i]), "curv_sd": sd}
    json.dump(summary, open(a.out.replace(".txt", ".json"), "w"))
    say("done")


if __name__ == "__main__":
    main()
