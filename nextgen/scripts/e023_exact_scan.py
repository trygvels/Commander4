#!/usr/bin/env python
"""E023 stage 1 (v2, after the stats review): exact-target block-shift scans of the collapsed hier
index likelihood on a simulated sky.

Target. As a function of the index fields alone, the per-pixel marginal that Commander4's
`likelihood: hier` samples, with the three amplitudes (CMB, sync, dust) integrated out under their
Gaussian priors and the amplitude hyperparameters theta at their joint MAP (profile), plus the
Laplace term of that profile where requested:

    L(beta) = max_theta { sum_p l_p(beta, theta) + log p(theta) }  [- 1/2 log det H(theta*)],
    l_p = -1/2 r_p^T C_p^-1 r_p - 1/2 log det C_p            (realization mode, r_p = d_p - M_p m_p)
    l_p = -1/2 r_p^T C_p^-1 r_p - 1/2 tr(C_p^-1 N_p) - 1/2 log det C_p
                                                             (expected mode, r_p = sky_p - M_p m_p)
    C_p = M_p T_p M_p^T + N_p,  m_p = (0, mu_s g_s,p, mu_d g_d,p),
    T_p = diag(cmb_sd^2, (s_s,p g_s,p)^2, (s_d,p g_d,p)^2),

over every pixel of the nside-64 sky (no mask), with hyperpriors mu ~ N(1, 0.5^2) and log s ~
N(log 0.3, 1^2) on every mu and log s entry (index_field.HIER_*), log s bounded below by
HIER_LOGS_MIN. The expected mode is the noise-averaged log-likelihood (the trace term is beta- and
theta-dependent and kept), so its peak is the pseudo-true bias without noise scatter; the
realization mode is the single-realization target, whose peak scatters by the block curvature sd.
The index prior is left out: the question is whether the data term is biased. H is the Hessian of
the profiled objective, by central differences of its analytic gradient.

Hyperparameter layouts: "global" (mu_c, log s_c per component), "perq" (log s_c per template-S/N
quartile), "perqmu" (mu_c and log s_c per quartile). The dust SED is one MBB at 545 GHz at the
fixed T or at the PySM T_d per pixel (a reference outside Commander4's model, whose T is scalar).
Templates: the data-built pseudo-templates on disk ("data") or the true amplitudes ("exact").

Reference ("truth" of a region index): the noise-free region fit, one (beta_s, beta_d) per
nside-16 region minimizing the profile chi-square of the noiseless sky with three free amplitudes
per pixel and the rms weights, under the cell's SED; computed per SED level, independent of the
templates and of the noise scale. S/N quartiles of the regions come from the Fisher sd computed
with the data-built template as the amplitude and the prior-mean indices (computable on real
data); 0 is the highest S/N.

Scans: one quartile of regions shifted by delta from the reference field, the other component at
its reference; fine scans over [-0.1, 0.1] for both components and every quartile, a wide sync
scan over [-3, 3] in the faintest quartile. Output: per scan, the argmax, the curvature sd at the
peak, the chi-square, log-determinant, trace and Laplace parts relative to delta 0, the profiled
hyperparameters at delta 0 and at the argmax, and (perq, perqmu) the profiled s_q and mu_q against
the realized sd and mean of a/g per quartile.
"""
import argparse, json, os, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import healpy as hp
from scipy.optimize import minimize

H_K = 4.799243073366221e-11   # h/k in K s, astropy value as in commander4.sky.diffuse_components
T_CMB = 2.7255
FREQS = np.array([23., 30., 33., 41., 44., 61., 70., 94., 100., 143., 217., 353.])
MU_PRIOR, LOGS_PRIOR, LOGS_MIN = (1.0, 0.5), (float(np.log(0.3)), 1.0), -8.0
CMB_SD = 100.0
PRIOR_MEAN = {"sync": -3.10, "dust": 1.60}


def cmb_mixing(nu):
    x = H_K * nu * 1e9 / T_CMB
    return x**2 * np.exp(x) / np.expm1(x)**2          # uK_RJ per uK_CMB


def sync_sed(bs):
    return (FREQS[None, :] / 30.0) ** bs[:, None]


def dust_sed(bd, T):
    x = H_K * FREQS[None, :] * 1e9 / T[:, None]
    xr = H_K * 545e9 / T[:, None]
    return (FREQS[None, :] / 545.0) ** (bd[:, None] + 1.0) * np.expm1(xr) / np.expm1(x)


class Hyper:
    """Layout of theta: mu_s[nmu], log s_s[ns], mu_d[nmu], log s_d[ns]; per-pixel expansion through
    the quartile index of each pixel."""

    def __init__(self, kind, qs, qd):
        self.kind, self.qs, self.qd = kind, qs, qd
        self.nmu = 4 if kind == "perqmu" else 1
        self.ns = 4 if kind in ("perq", "perqmu") else 1
        self.dim = 2 * (self.nmu + self.ns)
        n = qs.size
        self.imu = {"s": qs if self.nmu == 4 else np.zeros(n, int), "d": qd if self.nmu == 4 else np.zeros(n, int)}
        self.ils = {"s": qs if self.ns == 4 else np.zeros(n, int), "d": qd if self.ns == 4 else np.zeros(n, int)}
        o = 0
        self.sl = {}
        for name, k in (("mu_s", self.nmu), ("ls_s", self.ns), ("mu_d", self.nmu), ("ls_d", self.ns)):
            self.sl[name] = slice(o, o + k); o += k

    def theta0(self):
        return np.concatenate([[1.0] * self.nmu, [LOGS_PRIOR[0]] * self.ns, [1.0] * self.nmu, [LOGS_PRIOR[0]] * self.ns])

    def bounds(self):
        return [(None, None)] * self.nmu + [(LOGS_MIN, 3.0)] * self.ns + [(None, None)] * self.nmu + [(LOGS_MIN, 3.0)] * self.ns

    def expand(self, th):
        return (th[self.sl["mu_s"]][self.imu["s"]], th[self.sl["ls_s"]][self.ils["s"]],
                th[self.sl["mu_d"]][self.imu["d"]], th[self.sl["ls_d"]][self.ils["d"]])

    def reduce(self, dmu_s, dls_s, dmu_d, dls_d):
        g = np.empty(self.dim)
        g[self.sl["mu_s"]] = np.bincount(self.imu["s"], dmu_s, self.nmu)
        g[self.sl["ls_s"]] = np.bincount(self.ils["s"], dls_s, self.ns)
        g[self.sl["mu_d"]] = np.bincount(self.imu["d"], dmu_d, self.nmu)
        g[self.sl["ls_d"]] = np.bincount(self.ils["d"], dls_d, self.ns)
        return g

    def hyperprior(self, th):
        lp, g = 0.0, np.zeros(self.dim)
        for name in ("mu_s", "mu_d"):
            x = th[self.sl[name]]; lp += np.sum(-0.5 * ((x - MU_PRIOR[0]) / MU_PRIOR[1])**2)
            g[self.sl[name]] = -(x - MU_PRIOR[0]) / MU_PRIOR[1]**2
        for name in ("ls_s", "ls_d"):
            x = th[self.sl[name]]; lp += np.sum(-0.5 * ((x - LOGS_PRIOR[0]) / LOGS_PRIOR[1])**2)
            g[self.sl[name]] = -(x - LOGS_PRIOR[0]) / LOGS_PRIOR[1]**2
        return lp, g

    def fmt(self, th):
        f = lambda x: f"{x[0]:.4f}" if x.size == 1 else str(np.round(x, 4).tolist())
        return (f"mu_s {f(th[self.sl['mu_s']])} s_s {f(np.exp(th[self.sl['ls_s']]))} "
                f"mu_d {f(th[self.sl['mu_d']])} s_d {f(np.exp(th[self.sl['ls_d']]))}")


class Sky:
    def __init__(self, d, sky, rms):
        self.d, self.sky, self.rms2 = d, sky, rms**2
        self.n = d.shape[0]
        self.cm = np.broadcast_to(cmb_mixing(FREQS), (self.n, 12))
        self.gs = self.gd = None

    def terms(self, Ms, Md, th, hyper, expected, grad=False):
        """Per-pixel sums: (chi2 part, logdet part, trace part) and the gradient wrt theta."""
        mus, ls_s, mud, ls_d = hyper.expand(th)
        M = np.stack([self.cm, Ms, Md], -1)                                   # (n, 12, 3)
        m = np.stack([np.zeros(self.n), mus * self.gs, mud * self.gd], -1)
        Tv = np.stack([np.full(self.n, CMB_SD**2), (np.exp(ls_s) * self.gs)**2,
                       (np.exp(ls_d) * self.gd)**2], -1)
        C = np.einsum("nij,nj,nkj->nik", M, Tv, M) + self.rms2[:, None, :] * np.eye(12)[None]
        r = (self.sky if expected else self.d) - np.einsum("nij,nj->ni", M, m)
        L = np.linalg.cholesky(C)
        ld = -np.log(np.einsum("nii->ni", L)).sum()
        sol = np.linalg.solve(C, r[..., None])[..., 0]
        chi = -0.5 * np.einsum("ni,ni->", r, sol)
        tr = 0.0
        if expected or grad:
            CiM = np.linalg.solve(C, M)                                       # (n, 12, 3)
        if expected:
            Cinv_diag = np.einsum("nij,nij->ni", np.linalg.inv(L).transpose(0, 2, 1), np.linalg.inv(L).transpose(0, 2, 1))
            tr = -0.5 * np.einsum("ni,ni->", self.rms2, Cinv_diag)
        if not grad:
            return chi, ld, tr, None
        u = np.einsum("nij,ni->nj", M, sol)                                    # M^T C^-1 r
        w = np.einsum("nij,nij->nj", M, CiM)                                   # diag M^T C^-1 M
        dmu_s, dmu_d = u[:, 1] * self.gs, u[:, 2] * self.gd
        dls_s = Tv[:, 1] * (u[:, 1]**2 - w[:, 1])
        dls_d = Tv[:, 2] * (u[:, 2]**2 - w[:, 2])
        if expected:
            vNv = np.einsum("ni,nij,nij->nj", self.rms2, CiM, CiM)             # v^T N v per comp
            dls_s = dls_s + Tv[:, 1] * vNv[:, 1]
            dls_d = dls_d + Tv[:, 2] * vNv[:, 2]
        return chi, ld, tr, hyper.reduce(dmu_s, dls_s, dmu_d, dls_d)

    def objective(self, Ms, Md, hyper, expected):
        def f(th):
            chi, ld, tr, g = self.terms(Ms, Md, th, hyper, expected, grad=True)
            lp, gp = hyper.hyperprior(th)
            return -(chi + ld + tr + lp), -(g + gp)
        return f

    def profile(self, Ms, Md, th0, hyper, expected):
        f = self.objective(Ms, Md, hyper, expected)
        res = minimize(f, th0, jac=True, method="L-BFGS-B", bounds=hyper.bounds(),
                       options={"maxiter": 300, "ftol": 1e-13, "gtol": 1e-5})
        chi, ld, tr, _ = self.terms(Ms, Md, res.x, hyper, expected)
        lp, _ = hyper.hyperprior(res.x)
        return res.x, dict(chi=float(chi), ld=float(ld), tr=float(tr), lp=float(lp)), res.nfev

    def laplace(self, Ms, Md, th, hyper, expected, h=1e-3):
        """-1/2 log det H, H the Hessian of the profiled objective at th (central differences of
        the analytic gradient; symmetrized)."""
        f = self.objective(Ms, Md, hyper, expected)
        k = th.size; H = np.empty((k, k))
        for i in range(k):
            e = np.zeros(k); e[i] = h
            H[i] = (f(th + e)[1] - f(th - e)[1]) / (2 * h)
        H = 0.5 * (H + H.T)
        sign, logdet = np.linalg.slogdet(H)
        return float(-0.5 * logdet) if sign > 0 else float("nan")


def region_fit(sky, rms, T, cm, reg, nreg, bs0, bd0, bounds=((-4.0, -2.0), (1.1, 2.1))):
    """Noise-free region fit: per region the (beta_s, beta_d) minimizing the profile chi-square of
    the noiseless sky over three free amplitudes per pixel with rms weights."""
    w = 1.0 / rms**2
    out = np.empty((nreg, 2))
    order = np.argsort(reg, kind="stable"); counts = np.bincount(reg, minlength=nreg)
    starts = np.concatenate([[0], np.cumsum(counts)])
    for r in range(nreg):
        pix = order[starts[r]:starts[r + 1]]
        s_r, w_r, c_r, T_r = sky[pix], w[pix], cm[pix], T[pix]
        def chi2(b):
            if not (bounds[0][0] <= b[0] <= bounds[0][1] and bounds[1][0] <= b[1] <= bounds[1][1]):
                return 1e30
            M = np.stack([c_r, sync_sed(np.full(pix.size, b[0])), dust_sed(np.full(pix.size, b[1]), T_r)], -1)
            A = np.einsum("pbi,pb,pbj->pij", M, w_r, M); rhs = np.einsum("pbi,pb,pb->pi", M, w_r, s_r)
            amp = np.linalg.solve(A, rhs[..., None])[..., 0]
            res = s_r - np.einsum("pbi,pi->pb", M, amp)
            return float(np.sum(res**2 * w_r))
        x0 = np.array([np.clip(bs0[r], *bounds[0]), np.clip(bd0[r], *bounds[1])])
        res = minimize(chi2, x0, method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-12, "maxiter": 2000, "initial_simplex": np.array([x0, x0 + [0.02, 0], x0 + [0, 0.02]])})
        out[r] = res.x
    return out[:, 0], out[:, 1]


def fisher_sd_regions(c, a, beta, rms, T, reg, nreg):
    """Region Fisher sd of the index with the amplitude a known (as e020_analyze.index_fisher_weight),
    1/sqrt(sum over the region's pixels)."""
    nu_ref = 30.0 if c == "sync" else 545.0
    sed = sync_sed(beta) if c == "sync" else dust_sed(beta, T)
    wpix = np.sum((a[:, None] * sed * np.log(FREQS / nu_ref)[None, :] / rms) ** 2, axis=1)
    return 1.0 / np.sqrt(np.maximum(np.bincount(reg, wpix, nreg), 1e-300))


def quartiles(f):
    q = np.empty(f.size, int); order = np.argsort(f); n4 = f.size // 4
    for k in range(4):
        q[order[k * n4:(k + 1) * n4 if k < 3 else None]] = k
    return q


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells", nargs="+", default=["exact_trueT_global", "exact_fixT_global", "data_trueT_global",
                                                    "data_fixT_global", "data_fixT_perq", "data_trueT_perq",
                                                    "data_fixT_perqmu"])
    ap.add_argument("--real-cells", nargs="*", default=["data_fixT_global", "exact_trueT_global"],
                    help="cells also scanned in realization mode (no Laplace term)")
    ap.add_argument("--fine", type=float, nargs=3, default=[-0.10, 0.10, 0.01])
    ap.add_argument("--wide", type=float, nargs=3, default=[-3.0, 3.0, 0.2])
    ap.add_argument("--laplace-at", type=int, default=5, help="Laplace term at this many deltas per scan (0 = none)")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    D = a.dir
    t = np.load(f"{D}/truth.npz")
    nside = int(t["nside"]); npix = hp.nside2npix(nside)
    d = np.array([hp.read_map(f"{D}/band{int(f):03d}_map.fits", field=0) for f in FREQS]).T
    rms = np.array([hp.read_map(f"{D}/band{int(f):03d}_rms.fits", field=0) for f in FREQS]).T
    skymap = np.asarray(t["sky"], float)
    if skymap.shape != d.shape:
        skymap = skymap.T
    g_data = {"sync": np.maximum(hp.read_map(f"{D}/template_sync_30GHz.fits", field=0), 1e-3),
              "dust": np.maximum(hp.read_map(f"{D}/template_dust_545GHz.fits", field=0), 1e-3)}
    a_true = {"sync": np.asarray(t["a_sync"], float), "dust": np.asarray(t["a_dust"], float)}
    g_exact = {k: np.maximum(v, 1e-3 * np.mean(np.abs(v))) for k, v in a_true.items()}
    T_fixed = float(t["T_fixed"]); T_true = np.asarray(t["T_d"], float)
    nreg = hp.nside2npix(16)
    reg = hp.ud_grade(np.arange(nreg, dtype=float), nside).astype(int)
    sky = Sky(d, skymap, rms)
    out = open(a.out, "w")
    def say(s):
        print(s, flush=True); out.write(s + "\n"); out.flush()
    say(f"E023 exact-target scans (v2) on {D}: nside {nside}, noise scale {float(t['noise_scale'])}, "
        f"T_fixed {T_fixed:.3f} K; truth sd beta_s {np.asarray(t['beta_s']).std():.4f}, "
        f"beta_eff3 {np.asarray(t['beta_eff3']).std():.4f}, beta_d {np.asarray(t['beta_d']).std():.4f}")
    # quartiles from the template-based Fisher sd (real-data computable)
    qs16 = quartiles(fisher_sd_regions("sync", g_data["sync"], np.full(npix, PRIOR_MEAN["sync"]), rms, np.full(npix, T_fixed), reg, nreg))
    qd16 = quartiles(fisher_sd_regions("dust", g_data["dust"], np.full(npix, PRIOR_MEAN["dust"]), rms, np.full(npix, T_fixed), reg, nreg))
    qs64, qd64 = qs16[reg], qd16[reg]
    fs_true = fisher_sd_regions("sync", a_true["sync"], np.asarray(t["beta_s"], float), rms, T_true, reg, nreg)
    fd_true = fisher_sd_regions("dust", a_true["dust"], np.asarray(t["beta_d"], float), rms, T_true, reg, nreg)
    for c, q16, ft in (("sync", qs16, fs_true), ("dust", qd16, fd_true)):
        say(f"  {c} quartiles (template Fisher sd): median true-amplitude region Fisher sd per quartile "
            + " ".join(f"{np.median(ft[q16 == k]):.4f}" for k in range(4))
            + f"; block sd (quartile, amplitudes known) " + " ".join(f"{1/np.sqrt(np.sum(1/ft[q16 == k]**2)):.5f}" for k in range(4)))
    hi = np.abs(hp.pix2ang(nside, np.arange(npix), lonlat=True)[1]) > 20
    for name in ("sync", "dust"):
        r = g_data[name][hi] / a_true[name][hi]
        say(f"  data template {name}: g/true at |b|>20 16/50/84 pct {np.percentile(r,16):.3f}/"
            f"{np.percentile(r,50):.3f}/{np.percentile(r,84):.3f}")
    # realized a/g per quartile for the two template levels
    realized = {}
    for tmpl, g in (("data", g_data), ("exact", g_exact)):
        realized[tmpl] = {}
        for c, q64 in (("sync", qs64), ("dust", qd64)):
            r = a_true[c] / g[c]
            realized[tmpl][c] = [(float(r[q64 == k].mean()), float(r[q64 == k].std())) for k in range(4)]
            say(f"  realized a/g ({tmpl} template, {c}) mean/sd per quartile: "
                + " ".join(f"{m:.3f}/{s:.3f}" for m, s in realized[tmpl][c]))
    # references: noise-free region fits per SED level
    bs_mean = hp.ud_grade(np.asarray(t["beta_s"], float), 16)
    refs = {}
    for sed in ("fixT", "trueT"):
        T = np.full(npix, T_fixed) if sed == "fixT" else T_true
        bd_mean = hp.ud_grade(np.asarray(t["beta_eff3"] if sed == "fixT" else t["beta_d"], float), 16)
        t0 = time.time()
        bs_ref, bd_ref = region_fit(skymap, rms, T, np.asarray(sky.cm), reg, nreg, bs_mean, bd_mean)
        refs[sed] = (bs_ref, bd_ref)
        say(f"  reference ({sed}): region fit in {time.time()-t0:.0f} s; fit minus region mean: sync rms "
            f"{np.sqrt(np.mean((bs_ref-bs_mean)**2)):.4f} max {np.max(np.abs(bs_ref-bs_mean)):.4f}, dust rms "
            f"{np.sqrt(np.mean((bd_ref-bd_mean)**2)):.4f} max {np.max(np.abs(bd_ref-bd_mean)):.4f}; "
            f"at the bounds: {int(np.sum((bs_ref <= -3.999) | (bs_ref >= -2.001)))} sync, {int(np.sum((bd_ref <= 1.101) | (bd_ref >= 2.099)))} dust")
    np.savez(a.out.replace(".txt", "_refs.npz"), bs_fixT=refs["fixT"][0], bd_fixT=refs["fixT"][1],
             bs_trueT=refs["trueT"][0], bd_trueT=refs["trueT"][1], qs16=qs16, qd16=qd16)
    summary = {"quartiles": {"sync": qs16.tolist(), "dust": qd16.tolist()}, "realized": realized, "cells": {}}
    fine = np.arange(a.fine[0], a.fine[1] + 1e-9, a.fine[2])
    wide = np.arange(a.wide[0], a.wide[1] + 1e-9, a.wide[2])
    if a.quick:
        a.cells, a.real_cells, fine, wide = a.cells[:1], [], np.array([-0.01, 0.0, 0.01]), np.array([])
    for cell in a.cells:
        tmpl, sed, kind = cell.split("_")
        sky.gs, sky.gd = (g_data if tmpl == "data" else g_exact)["sync"], (g_data if tmpl == "data" else g_exact)["dust"]
        T = np.full(npix, T_fixed) if sed == "fixT" else T_true
        hyper = Hyper(kind, qs64, qd64)
        bs_ref, bd_ref = refs[sed]
        bs0, bd0 = bs_ref[reg], bd_ref[reg]
        for expected in ([True, False] if cell in a.real_cells else [True]):
            mode = "expected" if expected else "realization"
            key = f"{cell}:{mode}"
            t0 = time.time()
            th0, p0, nfev = sky.profile(sync_sed(bs0), dust_sed(bd0, T), hyper.theta0(), hyper, expected)
            tot0 = p0["chi"] + p0["ld"] + p0["tr"] + p0["lp"]
            lap0 = sky.laplace(sync_sed(bs0), dust_sed(bd0, T), th0, hyper, expected) if a.laplace_at else 0.0
            say(f"\n== cell {cell} [{mode}]: profiled at the reference field ({nfev} evals, {time.time()-t0:.1f} s): "
                f"{hyper.fmt(th0)}; chi2 {p0['chi']:.1f}, logdet {p0['ld']:.1f}, trace {p0['tr']:.1f}, "
                f"hyperprior {p0['lp']:.2f}, Laplace {lap0:.2f}")
            if kind != "global":
                for c in ("s", "d"):
                    comp = "sync" if c == "s" else "dust"
                    sq = np.exp(th0[hyper.sl[f"ls_{c}"]]); mq = th0[hyper.sl[f"mu_{c}"]]
                    say(f"   {comp}: profiled s_q {np.round(sq, 3).tolist()} against realized sd(a/g) "
                        f"{[round(x[1], 3) for x in realized[tmpl][comp]]}; mu_q {np.round(mq, 3).tolist()} against mean(a/g) "
                        f"{[round(x[0], 3) for x in realized[tmpl][comp]]}")
            summary["cells"][key] = {"theta0": th0.tolist(), "parts0": p0, "laplace0": lap0, "scans": {}}
            scans = [(f"{c}_Q{q+1}", c, q, fine) for c in ("s", "d") for q in range(4)]
            if wide.size:
                scans.append(("s_Q4_wide", "s", 3, wide))
            for label, comp, q, deltas in scans:
                sel = (qs64 if comp == "s" else qd64) == q
                th = th0.copy(); rows = []
                lap_idx = set(np.linspace(0, len(deltas) - 1, a.laplace_at).round().astype(int).tolist()) if a.laplace_at else set()
                t1 = time.time()
                for k, dl in enumerate(deltas):
                    bs = bs0 + dl * sel if comp == "s" else bs0
                    bd = bd0 + dl * sel if comp == "d" else bd0
                    Ms, Md = sync_sed(bs), dust_sed(bd, T)
                    th, p, nfev = sky.profile(Ms, Md, th, hyper, expected)
                    lap = sky.laplace(Ms, Md, th, hyper, expected) - lap0 if k in lap_idx else float("nan")
                    rows.append(dict(delta=float(dl), chi=p["chi"] - p0["chi"], ld=p["ld"] - p0["ld"], tr=p["tr"] - p0["tr"],
                                     lp=p["lp"] - p0["lp"], total=p["chi"] + p["ld"] + p["tr"] + p["lp"] - tot0, laplace=lap, theta=th.copy()))
                tot = np.array([r["total"] for r in rows]); i = int(np.argmax(tot))
                sd = float("nan")
                if 0 < i < len(rows) - 1:
                    h = deltas[1] - deltas[0]; curv = (tot[i + 1] - 2 * tot[i] + tot[i - 1]) / h**2
                    sd = float(np.sqrt(-1.0 / curv)) if curv < 0 else float("nan")
                laps = np.array([r["laplace"] for r in rows]); lapv = laps[np.isfinite(laps)]
                lap_note = f"Laplace term over the scan min {lapv.min():+.2f} max {lapv.max():+.2f}" if lapv.size else "no Laplace"
                tot_lap = tot.copy()
                if lapv.size:
                    tot_lap[np.isfinite(laps)] += laps[np.isfinite(laps)]
                    i_lap = int(np.argmax(tot_lap[np.isfinite(laps)])); i_lap = np.where(np.isfinite(laps))[0][i_lap]
                    lap_note += f"; argmax with Laplace among those deltas {deltas[i_lap]:+.3f}"
                say(f"-- {cell} [{mode}] {label}: {int(sel.sum())} pixels; argmax delta {deltas[i]:+.3f} (curvature sd "
                    f"{sd:.4f}), gain over delta 0 {tot[i] - tot[np.argmin(abs(deltas))]:.2f}; {lap_note}; "
                    f"hyperparameters at argmax {hyper.fmt(rows[i]['theta'])}; {time.time()-t1:.0f} s")
                for r in rows:
                    say(f"   {r['delta']:+.3f} chi2 {r['chi']:12.2f} logdet {r['ld']:12.2f} trace {r['tr']:10.2f} hyper {r['lp']:8.2f} "
                        f"total {r['total']:12.2f}" + (f" laplace {r['laplace']:+8.2f}" if np.isfinite(r["laplace"]) else ""))
                summary["cells"][key]["scans"][label] = {k: [r[k] for r in rows] for k in ("delta", "chi", "ld", "tr", "lp", "total", "laplace")}
                summary["cells"][key]["scans"][label].update(argmax=float(deltas[i]), curv_sd=sd, theta_argmax=rows[i]["theta"].tolist())
                json.dump(summary, open(a.out.replace(".txt", ".json"), "w"))
    say("done")


if __name__ == "__main__":
    main()
