"""E022 analysis: conditional, collapsed and Gibbs-chain posteriors of the per-region polarized
synchrotron index on the BP XV low-resolution surrogate (see E022 Target and priors).

Targets per region r (separable; 1-D in beta_r):
  T1 conditional at fixed amplitudes a: -chi2_r(beta)/2 + log prior, evaluated at amplitudes
     drawn from P(a | beta*, d) (beta* = T2 mean with prior for Spur and Plane, the prior mean
     for the prior-only regions), averaged over draws; the same at the BP release posterior-mean
     amplitudes is reported as a consistency diagnostic of the surrogate only
  T2 collapsed, flat amplitude prior: sum_{p,s} [ b^T M^-1 b / 2 - log det M / 2 ] + log prior,
     M = A^T N^-1 A, b = A^T N^-1 d (3 x 3 per pixel and Stokes); profile variant drops log det
  T3 two-block Gibbs: a | beta, d (per pixel 3-D Gaussian) and beta_r | a, d (exact inverse-CDF
     draw of T1 on the grid); High Latitude and Galactic Center from the prior alone as in BP XV.
Also T1 with the dense WMAP QU inverse covariances (other regions fixed at their T1 modes).
"""
import argparse, json, os, time
import numpy as np
import healpy as hp
from astropy.io import fits

T_CMB = 2.7255
H_OVER_K = 0.0479924
NU_S, NU_D = 30.0, 353.0
PRIOR_MEAN, PRIOR_SD = -3.3, 0.1
REGIONS = ("High Latitude", "Spur", "Galactic Center", "Plane")
PRIOR_ONLY = (0, 2)


def cmb_to_rj(nu):
    x = H_OVER_K * nu / T_CMB
    return x**2 * np.exp(x) / np.expm1(x) ** 2


def mbb_rj(nu, beta_d, T):
    """Modified blackbody in RJ units relative to NU_D: (nu/nu0)^(beta+1) (e^{x0}-1)/(e^x-1)."""
    x = H_OVER_K * nu / T; x0 = H_OVER_K * NU_D / T
    return (nu / NU_D) ** (beta_d + 1.0) * np.expm1(x0) / np.expm1(x)


class Surrogate:
    def __init__(self, path, lfi_noise="bp", cmb_prior_sd=0.0):
        z = np.load(os.path.join(path, "surrogate.npz"), allow_pickle=True)
        self.d = z["d"]; self.sigma = z["sigma"].copy(); self.nu = z["freqs"]; self.bands = [str(b) for b in z["bands"]]
        if lfi_noise == "white":
            for ib, b in enumerate(self.bands):
                if b in ("030", "044", "070"): self.sigma[ib] = z["sigma_white"][ib]
        # Gaussian per-pixel prior N(0, sd^2) on the CMB amplitude (uK_CMB); 0 means flat
        self.prior_prec = np.zeros(3); self.cmb_prior_sd = float(cmb_prior_sd)
        if cmb_prior_sd > 0: self.prior_prec[2] = 1.0 / cmb_prior_sd**2
        self.amp = z["amp"]; self.T_d = z["T_d"]; self.beta_d = float(z["beta_d"]); self.region = z["region"]
        self.nb, _, self.npix = self.d.shape
        self.w = 1.0 / self.sigma**2
        ok = np.all(np.isfinite(self.d), axis=(0, 1)) & np.all(np.isfinite(self.sigma), axis=(0, 1)) & (self.sigma.min(axis=(0, 1)) > 0)
        self.ok = ok
        # pixels without valid data get zero data and negligible weight so that every 3 x 3 system
        # stays positive definite; they are excluded from every region sum through `ok`
        self.w[:, :, ~ok] = 1e-12; self.d[:, :, ~ok] = 0.0
        # fixed mixing columns: dust [nb, npix], cmb [nb]
        self.A_d = np.array([mbb_rj(nu, self.beta_d, self.T_d) for nu in self.nu])
        self.A_c = cmb_to_rj(self.nu)
        self.path = path

    def f_s(self, beta):
        return (self.nu / NU_S) ** beta  # [nb]

    def mixing(self, beta_pix):
        """A[nb, npix, 3] for a per-pixel synchrotron index."""
        A = np.empty((self.nb, self.npix, 3))
        A[:, :, 0] = (self.nu[:, None] / NU_S) ** beta_pix[None, :]
        A[:, :, 1] = self.A_d
        A[:, :, 2] = self.A_c[:, None]
        return A

    # --- T1 sufficient statistics: chi2_r(beta) = C0 - 2 sum_nu f_nu S1_nu + sum_nu f_nu^2 S2_nu
    def cond_stats(self, amp, r):
        m = (self.region == r) & self.ok
        a_s = amp[0][:, m]                                   # [2, n]
        res = self.d[:, :, m] - self.A_d[:, None, m] * amp[1][None, :, m] - self.A_c[:, None, None] * amp[2][None, :, m]
        w = self.w[:, :, m]
        S1 = np.einsum("bsp,sp,bsp->b", w, a_s, res); S2 = np.einsum("bsp,sp->b", w, a_s**2)
        C0 = np.einsum("bsp,bsp->", w, res**2)
        return C0, S1, S2

    def cond_logp(self, grid, amp, r, prior=True):
        C0, S1, S2 = self.cond_stats(amp, r)
        f = self.f_s(grid[:, None])                           # [ng, nb]
        chi2 = C0 - 2 * f @ S1 + (f**2) @ S2
        lp = -0.5 * chi2
        if prior: lp = lp - 0.5 * ((grid - PRIOR_MEAN) / PRIOR_SD) ** 2
        return lp

    # --- T2 collapsed
    def collapsed_logp(self, grid, r, prior=True, logdet=True):
        m = (self.region == r) & self.ok
        d = self.d[:, :, m]; w = self.w[:, :, m]; Ad = self.A_d[:, m]; Ac = self.A_c
        out = np.empty(grid.size)
        for i, beta in enumerate(grid):
            fs = self.f_s(beta)
            A = np.empty((self.nb, m.sum(), 3)); A[:, :, 0] = fs[:, None]; A[:, :, 1] = Ad; A[:, :, 2] = Ac[:, None]
            # M[s,p,i,j] = sum_b w[b,s,p] A[b,p,i] A[b,p,j]; b[s,p,i] = sum_b w A d
            M = np.einsum("bsp,bpi,bpj->spij", w, A, A) + np.diag(self.prior_prec); bv = np.einsum("bsp,bpi,bsp->spi", w, A, d)
            L = np.linalg.cholesky(M); y = np.linalg.solve(L, bv[..., None])[..., 0]
            q = 0.5 * np.sum(y**2)
            ld = np.sum(np.log(np.diagonal(L, axis1=-2, axis2=-1))) if logdet else 0.0
            out[i] = q - ld
        if prior: out = out - 0.5 * ((grid - PRIOR_MEAN) / PRIOR_SD) ** 2
        return out

    # --- T3 amplitude draw given per-region beta
    def draw_amp(self, beta_r, rng):
        beta_pix = beta_r[self.region]
        A = self.mixing(beta_pix)                             # [nb, npix, 3]
        M = np.einsum("bsp,bpi,bpj->spij", self.w, A, A) + np.diag(self.prior_prec); bv = np.einsum("bsp,bpi,bsp->spi", self.w, A, self.d)
        L = np.linalg.cholesky(M)
        mean = np.linalg.solve(M, bv[..., None])[..., 0]
        z = rng.standard_normal(mean.shape)
        # sample: mean + L^-T z
        dev = np.linalg.solve(np.swapaxes(L, -1, -2), z[..., None])[..., 0]
        a = mean + dev                                        # [2, npix, 3]
        return np.transpose(a, (2, 0, 1))                     # [3, 2, npix]


def moments(grid, lp):
    p = np.exp(lp - lp.max()); p /= p.sum()
    mean = np.sum(p * grid); sd = np.sqrt(np.sum(p * (grid - mean) ** 2))
    return float(mean), float(sd), p


def draw_from_grid(grid, lp, rng):
    p = np.exp(lp - lp.max()); c = np.cumsum(p) - 0.5 * p; c /= np.sum(p)  # cell-centred CDF
    return float(np.interp(rng.random(), c, grid))


def tau_ips(x):
    """Integrated autocorrelation time by Geyer's initial positive sequence; returns (tau, err)."""
    x = np.asarray(x, dtype=float); n = x.size; xc = x - x.mean()
    f = np.fft.rfft(xc, 2 * n); acf = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); acf /= acf[0]
    tau = -1.0; K = 0
    for k in range(0, n - 1, 2):
        g = acf[k] + acf[k + 1]
        if g <= 0: break
        tau += 2 * g; K = k + 1
    tau = max(tau, 1.0)
    return float(tau), float(tau * np.sqrt(2 * (2 * K + 1) / n))  # Sokal error with the truncation lag


def dense_wmap_cond(sur, amp, beta_mode, grid, r):
    """T1 for region r with the dense WMAP N^-1 on the WMAP bands (others diagonal), other
    regions fixed at beta_mode; returns the log posterior on the grid (with prior)."""
    lp = np.zeros(grid.size)
    C0, S1, S2 = sur.cond_stats(amp, r)
    beta_pix = beta_mode[sur.region].copy()
    m_r = (sur.region == r) & sur.ok
    for ib, b in enumerate(sur.bands):
        f0 = (sur.nu[ib] / NU_S) ** beta_mode[r]
        fgrid = (sur.nu[ib] / NU_S) ** grid
        if b.startswith("W"):
            z = np.load(os.path.join(sur.path, f"band_{b}.npz"), allow_pickle=True)
            with fits.open(str(z["ninv_file"])) as h:
                ninv = np.asarray(h[0].data if h[0].data is not None else h[1].data, dtype=np.float64)
            nest = bool(z["ninv_nested"])
            # map our RING vectors into the N^-1 ordering (mK^-2 -> uK_RJ^-2 handled by scaling)
            scale = 1e3 * cmb_to_rj(sur.nu[ib])                 # map factor mK_CMB -> uK_RJ
            ninv = ninv / scale**2
            def to_n(v):
                return hp.reorder(v, r2n=True) if nest else v
            A0 = sur.mixing(beta_pix)[ib]                        # [npix, 3]
            res0 = np.stack([sur.d[ib, s] - np.einsum("pi,ip->p", A0, amp[:, s, :]) for s in range(2)])
            v = np.stack([np.where(m_r, amp[0, s], 0.0) for s in range(2)])
            res0[:, ~sur.ok] = 0.0; v[:, ~sur.ok] = 0.0
            R = np.concatenate([to_n(res0[0]), to_n(res0[1])]); V = np.concatenate([to_n(v[0]), to_n(v[1])])
            c0 = R @ ninv @ R; c1 = V @ ninv @ R; c2 = V @ ninv @ V
            df = fgrid - f0
            lp += -0.5 * (c0 - 2 * df * c1 + df**2 * c2)
        else:
            lp += -0.5 * (-2 * fgrid * S1[ib] + fgrid**2 * S2[ib])  # this band's diagonal share
    return lp - 0.5 * ((grid - PRIOR_MEAN) / PRIOR_SD) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.expanduser("~/data/bpxv_surrogate"))
    ap.add_argument("--iters", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=2022)
    ap.add_argument("--out", default="nextgen/experiments/E022_bp_xv_polarized_index_gibbs_coupling")
    ap.add_argument("--no-dense", action="store_true")
    ap.add_argument("--lfi-noise", choices=("bp", "white"), default="bp", help="LFI noise sd: BP nside-16 simulated rms, or white-noise propagation")
    ap.add_argument("--cmb-prior-sd", type=float, default=0.0, help="Gaussian per-pixel prior sd on the CMB Q/U amplitude in uK_CMB (0: flat)")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sur = Surrogate(a.dir, lfi_noise=a.lfi_noise, cmb_prior_sd=a.cmb_prior_sd)
    tag = a.tag or f"{a.lfi_noise}_cmb{a.cmb_prior_sd:g}"
    print(f"variant {tag}: LFI noise {a.lfi_noise}, CMB prior sd {a.cmb_prior_sd:g} uK_CMB")
    print(f"bands {sur.bands}, pixels used {sur.ok.sum()} of {sur.npix}; region sizes",
          [int(((sur.region == r) & sur.ok).sum()) for r in range(4)])
    grid = np.linspace(-4.5, -1.5, 2001)
    amp0 = sur.amp
    res = {"bands": sur.bands, "regions": {}, "variant": {"lfi_noise": a.lfi_noise, "cmb_prior_sd": a.cmb_prior_sd, "iters": a.iters}}
    # T2 first (independent of amplitudes), then T1 at amplitude draws consistent with the data
    for r, name in enumerate(REGIONS):
        t0 = time.time()
        m_np = sur.collapsed_logp(grid, r, prior=False); m_p = sur.collapsed_logp(grid, r, prior=True)
        pr_np = sur.collapsed_logp(grid, r, prior=False, logdet=False); pr_p = sur.collapsed_logp(grid, r, prior=True, logdet=False)
        c_bp = sur.cond_logp(grid, amp0, r, prior=True)
        rr = {}
        for key, lp in (("marg", m_np), ("marg_prior", m_p), ("prof", pr_np), ("prof_prior", pr_p), ("cond_bpmean_prior", c_bp)):
            mean, sd, _ = moments(grid, lp); rr[key] = {"mean": mean, "sd": sd, "at_edge": bool(abs(mean - grid[0]) < 0.01 or abs(mean - grid[-1]) < 0.01)}
        rr["npix"] = int(((sur.region == r) & sur.ok).sum())
        res["regions"][name] = rr
        print(f"{name:16s} n={rr['npix']:4d} marg {rr['marg']['mean']:.4f}+-{rr['marg']['sd']:.4f} (prior {rr['marg_prior']['mean']:.4f}+-{rr['marg_prior']['sd']:.4f}) | "
              f"prof {rr['prof']['mean']:.4f}+-{rr['prof']['sd']:.4f} (prior {rr['prof_prior']['mean']:.4f}+-{rr['prof_prior']['sd']:.4f}) | "
              f"cond at BP means {rr['cond_bpmean_prior']['mean']:.4f}+-{rr['cond_bpmean_prior']['sd']:.5f}{' EDGE' if rr['cond_bpmean_prior']['at_edge'] else ''} [{time.time()-t0:.1f}s]")
    beta_star = np.array([PRIOR_MEAN, res["regions"]["Spur"]["marg_prior"]["mean"], PRIOR_MEAN, res["regions"]["Plane"]["marg_prior"]["mean"]])
    rng = np.random.default_rng(a.seed + 1)
    ndraw = 20; sds = np.zeros((ndraw, 4, 2)); means = np.zeros((ndraw, 4, 2)); dense_sd = np.zeros((min(ndraw, 5), 4))
    for k in range(ndraw):
        amp = sur.draw_amp(beta_star, rng)
        for r in range(4):
            for j, pr in enumerate((False, True)):
                means[k, r, j], sds[k, r, j], _ = moments(grid, sur.cond_logp(grid, amp, r, prior=pr))
            if not a.no_dense and k < dense_sd.shape[0] and any(b.startswith("W") for b in sur.bands):
                _, dense_sd[k, r], _ = moments(grid, dense_wmap_cond(sur, amp, beta_star, grid, r))
    print("T1 at amplitude draws from P(a | beta*, d), beta* =", np.round(beta_star, 4))
    for r, name in enumerate(REGIONS):
        rr = res["regions"][name]
        rr["cond"] = {"mean": float(means[:, r, 0].mean()), "sd": float(sds[:, r, 0].mean()), "sd_spread": float(sds[:, r, 0].std()), "mean_spread": float(means[:, r, 0].std())}
        rr["cond_prior"] = {"mean": float(means[:, r, 1].mean()), "sd": float(sds[:, r, 1].mean()), "sd_spread": float(sds[:, r, 1].std()), "mean_spread": float(means[:, r, 1].std())}
        rr["ratio_prior"] = rr["marg_prior"]["sd"] / rr["cond_prior"]["sd"]; rr["ratio_noprior"] = rr["marg"]["sd"] / rr["cond"]["sd"]
        rr["tau_gauss_prior"] = 2 * rr["ratio_prior"] ** 2 - 1
        rr["prior_to_post_width"] = PRIOR_SD / rr["marg_prior"]["sd"]
        if not a.no_dense: rr["cond_prior_denseWMAP"] = {"sd": float(dense_sd[:, r].mean()), "sd_spread": float(dense_sd[:, r].std())}
        print(f"  {name:16s} cond {rr['cond']['mean']:.4f}+-{rr['cond']['sd']:.4f} (sd spread over draws {rr['cond']['sd_spread']:.4f}, mean spread {rr['cond']['mean_spread']:.4f}); "
              f"with prior sd {rr['cond_prior']['sd']:.4f}" + (f"; dense-WMAP sd {rr['cond_prior_denseWMAP']['sd']:.4f}" if not a.no_dense else "") +
              f" | ratio marg/cond {rr['ratio_prior']:.2f} (no prior {rr['ratio_noprior']:.2f}) -> Gaussian tau {rr['tau_gauss_prior']:.1f}; prior/posterior width {rr['prior_to_post_width']:.2f}")

    # T3 Gibbs: BP XV scheme (High Latitude and Galactic Center from the prior), then all four
    # regions from the likelihood (what the flat-prior marginal does to the prior-dominated regions)
    def run_chain(prior_only, ctag, start):
        rng = np.random.default_rng(a.seed)
        beta = np.array(start, dtype=float); chain = np.empty((a.iters, 4)); cvar = np.full((a.iters, 4), np.nan); t0 = time.time()
        for it in range(a.iters):
            amp = sur.draw_amp(beta, rng)
            for r in range(4):
                if r in prior_only: beta[r] = PRIOR_MEAN + PRIOR_SD * rng.standard_normal()
                else:
                    lp = sur.cond_logp(grid, amp, r, prior=True); beta[r] = draw_from_grid(grid, lp, rng)
                    cvar[it, r] = moments(grid, lp)[1] ** 2
            chain[it] = beta
        np.save(os.path.join(a.out, f"gibbs_chain_{ctag}_{tag}.npy"), chain)
        burn = a.iters // 10; out = {}
        for r, name in enumerate(REGIONS):
            if r in prior_only: continue
            x = chain[burn:, r]; tau, err = tau_ips(x); n = x.size
            xc = x - x.mean(); lag1 = float(np.sum(xc[1:] * xc[:-1]) / np.sum(xc**2))
            var_T2 = float(x.var(ddof=1)); e_cvar = float(np.nanmean(cvar[burn:, r]))
            rr = res["regions"][name]
            out[name] = {"mean": float(x.mean()), "mean_mc_err": float(x.std(ddof=1) * np.sqrt(tau / n)), "sd": float(x.std(ddof=1)),
                         "sd_rel_mc_err": float(np.sqrt(tau / (2 * n))), "tau": tau, "tau_err": err, "n": int(n), "ess": float(n / tau),
                         "lag1": lag1, "lag1_identity": 1.0 - e_cvar / var_T2, "mean_cond_sd_along_chain": float(np.sqrt(e_cvar)),
                         "tau_from_lag1": float((1 + lag1) / (1 - lag1)),
                         "mean_minus_T2_in_mc_err": float((x.mean() - rr["marg_prior"]["mean"]) / (x.std(ddof=1) * np.sqrt(tau / n))),
                         "at_edge": bool(np.mean(x > grid[-1] - 0.02) > 0.05 or np.mean(x < grid[0] + 0.02) > 0.05)}
            o = out[name]
            print(f"  [{ctag}] {name:16s} Gibbs mean {o['mean']:.4f}+-{o['mean_mc_err']:.4f} sd {o['sd']:.4f} (rel err {o['sd_rel_mc_err']:.2f}) vs T2 {rr['marg_prior']['mean']:.4f}+-{rr['marg_prior']['sd']:.4f} "
                  f"(mean off by {o['mean_minus_T2_in_mc_err']:.1f} MC err); tau {tau:.1f}+-{err:.1f} (Gaussian {rr['tau_gauss_prior']:.1f}, from lag-1 {o['tau_from_lag1']:.1f}); "
                  f"lag-1 {lag1:.3f} vs identity 1 - E[var_cond]/var {o['lag1_identity']:.3f} (E[sd_cond] along chain {o['mean_cond_sd_along_chain']:.4f}); ESS {n/tau:.0f}{' EDGE' if o['at_edge'] else ''}")
            if a.iters < 200 * rr["tau_gauss_prior"]: print(f"    WARNING: {a.iters} iterations is below 200 tau_gauss = {200*rr['tau_gauss_prior']:.0f}")
        print(f"  [{ctag}] {a.iters} iterations in {time.time()-t0:.0f}s")
        return out
    res["gibbs_bp_scheme"] = run_chain(PRIOR_ONLY, "bp", [PRIOR_MEAN, -3.17, PRIOR_MEAN, -3.0])
    res["gibbs_all_regions"] = run_chain((), "all", [PRIOR_MEAN, -3.17, PRIOR_MEAN, -3.0])
    for name in REGIONS:
        if name in res["gibbs_bp_scheme"]: res["regions"][name]["gibbs"] = res["gibbs_bp_scheme"][name]
    with open(os.path.join(a.out, f"result_{tag}.json"), "w") as f: json.dump(res, f, indent=1)
    print("written", os.path.join(a.out, f"result_{tag}.json"))


if __name__ == "__main__":
    main()
