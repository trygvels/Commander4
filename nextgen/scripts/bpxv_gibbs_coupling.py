"""E022 analysis: conditional, collapsed and Gibbs-chain posteriors of the per-region polarized
synchrotron index on the BP XV low-resolution surrogate (see E022 Target and priors).

Targets per region r (separable; 1-D in beta_r):
  T1 conditional at fixed amplitudes a (BP release means): -chi2_r(beta)/2 + log prior
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
    def __init__(self, path):
        z = np.load(os.path.join(path, "surrogate.npz"), allow_pickle=True)
        self.d = z["d"]; self.sigma = z["sigma"]; self.nu = z["freqs"]; self.bands = [str(b) for b in z["bands"]]
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
            M = np.einsum("bsp,bpi,bpj->spij", w, A, A); bv = np.einsum("bsp,bpi,bsp->spi", w, A, d)
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
        M = np.einsum("bsp,bpi,bpj->spij", self.w, A, A); bv = np.einsum("bsp,bpi,bsp->spi", self.w, A, self.d)
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
    p = np.exp(lp - lp.max()); c = np.cumsum(p); c /= c[-1]
    return float(np.interp(rng.random(), c, grid))


def tau_ips(x):
    """Integrated autocorrelation time by Geyer's initial positive sequence; returns (tau, err)."""
    x = np.asarray(x, dtype=float); n = x.size; xc = x - x.mean()
    f = np.fft.rfft(xc, 2 * n); acf = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); acf /= acf[0]
    tau = -1.0
    for k in range(0, n - 1, 2):
        g = acf[k] + acf[k + 1]
        if g <= 0: break
        tau += 2 * g
    tau = max(tau, 1.0)
    return float(tau), float(tau * np.sqrt(2 * (2 * tau + 1) / n))  # Sokal-type error


def dense_wmap_cond(sur, amp, beta_mode, grid, r):
    """T1 for region r with the dense WMAP N^-1 on the WMAP bands (others diagonal), other
    regions fixed at beta_mode; returns the log posterior on the grid (with prior)."""
    lp = np.zeros(grid.size)
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
            res0 = np.where(np.isfinite(res0), res0, 0.0); v = np.where(np.isfinite(v), v, 0.0)
            R = np.concatenate([to_n(res0[0]), to_n(res0[1])]); V = np.concatenate([to_n(v[0]), to_n(v[1])])
            c0 = R @ ninv @ R; c1 = V @ ninv @ R; c2 = V @ ninv @ V
            df = fgrid - f0
            lp += -0.5 * (c0 - 2 * df * c1 + df**2 * c2)
        else:
            C0, S1, S2 = sur.cond_stats(amp, r)  # all bands diagonal; take this band's share
            lp += -0.5 * (-2 * fgrid * S1[ib] + fgrid**2 * S2[ib])
    return lp - 0.5 * ((grid - PRIOR_MEAN) / PRIOR_SD) ** 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=os.path.expanduser("~/data/bpxv_surrogate"))
    ap.add_argument("--iters", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=2022)
    ap.add_argument("--out", default="nextgen/experiments/E022_bp_xv_polarized_index_gibbs_coupling")
    ap.add_argument("--no-dense", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sur = Surrogate(a.dir)
    print(f"bands {sur.bands}, pixels used {sur.ok.sum()} of {sur.npix}; region sizes",
          [int(((sur.region == r) & sur.ok).sum()) for r in range(4)])
    grid = np.linspace(-4.5, -1.5, 2001)
    amp0 = sur.amp
    res = {"bands": sur.bands, "regions": {}}
    beta_mode = np.zeros(4)
    for r, name in enumerate(REGIONS):
        t0 = time.time()
        c_np = sur.cond_logp(grid, amp0, r, prior=False); c_p = sur.cond_logp(grid, amp0, r, prior=True)
        m_np = sur.collapsed_logp(grid, r, prior=False); m_p = sur.collapsed_logp(grid, r, prior=True)
        pr_np = sur.collapsed_logp(grid, r, prior=False, logdet=False); pr_p = sur.collapsed_logp(grid, r, prior=True, logdet=False)
        rr = {}
        for key, lp in (("cond", c_np), ("cond_prior", c_p), ("marg", m_np), ("marg_prior", m_p), ("prof", pr_np), ("prof_prior", pr_p)):
            mean, sd, _ = moments(grid, lp); rr[key] = {"mean": mean, "sd": sd}
        beta_mode[r] = grid[np.argmax(c_p)]
        rr["ratio_prior"] = rr["marg_prior"]["sd"] / rr["cond_prior"]["sd"]
        rr["ratio_noprior"] = rr["marg"]["sd"] / rr["cond"]["sd"]
        rr["tau_gauss_prior"] = 2 * rr["ratio_prior"] ** 2 - 1
        rr["prior_to_post_width"] = PRIOR_SD / rr["marg_prior"]["sd"]
        rr["npix"] = int(((sur.region == r) & sur.ok).sum())
        res["regions"][name] = rr
        print(f"{name:16s} n={rr['npix']:4d} cond {rr['cond']['mean']:.4f}+-{rr['cond']['sd']:.5f} (prior {rr['cond_prior']['sd']:.5f}) | "
              f"marg {rr['marg']['mean']:.4f}+-{rr['marg']['sd']:.4f} (prior {rr['marg_prior']['mean']:.4f}+-{rr['marg_prior']['sd']:.4f}) | "
              f"prof {rr['prof']['mean']:.4f}+-{rr['prof']['sd']:.4f} | ratio {rr['ratio_prior']:.1f} -> tau_gauss {rr['tau_gauss_prior']:.0f} [{time.time()-t0:.1f}s]")
    np.savez(os.path.join(a.out, "grids.npz"), grid=grid)

    if not a.no_dense and any(b.startswith("W") for b in sur.bands):
        print("T1 with dense WMAP N^-1 (other regions at their conditional modes):")
        for r, name in enumerate(REGIONS):
            lp = dense_wmap_cond(sur, amp0, beta_mode, grid, r)
            mean, sd, _ = moments(grid, lp)
            res["regions"][name]["cond_prior_denseWMAP"] = {"mean": mean, "sd": sd}
            print(f"  {name:16s} {mean:.4f}+-{sd:.5f} (diagonal: {res['regions'][name]['cond_prior']['sd']:.5f})")

    # T3 Gibbs
    rng = np.random.default_rng(a.seed)
    beta = np.array([PRIOR_MEAN, -3.17, PRIOR_MEAN, -3.0])
    chain = np.empty((a.iters, 4)); t0 = time.time()
    for it in range(a.iters):
        amp = sur.draw_amp(beta, rng)
        for r in range(4):
            if r in PRIOR_ONLY:
                beta[r] = PRIOR_MEAN + PRIOR_SD * rng.standard_normal()
            else:
                beta[r] = draw_from_grid(grid, sur.cond_logp(grid, amp, r, prior=True), rng)
        chain[it] = beta
        if it % 500 == 0: print(f"  gibbs it {it} beta {np.round(beta, 3)} [{time.time()-t0:.0f}s]")
    np.save(os.path.join(a.out, "gibbs_chain.npy"), chain)
    burn = a.iters // 10
    for r, name in enumerate(REGIONS):
        if r in PRIOR_ONLY: continue
        x = chain[burn:, r]; tau, err = tau_ips(x)
        rr = res["regions"][name]
        rr["gibbs"] = {"mean": float(x.mean()), "sd": float(x.std(ddof=1)), "tau": tau, "tau_err": err, "n": int(x.size),
                       "ess": float(x.size / tau)}
        print(f"{name:16s} Gibbs: mean {x.mean():.4f} sd {x.std(ddof=1):.4f} (T2 with prior: {rr['marg_prior']['mean']:.4f}+-{rr['marg_prior']['sd']:.4f}); "
              f"tau {tau:.1f} +- {err:.1f} (Gaussian prediction {rr['tau_gauss_prior']:.1f}); ESS {x.size/tau:.0f} of {x.size}")
    with open(os.path.join(a.out, "result.json"), "w") as f: json.dump(res, f, indent=1)
    print("written", os.path.join(a.out, "result.json"))


if __name__ == "__main__":
    main()
