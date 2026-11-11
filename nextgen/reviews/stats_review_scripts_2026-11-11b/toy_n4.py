"""A small HEALPix toy for the simulations of this review: nside 4 (384 elements per band), five delta bands,
three latitude regions, dense noise with a large-scale excess in the four low bands. Reviewer's own construction."""
import sys

import healpy as hp
import numpy as np

sys.path.insert(0, "nextgen")
sys.path.insert(0, "nextgen/reviews/stats_review_scripts_2026-11-11b")
from lowres.likelihood import BandData  # noqa: E402
from lowres.window import WindowMatcher  # noqa: E402
import rv2_common as rv  # noqa: E402

NSIDE = 4
NPIX = hp.nside2npix(NSIDE)
N = 2 * NPIX
FREQS = rv.FREQS
BETA_TRUE = np.array([-3.10, -3.25, -2.95])


def regions():
    lat = np.abs(hp.pix2ang(NSIDE, np.arange(NPIX), nest=True, lonlat=True)[1])
    return np.tile(np.digitize(lat, [20.0, 45.0])[::1] * 0 + (2 - np.digitize(lat, [20.0, 45.0])), 2)   # 0 high, 1 mid, 2 plane


class Toy:
    def __init__(self, seed=5, sigma=(0.20, 0.20, 0.20, 0.20, 0.20), red=1.0, amp=(3.0, 5.0, 10.0), target=60.0):
        self.rng0 = np.random.default_rng(seed)
        self.WM = WindowMatcher(NSIDE, target_arcmin=target)
        self.region = regions()
        Y, ell = self.WM.Y, self.WM.ell
        self.covs, self.chol, self.ninv = [], [], []
        for i, s in enumerate(sigma):
            if i < 4 and red > 0:
                c = s ** 2 * (np.eye(N) + red * (Y * (ell / 2.0) ** -2.0) @ Y.T * 4.0)     # white plus a red large-scale term
            else:
                c = s ** 2 * np.eye(N)
            c = 0.5 * (c + c.T)
            self.covs.append(c)
            self.chol.append(np.linalg.cholesky(c))
            ni = np.linalg.inv(c)
            self.ninv.append(0.5 * (ni + ni.T))
        self.a_s = self.rng0.standard_normal(N) * np.asarray(amp)[self.region]
        self.a_d = self.rng0.standard_normal(N) * 4.0 * (1.0 + 2.0 * (self.region == 2))
        self.s_true = [rv.sync_fun(f)(BETA_TRUE)[0][self.region] for f in FREQS]
        self.fd = [rv.dust_coef(f) for f in FREQS]

    def signal(self, a_s=None, s_cols=None):
        a_s = self.a_s if a_s is None else a_s
        s_cols = self.s_true if s_cols is None else s_cols
        return [s * a_s + f * self.a_d for s, f in zip(s_cols, self.fd)]

    def noise(self, rng):
        return [L @ rng.standard_normal(N) for L in self.chol]

    def bands(self, dvecs, ninvs=None):
        ninvs = self.ninv if ninvs is None else ninvs
        return [BandData(f"{int(f)}", d, ni, rv.sync_fun(f), fd) for f, d, ni, fd in zip(FREQS, dvecs, ninvs, self.fd)]


def zstats(est, sd, truth):
    z = (est - truth) / sd
    m = z.shape[0]
    return z.mean(0), z.std(0, ddof=1), np.mean(np.abs(z) < 1, axis=0), 1 / np.sqrt(m), 1 / np.sqrt(2 * (m - 1))
