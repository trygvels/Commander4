"""Removing elements from a band with a dense noise covariance, exactly.

Dropping elements m from the data is marginalizing over them: the kept elements k have covariance N_kk, whose
inverse is the Schur complement of the full inverse covariance,

    (N_kk)^-1 = Ninv_kk - Ninv_km (Ninv_mm)^-1 Ninv_mk .

For a band that has been reduced to its large-scale coefficients (a singular Ninv, `e042e_lowpass.py`) the block
Ninv_mm can be singular and its pseudo-inverse is used: directions of the masked elements that the band does not
constrain carry nothing to remove. The CMB prior on the kept elements is the inverse of the kept block of its
covariance.
"""
import numpy as np


def marginalize(d, ninv, keep, pinv=False):
    """Data and inverse covariance of the kept elements. `ninv` dense (n, n) or diagonal (n,)."""
    keep = np.asarray(keep, dtype=bool)
    if ninv.ndim == 1:
        return d[keep], ninv[keep]
    m = ~keep
    if not m.any():
        return d, ninv
    kk = ninv[np.ix_(keep, keep)]
    km = ninv[np.ix_(keep, m)]
    mm = ninv[np.ix_(m, m)]
    if pinv:
        out = kk - km @ np.linalg.pinv(0.5 * (mm + mm.T), rcond=1e-10, hermitian=True) @ km.T
    else:
        out = kk - km @ np.linalg.solve(mm, km.T)
    return d[keep], 0.5 * (out + out.T)


def restrict_prior(cov, keep):
    """Prior precision of the kept elements from the full prior covariance."""
    keep = np.asarray(keep, dtype=bool)
    sub = cov[np.ix_(keep, keep)]
    return np.linalg.inv(0.5 * (sub + sub.T))
