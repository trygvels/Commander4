"""Band-averaged spectra in thermodynamic units.

Maps calibrated on the CMB dipole are in K_CMB. A component with brightness-temperature (RJ) spectrum
f(nu) and unit amplitude at nu_ref is observed through a bandpass r(nu), given as the response to
antenna temperature, as

    m = int r(nu) f(nu) dnu / int r(nu) g(nu) dnu,      g(nu) = x^2 e^x / (e^x - 1)^2, x = h nu / k T_cmb,

where g converts a thermodynamic temperature difference to antenna temperature. The CMB has f = g, so
its column is 1 in every band. A delta bandpass at nu_0 gives m = f(nu_0) / g(nu_0).

This is the WMAP convention (Jarosik et al. 2003; Bennett et al. 2013 App. E). Commander3 writes
the same integral as sum(tau f / g) / sum(tau), which differs at second order in the bandwidth and by
less than 1e-3 for the WMAP K to V bands.
"""
import glob
import os

import numpy as np

H_OVER_K = 0.0479924466       # K / GHz
T_CMB = 2.7255


def g_rj(nu):
    """dT_antenna / dT_thermodynamic at nu (GHz)."""
    x = H_OVER_K * np.asarray(nu, dtype=float) / T_CMB
    return x ** 2 * np.exp(x) / np.expm1(x) ** 2


def mbb_rj(nu, beta, T, nu_ref):
    """Modified blackbody in RJ brightness temperature, unity at nu_ref."""
    nu = np.asarray(nu, dtype=float)
    x, x0 = H_OVER_K * nu / T, H_OVER_K * nu_ref / T
    return (nu / nu_ref) ** (beta + 1.0) * np.expm1(x0) / np.expm1(x)


class Band:
    """A frequency band: bandpass (nu in GHz, response to antenna temperature) and a label."""

    def __init__(self, name, nu, resp):
        self.name = name
        nu, resp = np.asarray(nu, dtype=float), np.asarray(resp, dtype=float)
        if nu.size == 1:
            self.nu, self.w = nu, np.ones(1)
        else:
            # trapezoid weights so that sums are integrals
            dnu = np.gradient(nu)
            self.nu, self.w = nu, resp * dnu
        self._norm = float(np.sum(self.w * g_rj(self.nu)))

    @classmethod
    def delta(cls, name, nu0):
        return cls(name, [nu0], [1.0])

    def average(self, f):
        """Band-averaged thermodynamic response to the RJ spectrum f evaluated on self.nu."""
        return float(np.sum(self.w * f)) / self._norm

    def power_law(self, beta, nu_ref):
        """Response to (nu / nu_ref)^beta for an array of beta; also its derivative in beta."""
        beta = np.atleast_1d(np.asarray(beta, dtype=float))
        lr = np.log(self.nu / nu_ref)
        p = np.exp(beta[:, None] * lr[None, :]) * self.w[None, :]
        return p.sum(1) / self._norm, (p * lr[None, :]).sum(1) / self._norm

    def mbb(self, beta, T, nu_ref):
        return self.average(mbb_rj(self.nu, beta, T, nu_ref))

    def effective_frequency(self, beta):
        """nu_eff such that the delta-band response to a power law of index beta equals the band average."""
        s, _ = self.power_law(beta, 1.0)
        grid = np.linspace(self.nu.min(), self.nu.max(), 4001)
        return float(grid[np.argmin(np.abs(grid ** beta / g_rj(grid) - s[0]))])


WMAP_RADIOMETERS = {"K": ["K11", "K12"], "Ka": ["Ka11", "Ka12"], "Q": ["Q11", "Q12", "Q21", "Q22"],
                    "V": ["V11", "V12", "V21", "V22"], "W": ["W11", "W12", "W21", "W22", "W31", "W32", "W41", "W42"]}
WMAP_NOMINAL = {"K": 22.8, "Ka": 33.0, "Q": 40.7, "V": 60.8, "W": 93.5}


def wmap_band(band, bandpass_dir):
    """WMAP band from the nine-year per-radiometer profiles (`wmap_bandpass_<R>_v5.cbp`: frequency in GHz
    and the response of the two detectors in linear RJ units). Detectors are summed, which weights each by
    the area under its curve, as the LAMBDA description recommends."""
    nus, resp = None, None
    for rad in WMAP_RADIOMETERS[band]:
        fn = os.path.join(bandpass_dir, f"wmap_bandpass_{rad}_v5.cbp")
        a = np.loadtxt(fn)
        if nus is None:
            nus, resp = a[:, 0], np.zeros(a.shape[0])
        r = np.interp(nus, a[:, 0], a[:, 1] / np.trapezoid(a[:, 1], a[:, 0]), left=0.0, right=0.0) \
            + np.interp(nus, a[:, 0], a[:, 2] / np.trapezoid(a[:, 2], a[:, 0]), left=0.0, right=0.0)
        resp += r
    resp = np.clip(resp, 0.0, None)
    return Band(f"WMAP_{band}", nus, resp)
