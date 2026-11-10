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


LFI_RIMO = "/Users/trygvels/data/planck_lfi_rimo/LFI_RIMO_R3.31.fits"
HFI_RIMO = "/Users/trygvels/data/planck_hfi_rimo/HFI_RIMO_R3.00.fits"
C_CM_GHZ = 29.9792458           # GHz per cm^-1


def lfi_band(freq, rimo=LFI_RIMO):
    """Planck LFI band-averaged bandpass from the RIMO (`BANDPASS_030` etc.: frequency in GHz and the
    transmission). The LFI convention is the WMAP one: the calibrated signal of a source with RJ
    temperature spectrum f is int tau f dnu / int tau g dnu (Planck 2013 results II and V; Commander3
    `comm_bp_mod`, type `LFI`), so the transmission is the response to antenna temperature. Checked: the
    ratios C(alpha) / C(-1) of the colour corrections computed from these bandpasses reproduce the
    published table of the Planck Explanatory Supplement to 3e-4 at 30, 44 and 70 GHz."""
    from astropy.io import fits
    with fits.open(rimo) as h:
        t = h[f"BANDPASS_{int(freq):03d}"].data
        nu, tau = np.asarray(t["WAVENUMBER"], float), np.asarray(t["TRANSMISSION"], float)
    return Band(f"Planck_{int(freq):03d}", nu, np.clip(tau, 0.0, None))


def hfi_band(freq, rimo=HFI_RIMO, cut=1e-5, psb=False):
    """Planck HFI band-averaged bandpass from the RIMO (`BANDPASS_F353` etc., or `BANDPASS_PSB353` with
    psb=True, the polarization-sensitive bolometers only: wavenumber in cm^-1 and the transmission to
    intensity). The calibrated signal of a source with RJ temperature spectrum f is
    int tau nu^2 f dnu / int tau nu^2 g dnu (Planck 2013 results IX; Commander3 type `HFI_cmb`), so the
    response to antenna temperature is tau nu^2. Samples below `cut` of the peak are dropped, which
    removes the far wings where the measured transmission is noise. Checked against the published unit
    conversions (243.5, 371.1, 482.4, 288.2 MJy/sr per K_CMB at 100 to 353 GHz against 244.1, 371.7,
    483.7, 287.5) and the 353 GHz dust colour correction (0.9052 against 0.9036 for 18 K and index 1.6).
    A modified blackbody (1.53, 19.6 K) seen through the 353 GHz band is 1.137 (all bolometers) or 1.121
    (PSBs) times its delta-band value at 353.0 GHz, so a delta band misstates the dust amplitude by that
    factor."""
    from astropy.io import fits
    with fits.open(rimo) as h:
        t = h[f"BANDPASS_{'PSB' if psb else 'F'}{int(freq):03d}"].data
        nu, tau = np.asarray(t["WAVENUMBER"], float) * C_CM_GHZ, np.asarray(t["TRANSMISSION"], float)
    ok = (tau > cut * tau.max()) & (nu > 0.5 * freq) & (nu < 2.0 * freq)
    return Band(f"Planck_{int(freq):03d}", nu[ok], tau[ok] * nu[ok] ** 2)
