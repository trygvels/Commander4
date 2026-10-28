"""Q17 (E022): conditional width of the polarized synchrotron index per BP XV region at
BeyondPlanck's own resolution, from one released chain sample's synchrotron amplitude alms.
Information term I = sum over bands (LFI 30, 44, 70 at native nside through their beams; WMAP
Ka, Q, V at nside 16 with the diagonal of N) and Stokes and pixels in the region of
w (a_s f_nu ln(nu/30))^2; sd_cond = 1 / sqrt(I + 1/0.1^2). The same with the release
posterior-mean amplitude map (60 arcmin, no per-pixel noise) shows how much of the width is
set by the sample's small-scale amplitude noise. Usage: bp_sample_cond_width.py [sample]"""
import io, os, sys, time, urllib.request, h5py, numpy as np, healpy as hp
from astropy.io import fits
sys.path.insert(0, os.path.dirname(__file__))
from prep_bpxv_surrogate import cmb_to_rj, col, BP

URL = "https://sdc.uio.no/vol/cosmoglobe-data/BeyondPlanck/releases/v2/BP_c0001_v2.h5"
SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "000500"
LMAX_USE = 1500; NU_S = 30.0
BANDS = (("030", 512, 32.29, 28.4), ("044", 512, 27.0, 44.1), ("070", 1024, 13.25, 70.1))
WMAP = (("WKa", 33.0), ("WQ", 40.6), ("WV", 60.8))
SUR16 = os.path.expanduser("~/data/bpxv_surrogate")
CACHE = os.path.expanduser(f"~/data/beyondplanck/chain_beta_s/c0001_{SAMPLE}_synch_ampQU.npy")


class HttpFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url; self.pos = 0; self.fetched = 0
        self.size = int(urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60).headers["Content-Length"])
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]; return self.pos
    def readinto(self, b):
        if self.pos >= self.size: return 0
        end = min(self.pos + len(b), self.size) - 1
        for attempt in range(5):
            try:
                data = urllib.request.urlopen(urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"}), timeout=300).read(); break
            except Exception:
                if attempt == 4: raise
                time.sleep(10)
        self.fetched += len(data); b[:len(data)] = data; self.pos += len(data); return len(data)


def commander_to_healpy(packed, lmax_in, lmax_out):
    """Commander real-packed alms (index l^2+l+m real, l^2+l-m imaginary, both over sqrt 2 for
    m > 0) to healpy complex alms up to lmax_out."""
    out = np.zeros(hp.Alm.getsize(lmax_out), dtype=np.complex128)
    for l in range(lmax_out + 1):
        base = l * l + l
        out[hp.Alm.getidx(lmax_out, l, 0)] = packed[base]
        if l > 0:
            m = np.arange(1, l + 1)
            out[hp.Alm.getidx(lmax_out, l, m)] = (packed[base + m] + 1j * packed[base - m]) / np.sqrt(2.0)
    return out


t0 = time.time()
if os.path.exists(CACHE):
    qu = np.load(CACHE); print("cached sample alms", qu.shape)
else:
    hf = HttpFile(URL); f = h5py.File(io.BufferedReader(hf, buffer_size=1 << 20), "r")
    g = f[SAMPLE + "/synch"]; lmax = int(g["amp_lmax"][()]); print("amp_lmax", lmax, "shape", g["amp_alm"].shape)
    ncoef = (LMAX_USE + 1) ** 2
    qu = np.array(g["amp_alm"][1:3, :ncoef], dtype=np.float64); np.save(CACHE, qu)
    print(f"read E and B alms up to l {LMAX_USE}: {hf.fetched/1e6:.1f} MB in {time.time()-t0:.0f}s", flush=True)
alm_E = commander_to_healpy(qu[0], LMAX_USE, LMAX_USE); alm_B = commander_to_healpy(qu[1], LMAX_USE, LMAX_USE)
alm_T = np.zeros_like(alm_E)

# regions at nside 16 (surrogate labels), upgraded to each band's nside by nested parent lookup
reg16 = np.load(os.path.join(SUR16, "comps.npz"))["region"]
def upgrade(m_low, nside):
    fac = (nside // 16) ** 2; m_nest = hp.reorder(m_low, r2n=True)
    return hp.reorder(m_nest[np.arange(hp.nside2npix(nside)) // fac], n2r=True)

# release posterior-mean synch map (60 arcmin) for comparison
msyn = f"{BP}/BP_synch_IQU_n1024_v2.fits"; mean_q, mean_u = col(msyn, "Q_MEAN"), col(msyn, "U_MEAN")

NAMES = ("High Latitude", "Spur", "Galactic Center", "Plane")
info_sample = np.zeros(4); info_mean = np.zeros(4); per_band = {}
for name, nside, fwhm, nu in BANDS:
    fn = f"{BP}/BP_{name}_IQU_n{nside:04d}_v2.fits"
    rq, ru = col(fn, "Q_RMS") * cmb_to_rj(nu), col(fn, "U_RMS") * cmb_to_rj(nu)
    bl = hp.gauss_beam(np.radians(fwhm / 60), lmax=LMAX_USE, pol=True)
    pw = hp.pixwin(nside, pol=True, lmax=LMAX_USE)
    aE = hp.almxfl(alm_E, bl[:, 1] * pw[1]); aB = hp.almxfl(alm_B, bl[:, 2] * pw[1])
    _, aq, au = hp.alm2map([alm_T, aE, aB], nside, lmax=LMAX_USE, pol=True)
    # the mean map is already at 60 arcmin; take it as is at this nside (no extra beam)
    mq, mu = hp.ud_grade(mean_q, nside), hp.ud_grade(mean_u, nside)
    reg = upgrade(reg16, nside)
    f_nu = (nu / NU_S) ** -3.2; g = (f_nu * np.log(nu / NU_S)) ** 2
    for r in range(4):
        m = reg == r
        info_sample[r] += g * (np.sum(aq[m] ** 2 / rq[m] ** 2) + np.sum(au[m] ** 2 / ru[m] ** 2))
        info_mean[r] += g * (np.sum(mq[m] ** 2 / rq[m] ** 2) + np.sum(mu[m] ** 2 / ru[m] ** 2))
    per_band[name] = (info_sample.copy(), info_mean.copy())
    print(f"band {name}: sample amp Q rms in Spur {aq[reg==1].std():.2f} uK_RJ (mean map {mq[reg==1].std():.2f}); cumulative sd_cond Spur {1/np.sqrt(info_sample[1]+100):.5f} (mean map {1/np.sqrt(info_mean[1]+100):.5f}) [{time.time()-t0:.0f}s]", flush=True)
# WMAP at nside 16 with the diagonal of N (surrogate band files) and the sample amplitude smoothed to the WMAP beams
for name, nu in WMAP:
    z = np.load(os.path.join(SUR16, f"band_{name}.npz"))
    fwhm = {"WKa": 37.0, "WQ": 29.0, "WV": 20.0}[name]
    bl = hp.gauss_beam(np.radians(fwhm / 60), lmax=LMAX_USE, pol=True); pw = hp.pixwin(16, pol=True, lmax=47)
    aE = hp.almxfl(alm_E, bl[:, 1]); aB = hp.almxfl(alm_B, bl[:, 2])
    _, aq, au = hp.alm2map([alm_T, aE, aB], 16, lmax=LMAX_USE, pol=True)
    mq, mu = hp.ud_grade(mean_q, 16), hp.ud_grade(mean_u, 16)
    f_nu = (nu / NU_S) ** -3.2; g = (f_nu * np.log(nu / NU_S)) ** 2
    for r in range(4):
        m = reg16 == r
        info_sample[r] += g * (np.sum(aq[m] ** 2 / z["sq"][m] ** 2) + np.sum(au[m] ** 2 / z["su"][m] ** 2))
        info_mean[r] += g * (np.sum(mq[m] ** 2 / z["sq"][m] ** 2) + np.sum(mu[m] ** 2 / z["su"][m] ** 2))
print("\nconditional sd of beta_s at BP's resolution (information term, prior N(-3.3, 0.1) included), sample", SAMPLE)
print("| Region | sd_cond with the sample amplitudes | sd_cond with the release mean map | ratio |")
print("|---|---|---|---|")
for r in range(4):
    a, b = 1 / np.sqrt(info_sample[r] + 100), 1 / np.sqrt(info_mean[r] + 100)
    print(f"| {NAMES[r]} | {a:.5f} | {b:.5f} | {b/a:.1f} |")
