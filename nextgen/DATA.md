# Public data for testing and science

Sources found in the 2026-10-25 literature survey. Sizes and contents are from project pages or
papers as reported by the research agents; anything marked (verify) was not confirmed. Record what
you actually download (path, size, checksum, date) in the "Local copies" table.

## Cosmoglobe and BeyondPlanck products
- **Cosmoglobe DR1** (WMAP + LFI jointly from TOD): `https://sdc.uio.no/vol/cosmoglobe-data/cosmoglobe_DR1/`
  - 13 IQU frequency maps (WMAP K to W4, LFI 30/44/70), nside 512 (1024 at 70 GHz, inferred from file
    names), component maps, CMB map, synchrotron beta maps, half-mission maps.
  - Two main chain files of ~841 GB each. Do not download in full; read selected datasets remotely
    or ask for subsets. No TOD listed.
  - Product page: `https://cosmoglobe.uio.no/products/cosmoglobe-dr1.html`
- **BeyondPlanck v2**: maps, masks, likelihood; four chains of ~1.6 TB each. Same host.
  Page: `https://cosmoglobe.uio.no/products/beyondplanck.html`
- **Cosmoglobe DR2** (DIRBE): full and half-mission maps, beams, bandpasses, four dust components;
  chains "coming soon". Page: `https://cosmoglobe.uio.no/products/cosmoglobe-dr2.html` (host tsih3.uio.no).

## Planck
- **PR4 / NPIPE**: frequency maps, likelihoods, chains on the public portal
  `https://portal.nersc.gov/cfs/cmb/planck2020/` (also referenced as
  `https://portal.nersc.gov/project/cmb/planck2020/`). The full release with 600 simulations and TOD
  is at NERSC `/global/cfs/cdirs/cmb/data/planck2020` and needs an account. Simulation sizes (verify).
- **PR3 (2018) maps and component products**: Planck Legacy Archive, `https://pla.esac.esa.int/` (verify).

## WMAP
- **Final (9-year, DR5)** on LAMBDA: per-DA and per-year maps, calibrated and uncalibrated TOD, beams,
  far sidelobes, bandpasses, templates, likelihood:
  `https://lambda.gsfc.nasa.gov/product/wmap/dr5/m_products.html`

## Low-frequency and ground-based (degeneracy breakers; all verify)
- **QUIJOTE MFI wide survey** (11, 13, 17, 19 GHz, ~1 deg, ~29,000 deg^2, 5% calibration). Papers
  2301.05113, 2301.05117; data release location to verify.
- **C-BASS North** (4.76 GHz, 0.73 deg). Paper 2202.10411; public maps to verify.
- **S-PASS** (2.3 GHz, southern sky; Faraday rotation matters). Paper 1802.01145; maps to verify.
- **CLASS** 40 and 90 GHz large-scale polarization (2404.17567, 2501.11904); maps to verify.
- **Haslam 408 MHz**: widely used; note evidence of local flux-scale errors up to 1.6x (2409.06770).

## Simulations
- **PySM 3** sky models (`pip install pysm3`; templates download on demand): dust d0 to d12 (d9 to d11
  GNILC-based, d12 the 3D MKD model), synchrotron s1 to s7 (s5/s6 log-poltens, s7 with curvature),
  AME a1/a2, free-free f1, CO co1 to co3, CIB, SZ, radio sources. Docs:
  `https://pysm3.readthedocs.io/en/latest/models/`
- **PEGS skies** (three complexity levels, 2502.20452). Mapping to presets such as d9s4 / d10s5 /
  d12s7 is an inference (verify).
- **LiteBIRD**: first official simulation release, 500 nside-512 map realizations and one year of TOD
  for ~1/3 of detectors (2507.07122). Public access outside the collaboration (verify).
  `litebird_sim` is public on GitHub.
- **Commander4's own `simgen`** (in this repo) for TOD-level simulations with known truth.

## Local copies

| Dataset | Path | Size | Date | Notes |
|---|---|---|---|---|
| | | | | |
