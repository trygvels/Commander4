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
| Cosmoglobe DR1 frequency maps (WMAP K to W, LFI 30/44/70) and component maps (cmb, synch, dust, ame, freefree), `param_CG_v09.txt`, `synch_beta/` | `~/data/cosmoglobe_DR1/` | 8.1 GB | 2026-10-25 | Used in E014. WMAP in mK, LFI in uK. Dust map referenced at 70 GHz (F13). |
| DR1 maps, dipole and high-latitude monopole removed, plus rms maps and the nside-64 masks (dx12 common, b20) | `~/data/q1_inputs/` | 2.3 GB | 2026-10-25 | Made by `nextgen/scripts/prep_q1_inputs.py`; the E014 inputs. |
| BeyondPlanck v2 inputs: BP frequency and component maps, Haslam (`haslam_n0512_60arcmin_uK`, rms 1 percent + 0.8 K), WMAP9 pol maps, AME masks, `BP_mean_v2.h5`, bandpasses | `~/data/beyondplanck/` | 10 GB | 2026-10-25 | Haslam map and rms used in E014. |
| NPIPE A-split 353 and 857 GHz (UiO BeyondPlanck copies, `npipe6v20A_*_map_uK`, rms mcscaled/hrscaled) and their nside-512 preparations | `~/data/beyondplanck/npipe/`, `.../npipe/n512/` | 3.2 GB | 2026-10-25 | Header says K_CMB, values are uK_CMB; dipole kept; 857 in BP's band-integrated K_CMB convention. Fixed by `prep_pr4_maps.py --release npipe_uio` (E014 pass 2). Used in E014. |
| Planck PR3 100 GHz full mission (R3.01, nside 2048) and nside-512 preparation | `~/data/planck_pr3/`, `.../n512/` | 2.0 GB | 2026-10-25 | Used in E014. |
| Planck PR3 143 and 217 GHz (R3.01) | `~/data/planck_pr3/` | 62 MB each, TRUNCATED | 2026-10-25 | PLA dropped the transfers at about 60 MB (three attempts); unusable, not used. Retry from PLA or a mirror. |
| Planck PR4 545 and 857 GHz full mission (R4.00, nside 2048) and nside-512 preparations | `~/data/planck_pr4/`, `.../n512/` | 1.4 GB | 2026-10-25 | Downloaded from NERSC at 30 kB/s (other PR4 bands not fetched). Prepared but not used in pass 2 (replaced by the NPIPE A-split UiO files; the full-mission 545/857 remain available). |
| Haslam 408 MHz destriped, desourced (Remazeilles 2014) | `~/data/haslam/` | 12 MB | 2026-10-25 | Reference copy; E014 uses the BP-smoothed version above. |
| Cosmoglobe DR2 | `~/data/cosmoglobe_DR2/` | empty | | Not downloaded (needs the DR2 parameter file; see STATUS "Needs Trygve"). |
| Commander4 chain outputs (E014 passes) | `~/data/c4_chains/q1_n64/` | 10 GB | 2026-10-25/26 | ridge / resolve_chisq / conditional / marginal, suffixes nohfi, npipemask, npipetuned. |

## External low-frequency data: verified 2026-10-25

Every statement below comes from a page or file listing fetched on 2026-10-25, or from a paper
whose abstract (A) or full text (F) was read. Anything not shown by a fetched page is marked
(inference) or (not verified). Transfer speeds were not measured in this session. No em-dashes.

### Summary table

| Dataset | Public | Host | Format | Cite |
|---|---|---|---|---|
| QUIJOTE MFI wide survey DR1 (11/13/17/19 GHz) | yes | LAMBDA | HEALPix nside 512 RING, Galactic, mK_CMB, IQU + hits + weights + QU cov | 2301.05113 (A) |
| C-BASS North 4.76 GHz | no (not released) | none found | n/a | 2202.10411 (A) for preliminary-map science |
| S-PASS 2.3 GHz DR1 | yes | LAMBDA and INAF site | HEALPix nside 1024 RING, Galactic, K, separate I/Q/U + sensitivity | 1903.09420 (A) |
| CLASS 40 GHz DR1 (2016 to 2022) | yes | LAMBDA | HEALPix nside 128 RING, Celestial, K, V/Q/U + hits blocks | 2305.01045 (A), 2309.00675 (A), 2404.17567 (A) for combined |
| CLASS 90 GHz | spectra only, no maps | LAMBDA | EE C_l tables | 2501.11904 (A) |
| Haslam 408 MHz, Remazeilles 2014 | yes | LAMBDA | HEALPix nside 512 RING, Galactic, K | 1411.3628 (F) |
| WMAP 9-year per-DA I/Q/U | yes | LAMBDA | HEALPix nside 512 NESTED, Galactic, mK thermodynamic, + N_OBS and QQ/QU/UU weights | WMAP DR5 (page) |
| Planck PR4 / NPIPE frequency maps | yes | NERSC portal; CMB-S4 Globus mirror (direct HTTPS) | PLA-format R4.00 FITS, nside 1024 (LFI) / 2048 (HFI) | 2007.04997 (not re-read here) |

### 1. QUIJOTE MFI wide survey, Data Release 1 (public, LAMBDA)

- Release page: `https://lambda.gsfc.nasa.gov/product/quijote/quijote_mfi_data_get.html` (fetched).
  Description page: `https://lambda.gsfc.nasa.gov/product/quijote/quijote_mfi_data_info.html` (fetched).
  The IAC page `research.iac.es/proyecto/quijote/pages/publications/mfi-wide-survey-data-release.php`
  timed out twice on 2026-10-25; the LAMBDA copy is the one verified.
- Files (base `https://lambda.gsfc.nasa.gov`, all delivered 2023-01-01 unless stated):
  - Sky maps `/data/suborbital/QUIJOTE/QUIJOTE-MFI-Release1/quijote_mfi_skymap_{11,13,17,19}ghz_512_dr1.fits`,
    113.3 MB each; 1-degree smoothed `quijote_mfi_smth_skymap_{11,13,17,19}ghz_512_dr1.fits`, 63 MB each.
  - Half-mission null-test splits `.../nulltests/quijote_mfi_[smth_]skymap_XXghz_512_dr1_half{1,2}.fits`
    (16 files, same sizes).
  - Masks `.../masks/mask_quijote_ncp_lowdec_satband_nside512.fits` (default analysis mask) and
    `mask_quijote_satband_ncp86_nside512.fits`, 12.6 MB each.
  - RIMO `.../RIMO/rimo_quijote_mfi_{params,beamrp,beamtf,bandpass}_dr1.fits` (12 KB to 374 KB).
  - Point source catalogues, polarized component-separation maps at nside 64, Galactic-plane
    intensity component maps (v1.1, 2023-12-19), cPILC beta_s and moment maps at nside 32/64 (v1.2,
    2025-12-02).
  - Explanatory Supplement v1.2 (2025-10-22, 20 pages): `/data/suborbital/QUIJOTE/QUIJOTE-MFI-Release1/QUIJOTE_MFI_DR1_Explanatory_Supplement_v1.2.pdf` (read, pages 1 to 12).
- Format (Supplement Sect. 3 and Table 5): HEALPix RING, nside 512, Galactic, Stokes I/Q/U in
  mK_CMB, COSMO polarization convention ("to convert from COSMO to IAU, you have to swap the sign
  of the U Stokes parameter"). Nine columns: I_STOKES, Q_STOKES, U_STOKES, NHITS_I, NHITS_QU,
  WEI_I, WEI_Q, WEI_U (mK_CMB^-2), COV_QU (mK_CMB^2). Hits are 40 ms samples (integration time
  0.04 N_hits s). Weight and covariance maps "are renormalised using the factors described in
  Table 12 of Rubino-Martin et al. 2023 to also account for residual 1/f noise" so that they match
  the null-test one-point statistics. No separate rms map; rms = WEI^-1/2 (inference from the
  column definition).
- Beams (Supplement Table 3, from Rubino-Martin 2023 Table 3): FWHM 55.38 / 55.84 / 38.95 / 40.32
  arcmin at 11 / 13 / 17 / 19 GHz; effective frequencies 10.98 / 12.89 / 16.85 / 18.85 GHz for
  alpha = -1. The 17 and 19 GHz maps are fixed-weight combinations of horns 2 and 4 (weights
  Table 2: w_417 = 0.362 (I), 0.732 (QU); w_419 = 0.419 (I), 0.788 (QU)); the effective beam,
  B_l and bandpass are the same linear combinations (Eqs. 3, 5, 8). Colour-correction polynomials
  in Table 4 and in the `fastcc` code.
- Noise (Supplement Sect. 3): "sensitivities in polarization (Stokes Q and U maps) are within the
  range 35 to 40 uK per 1-degree beam. In intensity, they are within the range 60 to 200 uK per
  1-degree beam" with clear 1/f residuals on large scales in intensity, worst at 17 and 19 GHz;
  11 and 13 GHz share common-mode 1/f noise, as do 17 and 19 GHz.
- Calibration: abstract of 2301.05113 (A) gives "a conservative overall calibration uncertainty
  of 5%".
- Transfer-function caveat (Supplement Sect. 3.1, FDEC filtering): residual RFI is removed with a
  function of declination (median of all pixels at each declination, built with a |b| < 10 deg
  mask and masks of the 10 per cent brightest pixels), applied separately in I and in P (P in
  local equatorial coordinates). "When correlating the final QUIJOTE MFI wide-survey maps with
  external datasets, it should be taken into account that those large scale modes are not present
  in the QUIJOTE maps. For a direct comparison, the FDEC filtering should also be applied to the
  external datasets." Code: `https://github.com/jarubinomartin/sancho` (IDL and Python). The
  LAMBDA info page repeats: "before making comparisons between QUIJOTE and other data sets, please
  be sure to consult Section 3.1 on the filtering applied to the QUIJOTE maps." Implication for
  Commander4 (inference): QUIJOTE cannot be ingested as a plain frequency map in a joint
  pixel-domain fit; either the FDEC operator must be applied to the model sky at every likelihood
  evaluation, or a mask/large-scale marginalisation is needed.
- Masks: default mask excludes dec > 70 deg (NCP), dec < -12 deg, and the satellite band
  -12 < dec < 6 deg; the satband+ncp86 mask is already applied to all released maps.
- Cite: Rubino-Martin et al. 2023, arXiv:2301.05113 (MNRAS 519, 3383); acknowledgement text on the
  LAMBDA info page.

### 2. C-BASS North 4.76 GHz (not public)

- Searches run: web search "C-BASS north 4.76 GHz map public release download LAMBDA"; web search
  with site hints lambda.gsfc.nasa.gov and cbass.web.ox.ac.uk; arXiv API query
  `all:"C-BASS" AND all:northern AND all:survey` sorted by date (15 results, newest 2408.02554,
  none is a data-release paper); LAMBDA foreground index and polarization-surveys page
  (`/product/foreground/fg_pol_survey.html`, fetched: C-BASS row has only the Oxford website link,
  no data link); Oxford site `https://cbass.web.ox.ac.uk/` (fetched: no data access or release page).
- Papers still cite the northern map as "Taylor et al. in prep." (2606.21334 full text, Sect. on
  data, read via arXiv HTML; 2202.10411 abstract says "a preliminary version of the C-BASS
  intensity map").
- Calibration numbers quoted by users of the private map: 2606.21334 states "an absolute
  calibration uncertainty of 2.5%" from the beam model and separately propagates a 5% relative
  calibration uncertainty; search summaries mention ~3% in other papers (not verified).
- Conclusion: no public C-BASS North maps as of 2026-10-25. Access would be on request to the
  collaboration (inference; the Oxford page gives only a contact link). This contradicts
  `nextgen/research/02_cmb_end_to_end_landscape.md` line 165, which says low-frequency anchors
  including C-BASS "are now public".

### 3. S-PASS 2.3 GHz DR1 (public, LAMBDA and INAF)

- LAMBDA description `https://lambda.gsfc.nasa.gov/product/foreground/fg_s_pass_info.html` and
  download `https://lambda.gsfc.nasa.gov/product/foreground/fg_s_pass_get.html` (both fetched).
  INAF page `https://sites.google.com/inaf.it/spass` (fetched): "In this site we publicly release
  the DR1 maps"; DR2 (spectropolarimetric cube) "will be released at a later stage".
- Files (base `https://lambda.gsfc.nasa.gov/data/foregrounds/spass/`, delivered 2019-04-19):
  `spass_dr1_1902_healpix_Tb.{i,q,u}.fits` (49 MB each), `spass_dr1_1902_healpix_Tb_sens.fits`
  (49 MB, Q/U sensitivity in K), Aitoff versions (76 MB each), `rm.fits` and `rm_err.fits`
  (57 KB, 2020-06-24, rotation measure from S-PASS + WMAP + Planck).
- Format (info page): HEALPix RING, nside 1024, Galactic, units K "full-beam Brightness
  Temperature", 1.19 Jy/K, FWHM 8.9 arcmin, effective frequency 2303 MHz, bandwidth 168 MHz
  (2176 to 2216 and 2272 to 2400 MHz). IAU polarization-angle convention, header keyword
  POLCCONV = 'IAU'; "Converting ... means changing the sign of Stokes U" for WMAP/Planck/HEALPix.
- Noise and calibration (info page): Q/U mean sensitivity 0.81 mK per beam; Stokes I confusion
  limited at 9 mK (11 mJy/beam); Q/U zero offset set by the scanning strategy, overall mean signal
  recovered to better than 0.5%; Stokes I zero level set against literature absolute flux
  measurements with error 70 mK (83 mJy/beam). Coverage dec < -1 deg. The LAMBDA page does not
  state an absolute gain-calibration percentage; 2606.21334 assumes 5% for S-PASS.
- Caveats: depolarization strongest in the inner Galaxy (Sagittarius Arm); Faraday rotation at
  2.3 GHz is not discussed on the LAMBDA page but an RM map is provided (use it; inference).
- Cite: Carretti et al. 2019, MNRAS 489, 2330, arXiv:1903.09420; required acknowledgement "This
  work has made use of S-band Polarisation All Sky Survey (S-PASS) data."

### 4. CLASS 40 GHz DR1 and 90 GHz (40 GHz maps public; 90 GHz spectra only)

- Product table `https://lambda.gsfc.nasa.gov/product/class/class_prod_table.html` (fetched) lists
  three DOI releases: 2016 to 2022 40 GHz DR1 products; 2021 to 2024 90 GHz power spectra; 2016
  to 2018 circular polarization.
- 40 GHz maps `https://lambda.gsfc.nasa.gov/product/class/class_dr1_pol_maps_get.html` (fetched):
  `/data/suborbital/class/dr1/maps/class_dr1_40GHz_skymap[_m1|_m2]_n128.fits`, 6.8 MB each,
  2023-08-30. FITS header (fetched): columns V_POLARISATION, Q_POLARISATION, U_POLARISATION in K,
  then VV, VQ, VU, QQ, QU, UU in sec (hits blocks, "proportional to the inverse-variance assuming
  constant white noise"); NSIDE 128, RING, COORDSYS 'C' (Celestial), POLCCONV 'COSMO'. Info page:
  resolution 1.56 deg, coverage about 75% of sky, -76 < dec < +30 deg; no calibration
  uncertainty stated.
- Transfer function `class_dr1_map_xfer_func_get.html` (fetched):
  `/data/suborbital/class/dr1/transfer_function/class_dr1_40GHz_VEB_transfer_matrix_lmax383_{hits,uniform}.hdf5`,
  11 MB each. Info page: a 1152 x 1152 matrix F_ll' with nine 384 x 384 blocks, "Ordering is VV,
  EE, BB on the diagonals", multipole range 0 to 383, acts on C_l ("C_l,out = F_ll' C_l',in"),
  unbiased with a pseudo-C_l estimator for l > 4. The uniform-weight version is "useful for
  simulating CLASS filtering's effect on other experiments". No pixel-space or alm-space operator
  is released. Abstract of 2305.01045 (A): the transfer function recovers about 75% of EE/BB/VV
  power at l = 20 and about 45% at l = 10; white noise 110 uK arcmin; correlated noise rising as
  l^-2.2 with knee l about 16 (linear pol).
- Noise simulations `class_dr1_sims_get.html` (fetched): three tar.gz of 308 MB (coadd, m1, m2);
  count not stated on the page. Beam `class_dr1_beam_get.html`: `class_dr1_40GHz_br.csv` and
  `class_dr1_40GHz_bl.csv` (Moon-derived, with per-l fractional uncertainty); no FWHM stated.
  Masks `class_dr1_masks_get.html`: synchrotron masks s0 to s9 and dust masks d0 to d9 (16 MB
  files) plus a source mask, nside 128.
- Reobserved external maps `class_dr1_reobs_maps_get.html` (fetched): WMAP 9yr K, Ka, Q and Planck
  PR4 353 GHz passed through the CLASS pipeline, nside 128, 2.3 MB each. These are the practical
  route to compare a model sky with CLASS without a pixel-space filter operator (inference).
- Combined CLASS + WMAP Q products `class_dr1_combined_products_get.html` (fetched):
  `class_dr1_combined_map_QU_40GHz_n256.fits` (6.3 MB, nside 256 RING Celestial, 1 deg, 73% sky),
  200 combined noise sims (1.1 GB tar.gz), weights CSV (lmax = mmax = 383), and
  `class_dr1_combined_beta_s_nside32.fits` (7 deg, Galactic). Dust was pre-subtracted with
  beta_d = 1.53, T_d = 19.6 K using PR4 353 GHz. Cite Shi et al. 2024, arXiv:2404.17567.
- 90 GHz `class_90ghz_ee_spectra_2025_doi.html` (fetched): three files under 2 KB, unbinned EE for
  2 <= l <= 30 and binned EE for 31 <= l < 301; DOI 10.71662/9h5m-1369; no maps and no statement
  about future map release. Paper: 2501.11904 (A).
- Cite for 40 GHz: Li et al. 2023, arXiv:2305.01045 (ApJ 956, 77); Eimer et al. 2024,
  arXiv:2309.00675 (ApJ 963, 92).

### 5. Haslam 408 MHz, Remazeilles et al. 2014 version (public, LAMBDA)

- Pages fetched: `https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_info.html` and
  `fg_2014_haslam_408_get.html`. Files (base `https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/`,
  2014-10-15): `haslam408_dsds_Remazeilles2014.fits` (12.6 MB, desourced + destriped, the main
  product), `haslam408_ds_Remazeilles2014.fits` (12.6 MB, destriped only),
  `haslam408_dsds_Remazeilles2014_ns2048.fits` (402.7 MB, small scales added, for simulations
  only), `haslam408_bitmask_Remazeilles2014.fits` (6.3 MB; bit 0 Gaussian-fitted sources, bit 1
  inpainted).
- FITS header (fetched): TEMPERATURE in K, NSIDE 512, ORDERING RING, COORDSYS GALACTIC; header
  cites arXiv:1411.3628.
- Correct citation: Remazeilles, Dickinson, Banday, Bigot-Sazy, Ghosh 2015, MNRAS 451, 4311,
  arXiv:1411.3628 (verified via arXiv API; the LAMBDA info-page text shows the digits transposed
  as 1411.3268 while its link is correct).
- Uncertainty ingredients stated in the paper (F, full text read):
  - Effective beam 56.0 +/- 1.0 arcmin (abstract; Sect. 3.1 gives 56.0 +/- 0.6 from the full map),
    not the 51 arcmin usually quoted. Positional accuracy about 7 arcmin (Sect. 3.2).
  - Zero level: Haslam et al. 1982 quote +/- 3 K; the authors leave the zero level untouched,
    minimum map temperature 11.4 K (Sect. 4.4). Wehus et al. 2014 found a monopole of 8.9 +/- 1.3 K
    including CMB, "somewhat dependent on the foreground model" (Sect. 4.4).
  - Absolute calibration "often quoted as +/- 10%"; the authors argue the overall scale is likely
    better (about 5%) but that the main uncertainty is full-beam versus main-beam scale, which can
    be 30 to 50% for point sources; the map is on the full-beam scale, "the (approximately) correct
    scale for studying diffuse synchrotron emission" (Sect. 4.4).
  - Residual striping: rms reduced from 0.35 K to 0.10 K in a low-background band (Sect. 3.3);
    confusion noise about 0.1 K (Sect. 3.4.4).
  - The paper gives no single recommended per-pixel noise model. A working model consistent
    with it (inference): white rms of order 0.1 K to 1 K per 56 arcmin beam (striping plus
    confusion), a free monopole with prior sigma about 3 K, and a free gain with prior sigma 5 to
    10%. BeyondPlanck/Cosmoglobe inputs `haslam_n0512_60arcmin_uK` with `Haslam_rms_*` files
    (E013, E014) encode one such choice; compare before reuse.

### 6. WMAP 9-year (DR5) per-DA maps (public, LAMBDA)

- Index `https://lambda.gsfc.nasa.gov/product/wmap/dr5/m_products.html` (fetched). Per-DA listing
  pages (fetched): raw I/Q/U `maps_da_r9_iqu_9yr_get.html`; 1 deg smoothed I/Q/U
  `maps_da_smth_r9_iqu_9yr_get.html`. Other listings (links only, not fetched):
  `maps_da_r9_iqus_9yr_get.html` (adds spurious-polarization S map),
  `maps_forered_da_r9_iqu_9yr_get.html`, `maps_deconv_da_r9_i_9yr_get.html` and
  `maps_deconv_da_smth_r9_i_9yr_get.html` (deconvolved, intensity only; no deconvolved Q/U
  products exist), `maps_da_r4_quninv_1yr_get.html` (nside 16 Q/U inverse covariance).
- Files (base `https://lambda.gsfc.nasa.gov/data/map/dr5/skymaps/9yr/`, 2012-11-11):
  raw `raw/wmap_iqumap_r9_9yr_{K1,Ka1,Q1,Q2,V1,V2,W1,W2,W3,W4}_v5.fits`, archive
  `raw/wmap_iqumap_r9_nineyear_v5.tar.gz` (871.16 MB for all);
  smoothed `smoothed/wmap_smth_iqumap_r9_9yr_<DA>_v5.fits`, archive
  `smoothed/wmap_smth_iqumap_r9_nineyear_v5.tar.gz` (872.78 MB). The page lists 96.01 MB per DA.
- FITS header of the smoothed K1 map (fetched): extension 1 columns TEMPERATURE, Q_POLARISATION,
  U_POLARISATION in "mK,thermodynamic" and N_OBS in counts; extension 2 weight arrays N_OBS, QQ,
  QU, UU (counts, "must be divided by a scaling factor sigma_0^2 to provide a weight array in
  physical units (mK^-2)"); NSIDE 512, ORDERING NESTED. Pixel noise sigma = sigma0 / sqrt(N_obs).
- sigma0 per DA (mK, from `skymap_info.html`, fetched): I / Q,U: K1 1.429 / 1.435; Ka1 1.466 /
  1.472; Q1 2.245 / 2.254; Q2 2.131 / 2.140; V1 3.314 / 3.324; V2 2.949 / 2.958; W1 5.899 / 5.912;
  W2 6.562 / 6.577; W3 6.941 / 6.958; W4 6.773 / 6.795. CMB dipole removed from Stokes I.
  Coordinates Galactic. Caveat for the smoothed maps (inference): the N_OBS column after 1 deg
  smoothing does not give the pixel noise of the smoothed map; the sigma0/sqrt(N_obs) formula is
  for the unsmoothed Res 9 maps.
- No separate rms maps are released; N_OBS plus sigma0 is the noise product.

### 7. Planck PR4 (NPIPE) frequency maps: mirrors

- NERSC portal `https://portal.nersc.gov/cfs/cmb/planck2020/` (fetched, last updated 2025-01-07):
  `frequency_maps/{Light_maps,Multi-detector,Single-detector,Single-frequency}/`. Listings fetched:
  - `Single-frequency/`: 31 files, `{LFI,HFI}_SkyMap_<freq>_<nside>_R4.00_{full,full-ringhalf-1,full-ringhalf-2}.fits`
    (LFI 030/044/070 at 1024: 480 MB each, plus `year-1-3` and `year-2-4` for 030 and 044;
    HFI 100 to 353 at 2048: 1.9 GB each; 545/857: 576 MB each). Total about 30 GB.
  - `Light_maps/`: nine `*-field-IQU_*_R4.00_full.fits` (I/Q/U only; 144 MB LFI, 576 MB HFI,
    192 MB at 545/857 intensity only).
  - `Multi-detector/`: 69 files, detector-set maps 2020-06-16 and single-detector maps
    2025-01-07 (217 and 353 single detectors at nside 4096, 2.3 GB).
  The page states no mirror. Full release with 600 simulations and TOD needs a NERSC account.
- CMB-S4 mirror `https://data.cmb-s4.org/planck_pr4.html` (fetched): "The Planck PR4 data is
  publicly available." Categories: full sky maps (27 files, 10.4 GiB), half-ring maps (90 files,
  77.3 GiB), low resolution maps (63 files, 155.4 GiB), effective beams (1049 files, 312.7 MiB),
  single channel maps (171 files, 63.7 GiB). Files use NERSC internal names
  (`npipe6v20_<freq>_map.fits` 144 MiB LFI / 576 MiB HFI, `_hmap.fits`, `_wcov_mcscaled.fits`),
  not the PLA `R4.00` names. Direct HTTPS hrefs on
  `https://g-456d30.0ed28.75bc.data.globus.org/public/planck/planck_pr4/fullsky/<file>` (fetched
  from `planck_pr4-fullsky.html`), plus a Globus web-manager path `/public/planck/planck_pr4/`.
  The half-ring page returned 404 at the guessed URL; navigate from the index. Speed not measured.
- Hugging Face `https://huggingface.co/datasets/astro-legacy-archive/planck-frequency-maps-pr4`
  (fetched): the nine R4.00 full-mission maps as Parquet (8.73 GB total), columns I/Q/U/HITS and
  covariance preserved "exactly as released" per the card, SHA-256 of source FITS in metadata,
  license "unknown". Not a FITS mirror; use only for bulk fetch then verify against NERSC hashes
  (inference).
- Planck Legacy Archive: the PLA wiki (`wiki.cosmos.esa.int/planck-legacy-archive`) reset the
  connection three times; whether PR4 R4.00 maps are served from pla.esac.esa.int was not
  verified here. IRSA `https://irsa.ipac.caltech.edu/data/Planck/` (fetched) has release_1 to
  release_3 only, no PR4.
- Recommendation (inference): for PR4 full-mission I/Q/U, use the CMB-S4 Globus HTTPS links or
  Globus transfer; fall back to NERSC `Light_maps/` for the smaller IQU-only files.

### Corrections to existing notes

- `DATA.md` line 32 to 36: QUIJOTE, S-PASS and CLASS 40 GHz are public on LAMBDA with the URLs
  above; C-BASS North is not public; CLASS 90 GHz has spectra only.
- `research/02_cmb_end_to_end_landscape.md` line 165 ("Low-frequency anchors (C-BASS, ...) are
  now public") is wrong for C-BASS.
- `research/03_instrument_and_external_data.md` line 69: the C-BASS entry should note the map is
  private (Taylor et al. in prep.).

## Unit conventions for HFI 545 and 857 GHz (verified 2026-10-25, E014)

PR4 (NPIPE) delivers 545 and 857 GHz in K_CMB, converted from MJy/sr with the band-integrated HFI
factors (58.04 and 2.27 MJy/sr per K_CMB, Planck 2013 IX). The BeyondPlanck input files at UiO
(`npipe6v20A_857_map_uK.fits`, header says K_CMB, values in uK_CMB) are the PR4 A-split maps in the
same convention (slope 0.99 against the PR4 full map). Commander4 converts a K_CMB band back to
brightness with the monochromatic dB/dT at the band centre (1.44 MJy/sr per K_CMB at 857 GHz), so
feeding these maps as K_CMB makes 857 GHz 1.58x too faint (and 545 GHz about 1.1x). Convert them
back to MJy/sr first (`nextgen/scripts/prep_pr4_maps.py` does this for `--release pr4` and
`npipe_uio`); PR3 545/857 are already in MJy/sr. The same monochromatic conversion costs 2.5 percent
at 353 GHz (colour correction of the delta bandpass), left uncorrected in E014.

PR4 frequency maps keep the solar dipole (PR3 maps do not): subtract 3362.08 uK toward
(l, b) = (264.021, 48.253) before a maps-only run.

## Cosmoglobe DR1 thermal dust reference frequency (verified 2026-10-25, E014)

The DR1 dust component file `CG_dust_IQU_n1024_v1.fits` says `NU_REF_T = 545 GHz` (uK_RJ) in its
header, but its `I_MEAN`, `Q_MEAN` and `U_MEAN` are the dust amplitude at 70 GHz. The run's
parameter file `param_CG_v09.txt` sets `COMP_NU_REF_T03 = 857` and `COMP_NU_REF_P03 = 353` at line
1731 and overrides them at line 1774 with `COMP_NU_REF_T03 = 70 857` and `COMP_NU_REF_P03 = 70`;
Commander3 (`comm_param_mod.f90`, `get_parameter_hashtable(..., par_dp=...)`) reads one double from
the string, so the reference frequency became 70 GHz for I and P. Regressions of the NPIPE and PR4
857 GHz maps on `I_MEAN` at |b| > 30 deg give slope 12.1 to 12.3 against 11.3 predicted for a 70 GHz
reference with the DR1 per-pixel beta and T (1.11 for 545, 1.00 for 857); 353 GHz I gives 7.2 vs
7.3, 353 GHz Q/U give 8.3 vs 7.3. The paper (Watts et al. 2023, Fig. 29 caption) states the dust map
is shown at 70 GHz. The BeyondPlanck dust file `BP_dust_IQU_n1024_v2.fits` carries the same wrong
header (545) against a parameter file saying 857, so the release script writes 545 regardless. The
other DR1 component headers (synchrotron 0.408 / 30, free-free 40, AME 22 GHz) match the parameter
file. Anyone evaluating the DR1 dust SED from the header will be off by about 11x at 857 GHz.
