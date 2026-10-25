> Raw report from a research subagent, 2026-10-25. Most sources were read at abstract level.
> Treat every number as a lead: read the paper before citing it. Curated version: `../LITERATURE.md`.

## Commander4 literature review, 2022-2026 (up to 25 Oct 2026)

**Method and gaps:** Abstracts come from the arXiv API. Citation lists come from INSPIRE. Data pages were fetched directly. The web-search budget ran out partway through, so four things are unverified and should be treated as unknown:
- whether SO simulations are public;
- the current status of the CMB-S4 project;
- whether the LiteBIRD simulations can be accessed outside the collaboration;
- the sizes of the NPIPE simulation files.

### 1. Status of the Cosmoglobe programme

**DR1 (2023, A&A 679/686): WMAP and LFI TOD analysed jointly, plus HFI, Haslam and other maps as external constraints.**
- One full WMAP+LFI Gibbs sample costs 812 CPU-hours, 603 of them for WMAP low-level processing [2303.08095].
- Gain and transmission-imbalance modelling cuts WMAP polarization residuals.
- The W-band polarization map, which the WMAP team left out, now looks visually consistent with V-band.
- The K-band vs LFI 30 GHz discrepancy appears resolved.
- The CMB dipole comes out about 11 µK higher than WMAP9 and about 2.5σ above BeyondPlanck, but agrees with PR4.

**Other DR1 papers:**
- Birefringence from the reprocessed maps [2305.02268].
- A full-sky polarized synchrotron model [2310.13740]: amplitude rms 3.2 µK at 30 GHz, B/E = 0.40±0.02, β_s = -3.07±0.07 from T-T plots.
- A paper implementing the Racine et al. (2016) sampler for cosmological parameters in Commander [2306.15511]. It costs about 2000 CPU-hours per independent sample, and the authors say this needs to drop by 10x, for example through multigrid preconditioners or derivative-based samplers.

**DR2 (COBE-DIRBE): papers 2024-2026, data online 4 Jul 2025.**
- DIRBE TOD analysed jointly with FIRAS, HFI, WISE and Gaia [2408.10952].
- CIB monopole [2406.01491] and zodiacal light [2408.11004].
- Starlight [2601.07831] and dust/CII correlations [2601.07818].
- A four-component dust model (hot, cold, nearby, Hα-correlated) fitted to HFI [2601.10640] and to DIRBE [2601.07822].

**Newer Commander methods:**
- N+2 mapmaking [2509.22322].
- Bolometer transfer-function deconvolution inside maximum-likelihood mapmaking [2405.07729].
- Frequency-correlated noise in Commander2 [2609.04510].
- A low-frequency all-sky intensity model [2606.21334]: Commander fit to 35 maps including S-PASS, C-BASS, QUIJOTE and Cosmoglobe WMAP/LFI. It finds synchrotron at 4.76 GHz about 2x higher than the Planck 2015 prediction.
- LiteBIRD feasibility with Commander3 [2507.05324]: about 3000 CPU-hours per Gibbs sample, and a 3-year mission of 238 TB (70 TB compressed). Only 1/f noise was simulated so far.

**Commander4:**
- The Cosmoglobe Moriond 2026 overview [2604.16901] describes Commander4 as a multi-node rewrite with a Python front end and C++ core.
  - It lists ongoing work on AKARI, FIRAS and HFI reanalyses, and on SO, CHIPASS, Fermi-LAT and COMAP.
  - It does not mention a "DR3".
- The GitHub repo shows MPI-split TOD and component-separation ranks, a ducc0 dependency, global CG or per-pixel amplitude solves, and Metropolis-Hastings for spectral parameters. The page says nothing about GPU support.
- Eriksen's ICRANet 2024 talk says an ERC Advanced Grant funds it, with a target of about 100k cores and a first science application of re-analysing raw HFI TOD.
- At the Nordita workshop (Aug-Sep 2026), Lunde said Commander4 is less mature than Commander3 but already runs preliminary SO and LiteBIRD simulations.
- Galloway presented "OpenHFI", an end-to-end Bayesian HFI pipeline, at the same workshop. I found no OpenHFI paper.

### 2. BeyondPlanck results and disagreements

**Main results** [2011.05609, 2205.10104, 2208.14276]:
- τ = 0.066±0.013, with error bars about 30% larger than the official pipelines.
- Solar dipole 3362.7±1.4 µK.
- 44 GHz is fully used for the first time; it had failed null tests in the official LFI analysis.
- New effects found: gain over-smoothing, and non-1/f noise at 30 and 44 GHz.
- Gains differ from Planck DR4 by -0.67%, 0.12% and -0.03% at 30, 44 and 70 GHz [2011.08082].
- Strong CMB vs free-free degeneracy at intermediate scales [2208.14276].
- Synchrotron: β_s = -3.17±0.06 in the Spur and -3.03±0.07 in the Plane. The high-latitude region is prior-dominated and the Galactic centre is contaminated by systematics. Dust β_d = 1.64±0.03 [2011.08503].
- AME: the polarization-fraction limit depends on the β_s prior, from p<0.6% down to a nominal 2.5±1.0% detection [2201.03530].

**I found no dedicated rebuttal paper.** I scanned the INSPIRE citation lists of BeyondPlanck I, X, XII and XV and Cosmoglobe DR1. The disagreements in the literature are numerical and methodological, as follows.

**τ:**
- Other low-ℓ analyses come out lower:
  - NPIPE: 0.051±0.006 [2007.04997].
  - PR4 with LoLLiPoP: about 0.058 [2309.10034, 2504.13254].
  - de Belsunce: 0.0598±0.0059 [2205.13968].
  - CLASS: 0.053 (+0.018/-0.019) [2501.11904].
  - elica (LFI70 + SRoll2 HFI + WMAP): 0.0575 (+0.0048/-0.0058), which its authors present as robustly low [2603.22454].
- Sailer et al. argue the opposite. They get τ = 0.090±0.012 from DESI, high-ℓ CMB and lensing, which is in 3-5σ tension with Planck low-ℓ, and they call Planck low-ℓ polarization a "single-point failure" for modern cosmology [2504.16932].
- Hergt et al. find no robust stand-alone τ-driven discrepancy [2602.06115].

**Dipole:** Cosmoglobe DR1's own value is about 2.5σ above BeyondPlanck [2303.08095].

**Synchrotron spectral index:**
- Weiland et al. find steepening to about -3.4 near the south Galactic pole and -3.25±0.04 in the BICEP/Keck field [2203.11445].
- QUIJOTE XIX finds -3.05±0.01 in the plane and -3.13±0.02 at high latitude [2510.17761].
- CLASS+WMAP strongly rejects a uniform β_s (PTE<0.001) [2404.17567].
- The common thread is that these results ask for finer spatial structure than BeyondPlanck's four regions.

**Haslam anchor:**
- Bayesian evidence points to a flux-scale error of up to 1.6x in the Haslam map, in a small region [2409.06770].
- A Gibbs scheme to recalibrate it has been proposed [2509.11894].
- Hoerning's 2x higher synchrotron at 4.76 GHz [2606.21334] also challenges models anchored on Haslam.

### 3. Competing and complementary methods

**Parametric:**
- FGBuster-style maximum likelihood [0804.2645, 1609.03807].
- BICEP/Keck multicomponent likelihood: r<0.036 [2110.00483]. BK XX adds component-separated maps, with 84% correlation to the baseline [2509.21648].
- FURAX JAX/GPU parametric component separation with K-means clustering of spectral parameters [2604.08463]. It is about 100x faster than FGBuster's optimizer, and the 68% upper limit on r drops by about 30% on LiteBIRD-like simulations.

**Minimally informed:** the method [2308.03412], its pixel-domain implementation MICMAC [2405.18365], and forecasts [2502.15511].

**Blind and semi-blind:**
- SMICA [astro-ph/0211504] and GNILC [1103.1166].
- cMILC and ocMILC [2402.17579].
- MC-NILC for LiteBIRD gain requirements [2411.02080].
- NILC combined with simulation-based inference: unbiased on PySM d12, r = (1.09±0.27)e-2 for an input of 0.01 [2512.16869].

**Reported biases on r:**
- SO SAT pipelines: >2σ bias for complex foregrounds, fixed by marginalizing over residuals [2302.04276, 2407.06891].
- CMB-S4: 1-2σ biases with complex foregrounds; σ(r) of 3-5e-4 [2502.04300].
- Forecasts: SO expanded SATs σ_r = 1.2e-3 [2512.15833]. LiteBIRD after its reformation is a single telescope with 12 bands and a requirement of δr<0.002 [2607.25285].

**Differentiable tools:**
- FURAX [2603.19600], including the POMME half-wave-plate mapmaker [2606.31720].
- jax-healpy (docs page).
- s2fft: differentiable GPU spherical harmonic transforms [2311.14670].
- Flinch: field-level inference, where MicroCanonical Langevin beats HMC by orders of magnitude [2510.26691].

**Bayesian work outside Cosmoglobe:**
- de Belsunce et al. [2205.13968, 2207.04903].
- Gibbs-ILC [2209.07179].
- HMC two-population dust from HFI [2310.04410].
- Improved CMB Gibbs samplers (centred overrelaxation), 10-100x in effective sample size [2111.07664].
- An information-field-theory component-separation talk (Nordita, 2026); no paper found.

### 4. Public test data

**Cosmoglobe DR1** (sdc.uio.no/vol/cosmoglobe-data/cosmoglobe_DR1/):
- 13 IQU frequency maps (WMAP K to W4 and LFI). Nside 512, and 1024 for 70 GHz, inferred from file names.
- Component maps and a CMB map.
- Two main chain files of 841 GB each.
- Half-mission chains and maps.
- β_s maps.
- No TOD is listed.

**BeyondPlanck v2:** four chains of 1.6 TB each, plus maps, masks and likelihood. Same host.

**Cosmoglobe DR2:** DIRBE maps (full and half-mission), beams, bandpasses and the four dust components. Chains are "coming soon". Hosted under tsih3.uio.no.

**NPIPE / PR4:**
- Frequency maps, likelihoods and chains on the public NERSC portal: https://portal.nersc.gov/cfs/cmb/planck2020/
- The full release at NERSC /global/cfs/cdirs/cmb/data/planck2020 needs an account.
- The NPIPE paper states it includes 600 signal+noise+systematics simulations and raw and processed TOD.

**WMAP9 (LAMBDA DR5):** per-DA and per-year maps, uncalibrated and calibrated TOD, beams, far sidelobes, bandpasses, templates and likelihood.

**PySM 3:**
- Dust d0-d12; d9/d10/d11 are GNILC-based and d12 is the 3D six-layer MKD model.
- Synchrotron s1-s7; s5/s6 use logpoltens and s7 adds curvature.
- a1/a2, f1, co1-co3, CIB, SZ and radio sources.
- The 2025 PEGS paper defines three skies spanning a range of complexity [2502.20452]. The low/medium/high mapping to presets (d9s4 / d10s5 / d12s7) is my inference, not from the docs.

**LiteBIRD:** first official release of 500 Nside-512 maps and one year of TOD for about 1/3 of detectors [2507.07122]. litebird_sim is public; whether the data are public is unverified.

### 5. Open problems the field names

- **Foreground complexity:** along-line-of-sight 3D dust variation is a "significant challenge" [2411.11649]. The canonical modified-blackbody and power-law description "will reach its limits" for LiteBIRD [2607.21247, 2607.20080]. There is evidence for multiple dust populations [2310.04410, 2601.10640].
- **Spatially varying parameters:** the trade-off between clustering and model complexity [2604.08463]. β_s variation and curvature are poorly constrained without 10-20 GHz data [2203.11445, 2201.03530].
- **Instrument systematics:**
  - Bandpass knowledge, with a budget of Δr<6.5e-6 per effect [2506.22217].
  - Gain: parametric methods need instrument modelling more than NILC does [2411.02080].
  - Beam chromaticity [2503.01791].
  - Coupling of half-wave-plate non-idealities with detector non-linearity [2609.00841].
  - Frequency-correlated noise [2609.04510].
- **Validation:** independent pipelines must cross-validate any detection [2502.04300], and non-Gaussian residual tests are being proposed [2609.17531].
- **τ:** remains the least well-constrained ΛCDM parameter and is far from its cosmic-variance limit [2504.16932].
- **Compute:** a 10x reduction per independent sample is needed [2306.15511]; LiteBIRD-scale volumes need about 3000 CPU-hours per sample [2507.05324].

### Implications for a next-generation Commander

- **Close the 10x cost gap.** Candidate ingredients: multigrid or multilevel preconditioned CG; derivative-based samplers (HMC or MCLMC) for spectral and cosmological parameters; and overrelaxed or auxiliary-variable Gibbs. Make GPU spherical harmonic transforms swappable (ducc0 vs s2fft). Report the result as ESS per CPU-hour against Commander3 on DR1.
- **Replace fixed regions with adaptive spatial models of spectral parameters.** Options are hierarchical or GP priors, clustering, or moment-expansion components. Compare them with Bayesian evidence and posterior-predictive checks. This directly answers the β_s steepening seen by Weiland, QUIJOTE and CLASS.
- **Sample calibration and scale nuisances jointly.** That includes Haslam and other external-map calibration factors, bandpass shifts and beam chromaticity. Add frequency-correlated noise. Low-frequency anchors (C-BASS, S-PASS, QUIJOTE, CLASS 40 GHz) are now public and already used with Commander.  [CORRECTION 2026-10-25: no public C-BASS maps exist; see DATA.md verified section]
- **Publish an end-to-end τ.** The target is LFI+WMAP plus HFI with full error propagation, set explicitly against the τ≈0.058 vs 0.09 debate. This is the clearest high-impact result.
- **Benchmark against shared references.** Run r-bias tests on d9s4 → d10s5 → d12s7 and the PEGS skies, in the same SO and CMB-S4 configurations used by BBPower, FGBuster and NILC. Use the NPIPE 600 simulations and the LiteBIRD release for frequentist validation, keeping the Bayesian vs frequentist simulation distinction of BeyondPlanck IV [2209.04437].
- **Interoperate rather than compete.** FURAX operators, POMME, jax-healpy and litebird_sim inputs are natural plug-ins, and MICMAC or NILC make good cross-checks.
- **Start validation from public products.** Use the DR1 and BeyondPlanck chains for regression tests, WMAP TOD from LAMBDA, and NPIPE TOD at NERSC.

### References (arXiv abs/<ID>)
- **Cosmoglobe:**
  - DR1 and parameter sampling: 2303.08095, 2305.02268, 2310.13740, 2306.15511
  - DR2: 2408.10952, 2406.01491, 2408.11004, 2601.07831, 2601.07818, 2601.10640, 2601.07822
  - Overview, LiteBIRD feasibility, low-frequency model, open science, WMAP Q-band: 2604.16901, 2507.05324, 2606.21334, 2205.11262, 2202.11979
  - New methods and AME/PAH: 2509.22322, 2405.07729, 2609.04510, 2601.20545
- **BeyondPlanck:** 2011.05609, 2205.10104, 2208.14276, 2208.14293, 2011.08082, 2011.08503, 2201.03530, 2201.03509, 2201.08188, 2209.04437
- **τ:** 2007.04997, 2309.10034, 2504.13254, 2205.13968, 2207.04903, 2501.11904, 2603.22454, 2504.16932, 2602.06115
- **Synchrotron and Haslam:** 2203.11445, 2510.17761, 2404.17567, 2409.06770, 2509.11894
- **Methods:** 0804.2645, 1609.03807, 2110.00483, 2509.21648, 2604.08463, 2603.19600, 2606.31720, 2308.03412, 2405.18365, 2502.15511, astro-ph/0211504, 1103.1166, 2402.17579, 2411.02080, 2512.16869, 2311.14670, 2510.26691, 2209.07179, 2310.04410, 2111.07664, 1512.06619, 0709.1058, astro-ph/0310080, 1701.00274
- **Experiments, foregrounds and systematics:**
  - SO: 2302.04276, 2407.06891, 2512.15833, 2503.01791
  - CMB-S4: 2502.04300
  - LiteBIRD: 2607.25285, 2507.07122, 2507.04918, 2506.22217, 2609.00841, 2607.20080, 2607.21247
  - Foreground modelling and validation: 2411.11649, 2609.17531, 2502.20452
- **Pages:**
  - Cosmoglobe: github.com/Cosmoglobe/Commander4, github.com/Cosmoglobe, cosmoglobe.uio.no/products/{cosmoglobe-dr1,cosmoglobe-dr2,beyondplanck}.html, cosmoglobe.uio.no/news/cosmoglobe-dr2.html
  - Talks: indico.icranet.org/event/8/contributions/1947/, indico.fysik.su.se/event/9149/contributions/ (15460, 15464, 15578), indico.global/event/7910/contributions/130990/, isas.jaxa.jp/home/astro_colloq/20260311_abst.txt
  - Data: portal.nersc.gov/project/cmb/planck2020/, lambda.gsfc.nasa.gov/product/wmap/dr5/m_products.html
  - Tools: pysm3.readthedocs.io/en/latest/models/, jax-healpy.readthedocs.io, github.com/CMBSciPol/furax

**Possible memory note:** Commander4 is ERC-AdG-funded, has a Python front end with a C++ core and MPI, and has OpenHFI as its first target. The τ debate, with low-ℓ around 0.058 versus around 0.09 from DESI+CMB, is the clearest opening for a high-impact paper.
