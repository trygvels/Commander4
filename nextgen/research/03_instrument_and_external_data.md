> Raw report from a research subagent, 2026-10-25. Most sources were read at abstract level.
> Treat every number as a lead: read the paper before citing it. Curated version: `../LITERATURE.md`.

# Commander4 research brief: instrument degeneracies, local foregrounds, external data, τ, error propagation

Everything below comes from arXiv abstract or full-text pages and official pages that I fetched in this session. Items marked "(inference)" are my own synthesis.

## 1. Instrument-model degeneracies

**Gain vs polarized sky**
- **LFI (BeyondPlanck VII):** gain is split into three terms: an absolute term (orbital dipole only), a per-detector relative term, and a per-pointing-period term. The paper identifies a degeneracy between the dipole and the polarization quadrupole, plus anomalous gain jumps. Adding WMAP helps break them [2011.08082].
- **HFI (NPIPE):** the paper states that the Planck scan strategy "allows for apparent gain fluctuations that are compensated with scan-synchronous polarization residuals".
  - NPIPE breaks this with foreground polarization priors (30/217/353 GHz templates) during calibration.
  - The cost is a large-scale polarization transfer function, which the PR4 likelihoods correct [2007.04997, 2309.10034].
- **WMAP:** the polarized sky signal is strongly degenerate with the transmission-imbalance factors ("poorly measured modes"), worst in V and W.
  - Cosmoglobe DR1 fits the gains to the full sky signal jointly with LFI, which measures those modes well [2303.08095].
  - Results: lower large-scale polarization residuals, W-band visually consistent with V-band, and the K-band vs 30 GHz discrepancy resolved.
  - Still open: 1–2 µK large-scale temperature residuals, noisy W-band gains and baselines, and an excess χ² from the noise approximation. The authors call the maps "preliminary" and report no τ.
- **Dipole tension:** DR1 gives 3366.2±1.4 µK and BeyondPlanck gives 3362.7±1.4 µK, a 2.5σ difference. DR1 agrees with NPIPE (3366.6±2.7 µK).

**Correlated noise.** BeyondPlanck VI finds significant time variation in all noise PSD parameters, tracking housekeeping data. It also finds non-1/f excess at 30/44 GHz within 0.1–1 Hz, modelled with an added lognormal term [2011.06650]. BeyondPlanck I additionally reports "gain over-smoothing" [2011.05609].

**Bandpass**
- **BeyondPlanck X:** the nominal LFI bandpasses were cleaned (61 GHz ground-measurement systematic, standing-wave ripples, band edges), giving a net +0.6 GHz shift at 70 GHz [2201.03417].
  - Fitted 30 GHz shift: Δ₃₀ = 0.24±0.03 GHz.
  - Relative per-detector shifts come from spurious maps.
  - This is the first time beam-mismatch leakage was included for LFI.
- BeyondPlanck I notes that without priors, bandpass parameters and foreground spectral parameters have no stable joint solution.
- HFI 2018 and NPIPE fit bandpass-mismatch coefficients in temperature. NPIPE does not fit them in polarization; it applies a sky-model correction instead [1807.06207, 2007.04997].

**Sidelobes**
- BeyondPlanck VIII uses spin-harmonic convolution, 3–10× faster than LevelS. Its sidelobe rms maps reflect sky-model variations only, not beam uncertainty [2201.03478].
- NPIPE: the far-sidelobe template is "too faint and degenerate" to fit below 545 GHz.
- The official WMAP pipeline neglected far sidelobes; Cosmoglobe models them explicitly.
- LiteBIRD: the integrated far-sidelobe power must be known to about 10⁻⁴ to keep δr < 1.9×10⁻⁵ [2312.09001].

**ADC nonlinearity.** NPIPE: "residual ADC nonlinearity remains the primary systematic error in the large-scale HFI polarization data." SRoll2 fits a higher-order ADC model plus polarization efficiencies; in simulations the residuals fall "well below detector noise" [1901.11386].

**Transfer functions**
- NPIPE fits residual bolometer transfer-function templates.
- PR4 corrects an E-mode transfer function caused by degeneracies between signal and systematics templates.
- CLASS: time-domain filtering removes almost all power at the lowest ℓ and is corrected with a pixel-space transfer matrix. This is the main limit on CLASS τ [2501.11904].

## 2. Zodiacal light and local foregrounds

- **Planck 2013 XIV** [1303.5074]:
  - Fitted emissivities for the COBE (Kelsall 1998, "K98") model: diffuse cloud, three dust bands, circumsolar ring, trailing feature.
  - The cloud's emissivity falls with wavelength; the bands fall more slowly, which suggests different grains.
  - It also fitted far-sidelobe Galactic pickup, and found the zodi correction small compared with the CMB TT spectrum.
- **NPIPE:**
  - Fits six K98 templates per frequency, but only the time-varying part (for example, cloud emissivity 0.53 at 857 GHz and 0.0124 at 100 GHz).
  - The stationary part is deferred to component separation, so it stays degenerate with dust and CIB at map level (inference).
- **Cosmoglobe DR2:**
  - DIRBE time-ordered data fitted jointly with FIRAS, Gaia, HFI and WISE [2408.10952].
  - Zodi parameters re-derived (K98 form) jointly with dust and starlight; predictions within 3% of K98 at 12/25 µm.
  - Found a static solar-centric excess at 4.9–60 µm, and zero-levels about 30 kJy/sr lower than K98 at 1.25–3.5 µm [2408.11004].
  - CIB monopole at 2.2 µm: 10.2±1.2 nW m⁻² sr⁻¹, 74% below the official value. At 12 µm only an upper limit is possible (55 vs the legacy 468) [2406.01491].
- **Gjerløw, FIR2025 talk:**
  - The legacy DIRBE 100 µm map, widely used for dust modelling, has significant zodi contamination.
  - A three-component dust model (cold, hot, nearby) fits 100 GHz to 1 µm.
  - Intermediate HFI channels are not yet in the global fit.
- **Tools and plans:** ZodiPy [2205.12962]; a planned AKARI reanalysis focused on zodi [2604.16901].

## 3. External degeneracy breakers

| Data | Bands | Sky | Key results |
|---|---|---|---|
| QUIJOTE MFI | 11/13/17/19 GHz, 1° | ~29,000 deg², north | Polarization noise 35–40 µK per 1° beam; 5% calibration uncertainty; β(11–23 GHz) peaks at −3.09 (σ 0.14) [2301.05113]. With WMAP+Planck: significant spatial variation in β_s and a hint of curvature [2203.04861]. MFI2 forecast: ~10× tighter β_s, ~5× tighter curvature [2511.14572] |
| C-BASS | 4.76 GHz, 0.73° | North, δ > −15° | β ≈ −3.10 (4.76–22.8 GHz, intensity) [2202.10411]. 2026 Commander intensity model with S-PASS, QUIJOTE and Cosmoglobe WMAP/LFI: synchrotron at 4.76 GHz is ~2× the Planck 2015 prediction [2606.21334]. C-BASS+S-PASS polarization: the North Polar Spur/Loop I and Fan region must be modelled explicitly [2605.23489] |
| S-PASS | 2.3 GHz | Southern sky | β_s = −3.22±0.08 (2.3–33 GHz); B/E ≈ 0.5; dust–synchrotron correlation up to ~40%; no region has total polarized foreground below r ≈ 10⁻³ [1802.01145] |
| CLASS | 40, 90 GHz | Dec −74° to +28°; f_sky 37% at 90 GHz | 40 GHz + WMAP Q-band: 2–3× lower noise at ℓ 30–100; β ≈ −3.08, uniform index rejected (PTE < 0.001) [2404.17567]. 90 GHz × Planck: τ = 0.053 +0.018/−0.019 (ℓ 2–30), the first ground-based τ; CLASS × WMAP gives only τ < 0.070 [2501.11904] |
| SPIDER | 95/150 GHz (280 GHz on 2nd flight) | 4.8% | r < 0.11 (Feldman–Cousins) [2103.13334] |
| BICEP/Keck BK18 | 95/150/220 GHz | Small patch | r < 0.036, σ(r) = 0.009; 220 GHz dust S/N beats Planck 353 GHz in their field [2110.00483]. End-2025 combination: r < 0.034 [2512.10613] |
| LiteBIRD | 15 bands, 34–448 GHz | Full sky | 0.5° at 100 GHz, 2.2 µK·arcmin [2202.02773]. σ(τ) about 3× better than Planck, with bias if the reionization model is wrong (Giardiello, Nordita 2026 talk). Bandpass needs ≲1.5 GHz resolution [2506.22217] |
| Taurus (planned balloon) | 150/220/280/350 GHz | ~70% | Aims for τ within 2× the cosmic-variance limit [2407.01438] |

SPIDER and BICEP/Keck constrain the dust spectrum and r, but not τ directly (inference).

## 4. Status of τ

| Analysis | τ |
|---|---|
| Planck 2018 [1807.06209] | 0.054±0.007 |
| SRoll2 HFI, low-ℓ polarization only / combined [1908.09856] | 0.0566 +0.0053/−0.0062 / 0.059±0.006 |
| SRoll2 TTTEEE (momento) + CamSpec [2103.14378] | 0.0627 +0.0050/−0.0058 |
| NPIPE [2007.04997] | 0.051±0.006 |
| PR4 HiLLiPoP+LoLLiPoP [2309.10034] | 0.058±0.006 |
| PR4, model-averaged [2504.13254] | 0.0576±0.0060 |
| WMAP+LFI legacy maps [2005.05600] | 0.069 +0.011/−0.012 |
| BeyondPlanck, LFI+WMAP [2205.10104] | 0.066±0.013 (errors ~30% larger than official) |
| CLASS × Planck [2501.11904] | 0.053 +0.018/−0.019 |
| elica: LFI 70 + SRoll2 100/143 + WMAP, cross-spectra (2026) [2603.22454] | 0.0575 +0.0048/−0.0058 |

**Tension.** DESI BAO plus high-ℓ CMB and lensing, in ΛCDM, gives τ = 0.090±0.012. That is 3–5σ from the low-ℓ polarization result, which the authors call a "single-point failure" risk for most CMB analyses [2504.16932]. Two independent astrophysical estimates disagree with a high τ:
- 0.0552 (+0.0019/−0.0026 stat, +0.0075/−0.0049 sys) [2601.09644];
- 0.067±0.011 [2607.26373].

**What limits it:**
- HFI: residual ADC nonlinearity, and the gain/polarization degeneracy with its transfer function.
- WMAP: transmission imbalance.
- LFI: noise.
- Ground-based: filtering.
- Also: dependence on the assumed reionization-history prior [2504.13254].

## 5. Bayesian end-to-end vs Monte Carlo error propagation

**Bayesian end-to-end (Commander).**
- *For:*
  - A joint posterior marginalizes over non-linear instrument–sky degeneracies, giving larger and more complete errors [2205.10104].
  - Joint datasets break internal degeneracies.
  - Constrained realizations keep correlations specific to the actual data [2209.04437].
- *Against:*
  - Only effects that are explicitly parameterized get marginalized.
  - Results depend on priors (bandpass vs spectral parameters).
  - Slow convergence: about 2000 samples for the 225-mode low-resolution covariance, and about 4000 for τ.
  - Exact Cℓ/parameter sampling costs about 2000 CPU-hours per independent sample and needs roughly a 10× speedup [2306.15511].
  - LiteBIRD-scale runs cost about 3000 CPU-hours per sample, even with only 1/f noise simulated [2507.05324].

**Simulation-based (NPIPE's 600 simulations, SRoll2 simulations, SimBaL/momento/DELFI, elica's 500 simulations, simulation-based inference).**
- *For:*
  - Known inputs make them good for validation, and cosmic variance is included.
  - Cross-spectra remove noise bias [2603.22454].
  - Emulators can be trained from about 10 simulations (Campeti, 2026 talk), and NILC plus simulation-based inference is robust to model misspecification [2512.16869].
- *Against:*
  - Only fiducial systematics are included, and data-conditioned correlations are missing.
  - The transfer-function correction depends on the fiducial model (inference).

BeyondPlanck IV concludes the two approaches are complementary [2209.04437].

## Implications for Commander4 instrument and data modelling

- **HFI time-ordered data comes first for τ.** Build parametric ADC-nonlinearity, polarization-efficiency, bandpass-mismatch (temperature and polarization) and bolometer transfer-function templates into the Gibbs chain. HFI reanalysis is listed as a future Cosmoglobe project [2604.16901].
- **Make the gain/polarization degeneracy explicit.** Report the effective E-mode transfer function and how sensitive it is to the sky prior. Validate with random-input (frequentist) simulations.
- **Anchor absolute calibration across experiments.** Use the orbital dipole for LFI+WMAP+HFI, and resolve the 2.5σ dipole tension between BeyondPlanck and DR1.
- **Noise model.** Use per-scan, non-stationary PSDs with extra components (such as the lognormal term) and housekeeping-informed priors.
- **Bandpass.** Sample shifts and shapes with physically informed priors, and track their degeneracy with β_s, β_d and T_d.
- **Sidelobes.** Sample sidelobe amplitude or integrated power with GRASP-based priors. The 10⁻⁴ LiteBIRD requirement sets the target precision.
- **Zodi in the time domain.** Use ZodiPy-type per-component emissivities fitted jointly with dust. Include the stationary part and the solar-centric excess, and stop treating DIRBE 100 µm as clean.
- **Low-frequency external data.** Ingest QUIJOTE (with a 5% calibration prior), C-BASS, S-PASS (with a Faraday mask or model) and CLASS 40 GHz, each with its own transfer function. Allow curved or spatially varying synchrotron, and templates for local structures.
- **Ground-based maps.** Support filter/transfer-matrix-aware likelihoods; CLASS shows filtering is the dominant τ limit.
- **Sampling speed.** Address the ~4000-sample requirement for τ with derivative-based or multigrid samplers.
- **Robustness checks.** Ship both constrained and random-input simulation sets, and provide a cross-spectrum (elica-style) check alongside the pixel/Blackwell-Rao likelihoods.
- **Reionization models.** Sample flexible x_e(z) with flat-τ priors.

## Gaps

I hit the web search limit before I could check:
- C-BASS South status;
- BICEP/Keck results after BK18;
- SPIDER 280 GHz science results;
- an official LiteBIRD σ(τ) number;
- NPIPE §4.3 transfer-function values.

None of these were guessed.

## References

- 2011.05609 https://arxiv.org/abs/2011.05609 (BeyondPlanck I)
- 2011.06650 https://arxiv.org/abs/2011.06650 (noise)
- 2011.08082 https://arxiv.org/abs/2011.08082 (gain)
- 2201.03417 https://arxiv.org/abs/2201.03417 (bandpass)
- 2201.03478 https://arxiv.org/abs/2201.03478 (sidelobes)
- 2208.14276 https://arxiv.org/abs/2208.14276
- 2205.10104 https://arxiv.org/abs/2205.10104
- 2209.04437 https://arxiv.org/abs/2209.04437
- 2303.08095 https://arxiv.org/abs/2303.08095 (Cosmoglobe DR1 WMAP)
- 2310.13740 https://arxiv.org/abs/2310.13740
- 2306.15511 https://arxiv.org/abs/2306.15511
- 2507.05324 https://arxiv.org/abs/2507.05324
- 2604.16901 https://arxiv.org/html/2604.16901
- 2408.10952, 2408.11004, 2406.01491, 2205.12962 https://arxiv.org/abs/2408.11004
- Gjerløw FIR2025 talk: https://farinfrared.web.rug.nl/pres/Gjerlow-FIR2025.pdf
- 1303.5074 https://arxiv.org/abs/1303.5074 (Planck 2013 XIV, zodi)
- 2007.04997 https://arxiv.org/abs/2007.04997 (NPIPE)
- 2309.10034 https://arxiv.org/abs/2309.10034
- 2504.13254 https://arxiv.org/abs/2504.13254
- 1807.06207 https://arxiv.org/abs/1807.06207
- 1807.06209 https://arxiv.org/abs/1807.06209
- 1901.11386 https://arxiv.org/abs/1901.11386 (SRoll2)
- 1908.09856 https://arxiv.org/abs/1908.09856
- 2103.14378 https://arxiv.org/abs/2103.14378
- 2005.05600 https://arxiv.org/abs/2005.05600
- 2603.22454 https://arxiv.org/abs/2603.22454 (elica)
- 2504.16932 https://arxiv.org/abs/2504.16932
- 2601.09644 https://arxiv.org/abs/2601.09644
- 2607.26373 https://arxiv.org/abs/2607.26373
- 2501.11904 https://arxiv.org/abs/2501.11904 (CLASS τ)
- 2404.17567 https://arxiv.org/abs/2404.17567
- 2301.05113 https://arxiv.org/abs/2301.05113 (QUIJOTE)
- 2203.04861 https://arxiv.org/abs/2203.04861
- 2511.14572 https://arxiv.org/abs/2511.14572
- 2202.10411 https://arxiv.org/abs/2202.10411 (C-BASS)
- 2606.21334 https://arxiv.org/abs/2606.21334
- 2605.23489 https://arxiv.org/abs/2605.23489
- 1802.01145 https://arxiv.org/abs/1802.01145 (S-PASS)
- 2103.13334 https://arxiv.org/abs/2103.13334 (SPIDER)
- 2110.00483 https://arxiv.org/abs/2110.00483 (BK18)
- 2512.10613 https://arxiv.org/abs/2512.10613
- 2202.02773 https://arxiv.org/abs/2202.02773 (LiteBIRD)
- 2312.09001 https://arxiv.org/abs/2312.09001
- 2506.22217 https://arxiv.org/abs/2506.22217
- 2407.01438 https://arxiv.org/abs/2407.01438 (Taurus)
- 2512.16869 https://arxiv.org/abs/2512.16869
- Nordita 2026 talks: https://indico.fysik.su.se/event/9149/contributions/15464/ (Lunde, "Commander4: Status and future") and https://indico.fysik.su.se/event/9149/contributions/15439/ (Giardiello, LiteBIRD τ forecast)

**Possibly worth remembering (for /areas/commander-improvement.md):**
- Commander4 has a Python interface with C++ for the low-level code, and scales across many compute nodes.
- As of August 2026 it can run preliminary analyses of SO and LiteBIRD simulations (Lunde talk, Nordita).
- Cosmoglobe's planned projects include an HFI reanalysis, AKARI, FIRAS and CHIPASS.
