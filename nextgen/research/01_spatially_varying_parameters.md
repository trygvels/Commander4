> Raw report from a research subagent, 2026-10-25. Most sources were read at abstract level.
> Treat every number as a lead: read the paper before citing it. Curated version: `../LITERATURE.md`.

## Spatially varying spectral parameters for Commander4: literature review (2018 to 2026) and recommendations

Everything below comes from arXiv abstracts or full text I retrieved, plus two codebases I read (Cosmoglobe/Commander3 and trygvels/Commander4 HEAD). Where I could not verify something, I say so. The web-search budget ran out near the end, so a few sub-topics are thin; these are flagged.

### 1. Ways to let spectral parameters vary across the sky

**Fixed regions and multi-resolution pixels**
- BeyondPlanck XV (Svalheim+ 2023) fit β_s in 4 large regions, smoothed each proposal with a 10° beam to suppress edge effects, and fit β_d as one full-sky value.
  - Only the Spur and Plane regions were data-constrained (β_s = −3.17±0.06 and −3.03±0.07). High latitude was prior-dominated.
  - Widening the prior to σ = 0.2 let β_s reach about −2.7. The CMB then absorbed the excess foreground and biased τ.
- Commander3 already has two relevant modes in `comm_nonlin_smod.f90`:
  - An alm sampler for β with a Gaussian-beam prior per multipole (`ALMSAMP_PRIOR_FWHM`). The prior width at each ℓ is the `p_gauss` σ times exp(−ℓ(ℓ+1)σ_b²/2).
  - A pixel-region mode. It smooths the region map, clips it to the uniform bounds with the in-code comment "in case of ringing", then converts to alms.
  - So ringing was a real problem in practice. The prior hyperparameters were fixed, never sampled.
- Rizzieri+ 2025 (fgbuster) compared multi-patch, multi-resolution (nside [64, 8, 4] for β_d, T_d, β_s), latitude-adaptive and clustered subsets:
  - Statistical residuals depend only weakly on the foregrounds. Systematic residuals depend strongly on both the foregrounds and the β resolution.
  - The latitude-adaptive configuration optimized for LiteBIRD gives an r bias of 0.8×10⁻³ (d1s1) and 3×10⁻³ (d10s5). Applied to a PICO-like setup it gives about 2×10⁻².
  - 200 "ideal" clusters drawn from the true β maps give negligible bias.
- Errard & Stompor 2019: ignoring spatial variability biases r by O(0.01).

**Adaptive clustering**
- Puglisi+ 2022: spectral clustering that combines pixel proximity with measurement similarity weighted by the errors. It can also cluster on HI cloud counts (N_c).
- Carones+ 2023 (MC-NILC): clusters defined by the ratio of B-mode maps at two frequencies.
- Grumitt+ 2020: mean-shift clustering into 171 regions.
- Kabalan+ 2026 (FURAX, JAX on GPU): a K-means grid search that minimizes r̂ + σ(r) on simulations. It reduces the 68% upper limit on r by about 30% relative to the fixed multi-resolution configuration.
  - The chosen patch counts are about 7000 to 9700 for β_d, 300 to 2500 for T_d, and 150 to 300 for β_s.
  - The selection uses simulation truth, so it is not yet usable on real data.

**Hierarchical and Gaussian-random-field priors**
- Grumitt+ 2020 (NUTS): per-pixel β drawn from a Normal hyper-distribution per region, with sampled hypermean and hyperstd and half-normal hyperpriors, in a non-centred parameterization.
  - With input r = 5×10⁻³ (LiteBIRD + C-BASS): pooled regions gave r = (9.0±1.1)×10⁻³, the hierarchical model gave (5.2±1.0)×10⁻³.
  - The authors report funnel geometry, and the model could not run on LiteBIRD alone because of divergences.
- de Belsunce+ 2022: a Gaussian-kernel prior covariance on β with a fixed θ_corr = 5° (not inferred), solved by Newton-Raphson. They warn that results are "strongly prior dominated."
- I found no paper that infers a β power spectrum or correlation length as a free hyperparameter inside a CMB Gibbs sampler. That appears to be an open gap.

**Wavelet, sparsity and total-variation priors**
- Only blind or semi-blind methods appeared: L-GMCA (Bobin+ 2013), HGMCA (Wagner-Carena+ 2020) and needlet ILC variants.
- I found no total-variation or sparsity prior on β maps (not verified as absent).

### 2. Moment expansion

- Chluba+ 2017: averaging over the line of sight and the beam distorts SEDs. Expanding around pivot values gives derivative templates (∂_β f, ∂²_β f, cross terms) with "moment" amplitudes, and the distortion is scale dependent.
- Power-spectrum pipelines:
  - Azzoni+ 2021 (Simons Observatory): unbiased r with σ(r) ≈ 0.003, and r < 0.06 from BK data.
  - Vacher+ 2022 (LiteBIRD): moments in β alone are insufficient. Moments in β and T give unbiased r with σ(r) = 8.8×10⁻⁴, but correlations between moments and the CMB inflate σ(r).
- Map and needlet pipelines:
  - cMILC (Remazeilles+ 2021) and ocMILC (Carones & Remazeilles 2024) deproject moments, giving unbiased r for four PICO foreground models.
  - QUIJOTE XIX (Adak+ 2025) combined moments with cILC: β_s = −3.05±0.01 in the plane and −3.13±0.02 at high latitude.
  - The Delta-map method (Ichiki+ 2019) works to first order in δβ and gives r ≲ 4×10⁻⁴ in noiseless simulations. It can model decorrelation with separate Q and U parameters.
  - Hybrid map-Cℓ method (Azzoni+ 2023): r_FG ≲ 5×10⁻⁴ for a 20 to 30% increase in σ.
- Physics with spin moments:
  - Vacher+ 2023 and 2025 show frequency-dependent polarization angles and E/B ratios, and 3D dust complexity that different methods respond to differently.
  - Guillet+ 2025: Planck SRoll2 data prefer temperature-only fluctuations (σ_T ≈ 2.7 K). Polarized and intensity fluctuations correlate above 85%.
- Pipeline comparisons: complex foregrounds give 1 to 2σ (CMB-S4 2025) or more than 2σ (Simons Observatory, Wolz+ 2024) bias with the default pipelines. Marginalizing over foreground residuals fixes this at a moderate cost in σ(r).

### 3. Foreground models and data that break degeneracies

**Synchrotron**
- QUIJOTE MFI at 11 and 13 GHz (de la Hoz+ 2023):
  - Mean β_s = −3.08 with dispersion 0.13, mapped at 2° smoothing.
  - Uniform curvature c_s = −0.0797±0.0012, but the curved model is not significantly preferred.
  - The prior controls β_s outside the bright regions.
- MFI2 forecast (Almeida+ 2025): about 10× better β_s and 5× better curvature. Pixel-level curvature is detectable only if |C_s| ≳ 0.14 in the plane.
- S-PASS at 2.3 GHz with WMAP (Fuskeland+ 2021, 15° patches): β_s steepens from about −2.8 near the plane to about −3.3 at high latitude.
- Weiland+ 2022 (about 3.5° scales, 1.4 to 33 GHz): β steepens by 0.25 to 0.4 away from the plane, and β = −3.25±0.04 in the BK field.
- C-BASS at 5 GHz: decisive in the Grumitt simulations (above). Shaw+ 2026 finds local structures (North Polar Spur, Fan region) must be modelled explicitly.
- Polarized AME is degenerate with β_s:
  - BeyondPlanck XVI: p_AME < 0.6% if β_s is steeper than −3.1, but a nominal 2.5±1.0% if β_s = −3.0.
  - QUIJOTE XVIII: limits below about 1% in most of the bright regions studied.

**Dust**
- Planck 2018 XI: β_d = 1.53±0.02 with no detected decorrelation.
- SPIDER 2025: two sub-regions give β_d = 1.52 and 1.09, which differ at 3.9σ.
- Extending LiteBIRD's high-frequency limit from 400 to 600 GHz cuts σ(r) by 30 to 50% when dust β and T vary (Fuskeland+ 2023).
- Cosmoglobe DR2 VI and VII model dust as hot, cold, nearby and Hα-correlated components with uniform SEDs. These reproduce published β/T maps economically, so apparent spatial variation may be a mixture of uniform components.

**Summary of which data help where**: 2 to 20 GHz data (C-BASS, QUIJOTE, S-PASS) constrain β_s at 1 to 4° scales. Frequencies above 500 GHz and FIRAS/DIRBE constrain T_d. I did not retrieve results for GMIMS, LWA, BLAST or PIPER (not verified).

### 4. Machine-learning and data-driven priors

- Diffusion models as dust priors for posterior sampling: Heurtel-Depeiges+ 2023, and GDiff 2024, which samples cosmological parameters by Gibbs alternation.
- Scattering transforms: dust denoising (Delouis+ 2022), polarized dust separation (Tsouros+ 2026), dust/CIB separation (Sinha+ 2026).
- ForSE+ GANs, for simulations only.
- NILC combined with simulation-based inference for r (Duivenvoorden+ 2025).
- I found no ML prior used on spectral-index fields inside a Bayesian sampler, and no SED emulator used in one. Commander3 uses spline lookup tables limited to two dimensions. JAX autodiff (FURAX, about 100× faster than fgbuster's optimizer) is the practical route to gradients.

### 5. Optimal resolution and prior-volume effects

**Resolution**
- Systematic and statistical residuals move in opposite directions as patches get finer. Statistical residuals bottom out at the noise level (Rizzieri+ 2025).
- The optimum depends on the experiment and on the foreground model. Empirically β_d wants fine resolution, T_d coarser and β_s coarsest (Kabalan+ 2026).
- No published analytic optimum exists. A Wiener filter of the β field (hierarchical GRF prior) gives the minimum-mean-squared-error (MMSE) trade-off automatically (my inference).

**Prior volume and Jeffreys priors**
- Eriksen+ 2008: a uniform β prior produced a monopole bias of more than 4σ. The Jeffreys prior √F_ββ removes it. Planck 2015 X adopts it "to suppress prior volume effects."
- Stompor+ 2009: the spectral (profile) likelihood has the same maximum as the full joint likelihood.
- BeyondPlanck XIII, Eq. 23: the amplitude-marginal likelihood adds +½ ln|(AᵀN⁻¹A)⁻¹|. The ridge (profile) form omits that term.
- My own derivation: a joint Jeffreys prior on (a, β) contributes det(AᵀN⁻¹A)^{1/2}, which cancels the marginal's determinant. So profile likelihood × conditional Jeffreys for β is roughly "joint Jeffreys."

**Commander4 code check (HEAD of trygvels/Commander4, `compsep/mcmc.py`)**
- Each step makes a symmetric random-walk proposal, re-solves the amplitudes by CG, and accepts on the data χ² only.
- If the re-solve is a draw with fluctuations (`optimize: false`, the default), this is not an exact Metropolis-Hastings step for either the joint or the marginal target. The correct joint step uses the Eq. 23 marginal ratio, then draws a.
- If the re-solve is the MAP (`optimize: true`), it is the profile likelihood, but missing the aᵀS⁻¹a prior term.
- This is worth checking against your recent acceptance fix.

## Recommendations for Commander4, ranked

1. **Hierarchical GRF prior on each spectral field.** Write β(n) = β̄ + δβ(n), with δβ ~ GRF.
   - Use C_ℓ^β = A(ℓ/ℓ₀)^α, or a Matérn field with correlation length λ. Sample (β̄, A, α or λ) in the Gibbs chain with weakly informative hyperpriors (half-normal on A).
   - Use log T_d rather than T_d.
   - The posterior mean is then locally adaptive: data-dominated in the plane, shrunk to β̄ at high latitude. This is "as fine as S/N allows" without hand-drawn boundaries.
   - It generalizes Grumitt's per-region model and Commander3's fixed alm prior.
2. **Avoid ringing by defining the field in pixel space.** Use a sparse GMRF/SPDE precision on HEALPix neighbours at nside_β (for example 64 for β_d and T_d, 16 to 32 for β_s), not truncated alms. If you keep alms, set lmax well beyond the prior's cutoff.
3. **Use a gradient-based sampler** (HMC/NUTS or MALA) in whitened, non-centred coordinates, with JAX gradients of the SEDs. Base the likelihood on a low-resolution, amplitude-marginal (or profile) spectral likelihood with beam-aware smoothing (Rizzieri+ 2025). Alternate this with the CG amplitude solve, and update the hyperparameters with joint marginal moves.
4. **Fix the likelihood definition first.** Use the Eq. 23 marginal with a stochastic log-det (Hutchinson/Lanczos), or deliberately use the profile likelihood with an explicit Jeffreys term. Calibrate either choice with simulation-based calibration (SBC) coverage tests.
5. **Add first-order moment components** (∂_β f times an amplitude map) as linear CG components at full resolution with Gaussian priors. This absorbs variation below the β-field resolution and line-of-sight averaging. Use intensity β/T maps as the prior means, since polarized and intensity fluctuations correlate above 85%.
6. **Add data where it matters.** C-BASS, QUIJOTE MFI and Faraday-corrected S-PASS for β_s. HFI 545/857 and DIRBE for T_d. Keep synchrotron curvature as one global parameter, and treat polarized AME as a nuisance with a prior of p ≲ 1%.
7. **Validate on simulations.** Use PySM d1s1, d10s5 and d12 skies and report r bias and CMB residuals against the fixed nside [64, 8, 4] baseline. Clustering the posterior-mean β map is a fallback if fixed regions are needed for production runs.

## References
- Rizzieri+ 2025, arXiv:2510.08534, https://arxiv.org/abs/2510.08534
- Kabalan+ 2026, arXiv:2604.08463, https://arxiv.org/abs/2604.08463
- Puglisi+ 2022, arXiv:2109.11562, https://arxiv.org/abs/2109.11562
- Carones+ 2023, arXiv:2212.04456, https://arxiv.org/abs/2212.04456
- Grumitt+ 2020, arXiv:1910.14170, https://arxiv.org/abs/1910.14170
- de Belsunce+ 2022, arXiv:2205.13968, https://arxiv.org/abs/2205.13968
- Svalheim+ (BeyondPlanck XV), arXiv:2011.08503, https://arxiv.org/abs/2011.08503
- Andersen+ (BeyondPlanck XIII), arXiv:2201.08188, https://arxiv.org/abs/2201.08188
- Galloway+ (BeyondPlanck III), arXiv:2201.03509, https://arxiv.org/abs/2201.03509
- Watts+ (Cosmoglobe DR1 III), arXiv:2310.13740, https://arxiv.org/abs/2310.13740
- Errard & Stompor 2019, arXiv:1811.00479, https://arxiv.org/abs/1811.00479
- Rizzieri+ 2025 (beams), arXiv:2409.12944, https://arxiv.org/abs/2409.12944
- Wagner-Carena+ 2020, arXiv:1910.08077, https://arxiv.org/abs/1910.08077
- Bobin+ 2013, arXiv:1206.1773, https://arxiv.org/abs/1206.1773
- Chluba+ 2017, arXiv:1701.00274, https://arxiv.org/abs/1701.00274
- Azzoni+ 2021, arXiv:2011.11575, https://arxiv.org/abs/2011.11575
- Vacher+ 2022, arXiv:2111.07742, https://arxiv.org/abs/2111.07742
- Remazeilles+ 2021, arXiv:2006.08628, https://arxiv.org/abs/2006.08628
- Rotti & Chluba 2021, arXiv:2006.02458, https://arxiv.org/abs/2006.02458
- Carones & Remazeilles 2024, arXiv:2402.17579, https://arxiv.org/abs/2402.17579
- Vacher+ 2025, arXiv:2411.11649, https://arxiv.org/abs/2411.11649
- Vacher+ 2023, arXiv:2210.14768, https://arxiv.org/abs/2210.14768
- Guillet+ 2025, arXiv:2510.18305, https://arxiv.org/abs/2510.18305
- Azzoni+ 2023, arXiv:2210.14838, https://arxiv.org/abs/2210.14838
- Ichiki+ 2019, arXiv:1811.03886, https://arxiv.org/abs/1811.03886
- Minami & Ichiki 2023, arXiv:2212.01773, https://arxiv.org/abs/2212.01773
- Wolz+ 2024, arXiv:2407.06891, https://arxiv.org/abs/2407.06891
- CMB-S4 2025, arXiv:2502.04300, https://arxiv.org/abs/2502.04300
- Duivenvoorden+ 2025, arXiv:2512.16869, https://arxiv.org/abs/2512.16869
- Mustafa+ 2026, arXiv:2603.11026, https://arxiv.org/abs/2603.11026
- de la Hoz+ 2023 (QUIJOTE VIII), arXiv:2301.05117, https://arxiv.org/abs/2301.05117
- Adak+ 2025 (QUIJOTE XIX), arXiv:2510.17761, https://arxiv.org/abs/2510.17761
- Almeida+ 2025, arXiv:2511.14572, https://arxiv.org/abs/2511.14572
- González-González+ 2024 (QUIJOTE XVIII), arXiv:2409.03418, https://arxiv.org/abs/2409.03418
- Fuskeland+ 2021, arXiv:1909.05923, https://arxiv.org/abs/1909.05923
- Weiland+ 2022, arXiv:2203.11445, https://arxiv.org/abs/2203.11445
- Shaw+ 2026, arXiv:2605.23489, https://arxiv.org/abs/2605.23489
- Herman+ (BeyondPlanck XVI), arXiv:2201.03530, https://arxiv.org/abs/2201.03530
- SPIDER 2025, arXiv:2407.20982, https://arxiv.org/abs/2407.20982
- Planck 2018 XI, arXiv:1801.04945, https://arxiv.org/abs/1801.04945
- Fuskeland+ 2023, arXiv:2302.05228, https://arxiv.org/abs/2302.05228
- Cosmoglobe DR2 V, VI, VII: arXiv:2601.07818, 2601.10640, 2601.07822, https://arxiv.org/abs/2601.10640
- Heurtel-Depeiges+ 2023, arXiv:2310.16285, https://arxiv.org/abs/2310.16285
- GDiff 2024, arXiv:2402.19455, https://arxiv.org/abs/2402.19455
- Delouis+ 2022, arXiv:2207.12527, https://arxiv.org/abs/2207.12527
- Tsouros+ 2026, arXiv:2602.04528, https://arxiv.org/abs/2602.04528
- Sinha+ 2026, arXiv:2603.06110, https://arxiv.org/abs/2603.06110
- ForSE+ 2024, arXiv:2406.14519, https://arxiv.org/abs/2406.14519
- Stompor+ 2009, arXiv:0804.2645, https://arxiv.org/abs/0804.2645
- Eriksen+ 2008, arXiv:0709.1058, https://arxiv.org/abs/0709.1058
- Planck 2015 X, arXiv:1502.01588, https://arxiv.org/abs/1502.01588
- MICMAC: arXiv:2405.18365, 2308.03412, 2502.15511, https://arxiv.org/abs/2405.18365

**Code inspected** (clones in /tmp/claude-0/-home-claude/1540e426-15b0-5fbf-8cc9-6e636ee967c6/scratchpad):
- c3/Commander/commander3/src/comm_nonlin_smod.f90, lines 385 to 560
- c3/Commander/commander3/src/comm_diffuse_comp_smod.f90, lines 1237 to 1257
- c4/Commander4/src/commander4/compsep/mcmc.py, lines 138 to 243
