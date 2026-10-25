> Raw report from a research subagent, 2026-10-25. Most sources were read at abstract level.
> Treat every number as a lead: read the paper before citing it. Curated version: `../LITERATURE.md`.

## Adjacent-field methods that could improve Commander4: survey (Oct 2026)

Scope: I took every number below from the retrieved arXiv abstract or HTML pages, or from software docs. Most papers I read only at abstract level. Where a figure came from a full-text extraction, or I could not check it, I say so.

### 1. Samplers for very high-dimensional posteriors
- **MCLMC (2303.18221)**: on a φ⁴ lattice it was 12 to 32× faster than HMC. On hierarchical models it needed about 3× fewer gradients than NUTS.
- **MAMS (2503.01707, NeurIPS 2025)** is MCLMC with an MH correction, so it is asymptotically unbiased. It auto-tunes and reports "up to 7×" speedup over NUTS.
- **Field-level evidence**:
  - Bayer, Seljak and Modi (2307.09504): MCLMC was more than 10× more efficient than HMC at 2.6×10⁵ dimensions.
  - Simon, Lanusse and de Mattia (2504.20130): unadjusted MCLMC scales past 10⁶ dimensions and beats the adjusted samplers by more than an order of magnitude. Preconditioning the latent variables "matters." At 128³ it took under 4 h on 8 GPUs, against more than 80 h for NUTS.
  - **Flinch.jl (2510.26691)** is the closest analogue to Commander4. It does curved-sky masked CMB field-level inference with more than 10⁶ parameters at nside 512, using ducc SHTs with hand-written adjoints (a gradient costs about 2.1× a forward pass). MCLMC got roughly 80 to 100× more ESS per gradient than NUTS (7 h against about 6 days).
- **Mass-matrix adaptation (2603.18845, nutpie)**: minimising Fisher divergence beat the Stan/PyMC diagonal adaptation. The low-rank-plus-diagonal form gave a median 4× gain over 114 models.
- **Dimension-robust function-space MCMC**:
  - pCN has a dimension-independent spectral gap, while the random-walk gap collapses (Hairer, Stuart and Vollmer, 1112.1392).
  - Non-centred hierarchical parameterisations plus dimension-robust samplers keep hyperparameter mixing stable as resolution increases (Chen et al., 1803.03344).
  - Geometric MCMC in a data-informed subspace, with pCN on the complement, gave "up to two orders of magnitude" over pCN (Beskos et al., 1606.06351).
  - In Commander these matter only for the non-Gaussian fields, such as spectral-index maps with Gaussian priors. The amplitude draw is already exact.
- **RMHMC (0907.1100)**: the metric costs a lot per step. The practical version is a fixed Fisher or Wiener-type mass matrix.
- **flowMC (2211.06397)**: it works well at low dimension, for example Jim's gravitational-wave runs in about a minute (2302.05333). Grenioux et al. (2302.04763) find flow proposals reliable only "up to moderately high dimensions." It is unsuitable for the field blocks and plausible only for a handful of multimodal global parameters.

### 2. Gibbs-specific improvements
- **Rescaling moves**:
  - Jewell et al. (0807.0624): propose a new C_ℓ, rescale the sky deterministically, and get acceptance near 1 at low S/N.
  - Racine et al. (1512.06619): rescale only the fluctuation part, which avoids the high-S/N inefficiency.
  - Eskilt et al. (A&A 678, A169; 2306.15511) implemented this for cosmological parameters in Commander3: about 20 CPU-h per sample and a correlation length of about 100, so about 2000 CPU-h per independent sample. They say a further 10× reduction is needed and suggest multigrid preconditioners or derivative-based sampling.
- **Ducrocq, Chopin, Errard and Stompor (2111.07664)** is the most directly relevant paper here. ASIS (interweaving) gave about parity with centred Gibbs. "Centered overrelax" (auxiliary variables plus overrelaxation) and "Centered 1" gave large ESS/s gains over centred Gibbs:

  | Sampler | EE | BB |
  |---|---|---|
  | Centered overrelax, near full sky | 6.9× | 36× |
  | Centered overrelax, cut sky | 12.7× | 76× |
  | Centered 1, cut sky (BB only) | | 184× |

  These per-sampler figures come from a full-text extraction. The abstract confirms only the 1 to 2 orders of magnitude.
- **Partially collapsed Gibbs (van Dyk and Park, JASA 2008)** converges at least as fast as the parent Gibbs sampler. Trimming steps can break the stationary distribution, so step order matters.
- **Collapsing amplitudes for spectral indices**:
  - Stompor et al. (0804.2645) give an analytic spectral likelihood for β whose maximum equals the full-problem maximum likelihood.
  - Grumitt et al. (1910.14170): NUTS was about 10³× faster than MH for β, and a hierarchical per-pixel β model reduced the bias on r.
- **MUSE (2112.09354)**: an approximate marginal for global parameters, run on SPT-3G-sized data with 6×10⁶ latent dimensions. The implicit-differentiation version (2209.10512) was 5 to 397× faster than HMC on test problems.
- **Messenger fields** are covered in section 4.
- **Gibbs elsewhere**: radio cosmology is converging on the same Gibbs and constrained-realisation design. Examples are Hydra (2211.05088, 2403.13767), joint gain/sky Gibbs for MeerKLASS (2509.10992), and a joint sky, flux-scale and spectral-index Gibbs sampler (2509.11894).

### 3. Spatial statistics
- **SPDE (Lindgren, Rue and Lindström, JRSS-B 2011)** turns Matérn fields into sparse GMRFs on any triangulation, explicitly including the sphere (global temperature data, 15,182 vertices). Building the precision costs O(n). A sparse Cholesky costs O(n^{3/2}) in 2D. The 10-year review is 2111.01084.
- **NNGP (1406.7343)** has flops linear in n.
- **GraphGP (2606.11402)** runs Vecchia on GPU in JAX/CUDA to nearly 10⁹ parameters. It provides the log-determinant and kernel derivatives.
- **Uses in astrophysics**:
  - INLA for IFU galaxy maps, 721 cubes (1802.06280).
  - INLA for radiative-transfer images, about 5× faster (2211.02602).
  - A GP 3D dust map in spherical coordinates at Nside 256 with 661M degrees of freedom (iterative charted refinement plus MGVI, 2308.01295; NIFTy.re is 2402.16683).
  - Hierarchical GP kernels for 21-cm foregrounds cut residual standard deviation by up to 30% (2407.11296).
- I found **no SPDE, GMRF or INLA application on HEALPix for CMB foreground parameters**. That is a gap and an opportunity.

### 4. Linear algebra
- **Trace and log-determinant estimators**:
  - Hutchinson needs O(1/ε²) matrix-vector products; Hutch++ needs O(1/ε) and is optimal (2010.09649).
  - XTrace reports errors "orders of magnitude smaller" than Hutch++ for the same budget (2301.07825).
  - Stochastic Lanczos quadrature (Ubaru, Chen and Saad, SIMAX 2017) needs about 4√κ·log(K/ε) Lanczos steps and about 24·log(2/η)/ε² probes. It was demonstrated up to n ≈ 10⁶ (full-text extraction).
  - Dorn and Enßlin (1504.02661) give a stochastic log-determinant for implicit operators in the IFT context.
  - GPyTorch's BBMM reduces exact GP cost from O(n³) to O(n²) (1809.11165).
  - Preconditioning reduces log-determinant estimator variance and speeds training by up to 10× (Wenger et al., 2107.00243).
  - `matfree` provides differentiable Lanczos in JAX (2405.17277).
  - Because the step count scales as √κ, the existing CG preconditioner directly sets the cost of the log-determinant.
- **CMB preconditioners and solvers**:
  - Multilevel solver (Seljebotn et al., 1308.5299): 3 W-cycles reached under 1 µK error per pixel on Planck 143 GHz in about 6 min on 16 cores. CG never reached that accuracy.
  - Pseudo-inverse preconditioner (1710.00621): 2 to 3× faster than diagonal on the full sky.
  - Two-level map-making preconditioner (1408.3048): up to 5× fewer iterations and 4× less time.
  - Papež, Grigori and Stompor (1803.03462) show messenger methods are fixed-point iterations on a preconditioned system, so CG is "comparable and typically better." They also find PCG performance depends on the starting vector.
  - Conflict: the dual-messenger paper (1702.08852) claims 3 to 4× faster than PCG. The two sets of results disagree.
  - I could not retrieve the original Smith et al. 2007 multigrid C⁻¹ preconditioner.
- **ML-learned solvers**:
  - WienerNet is about 1000× faster than CG but approximate (flat sky, 1905.05846).
  - DeepWiener extends this to polarisation and reports only qualitative accuracy (2412.10580).
  - Learned preconditioners that keep CG exactness: NeuralIF (2305.16368) and UNO-CG (2508.02681). No figures were in the abstracts.
- **Recycling across a sequence of systems**: Krylov recycling and deflation for sequences of systems inside MCMC exists (py-deflation code; survey 2001.10347). It fits Gibbs well, because the operator changes slowly between iterations.

### 5. Transferable practice from other fields
- **Large-scale structure**:
  - BORG combines block sampling with HMC (1806.11117).
  - Manticore validates with posterior predictive checks on power spectra, bispectra and initial-condition Gaussianity (2505.10682).
  - The JAX benchmark (2504.20130) and Flinch both show MCLMC plus careful preconditioning working at 10⁶ dimensions or more.
- **Pulsar timing arrays**:
  - NANOGrav's Discovery is a JAX/GPU likelihood built from composable Woodbury kernels. It is a pattern for 1/f noise and template marginalisation. No benchmarks are published.
  - Variational inference plus flows analyses the NANOGrav 15-year data in tens of minutes, against hours to days (2405.08857).
- **Gravitational waves**: Jim (JAX plus flowMC) takes about 1 min of sampling for binary black holes (2302.05333) and 21 to 33 min end to end for binary neutron stars (2404.11397).
- **21-cm**: GPR foreground removal is now improved with learned VAE kernels (2311.16633) and hierarchical kernels (2407.11296).
- **Differentiable SHTs and CMB JAX tooling**:
  - s2fft (2311.14670): up to 400× faster than C codes and supports HEALPix. Machine precision holds only for sampling-theorem grids, which HEALPix is not.
  - cuHPX (2510.01785): over 20× faster, with analytic gradients.
  - jax-healpy (Chanial et al.).
  - Furax (2603.19600): JAX operators for map-making and component separation.
  - GPU component separation (2604.08463): about 100× faster than fgbuster's TNC optimiser.
- **Main lesson**: write the SHT adjoints by hand, as Flinch does, rather than differentiating through the SHT code.

### 6. Validation practice
- **Simulation-based calibration and coverage**:
  - SBC (Talts et al., 1804.06788).
  - Modrák et al. (2211.02383) add data-dependent test quantities, especially the joint likelihood. These catch a posterior that equals the prior, which plain SBC misses.
  - SBC for Bayes factors needs several hundred simulations (2508.11814).
  - TARP (2302.03026) tests coverage in high dimension without evaluating densities.
- **Large collaborations**:
  - Planck FFP8 (1509.06348): about 10⁴ Monte Carlo realisations at about 25M CPU-h.
  - BeyondPlanck IV (2209.04437) separates Bayesian simulations (constrained) from frequentist ones (random), and validated gain, noise and map-making on 30 GHz simulations.
  - BeyondPlanck reports a χ² probability-to-exceed (32%) as a posterior predictive check.
  - A 2026 field-level BAO paper checks coverage with 1000 mocks (2603.15732).
- **Cost of a full Gibbs sample**:

  | Configuration | CPU-h per sample |
  |---|---|
  | Commander3, BeyondPlanck (2201.03509) | 170 |
  | Cosmoglobe DR1, WMAP plus LFI (2303.08095) | 812 |
  | LiteBIRD (2507.05324) | about 3000 |

  Full-resolution SBC is therefore infeasible. It has to be done at reduced resolution.

### Ranked opportunities for Commander4
1. **Faster CG solves.** Warm starts, Krylov recycling across Gibbs steps, and better preconditioners (multilevel or pseudo-inverse, or a learned preconditioner inside CG).
   - Gain: 2 to 5× from preconditioners alone (published); possibly much more with multilevel.
   - Cost: moderate.
   - Risk: low, because the solves stay exact.
2. **Collapsed and joint moves for β against amplitudes.** Marginalise amplitudes analytically in the β step, or use joint (β, a) proposals, following partially collapsed Gibbs.
   - Gain: potentially orders of magnitude in β correlation length (by analogy with Grumitt's 10³ and the rescaling moves).
   - Cost: low for per-pixel β; needs a log-determinant once there are spatial priors.
   - Risk: step-order validity.
3. **Auxiliary-variable and overrelaxed C_ℓ moves** (Ducrocq et al.) plus Racine rescaling.
   - Gain: 10 to 100× ESS/s at low S/N (on C_ℓ-only problems).
   - Cost: low to moderate.
   - Risk: untested together with foreground and instrument steps.
4. **MAMS or MCLMC instead of NUTS**, with Fisher or low-rank preconditioning and tuning carried over between Gibbs iterations.
   - Gain: 3 to 100× ESS per gradient (problem-dependent).
   - Cost: low if the gradients are in JAX (BlackJAX).
   - Risk: unadjusted MCLMC is biased, so use MAMS. Re-tuning on every Gibbs conditional costs extra.
5. **Stochastic log-determinants** (Hutch++, XTrace or SLQ with the CG preconditioner) for marginal likelihoods, model comparison, and collapsed hyperparameter steps.
   - Gain: enables things Commander cannot do now.
   - Cost: probes × Lanczos steps × SHTs.
   - Risk: noise in the log-determinant biases MH unless handled carefully.
6. **GMRF/SPDE priors for spectral-index maps on HEALPix.**
   - Gain: cheap exact log-determinants and non-stationarity.
   - Cost: moderate.
   - Risk: no CMB precedent.
7. **GPU SHTs** (s2fft, cuHPX, ducc with custom adjoints).
   - Gain: about 10 to 400× on the SHT kernels (vendor claims).
   - Cost: high.
   - Risk: HEALPix quadrature precision and multi-node scaling.
8. **MUSE** as a fast cross-check on global parameters.
   - Gain: large speedup.
   - Cost: moderate.
   - Risk: approximate (Gaussianised).
9. **Validation suite**: low-resolution SBC with joint-likelihood test quantities, TARP, and posterior predictive χ² maps.
   - Gain: credibility of the error propagation.
   - Cost: hundreds of low-resolution runs.
   - Risk: low.

flowMC is not recommended except for a few multimodal instrument parameters. Two items need a full-text read before relying on them: the SLQ bounds and the per-sampler ratios in Ducrocq et al.

### References
- 2303.18221 https://arxiv.org/abs/2303.18221
- 2503.01707 https://arxiv.org/abs/2503.01707
- 2307.09504 https://arxiv.org/abs/2307.09504
- 2504.20130 https://arxiv.org/abs/2504.20130
- 2510.26691 https://arxiv.org/abs/2510.26691
- 2603.18845 https://arxiv.org/abs/2603.18845
- 1112.1392 https://arxiv.org/abs/1112.1392
- 1202.0709 https://arxiv.org/abs/1202.0709
- 1803.03344 https://arxiv.org/abs/1803.03344
- 1606.06351 https://arxiv.org/abs/1606.06351
- 0907.1100 https://arxiv.org/abs/0907.1100
- 2211.06397 https://arxiv.org/abs/2211.06397
- 2302.04763 https://arxiv.org/abs/2302.04763
- 0807.0624 https://arxiv.org/abs/0807.0624
- 1512.06619 https://arxiv.org/abs/1512.06619
- 2306.15511 https://arxiv.org/abs/2306.15511
- 2111.07664 https://arxiv.org/abs/2111.07664
- van Dyk & Park, JASA 2008: https://www.ma.imperial.ac.uk/~dvandyk/Research/08-jasa-pcg.pdf
- 0804.2645 https://arxiv.org/abs/0804.2645
- 1910.14170 https://arxiv.org/abs/1910.14170
- 2112.09354 https://arxiv.org/abs/2112.09354
- 2209.10512 https://arxiv.org/abs/2209.10512
- 2211.05088 https://arxiv.org/abs/2211.05088
- 2403.13767 https://arxiv.org/abs/2403.13767
- 2509.10992 https://arxiv.org/abs/2509.10992
- 2509.11894 https://arxiv.org/abs/2509.11894
- Lindgren, Rue & Lindström 2011: https://joeid.folk.ntnu.no/StatGr/Lindgren_etal.pdf
- 2111.01084 https://arxiv.org/abs/2111.01084
- 1406.7343 https://arxiv.org/abs/1406.7343
- 2606.11402 https://arxiv.org/abs/2606.11402
- 1802.06280 https://arxiv.org/abs/1802.06280
- 2211.02602 https://arxiv.org/abs/2211.02602
- 2308.01295 https://arxiv.org/abs/2308.01295
- 2402.16683 https://arxiv.org/abs/2402.16683
- 2407.11296 https://arxiv.org/abs/2407.11296
- 2010.09649 https://arxiv.org/abs/2010.09649
- 2301.07825 https://arxiv.org/abs/2301.07825
- Ubaru, Chen & Saad (SLQ): https://www-users.cse.umn.edu/~saad/PDF/ys-2016-04.pdf
- 1504.02661 https://arxiv.org/abs/1504.02661
- 1809.11165 https://arxiv.org/abs/1809.11165
- 2107.00243 https://arxiv.org/abs/2107.00243
- 2405.17277 https://arxiv.org/abs/2405.17277
- 1308.5299 https://arxiv.org/abs/1308.5299
- 1710.00621 https://arxiv.org/abs/1710.00621
- 1408.3048 https://arxiv.org/abs/1408.3048
- 1803.03462 https://arxiv.org/abs/1803.03462
- 1702.08852 https://arxiv.org/abs/1702.08852
- 1906.10704 https://arxiv.org/abs/1906.10704
- 1905.05846 https://arxiv.org/abs/1905.05846
- 2412.10580 https://arxiv.org/abs/2412.10580
- 2305.16368 https://arxiv.org/abs/2305.16368
- 2508.02681 https://arxiv.org/abs/2508.02681
- 2001.10347 https://arxiv.org/abs/2001.10347
- 1806.11117 https://arxiv.org/abs/1806.11117
- 2505.10682 https://arxiv.org/abs/2505.10682
- Discovery: https://github.com/nanograv/discovery
- 2405.08857 https://arxiv.org/abs/2405.08857
- 2302.05333 https://arxiv.org/abs/2302.05333
- 2404.11397 https://arxiv.org/abs/2404.11397
- 2311.16633 https://arxiv.org/abs/2311.16633
- 2311.14670 https://arxiv.org/abs/2311.14670
- 2510.01785 https://arxiv.org/abs/2510.01785
- jax-healpy: https://jax-healpy.readthedocs.io/
- 2603.19600 https://arxiv.org/abs/2603.19600
- 2604.08463 https://arxiv.org/abs/2604.08463
- 1804.06788 https://arxiv.org/abs/1804.06788
- 2211.02383 https://arxiv.org/abs/2211.02383
- 2508.11814 https://arxiv.org/abs/2508.11814
- 2302.03026 https://arxiv.org/abs/2302.03026
- 2209.04437 https://arxiv.org/abs/2209.04437
- 1509.06348 https://arxiv.org/abs/1509.06348
- 2603.15732 https://arxiv.org/abs/2603.15732
- 2201.03509 https://arxiv.org/abs/2201.03509
- 2303.08095 https://arxiv.org/abs/2303.08095
- 2507.05324 https://arxiv.org/abs/2507.05324
- 2604.16901 https://arxiv.org/abs/2604.16901 (Cosmoglobe overview: Commander4 is a multi-node rewrite with a Python interface over C++)
