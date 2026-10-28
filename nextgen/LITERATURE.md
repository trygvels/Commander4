# Curated literature

Compiled from the four research reports of 2026-10-25 (full reports in `research/`). Most entries
were read at abstract level by research agents; read the paper before citing a number. arXiv IDs
link as `https://arxiv.org/abs/<id>`.

## Cosmoglobe, BeyondPlanck, Commander
| arXiv | Takeaway |
|---|---|
| 2011.05609 | BeyondPlanck I: framework; tau = 0.066 +- 0.013 (2205.10104) with ~30% larger errors. Read in full to Sect. 9.1 (2026-10-27): 4 chains x 1000 samples, burn-in 200, 3200 retained, 1.3 wall h per sample, 670 000 CPU h; no Gelman-Rubin; intensity beta marginal over a (Eq. 118), polarization conditional on a (Eq. 119, Sect. 8.2) |
| 2011.08503 | BeyondPlanck XV (journal XIV): beta_s in 4 regions; high latitude prior-dominated; widening the prior biases tau. Read in full (2026-10-27): Q/U chi^2 over LFI 30/44/70, WMAP Ka/Q/V, HFI 353; conditional on a (Eq. 26), region map smoothed 10 deg; priors N(-3.3, 0.1), Spur and Plane data-constrained (-3.17 +- 0.06, -3.03 +- 0.07), High Latitude and Galactic Centre prior draws; inner Metropolis steps per Gibbs iteration set so the lag-n correlation is below 0.1 (n not stated); 2 chains, about 500 ensemble samples, beta_s Gibbs correlation length 30 to 50 (Sect. 5.3.2); prior shift 0.5 moves the posterior by at most 0.1 |
| 2201.03530 | BeyondPlanck XVI: polarized AME limit depends on the beta_s prior (<0.6% vs 2.5 +- 1.0%). Read in full (2026-10-28, `research/2026-10-28_beyondplanck_xvi.md`): Q/U only, nside 64, 1 deg; LFI 30/44/70 (1000 BP map realizations cycled per Gibbs step), WMAP K/Ka/Q/V, DR4 353 GHz as fixed dust template, f_sky 0.91; beta_s per pixel common to Q and U, Metropolis conditional on amplitudes with proposal width equal to the prior width, 'repeated until convergence' (Sect. 4.3), inner steps, chains, burn-in, acceptance, correlation length and beta_s uncertainties not stated; 500 Gibbs samples per chain (Sect. 5); beta_d drawn from its prior 1.62 +- 0.04, T_d fixed 19.6 K; beta_s prior-dominated (Sect. 5.2); AME limits 1.0 vs 0.6 percent and 2.5 vs 3.5 percent quoted inconsistently between Sect. 5.3, the abstract and Sect. 6; monopoles never mentioned |
| 2201.08188 | BeyondPlanck XIII: amplitude-marginal likelihood (their Eq. 23, with +1/2 ln det (A^T N^-1 A)^-1, i.e. the correct flat-prior sign; the code sign is the opposite, F10). Read in full (2026-10-27): intensity beta_s likelihood disabled (prior only), beta_d prior N(1.56, 0.03), nu_p the only data-constrained intensity spectral parameter, full sky, band monopoles re-estimated per proposal (marginalizing them widens the log-likelihood by 1.5 to 2); 2 chains x 750, burn-in 10, visual diagnostics only; band monopoles mutually correlated and anticorrelated with component monopoles (Sect. 6.2, Table 3) |
| 2201.03509 | BeyondPlanck III: computational infrastructure; 170 CPU-h per Gibbs sample. Read in full (2026-10-27): TOD 96.9 h, amplitudes 23.9 h, spectral indices 40.3 h per sample (Table 2); one chain per 128-core node; CG capped at 100 iterations; no Markov convergence diagnostics; 300 000 CPU h total (BP I says 670 000) |
| 2209.04437 | BeyondPlanck IV: Bayesian (constrained) vs frequentist simulations |
| 2303.08095 | Cosmoglobe DR1 WMAP+LFI; 812 CPU-h per sample; dipole 2.5 sigma above BeyondPlanck |
| 2310.13740 | DR1 polarized synchrotron model; beta_s = -3.07 +- 0.07 from T-T plots |
| 2306.15511 | Cosmological-parameter sampling in Commander; ~2000 CPU-h per independent sample, 10x needed |
| 2408.10952, 2408.11004, 2406.01491 | DR2 (DIRBE): global fit, zodi, CIB monopole |
| 2601.10640, 2601.07822 | DR2 four-component dust model |
| 2604.16901 | Cosmoglobe overview 2026: Commander4 multi-node rewrite, planned HFI/AKARI/FIRAS work |
| 2507.05324 | LiteBIRD with Commander3: ~3000 CPU-h per sample |
| 2606.21334 | Low-frequency intensity model with S-PASS, C-BASS, QUIJOTE; synchrotron at 4.76 GHz ~2x Planck 2015 |

## Spatially varying spectral parameters
| arXiv | Takeaway |
|---|---|
| 1910.14170 | Grumitt+ 2020: NUTS ~10^3x faster than MH; hierarchical per-pixel beta removes r bias; funnel divergences |
| 2510.08534 | Rizzieri+ 2025: multi-patch vs multi-resolution vs adaptive; systematic residuals depend strongly on resolution |
| 2604.08463 | Kabalan+ 2026 (FURAX): GPU K-means patches, ~30% lower r upper limit; selection uses sim truth |
| 2109.11562 | Puglisi+ 2022: spectral clustering of sky regions |
| 2205.13968 | de Belsunce+ 2022: Gaussian-kernel beta prior with fixed 5 deg correlation; strongly prior dominated |
| 1811.00479 | Errard & Stompor 2019: ignoring spatial variability biases r by O(0.01) |
| 0804.2645 | Stompor+ 2009: spectral (profile) likelihood |
| 0709.1058 | Eriksen+ 2008: Jeffreys prior removes >4 sigma prior-volume bias |

## Moments and foreground models
| arXiv | Takeaway |
|---|---|
| 1701.00274 | Chluba+ 2017: moment expansion of averaged SEDs |
| 2111.07742 | Vacher+ 2022: beta-only moments insufficient for LiteBIRD; beta and T needed |
| 2006.08628, 2402.17579 | cMILC / ocMILC moment deprojection |
| 2411.11649 | Vacher+ 2025: 3D dust complexity |
| 2510.18305 | Guillet+ 2025: SRoll2 prefers temperature-only fluctuations |
| 2302.05228 | Fuskeland+ 2023: extending to 600 GHz cuts sigma(r) by 30 to 50% |
| 2203.11445 | Weiland+ 2022: beta_s steepens by 0.25 to 0.4 off the plane |
| 2510.17761 | QUIJOTE XIX: beta_s = -3.05 (plane) vs -3.13 (high latitude) |
| 2404.17567 | CLASS 40 GHz: uniform beta_s rejected, PTE < 0.001 |
| 2511.14572 | QUIJOTE MFI2 forecast: ~10x tighter beta_s |
| 2605.23489 | C-BASS + S-PASS: local structures must be modelled |

## Samplers, Gibbs and linear algebra
| arXiv | Takeaway |
|---|---|
| 2303.18221, 2503.01707 | MCLMC; MAMS (adjusted), up to 7x over NUTS |
| 2510.26691 | Flinch: CMB field-level at nside 512, MCLMC 80 to 100x ESS per gradient vs NUTS; hand-written SHT adjoints |
| 2504.20130 | Field-level LSS: unadjusted MCLMC beyond 10^6 dims; preconditioning matters |
| 2603.18845 | nutpie mass-matrix adaptation via Fisher divergence; low-rank+diag ~4x |
| 1112.1392, 1803.03344 | pCN dimension-robust; non-centred hierarchical parameterizations |
| 0807.0624, 1512.06619 | Jewell / Racine rescaling moves for C_l and sky |
| 2111.07664 | Ducrocq+ 2021: auxiliary-variable / overrelaxed Gibbs, 1 to 2 orders of magnitude on C_l |
| 2112.09354, 2209.10512 | MUSE approximate marginal for global parameters |
| 2010.09649, 2301.07825 | Hutch++, XTrace trace estimation |
| 1504.02661 | Stochastic log-det for implicit operators (IFT) |
| 1308.5299, 1710.00621, 1408.3048 | Multilevel and pseudo-inverse CMB preconditioners, 2 to 5x |
| 1803.03462 | Messenger fields are preconditioned fixed-point iterations; CG comparable or better |
| 2001.10347 | Krylov recycling and deflation for sequences of systems |
| Lindgren, Rue, Lindstrom 2011 (JRSS-B) | SPDE approach: Matern fields as sparse GMRFs, including on the sphere |
| 2111.01084 | SPDE 10-year review |
| 2308.01295, 2402.16683 | GP 3D dust map at nside 256 (NIFTy.re) |

## Differentiable tools
| arXiv / link | Takeaway |
|---|---|
| 2603.19600 | FURAX: JAX operators for map-making and component separation |
| 2311.14670 | s2fft: differentiable GPU SHTs (exact only on sampling-theorem grids) |
| 2510.01785 | cuHPX: GPU HEALPix SHTs with analytic gradients |
| jax-healpy.readthedocs.io | JAX HEALPix utilities |

## Validation
| arXiv | Takeaway |
|---|---|
| 1804.06788 | Simulation-based calibration (Talts+ 2018) |
| 2211.02383 | SBC with data-dependent test quantities (Modrak+ 2022) |
| 2302.03026 | TARP coverage in high dimension |

## tau and instrument systematics
| arXiv | Takeaway |
|---|---|
| 2504.16932 | DESI+CMB gives tau = 0.090 +- 0.012, 3 to 5 sigma from low-l |
| 2309.10034, 2504.13254 | PR4 low-l tau ~0.058 |
| 2603.22454 | elica: tau = 0.0575 (+0.0048 / -0.0058) |
| 2501.11904 | CLASS x Planck: tau = 0.053 (+0.018 / -0.019) |
| 2007.04997 | NPIPE: ADC nonlinearity is the primary large-scale HFI polarization systematic |
| 1901.11386 | SRoll2: ADC nonlinearity and polarization efficiency fitting |
| 2011.08082, 2011.06650, 2201.03417, 2201.03478 | BeyondPlanck gain, noise, bandpass, sidelobes |
| 2312.09001, 2506.22217 | LiteBIRD sidelobe and bandpass requirements |

Full-read report with section references for the four BeyondPlanck papers above: `research/2026-10-27_beyondplanck_chain_lengths_and_index_sampling.md` (BP I Sects. 9.2 to 11 not reached: no accessible full text; BP XVI not yet read).
