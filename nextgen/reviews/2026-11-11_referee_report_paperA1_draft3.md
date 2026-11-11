# Referee report on Paper A1, draft 3

- **Written:** 2026-11-11. Reading started between 08:04 and 08:07 CET (my first listing shows `sec_impl.tex`
  before its 08:07 edit); report written at 08:39 CET (`date`).
- **Referee:** independent agent, hostile brief. I did not see how the paper was written. I did not trust the
  brief, the response to the first referee or the stats review: what I use from them I checked.
- **Draft reviewed:** `nextgen/papers/paperA1/*.tex` at commit fe6fd74 (08:16). The paper changed twice while I
  read it (`sec_impl.tex` at 08:07 and 08:16; nothing else). Line numbers refer to that state. md5 prefixes:
  abstract 310009f0, sec_data 9200d9d3, sec_concl c1a590a3, sec_impl 6b223fad, gen_moves_text 734b6a0b,
  gen_tab_e040 36a63c08, gen_tab_cmbprior e44cf634, gen_tab_moves 4c857c5a. `main.pdf` has 28 pages.
- **Evidence read:** records E040 and E043 in full; E039 (criteria, addendum B, outputs), E022 (setup and
  result table), E041 (ladder output), E042 (result files); the first referee report, the response to it, the stats review of the three additions; the
  three citation files and the two prior-art notes of 2026-11-09 where they bear on a reference; FINDINGS F52
  to F54; `lowres/block.py`, `lowres/window.py`, `lowres/data.py` (regions), `paperA1_make_tables.py`,
  `e039d_cmbprior_sim.py`, `e040_wmap_region_indices.py` (prediction), `e043_compare.py`; the `.npz` and text
  outputs named below; `figures/fig_blocks.pdf`.
- **Own checks:** `nextgen/reviews/referee_scripts_2026-11-11b/` (ten scripts with their `.out`, and the helper
  `list_npz.py`). They read saved `.npz` and text files; one toy simulates block scores directly (40 blocks,
  two parameters, numpy Cholesky of 2 by 2 matrices only). No fit was run, nothing at nside 8 or 16 was
  recomputed, no covariance matrix was loaded, no scipy Cholesky call was made. A few seconds each, under
  200 MB.
- **Web:** arXiv abstract pages and one arXiv HTML page, for reference records and for one passage of
  Svalheim et al. No personal or project detail in any query.

Convention: "the paper says" quotes the draft; "the record says" quotes an experiment record; "(my check)"
is a calculation of mine with the script named; "(from memory)" is unverified.

## 1. Verdict in ten lines

1. **Fatal objection: yes, to one claim, not to the paper.** The claim: "The scatter of the score over patches
   of sky shows that the high-latitude region and the Galactic plane are not one index each"
   (`abstract.tex:31-32`, `sec_data.tex:290`, `sec_concl.tex:26-28`). For High Latitude the evidence does
   not support it. 65 percent of the block variance behind the ratio 3.5 involves the Galactic Center's
   scores (25 percent alone, 38 percent as a cross term), carried into the High Latitude index by the
   coupling of the two indices; one nside-2 patch (l = 337.5, b = 0) that both regions share gives 80 percent
   of the total. With the region's own scores the ratio is 2.1, and 4 of the 60 null simulations lie above
   that (my check).
2. For the Plane the scatter is real, but 58 percent of the block variance is one patch, the one that
   contains Tau A (score 25 model standard deviations from zero; residuals of 20 to 30 uK in Ka, Q and V in
   one pixel). Deleting that patch moves the Plane index by -0.074 to first order, from -2.96 to -3.04.
   Svalheim et al. (2023), whose regions the paper says it uses, state that they mask Tau A and the
   brightest 5 percent of the plane for the index. The paper applies no mask.
3. The region called High Latitude is mostly the inner Galactic plane: 62 percent of its information is in
   patches centred on b = 0, 55 percent of its squared polarized amplitude within 10 degrees of the plane.
   The regions are a surrogate built in E022 from the released index map, with a Galactic Center cut at 25
   degrees chosen there. They are not "the four of Svalheim et al." and the paper shows no map of them.
4. The methods results of Sects. 3 to 5 stand. The corrections the first referee asked for were made and
   nothing there was changed for the worse.
5. The estimator comparison stands as a measurement (-3.08, -3.09, -2.92 on identical maps). Two repairs:
   12.6 is the distance of the 30 GHz optimum from the profile optimum; the change of convention is 0.167,
   which is 13.4 model standard deviations (3.8 block errors). And "as the identity predicts to 9 percent"
   is a Newton step with the Hessian of the fitted marginal on an identity that is exact.
6. The white-prior sentence (`abstract.tex:29-30`) is not worded for what was shown: nside 8, a prior with
   1.9 times the simulated variance (not "the right total variance"), three of four regions, a Monte Carlo
   error of 10 points, the matched prior itself up to 14 percent off on the same skies. A paired calculation
   from the same file gives 15 to 23 percent. The stats review found the opposite sign in another setup.
7. The memory-fault account is supported after the author's own corrections of 08:07 and 08:16. Left over:
   the nside-16 check is of the profile fit only (the flat-prior rows of the central table were computed
   on the old path and not rerun), "every one of them" is not true of the 48-region stage, and "of this
   size" in the abstract reads as if the real-data fits had been exposed.
8. Numbers: 26 mismatches or misattributions, listed in Sect. 3. None changes a conclusion of Sects. 3 to 5.
   Five are in sentences of the abstract or the conclusions (rows 1, 2, 5, 13, 14).
9. The pending items cannot settle what they are attached to as designed: 12 realizations bound a p-value
   at 1 in 13, the null has the model true in every region (so it cannot separate items 1 and 2 above from
   an index spread), and 6 paired differences fix a scatter to 30 percent.
10. Strongest points: the tempering law; the reference-frequency identity shown on identical maps; the
    plug-in trap with its coverage ladder; record keeping that let me find all of the above from saved
    files without running a fit. Weakest: Sect. 6.3 from "Errors that do not assume the model" onward and
    Sect. 6.4, which interpret misfit without having looked at where on the sky it sits.

## 2. Claims of the abstract and the conclusions

"Qualified" means supported with a qualification the text lacks.

| # | Claim (where) | Evidence | Verdict |
|---|---|---|---|
| 1 | Re-solve move: pseudo-marginal alone; inside a Gibbs sweep a tempered profile posterior with T = 1/2 + 1/(4 c_k); width 0.79 k^(1/4); confirmed within 0.4 percent (`abstract.tex:5-12`) | Sect. 3, Table 1; re-derived by the first referee; text unchanged in substance since | Supported. Table 1 caption has a slip (Sect. 3, row 12) |
| 2 | Moving the reference multiplies the index prior by exp(n beta ln(nu'/nu)); shift ln(nu_ref/nu_piv) sigma_1^2 (`abstract.tex:12-16`) | Sect. 4; identity to nine digits in a unit test | Supported |
| 3 | The shift "has opposite signs for references of 23 and 30 GHz with the WMAP and Planck LFI bands" (`abstract.tex:16-18`, `sec_intro.tex:30-32`) | Table `tab:flat` block 1 and block 3 only. With a free CMB both conventions pull flat: +0.006 and +0.167, +0.058 and +0.859 (Table `tab:flat`, block 2); +0.52 and +1.32 on the maps (`s4_bpxv_n16_wm_cmbflat.npz`) | Qualified: true when the CMB is constrained or absent. With a free CMB, the configuration of Sect. 6.6, the sign is the same |
| 4 | Profile score unbiased; curvature too sharp by sqrt(1 + T/J); plug-in sandwich biased the wrong way; debiased form nominal over four decades (`abstract.tex:18-21`) | E041: T/J from 4e-5 to 1.56 in the measured cells (`E041/ladder.txt`) | Supported. The classical errors-in-variables origin is cited through Engel (1985) only |
| 5 | Zero-mean prior: same estimator for one component and uniform noise; otherwise biased by up to 0.6 sd (`abstract.tex:21-23`) | E041 (-0.58 to +0.17) | Supported on the ladder. On the maps it is 0.59 to 0.79 sd from the profile in all four regions (my check), outside the 0.5 written as R3; the paper gives "within 0.8" and not that R3 failed |
| 6 | "the high-latitude polarized synchrotron index is -3.08 ... -3.09 ... -2.92" (`abstract.tex:24-26`) | `s4_bpxv_n16_wm_lcdm.npz`: -3.0825, -3.0923, -2.9249 | Numbers supported. "High-latitude" is not: see 4.5 |
| 7 | "12.6 statistical standard deviations from a change of convention" (`abstract.tex:27`); "by 0.17 ..., 12.6 times its statistical error" (`sec_concl.tex:14-15`) | (flat30 - profile) / sd = 0.1576 / 0.01248 = 12.63; (flat30 - flat23) / sd = 0.1674 / 0.01248 = 13.41 | Qualified. The change of convention is 13.4; 12.6 is measured from the profile optimum. Both are in units of the model error, which the same abstract says is 3.5 times too small here: 3.8 and 3.6 block errors |
| 8 | "as the identity predicts to 9 percent" (`abstract.tex:27-28`) | measured 0.1674 against H^-1 n ln(30/23) = 0.1844: 9.2 percent; H is the Hessian of the fitted marginal at its 23 GHz optimum (`e040_wmap_region_indices.py:141-143`) | Qualified. The identity is exact; the 9 percent is the error of a first-order step that uses the curvature of the function whose two optima are compared. The prediction written before the data (J_c^-1) is off by 11 percent here and 18 percent with a free CMB |
| 9 | "With a free CMB amplitude the same change of convention moves it by 0.8, again as predicted" (`abstract.tex:28`) | 0.7993 against 0.8740 (8.5 percent) | Qualified as row 8. Not said in the abstract: both conventions are then far from the profile optimum (+0.52 and +1.32) |
| 10 | "In every case that we could check ... the first-order identity predicts the move to 10 percent" (`sec_concl.tex:17-18`) | 12 ratios from 0.898 to 0.999, all below one | Supported to 10.2 percent (Spur, flat CMB prior). Systematically high, which is second order and worth saying |
| 11 | "A white prior on the CMB with the right total variance leaves the index unbiased on simulated LCDM skies and makes its error bar too small by a quarter to a third" (`abstract.tex:29-30`, `sec_concl.tex:20-21`, `sec_data.tex:200-202`) | `E040/cmbprior_sim_n8.txt`, 60 realizations at nside 8 | **Not supported as worded.** See 4.6 |
| 12 | "The scatter of the score over patches of sky shows that the high-latitude region and the Galactic plane are not one index each" (`abstract.tex:31-32`) | block over sandwich 3.51 and 7.52 against a null at nside 8 with largest values 2.82 and 2.82 | **High Latitude: not supported. Plane: qualified** (one patch with Tau A carries 58 percent). See 4.2 to 4.4 |
| 13 | "and gives errors of 0.04 and 0.10 where the model gives 0.012 and 0.013" (`abstract.tex:32-33`) | `block_v3_maps_lcdm.npz`: 0.0437, 0.0999; 0.0125, 0.0133 | Numbers supported. 0.044 is the robust error of the High Latitude index in the joint four-region fit and is mostly the Galactic Center's misfit; the region's own scores give 0.026. 0.100 is 0.063 without the Tau A patch's term. In the paper's own simulations this form is too large at High Latitude (rms of z 0.82, 0.59, 0.90) |
| 14 | "the index also moves by 0.02 to 0.11 with the angular scales that inform it" (`abstract.tex:33`), "in every region and by several times what the model allows for nested fits" (`sec_concl.tex:29-30`) | `lowpass*_lcdm.npz`: 0.017 to 0.110; 2.6 to 10.7 nested sd | Supported as a size. The nested yardstick assumes one index per region, which the paper rejects two paragraphs earlier. In the Plane the shift (-0.055 to -0.110) has the sign and size of removing the Tau A patch (-0.074); not tested |
| 15 | "by 0.01 to 0.07 with the prior on the CMB" | 0.006 to 0.073 | Supported |
| 16 | "by up to 0.03 for 0.05 in the fixed dust index" | 0.0274 | Supported |
| 17 | "the official Planck 30 GHz polarization maps do not fit one power law with WMAP K band" (`abstract.tex:35`) | reduced chi-square 2.64, 4.19, 2.32 (PR4, LCDM prior, three noise arms), 5.97 (2018 map, white prior) | Supported. About 70 percent of the excess chi-square is in the Galactic Center and Plane regions, 13 percent of the elements, unmasked (per element and band 8.6 and 1.7 against 0.52 and 0.50 elsewhere, `l_K30_wm.txt`). The sign against Weiland et al. for the 2018 map is still not discussed |
| 18 | "a silent memory fault ... that sat under every low-resolution fit of this size on our machine ... and that no number changed" (`abstract.tex:35-36`) | E043; my comparison of the files | Qualified. See 4.7 |
| 19 | "for a band set without K band it is more than a unit" (`sec_concl.tex:16`) | -4.04 against -1.60 (bound) | Supported, with a flat CMB prior. B2 is pending |
| 20 | "half a dozen patches carry the information of each region, and with those errors the regions do not differ" (`sec_concl.tex:28-29`) | n_eff 5.7, 8.8, 2.5, 5.1 | Supported for three regions. No block error is quoted for the Galactic Center, so "the four regions do not differ" (`sec_data.tex:294`) has no error for one of them |
| 21 | "The spread of published values ... is of the size of these effects and of the convention effect together" (`sec_concl.tex:33-34`) | none | Opinion; fine if left as one |
| 22 | Claims of Sects. 3 to 5 accepted by the first referee | diff of 9ed43ba against HEAD for `sec_model`, `sec_legacy`, `sec_flat`, `sec_profile`, `app` | Corrections 11 to 21 of the first report are all in the text. Nothing worsened. Open from that report: Fuller (1987) for the errors-in-variables variance; a reference for the 2018 30 GHz map |

## 3. Number mismatches

Recomputed with `check_numbers_sec6.py`, `check_block.py` (outputs beside them) and by reading the text files.
Everything in the three generated tables matches its `.npz` to the printed digit. The slips are in the
hand-written prose and in which configuration a number belongs to.

| # | Where | Paper | Source | Correct |
|---|---|---|---|---|
| 1 | `abstract.tex:27` | "12.6 ... from a change of convention" | `s4_bpxv_n16_wm_lcdm.npz` | 13.4 for the change 23 to 30 GHz; 12.6 is flat30 minus profile |
| 2 | `sec_concl.tex:14-15` | "moves ... by 0.17 ..., 12.6 times its statistical error" | same | 0.167 is 13.4 times 0.0125 |
| 3 | `sec_data.tex:205` | Spur with a flat CMB prior "-3.06 +- 0.04" | `s4_bpxv_n16_wm_cmbflat.npz`: -3.0548 +- 0.0358 | -3.05 (rounded twice, through -3.055). STATUS.md has the same slip |
| 4 | `sec_data.tex:196` | flat prior "calibrated at 2.2 to 2.5 times the error" | `cmbprior_sim_n8.txt`: 0.0357, 0.0381, 0.0266, 0.0356 over 0.0159, 0.0181, 0.0108, 0.0161 | 2.1 to 2.5 (Spur 2.10) |
| 5 | `sec_data.tex:218-219`, `sec_concl.tex:17-18` | "within 10 percent in every region and for every prior" | Spur, flat prior: 0.3799 / 0.4230 = 0.898 | 10.2 percent; "to about 10 percent" |
| 6 | `sec_data.tex:215` | pivot "19.7 to 19.9 GHz, to first order" | measured offset at 23 GHz with the predicted slope | with the measured difference of the two optima: 19.4, 19.5, 19.9, 19.6. The first-order pivot is known to 0.3 GHz at shifts of 0.5 to 1.3. "about 19.5" |
| 7 | `sec_data.tex:263-265` | scatter "within 10 percent" of sqrt(spread^2 / n_eff + sandwich^2) | `block_v2_sim_spread01.txt`: 0.0371, 0.0405, 0.0585, 0.0445 against 0.0379, 0.0406, 0.0629, 0.0509 | Plane 13 percent low |
| 8 | `sec_data.tex:366` | removing the 20 K-band modes "changes the K, 30 index by 0.03" | `E042/modes_ctrl.txt` against `s2_K30_wm.npz` | +0.033, -0.009, +0.092, +0.006. 0.03 is High Latitude; the Galactic Center moves by 0.09 |
| 9 | `sec_data.tex:361-363` | residual "16 to 43 times" its noise model "in B as much as in E"; control "2 to 4" | `E042/predict30.txt`, `predictKa.txt` | these are the values after removing a fitted amplitude slope per region. Before: 81 and 40 at l = 2 to 3, 32 and 21 at l = 4 to 6 (E twice B at the lowest multipoles); control 1.5 to 3.3. Say which |
| 10 | `sec_data.tex:88-90, 109` (Table `tab:e039`) | "predicted, coupled", "Eq. (coupled)" | `E039/sim_n8_bpxv.txt` summary: "coupled J^-1 (n ln(30/23))" | the row is the J form, not the H form that Eq. (coupled) states |
| 11 | `sec_data.tex:152-153` (caption of Table `tab:e040`) | "three free amplitudes on five bands leave an expectation of 0.6" | Table `tab:cmbprior`: 0.6, 0.6, 0.4 | three free amplitudes leave 0.4. 0.6 is with the CMB amplitude held by its prior |
| 12 | `sec_legacy.tex:131` (caption of Table 1) | "The grid and the simulation sit 0.1 to 0.2 percent below the law" | Table 1 | the grid is 0.32 percent below at k = 10; the simulation is above the law at k = 1e4 and 1e5 |
| 13 | `abstract.tex:29`, `sec_concl.tex:20`, `sec_data.tex:200` | white prior "with the right total variance" | `cmbprior_sim_n8.txt:1`: simulated CMB rms 0.145 uK; prior 0.2 uK (`e039d_cmbprior_sim.py:73`) | prior variance 1.9 times the simulated one. At nside 16 the ratio is 1.4 in band and 1.0 against the pixel-averaged sky; the simulation is not that case |
| 14 | `abstract.tex:30`, `sec_data.tex:192-193` | error "too small by a quarter to a third", "24 to 32 percent" | sd of z 1.315, 1.241, 0.961, 1.297; matched prior on the same skies 1.082, 0.980, 0.930, 1.141; standard error 0.09 to 0.12 | relative to the matched prior 22, 27, 3, 14 percent. From quoted LCDM error and paired scatter, sqrt(0.0159^2 + 0.0073^2) / 0.0152 and so on: 15, 23, 4, 21 percent (my check) |
| 15 | `sec_data.tex:193-194` | white-prior intervals "cover 48 to 65 percent" | white 0.550, 0.483, 0.717, 0.650; matched 0.600, 0.683, 0.700, 0.617 | in the Plane the white prior covers more than the matched one; binomial error 0.06. Only the Spur is clear |
| 16 | `sec_data.tex:257` | "With independent patches (white CMB prior, nside 8)" | WMAP covariances are dense; stats review item 3 reads the same simulation as showing the coupling at High Latitude | the patches are not independent in that simulation |
| 17 | `sec_data.tex:206-207` | flat minus LCDM "1.8, 0.7 and -0.7" and "3.4" of the simulated scatter | yardstick is the nside-8 paired scatter 0.031, 0.037, 0.021, 0.034 | at nside 16 the paper's own nested rule gives sqrt(sd_flat^2 - sd_lcdm^2) = 0.026, 0.032, 0.018, 0.027: moves of 2.1, 0.8, 4.1, -0.9. The rule reproduces the nside-8 paired scatter to 15 percent (0.032, 0.034, 0.024, 0.032) |
| 18 | `sec_data.tex:196-199` | white minus LCDM "unremarkable", "4.4 standard deviations" | nside-8 paired scatter for an nside-16 move | marked pending by the author; see Sect. 4.8 for what it has to show |
| 19 | `sec_data.tex:45` | "the residual bias is 0.002 or less" | `E038/window_wmap_matched.txt`: +0.0002, +0.0018, +0.0023, +0.0023 | 0.0023; "below 0.0025" as the E040 record has it |
| 20 | `sec_data.tex:409` | "Two published analyses in this configuration" | round-2 citation file, Watts et al. Sect. 2.4: "the introduction of the K-band data is not sufficient to break it" | Watts et al. include K band. One analysis is in this configuration |
| 21 | `sec_data.tex:245` | block error "uncertain by a third or so" | null simulations: sd of ln(ratio) 0.27, 0.21, 0.66, 0.38; review: 2 to 5 degrees of freedom when a spread dominates | a third to a half |
| 22 | `app.tex:46-47` | "We print this sensitivity with every fit that uses such a prior" | `s4_bpxv_n16_wm_lcdm.txt` has no such line; stats review, code defect 5 | not printed for any LCDM-prior fit, that is for the central table |
| 23 | `sec_impl.tex:54-55` | exposed simulations: "every one of them was run again ... and compared" | E043, V5: 29 of 95 realizations of the 48-region stage compared, 9 entries beyond the tolerance set in advance | "every one that this paper uses" |
| 24 | `gen_moves_text.tex:18-19` | low-pass shifts "the same within 0.025 with the white CMB prior" | Spur: 0.008, 0.016, 0.012 (white) against 0.017, 0.031, 0.033 (LCDM) | true as written and hides a factor two in the Spur, where the white-prior shift at L = 20 is 2.1 nested sd |
| 25 | `sec_data.tex:65`, captions of Tables `tab:e040` and `tab:bpset` | "the four regions of Svalheim et al. (2023)" | E022 record lines 91 to 98; `lowres/data.py:161-165` | a surrogate: nearest-plateau assignment on the released index map, Galactic Center defined as prior-only pixels within 25 degrees of the centre, no processing mask |
| 26 | `refs.tex:27` | "Nicholls, G. K., Fox, C., & Watt, A. M." | arXiv:1205.6857 | third author Alexis Muir Watt |

Also, not in the paper: `reviews/2026-11-11_response_to_referee_draft1.md:20` says "done, B2 computed" and names
`E040/s3_bpset_cmb02_wm.txt` and `s3_bpset_lcdm_wm.txt`. Neither existed at 08:33; the second was started at
08:36 and was empty at 08:40. Line 9 of the same file says the fits are running, and the paper marks B2
pending. The table row is wrong.

Hard-coded numbers in the "generated" prose (`paperA1_make_tables.py:214-283`): the control shifts, 0.95 to
1.03, 0.9 to 1.1, 0.025, the white-prior band-subset shifts, 0.005 and 0.010, 0.002 to 0.015, -3.08 and
-2.96, 0.01 to 0.02, 0.04 to 0.10, 22.4 and 22.8. All agree with their sources today; they will not follow
a rerun.

## 4. Methods check

### 4.1 Eq. (block): what it estimates

Re-derived for one index and independent blocks, with h_B = D_B / D and Var(G_B) = D_B. At the optimum the
block score is G_B - h_B sum G, so

- E[sum G_B^2] / D^2 = (1 - sum h^2) / D: the plain cluster sandwich is low by 1 - 1/n_eff.
- E[sum G_B^2 / (1 - h_B)^2] / D^2 = sum h / (1 - h) / D: the paper's form is high by that factor.
- Deleting block B moves the estimate by G_B / (D (1 - h_B)), so the paper's form is the sum of squared
  delete-one moves, G / (G - 1) times the jackknife variance.
- With independent offsets of variance tau^2 the estimate scatters by tau^2 / n_eff + 1 / D about the
  population mean.

`toy_block_leakage.py` confirms all three expectations to three digits (n_eff 6.2: 0.838, 1.200, ratio to the
jackknife 1.0256). They agree with `lowres/block.py:12-15` and with the stats review. The statements of
`sec_data.tex:235-241` are correct and no stronger than Bell and McCaffrey (2002) and Cameron and Miller
(2015) support. Three things the paper does not say:

1. Eq. (block) is written with a scalar h_B and an outer product G_B G_B^T. The code uses one leverage per
   region and block, from the diagonal of the block Jacobian, with negative D_B set to zero (two at High
   Latitude), and the full 4 by 4 inverse of J_c on both sides (`block.py:74-96`).
2. The form is biased high under the model before any coupling. At the leverage of the maps sqrt(sum h / (1 -
   h)) is 1.12, 1.07, 1.31, 1.13 (my check). The paper reports null medians of 1.2 to 1.3 and attributes
   them to the prior ("because that prior couples the patches", `sec_data.tex:287`). The High Latitude median
   was 1.25 with the white prior and 1.28 with the LCDM prior; the two simulations also differ in
   covariances, so the attribution is not isolated.
3. Because the bread is the full J_c^-1, the block variance of one region contains the block scores of the
   others. This is the next point.

### 4.2 High Latitude: the ratio 3.5 is mostly the Galactic Center

`check_block_crossregion.py`, from the saved G, D (`block_v3_maps_lcdm.npz`) and J, T
(`s4_bpxv_n16_wm_lcdm.npz`). The stored values are reproduced exactly (0.0437, 0.0377, 0.0789, 0.0999).
Writing the High Latitude variance as sum over region pairs of J_c^-1[HL, s] M[s, t] J_c^-1[t, HL]:

| Term | Share of the High Latitude block variance |
|---|---|
| High Latitude scores alone | 0.35 |
| Galactic Center scores alone | 0.25 |
| cross term High Latitude with Galactic Center | 0.38 |
| Plane, Spur and their cross terms | 0.02 |

- In units of the model covariance the middle matrix M is 4.1 for High Latitude, 3.5 for the Spur, 72 for the
  Galactic Center and 55 for the Plane. The indices of High Latitude and Galactic Center are coupled with a
  correlation of 0.20 in J_c^-1. 0.20^2 times 72 is 3: the Galactic Center's misfit alone adds three model
  variances to the High Latitude index.
- By patch: block 18 (nside-2 nested pixel 18, centre l = 337.5, b = 0) gives 0.80 of the High Latitude
  block variance. It holds 37 percent of the High Latitude information and 40 percent of the Galactic
  Center's. Its High Latitude score is 1.6 model standard deviations from zero. Its Galactic Center score is
  7.0. The two terms add with the same sign (-0.020 and -0.019 in index units before squaring).
- With the region's own scores only, the High Latitude block error is 0.026 and the ratio to the sandwich
  error 2.08 (CR0: 1.61). With the other regions' scores replaced by their model covariance: 0.026, ratio
  2.06. Against the 60 null simulations at nside 8 with the LCDM prior (`sim_n8_lcdm_blocknull.npz`), 4 lie
  above 2.08 (3 above 1.61 for CR0). Against the white-prior null, 1 of 60. That is a p-value of 0.02 to
  0.07 at the wrong resolution, not "outside the range".
- 59 percent of even the own-score variance is the same block 18, a score of 1.6 sd multiplied by 1 / (1 -
  0.37)^2 = 2.5. Without that term the own-score block error is 0.016 against a sandwich error of 0.012.
- The toy confirms the mechanism: two indices coupled at 0.2, region A true to the model, region B with
  scores over-dispersed by 70 in variance. The block error of A is 2.0 times its model error in the median,
  and even its own-score ratio has a 97.5th percentile of 2.5 (1.7 to 1.8 when B fits).
- The chi-square of the High Latitude region is 0.584 against 0.6. Nothing else in the paper says that this
  region misfits at all scales.

What stands: 0.044 is a legitimate robust error for the High Latitude index of this joint fit. What it
measures is that the index is tied to a region that fails. It is an argument for decoupling or masking the
Galactic Center, not evidence about the High Latitude region. The same holds in part for the Spur: 27 percent
of its block variance is the Plane's scores (24 percent from block 24 alone, the Tau A patch); its
own-score ratio is 1.90, not 2.29. The sentences that have to go: `abstract.tex:31-32` for High Latitude; `sec_data.tex:288-293` ("outside
the range at High Latitude", "So the High Latitude region ... not one index", "a spread of that size is the
scale of what the High Latitude ratio asks for"); `sec_concl.tex:26-28`; the caption of Fig. `fig:blocks`,
"What the block errors ... are made of", which is not true of the light band in the High Latitude panel;
the title of F54. F54 itself ends "Any misfit scatters the block scores, so a spread of the index is one
reading, not a measurement" and the E040 record says the same at line 733. The paper is stronger than its
own record.

### 4.3 Plane: one patch with Tau A

`check_block_positions.py`, `check_tauA.py`, `check_tauA_chi2.py`, `check_block_influence.py`.

- Block 24 (centre l = 180, b = -19.5) carries 4.6 percent of the Plane's information and 58 percent of its
  block variance. Its score is 25 model standard deviations from zero and its one-step index is 1.5 flatter
  than the region's. Fig. `fig:blocks` shows it as an arrow at the top of the Plane panel.
- Tau A (l = 184.56, b = -5.78) is in that block, in nside-16 nested pixel 1591. That pixel has the largest
  fitted polarized amplitude of the 323 Plane pixels (175 uK against a median of 38) and residuals in Q of
  +13.6 (K), -29.7 (Ka), -23.7 (Q band), -19.9 (V) uK where the rms residual of the region is 1.0 to 2.6 uK.
  The neighbour, pixel 1597, has 4 to 13 uK. One pixel in 323 holds 18 percent of the region's unweighted
  residual sum of squares.
- With a diagonal noise proxy the two pixels carry a chi-square of about 440. The Plane's excess over the
  expectation is 549 (0.770 against 0.6 on 3230 element-bands). So "30 percent off in the Plane"
  (`sec_data.tex:223-224`) is mostly Tau A. I did not load the covariances; the proxy is rough.
- Deleting block 24 moves the Plane index by G_B / (D - D_B) = -0.074, 5.6 sandwich errors, to about -3.04
  (first order, not a refit). The gap between the Plane and the other regions is 0.12. The low-pass shifts
  in the Plane are -0.055, -0.084, -0.110. A compact source has a white harmonic spectrum, so an index read
  from l up to L loses its pull as L falls: the sign and the size fit. With a flat CMB prior the Plane fits
  (0.404 against 0.4), as it would if a free achromatic amplitude absorbs a flat-spectrum excess in one
  pixel. These three are inferences, not tests.
- Without the term of block 24 the Plane's block error is 0.063; without the two largest terms 0.039, still
  2.9 times the model error, with five more patches at 2.5 to 5.2 sd (blocks 21, 4, 26, 39 and 22, the outer
  plane and the Fan). So the Plane does scatter beyond the model without Tau A. How far beyond is unknown until
  the fit is redone without it.
- Svalheim et al. (2023), Sect. 2.1.1, as returned by an HTML reading of arXiv v2 (not the PDF): "it is
  important to impose a processing mask during the posterior evaluation"; "we construct a dedicated spectral
  index processing mask for synchrotron emission"; "thresholding on the strongest 5 % of the signal to
  remove the brightest part of the Galactic plane"; "to remove strong sources such as Tau A and their
  surrounding area". Caption of their Fig. 1: "The mask is chosen to reduce temperature-to-polarization
  leakage bias, and to exclude bright point sources." The E022 record, line 97: "No mask beyond non-finite
  pixels (BP XV's processing mask is not on disk; the regions already separate the plane)". The paper says
  "No mask is applied" (`sec_data.tex:66`) and does not say that the source analysis applied one.
- The paper names "compact structure in the plane" as a candidate for the scale dependence
  (`gen_moves_text.tex:27-28`) and then says "We do not know the cause". One look at the per-patch scores
  it already plots would have named the first suspect.

### 4.4 Is the null reference adequate

No, for five reasons, in order of weight.

1. It has the model true in all four regions at once. On the maps the score variance of the Galactic Center
   and of the Plane is 33 and 42 times the model's (plain form; 72 and 55 with the leverage correction). A null that is to license "High Latitude is not one index" must
   carry that failure, or the test must use the region's own scores (4.2).
2. It is at nside 8, the maps at nside 16 (marked pending).
3. The simulated synchrotron sky is the DR1 map smoothed to 5 degrees (`e039d_cmbprior_sim.py`,
   `truth_sync_30`): no compact sources, less concentrated. Largest leverage at High Latitude 0.32 in the
   simulations against 0.37 on the maps; n_eff 7.7 against 5.7. The paper does not say the sky is smoothed.
4. The noise is drawn from the covariances that the fit uses. If the released matrices understate the
   large-scale noise, as the paper quotes Watts et al. for, the block scores scatter with no index spread at
   all, and this null cannot see it.
5. 60 realizations: "0 of 60" is p below 0.016.

Readings of a large block error other than an index spread of the diffuse emission, none excluded by the
paper: misfit of a coupled region (shown here for High Latitude and the Spur); a compact source with its own
spectrum (shown for the Plane); understated large-scale noise; dust that does not follow one modified
blackbody; window residuals at bright compact structure (Gaussian beams to multipole 47); temperature to
polarization leakage in the bright plane; coupling of patches by covariance and prior (the stats review).

### 4.5 The regions

`check_region_map.py` prints a coarse map from the saved region array of the central fit.

- High Latitude: 2194 pixels, 140 of them within 10 degrees of the plane, 131 of these at l = 25 to 68 and
  l = 298 to 335 (the inner Galaxy outside the 25 degree cut), 9 at l = 200 to 245.
  Those 140 hold 55 percent of the region's squared fitted amplitude; pixels above 30 degrees latitude hold
  27 percent. The five brightest High Latitude pixels are at b = 0 to 2.4, l = 323 to 335. By information: 62
  percent in the eight patches centred on b = 0, 14 percent in patches centred at 42 degrees and above.
- The reason is in the E022 record, lines 91 to 96: regions were made "from QU_BETA_MEAN and QU_BETA_STDDEV
  of the BP synchrotron map at nside 16 ... nearest-plateau assignment ... Galactic Center = prior-only
  pixels within 25 degrees of the centre; the rest of the prior-only pixels are High Latitude". The two
  prior-only regions cannot be told apart in the released map, so the boundary is E022's choice. Every
  inner-plane pixel beyond 25 degrees that carries the prior value became High Latitude.
- Spur: 98 of 464 pixels within 10 degrees of the plane, as strips along the edge of the Plane region. The
  released map is smoothed by 10 degrees, so between the Plane plateau (-3.0) and the prior plateau (-3.3) it
  passes through the Spur's value (-3.17), and nearest-plateau assignment labels the edge "Spur". At least
  19 percent of the Spur's information is south of the plane (blocks 16 and 45). `sec_concl.tex:24` calls it
  "the North Polar Spur".
- Consequences. "The high-latitude polarized synchrotron index" (`abstract.tex:25`, `sec_concl.tex:14, 24`)
  is the index of the inner Galactic plane at 25 to 70 degrees from the centre more than of anything else.
  That High Latitude, Spur and Galactic Center all give -3.08 is then no surprise. The comparison of Sect.
  6.6 with what the published analyses report "at high latitude" compares different skies. The differential
  comparison of estimators is not affected: every row uses the same regions.
- I could not see the published figure of the regions. Whether their High Latitude region contains these
  pixels (before their mask) is not something I verified; that the paper's boundary is not theirs is in the
  project's own record.

### 4.6 The LCDM prior and the white-prior statement

- `window.py:51-70` against `sec_data.tex:51-58`: CAMB with H0 67.36, ombh2 0.02237, omch2 0.1200, tau 0.0544,
  As 2.1e-9, ns 0.9649; total EE and BB; Gaussian window of 60 arcmin and the polarization pixel window,
  squared; multipoles 2 to 3 nside - 1 = 47; floor 0.01 uK. The description is accurate. Planck 2018 VI gives
  tau = 0.054 +- 0.007 (arXiv abstract). The floor arm (0.11 uK: indices by 0.0006 or less, errors by 0.3
  percent) matches `block_v3_maps_lcdm_floor011.txt`.
- The statement of `abstract.tex:29-30` rests on `cmbprior_sim_n8.txt`. Problems:
  1. The simulated CMB has an rms of 0.145 uK per element and the white prior 0.2 uK (line 1 of the file and
     `e039d_cmbprior_sim.py:73`). The docstring of the script says "right total variance"; at nside 8 it is
     1.9 times too large. The case the sentence describes was not simulated.
  2. The stats review ran the case the sentence describes, a white prior of the same rms on LCDM skies
     (its Sect. 1, test t1d arm B), and found errors too large, sd of z 0.93, 0.86, 0.93, with the comment
     "the sign depends on the configuration, the miscalibration does not". The paper cites the review for
     other points and not for this.
  3. Rows 14 and 15 of Sect. 3: the excess over the matched prior on the same skies is 14 to 27 percent with
     a standard error near 10 points; the paired calculation gives 15 to 23 percent; the fourth region shows
     nothing; coverage separates the two priors in the Spur only.
  4. "Unbiased": the largest offset is -0.0048 +- 0.0027 in the Spur, -0.28 quoted sd. The paired mean of
     white minus LCDM is negative in all four regions (-0.0009, -0.0020, -0.0002, -0.0010, each about one
     standard error), the direction of the move on the maps. "Unbiased to 0.3 standard deviations at nside
     8" is what was shown. The appendix says such a bias grows as the square root of the number of
     elements, so nside 8 does not carry to nside 16 by itself.
  What can be said: in one nside-8 simulation a white prior of 0.2 uK on LCDM skies left the index unbiased
  to 0.3 sd and understated its error by about a fifth in three of four regions; the sign and size depend on
  the configuration. The exact answer for the map configuration needs no simulation: the covariance of the
  white-prior estimator under the LCDM covariance is J_w^-1 Var_true(score_w) J_w^-1 with matrices the fit
  already forms.

### 4.7 The memory fault

- My comparison (`check_e043_compare.py`): E039 stage 1 against its rerun, 219 of 221 lines identical after
  removing timings; the two others differ in the optimizer's iteration count only (6 against 7, 5 against 6)
  with every number equal. Plug-in addendum 52 of 52 lines. Nside-8 real-data arm 14 of 14. Old path against
  new path at nside 16: 12 of 12. The reruns are separate executions (the iteration counts show it), from the
  same seeds. A corrupted realization would have shown if it changed any printed number by 1e-4 in an
  index, 1/150 of the error. The comparison is adequate for the claim "no printed number changed".
- Exposure by code: the only calls on the faulty path were `cho_factor(lower=True)` in `likelihood.py` (order
  = elements times components) and `sl.cholesky(lower=True)` for the noise draws of the E039 scripts (commit
  e9589fb). Everything else uses `np.linalg` (LU or the sound triangle). The first reviewer's scripts of
  2026-11-10 work at 48 to 192 elements as far as I read them, below the threshold. The account in the E043 exposure
  table is consistent with the code.
- The author replaced "none of several hundred factorizations" by "about a dozen" at 08:07 and added 40 more
  factorizations at orders 12288 and 18432 at 08:16, with a positive control. That is now sound.
- Left over. (i) "its output is reproduced in every printed digit" (`sec_impl.tex:52-53`): the old-path and
  new-path reruns are profile-only fits of the white-prior baseline (12 lines against the original 14).
  The flat-prior and zero-mean rows at nside 16, which carry the 12.6, were not rerun on either path under
  the canary, and the LCDM run of the central table (02:46 to 03:23) was on the old path. Its profile row is
  reproduced on the sound path by `block_v3_maps_lcdm.npz` to six decimals (my check). The other rows rest
  on "no overrun at this order" alone. (ii) Row 23 of Sect. 3. (iii) `abstract.tex:35-36`, "sat under every
  low-resolution fit of this size": no size has been named in the abstract, and the fits the abstract
  reports were not exposed. Say "under every nside-8 simulation of this work".

### 4.8 The pivots, the nested rule, and the pending items

- Pivots of Sect. 6.3. With the LCDM prior, 23 exp(-(flat23 - profile) / (flat30 - flat23) ln(30/23)) =
  23.36, 23.34, 23.34, 23.41, and 23.33 to 23.39 with the predicted slope: "23.3 to 23.4" is right. With the
  flat prior see row 6 of Sect. 3. `sec_flat.tex:175-176` ("about 20 GHz") is the safer wording.
- Nested rule. Var(low - all) = Var(low) - Var(all) needs the all-scales estimator to be efficient under the
  model. E039 addendum B checks it for the large-scale estimator at nside 8, L = 12, white prior: 0.0092,
  0.0097, 0.0057, 0.0058 against 0.0104, 0.0092, 0.0061, 0.0052 (ratios 0.89 to 1.12; "0.9 to 1.1" is right).
  The stats review adds nside 4 with the LCDM prior. The paper also uses the rule for K, Ka against four
  bands (`gen_moves_text.tex:36-38`) and between two cuts (lines 20 to 21); neither use is simulated, and the
  text does not claim it is. The arithmetic is right in every case I recomputed (2.5, 0.9, 5.2, -0.2; 1.3,
  1.8, 3.4, 0.8; 2.8, 2.0, 5.4, 6.7). The objection is to the reading: "times what the model allows" uses
  sandwich errors of a model the paper has just rejected. Against block errors the High Latitude and Spur
  shifts at L = 20 are 0.5 and 0.8.
- What each pending item would have to show.
  1. Block null at nside 16 (`sim_n16_lcdm.txt`, 12 realizations). No realization had finished when the
     output was moved to `E040/failed/sim_n16_lcdm_first_attempt_too_much_memory.txt` at 08:36; a staged
     version (`e040j_sim_n16_stages.py`) was being tried at 08:40. For the Plane: all 12
     ratios below about 3. For High Latitude it is the wrong test (4.4, point 1). If it is kept, compare the
     own-score ratio (2.08), print the leverage of the simulated skies, and note that 12 values bound p at
     0.08; 4 of 60 at nside 8 exceed 2.08, so about one of 12 is expected to.
  2. Scatter between CMB priors at nside 16 (6 realizations with all priors). The paired sd of white minus
     LCDM would have to be at least about 0.006, 0.008 and 0.003 for "unremarkable" to stand at High
     Latitude, in the Spur and in the Galactic Center, and below 0.010 for the Plane to stay beyond 3 sd.
     Six pairs fix a standard deviation to 30 percent. They cannot show an sd of z of 1.25 against 1.0. The
     exact calculation of 4.6 can.
  3. B2. As written in the record: the flat-prior optimum at a 30 GHz reference within 0.05 of the profile
     optimum in the Spur and the Plane with the 0.2 uK CMB prior. If it holds, Sect. 6.6 stands as a
     statement about a free CMB amplitude. If the shift stays at tenths, `sec_flat.tex:173-175` (pivot of
     29.5 to 30.2 GHz with a constrained CMB) is false for these maps and has to be limited to E022. The
     LCDM-prior version should be reported with it, since that is the paper's prior.

## 5. Pre-registration check

Times from `git log` and file times.

| Statement in the paper | Record | Verdict |
|---|---|---|
| Expectations for the baseline "written down, with numbers, before the first fit to real data" (`sec_data.tex:68-70`) | R1 to R9 committed 2026-11-10 20:35; first real-data fit 22:04 | True for the white-prior baseline. The central table is now the LCDM fit, whose expectations (C1 to C3) were committed at 02:43 on 2026-11-11, after the first referee report and 4.6 hours of real-data results. The paragraph should say so |
| "had fixed in advance that the LCDM covariance would replace it if any index moved by more than one standard deviation" (`sec_data.tex:185-186`) | "If C1 fails ..." committed 02:43; LCDM fit started 02:46 | True, by three minutes, and after the referee had shown the white prior wrong by factors of 200. Give the time |
| "The rule for choosing between this form and two others ... was written before either test was read" (`sec_data.tex:265-267`) | rule committed 04:01 in the entry "Stage 5 results, first part (read 03:52 and 04:00)"; that entry reports the model-true simulation for the uncorrected form (`block_errors_sim_spread0.txt`, 03:59:32) and quotes n_eff from "first realizations of the leverage-corrected script" | Overstated. The rule was written after the first test had been read for the plain form and found low, and after the map ratios had been read. And "in doubt the more conservative" selects the paper's form by construction, since it is never below the other two. The summed distances were 1.35, 1.03, 0.91 with entries known to 0.1 to 0.2 |
| Stage 5 (block errors) pre-registered 03:20 | committed 03:21; first map result 03:52 | True |
| Null test G1 to G3 "written before it is read" | committed 06:36; run 06:33 to 07:10; a look at two realizations is disclosed | True |
| Addendum B (large-scale estimator) pre-registered 03:20 | committed 03:21; result 06:56 | True |
| C4, stage 4b, stage 4c | committed 03:57, 04:01, 06:36; outputs 04:54, 05:37 on, 07:47 on | True |
| "Six criteria were fixed in advance", two failed (`sec_data.tex:81, 126-128`) | V1 to V6 | True |

Failed expectations that the paper does not report as such:

- R7 and C2: a flat CMB prior was to change the profile index by less than 1 sd (R7) or 1.5 sd of the
  flat-prior error (C2). Measured against the central fit: +4.4, +1.5, +8.2, -1.8 baseline sd; the record
  marks C2 failed in two regions. The paper presents the same numbers against a yardstick chosen
  afterwards, under which "three regions are within the scatter" (`sec_data.tex:206-207`).
- R3 with the LCDM prior (zero-mean prior within 0.5 sd): 0.62, 0.79, 0.69, 0.59. Not marked in the record
  either.
- D1 (ratio of block to sandwich error is 1 under the model) and D2 (1.2 to 2 at High Latitude, 1.5 to 3 in
  the Plane on the maps): both failed; the paper reports the null medians and the map ratios without the
  expectations.
- C4c: the author expected one of two readings and found neither. E2 (two regions), F2 (one region), B1
  (two upper bounds), S3 (Plane), R4 (which band), R8, R12 (Spur) are in the records only. That is
  acceptable; the sentence "the records ... record the expectations that failed" is true.

## 6. References

Not covered by the three citation-verification files (ten entries). The two prior-art notes of 2026-11-09
cover them at summary level. My check is of the record only.

| Entry | Checked | Result |
|---|---|---|
| Alquier et al. 2016 | arXiv:1403.5496 | title and four authors as cited; journal line not on the arXiv page (Stat. Comput. 26, 29, from memory) |
| Andersen et al. 2023 | arXiv:2201.08188 | "BeyondPlanck XIII. Intensity foreground sampling, degeneracies, and priors", A&A 675, A13. The statement about marginal and ridge samplers is not in the abstract; it rests on the novelty note |
| Ceperley and Dewing 1999 | web search, arXiv physics/9812035 | J. Chem. Phys. 110, 9812 to 9820; exact for Gaussian noise of known variance, as cited |
| Dickinson et al. 2009 | arXiv:0903.4311 | ApJ 705, 1607. The Jeffreys prior is not in the abstract; novelty note quotes Sect. III.1 |
| Jew and Grumitt 2020 | arXiv:1907.11426 | DOI 10.1093/mnras/staa1233; "objective reference prior" in the abstract. By the novelty note their reference prior is on the amplitudes, so "priors on the indices" (`sec_flat.tex:12-14`) is loose for this one. Their abstract also reports an intrinsic scatter of the index, which is prior art for Sect. 6.3 |
| Lewis et al. 2000 | arXiv:astro-ph/9911177 | ApJ 538, 473 |
| Medina-Aguayo et al. 2016 | arXiv:1503.07066 | title and authors as cited; journal line from memory |
| Nicholls et al. 2012 | arXiv:1205.6857 | title as cited; third author is Muir Watt. "Over-dispersed" is not in the abstract; the prior-art note has it from their Fig. 2 at summary level |
| Planck 2018 VI | arXiv:1807.06209 | A&A 641, A6; tau 0.054 +- 0.007 |
| Planck 2015 X | arXiv:1502.01588 | title as cited; A&A 594, A10 from memory |

Statements against the verification files:

- Liang and Zeger (1986): the file read the abstract only and supports "standardly attributed". The paper
  says "in common use". Fine.
- Bell and McCaffrey (2002), Cameron and Miller (2015): bias of the plain form with few or unequal
  clusters; the (1 - h)^-1 form equal to G / (G - 1) times the jackknife in linear regression; t quantiles.
  All as the file has them. The form without the factor is Mancl and DeRouen (2001), not cited.
- Carter et al. (2017): the paper says "an analogue, in shares of the information". The file's flag 1 (their
  quantity is a share of the variance) is respected.
- Fuskeland et al. (2014): -2.98 +- 0.01 is their weighted mean over all Galactic plane regions, inner
  Galaxy included, with a longitude dependence of 0.17 amplitude. The paper compares it with its outer-plane
  region including Tau A (`gen_moves_text.tex:29-32`) and concludes "that favours the sky". The comparison
  is between different regions; the conclusion does not follow.
- Svalheim et al. (2023): quotation and prior as verified. Stronger than verified: "a synchrotron amplitude
  defined at 30 GHz, and a CMB amplitude without a prior \citep{svalheim2023}" (`sec_data.tex:377-378`). The
  reference frequency is from a figure caption and the CMB prior from the released parameter file (E022),
  not from the paper. Missing: their processing mask (4.3).
- Watts et al. (2023): the three statements match the file. Wrong: row 20 of Sect. 3.
- Weiland et al. (2022): numbers as verified. Not said: they attribute the discrepancy primarily to the 30
  GHz map; and the paper's own fit with the same 2018 map gives K to 30 flatter than K to Ka in its Plane
  region (-2.71 against -2.98), the opposite sign to their -3.27 against -2.90 in the Fan. The first referee
  raised this; the draft still passes over it.

## 7. What a hostile expert would say

1. **Scope.** "A methods paper whose data section is half its length, has three red items and declines to
   quote an index." The paper answers in part (`sec_data.tex:296-297`: the comparison is differential). It
   does not answer why Sects. 6.4 to 6.6 are then here. One sentence will not do it. Tables `tab:e040` and
   `tab:cmbprior` are the demonstration; the rest belongs to the paper that measures the index.
2. **Regions.** "Four regions from one earlier analysis" understates it: they are a reconstruction, the
   Galactic Center boundary is the author's, the High Latitude region is mostly the inner plane, and there
   is no figure. Not answered. One check: a map of the regions with the information share per patch beside
   it. One sentence: how they were made.
3. **Masks, leakage, bandpass mismatch.** Not answered, and no longer hypothetical: the Tau A patch carries
   most of the Plane's excess chi-square and about 0.07 of the 0.12 by which its index differs from the
   other regions (first order). The WMAP team's polarization mask and the processing mask
   of Svalheim et al. exist for this. One check: the central table with Tau A and its neighbours removed
   (rows and columns of N^-1 deleted, or a free amplitude per band there), and with the brightest 5 percent
   of the plane removed.
4. **Nine-year covariances as released.** Answered in part (`sec_data.tex:19-21`, the 20-mode test). The test
   is for K band, with the white prior, on the K, Ka fit, and by the paper's words tests "the modes the
   matrix knows about". The block-error null draws its noise from the same matrices. One sentence in Sect.
   6.3 would say that understated large-scale noise scatters block scores like an index spread does.
5. **Five to nine effective patches.** Answered (`sec_data.tex:241-247`, `sec_concl.tex:46-48`). Not said:
   one patch holds 37 percent at High Latitude and 36 percent in the Plane; with 2 to 5 degrees of freedom
   the 95 percent multiplier is 2.6 to 4.3, so "0.04" is a statement about 68 percent only.
6. **Novelty.** (i) Tempering law: Ball et al. have the leading term for the sign rule, and the paper says
   so. New are the constant 1/2, the exact c_k and the multi-move law. A hostile reader asks whether any
   published number used the move; the paper says it does not know. Answered honestly. (ii) Reference
   frequency: "a one-line Jacobian". The paper claims only the closed-form shift and the pivot, hedged
   ("as far as we can find"). Acceptable. (iii) Sandwich and debiasing: these are the variance and the
   moment correction of the functional errors-in-variables slope in a several-band dress (first referee,
   from memory, still unchecked). Citing Engel (1985), a Dutch methods journal, will not satisfy a
   statistician. One sentence and one reference (Fuller 1987, read first) would.
7. **The headline number.** "12.6 standard deviations in units of an error you say is 3.5 times too small."
   Not answered. Give both units.
8. **"As predicted."** "You differentiate an exact identity and find that Newton's method works to 9
   percent." Not answered. Say what is exact and what is first order, and quote the prediction that used no
   flat-prior fit (J_c^-1: 11 percent, 18 percent with a free CMB).
9. **The memory fault in the abstract.** "A pre-release operating system bug in a methods paper on spectral
   indices." It earns a paragraph in Sect. 8. In the abstract it costs two lines and invites the question
   which numbers were at risk; the true answer is none that the abstract reports.
10. **K band against 30 GHz.** "Your fit with the 2018 map has the opposite sign to Weiland et al. with the
    same map. One of you has an error, and yours has a reduced chi-square of 6 and no mask." Not answered.

## 8. Required changes

**Fatal**

1. `abstract.tex:30-33`, `sec_data.tex:283-294`, `sec_concl.tex:26-28`, caption of Fig. `fig:blocks`, title of
   F54. Withdraw "the high-latitude region ... [is] not one index". State that the block error of 0.044 is
   the error of the High Latitude index in the joint fit, that 65 percent of its variance involves the
   Galactic Center's scores, that the region's own scores give 0.026 (ratio 2.1; 4 of 60 null simulations
   above), and that no misfit of that region is shown. Or establish the claim: fit High Latitude with the
   Galactic Center decoupled or masked and compare own-score ratios with a null that has the same leverage.

**Major**

2. `sec_data.tex:65-66`, captions of Tables `tab:e040`, `tab:cmbprior`, `tab:bpset`, `abstract.tex:25`,
   `sec_concl.tex:14, 24`. Say what the regions are (E022 surrogate, nearest plateau, 25 degree cut, no
   mask), show them, and stop calling the first one "high latitude" while 62 percent of its information is
   in patches centred on the plane. Better: redo the central table on regions that are what they are named.
3. `sec_data.tex:66`, `sec_data.tex:221-225`, `gen_moves_text.tex:23-33` (template in
   `paperA1_make_tables.py`), `sec_concl.tex:29-33`. Name Tau A. Give the central table with it removed, or
   at the least the first-order numbers: Plane index -0.074, block error 0.063, about 440 of the 549 excess
   chi-square. Say that Svalheim et al. mask it. Remove "We do not know the cause" for the Plane until the
   masked fit has been looked at. Remove "that favours the sky" (Fuskeland's regions are not these).
4. `abstract.tex:29-30`, `sec_concl.tex:19-21`, `sec_data.tex:191-202`. Reword to what was simulated: nside 8,
   prior of 0.2 uK on a CMB of 0.145 uK, unbiased to 0.3 sd, error understated by about a fifth in three of
   four regions, sign found to depend on the configuration by the independent review. Give the matched
   prior's sd of z and coverage beside the white prior's. Or replace the simulation by the exact
   covariance of the white-prior estimator at the map configuration.
5. `abstract.tex:27`, `sec_concl.tex:14-15`. 13.4 for the change of convention, or say that 12.6 is measured
   from the profile optimum. Add the size in block errors.
6. `abstract.tex:27-28`, `sec_data.tex:169-171, 217-219`, `sec_concl.tex:17-18`. "The identity" is exact. Say
   that the first-order response, with the curvature of the fitted marginal, gives the move to 5 to 10
   percent and is high in all 12 cases; give the J_c prediction as the one made without a flat-prior fit.
7. `sec_data.tex:249-267`. Remove "With independent patches". Say that the simulated sky is smoothed to 5
   degrees. Say that the form is too large at High Latitude in all three simulations (0.82, 0.59, 0.90).
   Correct the sentence on when the rule was written (Sect. 5). State the other readings of a large block
   error (4.4).
8. The three pending items. Redesign as in 4.8 before the red text is replaced by numbers.
9. `sec_data.tex:204-209`. Use the nside-16 yardstick (2.1, 0.8, 4.1, -0.9) or mark the nside-8 one as such,
   and say that R7 and C2 failed.
10. `sec_data.tex:409-414`. Watts et al. include K band: one analysis is in this configuration. Add that
    Svalheim et al. mask the brightest plane and Tau A for the index. Attribute the 30 GHz reference and the
    free CMB amplitude to the figure caption and the parameter file.
11. `sec_data.tex:341-348`. Say that the fit with the 2018 map has the opposite sign to Weiland et al. in
    the outer plane, and that the misfit of K with 30 GHz sits mostly in the two unmasked bright regions.
12. `abstract.tex:16-18`, `sec_intro.tex:30-32`. Add "when the CMB amplitude is constrained"; with a free
    CMB both conventions pull the same way.

**Minor**

13. `sec_data.tex:205`: -3.05. `sec_data.tex:196`: 2.1 to 2.5. `sec_data.tex:215`: about 19.5 GHz.
    `sec_data.tex:263-265`: within 13 percent. `sec_data.tex:45`: 0.0023. `sec_data.tex:245`: a third to a
    half. `sec_data.tex:193-194`: give the matched coverage.
14. `sec_data.tex:361-366`: say that 16 to 43 is after removing a slope per region; give all four regions for
    the 20-mode test (the Galactic Center moves by 0.09).
15. Table `tab:e039`, last row: label it J^-1, or recompute with H.
16. Caption of Table `tab:e040`: two free amplitudes and a CMB amplitude held by its prior leave 0.6.
17. Eq. (block): one leverage per region, negative D_B set to zero, full J_c^-1; say that the form is high by
    sum h / (1 - h) under the model (1.12, 1.07, 1.31, 1.13 here).
18. `sec_data.tex:68-73, 185-186`: the baseline that was pre-registered before real data is the white-prior
    one; give the time of the replacement rule.
19. `sec_impl.tex:52-55`: the nside-16 rerun is of the profile fit; "every one that this paper uses".
    `abstract.tex:35-36`: "under every nside-8 simulation", or move the sentence out of the abstract.
20. `app.tex:46-47`: the sensitivity is not printed for LCDM-prior fits; print it or drop the sentence.
21. `sec_legacy.tex:131`: 0.2 to 0.3 percent for the grid; the simulation scatters about the law.
22. `gen_moves_text.tex:18-19`: give the Spur's white-prior shifts. Line 41: "where the model does not fit
    in any case" is said of High Latitude too, whose chi-square is 0.58.
23. `sec_concl.tex:24`: "Spur", not "North Polar Spur". `sec_data.tex:294`: three regions, since the fourth
    has no block error.
24. `refs.tex`: Muir Watt; add Mancl and DeRouen (2001), Fuller (1987) once read, Planck 2018 II for the 2018
    30 GHz map. `sec_flat.tex:12-14`: Jew and Grumitt's reference prior is not on the index.
25. `paperA1_make_tables.py`: compute the numbers now typed into the prose template, or mark them.
26. `reviews/2026-11-11_response_to_referee_draft1.md:20`: B2 is not computed.
27. FINDINGS F53 and F54, STATUS of 07:10: carry the same statements as rows 3, 13, 14 of Sect. 3 and item 1.

## 9. What I did not check

- I ran no fit. The own-score ratios, the decomposition by region and the delete-one shifts come from the
  saved G, D, J and T of one fit and are first order. The statement that the Plane index is -3.04 without
  the Tau A patch is an extrapolation, not a result.
- I did not load any covariance matrix. The chi-square attributed to the Tau A pixels uses a diagonal proxy.
- I identified Tau A by position and by the size of the residual. I did not check its flux against the
  fitted amplitude beyond order of magnitude (from memory), nor whether beam differences or leakage, and not
  its spectrum, make the residual.
- I did not see the figure of the regions in Svalheim et al. Their mask is quoted from an HTML reading of
  arXiv v2 by a small model, not from the PDF. Whether their Galactic Center region is larger than 25
  degrees is unverified.
- The null distribution of the own-score ratio was not simulated. I compared the own-score ratio of the maps
  with the stored full-matrix ratios of the null simulations, which agree with own-score ratios under the
  model to a few percent in my toy and in the "model-others" variant, not exactly.
- Sects. 3 to 5 and the appendix: I read the diff against draft 1 and did not re-derive them. I relied on the
  first referee's derivations for everything unchanged.
- E036, E038 (except the window residual), E041 (except the range of T/J) and the first half of E042 were
  not rechecked. The LFI noise model was not examined.
- Figures other than `fig_blocks.pdf` were not opened. The PDF text was not compared with the sources.
- The unit tests were not run. When I finished (08:40) the nside-16 simulation had failed for memory and was
  being restarted in stages, the LCDM-prior fit of the band set without K band had just been started and the
  white-prior one (B2) did not exist. `lowres/block.py` had an uncommitted change (memory handling only, by
  `git diff`); my line numbers are those of the committed file.
- For Dickinson et al., Planck 2015 X, Andersen et al. and Nicholls et al. I verified the record, not the
  sentence the paper attaches to it.
- The paper was edited twice while I read. I checked that only `sec_impl.tex` changed; anything after
  commit fe6fd74 is not covered.
