# Response to the third referee report on Paper A1 (draft 4)

Report: `2026-11-11_referee_report_paperA1_draft4.md`. Verdict: no fatal objection; 12 major and 19 minor items.
I checked the referee's numbers against the files before acting (the chi-square split, the scan curvature and
slope, the region counts, the one-element toy); they hold. Items that needed new computation were pre-registered
as H12 to H16 in the E044 record before any of the runs existed.

State: written 2026-11-11 13:40 with the text changes made; the rows marked "pending" are filled when the runs
are read (the paper shows them in red until then).

## Major

| # | Item | What was done |
|---|---|---|
| 1 | "No region is misfitted" needs its other half: the Plane has 993 for 862 along the CMB spectrum | Confirmed (0.970, 0.927, 1.012, 1.153; the Plane +3.2 sd). New Table `tab:chi` (generated) gives the split by region and direction with the K, Ka row; a paragraph "Goodness of fit" replaces the old sentences; abstract, Sect. 6.3 and conclusions now say that the totals show only that no region is above its expectation and that the model fails in the Plane along the CMB spectrum. The sampling error of the split is simulated (H14, `e044_chi2_split_null.py`, 40 realizations at nside 8 with both fits): pending. I went one step further, as the last report taught: asked where on the sky the excess sits. From the saved residuals with diagonal weights, twelve of the 435 Plane pixels carry about 142 of the excess of 151, and the largest hold W51 (rank 24 of 3072 pixels in the DR1 free-free model), Cygnus X (ranks 6 and 38) and W43 (rank 3). Reading: leakage of intensity into polarization at the brightest thermal sources, which a mask on polarized amplitude does not remove. Test pre-registered as H16 (the brightest 2 percent in free-free masked as well; LCDM and flat fits): pending. |
| 2 | The chi-square reading: percent out of the abstract or conditional; K, Ka value; significance of 1.008 against 0.956; literature that the note had and the paper did not; "without comment" wrong for Lattanzi et al. | The abstract gives the significance only. The body gives the K, Ka value (0.962 +- 0.019), the 2.3 sd between the two directions, says the test was not planned, cites Jarosik et al. 2007 (noise of the degraded maps not precisely described; r4 noise 1.5 to 5.7 percent above the r9 value, the other sign), Hinshaw et al. 2009 (added to the references, record checked), Lattanzi et al. and Planck 2018 V with "call them consistent with unity", and footnote 5 of Svalheim et al. 2023; says that a model-free test exists and was not made. |
| 3 | The curvature at the profile optimum is positive; -366 belongs to the bound | Confirmed from the scan (+109 to +123; J_c = 292). The paragraph is rewritten with the slope of the volume term (+399 per unit index, a pivot of 20.1 GHz) and the near-degeneracy of the synchrotron and CMB columns at an index of -0.08, which I reproduced with one element (slope 0.114 per amplitude, gain 0.28, pivot 20.5 GHz, maximum 3.8 per amplitude above -3.04). Dunkley et al. 2009 are cited for the mechanism, from the passage in the verification file. The generator no longer writes a Hessian or a first-order quantity for a fit with an optimum on a bound. Record R5, H11 and R7 corrected with a dated note. |
| 4 | Pivot "16.9 to 40.2 GHz" | Replaced by the scan value (20.1 GHz at High Latitude); `flat.pivot.*` is no longer defined for a fit on the bound. |
| 5 | Scale and band dependence: sensitivity of the test; clear in the Plane only; mask and regions changed together; one pattern, not two | The text gives two expected sd at l <= 20 (0.036, 0.059, 0.069, 0.063), says the Plane's change is gone and that the Spur and Intermediate changes are as large as before inside a larger noise, and treats the K, Ka arm and the flat-CMB arm as one observation at two sd. The present regions without the mask, all scales and l <= 20 (H12): pending. |
| 6 | "Moves by up to 0.11 with the prior on the CMB" beside "not significant" | Chosen: it is not a sensitivity. The abstract and the conclusions no longer list it as a move; it appears once, as part of the one pattern of item 5. |
| 7 | Zero-mean sentences worded for the map configuration | "In simulations at half the resolution of the map fit" in Sect. 5, Sect. 6.3 and the conclusions; "and more on the maps" replaced by the difference between the two estimators; Sect. 6.3 says what H9 decided. H9 was read at 14:06: the two realizations at nside 16 give -4.1, -2.2, -1.8, -0.7 sandwich errors, as large as the difference on the maps (-3.6, -1.6, -1.5, -0.2). The difference is the estimator's bias; the abstract now says "biased by up to four standard deviations in the configuration below", and the 3.6 sd is no longer listed among the things not understood. |
| 8 | The mask as a model choice | `--bright-fraction` added to the fit driver; thresholds of 3 and 10 percent (H13) and no mask (H12) are rows of Table `tab:moves`: pending. |
| 9 | "Test the sensitivity to the modes they flag": no such test | Clause deleted. |
| 10 | Stale reason for not quoting an index; caption | Sect. 6 opening, the closing paragraph of Sect. 6.4 and the conclusions give the present reasons (one mask and one set of regions; the Plane along the CMB spectrum; the 3.6 sd between two estimators; the chi-square). Caption: "The index of the one-power-law model fitted to ...". |
| 11 | Six realizations at nside 16 cannot place the map values | Forty more at nside 16 with block errors (H15): pending. Until then the text says "at half the resolution". The six that exist: medians 1.47, 1.57, 1.15, 0.97; none of six above the map value at Intermediate latitude (1.44) and in the Plane (1.28). |
| 12 | Name what was found after the fact; report the failure of H4; correct the earlier response | The paragraph "What was planned and what was not" names the failed expectations that bear on a conclusion (H4, the accuracy clauses of H1, the Spur's shift, the sign in H7) and what was found after the fact (the chi-square and its split, the CMB-prior sequence, the zero-mean disagreement, the dust index). The response to the second report said the paper stated the failure of H4; it did not until now. |

## Minor

| # | Item | What was done |
|---|---|---|
| 13 | "Exact first-order" | "First-order" everywhere; the appendix names both dropped terms and their order. |
| 14 | "Too small by 9 to 22 percent" | "The true error is 9 to 22 percent larger than the quoted one" (the quoted one is 8 to 18 percent too small; both are generated). |
| 15 | Number corrections | 0.23 to 0.60 (three interior regions); the floor from the generated value (0.002); 30 pixels for the Spur at nside 8; J_c = 292; "in three of them"; probabilities marked two-sided; "scatter of 0.009". |
| 16 | 1500 fits against a committed test of 500 | "In 1500 simulated fits at the configuration of its unit test"; bias printed as +0.020. |
| 17 | The nested yardstick is approximate | Stated in the text and in the paragraph on the white prior (relative order T/J; the LCDM covariance as the truth). |
| 18 | "Outside the North Polar Spur"; the Plane keeps 435 of 576; the threshold passes through the Spur | In the abstract and in "Regions and mask" (median kept amplitude in the Spur 22 uK against 28.6, verified). |
| 19 | What the BeyondPlanck mask has that ours lacks | Stated: ours has the first step, on another model, and a disc in place of the steps that rest on their chi-square. (Their description was read in HTML only.) |
| 20 | Spur block ratio quoted where none is tabulated | The sentence now gives three ratios and says why the Spur's is not tabulated. |
| 21 | Fig. `fig:blocks` drew 4 of 38 High Latitude patches | Redrawn: every patch, as its score in units of its model error against its share of the information. It shows one High Latitude patch at -4.7 with 5.5 percent of the information (l 276 to 352, b +45 to +87). |
| 22 | "At the edge of what its calculation allows" | "1.7 to 1.95 of the standard deviation expected for it". |
| 23 | Label the dust-index error | Done (Delta chi-square of one, three values, fixed temperature, power-law synchrotron). |
| 24 | "Together" | Removed. |
| 25 | BeyondPlanck XIV and XV | In the reference entry. |
| 26 | Abstract of 650 words | Rewritten at about 380 words; the history of the earlier regions is out of it. |
| 27, 28 | Undefined terms and a sentence that does not parse | Gone with the rewrite. |
| 29 | Generator | `lp*.chi.ratio`, `flat.pivot.*`, `flat.ratioH/J`, the Hessian keys and the single-set null keys removed; the nside-16 null pools both sets. The expectation per element and band uses the observed penalty (55); stated in a comment. |
| 30 | Record slips | Corrected with dates (128 bright pixels in the Plane; 30 pixels; 8.2; R5 and H11). |

## From "what a hostile expert would say"

- 7.2, the mask: threshold varied (H13), no mask (H12), thermal sources (H16); the text says what the threshold
  selects on.
- 7.3, the chi-square: the paragraph now reads as the referee's one sentence asks.
- 7.5: one clause on what the window does above multipole 47.
- 7.6: the promised test of the imbalance modes is deleted; leakage is now a tested hypothesis (H16), not only a
  word in a list.
- 7.7, K band: one sentence in Sect. 6.4 says that everything leans on K and that no arm removes it.
- 7.8: the runaway with a free CMB amplitude is tied to Dunkley et al. 2009.
- 7.9: "an independent referee" is gone; Sect. 6.1 says "an internal review" once.

## Not done

- A model-free test of the noise covariances (year or DA differences with the single-year matrices): the files
  are not on disk; Paper B (Q28).
- Variation of the regions themselves: stated in the paper as not done.
- The dropped terms of the mismatch formula on the maps (estimated by the referee at about a point of the 22).
