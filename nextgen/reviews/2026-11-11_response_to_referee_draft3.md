# Response to the second referee report on Paper A1 (draft 3)

Report: `2026-11-11_referee_report_paperA1_draft3.md`. The referee's central findings were checked and hold
(the check of the regions is in the E044 record, with numbers). The response is not a rewording: the real-data
section was redone as a new pre-registered experiment (E044, D035) on regions defined by geometry with a
processing mask, and the sections on Planck 30 GHz and on the band set without K band were removed from this
paper. Every number of the new real-data section is written by `nextgen/scripts/paperA1_make_e044.py` from
result files, as a table or as a named macro in the prose (`gen_numbers.tex`).

**State (2026-11-11 12:10):** all fifteen fits of E044 at nside 16 are read (R1 to R6 of the record), with the nulls
at nside 8 (120 realizations, and a null of the chi-square itself). Still running, and shown in red in the paper
until read: a scan of the flat-prior marginal (H11), six simulations at nside 16 (block-error null, low-pass
estimator, zero-mean prior H9, chi-square), and a third null set of 240 realizations at nside 8 (H10).

## Fatal

| # | Item | What was done |
|---|---|---|
| 1 | "The high-latitude region is not one index" is not supported: 65 percent of that block variance involves the Galactic Center's scores | Withdrawn (FINDINGS, amendment of 08:50; F54 withdrawn for that region). The claim is not in the paper. Block errors are now computed from each region's own scores (`lowres/block.py`, `CR3_own`, with a unit test in which a misfitting neighbour is ignored), and are read against a null simulated with the regions, mask, covariances and prior of the fit. On the new regions with the mask the own-score block errors are 1.63, 1.00, 1.44 and 1.28 times the sandwich errors and all four are inside the null: the paper now says that the block scores give no evidence that any region is more than one index. |

## Major

| # | Item | What was done |
|---|---|---|
| 2 | The regions are a reconstruction and are not what they are named; show them | Confirmed: 140 of the 2194 "High Latitude" pixels lie within 10 degrees of the plane and hold 55 percent of the region's squared amplitude. New regions by geometry (`lowres/regions.py`): |b| > 30, the northern Spur, 10 to 30, |b| < 10 degrees, with a map (Fig. `fig:regions`) and a paragraph that says how they are made and that the Spur outline and the mask recipe come from Svalheim et al. as published, not from their files. The central table, the CMB-prior table and the table of analysis choices are redone on them. |
| 3 | Tau A; a mask; "we do not know the cause"; Fuskeland | A processing mask (brightest 5 percent of the DR1 polarized synchrotron model, and a 7.5 degree disc around Tau A), applied exactly through the Schur complement of each inverse covariance (`lowres/mask.py`, unit-tested). The referee's first-order prediction was tested as pre-registered expectation H6: on the old regions with only the Tau A disc masked, the old Plane index moves by -0.084 (predicted -0.05 to -0.10; the referee's estimate -0.074) and its chi-square falls from 0.770 to 0.596. The paper reports this. The sentence comparing with Fuskeland et al. is gone. |
| 4 | The white-prior statement was worded beyond its simulation | Replaced, as the referee suggested, by the exact first-order bias and covariance of a fit whose CMB covariance is not the truth (`lowres/mismatch.py`, Eq. `eq:mismatch`, appendix), tested against 1500 simulated fits, and evaluated at the configuration of the map fit. Result (R3): the white-prior error is too small by 22, 21, 9 and 10 percent, the bias is below 0.1 of the quoted error, and the white-prior fit quotes a smaller error than the LCDM fit while its true error is 5 to 10 percent larger. The measured differences between the two fits are 1.74, 0.26, 1.87 and 0.05 of their expected standard deviation; the paper says that two are near two with the same sign. The earlier nside-8 simulation is cited as agreeing in kind, with its limits; the review's finding of the opposite sign in another configuration is stated. |
| 5 | 12.6 against 13.4; block-error units | The paper gives the move between the two conventions (flat prior at 30 GHz minus at 23 GHz) and its size in sandwich errors, per region, in the table and the text; no number is measured from the profile optimum under that name. |
| 6 | "The identity" is exact; what was tested is first order | Rewritten: the identity is exact; the table gives two first-order predictions (with the curvature of the fitted marginal, and with J_c, which needs no flat-prior fit), the text gives measured over predicted for both and says that both are always high. |
| 7 | Block-error paragraph: "independent patches", smoothed simulation sky, conservative at High Latitude, when the rule was written, other readings | The paragraph is rewritten. It no longer says the patches are independent; it says the simulated sky is the DR1 model smoothed to 5 degrees; it lists the other readings of a large block error; the sentence about when the selection rule was written is removed (the choice of form is described in the record with its true order of events). |
| 8 | The pending items could not settle what they were attached to | Redesigned: the null is simulated with the regions, mask and prior of the fit and compared with own-score ratios (120 realizations at nside 8 in two independent sets, and a few at nside 16); the scatter between CMB priors is computed exactly instead of from six pairs; B2 is out of this paper. |
| 9 | Flat against LCDM prior judged with an nside-8 yardstick; R7 and C2 failed | The flat-prior fit is compared with the LCDM fit in units of sqrt(sd_flat^2 - sd_lcdm^2) at nside 16, as the referee proposed: +1.95, +0.33, +1.86, +0.02 (R5). The paper reports that the two regions near two are the same two as for the white prior, with the other sign, and offers no cause. The failed expectations are in the E040 record; for E044 the paper states the expectations that failed where it matters (the Spur's reference shift, the block errors of two regions, the zero-mean prior, the sign of the Plane's change at large scales). Correction of 13:40 after the third report: draft 4 did not state the failure of the block-error expectation (H4); draft 5 does, in the paragraph "What was planned and what was not". |
| 10, 11 | Statements about Watts et al., Svalheim et al., Weiland et al. and the 30 GHz fits | The two subsections are removed from this paper (D035). The material and the referee's points go to Paper B with E042. |
| 12 | "Opposite signs for 23 and 30 GHz" needs "when the CMB amplitude is constrained" | Added in the abstract and the introduction. |

## Minor

| # | Item | What was done |
|---|---|---|
| 13 | Number slips in hand-written prose | The prose of the real-data section has no hand-typed result numbers any more. The remaining hand-typed numbers are from records of other experiments and were rechecked: "below 0.0025", the validation table, the unit-test numbers of Eq. `eq:mismatch`. |
| 14 | 30 GHz residual numbers | Section removed. |
| 15 | Table `tab:e039`, last row | Labelled as the J^-1 form, in the caption and the text. |
| 16 | "Three free amplitudes leave 0.6" | Caption corrected: two free amplitudes and a CMB amplitude held by its prior. |
| 17 | Eq. (block) against the code | The equation is now the own-score form that the paper uses; the text says that negative patch informations are set to zero, that the joint form exists and what it does, and that the form is biased high by sum h / (1 - h). |
| 18 | Which baseline was pre-registered before real data | The "Expectations" paragraph is rewritten: the present fits are not the first look at these maps; the LCDM rule was written three minutes before the first LCDM fit, after the first referee report. |
| 19 | Memory-fault wording | Out of the abstract. The implementation section says "every one of them that this paper uses", and that all real-data fits reported were run after the fix under the canary. |
| 20 | "We print this sensitivity with every fit" | Removed; the appendix now derives Eq. `eq:mismatch`. |
| 21 | Caption of Table 1 | Corrected (grid 0.2 to 0.3 percent below; the simulation scatters about the law). |
| 22 | Low-pass shifts of the Spur; "does not fit in any case" | The subsection is rewritten from the E044 arms (R6). With the mask the index from multipoles up to 30, 20 and 12 differs from the all-scale index by at most 1.4 of the standard deviation that the model predicts for the difference; on the old regions it was 2.8 to 10.4. The table gives every such arm in those units. |
| 23 | "North Polar Spur"; "four regions do not differ" | Rewritten with the new regions. |
| 24 | References | Muir Watt corrected. Fuller (1987) cited as the textbook for the errors-in-variables variance, Nakamura (1990) and Akritas and Bershady (1996) for the debiasing (records verified, `research/2026-11-11_citation_verification_errors_in_variables.md`; Fuller's text was not read and no formula is attributed to it). Mancl and DeRouen is not cited because it was not verified; the form is described without the name. "Reference priors on the indices" reworded. |
| 25 | Hard-coded numbers in the generated prose | The prose template is gone; numbers are macros written by the generator. |
| 26 | The response to the first referee said B2 was computed | That file is corrected below this one; B2 is out of this paper. |
| 27 | FINDINGS and STATUS carried the same statements | FINDINGS amended (08:50); STATUS entry of 09:05. |

## From "what a hostile expert would say"

- Scope: Sects. 6.5 and 6.6 are removed; Sect. 6.4 is one table and a short text.
- Masks and leakage: a mask is applied; temperature-to-polarization leakage and bandpass mismatch in the WMAP
  maps are named among the readings of misfit and are not modelled (stated in "What we did not do").
- Covariances as released: stated; the null draws its noise from the same matrices, which the block-error
  paragraph now says.
- Novelty of the sandwich and the debiasing: one sentence places both in the measurement-error literature.
- The memory fault: out of the abstract.

## What the referee could not check and E044 now has as results

- The Plane without Tau A is a fit, not an extrapolation (R2 of E044).
- Own-score block errors against a null of the same geometry are simulated.

## What E044 found that was not asked for

Each is in the record with the time it was read, and in the paper.

- **The zero-mean prior with a free scale is biased at this configuration** (R1 and the nside-8 simulation): 3.6,
  1.6, 1.5 and 0.2 sandwich errors steeper than the profile optimum on the maps; 1.05, 0.54, 0.44, 0.22 in 30
  simulations in which the profile optimum is unbiased. The abstract and the recommendations were changed. A
  simulation at nside 16 (H9) tests whether the simulated bias reaches the size seen on the maps.
- **The chi-square is low** (R4, R5): 0.973 of its expectation with the LCDM prior (2.5 sd), 0.956 with a flat CMB
  prior where the count of degrees of freedom is exact (3.4 sd); a simulated null gives 1.002 +- 0.003. I had
  written "the model fits" on a comparison with 0.6 per element and band and had not computed the significance.
  The deficit is not uniform: the one degree of freedom per element along the CMB spectrum carries 1.008 +-
  0.019. The literature has no statement that these matrices overstate the noise and two uncommented V-band
  values of 0.96 to 0.98 (`research/2026-11-11_wmap_lowres_covariance_normalization.md`).
- **With a free CMB amplitude and flat amplitude priors the High Latitude index has no optimum inside the allowed
  range**, in either convention (R5). The abstract had a template sentence for "the move at high latitude"; there
  is nothing to move. A scan is pre-registered (H11) so that the statement does not rest on an optimizer.
- **The dependence of the index on angular scale went away with the mask** (R6), and with it the reason the
  earlier drafts gave for not quoting an index. Sect. 6.4 is now "What moves the index": the paper gives the
  profile indices as the answer of this model on these regions with their sensitivities (CMB prior up to 0.11,
  bandpasses 0.05, dust index 0.02 in the plane), and says that they are not a recommended value.
- **The dust index that the three fits prefer** is 1.57 +- 0.02 (a parabola through three values of chi2_min; not
  planned), against 1.53 +- 0.02 from Planck.
