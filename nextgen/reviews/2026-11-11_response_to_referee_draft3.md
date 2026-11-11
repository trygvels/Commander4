# Response to the second referee report on Paper A1 (draft 3)

Report: `2026-11-11_referee_report_paperA1_draft3.md`. The referee's central findings were checked and hold
(the check of the regions is in the E044 record, with numbers). The response is not a rewording: the real-data
section was redone as a new pre-registered experiment (E044, D035) on regions defined by geometry with a
processing mask, and the sections on Planck 30 GHz and on the band set without K band were removed from this
paper. Every number of the new real-data section is written by `nextgen/scripts/paperA1_make_e044.py` from
result files, as a table or as a named macro in the prose (`gen_numbers.tex`).

**State at the time of writing (2026-11-11 10:00):** the central fit and the Tau A test are read; the fits with
the white and the flat CMB prior, the arms and the nside-16 null are running, and the paper shows their numbers
in red until they exist. This file is completed when they are read.

## Fatal

| # | Item | What was done |
|---|---|---|
| 1 | "The high-latitude region is not one index" is not supported: 65 percent of that block variance involves the Galactic Center's scores | Withdrawn (FINDINGS, amendment of 08:50; F54 withdrawn for that region). The claim is not in the paper. Block errors are now computed from each region's own scores (`lowres/block.py`, `CR3_own`, with a unit test in which a misfitting neighbour is ignored), and are read against a null simulated with the regions, mask, covariances and prior of the fit. On the new regions with the mask the own-score block errors are 1.63, 1.00, 1.44 and 1.28 times the sandwich errors and all four are inside the null: the paper now says that the block scores give no evidence that any region is more than one index. |

## Major

| # | Item | What was done |
|---|---|---|
| 2 | The regions are a reconstruction and are not what they are named; show them | Confirmed: 140 of the 2194 "High Latitude" pixels lie within 10 degrees of the plane and hold 55 percent of the region's squared amplitude. New regions by geometry (`lowres/regions.py`): |b| > 30, the northern Spur, 10 to 30, |b| < 10 degrees, with a map (Fig. `fig:regions`) and a paragraph that says how they are made and that the Spur outline and the mask recipe come from Svalheim et al. as published, not from their files. The central table, the CMB-prior table and the table of analysis choices are redone on them. |
| 3 | Tau A; a mask; "we do not know the cause"; Fuskeland | A processing mask (brightest 5 percent of the DR1 polarized synchrotron model, and a 7.5 degree disc around Tau A), applied exactly through the Schur complement of each inverse covariance (`lowres/mask.py`, unit-tested). The referee's first-order prediction was tested as pre-registered expectation H6: on the old regions with only the Tau A disc masked, the old Plane index moves by -0.084 (predicted -0.05 to -0.10; the referee's estimate -0.074) and its chi-square falls from 0.770 to 0.596. The paper reports this. The sentence comparing with Fuskeland et al. is gone. |
| 4 | The white-prior statement was worded beyond its simulation | Replaced, as the referee suggested, by the exact first-order bias and covariance of a fit whose CMB covariance is not the truth (`lowres/mismatch.py`, Eq. `eq:mismatch`, appendix), tested against 1500 simulated fits, and evaluated at the configuration of the map fit. The earlier nside-8 simulation is cited as agreeing in kind, with its limits; the review's finding of the opposite sign in another configuration is stated. |
| 5 | 12.6 against 13.4; block-error units | The paper gives the move between the two conventions (flat prior at 30 GHz minus at 23 GHz) and its size in sandwich errors, per region, in the table and the text; no number is measured from the profile optimum under that name. |
| 6 | "The identity" is exact; what was tested is first order | Rewritten: the identity is exact; the table gives two first-order predictions (with the curvature of the fitted marginal, and with J_c, which needs no flat-prior fit), the text gives measured over predicted for both and says that both are always high. |
| 7 | Block-error paragraph: "independent patches", smoothed simulation sky, conservative at High Latitude, when the rule was written, other readings | The paragraph is rewritten. It no longer says the patches are independent; it says the simulated sky is the DR1 model smoothed to 5 degrees; it lists the other readings of a large block error; the sentence about when the selection rule was written is removed (the choice of form is described in the record with its true order of events). |
| 8 | The pending items could not settle what they were attached to | Redesigned: the null is simulated with the regions, mask and prior of the fit and compared with own-score ratios (120 realizations at nside 8 in two independent sets, and a few at nside 16); the scatter between CMB priors is computed exactly instead of from six pairs; B2 is out of this paper. |
| 9 | Flat against LCDM prior judged with an nside-8 yardstick; R7 and C2 failed | The flat-prior fit is compared with the LCDM fit in units of sqrt(sd_flat^2 - sd_lcdm^2) at nside 16, as the referee proposed. The failed expectations are in the E040 record; for E044 the paper states the expectations that failed where it matters (the Spur's reference shift, the block errors of two regions, the zero-mean prior). |
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
| 22 | Low-pass shifts of the Spur; "does not fit in any case" | The subsection is rewritten from the E044 arms. |
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
