# Third referee report on Paper A1 (draft 4)

Date: 2026-11-11, read between 12:40 and 13:10 (clock). Paper at commit 2f5c999 (`nextgen/papers/paperA1/`,
`main.pdf` of 12:37, 29 pages); repository head 9fed07f. Referee: independent agent, adversarial brief. I read the
paper, the E044 record, the result files, the generator and the code; I recomputed the numbers of Sect. 6 from the
`.npz` files; I re-derived the masking, the prior-mismatch formula, the count of degrees of freedom and the
reference-frequency identity and tested them at toy size; I ran no fit at nside 8 or 16. Scripts:
`nextgen/reviews/referee_scripts_2026-11-11c/` (`list_macros.py`, `recompute_core.py`, `recompute_sims.py`,
`scan_curvature.py`, `chi2_by_region_direction.py`, `check_regions.py`, `rederive_small.py`, `volume_term_toy.py`).

Convention: "sd" is a standard deviation; "nested sd" is sqrt(sd_arm^2 - sd_base^2); line numbers are those of the
files at 2f5c999.

## 1. Verdict

1. Fatal objection: no. Every central claim of the abstract has evidence behind it; Sect. 8 has 12 major changes.
2. The differential result reproduces from the result files to every digit: the convention moves the index by 0.07
   to 0.48, 4.8 to 16.4 sd. The derivations I redid (masking, degrees of freedom, identity, mismatch) hold.
3. Nearest to fatal: "no region is misfitted" beside "noise overstated by 1.3 to 2.2 percent". The paper's own two
   fits give the Plane 993 for about 862 along the CMB spectrum (+3.1 sd; +4.2 sd under the noise reading).
4. One explanatory sentence is false: the curvature of the flat-prior marginal "at the profile optimum" is not
   negative. -366 is its value on the bound; at the profile optimum the scan gives about +115.
5. One generated number is nonsense in print: the pivot "16.9 to 40.2 GHz, below the lowest band" (`sec_flat.tex:178`).
6. "The dependence on angular scale went away with the mask" is shown for the Plane only: the test is 2 to 5 times
   less sensitive than before, and regions and mask were changed together.
7. The zero-mean sentences are worded for the map configuration and rest on nside 8; "more on the maps" is not shown.
8. The chi-square reading is an unplanned test with one reading in the abstract; K and Ka alone give 0.962 +- 0.019
   along a CMB-spectrum direction, which the paper omits; "exact first-order" is not exact (the appendix drops a term).
9. Pre-registration: H1 to H11 were committed before their results; H4's failure and four after-the-fact findings
   are not labelled in the paper. The six nside-16 null realizations cannot carry "agree".
10. Numbers: of 122 macros in the text 118 are what the text says; four carry an error of definition; five
    hand-typed statements are wrong, stale or loose (Sect. 3).

## 2. Claims

Verdict column: S supported; Q supported with a qualification the text lacks; N not supported.

| # | Claim (where) | Evidence I checked | Verdict |
|---|---|---|---|
| 1 | Tempering law, T = 1/2 + 1/(4 c_k), width 0.79 k^(1/4), confirmed to Monte Carlo errors of 0.4 percent (`abstract.tex:5-12`) | Unchanged since draft 3 except the caption of Table 1, which is now right: grid 0.15 to 0.32 percent below the law, simulation within 1.3 of its errors (0.32 to 0.36 percent). T(10) = 2.532 from c_10 recomputed. Not re-derived by me | S |
| 2 | Identity for the reference frequency; shift; opposite signs for 23 and 30 GHz "when the CMB amplitude is constrained" (`abstract.tex:12-18`) | Identity re-derived (4.4) and checked to 1e-11 with a dense CMB prior and with a flat one. On the maps with the LCDM prior the 23 GHz optimum is 0.9 sd steeper and the 30 GHz one 15.5 sd flatter than the profile optimum | S |
| 3 | Profile score unbiased, curvature too sharp, plug-in biased the wrong way, debiased form nominal over four decades (`abstract.tex:18-21`) | Unchanged. E041 has T/J from 2e-4 to 1.56 in the measured cells. Not re-derived | S |
| 4 (b) | Zero-mean prior "is biased ... by one standard deviation in simulations of the configuration below" (`abstract.tex:22-23`); "one sd in Sect. 6.3" (`sec_profile.tex:139-141`) | `sim_n8_zm.npz`: paired zero-mean minus profile -1.05, -0.54, -0.44, -0.22 sd, scatter 0.012 to 0.020, 30 realizations at nside 8 of one sky. Profile unbiased in the pooled 360 (+0.05, -0.11, 0.00, -0.05, each +- 0.05); in these 30 alone it is +0.26 +- 0.19, -0.36 +- 0.16, +0.06, 0.00 | Q: nside 8 is not "the configuration below" (a quarter of the amplitudes at twice the signal to noise); one sd is the High Latitude value only |
| 5 | Indices at \|b\| > 30 deg: -3.15 profile, -3.17 flat at 23 GHz, -2.70 flat at 30 GHz (`abstract.tex:26-28`) | `lcdm_all.npz`: -3.1469, -3.1744, -2.6960 | Q: the region is \|b\| > 30 deg outside the Spur (1487 of 1600 pixels; the 113 excluded are the North Polar Spur) |
| 6 (f) | Convention moves the index by 0.07 to 0.48, 4.8 to 16.4 sd; first order "predicts every move, too high by 8 to 33 percent" (`abstract.tex:29-31`, `sec_concl.tex:15-22`) | Recomputed: 0.478, 0.154, 0.204, 0.072; 16.4, 5.6, 10.5, 4.8 sd; prediction over measured minus one 33, 19, 16, 8 percent (curvature of the marginal) and 45, 27, 23, 11 percent (J_c) | S for the numbers. "Predicts" is generous at 33 percent. The pre-registered accuracy clauses of H1 failed in two places and the paper does not say so (Sect. 5) |
| 7 (a) | "With a free CMB amplitude the flat-prior marginal has no optimum at high latitude in either convention: the index runs to the edge of the allowed range, -1.6" (`abstract.tex:31-33`, `sec_concl.tex:16-19`) | `scan_flat.npz`: both marginals increase at every one of 20 steps from -4.6 to -1.6 in both slices; gain from the profile optimum to the edge 463 and 1601 (slice a), 593 and 2081 (slice b). At -1.6 the slope is still +338 per unit index (23 GHz) and the curvature has turned negative | Q: true inside the allowed range, shown on two one-dimensional slices; the range ends where the function is still accelerating upward, so -1.6 is the bound and nothing else. The explanation in the body is wrong (row 8) |
| 8 (a) | "At the profile optimum the curvature of the flat-prior marginal is negative in that direction (-366 against an information J of 316): the volume term outweighs the data" (`sec_data.tex:313-316`; record R5 and H11) | `e044_fit.py:179` evaluates `marginal_hessian` at the flat-prior optimum, which for this fit is on the bound (-1.6). The scan gives the curvature of minus ln L along the High Latitude index at the profile optimum: +109 (grid), +115 (local fit), +98 in the other slice. Profile curvature there +277 to +290 (J_c = 292); volume term about -175. The curvature of the marginal is positive from about -4.2 to about -2.0 and negative only above. 316 is J at the fitted amplitudes, not J_c | N. The volume term's curvature is negative (true) but does not outweigh the data at the profile optimum. What makes the function monotonic is the slope of the volume term, +396 per unit index there (0.133 per amplitude), against a net curvature that never exceeds 138. See 4.6 |
| 9 (a) | Pivot with a flat CMB prior "16.9 to 40.2 GHz, below the lowest band" (`sec_flat.tex:177-178`) | Generator line 188 applies a first-order formula to a fit on the bound with an indefinite Hessian: 27.3, 40.2, 16.9, 17.9 GHz | N. From the scan: 23 exp(-396/2974) = 20.1 GHz at High Latitude (19.8 from the other slice). Draft 3 had "about 20 GHz", which was right |
| 10 (a) | LCDM-prior scan has interior maxima near -3.17 and -2.74 (`sec_data.tex:323-325`) | `scan_lcdm.npz`: grid maxima at -3.10 and -2.80; three-point vertices -3.172 and -2.745; joint optima -3.174 and -2.696. Same slice definition as the flat-CMB scan | S (three-point parabola on a 0.15 grid; "near" is the right word) |
| 11 (b) | "The zero-mean prior lands 3.6 standard deviations from the profile optimum at high latitude" (`abstract.tex:33-34`); "it reaches one standard deviation in simulations and more on the maps" (`sec_profile.tex:116-117`); "in simulations it is the one that is biased" (`sec_concl.tex:25-27`) | Maps: -0.104, -0.043, -0.030, -0.004; -3.56, -1.56, -1.54, -0.24 sd | S as a difference between two estimators. Q for `sec_profile.tex:116-117`: that the difference on the maps is a bias of the zero-mean prior is what H9 is meant to test; if H9 fails, 2.5 sd of it is unexplained and is evidence against the model at High Latitude, for either estimator |
| 12 (c) | White CMB prior "makes the error bar too small by 9 to 22 percent, by an exact first-order calculation" (`abstract.tex:34-36`) | `white.npz`: true over quoted 1.218, 1.214, 1.086, 1.098; I reproduce cov_true = J_c^-1 (J + dK) J_c^-1 and the bias from the stored J, T, dK, eg. Formula re-derived (4.3) | Q: (i) not exact: two terms of first order in the mismatch are dropped (the appendix says so for one); (ii) "22 percent" is true over quoted minus one, so the quoted error is 8 to 18 percent too small; (iii) the truth is the paper's own LCDM covariance, which by `sec_data.tex:80-82` lacks 0.11 uK of the real pixel variance; (iv) dK is evaluated at fitted amplitudes without the debiasing that J gets |
| 13 (c) | Expected sd of white minus LCDM is sqrt(sd_true^2 - sd_lcdm^2); measured moves less the bias are -1.74, +0.26, -1.87, -0.05 of it (`sec_data.tex:280-284`) | Recomputed. Derivation in 4.5: Cov(white, LCDM) = J_L^-1 for the terms linear in the noise, so the rule holds to relative order T/J | Q: with J_c^-1 in place of the sandwich for the LCDM fit the first number is -1.57; state the rule as approximate |
| 14 (d) | "With the mask no region is misfitted" (`abstract.tex:36`, `sec_data.tex:193-194`, `sec_concl.tex:34-35`) | Per element and band 0.577, 0.591, 0.572, 0.614 against 0.598: -2.3, -0.2, -2.2, +0.9 sd | Q, and the largest gap in the paper: split as the paper splits the total (LCDM fit minus flat-CMB fit), the Plane carries 993 for about 862 along the CMB spectrum, ratio 1.153 +- 0.048, +3.1 sd (at least +2.95); the other three are 0.970, 0.927, 1.012. Under the paper's own reading of the deficit (all noise variance 0.956 of nominal) the Plane is +4.2 sd in that direction and +2.6 sd in total. See 4.7 |
| 15 (d) | "errors from the scatter of the score over patches of sky agree with the model's" (`abstract.tex:36-37`) | Ratios 1.63, 1.00, 1.44, 1.28; 360 nulls at nside 8: medians 1.34, 1.28, 1.13, 1.07; 51, 244, 50, 78 above the map values (14, 68, 14, 22 percent; Fisher p = 0.17) | Q: the null is at nside 8; the conclusions add "within what a dozen patches can say", the abstract does not; no power is quoted |
| 16 (d) | "the index does not depend on the angular scales or the band subset that inform it beyond the noise" (`abstract.tex:37-38`); largest of 24 is 2.4 (`sec_concl.tex:36-37`) | 24 values recomputed; largest 2.40; sum of squares 31.7; five exceed 1.85 where 1.5 are expected; the 12 low-pass values sum to 9.8; the K, Ka arm alone sums to 10.1 for 4 (p = 0.04 if independent) | Angular scale: Q (about right, but the test is weak, row 17). Band subset: Q, too strong: +2.4 and +2.0 nested sd in two regions with one sign, and the flat-CMB arm moves the same two kinds of region the same way for the same reason (less weight on Q and V) |
| 17 (d) | "The dependence on angular scale went away with the emission that the mask removes" (`sec_data.tex:460-466`); abstract's before and after (`abstract.tex:41-44`) | At l <= 20: old regions -0.024, -0.031, -0.055, -0.084 with nested sd 0.0085, 0.0094, 0.0066, 0.0081; new regions +0.004, -0.028, -0.047, +0.011 with nested sd 0.018, 0.029, 0.034, 0.032 | N as a causal statement for three of four regions. The Spur and Intermediate changes are as large as before and of the same sign; they are "inside the noise" because the noise of the test grew 2 to 5 times. Only the Plane changed clearly. Regions and mask were changed together; no fit separates them |
| 18 (d) | Units (sd_arm^2 - sd_base^2)^(1/2) for the low-pass, K and Ka, and flat-CMB arms (Table `tab:moves`) | Right yardstick for every one of these arms to relative order T/J (4.5): each uses less of the same data under the same model, and the base is the efficient fit if the CMB is LCDM | S, with the rule stated as approximate and conditional on the LCDM prior being the truth |
| 19 | "it moves by up to 0.11 with the prior on the CMB" (`abstract.tex:38-39`, `sec_concl.tex:38`, `sec_data.tex:486-487`) | +0.106 at High Latitude, which the paper's own Table 7 puts at 1.95 nested sd and Sect. 6.4 calls not significant | Q: the same number is used as a sensitivity in the abstract and as noise in Sect. 6.4. It cannot be both |
| 20 (e) | "The chi-square is 2.5 to 3.4 standard deviations low" (`abstract.tex:39-40`) | 16953.8 for 17420 (0.9732, -2.50 sd); 11101.4 for 11612 (0.9560, -3.35 sd). Count of degrees of freedom in `e044_fit.py:131` re-derived and confirmed by simulation (4.2). Null at nside 8: 1.0019 +- 0.0027 | S. Two readings of one data set, not a range. Two-sided probabilities 0.012 and 0.0008 |
| 21 (e) | "which as a statement about the released covariances is noise overstated by 1.3 to 2.2 percent" (`abstract.tex:40-41`) | 1 - sqrt(ratio) = 1.35 and 2.22 percent | Q: one reading of an unplanned test, in the abstract without "if". The paper's own numbers say the deficit is not uniform over directions (row 22) or regions (row 14); K and Ka alone give 0.962 +- 0.019 along a CMB-spectrum direction, which the paper does not report |
| 22 (e) | Deficit absent along the CMB spectrum: 1.008 +- 0.019 against 0.956 +- 0.013 (`sec_data.tex:300-307`) | Recomputed: 5852.4 / 5808 = 1.0076; difference 0.052 +- 0.023, 2.3 sd (2.0 sd with the 0.11 uK floor, 1.0016). Nested difference is chi-square with n degrees of freedom and independent of the flat chi-square (4.2) | Q: "none" overreads 2.3 sd; by region it is 0.970, 0.927, 1.012, 1.153; in the K, Ka, 353 fit the one degree of freedom per element, also a CMB-spectrum direction, has 0.962 +- 0.019 (-2.05 sd) |
| 23 (e) | Two analyses "print reduced chi-square values of 0.96 to 0.98 for V band, and near or above one for Ka and Q, without comment" (`sec_data.tex:204-206`) | Lattanzi et al. Table 1 read in the PDF: V 0.9627, 0.977, 0.9606, 0.974; Ka 1.026, 1.045, 1.012, 1.029; Q 0.9830, 1.004, 0.9778, 0.995. Planck 2018 V Table 12 (HTML through a small model): 1.004, 1.048, 0.970 | Q: Lattanzi et al. do comment: "The reduced chi^2 values we find are consistent with unity for all bands, in all masks considered" (their p. 4). See Sect. 6 |
| 24 (g) | Dust index 1.57 +- 0.02 from a parabola through three chi-square values (`sec_data.tex:471-476`) | 16977.09, 16953.83, 16949.42; vertex 1.5667, 1 / sqrt(a) = 0.0163 | S; it is labelled as a by-product and as a three-point parabola. Add that the error is the Delta chi-square = 1 width, conditional on 19.6 K and on a synchrotron power law, and that the minimum lies between two of three points with no check of the shape |
| 25 (h) | The indices are given "as the answer of this model on these regions, with what they are sensitive to, and not as a recommended value" (`sec_data.tex:488-490`); "We have not measured the polarized synchrotron index" (`sec_concl.tex:32`) | Read against the abstract, the caption of Table `tab:same`, `sec_data.tex:4-9` | Q: defensible, but the reason given at `sec_data.tex:6-9` is stale and the paper measures in three places. See 7.1 |
| 26 | "The spread of published values for this sky, -2.9 to -3.3, is of the size of these effects and of the convention effect together" (`sec_concl.tex:47-49`) | No published value is traced to a convention or to one of these effects | Q: suggestive, not shown. "Is of the size of" is all that can be said, and "together" claims more |
| 27 (i) | Sects. 3 to 5 against draft 3 | `git diff 56c9ed6 HEAD`: five edits. Four are improvements. One is a regression (row 9). One overreaches (row 11, `sec_profile.tex:116-117`). Table `tab:targets` agrees with Sect. 6 for the simulations | S apart from rows 9 and 11 |

### What each pending item has to show

1. **Block-error null at nside 16, 6 realizations** (`sec_data.tex:408`). It can show only that the null does not
   change with resolution: medians of the ratio within about 0.2 of the nside-8 medians (1.34, 1.28, 1.13, 1.07), and
   24 standardized errors with mean inside +-0.4 and rms inside 0.72 to 1.28. It cannot place the map values: under
   the nside-8 null, "0 of 6 above" has a probability of 0.40 at High Latitude. If the nside-16 ratios sit lower or
   tighter than at nside 8, "agree" is unsupported and a null of 40 or more is needed. State of the file when I
   finished (4 of 6, read 13:05, not a result): High Latitude 1.45, 1.32, 1.20, 1.54 against 1.63 on the maps;
   Spur 0.72, 1.75, 1.63, 1.51 against 1.00; Intermediate 1.14, 1.41, 1.12, 1.16 against 1.44; Plane 0.87, 1.09,
   0.96, 0.97 against 1.28. No simulated value reaches the map value in three regions, and the Plane's four values
   are tighter than the nside-8 null. Six will not settle this.
2. **Low-pass estimator at nside 16, 4 realizations, l <= 20** (`sec_data.tex:448`). Sixteen values of (low-pass minus
   all scales) / nested sd with mean inside +-0.5 and rms inside 0.66 to 1.34, and the same for the low-pass
   estimate against the truth. It covers l <= 20 only; the sentence must not extend it to 30 and 12. If the rms is
   near 0.66 the Intermediate value of -1.36 becomes -2.1.
3. **Zero-mean prior at nside 16, 2 realizations** (`sec_data.tex:241-242`, H9). Both realizations inside -2.5 to
   -4.5 (High Latitude), -1.0 to -2.2, -1.0 to -2.0, -0.1 to -0.5. Then `abstract.tex:22-23` can say "three
   standard deviations at the resolution of the map fit". If they give about -1, the abstract, `sec_profile.tex:116-117`
   and `sec_concl.tex:24-27` must say that the simulated bias explains a third of the difference on the maps, and
   Sect. 6.3 must list the remaining 2.5 sd as a disagreement between two estimators that the model does not
   predict. Either way the simulated sky is the DR1 model smoothed to 5 degrees; say that the bias "is set by the sky".

## 3. Numbers

I recomputed every `\N{}` macro that the text uses (123 keys) and the three generated tables from the `.npz` files.
All table entries agree. Mismatches and errors of definition:

| # | File and line | Printed | Correct | Cause |
|---|---|---|---|---|
| 1 | `sec_flat.tex:177-178` | pivot "16.9 to 40.2 GHz, below the lowest band" | about 20 GHz at High Latitude (20.1 from slice a of `scan_flat.npz`, 19.8 from slice b); undefined by this formula for the other regions | `paperA1_make_e044.py:188` uses 23 / exp((flat23 - beta) / hn) on a fit whose optimum is on the bound. The values are 27.3, 40.2, 16.9, 17.9, and 40.2 is not below 22.4 |
| 2 | `sec_data.tex:313-315` | curvature "at the profile optimum" -366 | about +115 at the profile optimum; -366 is at the flat-prior optimum on the bound | `e044_fit.py:179`; generator key `flat.Hdiag` (line 192) is the Hessian at `flat23` |
| 3 | `sec_data.tex:314-315` | "an information J of 316" | J_c = 292; 316 is J at the fitted amplitudes, which by Sect. 5.2 estimates J + T | generator key `flat.Jdiag` |
| 4 | `sec_data.tex:326-327` | "The other three ... are 0.23 to 1.44 flatter" | 0.23 to 0.60 (0.344, 0.603, 0.229); 1.44 is the region on the bound | `paperA1_make_e044.py:190`: the range is taken over all four regions |
| 5 | `sec_data.tex:82-83` | with the floor at 0.11 uK "the four indices change by 0.0006 or less" | 0.0020 or less: -0.0020, -0.0005, -0.0008, -0.0004 (`floor011.npz`); errors +0.24 to +0.32 percent | hand-typed from E040 (old regions); contradicts Table `tab:moves`, which prints -0.002 |
| 6 | `sec_data.tex:405` | "The Spur, with 66 pixels at this resolution" | 30 pixels, 60 elements at nside 8 (`named_regions(8)`; `sim_n8_chi2null.npz` has nper 778, 60, 366, 240) | 266 / 4; the same slip is in the record, line 474 |
| 7 | `abstract.tex:35`, `sec_data.tex:278`, `sec_concl.tex:28-29` | error "too small by 9 to 22 percent" | the true error is 9 to 22 percent larger than the quoted one; the quoted one is 8 to 18 percent too small | `paperA1_make_e044.py:203`: 100 (true / quoted - 1) |
| 8 | `sec_data.tex:270-272` | "A unit test with 1500 simulated fits"; measured bias "+0.021 +- 0.002" in the second region | the committed test (`lowres/tests/test_block_prior_lowpass.py:279-307`) runs 500 fits with 15 percent tolerances. My 1500 fits at its configuration: scatter 0.1133, 0.0889; bias +0.0164 +- 0.0029, +0.0205 +- 0.0023; predicted sd 0.1128, 0.0907; bias +0.0211, +0.0199 | wording; 0.0205 prints as +0.020 |
| 9 | `sec_data.tex:240-241` | maps show the difference "about three times larger in index" | 3.5, 2.9 and 3.3 times in three regions; in the Plane 0.0035 against 0.004, the same | wording |
| 10 | `sec_data.tex:198-199, 300` | "about one percent", "one in a thousand" | two-sided 0.012 and 0.0008; one-sided 0.006 and 0.0004 | say which |
| 11 | `sec_data.tex:395` | block ratio 1.0 for the Spur | Table `tab:same` prints no block error for the Spur | consistency |
| 12 | `sec_data.tex:147-149` | reduced chi-square "0.998 +- 0.009" (12 realizations) | 0.0091 is the scatter over 12, not the error of the mean (0.0026) | E039 record |

Checked and right, among others: all four generated tables; 5808 of 6144 and 2974, 266, 1698, 870; chi-square
16954 for 17420, 2.5 sd, 1.3 percent; 3.4 sd, 2.2 percent; 0.598 and the four sampling errors; pivot 23.3 with the
LCDM prior (23.26 to 23.33; the scan gives 23.3); 0.75, 0.84, 0.86, 0.93 and 0.69, 0.78, 0.81, 0.90; zero-mean
0.104, 0.043, 0.030, 0.004 and 3.6, 1.6, 1.5, 0.2; mismatch bias, ratios, pair sd, correlation 0.2; flat-CMB sd
ratios 2.12 to 2.36 and moves +1.95, +0.33, +1.86, +0.02; all arm moves and their units; dust 23.3, 4.4, 1.57 +-
0.02, 0.002 to 0.008; old low-pass 0.024 to 0.084 and 2.8 to 10.4; Tau A -2.964 to -3.048, 0.596; the 360 nulls
(medians, 97.5th percentiles, counts above, means and sd of z, coverage); the 40 chi-square nulls; the scan gains
463, 1601, 445 and vertices -3.17, -2.74. Table `tab:e039` and the text of Sect. 6.2 agree with the E039 record
line by line. Sect. 6.1: 1e-9, 2.0 to 2.2 uK, 65 to 164, 0.005 to 0.011, 0.0068 (0.00678 from the files; the E040
record prints 0.0067 from rounded values), 0.0025, 0.0085, 0.0035, 1.12 to 1.14, 0.13 percent, 0.17 uK, factor 200
agree with E038, E039, E040. Old regions: joint block ratios 3.5, 2.3, 8.9, 7.5 ("2.3 to 9") and 1.78, 0.77.

Generator, errors of definition that do not reach the text: `lp12/lp20/lp30/lp20_KKa.chi.ratio` and `.chi.nsig`
use the nominal 17420 for fits whose data were reduced to a few hundred modes (0.088 for l <= 20);
`chi.exp` uses the observed prior penalty (55.5) as its expectation (effect below 0.001); `null8.*` pools one set
only; `flat.pivot.*` and `flat.ratioH/J` are undefined for a fit on the bound and should not be written at all.
Record: line 49 "141 of them within 10 degrees of the plane" is 128 of the 154 bright pixels (141 is the number of
Plane pixels masked in all, 128 bright and 13 of the Tau A disc); R6 "8.1" is 8.2.

## 4. Methods

**4.1 Masking.** Dropping elements m is marginalizing over them; the kept block of N has inverse Ninv_kk - Ninv_km
Ninv_mm^-1 Ninv_mk. `lowres/mask.py:16-31` is that expression. Toy check: equal to the inverse of the kept block to
8e-16, and the chi-square of the kept data equals the full chi-square minimized over the masked data. The prior on
the kept elements is the kept block of the prior covariance, inverted (`restrict_prior`, and `e044_fit.py:108-111`).
In `e044_fit.py:90-99` the window operator is applied first and the Schur complement second, for every band,
353 GHz included; region and block arrays are cut with the same `keep`. The order is the right one: the window
is not local.

**4.2 Degrees of freedom.** Write the prior as n extra data, 0 = c + e with Cov(e) = S. Then chi2_min with the
penalty is the residual of a linear fit of bands * n + n data with 3n amplitudes: bands * n - 2n at fixed index, less
one per index. `e044_fit.py:131` computes 5n - 3n - 4 + n = 17420 and 11612 without the prior. At fixed index the
difference of the two minima is c_hat^T (Cov(c_hat) + S)^-1 c_hat, chi-square with n degrees of freedom and
independent of the flat-prior residual. Toy (96 elements, 6 bands, dense prior, 400 draws): 381.3 +- 1.4 for 382;
285.9 +- 1.2 for 286; difference 95.4 +- 0.7 for 96 with variance 209 for 192; correlation with the flat chi-square
+0.02. The paper's statements hold when the CMB is drawn from the prior. With the 0.11 uK that the prior lacks the
expectation along the CMB spectrum rises by less than one percent (my estimate from a noise of about 2 uK^2 in
that direction, not a calculation on the matrices); the fit with the floor at 0.11 uK gives 1.0016.

**4.3 Eq. (mismatch).** With Sigma_w = N + U S_w U^T and P the projector that removes the synchrotron and dust
columns in the Sigma_w metric, g_r = -(M_r a)^T P eps + (1/2) eps^T P_r eps, eps ~ N(0, Sigma_t). Then
E[g_r] = (1/2) tr(P_r U dS U^T), which equals the appendix's -tr[dS (F_r Z)^T N^-1 R C], and
Cov(g) = J + q^T dS q + T + tr(P_r Sigma_w P_s U dS U^T) + (1/2) tr(P_r U dS U^T P_s U dS U^T), and the expected
curvature is J + (1/2) tr(P_rs U dS U^T). The paper keeps J + T + q^T dS q and J. So the term q^T dS q is exact in
dS, and two terms of first order in dS are dropped: the change of the noise-noise term and the change of the
expected curvature. `app.tex:67-68` admits the first. Both are of the order of T, not of J. At the toy
configuration of the unit test (T/J = 0.012): paper's form 0.1118, 0.0901; with both terms 0.1139, 0.0920; 1500
fits 0.1133 +- 0.0021, 0.0889 +- 0.0016. There the noise-noise term quadruples under the mismatch (1.6 to 6.6) and
it moves the sd by 1 to 2 percent. On the maps T/J is 0.04 at High Latitude and dK/J is 0.43, so I expect the
dropped terms to be worth about a point of the 22; I could not compute them without the matrices. "Exact
first-order" (`abstract.tex:35-36`, `e044_fit.py:15-16`, response to the referee) should be "first-order".

**4.4 The identity.** A change of reference multiplies the synchrotron column of every element of region r by
(nu / nu')^beta_r outside the band integral, so A' = D A D, chi2_min is unchanged and
-(1/2) ln det A' = -(1/2) ln det A + sum_r n_r beta_r ln(nu'/nu). A prior on the CMB block is untouched by D.
`e044_flat_scan.py:87, 96` uses this with the right sign; toy check 7e-12 and 2e-11; the file's own direct
evaluations agree to 4e-9.

**4.5 The nested yardstick.** For two estimators that solve profile scores with different treatments of the CMB,
the terms linear in the noise have Cov(s_arm, s_base) = (M_r a)^T P_arm Sigma_L P_L (M_s a) = J_arm, because
Sigma_L P_L = 1 - Pi_L and P_arm M = 0. So Cov(beta_arm, beta_base) = J_L^-1 and Var(difference) = Var_arm -
J_L^-1, to relative order T/J. This covers the white prior, the flat prior, K and Ka alone and the low-pass arms.
The paper subtracts the sandwich variance of the base, which is right if the noise-noise parts of the two scores
are fully correlated. The difference matters only for the white prior at High Latitude (-1.57 to -1.74).

**4.6 Why the marginal has no optimum.** One element, delta bands at the K to V effective frequencies, free
synchrotron and CMB amplitudes (`volume_term_toy.py`): the volume term -(1/2) ln det A rises with the index all the
way to -0.08, where the synchrotron column is nearly degenerate with the CMB column, and there it is 3.9 per
amplitude above its value at -3.04. At -3.04 its slope is +0.124 per amplitude and its curvature -0.065; the scan
of the real fit gives +0.133 and -0.059 (396 and -175 for 2974 amplitudes), and the toy pivot is 20.3 GHz against
20.1. Between -3.04 and -1.6 the toy gains 0.30 per amplitude, 892 for 2974; the scan has 908. So the mechanism is
the pull toward the degeneracy with a free CMB amplitude. Its slope is about 400 per unit index at the profile
optimum. The data alone (curvature 290) would stop that 1.4 units away, which is the edge; the volume term's own
curvature (-175 there, -813 at -1.75) lowers the net curvature to 115 and then below zero near -2.0, so nothing
stops it inside the range. The paper's sentence has the sign of the volume curvature right and everything else
wrong: the net curvature at the profile optimum is positive, and it is the slope that the data cannot balance. A
CMB prior removes the degeneracy, which is why the LCDM scan has a maximum. The toy has no dust and white noise.

**4.7 Chi-square by region and direction** (`chi2_by_region_direction.py`; from `chi2_region` of `lcdm_all.npz` and
`flat.npz`, both of which sum to their totals).

| Region | Two directions outside the model (flat-CMB fit), chi2 / 2n | Along the CMB spectrum (LCDM minus flat), chi2 / expected | z | z if the noise variance is 0.956 of nominal |
|---|---|---|---|---|
| High Latitude | 0.963 +- 0.018 | 2857.6 / 2945.6 = 0.970 +- 0.026 | -1.1 | +0.6 |
| Spur | 1.020 +- 0.061 | 244.2 / 263.5 = 0.927 +- 0.087 | -0.8 | -0.4 |
| Intermediate | 0.929 +- 0.024 | 1702.0 / 1681.8 = 1.012 +- 0.034 | +0.3 | +1.7 |
| Plane | 0.963 +- 0.034 | 993.2 / 861.7 = 1.153 +- 0.048 | +3.1 | +4.2 |

The Plane is as low as High Latitude in the two directions that the model does not span, and 15 percent high in the
third. Limits of this split: the parts of a dense quadratic form by region are treated as chi-square (the nside-8
null confirms the sampling error for the LCDM fit: 0.033 against 0.032 in the Plane; nothing confirms it for the
difference); the penalty of 55 is spread evenly (with none of it in the Plane the excess is 2.95 sd); it is a look
after the fact over four regions. In the K, Ka, 353 fit the one degree of freedom per element (the CMB spectrum as
K and Ka see it) has 0.957, 1.012, 0.941, 0.997: the Plane is not high there, so the excess is in Q and V against
the K to Ka power law. The K, Ka arm of Table
`tab:moves` says the same thing in the index (+0.025, +2.0 nested sd). Candidates: leakage of intensity, an index
that varies from pixel to pixel, curvature. The block test did not see it (1.28, 78th percentile of the null).

**4.8 Regions against the text.** `lowres/regions.py` and `check_regions.py`: \|b\| > 30, 10 <= \|b\| <= 30, \|b\| < 10 on
pixel centres; Spur = surrogate Spur pixels with b > 20 and -45 < l < 65, 148 of 464, spanning l from -20 to 48 and
b from 22 to 72; mask = 154 pixels at or above the 95th percentile of the smoothed DR1 amplitude (28.6 uK_RJ; 128 in
the Plane, 15 in the Spur, 11 at Intermediate latitude) and 14 pixels of the Tau A disc; the arrays equal those
stored in `lcdm_all.npz`. The description in `sec_data.tex:57-64` is accurate. Two things it does not say: the
Plane keeps 435 of 576 pixels, so it is the plane outside its brightest quarter; and the Spur's kept pixels have a
median amplitude of 22 uK against the threshold of 28.6, so the threshold cuts through that region. Selection on
the data: the threshold is applied to a 5 degree average at a signal to noise of order 100, so selection on noise
is negligible; no definition uses an index.

**4.9 Block errors.** Eq. (block) is what `lowres/block.py:99-113` computes; I reproduce the four own-score values
and effective numbers from the stored G, D, J, T.

## 5. Pre-registration

- H1 to H8: the text between "## Expectations" and "## Command" at be81312 (09:00:35) is unchanged at head (diff is
  additions only). The error file of the first nside-16 job is stamped 09:00:37 and its result 09:43:06.
- H9: committed 4a2c8a7 at 10:13:23; unchanged since. No nside-16 simulation existed before 12:35. Clean.
- H10: committed 0c11355 at 10:24:38; unchanged since. The third null set was launched before the commit (draw at
  10:23:43, fits at 10:24:34); the first of 240 results could exist about 30 s after the commit and the file was
  read at 12:36. Clean, by seconds.
- H11: committed d0d3d7a at 10:57:32, amended 0df685b at 10:58:52; the scan started at 12:08:15. Clean. Its
  expectations hold. Its premise ("the curvature of the marginal at the profile optimum is negative there") is
  contradicted by the scan it ordered, R7 did not notice, and the paper repeats the premise (row 8).
- Failed expectations that the paper does not report: H4 for Intermediate and Plane (written 1.5 to 3 and 2 to 4,
  found 1.44 and 1.28; the response to the referee says the paper states this, and it does not); the two accuracy
  clauses of H1 (within 15 percent below a move of 0.3: Spur 16 percent; J_c "worse by up to 25 percent": 31
  percent low at High Latitude); H6 for the old Spur (-0.019 where "less than 0.01" was written); H8 on where the
  dust index matters most. `sec_data.tex:98-99` says the text reports failures "where it matters"; H4 matters,
  because the author expected two regions to fail the one-index model and they did not.
- Found after the fact and not labelled so in the paper: the total chi-square and its significance (R4: "A check I
  should have made"; H3 allowed 10 percent, which a 2.7 percent deficit passes); the split along the CMB spectrum
  (noted 11:40); the "sequence in the CMB prior"; the disagreement of the zero-mean prior (R1: "Not expected").
  Only the dust index is labelled (`sec_data.tex:475-476`). `sec_data.tex:92-94` reads as if every statement of
  the section had an expectation.
- H7 is reported fairly: the number is given with its sign and the reader can see that it missed the interval.

## 6. References

- Lattanzi et al. 2017: authors, JCAP 02 (2017) 041, arXiv:1611.01123 confirmed (arXiv record and DOI). Table 1
  confirmed in the PDF. The sentence at `sec_data.tex:204-206` is stronger than the source in one word: they do
  comment ("consistent with unity for all bands, in all masks considered"). Their Q band is below one in three of
  four entries (0.978 to 1.004), so "near or above one" is fair for Ka and loose for Q. Their fit is on the sky
  outside the WMAP processing mask with one scaling per band; it is not a per-pixel fit.
- Planck Collaboration V 2020: A&A 641, A5, arXiv:1907.12875 confirmed. Table 12 values 1.004, 1.048, 0.970
  confirmed in the ar5iv HTML through a small model, not in the PDF. "With the same matrices" is loose: that paper
  smooths the maps with a cosine window and says it adopts the official matrices smoothed in harmonic space.
- The research note found more than the paper uses, and the omissions lean one way. Not cited: (i) Jarosik et al.
  2007, Sect. 3.6, that the degraded maps have "noise properties not precisely described by the noise matrix", and
  their Table 5, where the r4 noise is 1.5 to 5.7 percent above the r9 value; (ii) the WMAP team's own pixel-space
  values, all at or above one (three-year, outside P06: 1.000 to 1.018 for Ka, Q, V in the pixel model, 1.016 to
  1.096 for the cleaned maps; nine-year null maps 1.009 +- 0.040; all as transcribed in the research note, not
  read by me); (iii) Hinshaw et al. 2009,
  "a slight tendency to overestimate the uncertainties at higher l", which supports the reading; (iv) Svalheim et
  al. 2023, footnote 5 of Sect. 5.1, a "uniform and slightly negative bias at high Galactic latitudes" in the
  chi-square, attributed there to the count of degrees of freedom. (iv) is the closest precedent, has the same
  sign and place, and is the author's own paper. "We found no statement of the accuracy of the normalization" is
  true as worded and incomplete without (i).
- Svalheim et al. 2023: the entry says "BeyondPlanck XIV" with A&A 675, A14, which is the published title; arXiv
  v2 is titled "BeyondPlanck XV". Say both, since the repository calls it XV throughout. The mask sentence
  (`sec_data.tex:65-66`): the arXiv HTML (through a small model, not the PDF) describes three steps: the strongest
  5 percent of the synchrotron amplitude smoothed to 4 degrees, a smoothed chi-square map thresholded at its
  strongest 10 percent "to remove strong sources such as Tau A", and an unsmoothed chi-square map for compact
  sources. The paper's mask has the first step (at 5 degrees, on another model) and a disc in place of the other
  two. "Ours follows that description" should say what it leaves out.
- Jarosik et al. 2007, 2011 for the matrices: "de-weight the poorly constrained transmission-imbalance modes" is
  the neutral wording that the verification file asks for (the Supplement and the LAMBDA page disagree on
  projection against a one-sigma term) and agrees with the file check (no null space). Bennett et al. 2013 is
  cited for the maps only. Watts et al. 2023: the sentence matches the quoted passage. Nothing in Sect. 6
  attributed to these four is stronger than the verification files, except the two items above.
- `sec_data.tex:20-21`: "we use them as released and test the sensitivity to the modes they flag." No such test is
  in the paper or in E044 (R4: "Not checked: whether the deficit is in particular modes").

## 7. What a hostile expert would say

**7.1 Scope.** "A table captioned 'The polarized synchrotron index of the WMAP nine-year maps', an index in the
abstract, a dust index with an error, a statement about the WMAP noise, and a sentence that says nothing was
measured." The paper answers in part (`sec_data.tex:488-494`, `sec_concl.tex:32`). What it does not answer: the
reason it gives for not delivering an index (`sec_data.tex:6-9`, choices that "move the index by as much as its
error or more") is no longer what Table `tab:moves` shows. After the mask the only such rows are the flat CMB
prior, which the paper calls noise, and the delta bandpasses, which it calls inexcusable. The honest reasons are
others: one mask and one set of regions with no variation of either; a 3 sd excess along the CMB spectrum in the
Plane; two estimators 3.6 sd apart at High Latitude with a third of it explained; a chi-square that is not
understood. One paragraph with these, and a caption that says "index of the one-power-law model", would do.

**7.2 The mask.** "A recipe, not a file, written by an author of the paper whose mask it imitates; a threshold on a
model made from the same data; and no row in Table `tab:moves` for the threshold." The first is stated. The second needs
one sentence (4.8). The third is not answered and is the obvious check after a draft whose result changed with the
mask: thresholds of 3 and 10 percent, and the present regions without the mask (`--no-mask` exists). The released
BeyondPlanck mask and regions would remove the word "reconstructed" from the paper.

**7.3 The low chi-square.** "Over-parameterized?" Answered: exact count with a flat CMB prior, a simulated null,
and the toy of 4.2. "Would the WMAP team accept the reading?" They would say that the matrices are a unit-lag
filter product scaled by a sigma_0 measured at high resolution, that they wrote in 2007 that the degraded maps are
not precisely described by them, that their own values are at or above one, and that 2 to 4 percent in variance
is inside what the nine-year change of sigma_0 alone moved (1.4 to 2.9 percent, per the research note). They would
ask for the test that needs no sky model: year or DA differences with the single-year matrices, which the
Supplement lists as released. And they would point at 4.7: a deficit that is absent in one direction and an
excess in the Plane is a statement about the fit before it is one about their noise. One sentence: "not planned;
if it is noise, it is not uniform over band combinations or regions; a model-free test needs the single-year
maps". The abstract should not carry the percent.

**7.4 Four regions.** "The index is known to vary; one index over 48 percent of the sky is a choice, and your test
of it has a dozen patches." Partly answered (`sec_data.tex:367-379`). Missing: what spread between patches the
block test would have detected, which E040 stage 5 measured and the paper no longer quotes; and 4.7, which is the
more sensitive diagnostic and does show something in the Plane.

**7.5 Nside 16 and the 1 degree window.** Answered: bias 0.005 to 0.011 without the window matching, below 0.0025
with it, against errors of 0.015 to 0.029. A reader will still ask what the window means on 3.7 degree pixels;
one clause ("to multipole 47; above it the pixel average is untouched") is in `lowres/window.py` and not in the
text.

**7.6 Leakage and bandpass mismatch.** "WMAP solves for a spurious map per differencing assembly; did you use it?"
The paper says leakage is not modelled. One sentence would help: the released Q and U matrices are marginal over
the spurious map, so first-order bandpass mismatch is in the covariance (the verification file has this), and
what is not is beam-asymmetry leakage near bright intensity, which is where 4.7 points. The promised test of the
imbalance modes (`sec_data.tex:20-21`) does not exist.

**7.7 K band.** "Everything hangs on K: the pivot is at 23.3 GHz and K has the largest known systematics." The
K, Ka arm does not address this (it keeps K). The sections that did were removed. One sentence saying so.

**7.8 Novelty.** Tempering law: new in this setting, placed against Ball et al. and Bardenet et al.; fine.
Identity: "a flat prior is not invariant under a change of variable; one line." The paper's claim is the closed
form of the shift and its pivot, which is fair, and Sect. 4 cites Dunkley et al. for the effect. The runaway
with a free CMB amplitude is the extreme of what Dunkley et al. saw (-2.6 in simulations with uniform priors on
both amplitudes); say so, with the degeneracy of 4.6, or an expert will say it for you. Sandwich and debiasing:
now placed in the measurement-error literature; fine. Prior mismatch: standard misspecification algebra; the
contribution is the number at this configuration, so the word "exact" is not needed and is wrong.

**7.9 "An independent referee".** The phrase occurs five times in Sect. 6 and the conclusions, and the story of
the earlier regions is told in the abstract, the introduction, Sects. 6.1, 6.3, 6.4 and the conclusions. For a
journal: once, in Sect. 6.1, as an internal review.

## 8. Required changes

**Fatal**

None.

**Major**

1. `abstract.tex:36`, `sec_data.tex:193-199` and `300-307`, `sec_concl.tex:34-35`. "No region is misfitted" needs
   its other half. Give the split of 4.7 by region next to the total (Table `tab:cmbprior` already prints 0.614
   against 0.385 for the Plane and nobody is told to subtract them). Say that the Plane carries 15 percent more
   chi-square than expected along the CMB spectrum (3 sd; 4 sd under the noise reading), and that "no region has a
   chi-square above its expectation" is all that the totals show. Add the flat-CMB fit and the regional difference
   to the nside-8 chi-square null so that the sampling error of this split is measured.
2. `abstract.tex:39-41`. Remove the percent from the abstract or make it conditional and say the test was not
   planned. In `sec_data.tex:193-207` and `297-307`: add the K, Ka value (0.962 +- 0.019 along a CMB-spectrum
   direction), give the significance of 1.008 against 0.956 (2.3 sd), replace "none along the CMB spectrum" by
   the numbers, cite Jarosik et al. 2007 Sect. 3.6, Hinshaw et al. 2009 and footnote 5 of Svalheim et al. 2023,
   and correct "without comment" for Lattanzi et al.
3. `sec_data.tex:313-316`; record R5, H11 and R7; any finding that repeats it. The curvature at the profile
   optimum is positive (about +115 against J_c = 292); -366 belongs to the bound. Rewrite with the slope (+396,
   pivot 20 GHz) and the degeneracy of 4.6. In the generator, do not print a Hessian evaluated on a bound as a
   property of the profile optimum.
4. `sec_flat.tex:177-178`, `paperA1_make_e044.py:188`. Replace "16.9 to 40.2 GHz" by the value from the scan, about
   20 GHz at High Latitude, and do not define `flat.pivot.*` for a fit on the bound.
5. `abstract.tex:37-38` and `41-43`, `sec_data.tex:445-466` and `485-486`, `sec_concl.tex:36-37`. State the
   sensitivity of the scale test (two nested sd at l <= 20 are 0.036, 0.059, 0.069, 0.063, against old effects of
   0.024 to 0.084) and say that the change is clear in the Plane only. "Went away with the emission that the mask
   removes" needs the fit that separates mask from regions: the present regions with `--no-mask`, all scales and
   l <= 20. For the band subset say "2.4 and 2.0 nested sd in two regions, same sign", and say that the K, Ka arm
   and the flat-CMB arm are one pattern (the index from the K to Ka lever arm is flatter), not two.
6. `abstract.tex:38-39`, `sec_concl.tex:38`, `sec_data.tex:486-487`. "Moves by up to 0.11 with the prior on the
   CMB" beside "1.95 nested sd, not significant". Choose: either a sensitivity, and then the K, Ka change is one
   too, or noise, and then it does not belong in the abstract as a move.
7. `abstract.tex:22-23`, `sec_profile.tex:116-117` and `139-141`, `sec_concl.tex:24-27`. Say "in simulations at
   half the resolution of the map fit" until H9 is read; delete "and more on the maps" or replace it by "and the
   two estimators differ by 3.6 sd on the maps"; say in Sect. 6.3 what H9 decides (Sect. 2, pending item 3).
8. Table `tab:moves` and `sec_data.tex:468-481`. Add the mask as a model choice: bright fraction 0.03 and 0.10,
   and no mask. Without it the paper has varied everything except the choice that changed its result.
9. `sec_data.tex:20-21`. Do the test of the flagged modes (deproject the released imbalance templates) or delete
   the clause.
10. `sec_data.tex:4-9` and the caption of Table `tab:same`. Replace the stale reason by the present ones (7.1);
    caption: "Index of the one-power-law model fitted to ...".
11. `sec_data.tex:395-411` and `abstract.tex:36-37`. The nside-16 null of six realizations cannot carry "agree"
    (Sect. 2, pending item 1). Either run enough to place 1.63, 1.44 and 1.28 (40 or more), or write "inside a
    null simulated at half the resolution; six realizations at full resolution do not contradict it", if they do
    not.
12. `sec_data.tex:92-99`. Name what was found after the fact (total chi-square, its split, the CMB-prior
    sequence, the zero-mean disagreement) and report the failure of H4. Correct the response to the referee,
    which says the paper already does.

**Minor**

13. `abstract.tex:35-36`, `e044_fit.py:15-16`: "first-order", not "exact first-order"; `app.tex:67-68`: name both
    dropped terms and their order (T/J).
14. `abstract.tex:35`, `sec_data.tex:278`, `sec_concl.tex:28-29`: "the true error is 9 to 22 percent larger than
    the quoted one".
15. `sec_data.tex:326-327`: 0.23 to 0.60. `sec_data.tex:82-83`: 0.0020 or less. `sec_data.tex:405`: 30 pixels.
    `sec_data.tex:314-315`: J_c = 292. `sec_data.tex:240-241`: three regions. `sec_data.tex:198-199, 300`:
    two-sided. `sec_data.tex:147-149`: scatter, not error.
16. `sec_data.tex:270`: "1500 simulated fits at the configuration of the unit test", and raise the committed test
    to match or quote its 500.
17. `sec_data.tex:280-284` and the caption of Table `tab:moves`: the yardstick holds to relative order T/J and
    assumes the LCDM covariance is the truth; -1.57 to -1.74.
18. `abstract.tex:27`: "at \|b\| > 30 deg outside the North Polar Spur". `sec_data.tex:57-64`: the Plane keeps 435 of
    576 pixels; the threshold passes through the Spur.
19. `sec_data.tex:65-66`: what the BeyondPlanck mask has that this one lacks (4 degrees, two chi-square steps).
20. `sec_data.tex:395` against Table `tab:same`: a Spur ratio is quoted where no Spur block error is given.
21. Fig. `fig:blocks`: 4 of 38 High Latitude patches are drawn. The panel does not show what the ratio of 1.63 is
    made of. Raise the cut or plot scores in units of their model error.
22. `sec_data.tex:296-297`: "at the edge of what its calculation allows" should be "1.7 to 1.9 of the sd expected
    for it".
23. `sec_data.tex:471-476`: label the error of the dust index (Delta chi-square = 1, three points, fixed
    temperature, power-law synchrotron).
24. `sec_concl.tex:47-49`: drop "together" or trace one published value.
25. `refs.tex:39`: "published as BeyondPlanck XIV; arXiv v2 is titled BeyondPlanck XV".
26. Abstract: about 650 words. Cut to 250: one sentence per method, then the three data sentences that are
    differential (rows 6, 7, 12), then one on what the fit passes and one on what it does not.
27. `abstract.tex:36-38`: "errors from the scatter of the score over patches of sky agree with the model's" and
    "the angular scales or the band subset that inform it" use terms the abstract has not defined. Suggested:
    "block errors from 48 patches agree with the model errors within their simulated scatter, and restricting the
    fit to multipoles below 30, 20 or 12 or to K and Ka changes the index by at most 2.4 of the sd expected for
    such a change".
28. `abstract.tex:40-41`: "which as a statement about the released covariances is noise overstated by" does not
    parse at first reading.
29. Generator: remove `lp*.chi.ratio`, `flat.pivot.*`, `flat.ratioH`, `flat.ratioJ`; compute the penalty's
    expectation or say that the observed one is used.
30. Record: line 49 (128, not 141), line 474 (30 pixels), R6 (8.2), R5 and H11 (item 3).
31. No em-dash found in the sources.

## 9. What I did not check

- I ran no fit and loaded no covariance. Everything about the maps comes from the saved `.npz` and text files. The
  split of 4.7 uses stored per-region sums; its sampling error is the paper's own approximation and is not
  simulated for the difference of two fits.
- The dropped terms of Eq. (mismatch) are computed at the toy configuration only. Their size on the maps is an
  estimate.
- Sects. 3 to 5 and the appendix up to line 47: I read the diff against draft 3 and checked Table 1 for internal
  consistency. I did not re-derive the tempering law, the record process, or the sandwich.
- E035, E036, E041, E042, E043 and the figures `fig_temper`, `fig_ladder` were not rechecked. Of E038, E039, E040
  only the numbers that Sect. 6 quotes. `sec_impl.tex` was read, not checked (its "0.25 and 0.77 MB" against the
  record's "352 kB to 768 kB" is not resolved by me).
- The low-pass construction with a singular inverse covariance and the pseudo-inverse in the mask: I relied on the
  statistical review of 2026-11-11 and did not re-derive it.
- The bandpass readers, the window operator and `lcdm_cmb_precision` were read, not tested. I did not check the
  CAMB spectra or the 353 GHz noise treatment.
- Planck 2018 V Table 12 and the mask of Svalheim et al. were read in HTML through a small model, not in the PDFs.
  Lattanzi et al. Table 1 and p. 4 were read in the PDF. Jarosik et al., Hinshaw et al. and the Supplement are
  quoted from the research note, not from the sources.
- The nside-16 simulations were running while I read; the four realizations quoted in Sect. 2 are a snapshot of
  13:05, not a result. Whether the profile optimum is unbiased at nside 16 on these regions is not established by
  anything I saw.
- The unit tests of `nextgen/lowres/tests` were not run (a long job held the memory); I ran my own toy versions.
- The compiled PDF: pages 14, 20 and 21 were looked at; the text of the whole was searched for unresolved macros
  (none) and for the three red items (present). Other pages were not compared with the sources.
- Anything committed after 9fed07f.
