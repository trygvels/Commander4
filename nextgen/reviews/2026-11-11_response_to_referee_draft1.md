# Response to the referee report on Paper A1, draft 1

Report: `2026-11-11_referee_report_paperA1_draft1.md`. This file lists, item by item, what was done and
where. Draft 2 (commit 5ae2489) made the text changes; draft 3 (this night, 2026-11-11 03:00 to 07:00)
added the computations that items 2, 7, 8 and 9 asked for, and one item nobody asked for (E043). Every
number in the three tables of Sect. 6 is now written by `nextgen/scripts/paperA1_make_tables.py` from
result files; none is typed.

**State at 10:00:** superseded in part by the second referee read and `2026-11-11_response_to_referee_draft3.md`: the real-data section was redone on
new regions with a mask (E044), so every statement below about Sect. 6 describes draft 3, not the current draft.

## Required items

| # | Item | Status | What was done | Where |
|---|---|---|---|---|
| 1 | Withdraw the quoted indices or replace them by what the cuts support | done | No index is quoted as a measurement anywhere. The paper was reduced to the methods and a differential comparison of estimators on the same maps (D032). | abstract, `sec_data.tex` intro and `sec:moves`, `sec_concl.tex` |
| 2 | Validate the low-pass anchor estimator on simulations before quoting numbers from it | done at nside 8, spot check at nside 16 running | E039 addendum B, 100 realizations at nside 8 (L = 12): mean z within one standard error of zero, sd of z 0.95 to 1.03, coverage 0.66 to 0.75, and the difference from the all-scales fit scatters as sqrt(sd_low^2 - sd_all^2) (0.89 to 1.12 of it). The independent stats review re-derived the construction and confirmed it on 400 realizations of its own. The low-pass rows of Table `tab:moves` are sizes of an effect; on the maps they are 3 to 10 times what the model allows for nested fits, with either CMB prior. | `E039_lowres_dense_validation.md` addendum B; `E039/lowpass_n8_l12.txt`; `reviews/2026-11-11_stats_review_prior_block_lowpass.md`; nside 16: `E040/sim_n16_lcdm.txt` (running) |
| 3 | Demonstrate the aliasing mechanism or call the scale dependence unexplained | done (the second) | "The cause is aliasing" removed. The text names two candidates and says which way Fuskeland et al. (2014) point. FINDINGS F49 amended and regraded. | `sec:moves`, paragraph "Angular scale"; FINDINGS amendment of 2026-11-11 |
| 4 | Repair the BeyondPlanck and Cosmoglobe reading | removed from this paper | Draft 2 repaired the text; B2 was queued and never completed. After the second referee read the whole subsection on the band set without K band, and the one on Planck 30 GHz, were taken out of Paper A1 (D035) and go to Paper B with E042. (An earlier version of this row said "done, B2 computed" and named two output files that did not exist; the second referee caught it.) | `DECISIONS.md` D035 |
| 5 | Correct the Commander3 statement to agree with E009 and F10 | done | Rewritten from the records. | `sec_legacy.tex`, last paragraph |
| 6 | Remove "differing by 0.5 between the 2018 and PR4 releases" | done | Removed from the abstract; the two fits are reported as bad fits on either side, not as a measured difference. | abstract, `sec:30` |
| 7 | LCDM covariance or flat prior on the CMB for every quoted fit, the other beside it, sensitivity printed | done | The LCDM covariance is the central choice (D034, by a rule written before the fit; my criterion C1 failed: the Plane moved by 2.3 sd). Table `tab:cmbprior` shows LCDM, white and flat. A simulation of the three priors on 60 LCDM skies: matched prior calibrated; white prior unbiased with errors 24 to 32 percent too small; flat prior calibrated at 2.2 to 2.5 times the error. The arms of Table `tab:moves` are repeated with the LCDM prior. The review found that the prior stops at multipole 47; with the power above that as its floor the indices change by 0.0006 or less. | `sec:same`; E040 stage 4, C4, 4b; FINDINGS F53 |
| 8 | Errors that survive misfit, chi-square per region, no error bar for the Galactic Center | done, nside-16 reference running | Block errors from the scatter of the score over 48 sky patches with a leverage correction (`lowres/block.py`, unit-tested), checked on simulations with and without a true spread, reviewed independently. The review showed that the reference for the ratio under the model must come from the configuration of the fit; with the LCDM-prior reference at nside 8 (60 realizations) the High Latitude region and the Plane are outside it and the Spur is inside, and the paper says so. Chi-square per region is in the tables; the Galactic Center has no block error and its model error is in parentheses. A figure shows the per-patch deviations. | `sec:same`, Eq. `eq:block`, Fig. `fig:blocks`; E040 stage 5; FINDINGS F54 |
| 9 | Say why the released PR4 covariances are not used; include the noise-scale arms | done | The sentence is in `sec:30`. K band with PR4 30 GHz under the LCDM prior: reduced chi-square 2.64 with the noise model, 4.19 with the white-noise columns only, 2.32 with the correlated term doubled; the high-latitude index stays at -3.40 to -3.44 against -3.06 for K with Ka. The misfit does not rest on the noise model. | `sec:30`; E040 stage 4c |
| 10 | Correct the statement on when the expectations were written | done | The paragraph "Expectations" now separates what was written before the first real-data fit from what was written after earlier real-data results. | `sec_data.tex`, "Expectations" |
| 11 | Sign in the appendix, prior precision in A, projector name | done | | `app.tex` |
| 12 | "exactly an exponential prior" | done | Reworded. | abstract |
| 13 | 13.1 and 13.6 are from a toy with two indices and the recorded steps | done | | `sec_legacy.tex` |
| 14 | "to 0.2 percent" in the abstract | done | Replaced by what Table 1 shows. | abstract |
| 15 | Number mismatches (rows 14, 15, 20, 24, 26, 36, 39, 42, 45, 47) | done, and superseded for Sect. 6 | The Sect. 6 tables are generated from the result files. The generation found three more rounding slips of mine in hand-typed numbers (Spur -3.097 for -3.0965, a reference shift of 0.037 for 0.0365, a ratio 2.8 for 2.75), which is why they are generated now. | `paperA1_make_tables.py`, `gen_tab_*.tex` |
| 16 | Bardenet et al.; Alonso et al. | done | | `sec_intro.tex`, `sec_legacy.tex` |
| 17 | References | done | BeyondPlanck XIII, Planck 2015 X, Dickinson et al. 2009 corrected. New in draft 3: Liang and Zeger 1986, Bell and McCaffrey 2002, Cameron and Miller 2015, Carter et al. 2017 (verified, `research/2026-11-11_citation_verification_cluster_robust.md`), Lewis et al. 2000, Planck 2018 VI. | `refs.tex` |
| 18 | Conditions of the zero-mean identity | done | | `sec_profile.tex` |
| 19 | Boxed shift formula is for one component | done | | `sec_flat.tex` |
| 20 | Forward reference in `sec_profile.tex` | done | | `sec_profile.tex` |
| 21 | Wording items | done | | several |

## Suggested items

- Split the paper: done (D032).
- Cite E022 for the band-set reading: done in `sec:bpset`.
- Say that the estimator comparison is differential: done, and the table is now for the adopted version.
- A table of every analysis choice that moves the index: `tab:moves`, generated.
- The Commander3 `marginal` sign error: left out of this paper; it is in FINDINGS (F10) and for Trygve to
  decide.

## What is new in draft 3 and was not asked for

- **E043.** A heap overflow in Apple Accelerate's Cholesky (upper-triangle option, matrix orders 1500 to
  6000, macOS 27.0 pre-release) under every nside-8 fit. Found from three crashes, measured with a malloc
  canary, fixed by using the other triangle. Every exposed result that the paper uses was rerun and
  compared as printed: no number changed. The nside-16 fits were not exposed. A paragraph in `sec_impl.tex`.
- **An independent statistical review of the three new pieces** (LCDM prior, block errors, large-scale estimator): no error in algebra or code; six
  statements of mine were worded too strongly and are corrected in the text; eight unit tests added.
- **The white CMB prior makes the error bar too small.** On LCDM skies the white-prior fit is unbiased and
  its error is 24 to 32 percent too small. That applies to the errors quoted in drafts 1 and 2.
- **Flat amplitude priors with a free CMB on WMAP alone.** The High Latitude optimum is -2.51 (reference 23
  GHz) and -1.71 (30 GHz) where the profile likelihood gives -3.03; the first-order prediction of the
  difference holds within 10 percent for shifts up to 0.8.

## Still open, and said so in the paper

- The cause of the scale dependence of the index.
- Simulation-based calibration of the ridge and marginal modes inside Commander4 itself (not run).
- The forward test of the legacy widths with amplitudes stored as maps (not run).
- An isotropic, empirical noise model for 30 GHz.
