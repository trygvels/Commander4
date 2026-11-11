# Response to the referee report on Paper A1, draft 1

Report: `2026-11-11_referee_report_paperA1_draft1.md`. This file lists, item by item, what was done and
where. Draft 2 (commit 5ae2489) made the text changes; draft 3 (this night, 2026-11-11 03:00 to 07:00)
added the computations that items 2, 7, 8 and 9 asked for, and one item nobody asked for (E043). Every
number in the three tables of Sect. 6 is now written by `nextgen/scripts/paperA1_make_tables.py` from
result files; none is typed.

**Not final (05:45):** the computations behind items 2, 4 (B2) and 9 are queued or running; their rows say what
was asked of them, and the outcomes are added when they are read. The document is not handed to a referee before that.

## Required items

| # | Item | Status | What was done | Where |
|---|---|---|---|---|
| 1 | Withdraw the quoted indices or replace them by what the cuts support | done | No index is quoted as a measurement anywhere. The paper was reduced to the methods and a differential comparison of estimators on the same maps (D032). | abstract, `sec_data.tex` intro and `sec:moves`, `sec_concl.tex` |
| 2 | Validate the low-pass anchor estimator on simulations before quoting numbers from it | computed | E039 addendum B: 100 realizations at nside 8 (L = 12), and a spot check at nside 16 with the LCDM prior (L = 20). Pre-registered B1 to B4. The low-pass rows of Table `tab:moves` remain sizes of an effect, not corrections. | `E039_lowres_dense_validation.md` addendum B; `E039/lowpass_n8_l12.txt`, `E039/lowpass_n16_l20.txt` |
| 3 | Demonstrate the aliasing mechanism or call the scale dependence unexplained | done (the second) | "The cause is aliasing" removed. The text names two candidates and says which way Fuskeland et al. (2014) point. FINDINGS F49 amended and regraded. | `sec:moves`, paragraph "Angular scale"; FINDINGS amendment of 2026-11-11 |
| 4 | Repair the BeyondPlanck and Cosmoglobe reading | done, B2 computed | The contradiction with Sect. 4.3 is removed; the published high-latitude values are called prior draws; "fitted freely" is attributed to Watts et al. with calibration; the Spur and Plane counter-evidence is stated; "does not need unmodelled systematics" is deleted; the reading is out of the abstract. B2 (white CMB prior) and the LCDM-prior version are run. | `sec:bpset`; `E040/s3_bpset_cmb02_wm.txt`, `E040/s3_bpset_lcdm_wm.txt` |
| 5 | Correct the Commander3 statement to agree with E009 and F10 | done | Rewritten from the records. | `sec_legacy.tex`, last paragraph |
| 6 | Remove "differing by 0.5 between the 2018 and PR4 releases" | done | Removed from the abstract; the two fits are reported as bad fits on either side, not as a measured difference. | abstract, `sec:30` |
| 7 | LCDM covariance or flat prior on the CMB for every quoted fit, the other beside it, sensitivity printed | done | The LCDM covariance is the central choice (D034, by a rule written before the fit; my criterion C1 failed). Table `tab:cmbprior` shows LCDM, white and flat. New: a simulation of the three priors on LCDM skies (white prior unbiased, errors 24 to 32 percent too small). Arms repeated with the LCDM prior. | `sec:same`; E040 stage 4, C4, 4b; FINDINGS F53 |
| 8 | Errors that survive misfit, chi-square per region, no error bar for the Galactic Center | done | Block errors from the scatter of the score over 48 sky patches with a leverage correction, checked on simulations with and without a true spread; chi-square per region in the table; the Galactic Center has no block error and its model error is in parentheses. | `sec:same`, Eq. `eq:block`; E040 stage 5; FINDINGS F54 |
| 9 | Say why the released PR4 covariances are not used; include the noise-scale arms | done | The sentence is in `sec:30`. Arms: white-noise columns only, and the power-law term doubled, for K with 30 GHz under the LCDM prior. | `sec:30`; E040 stage 4c |
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
