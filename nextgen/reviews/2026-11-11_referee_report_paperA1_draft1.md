# Referee report on Paper A1, draft 1

- **Written:** 2026-11-11, reading started 02:07 CET, report finished 02:37 CET (times from `date`)
- **Referee:** independent agent, hostile brief. I had not seen how the paper was written.
- **Draft reviewed:** `nextgen/papers/paperA1/*.tex` as last modified 2026-11-11 02:06 (commit 9ed43ba),
  and `main.pdf` of 02:06 (23 pages). Line numbers below refer to those files.
- **Evidence read:** records E035, E038, E039, E040, E041, E042 in full; E036 in its result tables; the
  parts of E009, E010, E014, E016, E017, E020, E021, E022, E032 that the paper cites; raw outputs
  named in Sect. 3; `reviews/2026-11-09_stats_review_legacy_move.md` (Sects. 2 to 5),
  `reviews/2026-11-10_stats_review_lowres.md` (verdict and derivations); both citation reports
  (summary tables and the items named in Sect. 4); the two prior-art notes of 2026-11-09 in part;
  `nextgen/lowres/likelihood.py` and `nextgen/scripts/e042e_lowpass.py`, `e042c_pr3_arm.py`.
- **Own checks:** `nextgen/reviews/referee_scripts_2026-11-11/` (`check_sec3_temper.py`,
  `check_sec45_app.py`, `check_cmb_prior_shape.py`, each with its `.out`; a few seconds each, one or
  two threads, well under 1 GB).
- **One result here is newer than the draft.** `E042/lowpass12.txt` finished at 02:25, after the PDF
  was compiled. I read it only once its final table and `.npz` were on disk and the queue had moved
  on to L = 30. I did not read `lowpass30.txt`, which was still being written.

Convention: "the paper says" quotes the draft; "the record says" quotes an experiment record; "(my
check)" is a calculation of mine; "(from memory)" is a statement I could not verify here and that
the authors should check before relying on it.

## 1. Verdict in ten lines

1. **Fatal objection: yes for one central claim, no for the methods.** Sects. 3 to 5 are right as far
   as I could re-derive and reproduce them.
2. What fails is the quoted WMAP index, "-3.11 +- 0.01 at high latitude, -3.11 +- 0.02 in the North
   Polar Spur and -3.06 +- 0.01 in the Galactic plane" (`abstract.tex:34-35`, Eq. `eq:result`). It
   comes from an estimator built after a failed control, never run on a simulation with known truth,
   and adopted before the stability test that the record made the condition for adopting it (E042,
   expectation Z4).
3. That test has now come in and fails. From L = 20 to L = 12 the K to V index moves by -0.029,
   +0.004, -0.040 and -0.019 (2.8, 0.6, 5.7 and 3.8 nested standard deviations; 3.1 in the Galactic
   Center by the record's own looser rule). All-scale, L = 20, L = 12 give -3.094, -3.113, -3.142 at
   high latitude. The index keeps steepening as small scales are removed. There is no stable value
   at the quoted precision.
4. The cause named for the effect, "The cause is aliasing" (`sec_data.tex:243`), is a hypothesis. The
   only simulation of sub-pixel structure in the records leaves 0.002 after window matching, and
   Fuskeland et al. (2014), with beam-matched 1 degree maps, find -2.98 +- 0.01 in the plane, which
   is the all-scale K, Ka value (-2.980) that the paper now calls biased.
5. Strongest point 1: the tempering law T = 1/2 + 1/(4 c_k), its multi-move form and the exact
   controls. I re-derived all of it and reproduced every number I tried.
6. Strongest point 2: the reference-frequency identity and its demonstration on identical maps
   (-3.10 against -2.94 at high latitude, predicted to 8 percent). This comparison needs no
   assumption about the sky and survives every objection below.
7. Strongest point 3: the plug-in trap E[J(a_hat)] = J + T with the coverage ladder, and the record
   keeping: expectations written before runs, failed predictions reported, independent
   re-derivations.
8. Weakest point 1: lines 2 to 4, made worse by a white CMB prior that is wrong by factors of 30 to
   200 at multipoles 2 to 5 (my check), exactly where the low-pass fit takes its information.
9. Weakest point 2: the reading of the BeyondPlanck and Cosmoglobe flat index. It contradicts the
   paper's own Sect. 4.3 and Sect. 3.5, rests on the pending run B2, misattributes a quotation, and
   does not explain why the Spur and Plane of those analyses were not pulled.
10. Weakest point 3: the 30 GHz statements. "Differing by 0.5 between the 2018 and PR4 releases"
    comes from fits with reduced chi-square 2.6 and 6 under an empirical isotropic noise model, and
    at high latitude it contradicts the paper's own cross-spectrum ratios, by which the two releases
    agree to 1 to 3 percent.

## 2. Check of the mathematics

Scripts and outputs are in `referee_scripts_2026-11-11/`. Where I write "verified" I re-derived the
result by hand and, where a number is attached, reproduced it numerically.

### Sect. 2 (model and targets)

- Eq. `eq:chisplit`, `eq:Lflat`, `eq:Lprof`: verified. The statement that the profile likelihood is
  the marginal under p(a | beta) proportional to det A^(1/2) is correct.
- Eq. `eq:JT` and `eq:uV`: verified by hand for one pixel and one component (score b c / A with
  variance u V + V; expected curvature u V), and by Monte Carlo for two regions, three components and
  dense noise (`check_sec45_app.out`, block D: Var(s) / (J + T) = 0.995 and 0.998; expected curvature
  0.0045 +- 0.0002 against J = 0.0046 where J + T = 0.039; mean score zero within errors).
- Gap, `sec_model.tex:59-60`: "the marginal keeps the form of Eq. (Lflat) with A -> A + S^-1". This
  drops the factor det S^(-1/2), which is constant only when S is fixed. With a free scale tau it is
  the term sum ln tau that the code carries. Say so.
- Gap, `sec_model.tex:62-63`: "tau_r maximized or marginalized under a prior flat in ln tau_r". The
  marginal likelihood tends to a positive constant as tau -> 0, so a prior flat in ln tau gives an
  improper posterior. Either give the cutoff or drop the marginalized option. This also bears on
  `sec_profile.tex:107-108`, "What the prior offers in exchange is a proper posterior".
- Unstated condition: Eq. `eq:JT` is derived for flat priors on every amplitude. With the CMB prior
  used on the data, A contains the prior precision and the statements hold only on average over CMB
  skies drawn from that prior. This is said in the appendix and nowhere in Sect. 2 or Sect. 6.

### Sect. 3 (the re-solve move)

- Extended target, Eq. `eq:extended`: verified. With proposal density det A(beta')^(1/2)
  exp(-X'/2) the Hastings ratio reduces to the chi-square difference only for a target with
  exp(-X), whose index marginal is the profile posterior and whose amplitude conditional has
  covariance (2A)^-1. X is Gamma(k/2, 1) under it.
- Inconsistency, `sec_legacy.tex:13`: step 3 accepts on the chi-square alone, but D in the next
  subsection contains the prior ratio. One of the two is incomplete.
- Eq. `eq:Tk`: verified. With alpha(D) = E min(1, exp(-(D + Z)/2)), alpha(0) = 1/2 + c_k and
  alpha'(0) = -c_k / 2; the small-step diffusion has d ln pi / d beta = 4 alpha'(0) U' / alpha(0)
  with U = chi2_min / 2 - ln p, hence 1/T = 2 c_k / (1/2 + c_k). The closed form
  c_k = Gamma(k) / (2^k Gamma(k/2) Gamma(k/2 + 1)) = 2^-(k+1) C(k, k/2) follows from a one-line
  integral; I checked it by quadrature for k = 1 to 1000 (ten digits) and reproduced the T column
  and the leading-term column of Table 1 exactly (`check_sec3_temper.out`, block a). T(0) = 1.
- Eq. `eq:Tlarge`: verified. The remainder is sqrt(pi k / 8) / (4k), which is O(k^-1/2).
- Eq. `eq:Tsigma`: verified, including sigma_E = sd(Z)/2 = sd(chi-square at a fresh draw)/sqrt 2.
  Unstated condition: T = 1/2 + 1/(8 f_Z(0)) needs the density of Z to be flat over a few units
  around zero, that is sd(Z) much larger than 1.
- A toy chain of my own (k = 1000, 4000 chains, step 0.1 sqrt T) gives a variance of 20.23 +- 0.11
  against T = 20.32 and an acceptance of 0.511 against 0.513 (block d).
- Eq. `eq:Rn`: verified. I re-derived it: the number of accepted moves per sweep is H_(n+1) - 1,
  the drift of proposal j is proportional to the expected density at the minimum of j draws, and
  the first-order drift does not depend on how earlier moves changed the retained X because each
  step is independent of the history. R_2 = 5/6 exactly. Quadrature gives 0.8333, 0.6465, 0.5428,
  0.4666, 0.3567 for n = 2, 5, 10, 20, 100 (block b), as in the paper.
- Loose statement, `sec_legacy.tex:104-105`: "more index moves narrow the law only as (ln n)^-1/4
  in width". That is the asymptote. Between n = 100 and 1000 the local exponent of R_n in ln n is
  -0.63, so -0.32 in width (my check). Say "asymptotically".
- Acceptance of the half-covariance repair: verified. The move at fixed index is an independence
  sampler with weight exp(-x/2); its stationary acceptance is 2 P(2 G2 < G1) = 2 I_(1/3)(k/2, k/2).
  I get 6.18e-4 at k = 100 and a decay rate of 0.0589 per amplitude (block c).
- Conditions that should be listed next to Eq. `eq:Tk`, because the reader of the abstract will not
  find them: exact conditional draws; the same k amplitudes redrawn by the Gibbs step and entering
  the acceptance; a symmetric proposal of fixed covariance; steps small against the stationary
  width. Outside the last one the variance is within 8 percent of T for k of 1e4 and above and
  within 9.3 percent at k = 100 (stats review, Sect. 5b), which the paper states as "within 8
  percent ... at any step size".
- **Error of description, `sec_legacy.tex:159`:** "The law with exact stored draws predicts 13.1 to
  13.6 at these settings." Eq. `eq:Rn` at k = 190825 and n = 20 gives sqrt(R_20 T_1) = 11.31 and an
  acceptance of 0.132 per proposal (block e). The 13.1 and 13.6 come from a toy chain in two
  dimensions with the recorded proposal widths of 5.6 and 93 profile widths
  (`E035/toy_e014_recorded_proposals.txt`). The paper never tells the reader that the real chain
  proposes two indices jointly with steps that large, so the number cannot be reproduced from the
  equations in the paper.

### Sect. 4 (flat prior and reference frequency)

- Eq. `eq:tilt`: verified by hand and numerically for three components, dense noise, band-averaged
  spectra and two regions (block A: the difference of ln L_flat between references of 30 and 23 GHz
  equals sum n_r beta_r ln(30/23) to ten digits; chi2_min unchanged to 1e-12).
- Wording, `abstract.tex:11-12`: "A flat amplitude prior is exactly an exponential prior on the
  index with rate n ln nu_ref." What is exact is the dependence on the reference frequency: moving
  the reference multiplies the index prior by an exponential. The flat-prior marginal is not the
  profile likelihood times an exponential; the remaining factor, sum_p (1/2) ln sum_b w nu_b^(2 beta),
  is curved (that curvature is the 2T of the next equation). Reword.
- Eq. `eq:response`: verified. The Hessian of (1/2) ln det A is 2T for one component with diagonal
  noise (block B, ratio 2.0000). It is 2.03 T for one component with dense noise, and in my
  two-component example with a free CMB amplitude it is -4.2 T, that is negative.
- Consequence that the paper does not draw: the boxed Eq. `eq:shift`, with its factor
  1 / (1 + 2T/J), holds for one component only. With a free CMB the curvature of the volume term
  can be negative, H can be smaller than J and need not be positive, and the flat-prior marginal
  then has no interior optimum. That is what Table `tab:bpset` shows (optimum on the bound). The
  text applies Eq. `eq:shift` to three-component fits (`sec_flat.tex:163-166`) and should say that
  only the high signal-to-noise limit is used there.
- First-order accuracy (block C, noise-free data, one component): the exact shift is within 1 to 6
  percent of g / (J + 2T) for shifts up to 0.13, and off by 15 to 30 percent for shifts of 0.35 to
  0.74. This agrees with `sec_flat.tex:99-100`.
- `sec_flat.tex:69-71`: "four times as many amplitudes for the same region information is four times
  the shift". True for T much smaller than J and J unchanged. The paper's own data give a factor
  2.5 between nside 8 and 16 (0.062 and 0.158), because J changes too. Say so where the data
  number is quoted.
- Pivot under projected weights (`sec_flat.tex:76-84`): verified. det A factorizes into a term
  independent of beta times m^T W~ m, and g is ln nu_ref minus a generalized mean of ln nu under W~,
  whose off-diagonal elements are negative for a positive spectrum F.
- Arithmetic checked: 128 ln(30/23) = 34.0; ln(30/23) x 128 x 0.059^2 = 0.118; the "predicted"
  column of Table `tab:flat` from the class medians and a pivot of 26.4 GHz (all eight entries);
  the implied pivots 26.04 GHz (twice) and 29.5 to 30.2 GHz from the E022 numbers.

### Sect. 5 (profile error, zero-mean prior)

- Eq. `eq:sandwich`, `eq:sdz`: verified (block D, first-order estimator variance over the sandwich:
  0.995 and 0.998).
- E[J(a_hat)] = J + T and Eq. `eq:debiased`: verified (block D: 0.0392 +- 0.0001 against J + T =
  0.0393 at T/J = 7.6). The claim that the plug-in sandwich is narrower than the curvature width is
  right: (J + 2T) / (J + T)^2 < 1/J. The dotted curve of Fig. `fig:ladder`,
  (1 + t)^(3/2) / (1 + 2t)^(1/2), is the correct expectation.
- Zero-mean prior identity (`sec_profile.tex:86-97`): verified by hand (Woodbury, then the
  substitution x = tau^2 A / (1 + tau^2 A)) and numerically (block E: agreement to 2e-13; same
  optimum to 1e-6; curvature ratio 0.8666 against 1 - n/q = 0.8666).
- Unstated conditions of that identity: one component (b_p and A are scalars); tau maximized, not
  marginalized; and "with its sandwich width" uses that the observed profile curvature is J on
  average. The model paragraph offers marginalization as an equal option, for which the identity
  is approximate.
- `sec_profile.tex:72-74`: "The coupled case with several regions, three components and the dense
  WMAP covariances is tested in Sect. 6." It is not yet. Table `tab:e039` uses "the sandwich error at
  the truth", and the test with errors from each realization is the red pending item
  (`E039/plugin_n8_bpxv.txt` is empty). On the real maps T/J_c is at most 0.013, so the debiased
  sandwich is never exercised on data. The sentence promises a test the paper does not contain.

### Appendix

- Score, its mean, the two variance terms, H_r and the trace identity: verified line by line.
- **Sign error, `app.tex:40-43`.** With the paper's definition s_r = -(1/2) d chi2_min / d beta_r, the
  expected score under a prior variance that is too small by kappa^2 is -(kappa^2 - 1) t, not
  +(kappa^2 - 1) t, for t_r = tr(Lam A^-1 M_r^T B^T) (block F: exact expectation +3.615 at kappa^2 = 4
  where (kappa^2 - 1) t = -3.615). The bias formula in the next line, -(kappa^2 - 1) J^-1 t, is
  right. The slip comes from the module docstring and the stats review, where "score" means the
  gradient of the minimized objective.
- Unstated: A and B in t_r must contain the prior precision. The symbol H is the projector here and
  the curvature in Sect. 4.
- The scalar kappa^2 does not describe the real case. The prior is white, the CMB is not. See
  objection 5.6: per multipole, kappa^2 runs from 200 at l = 2 to 0.35 at l = 16, and is zero for
  every B mode.

## 3. Check of numbers against records

"Match" means the paper's number equals the source after rounding.

| # | Where | Number in paper | Source | Number in source | Verdict |
|---|---|---|---|---|---|
| 1 | abstract, intro | 0.79 k^(1/4); factor 16 at k = 1.9e5 | derivation; E035 | 0.7916; 16.5 | match |
| 2 | abstract:11 | law confirmed "on toy chains to 0.2 percent" | stats review 5a; E035 | exact-kernel grid 0.1 to 0.3 percent below T; chains +-0.3 to 0.4 percent, one 0.46 percent off; E035 itself says 0.7 percent in width | **mismatch (wording)** |
| 3 | Table 1, k = 100 | 6.782, 6.267, 6.769, 6.751 +- 0.024, 0.538 | stats review 5a; my check | same; T and lead reproduced | match |
| 4 | Table 1, k = 2e5 | 280.75, 280.25, 280.31, 280.22 +- 1.01 | stats review 5a; my check | same | match |
| 5 | Table 1 caption | controls 1.000 +- 0.001 | `E035/toy.txt` | 1.000, 0.999, 1.001; 1.000, 1.000, 1.001 | match |
| 6 | Sect. 3.4 | R_n = 0.833, 0.647, 0.543, 0.467, 0.357 | `E035/fig_temper_toy.txt`; my check | 0.8333, 0.6465, 0.5428, 0.4666, 0.3567 | match |
| 7 | Sect. 3.4 | 6e-4 at k = 100; exp(-0.059 k) | stats review 4; my check | 6.18e-4; 0.0589 | match |
| 8 | Sect. 3.5 | 0.918 T; within 8 percent | stats review 5b | 0.918 (k = 1e4), 0.907 (k = 100) | match for k of 1e4 and above |
| 9 | Sect. 3.5 | measured 3.9, 5.7, 3.2, 4.4; toy 4.0 to 4.1, 5.7, 3.0, 4.3 | `E035/toy_at_earlier_settings.txt` | 4.03 to 4.08, 5.66 to 5.72, 3.04, 4.28 | match |
| 10 | Sect. 3.5 | widths 17.0, 16.4; 16.4, 16.6 | E035 addendum | 17.0, 16.4; 16.40, 16.59 | match |
| 11 | Sect. 3.5 | "The law ... predicts 13.1 to 13.6" | `E035/toy_e014_recorded_proposals.txt`; my check of Eq. Rn | toy 13.09 and 13.56; Eq. Rn gives 11.3 | number matches, **description wrong** |
| 12 | Sect. 3.5 | 682 +- 9; 437; 2.44 +- 0.06; 15.9, 16.6; 0.075; 0.074, 0.077; +3.1, +3.0 +- 0.5 | E035 addendum | 682 +- 9; 437; 2.44 +- 0.06; 15.95, 16.57; 0.0751; 0.0735, 0.0765; +3.1, +3.0 +- 0.5 | match |
| 13 | Table 2 | shifts -0.005, +0.007, +0.037; -0.061, +0.062, +0.341; -0.40, +0.34 | E036 result table | same | match |
| 14 | Sect. 4.3 | 0.118 against a measured 0.123 | E036 | 0.118; text says 0.122, table difference 0.123 | trivial mismatch |
| 15 | sec_flat:144 | "regions measured to 0.14" | Table 2 of the paper; E036 | class median 0.161 (0.140 is the noise x 2 arm) | **mismatch** |
| 16 | Table 2, lower blocks | +0.006, +0.167, +0.058, +0.859; -0.060, +0.063; factor 1.63 | E036 | same | match |
| 17 | Sect. 4.3 | +0.271, +0.172; 1.913, 1.216; 26.04 GHz | E022 table; E036 lines 419-420 | same | match |
| 18 | Sect. 4.2 | 1.01 and 1.03; 5 and 16 percent | E039 review item 1 | 1.008, 1.031; 0.955, 0.842 | match |
| 19 | Sect. 5.2 | 0.054; 1.59 at 1.56; 49 percent; 38 percent; 0.95 to 1.05; 0.635 | E041; `E041/ladder.txt` | 0.054; 1.588 at 1.555; 0.485; 0.3825; 0.946 to 1.052; 0.635 | match |
| 20 | sec_profile:101-104 | zero-mean bias "between -0.6 and +0.3"; coverage 0.59 to 0.73 | `E041/ladder.txt` | measured cells: -0.58 to +0.17 (the +0.27 is in a cell with sd 1.15, not measured); 0.5875 to 0.725 | **mismatch on +0.3** |
| 21 | Table 3 | all 44 entries | E039 result table | same (0.0680 is 0.0687 - 0.0008 = 0.0679) | match |
| 22 | Sect. 6.2 | 1.0006, 0.022; +0.03, 0.99, 0.998 +- 0.009 | E039 | 1.0006, 0.022; +0.03, 0.989, 0.9977, 0.0091 | match |
| 23 | Sect. 6.1 | 1e-9; 2.0 to 2.2 uK; 65 to 164 | E038 | 8.8e-10; 2.01 to 2.24; 65, 92, 101, 164 | match |
| 24 | sec_data:20-22 | chi2 1.77, 1.40, 1.21; excess "rising to between 3 and 7 at multipoles 2 to 6" | E038 (b) | 1.771, 1.401, 1.209; 3.3 to 6.9 (30 GHz), 2.1 to 3.7 (44), 1.2 to 3.7 (70); record text "2 to 7" | **mismatch on "3 and 7"** |
| 25 | Sect. 6.1 | (1.46, 1.16, 1.27), (1.28, 0.45, 1.25), (1.17, 0.16, 1.69); 22 sub-tests, none beyond 3 sigma | E038 | same | match |
| 26 | sec_data:35-36 | effective frequencies 22.44, 32.70, 40.34, 60.03 | headers of `E040/bpxv_n16_wm.txt` and every other fit | 22.41, 32.66, 40.30, 59.99 (drift-corrected) | **mismatch** |
| 27 | Sect. 6.1 | 28.09, 43.80, 69.41; factor 1.12 | `E040/s2_noK_lfi_wm.txt`; E040 amendment | 28.09, 43.80, 69.41; 1.121 | match |
| 28 | sec_data:45-47 | bias 0.005 to 0.011 and 0.009 to 0.020, "an upper estimate" | E038 (c), (c2) | 0.0048 to 0.0111; 0.0091 to 0.0196; "the upper-estimate caveat of the first results does not apply" | numbers match, **caveat contradicted** |
| 29 | Table 4 | all rows; T/J_c 0.013, 0.004, 0.0003, 0.002 | `E040/bpxv_n16_wm.txt` | -3.0940, -3.0965, -3.0867, -2.9950; 0.0121, 0.0149, 0.0088, 0.0120; shifts and predictions as printed; 0.0129, 0.0044, 0.0003, 0.0022 | match |
| 30 | Sect. 6.3 | 0.92 to 0.97; twelve sd; 0.062 at nside 8 | E040 | 0.92, 0.95, 0.95, 0.97; 12.4; +0.0618 | match (but see below) |
| 31 | Sect. 6.3 | window 0.007, 0.003, 0.009, 0.004 twice; residual 0.002 | E040; E038 (c2) | 0.0067, 0.0025, 0.0085, 0.0035; 0.0070, 0.0030, 0.0088, 0.0041; 0.0023 | match |
| 32 | Sect. 6.3 | nside 8: 0.23 flatter, 0.07 steeper; +0.22, -0.02 | E040 R8 | +0.232, -0.067; +0.217, -0.019 | match |
| 33 | Table 5, all scales | five rows | E040 stage 2; `E042/pr3_arm.txt` | same | match |
| 34 | Table 5, low-pass | four rows; Eq. result | `E042/lowpass20.txt` | -3.1132 +- 0.0144, -3.1120 +- 0.0166, -3.1401 +- 0.0109, -3.0613 +- 0.0136, and the other three rows | match |
| 35 | Sect. 6.4 | +0.08, +0.16, +0.11, +0.37; 3 to 14 percent; 0.32 to 0.33, 0.35 | E042 control and ratio table | +0.083, +0.162, +0.109, +0.365; 3, 6, 4, 14; 0.321 to 0.332, 0.353 | match |
| 36 | sec_data:254 | "errors grow by 10 to 20 percent only" | `lowpass20.txt` against `bpxv_n16_wm.txt` | 19, 11, 24, 13 percent | small mismatch |
| 37 | Sect. 6.4 | 1.9, 2.0, 3.6, 1.7 | my arithmetic from `lowpass20.txt` | 1.90, 1.98, 3.60, 1.65 | match |
| 38 | Sect. 6.5 | 0.50; 11 percent; 0.16, 0.44, 0.18; chi2 6 against 2.6 | `lowpass20.txt`; `pr3_arm.txt`; E040 | 0.499; 11.2; 0.159, 0.445, 0.174; 5.97, 2.64 | match |
| 39 | sec_data:305-307 | 30 GHz fainter by 11, 8, 2, 2 percent; "Ka control is brighter ... by 1 to 2 percent" | `E042/predict30.txt`, `predictKa.txt` | -0.114, -0.077, -0.024, -0.024; control +0.019, +0.014, +0.037, +0.010 | **mismatch: control is 1 to 4 percent** |
| 40 | Sect. 6.5 | 16 to 43 times; 2 to 4 for the control | same files | 33/43, 16/19; 3.7/2.2, 3.5/2.3 | match |
| 41 | Sect. 6.5 (i) | 20 modes; 0.03; less than 0.02 | E042 (c) | 20; +0.034; at most 0.015 | match |
| 42 | sec_data:317-318 | PR4 "7 to 15 percent fainter" with slopes 0.85 and 0.96 | E042 exploratory tables | slopes 0.854, 0.841, 0.955, so 4.5 to 16 percent; cross spectra 7 to 10 percent | **internally inconsistent** |
| 43 | Table 6 | all entries | `E040/s3_bpset_cmbflat_wm.txt` via E040 | -2.695 +- 0.040, -2.921 +- 0.050, -2.819 +- 0.025, -2.780 +- 0.042; -4.036, -3.432, -2.885, -2.997; -1.600, -2.314, -2.752, -2.479 | match |
| 44 | Sect. 7 | -0.026, +0.029, 0.000; 600 replicates; 1.11 to 1.42; -3.017 | E036; E016; E032 | same; 300 + 300; 1.109 to 1.418; -3.0160 to -3.0184 | match |
| 45 | sec_limits:24 | dust index shifted "by 0.11 at four stated standard deviations" | E021 | +0.114 at 21 sd (exact scan); "+0.13 at 4 posterior sd" (follow-ups) | **mismatch: two record numbers merged** |
| 46 | Sect. 8 | 0.03, 0.94, 0.97; -0.06 +- 0.04, +0.06 +- 0.05, 0.96, 0.94; 180 | registry E006; E014 | same | match |
| 47 | sec_concl:19-20 | published errors "four to five times ours" | quoted values over `lowpass20.txt` errors | 0.06/0.0144 = 4.2; 0.04/0.0136 = 2.9; 0.06/0.0166 = 3.6 | **mismatch: three to four** |
| 48 | sec_concl:17-19 | "agree with ... Fuskeland et al. (2014) ... within their errors" | citation report 1, item 11a | plane -2.98 +- 0.01, high latitude -3.12 +- 0.04; paper's Plane -3.06 +- 0.01 is 4.7 combined sd away | **mismatch (substantive)** |

**Count.** 12 mismatches in numbers (rows 2, 14, 15, 20, 24, 26, 36, 39, 42, 45, 47, 48). Ten are
small. Two matter: row 26 (the frequencies printed are not the ones the fits used, and the drift
correction that makes the difference is not mentioned in the paper) and row 48.

**Numbers I could not trace to an experiment record.**
- "Delta bands at the nominal WMAP frequencies would steepen a K to Ka index by 0.055"
  (`sec_data.tex:37-38`). The delta-band arm is in the red pending list. I get about 0.057 from the
  ratio of logarithms, so it is plausible, but it is a calculation and not a result.
- "a Cholesky factorization takes 9 s and one profile evaluation with its gradient 12 s on ten
  cores" (`sec_impl.tex:26`).
- "published LFI colour-correction ratios to 3e-4" (`sec_impl.tex:30`): found only in a docstring
  of `lowres/spectra.py`, not in a record.
- The class medians 0.019 and 0.161 of Table 2: not found in the part of E036 I read. Not searched
  exhaustively.
- The values of Fuskeland et al. (2014) that the conclusions lean on are not given in the paper.
- "whose processing masks differ from band to band (Bennett et al. 2013)" (`sec_data.tex:188-189`):
  neither citation report verifies this sentence.

**Places where the paper's wording is stronger than, or contradicts, the record.**
1. `sec_legacy.tex:159`: "The law ... predicts 13.1 to 13.6" (it is a toy at recorded step sizes; the
   law gives 11.3).
2. `sec_legacy.tex:159-160`: "The difference has a specific cause." The record says "closed at the
   level of analysis", with the forward run not done and the re-implementation not executing
   Commander4. The record also has the synchrotron shift, predicted -0.1 and measured -2.7 +- 2.0,
   which the paper leaves out while quoting the dust shift.
3. `sec_legacy.tex:190-191`: "Commander3's index samplers are of the latter kinds". E009 item 4 and
   finding F10 say that the 2024 Commander3 module (`sample_specind_mh` and siblings) "redraws
   amplitudes by CG and accepts on chi^2, like Commander4", and E009 item 1 says the BeyondPlanck
   polarized sampler holds the amplitudes fixed, which is neither of "the latter kinds". The
   sentence is wrong on both counts by the project's own record.
4. `sec_flat.tex:166-169`: "with the CMB constrained ... which may be why the index posteriors of
   Svalheim et al. show no such pull". E022 reads the released BeyondPlanck parameter file as a flat
   CMB prior (`COMP_CL_TYPE01 = none`), and Sect. 6.6 of the paper is built on "an unconstrained
   CMB amplitude" in the same analysis.
5. `sec_data.tex:46`: "an upper estimate" (row 28).
6. `sec_data.tex:57-58`: "The expectations for every run below were written down before the first
   fit to real data." The first real-data fit started at 22:04 on 2026-11-10. The expectations for
   stage 3 (22:50), for E042 (22:56, 23:00, 00:30) and for the low-pass arm (01:17 on 2026-11-11)
   were written before their own fits and after earlier real-data results had been read. Tables 5
   (lower block) and 6 and all of Sect. 6.5 fall under this.
7. `sec_data.tex:243`: "The cause is aliasing." The record says "Mechanism, consistent with the last
   bin of the ratio table".
8. `sec_data.tex:168-170` uses the nside-8 real-map number 0.062 as evidence, sixteen lines before
   "We discard the nside-8 fits of the real maps" (`sec_data.tex:186`).
9. `sec_data.tex:340-343`: "Both report that the high-latitude index, fitted freely, goes to
   'unreasonably flat' values". See Sect. 4 below: the quotation is verified in Watts et al. only.
10. `sec_data.tex:380`: "the profile likelihood on the same maps does not show it". The profile
    optimum of that band set is -2.70 at high latitude, 0.4 flatter than the WMAP value, which the
    next sentence admits.
11. `sec_limits.tex:29-30`: "the field mean is recovered (-3.017 against -3.000)". E032 gives a
    posterior sd of 0.0031, so the offset is five posterior standard deviations.
12. `sec_impl.tex:11-13`: the E006 run "matched the exact profile posterior within 0.03". The
    registry adds that this target sat 1.9 and 3.7 sd from the input because of the noise treatment
    in the harness, and that the flat-prior marginal chain was stuck in a local mode. By working
    principle 4 the paper should report both halves.

## 4. Claims of novelty and of prior work

Read levels are those of the citation reports. Several papers were read only through an HTML
summarizer (Svalheim et al., Watts et al., NPIPE); the reports themselves ask the first author to
check those quotes against the PDFs.

- **Andrieu and Roberts (2009).** Represented correctly: MCWM refreshes both estimates, "pi(theta) is
  typically not the invariant distribution", existence of an invariant law "is not obvious in
  general". The claim "Here it can be written down" (`sec_legacy.tex:49`) is acceptable only with
  the words "in the small-step limit" attached at that point; the exact kernel is not reversible
  for any power of the target, as the paper says later. The prior-art note found no explicit law
  in Alquier et al. (2016), Medina-Aguayo et al. (2016), Nicholls et al. (2012) or Ceperley and
  Dewing (1999), mostly at summary level. None of these four is cited. Nicholls et al. state that
  the naive algorithm is "slightly over-dispersed" and cite Ball et al.; Ceperley and Dewing give
  the exact correction for Gaussian noise. A referee from statistics will ask for them.
- **Ball, Fink and Bowler (2003).** Confirmed: T ~= sqrt(pi/8) sigma for the sign rule with Gaussian
  errors, approximate. The equation number was checked in arXiv v1, not in the journal.
- **Bardenet, Doucet and Holmes (2017).** Overstated. The paper says they "recover" Ball's
  temperature "without assuming Gaussian noise". The report says: a Berry-Esseen bound, a logistic
  fit to the normal distribution function, a Barker target pi^beta with beta = 1.702 sqrt(t) / (n
  sigma), which is a temperature of 0.588 s against Ball's 0.627 s, in the subsampling setting, and
  described by its authors as "impractical". Same proportionality, different constant. Write "obtain
  a temperature proportional to the noise, with a slightly different constant".
- **Eriksen et al. (2008).** Confirmed (offset -1.8 +- 0.4, "more than 4 sigma", Sect. 4.2).
- **Dunkley et al. (2009).** Confirmed on every point the paper uses, including the pivot at 23 GHz,
  the centre at -2.6, the growth with the number of pixels, the Jeffreys remark and the division by
  L_0. The high and low latitude values in the conclusions are correct.
- **Alonso et al. (2017).** Confirmed with the caption value ("more than 3 sigma"; the text says 4).
  "which makes their target the profile likelihood" (`sec_flat.tex:11-12`) is the project's own
  inference; the report says "The paper does not say so". Mark it as an inference. Alonso et al.
  also say their volume prior is "equivalent to the widely used Jeffreys prior for the spectral
  parameters" up to a logarithmic term, which is prior art for the paper's remark that the profile
  is a marginal under a Jeffreys-type prior (`sec_model.tex:55-57`).
- **Stompor et al. (2009).** Confirmed (Eq. 9 and Eq. 10). They also write that the profile
  curvature gives "useful, if not strictly precise, estimates" of the uncertainty. Sect. 5 should
  cite that sentence: the paper quantifies a caveat they stated.
- **Weiland et al. (2022).** The two numbers are right (-3.27 and -2.90, Fan region, PR3 30 GHz,
  statistical errors). Not said in the paper: with the same PR3 map the paper's own fit puts the K
  to 30 GHz index on the other side of K to Ka in every region (-2.71 against -2.98 in the Plane),
  the opposite sign to Weiland et al. A cited result is contradicted without comment.
- **Svalheim et al. (2023).** Three problems. (a) The report could not read the PDF and found no
  "unreasonably flat" passage in it. What it found in Sect. 4.4: priors of -3.3 +- 0.1; "we only
  include the prior for the high-Galactic and Galactic Center regions"; "very strong degeneracies
  with respect to instrumental systematic effects"; "relative gain uncertainties and degeneracies
  with the CMB quadrupole"; a wider prior "permits synchrotron indices as shallow as beta_s = -2.7".
  So in the published chains the high-latitude index is a prior draw, not a sampled quantity (E010
  says the same from the parameter file). The paper's "sample the polarized synchrotron index with
  a Metropolis step conditional on the amplitudes" (`sec_data.tex:338-340`) is true of the Spur
  and the Plane only. (b) The 30 GHz reference rests on a figure caption and on a table of a
  companion paper. (c) The published paper is BeyondPlanck XIV; the records call it XV.
- **Watts et al. (2023).** The quotation is theirs (Sect. 2.4 of arXiv v1): "when fitting both the
  polarized synchrotron amplitude and calibration freely without priors, the synchrotron spectral
  index at high Galactic latitudes tend to be biased toward unreasonably flat values, beta_s <~
  -2.5, which was likely due to a low level of unmodeled systematics, for instance
  temperature-to-polarization leakage". The paper drops "and calibration", which is a degeneracy
  that a maps-only test cannot address, and attributes the sentence to both papers. On K band
  against 30 GHz the paper cites Watts et al. for the poorly measured K-band modes and omits their
  headline claim that the discrepancy is "finally resolved" by joint reprocessing and "consistent
  with instrumental noise at high Galactic latitudes". It also omits, where it matters most, their
  statement that the nine-year covariance treatment "does not provide a complete description of the
  full large-scale uncertainties": the paper's adopted estimator takes the index from multipoles 2
  to 20 of exactly those maps with exactly those matrices.
- **Huber (1967), White (1982).** Confirmed as the origin of the form A^-1 B A^-1. Both are for
  independent identically distributed data and fixed dimension, which is not this setting. The
  sentence that carries the weight, "This is the errors-in-variables model with known error
  covariance, for which the consistency of the slope is classical; the classical caveat is that the
  variance of the estimator is not the inverse Fisher information" (`sec_profile.tex:8-11`), has no
  citation at all. The report verified it only in Engel (1985), citing Patefield (1977, 1978);
  Fuller (1987) could not be accessed.
- **Novelty of Sect. 5.** (From memory, to be checked in Fuller 1987, Sects. 1.3 and 2.2, and Gleser
  1981.) The variance of the slope in the functional errors-in-variables model with known error
  variances has a second term that is the product of the two error variances, and the standard
  moment estimator subtracts the error variance from the observed sum of squares. These are the
  one-band-pair versions of T and of J_c = J(a_hat) - T. If so, the "trap" and the "debiased
  sandwich" are the multi-band form of textbook results and the abstract should not present item
  (3) as new. What may be new is the closed form for dense noise and several components, and the
  coverage ladder.
- **Novelty of Sect. 4.** "What has not been written down, as far as we can find, is how large the
  pull is in closed form and what sets its sign" is hedged properly. The project's novelty note
  lists as "referee must-cite" for this claim: Planck 2015 X ("we adopt a Jeffreys prior ... to
  suppress prior volume effects"), Dickinson et al. (2009), Armitage-Caplan et al. (2012), Jew and
  Grumitt (2019), Gratton (2008), and BeyondPlanck XIII (Andersen et al.), which describes the
  marginal and ridge index samplers of Commander with the remark that the two differ by "the last
  determinant term". None is cited. The omission of BeyondPlanck XIII is the serious one: it is the
  published description of the very options the paper recommends in "What to do" item (i).
- **Bennett et al. (2013), NPIPE, Planck 2018 XI.** The report's flags 3, 4, 5 and 10 against the
  draft are still open: the dense matrices are described in Jarosik et al. (2007, 2011) and on
  LAMBDA, not in Bennett et al.; NPIPE says "half-ring"; the PR4 353 GHz maps appear to include the
  spider-web bolometers, so "the band of the polarization-sensitive bolometers" is the PR3
  description (1.137 against 1.121, a 1.4 percent effect on the dust extrapolation); the reference
  list lacks arXiv numbers and end pages.
- **A question the report raised and the draft does not answer** (round 2, gap 6): NPIPE Sect. 4.2
  describes released nside-16 noise covariance matrices for PR4. The paper builds an empirical
  model from one half-ring difference and does not say why the released matrices are not used.

## 5. Objections a hostile referee would raise, and what would answer each

**5.1 The quoted errors do not mean what they say when the model misfits.** The sandwich is the
variance of the estimator under the model. The model fails: chi-square per element and band of 1.82
in the Galactic Center and 0.78 against 0.6 in the Plane; K with Ka against K to V differ by 1.9,
2.0, 3.6 and 1.7 nested sd with one sign (chi-square 23.5 for four numbers, my arithmetic); nside 8
against 16 differ by up to 21 sd, of which the degrade explains the Galactic Center and one third
of the Plane (-0.067 measured, -0.019 simulated); all scales against L = 20 against L = 12 drift by
0.02 to 0.09. Under misfit the fitted number is a noise-weighted average of whatever the sky does,
and its uncertainty is set by which weighting was chosen, not by the noise. The paper says so in
"What we did not do" (`sec_concl.tex:44-45`) and then prints +- 0.01 in the abstract. The Galactic
Center is said to have an error that "means little" (`sec_data.tex:160`) and is given as -3.14 +-
0.01 in Eq. `eq:result`. The goodness of fit of the adopted low-pass fit is not reported at all:
`lowpass20.txt` prints chi-square 3036.0; counting degrees of freedom as E040 does I get about 2618
(my estimate), a ratio of 1.16, which is five to six standard deviations high; for K with Ka 1143.6
over about 870.
*What would answer it:* errors from a block jackknife or bootstrap over nside-2 or nside-4
super-pixels inside each region, which carry the misfit; chi-square and degrees of freedom per
region for the adopted fit; no error bar for a region whose fit is rejected.

**5.2 The low-pass baseline was chosen after seeing the data, and its own acceptance test has
failed.** The sequence in E042: a negative control, read at 01:17, failed; a first low-pass version
was ill-posed; a second version with an anchor band ran from 01:40 to 01:56; its output at L = 20 was
adopted as "the values we quote" in a draft compiled at 02:06, with the stability test in red.
The direction and rough size of the change were known before the expectations Z1 to Z4 were
written, because the high-pass control had just been read. The record's own rule: "if they do, L =
20 becomes the baseline ... If they do not, there is no resolution-independent index for these
regions at this depth and the paper says so." `lowpass12.txt` (02:25) gives -3.142 +- 0.018,
-3.108 +- 0.018, -3.180 +- 0.013, -3.081 +- 0.015. Against L = 20 the Galactic Center moves by 3.1
sd of the L = 12 error, so Z4 fails as written; by the nested measure that the paper itself uses
for band subsets the moves are 2.8, 0.6, 5.7 and 3.8 sd. A drift that continues from 9 degrees to
15 degrees is hard to assign to aliasing at 3.7 degree pixels. Nor has the estimator ever been run
on a simulation with known truth: E039 validates the all-scale fit, at nside 8 with 200
realizations and at the working nside 16 with 12. Working principle 3 of this project forbids real
data before that.
*What would answer it:* the E039 harness with the low-pass estimator at nside 16, at least 100
realizations, at L = 12, 20, 30; the three real-data values side by side with nested differences;
and then either a plateau, or the sentence the record prescribes.

**5.3 "Beam aliasing" is asserted, not shown, and the evidence in hand leans against it.**
(a) E038 (c2) simulated a sky with structure to multipole 1024 through each band's beam and found a
residual of 0.0002 to 0.0023 after window matching, ten to thirty times smaller than the 0.02 to
0.07 now attributed to aliasing. The paper's reply, that this sky "had no compact structure", is a
conjecture. (b) The window operator uses Gaussian beams; the record itself notes that "K band has
less power at the pixel scale than the Gaussian 52.8 arcmin window gives it". The WMAP beams are
not Gaussian, and an error in the K-band window gives the same signature with no aliasing. (c) My
estimate with Gaussian windows: power aliased from multipole 100 into the last bin raises Ka over K
by 8 percent, from multipole 70 by 3 percent; the measured rise in that bin is 8 percent, which
would need nearly all the power at multipoles 35 to 47 to come from multipole 100 and above. (d)
Fuskeland et al. (2014) smooth K and Ka to a common 1 degree beam at nside 64, where this aliasing
does not exist, and find -2.98 +- 0.01 in the plane. That is the paper's all-scale K, Ka value
(-2.980 +- 0.016), not its low-pass value (-3.04). The regions and the method differ, but the
simplest reading is that the flatter index at small scales in bright regions is on the sky, and
that "index from multipoles 2 to 20" measures a different quantity rather than a debiased one.
(e) The paper raises a third mechanism itself, sub-pixel weights that differ between bands in the
released maps (`sec_data.tex:186-190`), and does not test it.
*What would answer it:* a noise-free simulation with a high-resolution template that has compact
structure, the true WMAP beam transfer functions and the WMAP degrade, showing a bias of this sign
and size; the predicted ordering of the effect with beam width (Ka, Q, V against K), which the
ratio table can already test; the same analysis on maps smoothed to a common beam before
degrading. Until then the abstract should say "the index depends on angular scale by 0.02 to 0.09"
and not name a cause.

**5.4 The reading of the BeyondPlanck and Cosmoglobe flat index overreaches.**
(a) Internal contradiction. `sec_flat.tex:166-169` explains the absence of a pull in the published
posteriors by a constrained CMB at a 30 GHz pivot. Sect. 6.6 explains the presence of a pull in the
same analysis by an unconstrained CMB. E022 says the released configuration has no CMB prior. One
of the two passages has to go.
(b) The test is not finished by its own rule. E040 stage 3: "If B1 and B2 hold, the paper says
...". B2 is red in the draft. The abstract and the conclusions already say it.
(c) The mechanism predicts too much. On the paper's maps the flat-prior optimum at a 30 GHz
reference is -2.31 in the Spur and -2.48 in the Plane (Table 6), and the shift grows with the
number of amplitudes, which is far larger in the published analyses (E022, in a surrogate with LFI
and 353 GHz only: at nside 64 the Spur is pinned at the edge). The published Spur and Plane, the two regions that were
actually sampled with a likelihood, are -3.17 +- 0.06 and -3.03 +- 0.07. A reading that explains
the high latitudes has to say why these two were spared: the prior of width 0.1, a chain that does
not mix (E014 and E022 both show conditional steps far narrower than the marginal), smoothing
priors on the amplitudes (E022 notes that the published synchrotron and dust amplitudes carry
them, so the amplitude priors were not flat), or something else. The paper is silent.
(d) The high-latitude index of the published chains is a prior draw (Sect. 4 above). The "fitted
freely" behaviour is known from one sentence in Watts et al. about earlier runs "fitting both the
polarized synchrotron amplitude and calibration freely". Calibration is outside a maps-only test.
(e) The profile optimum of the same band set is itself flat: -2.70 +- 0.04 with a free CMB, -2.99
with the prior, against -3.11 from WMAP. The paper concludes that "these seven maps are not
consistent with one power law, a modified blackbody and a CMB of the expected size". That is an
unmodelled systematic or model failure by any name, in the configuration under discussion. The
conclusion "offers a reading that does not need unmodelled systematics" (`sec_concl.tex:33-34`)
contradicts the paragraph it summarizes.
(f) Sect. 3.5 misstates what Commander3 does (Sect. 3 above, item 3).
*What would answer it:* run B2; run stage 3 with the N(-3.3, 0.1^2) prior and report all four
regions; use E022, which has the published maps at nside 16 and already found the high-latitude
flat-prior marginal at the grid edge with a flat CMB, as the primary evidence instead of PR4; state
(c) and (d) in the text; reduce the claim to "the volume term of this configuration has the sign of
the reported behaviour and is large", and take it out of the abstract unless (c) is resolved. The
first author has access to the chains and to the people; "can only be settled with them" is an
invitation a referee will hand back.

**5.5 The PR4 and PR3 statements are stated above their evidence.**
(a) "its index against K band differing by 0.5 between the 2018 and PR4 releases"
(`abstract.tex:33`) compares -2.92 and -3.43 at high latitude from all-scale fits with reduced
chi-square 5.97 and 2.64. The paper holds that all-scale fits are biased, and E042 says of these
two: "The like-for-like low-pass comparison is the arm with the released PR3 covariance", which is
pending. (b) The same subsection says that at high latitude "the two 30 GHz maps agree with each
other to 1 to 3 percent in amplitude per multipole bin" (`sec_data.tex:324`). A 0.5 difference in
index over 0.22 in ln nu is 11 percent in amplitude. Both cannot describe the same pair of maps as
a spectral statement; the record's explanation is that the PR3 fit is driven by large-scale
pattern, "not an amplitude ratio at all". Then it is not an index, and the abstract should not call
it one. Where the two diagnostics agree is the plane (PR4 fainter than PR3 by 7 to 10 percent in
cross spectra, 0.33 in the fitted index). (c) Three measures of the high-latitude deficit of PR4
are in the records and do not agree: 11 percent from the fit, 11.4 +- 0.3 percent from a regression
on a noisy regressor, and about 5 percent at multipoles 2 to 10 "a 3 sd hint" with consistency
above, from cross spectra. The paper quotes the first two as results and the third as a caveat.
(d) The sign against Weiland et al. (2022) for PR3 is opposite and not discussed. (e) The leakage
reading is labelled exploratory, correctly, and then enters the argument that the maps "cannot be
combined". (f) "the two most sensitive synchrotron channels there are" (`sec_data.tex:276-277`)
ignores S-PASS, C-BASS and QUIJOTE. (g) "Sky signal assigned to the wrong index would be mostly E"
(`sec_data.tex:310`): published synchrotron B to E power ratios are roughly 0.35 to 0.5 (from
memory), so "mostly" is too strong; the argument survives as "E above B", which is what the Ka
control shows.
*What would answer it:* the pending arms (released PR3 covariance at matched multipoles, amplitude
factor, angle, DR1); one agreed measure of the deficit with its error; the conclusion scoped to
"the official PR3 and PR4 30 GHz maps with the official nine-year K band", with the Cosmoglobe DR1
claim of resolution stated. The safe statement, already supported, is that K with 30 GHz fits one
power law with reduced chi-square 2.6 to 6 and leaves 16 to 43 times the noise power at multipoles
2 to 6 in both E and B. That is enough to exclude 30 GHz from the index and needs no number for the
release difference.

**5.6 The treatment of the CMB prior.** (a) On the data the "profile likelihood" is not the
prior-free estimator of Sect. 5: the CMB carries a white Gaussian prior of 0.2 uK per element. Its
unbiasedness is then conditional on that prior being right, as the appendix says. (b) It is not
right in shape. A white prior of 0.2 uK at nside 16 is C_l = 1.6e-4 uK^2 sr for E and for B at
every multipole. Against the LCDM EE spectrum with the working window (`check_cmb_prior_shape.out`,
CAMB, tau = 0.054) the ratio of true to prior variance is 199, 128, 66, 29 and 12 at l = 2 to 6,
1.0 at l = 10, 0.35 to 0.67 from l = 12 to 47, and zero for all B modes. Weighted by mode count
over the E modes it is 7.4 for l = 2 to 20 and 19 for l = 2 to 12. The low-pass estimator moves the
index information to exactly these multipoles. (c) The paper prints no sensitivity for the fit it
quotes. E040 has +0.15 sd per unit of variance mismatch at high latitude for the all-scale fit;
`e042e_lowpass.py` does not print the quantity. With mismatches of order 3 to 9 on average the
effect could be of the order of the statistical error, and the scalar formula of the appendix does
not apply when the mismatch runs over three decades. (d) Conclusion item (iv), "state the
sensitivity of the index to that prior's variance", is not followed by the paper for its own
headline. (e) For the band set without K the profile optimum moves by 0.30 (seven sd) between a
flat and a 0.2 uK prior; the paper reads this as the 30 GHz problem. A prior that suppresses real
large-scale E modes by two orders of magnitude is an alternative that has to be excluded first.
*What would answer it:* the LCDM pixel covariance as the CMB prior (everything is dense already),
with the flat-prior fit beside it; both for all scales and for each L.

**5.7 The empirical LFI noise model.** Three parameters, isotropic in the whitened frame, fitted to
one half-ring difference and then tested on the same difference. It cannot see anything common to
the two halves, which for PR4 includes the calibration and the template corrections. The 30, 44
and 70 GHz maps are not independent on large scales (E040 notes that PR4 calibrates 44 and 70 GHz
with a polarized prior built on 30 GHz). The released PR4 nside-16 covariance matrices are not
used and the paper does not say why. The noise-scale arms announced in E038 ("white covariance
times r, and the model with A doubled") are not in the draft. Every LFI number in Tables 5 and 6
and all of Sect. 6.5 depend on this model, and "approximate in its error" does not cover the
optimum, which also depends on the relative weight of large and small scales.
*What would answer it:* the released matrices, or the arms; a test of the model on a split it was
not fitted to (detector or year split).

**5.8 Validation at the wrong resolution and of the wrong quantity.** The 200-realization
validation is at nside 8, where the simulated sky is constant inside the coarse pixels by
construction; the real analysis is at nside 16 with 12 realizations (Galactic Center mean z of
-0.54 +- 0.27). The error that is validated uses J and T at the truth. The error that is quoted
uses fitted amplitudes, and its coupled test is pending.

**5.9 The headline uses the large scales of maps whose large-scale covariance is disputed.** Watts
et al. (2023) say the nine-year matrices do not fully describe the poorly measured modes. The
paper's test (i) removes the 20 modes the K-band matrix flags; the K, Ka index moves by up to
0.015, which is as large as the quoted error, and the test was done on the all-scale fit, not on
the adopted one. The same projection for Ka, Q and V is not done.

**5.10 Is Sect. 3 about anything that was published?** The paper analyses a move in "the
development version of Commander4" and says Commander3 is different. If no published number used
the move, a referee will ask why it leads the paper. The honest position is in E009: Commander3's
2024 module has the same structure and "which production runs (likely DR2-era) used them with CG
updates is unknown". Say that, or find out. The real-map confirmation is a fit after the fact: one
measured factor (2.44) taken from the log of the chain it then explains. The paper admits this in
the conclusions; Sect. 3.5 and the abstract ("use it to account for a real-map chain") do not.

**5.11 Title, abstract and conclusions against the body.** The title promises "the polarized
synchrotron index of the public WMAP and Planck maps". The body delivers a WMAP K to V index that
is not stable under its own cut, and a statement that Planck 30 GHz cannot be used. Of the last
ten lines of the abstract, the estimator comparison is solid; the band set reading, the aliasing
bias, the release difference and the three quoted indices each depend on a red pending item (B2;
L stability; released PR3 covariance; dust, CMB-prior and start-value arms). "(statistical)" after
+- 0.01 does not make it safe when the same paper shows shifts of 0.02 to 0.09 from analysis
choices. Conclusion item (v), "model the difference between the beams", is not what the paper did;
it excluded scales.

**5.12 One paper or two.** Two. The first: Sects. 2 to 5 with Tables 3 and 4 as the demonstration
on real maps that the choice of estimator and of reference frequency moves the answer by twelve
standard deviations. That paper is nearly ready, is not exposed to objections 5.1 to 5.9, and
should claim nothing about the true index. The second: the index itself, with the scale
dependence, the CMB prior, the dust model, the 30 GHz anatomy and the BeyondPlanck configuration,
each with simulations and the pending arms. Every real-data result of the present Sect. 6 was
produced in the four hours before the draft was compiled (first fit 22:04, PDF 02:06), and its
second half reads as a lab notebook: results, a failed control, a remedy and an exploratory aside
in sequence. Kept together, the weak half will set the referee's
opinion of the strong half.

## 6. Required before submission, and suggested changes

### Required

1. **Major.** Withdraw the quoted indices of Eq. `eq:result` from the abstract and conclusions, or
   replace them by what the three cuts support: report all scales, L = 30, 20, 12 side by side with
   nested differences, apply rule Z4 as written (it fails at L = 12), and give a systematic range
   that covers the drift. (5.1, 5.2)
2. **Major.** Validate the low-pass anchor estimator on simulations with known truth at nside 16
   before quoting any real-data number from it. (5.2, 5.8)
3. **Major.** Either demonstrate the aliasing mechanism by simulation with compact structure and the
   true WMAP beams, or remove "beam-aliasing" and "The cause is aliasing" and describe a scale
   dependence of unknown origin. Reconcile with Fuskeland et al. (2014), and correct the sentence
   that claims agreement with them. (5.3; Sect. 3 row 48)
4. **Major.** Repair the BeyondPlanck and Cosmoglobe reading: remove the contradiction between
   `sec_flat.tex:166-169` and Sect. 6.6; run B2 and report it; address why the Spur and Plane were
   not pulled; say that the published high-latitude index is a prior draw and that the "fitted
   freely" sentence is from Watts et al. and includes calibration; delete "does not need unmodelled
   systematics". If these cannot be done, remove the reading from the abstract. (5.4)
5. **Major.** Correct `sec_legacy.tex:188-191` on Commander3 to agree with E009 and F10. (Sect. 3,
   item 3; 5.10)
6. **Major.** Remove "differing by 0.5 between the 2018 and PR4 releases" from the abstract, or
   support it with the like-for-like arm and reconcile it with the 1 to 3 percent agreement at
   high latitude and with the sign of Weiland et al. Scope the 30 GHz conclusion to the official
   products. (5.5)
7. **Major.** Replace the white CMB prior by the LCDM covariance or a flat prior for every fit that
   is quoted, show the other beside it, and print the sensitivity for the quoted fit. (5.6)
8. **Major.** Give errors that survive misfit (block jackknife or bootstrap), the chi-square and
   degrees of freedom of the adopted fit per region, and no error bar for the Galactic Center.
   (5.1)
9. **Major.** State why the released PR4 covariance matrices are not used, or use them; include the
   noise-scale arms. (5.7)
10. **Major.** Correct `sec_data.tex:57-58` on when the expectations were written. (Sect. 3, item
    6)
11. **Minor.** Fix the sign in `app.tex:41`; say that A contains the prior precision; rename the
    projector.
12. **Minor.** Reword `abstract.tex:11-12` ("exactly an exponential prior").
13. **Minor.** `sec_legacy.tex:159`: say that 13.1 and 13.6 are from a toy with two indices and the
    recorded step sizes, give the 11.3 of Eq. `eq:Rn`, give the errors of the measured widths (8
    and 2 percent), report the synchrotron shift too, and call the factor 2.44 post hoc where it is
    used.
14. **Minor.** Replace "to 0.2 percent" in the abstract by what Table 1 shows; add the conditions of
    Eq. `eq:Tk`; add the prior ratio to step 3; "asymptotically" for (ln n)^-1/4.
15. **Minor.** Correct the numbers of Sect. 3 rows 14, 15, 20, 24, 26, 36, 39, 42, 45, 47; mention
    the drift correction; remove "an upper estimate"; stop using the discarded nside-8 number or
    say why it is still valid for a difference.
16. **Minor.** Bardenet et al.: "recover" overstates. Mark the statement about Alonso et al. as an
    inference. Attribute the quotation in Sect. 6.6 correctly.
17. **Minor.** References: BeyondPlanck XIII; Planck 2015 X; Dickinson et al. 2009;
    Armitage-Caplan et al. 2012; Jew and Grumitt 2019; Gratton 2008; Ceperley and Dewing 1999;
    Nicholls et al. 2012; Medina-Aguayo et al. 2016; Alquier et al. 2016; a reference for the
    errors-in-variables statement (Fuller 1987 or Patefield 1977, read first); Jarosik et al. 2007
    and 2011 for the matrices; Planck 2018 II for PR3. Complete the Huber and White entries. Check
    the Svalheim and Watts quotes against the PDFs, as both citation reports ask.
18. **Minor.** State the conditions of the zero-mean identity (one component, tau maximized); fix
    the improper log-flat prior on tau; restore det S in the Gaussian-prior marginal.
19. **Minor.** State that the boxed Eq. `eq:shift` is for one component and that with a free CMB
    the curvature of the volume term can be negative.
20. **Minor.** Remove the forward reference at `sec_profile.tex:72-74` or finish the test it
    promises. Remove the stale pending note in the Beams paragraph. "half-ring", and the 353 GHz
    bolometers for PR4.
21. **Minor.** "two most sensitive synchrotron channels there are"; "mostly E"; "the field mean is
    recovered"; the E006 sentence (give the offset of the target from the input).

### Suggested

- Split the paper as in 5.12.
- Cite E022 for the band set reading: it uses the published maps, and its flat-CMB case already
  puts the high-latitude flat-prior marginal at the grid edge. It is better evidence than PR4.
- Table 4 compares estimators on the all-scale fit, which the paper then calls biased. Repeat the
  comparison in whatever version is adopted, or say explicitly that the comparison is differential
  and unaffected.
- Fig. 4: draw the reference line at the adopted value, and add Fuskeland et al. (2014).
- Decide whether the Commander3 `marginal` sign error (E009 item 3, used for the AME peak frequency
  in BeyondPlanck) belongs in a paper about volume terms in this code family. If the paper speaks
  about Commander3 at all, a referee who knows the code will ask.
- Give the plain-language relation of the low-pass estimate to T-T plot estimates, which weight
  bright compact structure; they may differ for a physical reason.
- Add a short table of every analysis choice that moves the K to V index, with its size: estimator,
  reference frequency, window, scale cut, band subset, CMB prior, dust index, bandpass drift,
  poorly measured modes. Most of the numbers exist in the records.

## 7. What I did not check

- I did not run Commander4, the dense real-data fits, or any script under `nextgen/scripts/`. All
  real-data numbers were compared with the text outputs, not regenerated.
- I did not read `.npz` files, the figure scripts `paperA1_fig_*.py`, `fit.py`, `spectra.py`,
  `window.py`, `lfi_noise.py` or `data.py` beyond a search for file names. I did not check the
  window operator, the harmonic basis or the low-pass projection numerically.
- I did not verify any cited paper myself. Sect. 4 compares the draft with the two citation
  reports. Statements marked "(from memory)" (the errors-in-variables variance in Fuller, the
  synchrotron B to E ratio) are mine and unverified.
- I did not read E002, E005, E006, E012 or E017 as records; their numbers were checked through
  E035, the registry and searches. E020, E021 and E032 were searched for the quoted numbers only.
  The first 180 lines and the addenda of E036 beyond the tables were not read.
- The stats review of the legacy move was read in Sects. 2 to 5 only; the lowres review in its
  verdict and first derivations only.
- My count of degrees of freedom for the low-pass fits is an estimate from the script, not from a
  printed value. My aliasing estimate in 5.3 (c) uses Gaussian beams.
- I did not check the 48-region runs, `lowpass30`, or anything else still being written.
- I did not check the PDF layout beyond the four figure files; Figs. 1 and 3 were not opened.
- The red pending items were not counted as errors. Where a conclusion depends on one, it is
  named in 5.2, 5.4, 5.5 and 5.11.
