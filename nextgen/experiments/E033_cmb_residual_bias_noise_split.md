# E033: the GLS CMB residual metric split into foreground bias and noise

- **Date started / finished:** 2026-11-09 22:10 / 2026-11-09 22:21
- **Branch and commit:** `nextgen` at 52061ad
- **Open question addressed:** review item (referee objection 7; the CMB column of Paper A Table 3);
  Q3 (what an index field buys for the CMB)
- **Validation level:** diagnostic (no sampler or target change)

## Question
The "GLS CMB residual" reported for every chain since E020 is the rms of the flat-prior GLS CMB
estimate at the chain's index fields minus the true CMB. An inline check on 2026-11-09 21:11
found that at ladder noise scale 1 this number is the estimator's noise (93.13 uK at the true
SEDs) and that it falls when the assumed synchrotron index is steeper (90.95 uK at -0.30), so
that chain 3's "improvement" over chain 1 (90.90 against 96.42 uK) may carry no information about
foreground residuals. How large is the foreground-bias part of the metric for each scale 1 chain,
does it rank the chains differently from the total, and does any chain's index field give a
smaller CMB foreground residual than a constant index? The answer decides whether the CMB column
of Table 3 can stay, and gives the baseline row the referee asked for.

## Target and priors
No sampling. The estimator is linear in the data: per pixel c_hat = w(theta)^T d with
w^T = e_c^T (M^T N^-1 M)^-1 M^T N^-1, M the (bands, 3) mixing matrix of CMB, power-law
synchrotron and modified-blackbody dust at the assumed (beta_s, beta_d, T_d), N the diagonal
white noise covariance, and flat priors on the three amplitudes. With d = sky + n,

    c_hat - c = [w^T sky - c] + w^T n,

exactly. The first bracket is the foreground bias of the assumed SEDs (zero when they equal the
true per-pixel SEDs and the sky follows the model), the second is noise with variance
[(M^T N^-1 M)^-1]_cc per pixel, which depends on the assumed SEDs and not on the sky.

## Setup
`/Users/trygvels/data/e023_pysm_n1` (PySM d1s1c1, nside 64, 12 bands, noise scale 1, seed 2027):
`truth.npz` holds the noiseless band signal `sky`, the true dust band column `dust_bands`, the
true CMB and the true index and temperature maps; the band maps hold sky plus noise. Script
`nextgen/scripts/e033_cmb_bias_noise_split.py` (imports `gls_cmb` from `e020_analyze.py`, so the
estimator is the one behind every earlier number). Rows: the true per-pixel SEDs (control); the
true indices with the oracle region-mean T map; the truth's nside-16 region means; constant
fields at the truth means and at the prior means; the noise-free region-fit references
(`E027/refs_nside16.npz`); and the posterior-mean fields of the scale 1 chains (stage A'
`gmrf16_exact_mumap_nc_Tmap_mu01`, chains 1, 2, 3 of E030, E026 `gmrf16_cs_mumap_nc` and
`gmrf16_exact_mumap_nc` with the fixed T). Each row: bias rms (all sky, |b| > 20, and per
synchrotron S/N quartile), analytic noise rms, realized noise rms, total rms, and the closure
total against sqrt(bias^2 + noise^2). Mac mini, seconds per row.

## Reference
The control row: with the true per-pixel SEDs (true beta_s, the true dust band column) the bias
term must vanish to rounding, because the sky is CMB plus a per-pixel power law plus the dust
column. The total at the true SEDs must reproduce the E020 floor (93.13 uK at scale 1). The
analytic noise rms must agree with the realized one to its sampling error (0.30 uK for 49152
pixels).

Pre-registered before the script is run (2026-11-09 22:10), from the inline check and rough
leakage estimates:
- P1 (control): bias rms below 0.05 uK at the true per-pixel SEDs. If it fails, the sky has
  structure outside the per-pixel model and the floor has to be redefined first.
- P2 (closure): |total - sqrt(bias^2 + realized noise^2)| below 0.5 uK in every row, and the
  realized noise rms within 1 uK of the analytic one.
- P3 (the total metric is the noise term): the analytic noise rms at chain 3's fields is 5.0 to
  6.0 uK below that at chain 1's fields, which is the whole of the 5.5 uK difference in the total.
- P4 (size of the bias): all-sky bias rms below 5 uK for the exact-template chain (stage A'),
  5 to 15 uK for the cs-template chains 1, 2, 3, and chains 1 and 3 within 30 percent of each
  other.
- P5 (skill against a constant): the exact-template chain's all-sky bias is at most half that of
  the constant field at the truth means; chain 1 beats the constant; chain 3 does not beat it by
  more than 20 percent (its brightest regions sit -0.09 from the reference).
- P6 (high latitude): at |b| > 20 every bias is below 8 uK, and chain 1's is 1.5 to 2.5 times
  chain 3's (the ratio of the faint displacements, +0.30 against -0.16).
These are guesses from one-pixel leakage estimates (leak about 0.5 a_s(30 GHz) delta beta) and
are expected to fail in places; the point is to fix them before looking.

## Command
```bash
cd /Users/trygvels/code/Commander4 && source .venv/bin/activate
D=/Users/trygvels/data/e023_pysm_n1
python nextgen/scripts/e033_cmb_bias_noise_split.py --dir $D \
  --refs nextgen/experiments/E027/refs_nside16.npz --t-map $D/truth_Td_nside16.fits \
  --tags gmrf16_exact_mumap_nc_Tmap_mu01 gmrf16_cs_mumap_nc_Tmap_mu03 \
         gmrf16_cs_mumap_nc_Tmap_mu05 gmrf16_cs_mumap_nc_Tmap_mu03_ref23 \
  --tags-fixed-t gmrf16_cs_mumap_nc gmrf16_exact_mumap_nc \
  > nextgen/experiments/E033/split_n1.txt
```

## Result
Run 2026-11-09 22:10:43 to 22:10:50, rerun 22:20 with the band-power block added (the first 45
lines are unchanged). Output `E033/split_n1.txt`. The weights reproduce `gls_cmb` to 4e-12 uK.
All numbers are rms in uK_CMB at nside 64 and noise scale 1.

| assumed fields | bias, all sky | bias, abs(b) above 20 | bias per sync quartile, bright to faint | analytic noise | total |
|---|---|---|---|---|---|
| true per-pixel SEDs (control) | 0.00 | 0.00 | 0.00, 0.00, 0.00, 0.00 | 93.13 | 93.18 |
| true indices, region-mean T map | 1.73 | 0.22 | 3.32, 0.92, 0.32, 0.13 | 93.12 | 93.20 |
| truth nside-16 region means, T map | 5.95 | 0.42 | 11.83, 1.13, 0.63, 0.30 | 93.12 | 93.36 |
| noise-free region fit (trueT indices), T map | 4.43 | 0.45 | 8.74, 1.27, 0.68, 0.30 | 93.14 | 93.33 |
| constant at the truth means, T map | 16.31 | 0.95 | 32.34, 3.81, 1.76, 0.60 | 93.08 | 94.51 |
| constant at the truth means, fixed T | 19.48 | 0.83 | 38.76, 3.57, 1.40, 0.57 | 93.00 | 95.00 |
| constant at the prior means (-3.1, 1.6), T map | 10.50 | 1.65 | 20.50, 3.94, 1.93, 1.09 | 90.64 | 91.28 |
| truth region means, sync shifted -0.30 | 5.84 | 1.44 | 11.42, 1.99, 1.25, 0.81 | 90.93 | 91.19 |
| truth region means, sync shifted +0.30 | 8.47 | 1.68 | 16.70, 2.29, 1.45, 0.92 | 95.97 | 96.37 |
| stage A' chain (exact templates, sd 0.1, T map) | 5.60 | 0.84 | 10.74, 2.85, 1.30, 0.60 | 92.83 | 93.08 |
| E030 chain 1 (cs, sd 0.3, reference 30 / 545) | 9.83 | 3.34 | 17.29, 7.85, 4.41, 2.59 | 95.98 | 96.42 |
| E030 chain 2 (cs, sd 0.5) | 6.70 | 1.56 | 12.62, 3.80, 2.04, 1.28 | 93.78 | 94.04 |
| E030 chain 3 (cs, sd 0.3, reference 23 / 353) | 8.17 | 2.05 | 15.41, 4.70, 2.32, 1.53 | 90.49 | 90.90 |
| E026 cs chain, fixed T | 7.63 | 1.77 | 14.23, 4.69, 2.46, 1.50 | 94.56 | 94.90 |
| E026 exact chain, fixed T | 6.98 | 1.15 | 13.60, 2.77, 1.23, 0.73 | 93.10 | 93.46 |

Attribution, all-sky bias with one component's field replaced by the truth: chain synchrotron with
the true dust column 0.89 (stage A'), 4.02 (chain 1), 2.41 (chain 2), 2.09 (chain 3); true
synchrotron with the chain's dust 5.51, 8.33, 6.70, 7.19. The dust field carries most of the bias
in every chain. Each chain's field against a constant at that chain's own field means: 5.60
against 14.18, 9.83 against 15.49, 6.70 against 11.35, 8.17 against 10.63.

Pre-registered items:

- P1 pass: bias 0.00 at the true per-pixel SEDs, in both forms of the control.
- P2 pass: closure within 0.07 uK in every row; realized noise within 0.07 uK of the analytic.
- P3 pass: the analytic noise rms is 95.98 at chain 1's fields and 90.49 at chain 3's, a difference
  of 5.49 against 5.52 in the totals.
- P4 fails for the exact-template chain (5.60 against "below 5") and passes for the rest (9.83,
  6.70, 8.17, chains 1 and 3 within 20 percent). The threshold ignored the resolution floor: the
  truth's own region means give 5.95, so 5.60 is at the floor of an nside-16 field.
- P5 passes for stage A' (5.60 against half of 16.31) and for chain 1, and fails for chain 3: its
  bias is 50 percent below the constant's, not within 20 percent. The guess assumed the
  synchrotron displacement sets the all-sky bias; the bright quartile and the dust field set it.
- P6 pass: every bias at abs(b) above 20 is below 3.4 uK; chain 1 over chain 3 is 1.63.

Band powers (exploratory, added after the first read, not pre-registered): sqrt of the mean
D_l of the bias map at abs(b) above 20, in uK, for l = 2 to 10, 11 to 30, 31 to 100, and the ratio
of bias power to noise power in the first bin at this noise scale and at 0.03 times it (the bias
does not depend on the noise scale, the noise power scales with its square).

| assumed fields | 2 to 10 | 11 to 30 | 31 to 100 | bias/noise power, l 2 to 10, scale 1 | at scale 0.03 |
|---|---|---|---|---|---|
| truth region means, T map | 0.024 | 0.086 | 0.319 | 3.4e-5 | 0.038 |
| noise-free region fit, T map | 0.065 | 0.131 | 0.325 | 2.5e-4 | 0.28 |
| constant at the truth means, T map | 0.315 | 0.455 | 0.527 | 6.0e-3 | 6.7 |
| constant at the prior means | 0.708 | 0.590 | 0.509 | 3.2e-2 | 36 |
| stage A' chain | 0.289 | 0.393 | 0.457 | 5.1e-3 | 5.6 |
| chain 1 | 1.916 | 1.400 | 0.908 | 0.21 | 231 |
| chain 2 | 0.792 | 0.609 | 0.571 | 3.7e-2 | 42 |
| chain 3 | 0.934 | 0.715 | 0.583 | 5.6e-2 | 62 |
| E026 exact chain, fixed T | 0.599 | 0.478 | 0.481 | 2.2e-2 | 24 |

## Negative control
The control row doubles as the test of the split: a wrong sky file or a unit error shows as a
nonzero bias at the true SEDs.

## Conclusion
(read 2026-11-09 22:12 and 22:22)

1. The total flat-prior GLS CMB residual at noise scale 1 is the noise term. It is 93 uK of noise
   with 6 to 10 uK of bias in quadrature, and the noise term itself moves by 5.5 uK between the
   chains because a steeper assumed synchrotron index lowers the noise variance of the estimator
   (90.9 at a shift of -0.30, 96.0 at +0.30). Every statement in E026 to E030 and in Paper A
   Table 3 that ranks models by this total at scale 1 is withdrawn, including "chain 3 improves the
   CMB residual" (item (vi) of its pre-registration passes without weight) and the E029 attribution
   of chain 1's excess to its synchrotron field. The bias term is the metric from now on.
2. On the bias term the all-sky number is a Galactic-plane number: the bright quartile carries
   10 to 39 uK and the other three 0.1 to 8 uK. It is set by the dust field and by the nside-16
   resolution (floor 5.95 at the truth's own region means, 1.73 from the region-mean T map alone).
   Every chain beats a constant index there, by 30 to 60 percent.
3. Away from the plane the ranking reverses. At abs(b) above 20 the constant at the truth means has
   bias 0.95 uK; the exact-template chain matches it (0.84) and every cs-template chain is worse
   (1.56 to 3.34). In the largest scales, l = 2 to 10, chain 1 has six times the bias amplitude
   of the constant. At noise scale 1 the index fields of E026 to E030 do not improve the
   high-latitude CMB over a constant index; with a 30 percent template error they degrade it. This
   is the regime of E031 in which a region's index is not constrained by the data, so the field's
   error (0.18 to 0.43 per pixel) exceeds the true variation (0.063).
4. Size of the stake in intensity: at abs(b) above 20 a constant index leaves a bias of 0.3 to 0.5 uK
   per band-power bin, against a CMB signal near 30 uK in the same units. At WMAP/LFI depth that
   bias is 6.7 times the noise power at l = 2 to 10 and below the noise above l = 30; a field at
   the truth's region means removes it (ratio 0.04). The gain from index fields for the intensity
   CMB is real and small. It is a foreground-science gain first.

What this does not establish: anything at noise scale 0.03 for a sampled field (E032), and
anything for polarization, where the same bias is compared with a signal of order 0.1 to 1 uK.

Changes: F40 added; the CMB column of Paper A Table 3 is to be replaced by the bias term at
abs(b) above 20 and its low-l band power; `e029_cmb_attribution.py` numbers are superseded.

## Follow-ups
- The same split at noise scale 0.03 for the E032 chains, where the noise term is 2.79 uK and
  the bias term can be seen directly (pre-registered as item (x) of E032).
- The band-power comparison deserves a mode-decoupled estimator and a mask before it enters a
  paper; here it is a pseudo-spectrum over f_sky.
