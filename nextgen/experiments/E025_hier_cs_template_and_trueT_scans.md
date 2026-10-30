# E025: Stage 2 of Q16 by exact-target scans: component-separation-quality templates and the true-T dust SED under the hier prior

- **Date started / finished:** 2026-10-30 /
- **Branch and commit:** nextgen, e6f4a79
- **Open question addressed:** Q16 (D022); M2.5 / M2.6
- **Validation level:** L3 (misspecified simulation with known truth), exact target only

## Question
Is the collapsed hier index target "unbiased at this depth" (E023 rule: every quartile's
preferred shift below 0.02 sync and 0.01 dust at 0.03 x noise, inside the block sd at noise 1,
wide sync scan within 0.1) when the template error is a smooth multiplicative field with no
coupling to T_d or to the index (the error a component-separation template would carry), with
the dust SED at the fixed T and at the true T_d, and with a free template multiplier per
nside-16 region? If yes for some configuration, that configuration is the first stage-2
chain (gmrf16, noise-1 sky, flat CMB prior, dense metric); if no for all, the hier prior for
dust at high S/N needs sampled T_d or is replaced by a weak Gaussian prior (D022).

## Target and priors
As E023 (stats-reviewed there) with the CMB amplitude prior sd at 1e4 (D021, F26):
L(beta) = max over theta of the expected log-likelihood on the noiseless sky with the trace
term, hyperparameters profiled, no index prior; layouts "global" (mu_c, s_c) and "perregmu"
(mu_c per nside-16 region with the N(1, 0.5^2) hyperprior each, s_c per template-S/N quartile).
Templates: g_c = a_c exp(f_c) with f_c a Gaussian random field of rms 0.3 and Gaussian
correlation FWHM 10 deg (seed 11; `--cs-scatter 0.3 --cs-fwhm 10 --cs-seed 11`), floored
positive; so a | g = g exp(-f), a log-normal multiplicative error that is not the prior's law
(F23: a known misspecification, but smooth and uncoupled). Dust SED: one MBB at the fixed T or
at the PySM T_d,p (reference for sampled T). Reference, quartiles, scans, Laplace term (global
layout only) as E023.

## Setup
Skies `/Users/trygvels/data/e020_d1s1` (0.03 x) and `/Users/trygvels/data/e023_pysm_n1` (1 x);
`nextgen/scripts/e023_exact_scan.py` with the new "cs" template level; cells cs_fixT_global,
cs_trueT_global, cs_fixT_perregmu, cs_trueT_perregmu per sky, expected mode; both skies in
parallel on the Mac mini (free), about 4 h.

## Reference
As E023: the noise-free region fit per SED level; the E023 exact-template cells as the
unbiased reference (clean at both depths with the flat CMB prior at noise 1) and the E023
data-built cells as the biased one. Pre-registered reading: the E023 rule per cell; attribution
by the fixed-T against true-T pair (what sampled T buys) and by the global against perregmu
pair (what a template-error field buys); the first configuration that passes at both depths
goes to the stage-2 chain; "none passes" sends dust to sampled T_d (M2.6) or a weak Gaussian
prior.

## Command
```bash
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e020_d1s1 --out nextgen/experiments/E025/scan_n003.txt --cmb-sd 1e4 --cells cs_fixT_global cs_trueT_global cs_fixT_perregmu cs_trueT_perregmu --real-cells
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e023_pysm_n1 --out nextgen/experiments/E025/scan_n1.txt --cmb-sd 1e4 --cells cs_fixT_global cs_trueT_global cs_fixT_perregmu cs_trueT_perregmu --real-cells
```

## Result
Both skies 00:40 to 04:14 (`E025/scan_n003.txt`, `E025/scan_n1.txt`); cs template realized a/g
= exp(-f) mean / sd 1.007 / 0.309 (sync), 1.060 / 0.319 (dust), per quartile 0.99 to 1.11 /
0.30 to 0.36. Gains in nats (2 nats about 2 block sd); "edge" = the scan boundary at 0.1.

**0.03 x noise.**
- cs_fixT_global (profiled mu_s 0.991, s_s 0.288, mu_d 1.029, s_d 0.302): sync 0.000, -0.010,
  +0.020, -0.020 (gains 0, 3.7, 5.3, 6.3: 2.7 to 3.5 sd), dust 0.000, 0.000, +0.070, +0.090
  (gains 278 and 139 in Q3 and Q4); wide sync 0.000. A smooth, uncoupled 30 percent template
  error under the global law biases the two fainter dust quartiles by 0.07 to 0.09 and the sync
  at the 0.02 threshold: the F23 mechanism (a | g with g = a exp(f) is not the prior's law).
- cs_trueT_global: the same within 0.01 (sync Q4 -0.030, 12.8 nats; dust +0.070, +0.080).
- cs_fixT_perregmu (region means per quartile 0.99 to 1.08 with sd 0.29, tracking exp(-f) at the
  region scale; profiled s_q 0.055 to 0.076, the within-region residual of a 10-deg field):
  sync 0.000, -0.010, -0.010, -0.020 (gains 0, 1.5, 2.0, 2.6), dust 0.000, 0.000, +0.010, -0.020
  (gains 0, 0, 1.4, 2.3); wide 0.000.
- cs_trueT_perregmu: sync 0.000, -0.010, -0.010, -0.020 (2.2, 2.7, 3.9), dust 0.000, +0.010,
  +0.010, -0.030 (1.0, 1.7, 2.8); wide 0.000.

**Noise 1 (flat CMB prior).**
- cs_fixT_global (mu_s 0.829, s_s 0.231, mu_d 0.999, s_d 0.275): sync -0.050 (1.4 nats), edge
  (0.7), -0.040 (0.06), 0.000; dust +0.010 (4.2 nats, 2.9 sd), +0.100 edge (24), +0.100 edge (12),
  +0.100 edge (7.7); wide 0.000. Dust biased beyond 0.1 in Q2 to Q4; sync not significantly.
- cs_trueT_global: the same within 0.01 (dust 5.1, 21, 11, 7.1 nats).
- cs_fixT_perregmu (region means per quartile 0.96 to 1.06; s_q 0.055 to 0.235): sync -0.020
  (0.3), +0.090 (0.5), edge (1.1), edge (1.0); dust 0.000, -0.040 (1.1), +0.020 (0.1), edge
  (1.4); wide +0.2 (1.4). No gain above 1.5 nats.
- cs_trueT_perregmu: the same within 0.01 (gains 0.1 to 1.6).

Attribution: fixed against true T changes nothing in any cell at either depth; global against
per-region multiplier changes everything for dust (0.07 to 0.09 and beyond 0.1 to within 0.02
and to insignificance) and brings the sync to the threshold.

## Negative control
The E023 data-built cells (biased) and exact-template cells (clean) bracket the cs cells; the
E021 control sky has been run for the scan code.

## Conclusion
1. **A per-region template multiplier makes the hier prior usable with a component-separation-
   quality template.** With g = a exp(f), f smooth (10 deg) and 30 percent rms, uncoupled to the
   sky, the global hier law biases the fainter dust quartiles by 0.07 to 0.09 at 0.03 x noise and
   beyond 0.1 at realistic noise (the F23 law mismatch, now on a realistic error); with one free
   multiplier per nside-16 region (3072 per component, N(1, 0.5^2) hyperprior each) and a scatter
   per quartile, every shift is within 0.02 at 0.03 x noise (sync Q4 -0.02 at 2.3 to 3.9 sd,
   dust Q4 -0.02 to -0.03 at 2.1 to 2.4 sd; all others 0.00 to 0.01) and no scan at realistic noise
   has a gain above 1.6 nats (the largest shifts, 0.04 to 0.1, are in quartiles with block sds of
   0.03 to 0.13). Against the pre-registered rule: at 0.03 x the sync passes (0.02) and the dust
   fails by 0.01 to 0.02 in its faintest quartile (threshold 0.01); at realistic noise the rule
   ("inside the block sd") is met by 14 of 18 scans and missed at 1.0 to 1.7 sd by four, which is
   what unbiased scans do; so "consistent with unbiased at the resolution of this depth", not
   "unbiased" by the letter.
2. **Fixed T_d is not the problem.** In none of the eight cells does the true-T SED change a
   shift by more than 0.01: at these depths the dust temperature misfit is below everything
   else. Sampled T_d is deprioritized (D023).
3. **What the stage-2 chain needs:** the hier prior with a template-error field at the region
   scale, i.e. mu_c as a map at nside 16 (or a GMRF on log mu_c), not a scalar, in Commander4's
   `AmpPrior`; the per-quartile scatter is secondary (the profiled residual scatter is the
   within-region part of the error, 0.06 for a 10-deg field). The condition for a real template
   is that its within-region error be at the few-percent level; the region-scale error can be
   anything smooth.
4. Not shown: templates with errors at all scales (a real component separation), the posterior
   under the GMRF index prior (chains), the multiplier field's own prior (a flat N(1, 0.5^2) per
   region here; a GMRF would be the model), more than one error realization. F28.

## Follow-ups
- DONE 2026-10-30 (commit f6642a3, stats-reviewed): the region-scale template multiplier in
  Commander4's hier `AmpPrior` (`mu_nside`); the stage-2 chains run as E026. Was: then the stage-2 chain
  (gmrf16, noise-1 sky, cs template, flat CMB prior, dense metric) with the E021 metrics.
- A template with errors at all scales (true amplitude times a field with a power law to the
  pixel scale) to set the within-region error the method tolerates.
