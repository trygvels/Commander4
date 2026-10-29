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
(pending)

## Negative control
The E023 data-built cells (biased) and exact-template cells (clean) bracket the cs cells; the
E021 control sky has been run for the scan code.

## Conclusion
(pending)

## Follow-ups
