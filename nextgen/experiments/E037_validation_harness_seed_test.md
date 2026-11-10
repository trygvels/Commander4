# E037: is the offset of the exact ridge target from the input on the validation simulation the smoothed noise, or a bug?

- **Date started / finished:** 2026-11-09 23:40 / 2026-11-10 01:11
- **Branch and commit:** `nextgen` at 981deb6
- **Open question addressed:** Q7; F8; referee Sect. 3.2 and challenge report item G2
- **Validation level:** L1 (exact grid targets over many seeds, with a white-noise control)

## Question

On the E006 validation simulation the exact ridge target sits 1.9 and 3.7 posterior sd from the
input indices (F8). The record attributes this to smoothed noise being treated as independent
between pixels. That was one realization, and by working principle 8 an unexplained offset under
the harness that carries F2 is a bug until shown otherwise. Over many seeds, is the offset a
zero-mean scatter whose size matches the noise correlation, and does it vanish when the noise is
white?

What changes with the answer: if yes, F8 is confirmed and sized, and Q7 gets a number for the
harness configuration. If the mean offset is not zero, or the white-noise control does not have
unit scatter, there is a bug in the harness or in the ridge target, and F2's end-to-end evidence
has to be redone.

## Target and priors

The exact 2-d ridge (profile) and flat-prior marginal posteriors of (beta_s, beta_d) on a grid,
flat prior inside the harness bounds, per-pixel diagonal noise with the rms the harness writes.
Independent NumPy implementation taken from `aux/validate_mcmc_likelihood.py`
(`_exact_posterior`, `mixing`).

## Setup

The harness simulation, in memory, with the seed as a parameter: nside 16, seven bands 30 to 353
GHz, CMB, synchrotron (-3.0) and dust (1.6) with the harness amplitude model, rms 60, 60, 70, 24,
16, 24, 80 uK_RJ. Three variants per seed, same sky and same white noise draw:

- (h) the harness: signal and noise smoothed to 480 arcmin at lmax 47, rms map set to the sd of
  the smoothed noise, pixels treated as independent;
- (w) white control: no smoothing of signal or noise, rms as simulated;
- (c) the harness maps with the noise variance multiplied by N_eff, the number of pixels per
  smoothing kernel, 4 pi sigma_b^2 / Omega_pix = 10.8, which is the crude correction for counting
  each independent noise patch once.

40 seeds. z = (posterior mean - input) / posterior sd per index.

## Reference

**Predictions, written before the first run:**

- P1. Variant (w): mean z within +-0.4 and sd of z between 0.8 and 1.25 for both indices under
  the ridge target (the profile is calibrated at this S/N).
- P2. Variant (h): mean z within +-0.6 sd(z)/sqrt(40) times 2, that is no significant bias, and
  sd of z between 2.2 and 4.0, near sqrt(10.8) = 3.3. The single realization of F8 (1.9 and 3.7)
  is then typical.
- P3. Variant (c): sd of z between 0.7 and 1.3.
- If (w) fails, the harness or the target has a bug. If (h) shows a mean offset above 3 standard
  errors, the smoothing biases the index, which is a different statement from F8.

**Added control, predictions written 2026-11-10 00:57 after the 40-seed run was read and before this was
run.** The white-noise control failed P1 (sd of z 2.6 and 7.8). By the rule above that reads as a
bug in the harness or the target. There is a statistical explanation that must be excluded
first: the harness has a per-pixel amplitude S/N far below 1 over most of the sky, and the
profile's curvature understates the error by the sandwich factor sqrt(1 + 1/u) there (challenge
report Sect. A1b item 3; E036 at four times the noise). The discriminating run is the white
control with every rms divided by 10 and by 30.

- P4. At one tenth of the noise the white control has sd of z between 0.85 and 1.3 and mean z
  within +-0.5 for both indices, and at one thirtieth between 0.85 and 1.2. If so, the failure of
  P1 is the profile's width at low per-pixel S/N and not a bug.
- If sd of z stays above 1.5 at one thirtieth of the noise, the target or the simulation has a
  bug and F2's end-to-end evidence is re-examined.

## Command

```bash
cd /Users/trygvels/code/Commander4
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python nextgen/scripts/e037_harness_seed_test.py \
  --nseed 40 > nextgen/experiments/E037/seed_test.txt
for s in 0.1 0.0333; do OMP_NUM_THREADS=1 .venv/bin/python nextgen/scripts/e037_harness_seed_test.py --nseed 40 \
  --variants w --noise-scale $s > nextgen/experiments/E037/seed_test_white_noise_x$s.txt; done
```

## Result

**The 40-seed run** (`E037/seed_test.txt`, 2026-11-10 00:39 to 00:55, one thread).

| variant, target | beta_s: mean z (se), sd of z | beta_d: mean z (se), sd of z | median posterior sd, beta_s / beta_d |
|---|---|---|---|
| (h) harness, ridge | -0.18 (0.69), 4.39 | -0.51 (0.76), 4.83 | 0.127 / 0.098 |
| (w) white control, ridge | -0.99 (0.41), 2.59 | +2.08 (1.24), 7.82 | 0.256 / 0.156 |
| (c) harness with variance times 10.8, ridge | -0.16 (0.14), 0.90 | -0.03 (0.18), 1.13 | 0.374 / 0.270 |
| (h) harness, flat-prior marginal | +65 (12), 74 | +1698 (11), 67 | 0.011 / 0.0005 |
| (w) white control, flat-prior marginal | +285 (1.5), 9.7 | +1836 (2.7), 17 | 0.0035 / 0.0005 |

- **P2: no bias, as predicted, and a scatter above the predicted range.** In the harness the
  offset of the exact ridge target from the input has no significant mean (-0.18 +- 0.69 and
  -0.51 +- 0.76) and an sd of 4.4 and 4.8 posterior sd, above the 2.2 to 4.0 I allowed. The single
  realization behind F8 (1.9 and 3.7) is a typical draw.
- **P1 fails.** The white-noise control, which I predicted calibrated, is overconfident by 2.6
  (beta_s) and 7.8 (beta_d). So the smoothing is not the only cause of the F8 offset, and by the
  rule written above the outcome reads as a bug until the added control says otherwise.
- **P3 passes, with a caveat that removes its weight.** With the noise variance multiplied by
  the number of pixels per kernel the sd of z is 0.90 and 1.13. But the posterior is then broad
  (0.37 in a prior box of half-width 1), the synchrotron truth sits at the centre of that box,
  and a broad posterior reverting to the box centre reads as calibrated. The posterior mean
  scatters by 0.34 in (c) against 0.56 in (h) on the same maps, which a rescaling of the noise
  could not do without the box.
- In index units the ridge estimate of the full-sky synchrotron index scatters over seeds by
  0.56 (harness) and 0.66 (white) around -3.0, while its posterior sd is 0.13 and 0.26. The
  harness configuration determines its indices much less well than its posterior says.
- The flat-prior marginal is tens to thousands of posterior sd from the input in both variants.

**The added control** (`E037/seed_test_white_noise_x0.1.txt`, `..._x0.0333.txt`, finished 01:10; 40
seeds each, white-noise control only). Median per-pixel S/N of the synchrotron signal at 30 GHz
and of the dust signal at 353 GHz: 0.062 and 0.041 at the harness noise, 0.62 and 0.41 at one
tenth, 1.86 and 1.24 at one thirtieth.

| noise | ridge, beta_s: mean z (se), sd of z | ridge, beta_d | flat-prior marginal, mean z: beta_s / beta_d |
|---|---|---|---|
| harness level | -0.99 (0.41), 2.59 | +2.08 (1.24), 7.82 | +285 / +1836 |
| one tenth | -0.10 (0.13), 0.82 | -0.02 (0.17), 1.07 | -17 / +592 |
| one thirtieth | -0.05 (0.13), 0.80 | -0.08 (0.15), 0.97 | -10 / +20 |

- **P4 holds where it matters.** With the per-pixel S/N raised by 10 and by 30 the ridge target
  of the white control is calibrated: no bias, sd of z 0.82 and 1.07, then 0.80 and 0.97. The
  synchrotron values sit 0.03 and 0.05 below the lower edge of my band, with a standard error of
  0.09 on each, on the side of intervals slightly too wide. The bug threshold (sd of z above 1.5
  at one thirtieth) is far from met.
- The flat-prior marginal stays 10 to 20 posterior sd from the input even at a per-pixel S/N
  near 2.

## Negative control

Variant (w) at reduced noise is the control of the target and passes. Variant (c) is not a clean
control (see P3).

## Conclusion

1. **There is no bug behind F8.** The same target that is overconfident by 2.6 and 7.8 at the
   harness noise is calibrated when the per-pixel S/N is raised to 0.4 to 0.6 and above.
2. **F8's explanation was incomplete.** The offset of the exact ridge target from the input is a
   zero-mean scatter of 4.4 and 4.8 posterior sd in the harness. Smoothed noise treated as
   independent is one cause. The other, at least as large, is that the harness has a per-pixel
   amplitude S/N of 0.04 to 0.06, where the curvature of the profile likelihood understates the
   error: the white-noise control alone gives 2.6 and 7.8. This is the width failure that the
   challenge report writes as the factor sqrt(1 + 1/u) and that E036 saw at 1.4 in its faintest
   class.
3. **The harness tests a sampler against its target and nothing else.** F2's end-to-end evidence
   (chain against grid) is untouched. What does not hold is any reading of the harness target as
   a calibrated statement about the sky's index: at that noise its intervals are four to five
   times too narrow.
4. A coverage curve for the profile against per-pixel S/N now has three points from this
   experiment (7.8 and 2.6 at 0.04 to 0.06, about 1 at 0.4 to 0.6 and at 1.2 to 1.9) and five from
   E036 (1.0 to 1.4 between 1.1 and 4.5). It is the figure Paper A1 needs, with the sandwich
   factor overlaid.

What this does not establish: the split of the harness overconfidence between the smoothing and
the per-pixel S/N (variant (c) cannot give it, because the synchrotron truth sits at the centre
of the prior box), and the sandwich factor itself, which was not computed.

Changes: F8 amended; the review's assumption table; Q7 keeps its open part (the size of the
smoothing effect on real maps).

## Follow-ups

- The sandwich factor computed for the harness and for E036, against these sd of z.
- The harness with the truth away from the centre of its prior box, before variant (c) is used.
