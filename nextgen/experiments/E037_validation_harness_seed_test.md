# E037: is the offset of the exact ridge target from the input on the validation simulation the smoothed noise, or a bug?

- **Date started / finished:** 2026-11-09 23:40 / running
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

## Command

```bash
cd /Users/trygvels/code/Commander4
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python nextgen/scripts/e037_harness_seed_test.py \
  --nseed 40 > nextgen/experiments/E037/seed_test.txt
```

## Result

Pending.

## Negative control

Variant (w) is the control of the target, variant (c) of the explanation.

## Conclusion

Pending.

## Follow-ups
