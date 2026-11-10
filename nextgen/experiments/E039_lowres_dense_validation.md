# E039: region indices of polarized synchrotron with the released WMAP noise covariances, validation on simulations

- **Date started / finished:** 2026-11-10 20:32 / running
- **Branch and commit:** `nextgen` at 747ed50
- **Open question addressed:** Q21, Q20; Gates 1 and 3 of the challenge report, reduced to regions
- **Validation level:** L0 and L1 (unit tests against closed forms), then L4 in the sense of a
  realistic amplitude field that is drawn from no prior, with the real noise covariance

## Question

Paper A1's data section estimates one synchrotron index per sky region from WMAP K, Ka, Q, V and
Planck 353 GHz at nside 16, with the dense WMAP inverse covariances as released. Before any real
index is read: does the dense collapsed likelihood of `nextgen/lowres/` return unbiased region
indices with calibrated errors when the noise is drawn from those covariances, and do the three
targets behave as E036 found with white noise (profile calibrated once its width carries the
sandwich factor, flat-prior marginal shifted by the reference frequency, zero-mean prior equal to
the profile)?

What changes with the answer: a pass clears the estimator for the real maps (E040). A fail in
bias or coverage is found and fixed first.

## Target and priors

Per band b, in thermodynamic units, elements p (Q then U): d_b = s_b(beta_r(p)) a_s + f_b a_d + a_c
+ n_b, n_b ~ N(0, N_b). s_b is the band-averaged power law at nu_ref = 23 GHz (WMAP nine-year
bandpass profiles; delta band at 353 GHz), f_b the band-averaged modified blackbody (index 1.53,
19.6 K, reference 353 GHz), the CMB column is 1. N_b^-1: the released WMAP res 4 matrices; for 353
GHz the PR4 per-pixel Q/U covariance combined over sub-pixels. Priors: flat on a_d; N(0, 0.2^2
uK^2) per element on a_c; a_s flat (profile and flat-prior marginal) or N(0, tau_r^2) with tau_r
maximized jointly (zero-mean prior). beta_r flat on [-4.6, -1.6]. Objectives and the sandwich
covariance J^-1 (J + T) J^-1 are in `nextgen/lowres/likelihood.py`; derivation of T below.

Derivation of the score covariance. With P = N^-1 - N^-1 M A^-1 M^T N^-1 and H = M A^-1 M^T N^-1,
the profile is chi2_min = d^T P d and its score is s_r = -(1/2) d^T P_r d with P_r = dP / dbeta_r.
For d = M a + n: P M = 0 gives P_r M a = -P M_r a, so the mean of the score is zero and its
signal part has variance (M_r a)^T P N P (M_r' a) = (M_r a)^T P (M_r' a) = J_rr'. The noise part is
(1/2) tr(P_r N P_r' N) = (1/2) tr(H_r H_r'), and with H_r = (1 - H) M_r G + M A^-1 M_r^T N^-1 (1 - H),
G = A^-1 M^T N^-1, the two squared terms vanish (G (1 - H) = 0 and (1 - H) M = 0) and the cross
terms give T_rr' = tr(A^-1 M_r^T P M_r'). The expected curvature of chi2_min / 2 is J. For one
amplitude per element and diagonal noise this is J = sum u V and T = sum V, the form checked in
the E036 addendum.

## Setup

- Unit tests (`nextgen/lowres/tests/test_likelihood.py`, five tests, pass at 747ed50 plus the
  module): analytic gradients of the three objectives against finite differences in four
  configurations; equality with the per-element closed forms for diagonal noise; J and T against
  sum u V and sum V; the sandwich covariance against the scatter of 300 fits with dense noise,
  three components and low S/N; the reference-frequency identity of the flat-prior marginal.
- Simulation: truth synchrotron amplitudes from the Cosmoglobe DR1 polarized synchrotron map (5
  deg, scaled from 30 to 23 GHz with each region's true index), dust amplitudes from the PR4 353
  GHz map at the working resolution, CMB drawn N(0, 0.2^2). Noise drawn from the real covariances.
  Bands K, Ka, Q, V, 353.
- Region sets and true indices, chosen away from the centre of the allowed range: the four
  BeyondPlanck XV regions with -3.15 (High Latitude), -3.25 (Spur), -2.85 (Galactic Center), -2.95
  (Plane); and the 48 nside-2 pixels with -3.1 + 0.12 g, g a fixed standard normal draw.
- Bulk at nside 8 (data and covariances degraded as the real nside-8 analysis would be), a spot
  check at nside 16. Profile fits for every realization; flat-prior marginal at references of 23
  and 30 GHz and the zero-mean prior for a subset.

## Reference

The truth is known. z = (estimate - truth) / sd per region, with the naive sd from J^-1 and with
the sandwich sd.

**Predictions, written before the first simulation run:**

- V1. Profile, unbiased: mean z (sandwich) within +-0.2 in every BP XV region, and within +-0.3
  for the nside-2 regions whose sandwich sd is below 0.25, with 200 realizations.
- V2. Profile, width: sd of z with the sandwich sd between 0.9 and 1.1; with the naive sd between
  1.0 and 1.5 and never below the sandwich value by more than 0.03.
- V3. Coverage of +-1 sandwich sd between 0.64 and 0.72 over the regions with sd below 0.25.
- V4. Flat-prior marginal: optimum minus profile optimum has opposite signs at references of 23
  and 30 GHz or, if the pivot lies below 23 GHz as with a free CMB in E036, is larger at 30 than
  at 23 by ln(30 / 23) times the number of amplitudes times the squared profile sd, within 25
  percent, in the regions where that product is below 0.3.
- V5. Zero-mean prior with tau maximized: within 0.3 sandwich sd of the profile optimum in the BP
  XV regions.
- V6. The reduced chi-square at the profile optimum is 1.00 +- 0.03 (the model is the truth).

## Command

```bash
cd /Users/trygvels/code/Commander4
.venv/bin/python -m pytest -q -p no:cacheprovider nextgen/lowres/tests
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions bpxv --nreal 200 --full 40 \
  > nextgen/experiments/E039/sim_n8_bpxv.txt
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions hp2 --nreal 200 --full 0 \
  > nextgen/experiments/E039/sim_n8_hp2.txt
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 16 --regions bpxv --nreal 12 --full 0 \
  > nextgen/experiments/E039/sim_n16_bpxv.txt
```

## Result

Pending.

## Negative control

The flat-prior marginal (V4) is the known-biased target.

## Conclusion

Pending.

## Follow-ups
