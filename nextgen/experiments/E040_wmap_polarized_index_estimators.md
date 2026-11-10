# E040: polarized synchrotron index per region from WMAP K to V and Planck 353 GHz, several estimators on the same maps

- **Date started / finished:** 2026-11-10 20:35 (pre-registration; no real-data fit has been run) /
- **Branch and commit:** `nextgen` at 11b2957
- **Open question addressed:** Q21, Q20, Q17; the data section of Paper A1 (D030 item 2)
- **Validation level:** L5 (real data), after E039 has passed

## Question

On the released WMAP nine-year low-resolution Q/U maps with their dense inverse covariances, plus
Planck PR4 353 GHz as the dust band: what is the polarized synchrotron index per sky region under
each estimator (profile likelihood with naive and sandwich errors, flat-prior marginal stated at
23 and at 30 GHz, zero-mean amplitude prior), does the flat-prior marginal move with the reference
frequency by the amount the formula predicts, and how much do the answers move with the band set,
the CMB prior, the bandpass treatment and the dust index?

What changes with the answer: this is the data section of Paper A1. If the reference-frequency
shift appears on real maps at the predicted size, the paper's second claim is demonstrated on
data. If the estimators that are calibrated in simulation disagree with each other on real maps by
more than their errors, the model is wrong somewhere and that is what the section reports.

## Target and priors

As E039: bands K, Ka, Q, V (band-averaged spectra from the WMAP bandpass profiles, dense N^-1 as
released) and 353 GHz (PSB bandpass from the HFI RIMO, see the amendment below; PR4 Q/U covariance combined over sub-pixels); synchrotron power
law at nu_ref = 23 GHz with one index per region, flat on [-4.6, -1.6]; dust modified blackbody
with index 1.53 and 19.6 K, flat amplitude prior; CMB with N(0, 0.2^2 uK^2) per element.

## Setup

- Data: `~/data/wmap9_lowres/`, `~/data/wmap9_lowres_K/` (final nine-year release, files v5),
  `~/data/planck_pr4/HFI_SkyMap_353_2048_R4.00_full.fits`, bandpasses
  `~/data/wmap9_bandpass/`. All official releases (DATA.md, holdings table).
- Baseline: nside 16, the four BeyondPlanck XV regions, all five bands, all three estimators.
- Variations, profile with sandwich error only: nside 8; 48 nside-2 regions; three latitude
  bands; band sets without K, with K and Ka only, with K, Q and V only; flat CMB prior and prior sd
  of 0.1 and 0.4; delta bands at the nominal WMAP frequencies; dust index 1.48 and 1.58.
- No mask. The Galactic Center and Plane regions are kept as BeyondPlanck defined them.

## Reference

There is no truth. The references are (1) E039, which fixes what each estimator does when the
model is right, (2) the formula for the flat-prior shift, and (3) published values, which are
quoted from agent summaries and are not yet checked in the papers: BeyondPlanck XV, Spur -3.17 +-
0.06 and Plane -3.03 +- 0.07 from LFI, WMAP Ka to V and 353 GHz without K band; Dunkley et al.
2009, WMAP five-year K to V, mean -3.03 +- 0.04, low latitude -3.00 +- 0.04, high latitude -3.08
+- 0.06; Fuskeland et al. 2014, T-T plots, -2.98 in the plane to -3.12 at high latitude.

**Expectations, written before any real-data fit:**

- R1. Profile indices at nside 16, BP XV regions: Plane in [-3.10, -2.90], Galactic Center in
  [-3.20, -2.80], Spur in [-3.30, -3.00], High Latitude in [-3.35, -2.95]; the Plane flatter than
  the High Latitude region by 0.05 to 0.30.
- R2. The flat-prior marginal at a 30 GHz reference minus the one at 23 GHz equals the coupled
  prediction J^-1 (n ln(30 / 23)) within 25 percent in every BP XV region.
- R3. The zero-mean prior with a maximized scale is within 0.5 sandwich sd of the profile.
- R4. The four-region model does not fit: chi2_min over the nominal degrees of freedom is above
  1.05, because the index varies inside the regions. With 48 regions it is lower and still above
  1.02. The per-band chi-square is largest in K band.
- R5. Without K band the errors grow by a factor 3 to 8 and each index agrees with the full band
  set within 2.5 combined sd.
- R6. Delta bands at the nominal frequencies make every index steeper by 0.04 to 0.07.
- R7. A flat CMB prior changes the profile indices by less than 1 sandwich sd and increases the
  errors by a factor 1.2 to 2.5.
- R8. Nside 8 and nside 16 agree within 2 sd of the nside 8 error per region.
- R9. A dust index change of 0.05 moves the Plane and Galactic Center indices by less than 0.02
  and the other two by less than 0.005.

**Amendments before any real-data fit, 2026-11-10 21:01.** Two things found while the E039 simulations
were running, neither from real data.

1. *The 353 GHz band is not a delta band.* Checked against the HFI bandpass (RIMO R3.00; the
   reader reproduces the published unit conversions to 0.3 percent and the dust colour
   correction to 0.2 percent), a modified blackbody with index 1.53 and 19.6 K seen through the
   353 GHz band is 1.137 (all bolometers) or 1.121 (polarization-sensitive bolometers) times its
   delta-band value at 353.0 GHz. A delta band therefore overstates the dust amplitude at every
   other frequency by 12 to 14 percent, more than the dust-index variation of R9 (9 percent at V
   band). The baseline now uses the PSB bandpass; the delta band is a sensitivity arm
   (`--bp353 delta`), and all bolometers another. E039 is not affected: its simulations are
   generated and fitted with the same band.
   - R10. Delta band against PSB bandpass at 353 GHz: Plane and Galactic Center indices move by
     less than 0.03, the other two by less than 0.01, and the V-band chi-square per element
     changes by less than 0.02.
2. *The plug-in information is biased.* J evaluated at the fitted amplitudes has expectation J +
   T. The errors quoted for real data use J_c = J(a_hat) - T for the curvature and J_c^-1 J(a_hat)
   J_c^-1 for the sandwich (`fit.sandwich_plugin`; unit test with T / J of 0.1 to 0.5: plain
   plug-in sd of z above 1.08, debiased within 0.08 of 1). For the four BP XV regions T / J is
   below 0.02 in the E039 truth, so R1 to R9 do not change; for small faint regions it matters.
   E039 uses J and T at the truth and tests the theory, not the plug-in; the plug-in is tested in
   E041.

**Second amendment before any real-data fit, 2026-11-10 22:03** (from the independent review of the module
and the E038 beam check; details in the E039 and E038 records).

3. *WMAP bandpass drift.* The team's nine-year correction of the effective frequencies (0.13,
   0.12, 0.11 and 0.06 percent for K to V) is now applied; arm `--no-drift`.
   - R11. Without the drift correction the WMAP-only indices change by less than 0.004.
4. *Window matching.* With each band seen through its own beam a noise-free sky of index -3.1
   returns +0.005 to +0.011 for these bands (E038 c). `--match-window 60` brings every band to a
   common 60 arcmin window to multipole 47 with an invertible operator; the residual is below
   0.0025. By rule D6 of E038 the window-matched fit is the baseline and the unmatched fit an arm.
   - R12. The window-matched indices are lower (steeper) than the unmatched ones by 0.003 to 0.015
     in every region.
5. *Flat-prior shift.* R2 stays as written, against J_c^-1 (n ln(30 / 23)). The corrected
   prediction H^-1 (n ln(30 / 23)), H the Hessian of the marginal at its optimum, is printed
   beside it. With T / J below 0.02 the two differ by less than 4 percent.
6. *CMB prior.* The profile score is unbiased over CMB skies drawn from the prior, not for a
   variance that differs from it. Each run prints the first-order shift per unit of (true CMB
   variance over prior variance minus one). R7 is unchanged.
7. *Errors* are quoted only where T / J_c is below 1.

**Stage 2, Planck LFI and the K band against 30 GHz comparison; expectations written 2026-11-10 21:15,
before any fit that includes an LFI band.** Inputs: the PR4 30, 44 and 70 GHz Q/U maps with the
scale-dependent noise model of E038 (b2), band-averaged spectra from the LFI RIMO. What is already
known from E038 and may colour these expectations: on bright elements the unweighted ratio of the
30 GHz to the K-band map is 0.453 (Q) and 0.470 (U), which read as an index of -3.4 to -3.6
between 22.4 and 28.1 GHz.

- S1. {K, 30, 353}: the profile index is steeper than that of {K, Ka, Q, V, 353} by 0.1 to 0.5 in
  at least two of the four BP XV regions, and the difference exceeds 2 combined sd there.
- S2. {K, Ka, 353}: within 2 combined sd of {K, Ka, Q, V, 353} in every region.
- S3. {30, 44, 70, Ka, Q, V, 353} (no K band, the BeyondPlanck band set on different maps): Spur and
  Plane within 2 combined sd of the published -3.17 +- 0.06 and -3.03 +- 0.07.
- S4. Adding 30, 44 and 70 GHz to {K, Ka, Q, V, 353} moves each region index by less than 2 sd of
  the WMAP-only error and reduces the error by a factor 1.1 to 2.
- S5. With the white-noise columns in place of the noise model the LFI-dominated fits move by
  less than 1 sd and their chi-square per element at 30 GHz rises above 1.5.
- Arms that go with every LFI result: white covariance times the fitted r; the model with its
  power-law term doubled.
- If S1 holds, the discrepancy is the result, and the next step, stated now, is to ask which of
  three one-parameter descriptions removes it: a calibration factor on K band, a relative
  polarization angle between WMAP and LFI, or the removal of the K-band modes whose released
  variance is above ten times the median. None is fitted before S1 is read.

## Command

```bash
cd /Users/trygvels/code/Commander4; S=nextgen/scripts/e040_wmap_region_indices.py; E=nextgen/experiments/E040
.venv/bin/python -I $S --regions bpxv --estimators profile flat zm --out $E/bpxv_n16.npz > $E/bpxv_n16.txt
```

## Result

Pending. Not run before E039 is read.

## Negative control

The flat-prior marginal at two reference frequencies.

## Conclusion

Pending.

## Follow-ups
