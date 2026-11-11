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

**Stage 3, the BeyondPlanck band set and reference frequency; expectations written 2026-11-10 22:50, before
any fit with the flat-prior marginal on a band set that includes LFI.** What prompts it: the second
citation check (`../research/2026-11-10_citation_verification_paperA1_round2.md`, item 6c; round 1,
item 10) quotes Svalheim et al. (2023, Sect. 4.4) and Watts et al. (2023, Sect. 2.4): fitted
freely, the high-latitude polarized synchrotron index goes to "unreasonably flat values, beta_s
<~ -2.5", read there as unmodelled systematics and handled with a prior N(-3.3, 0.1^2); a wider
prior "permits synchrotron indices as shallow as beta_s = -2.7" and biases the CMB. Their index
step is a Metropolis step conditional on the amplitudes inside a Gibbs sweep, whose index marginal
is the flat-prior marginal (F15), with the synchrotron amplitude defined at 30 GHz. E022 put the
pivot of that band set at 26 GHz with a free CMB amplitude and at 30 GHz with a CMB prior of 0.19
uK. So the flat-prior volume term predicts a pull toward flat indices there when the CMB is free.
This stage tests that reading on PR4 and WMAP maps at nside 16. It does not reproduce their
pipeline, their maps or their resolution.

Runs: bands 30, 44, 70, Ka, Q, V, 353 (no K band), four BP XV regions, window-matched, LFI noise
model of E038; profile and flat-prior marginal at references of 23 and 30 GHz; once with a flat
CMB prior and once with the 0.2 uK prior. Script `nextgen/scripts/e040b_estimators.py` (a
memory-lean copy of the E040 script; same numerics).

- B1. Flat CMB prior, reference 30 GHz: the flat-prior optimum is flatter than the profile
  optimum in all four regions, by more than 0.3 in the High Latitude region (or it sits on the
  upper bound of the allowed range), by 0.1 to 0.5 in the Spur and by 0.05 to 0.3 in the Plane.
- B2. CMB prior of 0.2 uK, reference 30 GHz: the flat-prior optimum is within 0.05 of the profile
  optimum in the Spur and the Plane.
- B3. The profile index of the High Latitude region on this band set is between -3.4 and -2.9
  with either CMB prior, and its sandwich error is 1.5 to 4 times that of the band set with K.
- If B1 and B2 hold, the paper says that the flat pull reported for this configuration has the
  sign and the order of magnitude of the volume term of a flat amplitude prior stated at 30 GHz
  with an unconstrained CMB amplitude, and that at a finer pixelization it is larger in
  proportion to the number of amplitudes. It does not say that this is what happened in those
  analyses: that needs their chains. If B1 fails, the reading is withdrawn.

**Stage 4, after the referee report on draft 1 of the paper (`../reviews/2026-11-11_referee_report_paperA1_draft1.md`);
expectations written 2026-11-11 02:43 before these fits.** The referee's check: a white CMB prior of 0.2 uK
per element understates the LCDM E-mode variance by factors of 199, 128, 66, 29 and 12 at
multipoles 2 to 6 and overstates it above multipole 10 and for every B mode. The central table of
the paper (the four estimators on K, Ka, Q, V, 353, window-matched, nside 16) is repeated with
(i) the LCDM pixel covariance as the CMB prior (`lowres.window.lcdm_cmb_precision`: CAMB, Planck
2018 parameters, tau 0.0544, the 60 arcmin window and the pixel window, a 0.01 uK white floor;
checked: rms 0.17 uK per element, per-mode variance 7.95 uK^2 at multipole 2) and (ii) a flat CMB
prior, for which the profile score is unbiased for any CMB.

- C1. LCDM prior: each profile index is within 1 sandwich sd of the white-prior value, and the
  move of the flat-prior optimum from a 23 to a 30 GHz reference is within 15 percent of the
  white-prior move (+0.158, +0.075, +0.025, +0.037).
- C2. Flat CMB prior: each profile index is within 1.5 sd (of the flat-prior error) of the
  white-prior value, the errors are 1.2 to 2.5 times larger (R7), and the flat-prior optimum is
  flatter than the profile at both references (the pivot moves to or below the lowest band), by
  more at 30 GHz than at 23 GHz; the difference between the two references is larger than with
  the white prior.
- C3. In all three versions the difference between the references agrees with H^-1 (n ln(30 / 23))
  within 25 percent wherever it is below 0.3.
- If C1 fails, the white prior was shaping the baseline and the LCDM version replaces it in the
  paper's table.

## Command

```bash
cd /Users/trygvels/code/Commander4; S=nextgen/scripts/e040_wmap_region_indices.py; E=nextgen/experiments/E040
.venv/bin/python -I $S --regions bpxv --estimators profile flat zm --out $E/bpxv_n16.npz > $E/bpxv_n16.txt
```

## Result

E039 stage 1 (four regions) was read and recorded at 22:03 on 2026-11-10; the first real-data fit
started at 22:04. The 48-region arm waits for E039's 48-region simulation.

**Nside 8, four BP XV regions, K, Ka, Q, V, 353, window-matched, every estimator** (`E040/bpxv_n8_wm_all.txt`,
22:04 to 22:06; this is arm R8 and also the first execution of the script on real maps).

| region | elements | profile | sandwich sd | flat23 - profile | flat30 - profile | flat30 - flat23 | predicted J_c^-1 | predicted H^-1 | zero-mean - profile |
|---|---|---|---|---|---|---|---|---|---|
| High Latitude | 1124 | -3.0709 | 0.0145 | -0.0012 | +0.0605 | +0.0618 | +0.0652 | +0.0642 | +0.0004 |
| Spur | 212 | -3.1305 | 0.0155 | -0.0003 | +0.0247 | +0.0251 | +0.0258 | +0.0256 | +0.0004 |
| Galactic Center | 50 | -2.8552 | 0.0112 | -0.0002 | +0.0047 | +0.0049 | +0.0054 | +0.0050 | 0.0000 |
| Plane | 150 | -3.0615 | 0.0143 | -0.0002 | +0.0123 | +0.0125 | +0.0128 | +0.0127 | +0.0001 |

Chi-square at the optimum 7155.9 for 4604 nominal degrees of freedom, ratio 1.554. Per band and
element: K 0.31, Ka 1.89, Q 1.38, V 0.96, 353 GHz 0.00. By region, per element and band: High
Latitude 0.72, Spur 0.84, Galactic Center 5.78, Plane 0.80. T / J_c at most 0.005. CMB prior: per
unit of (true variance over prior variance minus one) the index moves by +0.16 sd (High Latitude)
and at most +0.06 sd elsewhere.

Read against the expectations (nside 8 only; the nside-16 baseline follows):

- **R2 holds on real maps.** The flat-prior optimum moves between a 23 and a 30 GHz reference by
  +0.062, +0.025, +0.005 and +0.013, which is 0.95, 0.97, 0.91 and 0.98 of the pre-registered
  prediction and 0.96, 0.98, 0.98 and 0.98 of the corrected one. For the High Latitude region
  that is 4.2 sandwich sd from a change of convention. At 23 GHz the flat prior is neutral for
  this band set (within 0.08 sd).
- **R3 holds**: the zero-mean prior is within 0.03 sd of the profile optimum.
- **R4: the misfit is there (1.55), the band is not the one I named.** I expected the largest
  per-band chi-square at K band; it is the smallest (0.31), because K band fixes the synchrotron
  amplitude, and Ka is the largest (1.89). The Galactic Center region carries a per-element
  chi-square seven times that of the other three.
- **R1, in part**: all four indices are inside their ranges. The Plane is flatter than the High
  Latitude region by 0.009 and not by the 0.05 to 0.30 I expected; the flat region is the
  Galactic Center (-2.855).

**Nside 16 baseline: four BP XV regions, K, Ka, Q, V, 353, window-matched, every estimator**
(`E040/bpxv_n16_wm.txt`, 22:06 to 22:44; 38 min because the script kept four likelihood objects
alive and the machine swapped).

| region | elements | profile | sandwich sd | flat23 - profile | flat30 - profile | flat30 - flat23 | predicted J_c^-1 | predicted H^-1 | zero-mean - profile |
|---|---|---|---|---|---|---|---|---|---|
| High Latitude | 4388 | -3.0940 | 0.0121 | -0.0077 | +0.1498 | +0.1575 | +0.1758 | +0.1721 | +0.0051 |
| Spur | 928 | -3.0965 | 0.0149 | -0.0034 | +0.0715 | +0.0749 | +0.0804 | +0.0789 | +0.0005 |
| Galactic Center | 182 | -3.0867 | 0.0088 | -0.0012 | +0.0242 | +0.0254 | +0.0273 | +0.0267 | +0.0090 |
| Plane | 646 | -2.9950 | 0.0120 | -0.0020 | +0.0345 | +0.0365 | +0.0379 | +0.0378 | +0.0029 |

Chi-square at the optimum 19792.7 for 18428 nominal degrees of freedom, ratio 1.074. Per band and
element: K 0.14, Ka 1.01, Q 1.06, V 0.99, 353 GHz 0.00. By region, per element and band (0.6 is
what three free amplitudes on five bands leave): High Latitude 0.58, Spur 0.59, Plane 0.78,
Galactic Center 1.82. T / J_c at most 0.013. CMB prior sensitivity +0.15 sd (High Latitude) per
unit of variance mismatch, at most +0.06 sd elsewhere.

- **R1 holds at nside 16**: Plane -2.995, Galactic Center -3.087, Spur -3.097, High Latitude -3.094,
  all inside their ranges, and the Plane is flatter than the High Latitude region by 0.099.
- **R2 holds.** Between a 23 and a 30 GHz reference the flat-prior optimum moves by +0.158, +0.075,
  +0.025 and +0.037: 0.90, 0.93, 0.93 and 0.96 of the pre-registered prediction and 0.92, 0.95,
  0.95 and 0.97 of the corrected one. In sandwich sd of the profile that is 13.0, 5.0, 2.9 and
  3.0. The same maps, regions, bands and noise model give a High Latitude index of -3.094
  (profile), -3.102 (flat prior at 23 GHz) or -2.944 (flat prior at 30 GHz). At nside 8, with a
  quarter of the amplitudes, the same move was +0.062.
- **R3 holds in three regions and fails in the Galactic Center**: zero-mean minus profile is
  +0.42, +0.03 and +0.24 sd, and +1.02 sd in the Galactic Center, the region that does not fit.
- **R4**: 1.074, above the 1.05 I expected. The misfit sits in the Galactic Center (three times
  the expected chi-square) and the Plane (1.3 times). The High Latitude and Spur regions fit at
  0.58 and 0.59 against an expected 0.6. K band is not the worst band, it is the best (0.14).
- **R8 fails, for a reason that is now understood.** Nside 8 against nside 16: High Latitude
  +0.023, Spur -0.034, Plane -0.067, Galactic Center +0.232 (1.6, 2.2, 4.7 and 21 sd of the
  nside-8 error). A noise-free sky with index -3.1 everywhere, degraded from nside 16 to nside 8
  with the real inverse covariances as the nside-8 arm does, returns +0.006, -0.018, -0.019 and
  +0.217 (`E038/window_wmap_matched_n8.txt`; the same with a common beam, so it is not the
  beams). The generalized-least-squares degrade assumes the signal is constant inside a coarse
  pixel, the dense weights differ between bands, and in the Galactic Center the signal is
  anything but constant. **Nside-8 fits of the real maps are not a valid cross-check and are
  discarded.** The E039 simulations are not affected: their truth is constant inside the nside-8
  pixels by construction. The same mechanism at a smaller scale, sub-pixel weights that differ
  between bands inside an nside-16 pixel, cannot be excluded for the released maps themselves and
  is the first suspect for the Galactic Center misfit.

**Window arm** (`E040/bpxv_n16.txt`, unmatched, profile only): -3.0873, -3.0940, -3.0782, -2.9915. The
window-matched indices are steeper by 0.0067, 0.0025, 0.0085 and 0.0035. **R12 holds in three
regions** and misses in the Spur by 0.0005 (0.0025 against a lower bound of 0.003). The noise-free
simulation of E038 (c) predicted that matching removes 0.0070, 0.0030, 0.0088 and 0.0041: the real
maps agree with it to 0.0006 in every region.

**Stage 2, first pairs** (window-matched, profile, debiased sandwich; `E040/bpxv_n16_KKa_wm.txt`,
`s2_K30_wm.txt`; 22:48 to 22:55).

| bands | High Latitude | Spur | Galactic Center | Plane | chi-square per element: K, second band |
|---|---|---|---|---|---|
| K, Ka, Q, V, 353 (baseline) | -3.094 +- 0.012 | -3.097 +- 0.015 | -3.087 +- 0.009 | -2.995 +- 0.012 | 0.14, Ka 1.01 |
| K, Ka, 353 | -3.062 +- 0.016 | -3.074 +- 0.020 | -3.046 +- 0.011 | -2.980 +- 0.016 | 0.10, Ka 1.02 |
| K, 30, 353 | -3.428 +- 0.019 | -3.180 +- 0.023 | -3.263 +- 0.013 | -3.035 +- 0.020 | 0.78, 30 GHz 1.85 |

- **S1 holds.** With K band and 30 GHz alone the index is steeper than the WMAP baseline by 0.334
  (High Latitude), 0.083 (Spur), 0.176 (Galactic Center) and 0.040 (Plane); two regions are in the
  0.1 to 0.5 band and beyond 2 combined sd (15 and 11.5), the Spur is at 3.0 sd and the Plane at
  1.7. Against K with Ka the gap is 0.365, 0.106, 0.217 and 0.055. In amplitude the High Latitude
  gap is 8.2 percent between two maps 0.2245 apart in ln nu.
- **S2 holds in three regions and fails in the Galactic Center** by the combined-sd measure I wrote
  down (K, Ka against the baseline: 1.6, 0.9, 2.8 and 0.8). The measure was too lenient: the two
  fits share K and Ka, so the proper sd of the difference is sqrt(sd_sub^2 - sd_full^2), and by
  that the High Latitude difference of +0.032 is 3.3 sd and the Galactic Center one 5.7 sd. Adding
  Q and V steepens every region by 0.015 to 0.041. Candidates: a spectrum that steepens with
  frequency, or the dust extrapolation from 353 GHz to V band; the dust-index arms bear on the
  second.
- The reduced chi-square of the K, 30 fit is 2.64, carried by the Galactic Center (8.6 per
  element and band) and the Plane (1.7); the High Latitude and Spur regions are at 0.52 and 0.50.
  At high latitude the two maps are each consistent with a power law and the power law is the
  wrong one by 0.37: this is what a multiplicative or a smooth large-scale difference between the
  maps does, and what a random mismatch does not do. E042 takes it apart.

**Stage 2, the band set without K band** (`E040/s2_noK_lfi_wm.txt`, 22:55 to 22:59; 30, 44, 70, Ka, Q, V,
353; profile, CMB prior 0.2 uK): High Latitude -2.994 +- 0.019, Spur -3.059 +- 0.024, Galactic
Center -3.048 +- 0.012, Plane -2.766 +- 0.018. Reduced chi-square 1.28; per element 0.53 (30 GHz),
1.40 (44), 1.27 (70), 1.12 (Ka), 1.02 (Q), 0.99 (V); by region 0.73, 0.74, 3.25, 1.68. CMB prior
sensitivity +0.54 sd per unit of variance mismatch at high latitude, four times that of the
WMAP-only set.

- **S3 holds for the Spur and fails for the Plane.** Against the published -3.17 +- 0.06 and
  -3.03 +- 0.07: the Spur is 1.7 combined sd away, the Plane 3.7 (0.26 flatter). The published
  analysis used its own 30 to 70 GHz maps, not PR4, and a prior N(-3.3, 0.1^2).
- **B3 of stage 3 holds**: the High Latitude profile index on this band set is -2.994, and its
  error is 1.56 times that of the set with K band.
- With 30 GHz as the lowest band the Plane is at -2.77; with K band in any combination it is
  between -2.98 and -3.04. Together with the K, 30 fit this says the PR4 30 GHz map is faint
  against both neighbours. Note that the PR4 44 and 70 GHz maps are calibrated with a polarized
  prior built on 30 GHz (NPIPE Sect. 4.3; second citation report, item 4b), so the three LFI
  bands are not independent witnesses on large scales.

**Stage 3, the BeyondPlanck band set with a flat CMB prior** (`E040/s3_bpset_cmbflat_wm.txt`, 23:43 to
00:15; 30, 44, 70, Ka, Q, V, 353, no K band, flat priors on every amplitude; a first attempt
crashed after 32 minutes at its last print statement on a bug in my memory-lean copy of the
script, which was fixed and exercised end to end at nside 8 before the rerun).

| region | profile | sandwich sd | flat prior at 23 GHz | flat prior at 30 GHz | flat30 - profile | flat30 - flat23 |
|---|---|---|---|---|---|---|
| High Latitude | -2.695 | 0.040 | -4.036 | -1.600 (upper bound) | +1.095 | +2.436 |
| Spur | -2.921 | 0.050 | -3.432 | -2.314 | +0.607 | +1.118 |
| Galactic Center | -2.819 | 0.025 | -2.885 | -2.752 | +0.067 | +0.132 |
| Plane | -2.780 | 0.042 | -2.997 | -2.479 | +0.302 | +0.518 |

Reduced chi-square 1.18; by region 0.58, 0.59, 2.13, 1.03. T / J_c 0.07 at high latitude.

- **B1 holds in direction in all four regions and in size at high latitude; the Spur and the Plane
  exceed my upper bounds** (+0.61 against 0.5 and +0.302 against 0.3). With every amplitude under a
  flat prior and the synchrotron amplitude defined at 30 GHz, the optimum of the index marginal
  for the High Latitude region is at the upper bound of the allowed range, 1.1 above the profile
  optimum; the Spur is at -2.31 and the Plane at -2.48. Stated at 23 GHz the same marginal puts
  the High Latitude region at -4.04. The pivot of this band set with a free CMB lies between the
  two, as E022 found (26 GHz).
- The first-order prediction is useless at this size (it gives +5.3 for a move that the bound
  stops at +2.4), as Fig. 3 of the paper already shows for shifts above 0.3.
- **The profile optimum itself depends on the CMB prior for this band set**: -2.695 +- 0.040 with a
  flat CMB prior against -2.994 +- 0.019 with 0.2 uK (High Latitude), -2.921 against -3.059
  (Spur), -2.780 against -2.766 (Plane). With a flat prior the optimum is unbiased for any
  achromatic component; with the 0.2 uK prior it is unbiased only if what is achromatic in the
  maps has that variance. A 7 sd move at high latitude says these seven maps contain more
  achromatic structure than a 0.2 uK CMB, which is another view of the 30 GHz problem of E042.
- What can and cannot be said about the published analyses: in this configuration (their band
  set and reference, PR4 and WMAP maps, nside 16) a conditional index step inside a Gibbs sweep
  with flat amplitude priors has a marginal whose optimum is "unreasonably flat" at high latitude
  by more than a unit of index, from the volume term alone, with no systematics invoked. That is
  the behaviour Svalheim et al. (2023) and Watts et al. (2023) report and constrain with a prior
  N(-3.3, 0.1^2). Whether it is what happened in their chains needs their chains: the maps, the
  pixelization, the CMB treatment and the mixing differ.

### Stage 4 result (i): LCDM CMB prior, read 2026-11-11 03:24

`E040/s4_bpxv_n16_wm_lcdm.txt` (K, Ka, Q, V, 353, window-matched, nside 16, every estimator; chi2_min 19767.7
against 19792.7 with the white prior, ratio to the nominal degrees of freedom 1.0727; T / J_c at most 0.014).
The run started at 02:46 on the Cholesky path that E043 later found faulty at smaller orders; order 18432 is
not exposed (E043, V1 and V2 check this fit's baseline directly).

| Region | Profile, white 0.2 uK | Profile, LCDM | Move | In sd | flat23 - profile | flat30 - profile | 30 - 23 | Predicted (H) | zero-mean - profile |
|---|---|---|---|---|---|---|---|---|---|
| High Latitude | -3.0940 +- 0.0121 | -3.0825 +- 0.0125 | +0.0115 | 0.9 | -0.0098 | +0.1576 | +0.1674 | +0.1844 | -0.0077 |
| Spur | -3.0965 +- 0.0149 | -3.0799 +- 0.0165 | +0.0166 | 1.0 | -0.0053 | +0.0896 | +0.0949 | +0.1016 | -0.0130 |
| Galactic Center | -3.0867 +- 0.0088 | -3.0807 +- 0.0088 | +0.0060 | 0.7 | -0.0016 | +0.0272 | +0.0287 | +0.0304 | +0.0061 |
| Plane | -2.9950 +- 0.0120 | -2.9641 +- 0.0133 | +0.0309 | 2.3 | -0.0032 | +0.0452 | +0.0485 | +0.0509 | -0.0079 |

Against the expectations:

- **C1 fails.** The Plane index moves by 2.3 sd and the Spur by 1.01 sd (limit 1 sd); the move between the
  references changes by +6, +27, +13 and +33 percent (limit 15 percent). The white prior was shaping the
  baseline, most where the CMB is least negligible against the noise per mode. As written in advance, the
  LCDM version replaces the white-prior version in the paper's table, and the white-prior numbers are kept
  as an arm.
- C3 holds for this version: the move between the references is 5 to 9 percent below H^-1 (n ln(30 / 23))
  in all four regions.
- The direction: all four indices are flatter with the LCDM prior, by 0.006 to 0.031. A prior that allows
  the CMB its large-scale E modes takes large-scale achromatic power out of what the white prior forced
  into synchrotron and noise.
- The zero-mean prior now sits below the profile in three regions (-0.008 to -0.013) where it sat above
  with the white prior (+0.001 to +0.005): 0.5 to 0.8 sd. The noise patterns of the bands differ, so the
  identity with the profile optimum does not apply (E041), and the sign is not predicted.
- Every all-scales number quoted from the white-prior baseline (band pairs, sensitivity arms, E042
  low-pass arms) is a white-prior number and is labelled so; the arms that enter the paper are repeated
  with the LCDM prior where they are quoted as indices and kept as differences where they are quoted as
  differences.

### Stage 4, check C4, pre-registered 2026-11-11 03:45: is the move between the CMB priors what a LCDM sky does?

C1 failed: with the LCDM prior the four indices are flatter than with the white prior by 0.012, 0.017, 0.006
and 0.031. Two readings. (a) The sky has a LCDM CMB, the white prior has the wrong shape, and the white-prior
fit is biased steep by about that much; then the LCDM version is the better number. (b) The LCDM prior
gives large-scale achromatic power somewhere to go, and what it absorbs in the real maps is not CMB (the
CMB E-mode signal at multipoles 2 to 6 is far below the WMAP noise per mode); then neither version is
"right" and the difference is a size of a systematic. A simulation separates them:
`nextgen/scripts/e039d_cmbprior_sim.py`, the E039 simulations with the CMB drawn from the LCDM covariance,
fitted with the LCDM prior, the white prior and a flat prior; nside 8 first (60 realizations), then nside 16
with a few realizations if the answer depends on resolution.

- C4a. LCDM prior on LCDM skies: unbiased (mean z within 2 standard errors of zero per region) and calibrated
  (sd of z 0.85 to 1.15).
- C4b. Flat CMB prior: unbiased, errors 1.2 to 2.5 times larger (R7).
- C4c. White prior on LCDM skies: I expect reading (b), held weakly: the LCDM E-mode power is small
  against the noise per mode at these multipoles, so the paired difference between the white-prior and
  LCDM-prior estimates should be below 0.003 in absolute value in every region, which is five to ten times
  smaller than the move on the real maps. If instead the simulation shows the real-map move (white steeper
  by 0.005 to 0.03), reading (a) holds.
- Consequence for the paper: under (a) the LCDM version is the central table and the white version is a
  cautionary arm. Under (b) the table shows both and the text says that the difference measures
  large-scale achromatic structure in the maps that is not LCDM, by region, and that the index is uncertain
  by that amount on top of everything else.

### Stage 5, pre-registered 2026-11-11 03:20: errors that do not assume the model fits (referee item 8)

The sandwich error assumes that the model is right inside each region: one index, the stated noise. The
per-region chi-square says it is close to right at high latitude and in the Spur (0.97 and 0.99 of the
nominal three degrees of freedom per element), less so in the Plane (1.29) and not in the Galactic
Center (3.0). `nextgen/scripts/e040d_block_errors.py` computes a block-robust error for the baseline fit:
the per-element score terms g_p = -sum_b (N_b^-1 r_b)_p (ds_b a_s)_p at the optimum are summed over the
48 nside-2 pixels (Q and U of a pixel in the same block), K_block = sum_B G_B G_B^T, and the covariance is
J_c^-1 K_block J_c^-1. Blocks are treated as independent and no small-sample factor n_B / (n_B - 1) is
applied (under 5 percent on the error for ten blocks or more).

Expectations, written before the script has run on the maps:

- D1. On simulations where the model is true the ratio of the block error to the sandwich error is 1,
  with scatter about 1 / sqrt(2 (n_B - 1)) for n_B blocks. Checked first, on the E039 nside-8 simulations
  (`e040e_block_errors_sim.py`); if the ratio is not 1 there the estimator is wrong and is not used.
- D2. On the maps I expect the ratio to exceed 1 in every region, because the index varies inside the
  regions: for a spread of 0.05 to 0.1 between blocks and the measured information per block, 1.2 to 2
  at High Latitude, 1 to 2 in the Spur, 1.5 to 3 in the Plane. For the Galactic Center (three or four
  blocks) no block error is quoted, whatever the number.
- D3. If the ratio is above 1.5 anywhere, the central table of Paper A1 carries both errors and says
  which statement each supports: the sandwich error is the error of the estimator under the model, the
  block error is the error of "the index of this region" when the region is not one index.
- This is the error of the estimate under misfit. It does not make a misfitting model right, and it is
  not a measurement of the spread of the index.

## Negative control

The flat-prior marginal at two reference frequencies.

## Conclusion

Pending.

## Follow-ups
