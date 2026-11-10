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
released) and 353 GHz (delta band, PR4 Q/U covariance combined over sub-pixels); synchrotron power
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
