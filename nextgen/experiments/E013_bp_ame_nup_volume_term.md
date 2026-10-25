# E013: How much did the sign-flipped `marginal` volume term move BeyondPlanck's AME peak frequency? (Q14)

- **Date started / finished:** 2026-10-25 / 2026-10-25
- **Branch and commit:** nextgen, 0c0c612
- **Open question addressed:** Q14
- **Validation level:** L5-type real-data check with an exact (grid) per-pixel reference; depends on
  E009 and on the independent stats-reviewer pass on the BP10 local sampler (workflow of 2026-10-25)

## Question
On the actual BeyondPlanck data and sampling mask, does the per-pixel posterior of the AME peak
frequency nu_p differ materially between the lnL the production code evaluated (ridge plus
+1/2 ln A), the flat-prior marginal (ridge minus 1/2 ln A), the ridge (profile) likelihood and the
conditional chi-square? If the shift is small relative to the quoted nu_p uncertainty, F10 is a
code-correctness note; if it is large, it is a result that changes published numbers and Trygve must
decide how to report it.

## Scope correction found while scoping
Cosmoglobe DR1 (`CG_param_c0001_v1.txt`, fetched 2026-10-25) modelled AME with the exponential
("exp", power-law-like) parameterization `ame_CG_exp.defaults`, index beta with
`COMP_BETA_INT_LNLTYPE = prior` and `pixreg` regions, i.e. prior draws, not the `marginal` lnL.
E009's "likely DR1" is therefore wrong for DR1; the `marginal` nu_p path is BeyondPlanck only
(`BP_param_c0001_v2.txt`: `COMP_NU_P_INT_LNLTYPE07 = marginal`, `COMP_NU_P_INT_LMAX07 = -1`,
smoothing scale 3, Gaussian prior N(21, 3) GHz, uniform [10, 70] GHz, nu range 15 to 100 GHz,
mask `mask_AME_nu_p_BP8_10deg.fits` / `mask_ame_nup_sampling_temp_bp10_...`). FINDINGS F10 and
OPEN_QUESTIONS Q14 are updated accordingly.

## Target and priors
Per index pixel q at the BP smoothing scale (nside and FWHM from
`COMP_NU_P_SMOOTHING_SCALE07 = 3` in the BP10 smoothing-scale table; to be read from the parameter
file and defaults), with every other component fixed at the BP posterior mean and subtracted from
the smoothed frequency maps, the residual r_{b,q} in band b is modelled as

    r_{b,q} = a_q f_b(nu_p,q) + n_{b,q},   n ~ N(0, sigma_{b,q}^2),

where f_b(nu_p) is Commander3's spinning-dust SED (the `spdust2_cnm.dat` template shifted in
log-frequency to peak at nu_p, converted to uK_RJ with the band's unit conversion, evaluated at
the band centre frequency; bandpass integration as in Commander3 if the bandpass files are
available, otherwise at the centre frequency, stated in the result). With A_q(nu_p) =
sum_b f_b^2 / sigma_b^2 and b_q(nu_p) = sum_b r_b f_b / sigma_b^2 (the sums over the bands in
[15, 100] GHz as BP used), the four per-pixel log densities are

    chisq:        -1/2 sum_b (r_b - a_q^fixed f_b)^2 / sigma_b^2                 (a_q at BP mean)
    ridge:        -1/2 chi2_min = -1/2 [ sum_b r_b^2/sigma_b^2 - b_q^2 / A_q ]
    marginal:      ridge - 1/2 ln A_q                                           (flat prior on a_q)
    code-marginal: ridge + 1/2 ln A_q                                          (what BP10 evaluated)

each plus the Gaussian prior log N(nu_p; 21, 3^2) GHz and the uniform bounds [10, 70] GHz.

## Setup
- Data: BeyondPlanck v2 frequency maps (BP 30/44/70 GHz), WMAP9 Ka, Q1, Q2, V1, V2 as used by BP
  (auxcmd3 inputs), Haslam 408 MHz (below the 15 GHz cut; excluded from the nu_p likelihood as in
  BP), and the HFI inputs that enter the nu_p band set only if any lie in [15, 100] GHz (none do).
  Component means from `BP_*_n1024_v2.fits` (cmb, synch, dust, freefree, ame; the AME file holds
  amplitude and nu_p maps) and the monopole/dipole values from `BP_param_c0001_v2.txt` or
  `BP_mean_v2.h5`.
- Smoothing to the BP local-sampler scale (nside and FWHM from the parameter file), noise rms
  propagated as Commander3 does at that scale (its smoothed rms maps: `BP7_*_rms_n0016_*` and
  `Haslam_rms_n0016_*` in auxcmd3 are these for the 10 deg and 60 arcmin scales; use the file
  matching smoothing scale 3 if present, else propagate white noise and state the approximation).
- Code: `nextgen/experiments/E013_bp_ame_nup_volume_term/run.py` (NumPy; grid evaluation per pixel,
  no MCMC). Machine: Mac mini.

## Reference
The per-pixel posterior for each lnL is one-dimensional and is evaluated exactly on a grid of 2000
nu_p values in [10, 70] GHz; posterior mean, sd and MAP per pixel follow by quadrature. No sampler
is involved, so there is no MC error. The comparison quantities are, over the BP nu_p sampling
mask: the map of (mean_code - mean_marginal), (mean_code - mean_ridge), in units of the per-pixel
posterior sd, and their distribution; the fraction of pixels whose posterior sits on a bound; and
the same for the BP posterior mean nu_p map (the published product) versus each of the four. The
conditioning on fixed other components is an approximation of the Gibbs step and is stated as such;
it isolates the effect of the lnL choice at the S/N of the real pixels, which is what Q14 asks.

## What the independent review changed in the setup
The stats-reviewer pass on the BP10 local sampler (two reviewers, 2026-10-25) established that BP10
sampled nu_p as ONE full-sky value (one pixel region), with the lnL summed over all unmasked
nside-16 pixels, and with the band monopoles profiled (refit deterministically per proposal,
`COMBINED_MONOPOLE_SAMPLING`). The released BP nu_p map is indeed constant (25.60 GHz). The
analysis therefore reports (i) the per-pixel posteriors as a diagnostic of where the volume term
matters and (ii) the full-sky single-parameter posterior, which is the production target. The
production mask (`mask_ame_nup_sampling_temp_bp10_n1024_udgraded_tr75ptc.fits`) is not public;
`mask_AME_nu_p_BP8_10deg.fits` (1204 of 3072 nside-16 pixels) is used instead.

## Command
```bash
python -I nextgen/experiments/E013_bp_ame_nup_volume_term/run.py --out nextgen/experiments/E013_bp_ame_nup_volume_term/out
```

## Result
Summary in `E013_bp_ame_nup_volume_term/summary.txt`, arrays in `out/e013_results.npz`.

Residual check: after subtracting the BP mean sky (without AME) and a per-band monopole/dipole,
the median |residual minus BP AME model| per noise sigma is 4.0 (30 GHz), 2.1 (Ka), 1.7 (44),
1.1 to 1.3 (Q, 70), 0.7 (V). The reconstruction is not perfect (delta bandpasses, Gaussian beams,
BP7 instead of BP10 smoothed rms, refit monopoles), so absolute nu_p values are not comparable with
BP's 25.60 GHz; the differences BETWEEN lnL variants on identical residuals are what this
experiment measures.

Per pixel (1204 mask pixels, 8 bands 28 to 70 GHz), posterior means relative to ridge, in ridge sd:

| variant | median shift | p16 to p84 | median sd [GHz] |
|---|---|---|---|
| code-marginal (ridge + 1/2 ln A) | +0.04 | +0.02 to +0.09 | 0.62 |
| marginal (ridge - 1/2 ln A) | -0.05 | -0.09 to -0.02 | 0.63 |
| chisq (amplitude fixed at BP mean) | -0.09 | -1.28 to +1.24 | 0.21 |

By AME peak S/N (max over bands of a_BP f / sigma): the code-marginal shift is +0.14 (S/N 3 to 10,
2 pixels), +0.11 (10 to 30, 173 px), +0.05 (30 to 100, 778 px), +0.02 (above 100, 242 px) ridge sd;
the correct marginal is the mirror image. No pixel posterior sits on a bound. The AME S/N in the BP
mask is high (median above 30), which is why the volume term is small: E009's toy showed the sign
only matters once many low-S/N pixels enter the sum.

Full-sky single nu_p (the production target), fine grid of 4001 points:

| variant | mean [GHz] | sd [GHz] | shift vs ridge |
|---|---|---|---|
| ridge | 26.658 | 0.0098 | 0 |
| marginal | 26.650 | 0.0096 | -0.76 sd |
| code-marginal | 26.665 | 0.0099 | +0.78 sd |
| chisq | 26.311 | 0.0037 | -35 sd |

The sign error moved the full-sky conditional posterior of nu_p by 0.008 GHz (0.8 of its
conditional sd of 0.01 GHz); the correct marginal would have moved it by the same amount the other
way. The conditional sd here is far smaller than any realistic nu_p uncertainty (the model
reconstruction offset alone is 1 GHz), so the effect on the published value is negligible. The
chisq variant (amplitude fixed) is 0.35 GHz away, which shows that conditioning on the amplitude
matters far more than the volume term for this parameter.

## Negative control
The code-marginal must reduce to the ridge and marginal results at high S/N (plane pixels); at low
S/N the three must separate in the direction predicted by E009's toy (code-marginal pulled towards
the bound where A is largest).

## Conclusion
Q14 for BeyondPlanck: the `marginal` sign error (F10) is a real code error but had a negligible
effect on the BP nu_p because the AME pixels in the sampling mask have high S/N (median above 30)
and nu_p was a single full-sky parameter: the shift is 0.8 of a 0.01 GHz conditional sd. It is a
correctness note for Paper A, not a correction to a published number. For DR1 the question does
not arise (exponential AME model with prior-only index steps). The DR2-era `sample_specind_mh`
usage remains unverified (no parameter files in the repository use MCMC groups).
Negative control: at S/N 3 to 10 the shift grows to 0.14 sd and the two signs mirror each other,
as E009 predicted; chisq separates strongly from the collapsed variants everywhere.

## Follow-ups
- If the effect is material, repeat with the amplitude redrawn (partially collapsed Gibbs) on a
  small region to confirm the per-pixel approximation does not change the picture.
