# E013: How much did the sign-flipped `marginal` volume term move BeyondPlanck's AME peak frequency? (Q14)

- **Date started / finished:** 2026-10-25 /
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

## Command
```bash
# filled in when run.py exists
```

## Result
(pending)

## Negative control
The code-marginal must reduce to the ridge and marginal results at high S/N (plane pixels); at low
S/N the three must separate in the direction predicted by E009's toy (code-marginal pulled towards
the bound where A is largest).

## Conclusion
(pending)

## Follow-ups
- If the effect is material, repeat with the amplitude redrawn (partially collapsed Gibbs) on a
  small region to confirm the per-pixel approximation does not change the picture.
