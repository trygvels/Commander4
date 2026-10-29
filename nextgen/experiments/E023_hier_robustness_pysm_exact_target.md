# E023: Robustness of the hierarchical index target to template and SED misspecification on the PySM sky, by exact-target scans

- **Date started / finished:** 2026-10-29 /
- **Branch and commit:** nextgen, fe50735
- **Open question addressed:** Q16 (robust amplitude-template prior and dust SED at high S/N); feeds M2.5 (PySM skies) and the real-data model (M3)
- **Validation level:** L3 (misspecified simulation with known truth), stage 1 exact target only

## Question
On the PySM d1s1c1 sky at realistic noise (E012 levels, noise scale 1) and at the E020 level
(0.03), how far does the exact collapsed hier target (the likelihood Commander4 samples under
`likelihood: hier`) move the index field away from the truth in each S/N quartile, which of the
two misspecifications (data-built template under one global scatter; fixed T_d) drives it, and
does a per-quartile scatter remove it? If the target is unbiased at noise scale 1 with the
present prior, Paper A's method can run on real intensity maps at that depth and the robust
prior is a refinement; if it is biased, the prior change is required before any real-data run
(Q16 "Decides").

## Target and priors
The stage-1 object is the collapsed per-pixel marginal that Commander4's `hier` likelihood
evaluates, as a function of the index fields alone, with the amplitudes integrated out and the
amplitude hyperparameters at their joint MAP (profile) for each index state:

  L(beta) = max over (mu_s, s_s, mu_d, s_d) of
            sum_p [ -1/2 r_p^T C_p^-1 r_p - 1/2 log det C_p ] + log p(mu_s) + log p(log s_s)
            + log p(mu_d) + log p(log s_d),

with, per pixel p and the 12 bands, M_p = [m_cmb, m_sync(beta_s,p), m_dust(beta_d,p)] (uK_RJ per
unit amplitude; CMB in uK_CMB, sync a power law referenced at 30 GHz, dust one MBB referenced at
545 GHz at the fixed T = 20.93 K, or at the PySM T_d,p where stated), prior means m_p = (0,
mu_s g_s,p, mu_d g_d,p), prior variances T_p = diag(100^2, (s_s g_s,p)^2, (s_d g_d,p)^2), C_p =
M_p T_p M_p^T + N_p (white, the rms maps), r_p = d_p - M_p m_p, and the E021 hyperpriors mu_c ~
N(1, 0.5^2), log s_c ~ N(log 0.3, 1^2). This is the E018/E021 target (stats-reviewed there)
minus the index prior and with the hyperparameters profiled instead of marginalized. The
hyperparameter posteriors have sd about 0.002 on 49,152 pixels (E021), so profile and marginal
agree to that order; the index prior is left out on purpose: the question is whether the data
term, which any index prior can only shrink, is biased.

Candidate prior (stage 1, option a of Q16): the global scatter s_c replaced by one scatter per
S/N quartile of the component, s_c,q, q = 1..4 (eight hyperparameters profiled with the same
hyperprior each), everything else unchanged.

Statistic: for each component c, each S/N quartile Q (the sync or dust Fisher-sd quartile of the
nside-16 regions, as in E021), the block shift delta that maximizes L(beta_truth + delta 1_Q)
over delta in [-0.1, 0.1] in steps of 0.01, plus a wide sync scan over [-3, +1] in the faintest
quartile for the E020 runaway; the curvature sd of the block at the peak; the split into the
chi-square and log-determinant parts. "Truth" for dust is beta_eff3 (the fixed-T effective index
of the nside-16 region, as in the E021 analyzer) when T is fixed, and the region-mean beta_d when
the dust SED carries the PySM T_d,p.

Factor cells (2 x 2 x 2):
- noise scale 1 (new sky, same seed 2027, `--sky pysm --noise-scale 1`) and 0.03 (the E020 sky);
- templates: data-built (the E020 pseudo-templates, g/true 0.41 to 1.66 for sync at |b| > 20)
  and the model-story templates (g = PySM amplitude, a = g r with r = max(1 + 0.3 eps, 0.2),
  for which the hier prior is exact per F23);
- dust SED: fixed T = 20.93 K and the PySM T_d,p (the no-SED-misfit reference).

## Setup
Skies at nside 64, 12 bands 23 to 353 GHz, white noise, from `nextgen/scripts/pysm_intensity_sim.py`:
`/Users/trygvels/data/e020_d1s1` (noise 0.03, E020) and `/Users/trygvels/data/e023_pysm_n1` (noise
scale 1, this experiment). Model-story templates for the PySM sky: built in the scan script from
`truth.npz` (a_sync, a_dust) with the same r law and seed as E021 attempt 2. Scan script:
`nextgen/scripts/e023_exact_scan.py` (numpy 12 x 12 marginal, the E021 `e021_iid_check.py` kernel
extended with per-pixel T_d and per-quartile scatter; hyperparameters profiled with L-BFGS from the
E021 realized values). Runs on the Mac mini alongside the E021 seed-2028 queue (one core).

## Reference
The truth index maps of the PySM sky (`truth.npz`: beta_s, beta_d, T_d, beta_eff3). Expected
behaviour where the model is exact: on the model-story, true-T cell the preferred shift must be
zero within the block curvature sd at both noise levels (this is the F23 model-story check on a
new sky, and the control that the scan code is right). The E020 chain gives the expected sign and
size at 0.03 with data-built templates and fixed T: sync shifts of order -1 or worse in the faint
half, dust +0.08 in the plane.

Pre-registered reading (stage 1):
- A cell is "unbiased at this depth" when every quartile's preferred shift is below 0.02 (sync)
  and 0.01 (dust), smaller than the E021 region-level noise error (0.033 / 0.018) at the same
  nside, and the wide sync scan peaks within 0.1 of the truth.
- Attribution: compare the data-built/fixed-T cell with the two single-misspecification cells;
  the misspecification whose removal brings the cell to "unbiased" is the driver.
- The per-quartile scatter is adopted for stage 2 only if it brings the data-built cell to
  "unbiased" at noise scale 1 without degrading the model-story cell.
- A shift that survives both the true-T SED and the model-story templates is a defect in the scan
  or the sky, not a prior result; stop and find it.

Stage 2 (chains, after the E021 seed-2028 queue): gmrf16 hier with the prior stage 1 selects on
the noise-1 PySM sky, E021 metrics and the Q16 acceptance (index RMSE within two bootstrap SE of
the E021 consistent-sky result per quartile; no region with |beta_s - truth| above 0.5). Written
up as a second stage of this record or as E024 if the prior changes the code.

## Command
```bash
python nextgen/scripts/pysm_intensity_sim.py --dir /Users/trygvels/data/e023_pysm_n1 --nside 64 --noise-scale 1 --seed 2027 --sky pysm --n-iter 1200 --n-warmup 400 --refresh-mass-at 50 200 --no-adapt-mass
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e020_d1s1 --out nextgen/experiments/E023/scan_n003.txt
python nextgen/scripts/e023_exact_scan.py --dir /Users/trygvels/data/e023_pysm_n1 --out nextgen/experiments/E023/scan_n1.txt
```

## Result
(pending)

## Negative control
The E020 cell (noise 0.03, data-built templates, fixed T) must reproduce F21: a sync runaway in
the faint quartiles (the wide scan's maximum far from the truth, below -5 or near -0.6) and a
positive dust shift in the brightest quartile. If the scan does not show it, the scan is wrong.

## Conclusion
(pending)

## Follow-ups
