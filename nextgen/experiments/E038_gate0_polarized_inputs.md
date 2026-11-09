# E038: Gate 0 of the polarized track, are the public low-resolution inputs usable as released?

- **Date started / finished:** 2026-11-09 23:48 / running (item a.1 done; a.2, b and c open)
- **Branch and commit:** `nextgen` at 0413181
- **Open question addressed:** Q21; Gate 0 of `../reviews/2026-11-09_way_forward_challenge.md` Sect. F
- **Validation level:** L0 (input checks, no inference)

## Question

The proposed K band against 30 GHz null test and everything after it assume that the WMAP
nine-year res 4 Q/U maps and inverse covariances are usable as released, that a trustworthy LFI
30 GHz low-resolution covariance exists publicly, and that WMAP and LFI can be put through one
window. Are they? What changes: if (b) fails, LFI enters with an inflated diagonal and "exact" is
dropped for it; if (c) fails, one window for both sets comes before any index.

## Target and priors

None. Deterministic checks.

## Setup

- (a.1) For each WMAP band: symmetry of N^-1, its spectrum and null space, the pixel rms implied
  by its inverse against the white-noise diagonal, how far the covariance is from white, and the
  chi-square of the map. Script `nextgen/scripts/gate0_wmap_ninv_check.py`. K band fetched
  2026-11-09 23:43 from LAMBDA (`data/map/dr5/skymaps/9yr/raw/` and `raw/ninv/`, 151 MB) into
  `~/data/wmap9_lowres_K/`; Ka, Q, V are in `~/data/wmap9_lowres/`.
- (a.2) E022's conditional check with the dense Ka, Q, V matrices reproduced within 5 percent.
- (b) An LFI 30 GHz dense low-resolution covariance: half-split null test with a PTE between 0.01
  and 0.99, diagonal within 20 percent of BeyondPlanck's simulated rms.
- (c) A beta = -3.1 sky pushed through the WMAP res 4 path and the LFI low-resolution path
  returns region indices within 0.01 with the noise off.

## Reference

I wrote no expectations down before a.1 was run. The only prior statement is the challenge
report's Gate 0 text: N^-1 symmetric to 1e-10 and positive semi-definite, with the processing
mask as its null space. The E031 white-noise estimate of the pixel rms (1.96 uK_RJ per nside-16
pixel for K, 1.99 for Ka) is the other number fixed in advance.

## Command

```bash
cd /Users/trygvels/code/Commander4
for b in "K 22.8 ~/data/wmap9_lowres_K" "Ka 33.0 ~/data/wmap9_lowres"; do set -- $b
  OMP_NUM_THREADS=1 .venv/bin/python -I nextgen/scripts/gate0_wmap_ninv_check.py --label $1 --nu $2 \
    --ninv $3/wmap_band_quninv_r4_9yr_$1_v5.fits --map $3/wmap_band_qumap_r4_9yr_$1_v5.fits \
    > nextgen/experiments/E038/wmap_$1_ninv_check.txt; done
```

## Result

**(a.1), K and Ka, 2026-11-09 23:46 to 23:48.**

| | K (22.8 GHz) | Ka (33.0 GHz) |
|---|---|---|
| Symmetry, max abs(N - N^T) / max abs(N) | 3.8e-10 | 1.0e-10 |
| Negative or null eigenvalues | none | none |
| Condition number | 425 | 298 |
| Pixel rms from the inverse, median (16 and 84 percent), uK_RJ | 2.01 (1.67, 2.42) | 2.04 (1.71, 2.45) |
| The same from the diagonal of N^-1 alone | 1.84 | 1.90 |
| Variance of the worst mode over the median mode | 92 | 65 |
| Second, fifth, tenth, fiftieth worst | 72, 25, 18, 4.9 | 39, 24, 16, 4.6 |
| Share of the total variance in the worst 10 and 100 modes | 0.047, 0.114 | 0.037, 0.101 |
| Mean multipole of the three worst modes (Q part) | 2.9, 6.2, 3.9 | 2.3, 3.9, 4.2 |
| Map chi-square for 6144 elements | 1.16e6 | 1.40e5 |

Against the Gate 0 text: symmetry holds to 1e-10 for Ka and to 4e-10 for K, which is the storage
precision. Both matrices are positive definite. There is no null space: the processing mask is
not projected out in these files, unlike what the Gate 0 text assumed. The pixel rms agrees with the E031 white-noise estimate
(1.96 and 1.99) to 3 percent. The noise is far from white on large scales: a few modes at
multipoles 2 to 6 carry 20 to 90 times the variance of the median mode. A white-noise Fisher
estimate of anything that lives at those multipoles, such as a latitude trend of the index or
the reionization signal, is optimistic by that order.

(a.2), (b), (c): not started.

## Negative control

None for a.1.

## Conclusion

Partial. The WMAP K and Ka inverse covariances are usable as released for a dense analysis. The
E031 white-noise numbers for the index per region hold at the pixel scale and are optimistic for
the largest scales.

## Follow-ups

- Q and V through the same check; a.2, b and c.
