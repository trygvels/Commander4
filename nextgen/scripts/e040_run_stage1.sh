#!/bin/sh
# E040 stage 1: WMAP K, Ka, Q, V and Planck 353 GHz. Run from the repository root, after E039 is read.
#   sh nextgen/scripts/e040_run_stage1.sh [threads]
# Every run writes <tag>.txt, <tag>.err and, where given, <tag>.npz under nextgen/experiments/E040/.
T=${1:-7}
export VECLIB_MAXIMUM_THREADS=$T OMP_NUM_THREADS=$T
S=nextgen/scripts/e040_wmap_region_indices.py
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { tag=$1; shift; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P $S "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
# baseline with every estimator, as pre-registered (no window matching) and window-matched (E038 c2)
run bpxv_n16            --regions bpxv --estimators profile flat zm
run bpxv_n16_wm         --regions bpxv --estimators profile flat zm --match-window 60
# profile-only arms, window-matched unless the arm is about the window
for W in "" "_wm"; do
  if [ -z "$W" ]; then M=""; else M="--match-window 60"; fi
  run bpxv_n16_noK$W      --regions bpxv --bands Ka Q V 353 $M
  run bpxv_n16_KKa$W      --regions bpxv --bands K Ka 353 $M
  run bpxv_n16_KQV$W      --regions bpxv --bands K Q V 353 $M
done
run bpxv_n16_cmbflat_wm   --regions bpxv --cmb-prior-sd 0 --match-window 60
run bpxv_n16_cmb01_wm     --regions bpxv --cmb-prior-sd 0.1 --match-window 60
run bpxv_n16_cmb04_wm     --regions bpxv --cmb-prior-sd 0.4 --match-window 60
run bpxv_n16_delta        --regions bpxv --delta-bands
run bpxv_n16_nodrift_wm   --regions bpxv --no-drift --match-window 60
run bpxv_n16_bd148_wm     --regions bpxv --dust-beta 1.48 --match-window 60
run bpxv_n16_bd158_wm     --regions bpxv --dust-beta 1.58 --match-window 60
run bpxv_n16_bp353delta_wm --regions bpxv --bp353 delta --match-window 60
run bpxv_n16_bp353full_wm --regions bpxv --bp353 full --match-window 60
run bpxv_n16_start27_wm   --regions bpxv --beta0 -2.7 --match-window 60
run bpxv_n16_start35_wm   --regions bpxv --beta0 -3.5 --match-window 60
run bpxv_n8_wm            --regions bpxv --nside 8 --match-window 60
run bpxv_n8               --regions bpxv --nside 8
run lat_n16_wm            --regions lat --match-window 60
