#!/bin/sh
# Queue 10: what the paper needs with the LCDM CMB prior as the central configuration (E040 stage 4, C1 failed).
# Sound Cholesky path, malloc canary, one dense job at a time.
E=nextgen/experiments/E040
X2=nextgen/experiments/E042
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I"
W="--match-window 60"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040f_block_errors_v2.py maps --cmb-prior-lcdm --out $E/block_v2_maps_lcdm.npz > $E/block_v2_maps_lcdm.txt 2> $E/block_v2_maps_lcdm.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040f_block_errors_v2.py maps --out $E/block_v2_maps_white.npz > $E/block_v2_maps_white.txt 2> $E/block_v2_maps_white.err
run e040b_estimators.py l_bpxv_n16_nowm        --regions bpxv --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_KKa_wm      --regions bpxv --bands K Ka 353 $W --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_nodrift_wm  --regions bpxv --no-drift $W --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_bd148_wm    --regions bpxv --dust-beta 1.48 $W --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_bd158_wm    --regions bpxv --dust-beta 1.58 $W --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_bp353full_wm --regions bpxv --bp353 full $W --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_start27_wm  --regions bpxv --beta0 -2.7 $W --cmb-prior-lcdm
run e040b_estimators.py l_bpxv_n16_delta_wm    --regions bpxv --delta-bands $W --cmb-prior-lcdm
for L in 20 12 30; do
  [ -s $X2/lowpass${L}_lcdm.txt ] && grep -q "^K,Ka,353" $X2/lowpass${L}_lcdm.txt || DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e042e_lowpass.py --lmax $L --sets K,Ka,Q,V,353 K,Ka,353 --cmb-prior-lcdm --out $X2/lowpass${L}_lcdm.npz > $X2/lowpass${L}_lcdm.txt 2> $X2/lowpass${L}_lcdm.err
done
# simulation spot checks at the working resolution (E040 C4 and E039 addendum B), before the no-K band set
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039d_cmbprior_sim.py --nside 16 --nreal 6 --priors lcdm white > $E/cmbprior_sim_n16.txt 2> $E/cmbprior_sim_n16.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039e_lowpass_sim_n16.py --lmax 20 --nreal 6 > nextgen/experiments/E039/lowpass_n16_l20.txt 2> nextgen/experiments/E039/lowpass_n16_l20.err
run e040b_estimators.py s3_bpset_lcdm_wm       --regions bpxv --bands 30 44 70 Ka Q V 353 $W --cmb-prior-lcdm --estimators profile flat
run e040b_estimators.py s3_bpset_cmb02_wm      --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
