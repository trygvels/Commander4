#!/bin/sh
# Queue 12 (replaces the tail of queue 10 and queue 11 after the independent review of 2026-11-11):
# block errors on the maps with diagnostics, the LCDM prior with the power above the band limit as its floor,
# K band with 30 GHz and its noise arms, the nside-16 simulation at the configuration of the map fit, the no-K band set.
E=nextgen/experiments/E040
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I"
W="--match-window 60"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040i_block_maps.py --prior lcdm --out $E/block_v3_maps_lcdm.npz > $E/block_v3_maps_lcdm.txt 2> $E/block_v3_maps_lcdm.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040i_block_maps.py --prior lcdm --floor-uk 0.11 --out $E/block_v3_maps_lcdm_floor011.npz > $E/block_v3_maps_lcdm_floor011.txt 2> $E/block_v3_maps_lcdm_floor011.err
run e040b_estimators.py l_K30_wm          --regions bpxv --bands K 30 353 $W --cmb-prior-lcdm
run e040b_estimators.py l_K30_wm_white    --regions bpxv --bands K 30 353 $W --cmb-prior-lcdm --lfi-noise white
run e040b_estimators.py l_K30_wm_a2       --regions bpxv --bands K 30 353 $W --cmb-prior-lcdm --lfi-a-scale 2.0
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040h_sim_n16.py --nreal 12 --full 6 --lmax 20 --out $E/sim_n16_lcdm.npz > $E/sim_n16_lcdm.txt 2> $E/sim_n16_lcdm.err
run e040b_estimators.py s3_bpset_lcdm_wm       --regions bpxv --bands 30 44 70 Ka Q V 353 $W --cmb-prior-lcdm --estimators profile flat
run e040b_estimators.py s3_bpset_cmb02_wm      --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
