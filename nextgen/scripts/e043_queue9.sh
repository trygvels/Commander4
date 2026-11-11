#!/bin/sh
# Queue 9 (after the E043 fix): every run on the sound Cholesky path and under the malloc canary. One dense job at a time.
# V1 and V2 of E043 first (the nside-16 baseline on the old path and on the new one), then what queues 7 and 8 still owed.
E=nextgen/experiments/E040
X=nextgen/experiments/E043
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
LOWRES_CHOLESKY_UPLO_U=1 DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040_wmap_region_indices.py --regions bpxv $W > $X/bpxv_n16_wm_oldpath.txt 2> $X/bpxv_n16_wm_oldpath.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040_wmap_region_indices.py --regions bpxv $W > $X/bpxv_n16_wm_newpath.txt 2> $X/bpxv_n16_wm_newpath.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040d_block_errors.py --out $E/block_errors_bpxv_n16_wm.npz > $E/block_errors_bpxv_n16_wm.txt 2> $E/block_errors_bpxv_n16_wm.err
run e040b_estimators.py s4_bpxv_n16_wm_cmbflat --regions bpxv $W --cmb-prior-sd 0 --estimators profile flat zm
run e040_wmap_region_indices.py bpxv_n16_nodrift_wm   --regions bpxv --no-drift $W
run e040_wmap_region_indices.py bpxv_n16_delta        --regions bpxv --delta-bands
run e040_wmap_region_indices.py bpxv_n16_bd148_wm     --regions bpxv --dust-beta 1.48 $W
run e040_wmap_region_indices.py bpxv_n16_bd158_wm     --regions bpxv --dust-beta 1.58 $W
run e040_wmap_region_indices.py bpxv_n16_bp353full_wm --regions bpxv --bp353 full $W
run e040_wmap_region_indices.py bpxv_n16_start27_wm   --regions bpxv --beta0 -2.7 $W
run e040b_estimators.py s3_bpset_cmb02_wm      --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
