#!/bin/sh
# E040 priority queue: the K band against 30 GHz pairs, then stage 3 (the BeyondPlanck band set and reference),
# then the master queue for everything else. One dense job at a time.
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
run e040_wmap_region_indices.py bpxv_n16_KKa_wm --regions bpxv --bands K Ka 353 $W
run e040_wmap_region_indices.py s2_K30_wm       --regions bpxv --bands K 30 353 $W
run e040_wmap_region_indices.py s2_noK_lfi_wm   --regions bpxv --bands 30 44 70 Ka Q V 353 $W
run e040b_estimators.py s3_bpset_cmbflat_wm     --regions bpxv --bands 30 44 70 Ka Q V 353 $W --cmb-prior-sd 0 --estimators profile flat
run e040b_estimators.py s3_bpset_cmb02_wm       --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
run e040_wmap_region_indices.py s2_all_wm       --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W
run e040_wmap_region_indices.py s2_30Ka_wm      --regions bpxv --bands 30 Ka 353 $W
sh nextgen/scripts/e040_master.sh > $E/master.log 2>&1
