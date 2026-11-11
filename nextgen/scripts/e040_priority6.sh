#!/bin/sh
# Priority queue, sixth version: the reworked low-pass arm first. One dense job at a time.
E=nextgen/experiments/E040
X=nextgen/experiments/E042
P=".venv/bin/python -I"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
lp() { L=$1; [ -s $X/lowpass$L.txt ] && grep -q "^30,Ka,353" $X/lowpass$L.txt || $P nextgen/scripts/e042e_lowpass.py --lmax $L --out $X/lowpass$L.npz > $X/lowpass$L.txt 2> $X/lowpass$L.err; }
W="--match-window 60"
lp 20
[ -s $X/pr3_arm.txt ] && grep -q "30(PR3) Ka 353" $X/pr3_arm.txt || $P nextgen/scripts/e042c_pr3_arm.py --out $X/pr3_arm.npz > $X/pr3_arm.txt 2> $X/pr3_arm.err
lp 12
lp 30
[ -s $X/nuisance.txt ] && grep -q "angle fixed" $X/nuisance.txt || $P nextgen/scripts/e042_k30_anatomy.py --cases amp_global amp_region angle --out $X/nuisance.npz > $X/nuisance.txt 2> $X/nuisance.err
run e040b_estimators.py s3_bpset_cmb02_wm       --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
run e040_wmap_region_indices.py s2_all_wm       --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W
run e040_wmap_region_indices.py s2_30Ka_wm      --regions bpxv --bands 30 Ka 353 $W
sh nextgen/scripts/e040_master.sh > $E/master.log 2>&1
