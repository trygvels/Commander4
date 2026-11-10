#!/bin/sh
# Priority queue, second version: stage 3 (BeyondPlanck band set), the E042 anatomy of K band against 30 GHz,
# the remaining band sets, then the master queue for everything else. One dense job at a time.
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
run e040_wmap_region_indices.py s2_noK_lfi_wm   --regions bpxv --bands 30 44 70 Ka Q V 353 $W
run e040b_estimators.py s3_bpset_cmbflat_wm     --regions bpxv --bands 30 44 70 Ka Q V 353 $W --cmb-prior-sd 0 --estimators profile flat
mkdir -p nextgen/experiments/E042
X=nextgen/experiments/E042
[ -s $X/base_cuts.txt ] && grep -q "K cut l<=20" $X/base_cuts.txt || $P nextgen/scripts/e042_k30_anatomy.py --cases base cut30 cutK --out $X/base_cuts.npz > $X/base_cuts.txt 2> $X/base_cuts.err
[ -s $X/modes_ctrl.txt ] && grep -q "worst K modes out" $X/modes_ctrl.txt || $P nextgen/scripts/e042_k30_anatomy.py --cases worstK ctrlKa --out $X/modes_ctrl.npz > $X/modes_ctrl.txt 2> $X/modes_ctrl.err
[ -s $X/nuisance.txt ] && grep -q "angle fixed" $X/nuisance.txt || $P nextgen/scripts/e042_k30_anatomy.py --cases amp_global amp_region angle --out $X/nuisance.npz > $X/nuisance.txt 2> $X/nuisance.err
run e040b_estimators.py s3_bpset_cmb02_wm       --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
run e040_wmap_region_indices.py s2_all_wm       --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W
run e040_wmap_region_indices.py s2_30Ka_wm      --regions bpxv --bands 30 Ka 353 $W
sh nextgen/scripts/e040_master.sh > $E/master.log 2>&1
