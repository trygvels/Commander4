#!/bin/sh
# E044: fits on the named regions with the mask. Sound Cholesky path, malloc canary, one heavy job at a time.
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e044_fit.py"
run() { tag=$1; shift; [ -s $E/$tag.txt ] && grep -q "^Plane" $E/$tag.txt && return 0; DYLD_INSERT_LIBRARIES=$C $P "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
run lcdm_all    --prior lcdm --estimators profile flat zm --block
run bpxv_tauA   --geometry bpxv --tau-a-only --prior lcdm --block
run white       --prior white --mismatch --estimators profile flat
run flat        --prior flat --estimators profile flat
run lp20        --prior lcdm --lowpass 20
run bd148       --prior lcdm --dust-beta 1.48
run bd158       --prior lcdm --dust-beta 1.58
run delta       --prior lcdm --delta-bands
run KKa         --prior lcdm --bands K Ka 353
run lp12        --prior lcdm --lowpass 12
run lp30        --prior lcdm --lowpass 30
run nodrift     --prior lcdm --no-drift
run bp353full   --prior lcdm --bp353 full
run lp20_KKa    --prior lcdm --lowpass 20 --bands K Ka 353
run floor011    --prior lcdm --floor-uk 0.11
