#!/bin/sh
# Queue 11: K band with PR4 30 GHz under the LCDM CMB prior, with the two noise-model arms the referee asked for
# (white-noise columns only; the power-law term of the noise model doubled). Sound path, canary, one job at a time.
E=nextgen/experiments/E040
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I"
W="--match-window 60"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
run e040b_estimators.py l_K30_wm          --regions bpxv --bands K 30 353 $W --cmb-prior-lcdm
run e040b_estimators.py l_K30_wm_white    --regions bpxv --bands K 30 353 $W --cmb-prior-lcdm --lfi-noise white
run e040b_estimators.py l_K30_wm_a2       --regions bpxv --bands K 30 353 $W --cmb-prior-lcdm --lfi-a-scale 2.0
