#!/bin/sh
# Queue after the referee report on draft 1: what the reduced Paper A1 needs. One dense job at a time.
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
[ -s nextgen/experiments/E039/plugin_n8_bpxv.txt ] && grep -q "all measured regions" nextgen/experiments/E039/plugin_n8_bpxv.txt || $P nextgen/scripts/e039b_plugin_check.py --nside 8 --regions bpxv --nreal 40 > nextgen/experiments/E039/plugin_n8_bpxv.txt 2> nextgen/experiments/E039/plugin_n8_bpxv.err
run e040b_estimators.py s4_bpxv_n16_wm_lcdm    --regions bpxv $W --cmb-prior-lcdm --estimators profile flat zm
run e040b_estimators.py s4_bpxv_n16_wm_cmbflat --regions bpxv $W --cmb-prior-sd 0 --estimators profile flat zm
run e040b_estimators.py s3_bpset_cmb02_wm      --regions bpxv --bands 30 44 70 Ka Q V 353 $W --estimators profile flat
