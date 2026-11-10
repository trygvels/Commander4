#!/bin/sh
# Coverage runs, one thread (a second thread is left for the reviewer's foreground checks).
cd /Users/trygvels/code/Commander4
export VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
D=nextgen/experiments/E039/stats_review_scripts
for cfg in "S1 2000" "S2 2000" "F1 800" "F2 800" "F3 800" "F4 800" "S4 2000" "S3 500"; do
  set -- $cfg
  .venv/bin/python $D/check_03_coverage.py $1 $2 2>&1 | grep -v -e RuntimeWarning -e "sign, logdet" -e "^  *[a-zA-Z_]* = " > $D/out_03_coverage_$1.txt
done
