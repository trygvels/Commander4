#!/bin/sh
# E044 R4: chi-square of the profile fit under the model at nside 16 (total, per region and per band), on the six
# realizations drawn by e044_queue_sim16.sh. Runs behind e044_queue_after.sh. One heavy job at a time.
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
DYLD_INSERT_LIBRARIES=$C .venv/bin/python -I nextgen/scripts/e044_chi2_null.py --dir $HOME/data/c4_tools/sim_n16_named --nreal 6 --out $E/sim_n16_chi2null.npz > $E/sim_n16_chi2null.txt 2> $E/sim_n16_chi2null.err
