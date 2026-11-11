#!/bin/sh
# E044 H9: the zero-mean prior against the profile optimum in simulations at nside 16 (regions, mask, covariances
# and prior of the map fit). Reads the realizations drawn by e044_queue_sim16.sh. One heavy job at a time.
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e040j_sim_n16_stages.py"
S=$HOME/data/c4_tools/sim_n16_named
DYLD_INSERT_LIBRARIES=$C $P zm --nreal 2 --dir $S --out $E/sim_n16_zm.npz > $E/sim_n16_zm.txt 2> $E/sim_n16_zm.err
