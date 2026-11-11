#!/bin/sh
# E044: simulations at nside 16 with the regions, mask, covariances and prior of the map fit (staged; one heavy job at a time).
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e040j_sim_n16_stages.py"
S=$HOME/data/c4_tools/sim_n16_named
DYLD_INSERT_LIBRARIES=$C $P draw --nreal 6 --lmax 20 --geometry named --dir $S > $E/sim_n16_draw.txt 2> $E/sim_n16_draw.err
DYLD_INSERT_LIBRARIES=$C $P blocknull --nreal 6 --dir $S --out $E/sim_n16_blocknull.npz > $E/sim_n16_blocknull.txt 2> $E/sim_n16_blocknull.err
DYLD_INSERT_LIBRARIES=$C $P lowpass --nreal 4 --dir $S --out $E/sim_n16_lowpass.npz > $E/sim_n16_lowpass.txt 2> $E/sim_n16_lowpass.err
