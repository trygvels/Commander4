#!/bin/sh
# Queue 13: the nside-16 simulation at the configuration of the map fit, in stages that fit the memory
# (e040j_sim_n16_stages.py replaces e040h_sim_n16.py, which needed 21 GB). Sound path, canary, one job at a time.
E=nextgen/experiments/E040
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e040j_sim_n16_stages.py"
S=$HOME/data/c4_tools/sim_n16
DYLD_INSERT_LIBRARIES=$C $P draw --nreal 12 --lmax 20 --dir $S > $E/sim_n16_draw.txt 2> $E/sim_n16_draw.err
DYLD_INSERT_LIBRARIES=$C $P blocknull --nreal 12 --dir $S --out $E/sim_n16_blocknull.npz > $E/sim_n16_blocknull.txt 2> $E/sim_n16_blocknull.err
DYLD_INSERT_LIBRARIES=$C $P white --nreal 6 --dir $S --out $E/sim_n16_white.npz > $E/sim_n16_white.txt 2> $E/sim_n16_white.err
DYLD_INSERT_LIBRARIES=$C $P lowpass --nreal 6 --dir $S --out $E/sim_n16_lowpass.npz > $E/sim_n16_lowpass.txt 2> $E/sim_n16_lowpass.err
