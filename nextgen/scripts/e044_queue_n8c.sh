#!/bin/sh
# E044 H10: a third, independent null set at nside 8 (240 realizations, seed 43) for the calibration of the sandwich
# error and the null of the block-error ratio. A light job (2 to 4.5 GB) that may run beside one nside-16 fit.
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e040j_sim_n16_stages.py"
S=$HOME/data/c4_tools/sim_n8_named_c
DYLD_INSERT_LIBRARIES=$C $P draw --nreal 240 --nside 8 --lmax 12 --seed 43 --geometry named --dir $S > $E/sim_n8_draw_c.txt 2> $E/sim_n8_draw_c.err
DYLD_INSERT_LIBRARIES=$C $P blocknull --nreal 240 --nside 8 --dir $S --out $E/sim_n8_blocknull_c.npz > $E/sim_n8_blocknull_c.txt 2> $E/sim_n8_blocknull_c.err
