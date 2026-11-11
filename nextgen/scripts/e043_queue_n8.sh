#!/bin/sh
# nside-8 chain after the E043 fix: every run on the sound Cholesky path and under the malloc canary.
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I"
X=nextgen/experiments/E043
E9=nextgen/experiments/E039
E0=nextgen/experiments/E040
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039b_plugin_check.py --nside 8 --regions bpxv --nreal 40 --no-drift > $X/plugin_n8_bpxv_rerun.txt 2> $X/plugin_n8_bpxv_rerun.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040e_block_errors_sim.py --nreal 60 > $E0/block_errors_sim_spread0.txt 2> $E0/block_errors_sim_spread0.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040e_block_errors_sim.py --nreal 60 --spread 0.1 > $E0/block_errors_sim_spread01.txt 2> $E0/block_errors_sim_spread01.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039c_lowpass_sim.py --lmax 12 --nreal 100 > $E9/lowpass_n8_l12.txt 2> $E9/lowpass_n8_l12.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039_lowres_validation.py --nside 8 --regions hp2 --nreal 95 --first 0 --full 0 --no-drift > $X/sim_n8_hp2_part1_rerun.txt 2> $X/sim_n8_hp2_part1_rerun.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039_lowres_validation.py --nside 8 --regions hp2 --nreal 105 --first 95 --full 0 --no-drift > $E9/sim_n8_hp2_part2.txt 2> $E9/sim_n8_hp2_part2.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e039b_plugin_check.py --nside 8 --regions hp2 --nreal 40 --no-drift > $E9/plugin_n8_hp2.txt 2> $E9/plugin_n8_hp2.err
