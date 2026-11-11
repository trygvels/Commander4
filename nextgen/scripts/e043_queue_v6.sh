#!/bin/sh
# E043 V6: the other nside-8 runs (E038 window checks at nside 8, the E040 nside-8 arm) on the sound path under the canary.
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I"
X=nextgen/experiments/E043
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e038_window_check.py --common-beam 60 --nside-fit 8 > $X/window_control_n8_rerun.txt 2> $X/window_control_n8_rerun.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e038_window_check.py --match-window 60 --nside-fit 8 > $X/window_wmap_matched_n8_rerun.txt 2> $X/window_wmap_matched_n8_rerun.err
DYLD_INSERT_LIBRARIES=$C $P nextgen/scripts/e040_wmap_region_indices.py --regions bpxv --nside 8 --match-window 60 --estimators profile flat zm > $X/bpxv_n8_wm_all_rerun.txt 2> $X/bpxv_n8_wm_all_rerun.err
