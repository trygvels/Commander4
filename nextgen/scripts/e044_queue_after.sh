#!/bin/sh
# E044, behind the fits of e044_queue.sh: the scans of the flat-prior marginal (H11), then the nside-16 simulations
# (e044_queue_sim16.sh) and the zero-mean simulation at nside 16 (e044_queue_sim16b.sh). One heavy job at a time.
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e044_flat_scan.py"
DYLD_INSERT_LIBRARIES=$C $P --prior flat --fit $E/flat.npz --out $E/scan_flat.npz > $E/scan_flat.txt 2> $E/scan_flat.err
DYLD_INSERT_LIBRARIES=$C $P --prior lcdm --fit $E/lcdm_all.npz --out $E/scan_lcdm.npz > $E/scan_lcdm.txt 2> $E/scan_lcdm.err
sh nextgen/scripts/e044_queue_sim16.sh > $E/queue_sim16.log 2>&1
sh nextgen/scripts/e044_queue_sim16b.sh > $E/queue_sim16b.log 2>&1
