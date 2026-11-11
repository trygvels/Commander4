#!/bin/sh
# E044 after the third referee read (H12, H13, H15, H16), behind e044_queue_after.sh. One heavy job at a time.
#   the present regions without the mask (all scales with block errors; l <= 20); mask thresholds of 3 and 10 percent;
#   the brightest 2 percent in free-free masked as well (LCDM and flat CMB prior); the chi-square null at nside 16
#   (e044_queue_after2.sh); forty realizations at nside 16 with block errors.
E=nextgen/experiments/E044
C=$HOME/data/c4_tools/libcanary.dylib
P=".venv/bin/python -I nextgen/scripts/e044_fit.py"
run() { tag=$1; shift; [ -s $E/$tag.txt ] && grep -q "^Plane" $E/$tag.txt && return 0; DYLD_INSERT_LIBRARIES=$C $P "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
run ffmask_lcdm --prior lcdm --mask-extra $E/mask_ff02.npy --block
run ffmask_flat --prior flat --mask-extra $E/mask_ff02.npy
run nomask      --prior lcdm --no-mask --block
run nomask_lp20 --prior lcdm --no-mask --lowpass 20
run bright03    --prior lcdm --bright-fraction 0.03
run bright10    --prior lcdm --bright-fraction 0.10
sh nextgen/scripts/e044_queue_after2.sh > $E/queue_after2.log 2>&1
Q=".venv/bin/python -I nextgen/scripts/e040j_sim_n16_stages.py"
S=$HOME/data/c4_tools/sim_n16_named_b
DYLD_INSERT_LIBRARIES=$C $Q draw --nreal 40 --lmax 20 --seed 45 --geometry named --dir $S > $E/sim_n16_draw_b.txt 2> $E/sim_n16_draw_b.err
DYLD_INSERT_LIBRARIES=$C $Q blocknull --nreal 40 --dir $S --out $E/sim_n16_blocknull_b.npz > $E/sim_n16_blocknull_b.txt 2> $E/sim_n16_blocknull_b.err
