#!/bin/sh
# Queue 8: the sensitivity arms of the WMAP-only fit that the reduced Paper A1 quotes, then the rest of E039.
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { script=$1; tag=$2; shift 2; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P nextgen/scripts/$script "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
run e040_wmap_region_indices.py bpxv_n16_nodrift_wm   --regions bpxv --no-drift $W
run e040_wmap_region_indices.py bpxv_n16_delta        --regions bpxv --delta-bands
run e040_wmap_region_indices.py bpxv_n16_bd148_wm     --regions bpxv --dust-beta 1.48 $W
run e040_wmap_region_indices.py bpxv_n16_bd158_wm     --regions bpxv --dust-beta 1.58 $W
run e040_wmap_region_indices.py bpxv_n16_bp353full_wm --regions bpxv --bp353 full $W
run e040_wmap_region_indices.py bpxv_n16_start27_wm   --regions bpxv --beta0 -2.7 $W
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions hp2 --nreal 105 --first 95 --full 0 > nextgen/experiments/E039/sim_n8_hp2_part2.txt 2> nextgen/experiments/E039/sim_n8_hp2_part2.err
.venv/bin/python -I nextgen/scripts/e039b_plugin_check.py --nside 8 --regions hp2 --nreal 40 > nextgen/experiments/E039/plugin_n8_hp2.txt 2> nextgen/experiments/E039/plugin_n8_hp2.err
