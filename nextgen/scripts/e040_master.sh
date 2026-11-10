#!/bin/sh
# E040 master queue, one dense job at a time, most informative first. Finished tags are skipped.
#   sh nextgen/scripts/e040_master.sh
S=nextgen/scripts/e040_wmap_region_indices.py
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { tag=$1; shift; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P $S "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
# 0. E042 (f): the Cosmoglobe DR1 reprocessed maps, white-noise weights
X=nextgen/experiments/E042
mkdir -p $X
[ -s $X/dr1_K30.txt ] && grep -q "^Plane" $X/dr1_K30.txt || $P nextgen/scripts/e040c_dr1_arm.py --bands K 30 353 --match-window 60 --out $X/dr1_K30.npz > $X/dr1_K30.txt 2> $X/dr1_K30.err
[ -s $X/dr1_KKa.txt ] && grep -q "^Plane" $X/dr1_KKa.txt || $P nextgen/scripts/e040c_dr1_arm.py --bands K Ka 353 --match-window 60 --out $X/dr1_KKa.npz > $X/dr1_KKa.txt 2> $X/dr1_KKa.err
[ -s $X/dr1_30Ka.txt ] && grep -q "^Plane" $X/dr1_30Ka.txt || $P nextgen/scripts/e040c_dr1_arm.py --bands 30 Ka 353 --match-window 60 --out $X/dr1_30Ka.npz > $X/dr1_30Ka.txt 2> $X/dr1_30Ka.err
[ -s $X/pr3_ncvm.txt ] && grep -q "30(PR4, same l) Ka 353" $X/pr3_ncvm.txt || $P nextgen/scripts/e042d_pr3_ncvm.py --lmax 24 --out $X/pr3_ncvm.npz > $X/pr3_ncvm.txt 2> $X/pr3_ncvm.err
# 1. the band pairs of the K band against 30 GHz comparison
run bpxv_n16_KKa_wm     --regions bpxv --bands K Ka 353 $W
run s2_K30_wm           --regions bpxv --bands K 30 353 $W
run s2_noK_lfi_wm       --regions bpxv --bands 30 44 70 Ka Q V 353 $W
run s2_all_wm           --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W
run s2_30Ka_wm          --regions bpxv --bands 30 Ka 353 $W
run bpxv_n16_noK_wm     --regions bpxv --bands Ka Q V 353 $W
# 2. sensitivity of the WMAP-only baseline
run bpxv_n16            --regions bpxv
run bpxv_n16_cmbflat_wm --regions bpxv --cmb-prior-sd 0 $W
run bpxv_n16_KQV_wm     --regions bpxv --bands K Q V 353 $W
run bpxv_n16_delta      --regions bpxv --delta-bands
run bpxv_n16_bd148_wm   --regions bpxv --dust-beta 1.48 $W
run bpxv_n16_bd158_wm   --regions bpxv --dust-beta 1.58 $W
run bpxv_n16_bp353delta_wm --regions bpxv --bp353 delta $W
run bpxv_n16_nodrift_wm --regions bpxv --no-drift $W
run bpxv_n16_cmb01_wm   --regions bpxv --cmb-prior-sd 0.1 $W
run bpxv_n16_cmb04_wm   --regions bpxv --cmb-prior-sd 0.4 $W
# 3. LFI noise arms and conventions for the K band against 30 GHz pair
run s2_K30_white_wm     --regions bpxv --bands K 30 353 $W --lfi-noise white
run s2_K30_A2_wm        --regions bpxv --bands K 30 353 $W --lfi-a-scale 2
run s2_K30_nodrift_wm   --regions bpxv --bands K 30 353 $W --no-drift
run s2_K30              --regions bpxv --bands K 30 353
run s2_all_white_wm     --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W --lfi-noise white
run s2_all_A2_wm        --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W --lfi-a-scale 2
# 4. the rest
run lat_n16_wm          --regions lat $W
run bpxv_n16_bp353full_wm --regions bpxv --bp353 full $W
run bpxv_n16_start27_wm --regions bpxv --beta0 -2.7 $W
run bpxv_n16_start35_wm --regions bpxv --beta0 -3.5 $W
run bpxv_n8             --regions bpxv --nside 8
run s2_K30_delta        --regions bpxv --bands K 30 353 --delta-bands
# 5. E039: plug-in addendum for the four regions, the rest of the 48-region simulation, its plug-in addendum
.venv/bin/python -I nextgen/scripts/e039b_plugin_check.py --nside 8 --regions bpxv --nreal 40 > nextgen/experiments/E039/plugin_n8_bpxv.txt 2> nextgen/experiments/E039/plugin_n8_bpxv.err
.venv/bin/python -I nextgen/scripts/e039_lowres_validation.py --nside 8 --regions hp2 --nreal 105 --first 95 --full 0 > nextgen/experiments/E039/sim_n8_hp2_part2.txt 2> nextgen/experiments/E039/sim_n8_hp2_part2.err
.venv/bin/python -I nextgen/scripts/e039b_plugin_check.py --nside 8 --regions hp2 --nreal 40 > nextgen/experiments/E039/plugin_n8_hp2.txt 2> nextgen/experiments/E039/plugin_n8_hp2.err
