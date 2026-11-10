#!/bin/sh
# E040 stage 2: Planck LFI and the K band against 30 GHz comparison. Run from the repository root after
# stage 1 has finished (one dense job at a time on this machine).
#   sh nextgen/scripts/e040_run_stage2.sh
S=nextgen/scripts/e040_wmap_region_indices.py
E=nextgen/experiments/E040
P=".venv/bin/python -I"
run() { tag=$1; shift; [ -s $E/$tag.txt ] && grep -q "^High Latitude\|^0 \|^region" $E/$tag.txt && return 0; $P $S "$@" --out $E/$tag.npz > $E/$tag.txt 2> $E/$tag.err; }
W="--match-window 60"
# the three pairs of the pre-registration and the full set
run s2_K30_wm          --regions bpxv --bands K 30 353 $W
run s2_noK_lfi_wm      --regions bpxv --bands 30 44 70 Ka Q V 353 $W
run s2_all_wm          --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W
run s2_30Ka_wm         --regions bpxv --bands 30 Ka 353 $W
# arms that go with every LFI result
run s2_K30_white_wm    --regions bpxv --bands K 30 353 $W --lfi-noise white
run s2_K30_whiter_wm   --regions bpxv --bands K 30 353 $W --lfi-noise white --lfi-noise-scale 1.21
run s2_K30_A2_wm       --regions bpxv --bands K 30 353 $W --lfi-a-scale 2
run s2_all_white_wm    --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W --lfi-noise white
run s2_all_A2_wm       --regions bpxv --bands K 30 Ka Q 44 V 70 353 $W --lfi-a-scale 2
run s2_K30_nodrift_wm  --regions bpxv --bands K 30 353 $W --no-drift
run s2_K30             --regions bpxv --bands K 30 353
run s2_K30_delta       --regions bpxv --bands K 30 353 --delta-bands
