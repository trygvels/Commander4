#!/bin/zsh
# Resumable fetch of the PR4 (R4.00) full-mission HFI maps 100 to 353 GHz from the NERSC portal.
BASE=https://portal.nersc.gov/cfs/cmb/planck2020/frequency_maps/Single-frequency
OUT=/Users/trygvels/data/planck_pr4
EXPECT=2013275520
for f in 100 143 217 353; do
  name=HFI_SkyMap_${f}_2048_R4.00_full.fits
  for attempt in $(seq 1 60); do
    size=$(stat -f %z $OUT/$name 2>/dev/null || echo 0)
    if [ "$size" = "$EXPECT" ]; then echo "$name complete ($(date +%H:%M:%S))"; break; fi
    echo "$name attempt $attempt from $size bytes ($(date +%H:%M:%S))"
    curl -s -S -C - --retry 0 -m 3600 -o $OUT/$name $BASE/$name
    sleep 5
  done
done
echo ALLDONE
