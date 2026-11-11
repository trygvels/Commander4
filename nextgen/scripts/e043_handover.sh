#!/bin/sh
# When the flat-CMB-prior fit of queue 9 has written its table, stop queue 9 (shell first, then the job it has
# just started) and run queue 10 in its place. Waits on the output file and on process IDs only.
F=nextgen/experiments/E040/s4_bpxv_n16_wm_cmbflat.txt
until grep -q "^Plane" $F 2>/dev/null || ! kill -0 21597 2>/dev/null; do sleep 10; done
if kill -0 21597 2>/dev/null; then
  kids=$(pgrep -P 21597)
  kill 21597
  sleep 1
  for k in $kids; do kill $k 2>/dev/null; done
  sleep 5
fi
rm -f nextgen/experiments/E040/bpxv_n16_nodrift_wm.txt nextgen/experiments/E040/bpxv_n16_nodrift_wm.err
sh nextgen/scripts/e043_queue10.sh > nextgen/experiments/E043/queue10.log 2>&1
