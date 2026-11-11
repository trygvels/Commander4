#!/bin/sh
# When the last low-pass arm of queue 10 has saved its results, stop queue 10 (shell first, then the job it has
# just started) and run queue 12 in its place. Waits on output files and on process IDs only.
F=nextgen/experiments/E042/lowpass30_lcdm.npz
until [ -s $F ] || ! kill -0 26584 2>/dev/null; do sleep 10; done
if kill -0 26584 2>/dev/null; then
  sleep 15
  kids=$(pgrep -P 26584)
  kill 26584
  sleep 1
  for k in $kids; do kill $k 2>/dev/null; done
  sleep 5
fi
rm -f nextgen/experiments/E040/cmbprior_sim_n16.txt nextgen/experiments/E040/cmbprior_sim_n16.err
sh nextgen/scripts/e043_queue12.sh > nextgen/experiments/E043/queue12.log 2>&1
