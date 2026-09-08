#!/bin/bash
# Usage: gen_hd_combo.sh <P> <SIMS_PER_CELL> <CHUNK> <OUT> [K]
P=$1; NSIMS=$2; CHUNK=$3; OUT=$4; K=${5:-2}; START=${6:-1}
: > "$OUT"
for cfg in C1 C2 C3; do
  for rho in 0 0.5 1.0 1.5 2.0 2.5; do
    nd=$(awk -v r=$rho 'BEGIN{print (r>0)?1:0}')
    s=$START; while [ $s -le $NSIMS ]; do e=$((s+CHUNK-1)); [ $e -gt $NSIMS ]&&e=$NSIMS
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$K" "$cfg" "$P" "$rho" "$nd" "$s" "$e" >> "$OUT"
      s=$((e+1)); done
  done
done
echo "wrote $OUT: $(wc -l < $OUT) chunk-tasks (K=$K, P=$P, $NSIMS sims/cell, $CHUNK/chunk)"
