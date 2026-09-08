#!/bin/bash
# Scan the existing negative-transfer chunk CSVs and emit a CHUNK=1 combo of the
# sim-ids that are MISSING from the target repeat count, for resubmission via
# face_chunk.cmd. Missing sims are caused by array tasks that hit the wall-clock
# limit and were killed before writing (no numerical failures); see the chunk
# logs (DUE TO TIME LIMIT). Same TSV columns as gen_hd_combo.sh:
#   K  CFG  P  RHO  NDEV  S0  S1   (here S0==S1, one sim per task).
# Usage: gen_missing_combo.sh <OUT.tsv>
# Scope: p in {50,100}, target 200 repeats/cell (the "didn't reach 200" cells).
set -euo pipefail
OUT=${1:?usage: gen_missing_combo.sh <OUT.tsv>}
CDIR=diagnosis/face_probe/validation/chunks
TARGET=200
: > "$OUT"
for P in 50 100; do
  for K in 2 4 8; do
    for CFG in C1 C2 C3; do
      for RHO in 0 0.5 1 1.5 2 2.5; do
        nd=$(awk -v r="$RHO" 'BEGIN{print (r>0)?1:0}')
        present=$(awk -F, 'FNR>1{gsub(/"/,"",$1);print $1}' \
                    "$CDIR"/negT_K${K}_${CFG}_p${P}_rho${RHO}_*.csv 2>/dev/null \
                  | sort -u)
        comm -23 <(seq 1 "$TARGET" | sort) <(printf '%s\n' "$present" | grep -v '^$' | sort) \
        | while read -r s; do
            [ -n "$s" ] && printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$K" "$CFG" "$P" "$RHO" "$nd" "$s" "$s" >> "$OUT"
          done
      done
    done
  done
done
echo "wrote $OUT: $(wc -l < "$OUT") missing-sim tasks (CHUNK=1, target $TARGET)"
