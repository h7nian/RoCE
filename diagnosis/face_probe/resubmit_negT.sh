#!/bin/bash
# Resubmit the negative-transfer chunks that were lost to the 45-minute wall-clock
# limit (pure timeouts -- the chunk logs show "DUE TO TIME LIMIT", with zero
# numerical/sim errors). Runs ONE sim per array task with a generous wall-time so
# the tasks cannot time out again, reusing the existing face_chunk.cmd runner and
# the combo produced by gen_missing_combo.sh. After the array finishes, re-run
# face_negT_aggregate.R to regenerate negtransfer_*.csv at the full 200 repeats.
#
# Usage: bash diagnosis/face_probe/resubmit_negT.sh
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"   # repo root = SLURM submit dir

COMBO=diagnosis/face_probe/validation/chunks/combo_resubmit_negT.tsv
[ -s "$COMBO" ] || { echo "missing/empty combo: $COMBO (run gen_missing_combo.sh first)"; exit 1; }
N=$(wc -l < "$COMBO")
echo "[resubmit_negT] submitting $N one-sim tasks from $COMBO"

sbatch --job-name=negT_resub \
       --time=08:00:00 \
       --cpus-per-task=2 \
       --array=1-"${N}"%50 \
       --export=ALL,COMBO_FILE="${COMBO}" \
       diagnosis/face_probe/face_chunk.cmd
