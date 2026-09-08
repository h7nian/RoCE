#!/bin/bash
#SBATCH --job-name=roce_legacy_qc
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_v1/logs/%j_legacy_qc.out
#SBATCH --error=results/direct_tate_v1/logs/%j_legacy_qc.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
INPUT_ROOT="${PROJECT_ROOT}/diagnosis/face_probe/validation"
OUTPUT_ROOT="${PROJECT_ROOT}/results/direct_tate_v1/legacy_p100_diagnostics"

module load R/4.2.2-gcc-8.2.0-vp7tyde

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" "${PROJECT_ROOT}/results/direct_tate_v1/logs"

Rscript scripts/slurm/diagnose_legacy_p100_summaries.R \
  "${INPUT_ROOT}" "${OUTPUT_ROOT}" 200
