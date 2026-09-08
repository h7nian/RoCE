#!/bin/bash
#SBATCH --job-name=roce_parallel_qc
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_parallel_qc.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_parallel_qc.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
if [[ "$#" -ne 3 ]]; then
  echo "usage: $0 REFERENCE.csv CANDIDATE.csv OUTPUT_DIR" >&2
  exit 1
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/compare_parallel_reproducibility.R "$@"
