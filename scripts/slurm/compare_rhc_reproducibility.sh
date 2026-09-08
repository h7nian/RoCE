#!/bin/bash
#SBATCH --job-name=roce_rhc_repro
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rhc_repro.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rhc_repro.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
REFERENCE_DIRECTORY="${1:?reference RHC directory is required}"
CANDIDATE_DIRECTORY="${2:?candidate RHC directory is required}"
OUTPUT_DIRECTORY="${3:?audit output directory is required}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/compare_rhc_reproducibility.R \
  "${REFERENCE_DIRECTORY}" "${CANDIDATE_DIRECTORY}" "${OUTPUT_DIRECTORY}"
