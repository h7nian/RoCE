#!/bin/bash
#SBATCH --job-name=roce_cv_tail
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=5
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_cv_tail.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_cv_tail.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB is required}"
OUTPUT_PATH="${ROCE_CV_TAIL_OUTPUT:?ROCE_CV_TAIL_OUTPUT is required}"
CONFIGURATION="${ROCE_CV_TAIL_CONFIG:-C3}"
SOURCE_COUNT="${ROCE_CV_TAIL_K:-4}"
SIMULATION_ID="${ROCE_CV_TAIL_SIM_ID:-1}"
TREATMENT_VALUE="${ROCE_CV_TAIL_A_VAL:-1}"
SOURCE_INDEX="${ROCE_CV_TAIL_SOURCE_INDEX:-1}"
OUTER_FOLD="${ROCE_CV_TAIL_OUTER_FOLD:-1}"
INNER_FOLD="${ROCE_CV_TAIL_INNER_FOLD:-2}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_NUISANCE_CV_THREADS=5
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
Rscript scripts/slurm/probe_nuisance_cv_tail.R \
  "${OUTPUT_PATH}" "${CONFIGURATION}" "${SOURCE_COUNT}" \
  "${SIMULATION_ID}" "${TREATMENT_VALUE}" "${SOURCE_INDEX}" \
  "${OUTER_FOLD}" "${INNER_FOLD}"
