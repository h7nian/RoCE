#!/bin/bash
# Submit the rho=0 n_folds coverage pilot as one Slurm array (one seed per task).
#
# usage: submit_pilot_nfolds.sh CONFIG K N_FOLDS SEED_START SEED_END [MAX_CONCURRENT]
# env:   ROCE_PILOT_LIB (default results/direct_tate_mc500_b5000/Rlib_pilot_nfolds_20260917)
#        ROCE_PILOT_ROOT (default results/direct_tate_mc500_b5000/pilot_nfolds_20260917)
#        ROCE_NUISANCE_CV_THREADS (default 5), ROCE_SLURM_PARTITION (default msismall)
#        ROCE_NUISANCE_SOLVER is passed through unchanged (set coordinate_descent
#        for a paired audit against the original solver)
set -euo pipefail

CONFIG="$1"; K="$2"; N_FOLDS="$3"; SEED_START="$4"; SEED_END="$5"
MAX_CONCURRENT="${6:-50}"
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
cd "${PROJECT_ROOT}"
PILOT_LIB="${ROCE_PILOT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_pilot_nfolds_20260917}"
PILOT_ROOT="${ROCE_PILOT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/pilot_nfolds_20260917}"
CV_THREADS="${ROCE_NUISANCE_CV_THREADS:-5}"
PARTITION="${ROCE_SLURM_PARTITION:-msismall}"
CPUS=$((2 * K * CV_THREADS))
OUT_DIR="${PILOT_ROOT}/raw"
LOG_DIR="${PILOT_ROOT}/logs"
mkdir -p "${OUT_DIR}" "${LOG_DIR}"

if [[ ! -d "${PILOT_LIB}/RoCE" ]]; then
  echo "pilot library missing: ${PILOT_LIB}/RoCE" >&2
  exit 1
fi
if [[ "${N_FOLDS}" -lt "${CV_THREADS}" ]]; then
  CV_THREADS="${N_FOLDS}"
  CPUS=$((2 * K * CV_THREADS))
fi

JOB_NAME="pilot_${CONFIG}_K${K}_f${N_FOLDS}"
JOB_ID="$(sbatch --parsable \
  --job-name="${JOB_NAME}" \
  --partition="${PARTITION}" \
  --array="${SEED_START}-${SEED_END}%${MAX_CONCURRENT}" \
  --cpus-per-task="${CPUS}" \
  --mem=16G \
  --time=12:00:00 \
  --output="${LOG_DIR}/%A_%a_${JOB_NAME}.out" \
  --error="${LOG_DIR}/%A_%a_${JOB_NAME}.err" \
  --export="ALL,ROCE_PROJECT_LIB=${PILOT_LIB},ROCE_NUISANCE_CV_THREADS=${CV_THREADS},OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,R_LIBS_USER=/users/0/zhan9381/Rlibs" \
  --wrap="module load R/4.2.2-gcc-8.2.0-vp7tyde && cd ${PROJECT_ROOT} && Rscript scripts/slurm/pilot_nfolds_task.R \${SLURM_ARRAY_TASK_ID} ${CONFIG} ${K} ${N_FOLDS} ${OUT_DIR}")"
echo "submitted ${JOB_NAME} seeds ${SEED_START}-${SEED_END} (%${MAX_CONCURRENT}, ${CPUS} CPUs, ${CV_THREADS} CV threads): ${JOB_ID}"
echo "${JOB_ID} ${JOB_NAME} ${SEED_START}-${SEED_END} cpus=${CPUS}" >> "${PILOT_ROOT}/submissions.log"
