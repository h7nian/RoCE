#!/bin/bash
# HISTORY: 2026-09-10 #0011 Precision of the weight rules at production scale
# Task:    500-replication precision study of weight rules A and B at one
#          (config, K, rho) cell, using the production simulation path.
#          Array task = one seed block. Measured cost is ~2.5-3 CPU-hours per
#          K = 4, p = 100 replicate, so a 100-replicate block on 20 CPUs does
#          not fit in 12 h; use ROCE_RULE_BLOCK=25 (about 4 h). Replicates are
#          checkpointed individually, so a resubmitted task resumes.
#          Env: ROCE_RULE_CONFIG, ROCE_RULE_K, ROCE_RULE_RHO, ROCE_RULE_BLOCK.

#SBATCH --job-name=roce_rule_precision
#SBATCH --time=12:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=20
#SBATCH --partition=msismall
#SBATCH --output=diagnosis/logs/%x-%A_%a.out
#SBATCH --error=diagnosis/logs/%x-%A_%a.err
#SBATCH --requeue

set -euo pipefail
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"
test -f DESCRIPTION || { echo "submit from the project root" >&2; exit 1; }

module load R/4.2.2-gcc-8.2.0-vp7tyde
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

CONFIG="${ROCE_RULE_CONFIG:?ROCE_RULE_CONFIG is required}"
K="${ROCE_RULE_K:?ROCE_RULE_K is required}"
RHO="${ROCE_RULE_RHO:?ROCE_RULE_RHO is required}"
BLOCK="${ROCE_RULE_BLOCK:-25}"
OFFSET=$(( (SLURM_ARRAY_TASK_ID - 1) * BLOCK ))

TASK="weight_rule_precision"
mkdir -p "diagnosis/logs" "diagnosis/out/${TASK}"
Rscript "diagnosis/${TASK}/${TASK}.R" "diagnosis/out/${TASK}" \
  "${CONFIG}" "${K}" "${RHO}" "${BLOCK}" "${OFFSET}"
