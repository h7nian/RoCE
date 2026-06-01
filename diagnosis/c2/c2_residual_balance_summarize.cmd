#!/bin/bash
#SBATCH --job-name=FACE-c2-rbsum
#SBATCH --output=diagnosis/c2/residual_balance_probe/log/%x_%j.out
#SBATCH --error=diagnosis/c2/residual_balance_probe/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_RBAL_OUTPUT_ROOT="${C2_RBAL_OUTPUT_ROOT:-diagnosis/c2/residual_balance_probe}"
mkdir -p "${C2_RBAL_OUTPUT_ROOT}/log" "${C2_RBAL_OUTPUT_ROOT}/summary"

echo "======================================================"
echo " FACE-HD C2 Residual Balance Summary"
echo "======================================================"
echo " Start: $(date)"
echo " Host:  $(hostname)"
echo " Job:   ${SLURM_JOB_ID}"
echo " Output root: ${C2_RBAL_OUTPUT_ROOT}"
echo "======================================================"

python3 diagnosis/c2/c2_residual_balance_summarize.py

echo "======================================================"
echo " Summary completed at $(date)"
echo "======================================================"
