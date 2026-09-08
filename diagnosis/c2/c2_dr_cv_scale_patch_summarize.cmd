#!/bin/bash
#SBATCH --job-name=FACE-c2-cvsum
#SBATCH --output=diagnosis/c2/dr_cv_scale_patch_probe/log/%x_%j.out
#SBATCH --error=diagnosis/c2/dr_cv_scale_patch_probe/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/c2/dr_cv_scale_patch_probe/log diagnosis/c2/dr_cv_scale_patch_probe/summary

export C2_RBAL_OUTPUT_ROOT="${C2_RBAL_OUTPUT_ROOT:-diagnosis/c2/dr_cv_scale_patch_probe}"

echo "======================================================"
echo " RoCE C2 DR-CV Scale Patch Summary"
echo "======================================================"
echo " Output root: ${C2_RBAL_OUTPUT_ROOT}"
echo " Start:       $(date)"
echo " Host:        $(hostname)"
echo " Job:         ${SLURM_JOB_ID}"
echo "======================================================"

python3 diagnosis/c2/c2_residual_balance_summarize.py

echo "======================================================"
echo " Summary completed at $(date)"
echo "======================================================"
