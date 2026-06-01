#!/bin/bash
#SBATCH --job-name=FACE-c2-nsum
#SBATCH --output=diagnosis/c2/nuisance_refit_probe/log/%x_%j.out
#SBATCH --error=diagnosis/c2/nuisance_refit_probe/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/c2/nuisance_refit_probe/log diagnosis/c2/nuisance_refit_probe/summary

echo "======================================================"
echo " FACE-HD C2 Nuisance Refit Summary"
echo "======================================================"
echo " Start: $(date)"
echo " Host:  $(hostname)"
echo " Job:   ${SLURM_JOB_ID}"
echo "======================================================"

python3 diagnosis/c2/c2_nuisance_refit_summarize.py

echo "======================================================"
echo " Summary completed at $(date)"
echo "======================================================"
