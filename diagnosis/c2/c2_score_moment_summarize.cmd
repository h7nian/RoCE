#!/bin/bash
#SBATCH --job-name=FACE-c2-score-sum
#SBATCH --output=diagnosis/c2/score_moment_audit/log/%x_%j.out
#SBATCH --error=diagnosis/c2/score_moment_audit/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2g
#SBATCH --time=00:10:00
#SBATCH -p msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_SCORE_AUDIT_OUTPUT_ROOT="${C2_SCORE_AUDIT_OUTPUT_ROOT:-diagnosis/c2/score_moment_audit}"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/summary"

echo "======================================================"
echo " FACE-HD C2 Score/Moment Summary"
echo "======================================================"
echo " Repo:        $(pwd)"
echo " Output root: ${C2_SCORE_AUDIT_OUTPUT_ROOT}"
echo " Start:       $(date)"
echo " Host:        $(hostname)"
echo " Job:         ${SLURM_JOB_ID:-NA}"
echo "======================================================"

python3 diagnosis/c2/c2_score_moment_summarize.py

echo "======================================================"
echo " Summary completed at $(date)"
echo "======================================================"
