#!/bin/bash
#SBATCH --job-name=FACE-pdf-extract
#SBATCH --output=diagnosis/method_alignment/audit_results/%x_%j.out
#SBATCH --error=diagnosis/method_alignment/audit_results/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=00:20:00
#SBATCH -p msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/method_alignment/audit_results

echo "======================================================"
echo " RoCE PDF Extraction"
echo "======================================================"
echo " Repo:  $(pwd)"
echo " Start: $(date)"
echo " Host:  $(hostname)"
echo " Job:   ${SLURM_JOB_ID:-NA}"
echo "======================================================"

python3 diagnosis/method_alignment/extract_method_pdfs.py --docs-dir docs

echo "======================================================"
echo " PDF extraction completed at $(date)"
echo "======================================================"
