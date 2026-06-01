#!/bin/bash
#SBATCH --job-name=FACE-c2-oracle-gsum
#SBATCH --output=diagnosis/c2/oracle_gamma_probe/log/%x_%j.out
#SBATCH --error=diagnosis/c2/oracle_gamma_probe/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_ORACLE_GAMMA_OUTPUT_ROOT="${C2_ORACLE_GAMMA_OUTPUT_ROOT:-diagnosis/c2/oracle_gamma_probe}"
mkdir -p "${C2_ORACLE_GAMMA_OUTPUT_ROOT}/log" "${C2_ORACLE_GAMMA_OUTPUT_ROOT}/summary"

python3 diagnosis/c2/c2_oracle_gamma_summarize.py
