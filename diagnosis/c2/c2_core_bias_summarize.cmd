#!/bin/bash
#SBATCH --job-name=FACE-c2-core-bias
#SBATCH --output=diagnosis/c2/core_bias_summary/log/%x_%j.out
#SBATCH --error=diagnosis/c2/core_bias_summary/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_CORE_BIAS_OUTPUT_ROOT="${C2_CORE_BIAS_OUTPUT_ROOT:-diagnosis/c2/core_bias_summary}"
mkdir -p "${C2_CORE_BIAS_OUTPUT_ROOT}/log" "${C2_CORE_BIAS_OUTPUT_ROOT}/summary"

python3 diagnosis/c2/c2_core_bias_summarize.py
