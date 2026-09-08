#!/bin/bash
# Shared launcher for bounded R diagnostics. Scientific arguments belong to
# the called script; this wrapper only fixes the R runtime and thread budget.
set -euo pipefail
if [[ "$#" -lt 1 || ! -f "$1" || "$1" != *.R ]]; then
  echo "usage: run_r_diagnostic.sh EXISTING_SCRIPT.R [ARGUMENTS...]" >&2
  exit 2
fi
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
command -v Rscript >/dev/null
exec Rscript "$@"
