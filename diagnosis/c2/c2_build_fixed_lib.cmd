#!/bin/bash
#SBATCH --job-name=FACE-build-fix
#SBATCH --output=diagnosis/c2/fix_validation/log/build_%j.out
#SBATCH --error=diagnosis/c2/fix_validation/log/build_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=00:40:00

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
LIB="${C2_FIX_LIB:?C2_FIX_LIB must be set (absolute path to isolated library)}"
mkdir -p "${LIB}" diagnosis/c2/fix_validation/log

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

# Dependencies (Rcpp, RcppEigen, ...) resolve from the user library; RoCE itself
# is installed fresh into the isolated LIB from the edited source (--preclean forces
# a clean recompile of the patched src/cv_utils.h + src/density_ratio.cpp).
export R_LIBS_USER="${HOME}/Rlibs"

echo "=== Building patched RoCE into ${LIB} at $(date) ==="
R CMD INSTALL --preclean --no-multiarch --library="${LIB}" .
echo "=== INSTALL done; verifying load ==="
Rscript -e ".libPaths(c('${LIB}', .libPaths())); suppressMessages(library(RoCE)); cat('RoCE loaded from fixed lib OK:', system.file(package='RoCE'), '\n')"
echo "=== BUILD COMPLETE $(date) ==="
