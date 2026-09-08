#!/bin/bash
#SBATCH --job-name=wald_build_smoke
#SBATCH --time=00:55:00
#SBATCH --mem=12g
#SBATCH --cpus-per-task=4
#SBATCH --output=diagnosis/face_probe/validation/wald_build_smoke_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p $HOME/Rlibs_waldon $HOME/Rlibs_waldoff

echo "===== INSTALL Wald-ON (lambda=0.5) -> ~/Rlibs_waldon ====="
ROCE_AGG_WALD_LAMBDA=0.5 R_LIBS=$HOME/Rlibs_waldon:$HOME/Rlibs \
  R CMD INSTALL --library=$HOME/Rlibs_waldon --no-multiarch --no-test-load . 2>&1 | tail -6 || { echo "[FAIL] waldon install"; exit 1; }

echo "===== INSTALL Wald-OFF (lambda=1e-6) -> ~/Rlibs_waldoff ====="
ROCE_AGG_WALD_LAMBDA=1e-6 R_LIBS=$HOME/Rlibs_waldoff:$HOME/Rlibs \
  R CMD INSTALL --library=$HOME/Rlibs_waldoff --no-multiarch --no-test-load . 2>&1 | tail -6 || { echo "[FAIL] waldoff install"; exit 1; }

echo "===== verify baked lambda in each lib (expect 0.5 and 1e-06) ====="
R_LIBS=$HOME/Rlibs_waldon:$HOME/Rlibs  Rscript -e 'library(RoCE); cat("  waldon  AGG_WALD_LAMBDA =", RoCE:::AGG_WALD_LAMBDA, "\n")'
R_LIBS=$HOME/Rlibs_waldoff:$HOME/Rlibs Rscript -e 'library(RoCE); cat("  waldoff AGG_WALD_LAMBDA =", RoCE:::AGG_WALD_LAMBDA, "\n")'

echo "===== SMOKE: n_sims=10, rho in {0,2.5}, C1 p=50 K=4 ====="
export ROCE_NSIMS=10 ROCE_PROBE_P=50 ROCE_PROBE_CONFIG=C1 ROCE_PROBE_K=4 ROCE_RHO_GRID="0,2.5"

echo "----- Wald-ON smoke -----"
ROCE_WALD_TAG=on_smoke R_LIBS=$HOME/Rlibs_waldon:$HOME/Rlibs \
  Rscript diagnosis/face_probe/face_negtransfer_waldprobe.R 2>&1 | grep -E 'AGG_WALD|one_round|target_only|bias|coverage|ci_width|DONE|RoCE from'

echo "----- Wald-OFF smoke -----"
ROCE_WALD_TAG=off_smoke R_LIBS=$HOME/Rlibs_waldoff:$HOME/Rlibs \
  Rscript diagnosis/face_probe/face_negtransfer_waldprobe.R 2>&1 | grep -E 'AGG_WALD|one_round|target_only|bias|coverage|ci_width|DONE|RoCE from'

echo "===== SMOKE COMPLETE ====="