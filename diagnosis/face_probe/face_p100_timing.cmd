#!/bin/bash
#SBATCH --job-name=p100_time
#SBATCH --time=16:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=1
#SBATCH --output=diagnosis/face_probe/validation/face_p100_timing_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript -e '
library(RoCE); t0<-Sys.time()
r<-run_single_simulation(sim_id=1, n_total=3000, K=2, p=100, config="C1", outcome_type="binary",
   estimand_type="superpopulation", n_folds=5, dgp_type="face", ate_deviation=1.0,
   n_deviated_sites=1, estimate_ate=TRUE, verbose=FALSE, n_cores_internal=1)
cat(sprintf("p=100 ONE sim, 1 core, n_folds=5, all baselines+ATE: %.1f min\n",
   as.numeric(difftime(Sys.time(),t0,units="mins"))))'
