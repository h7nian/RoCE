#!/bin/bash
#SBATCH --job-name=hd_timing
#SBATCH --time=02:00:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=4
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/roce_timing_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript -e '
library(RoCE)
t0<-Sys.time()
r<-run_single_simulation(sim_id=1, n_total=3000, K=2, p=50, config="C1", outcome_type="binary",
   estimand_type="superpopulation", n_folds=5, dgp_type="face", ate_deviation=1.0,
   n_deviated_sites=1, estimate_ate=TRUE, verbose=FALSE, n_cores_internal=1)
cat(sprintf("p=50 ONE sim (all baselines+ATE, n_folds=5): %.1f min, %d method-rows\n",
   as.numeric(difftime(Sys.time(),t0,units="mins")), nrow(r)))
cat("methods:", paste(unique(r$method),collapse=", "), "\n")
'
