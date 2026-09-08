#!/bin/bash
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
candidate_root=results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1
test ! -e "${candidate_root}/lib"
mkdir "${candidate_root}/lib"
R CMD INSTALL --preclean --library="${candidate_root}/lib" "${candidate_root}/RoCE"
Rscript -e '.libPaths(c("results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/lib",.libPaths())); library(RoCE); x<-testthat::test_file("diagnosis/tate_common_weight/test_outcome_cv_training_scale.R",reporter="summary"); d<-as.data.frame(x); stopifnot(sum(d$failed)==0,!any(d$error),!any(d$skipped)); cat(sum(d$passed),"targeted assertions passed\n")'
