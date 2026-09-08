#!/bin/bash
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROCE_TEST_INSTALLED=1 _R_CHECK_FORCE_SUGGESTS_=false
candidate_root=/users/0/zhan9381/FACE-HD/results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1
test ! -e "${candidate_root}/package_check"
mkdir "${candidate_root}/package_check"
cd "${candidate_root}/package_check"
R CMD build --no-build-vignettes --no-manual "${candidate_root}/RoCE"
R CMD check --no-manual --no-build-vignettes RoCE_0.1.0.tar.gz
rg -x 'Status: OK' RoCE.Rcheck/00check.log
