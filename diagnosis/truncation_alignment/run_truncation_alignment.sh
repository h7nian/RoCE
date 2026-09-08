#!/bin/bash
# HISTORY: 2026-09-07 #0003 Truncation alignment of the tilting calibrated loss
# Task:    Build an isolated candidate library from the workspace package plus
#          truncation_alignment.patch, run the score check against candidate and
#          baseline libraries, the C1 identity refit (needs the #0002 C1 cell), the
#          full testthat suite and R CMD check on the candidate. Writes
#          diagnosis/out/truncation_alignment/. Checks about 2 h, tests and
#          R CMD check about 4 h on 10 cores (submit with --time=10:00:00).
#          Usage: run_truncation_alignment.sh [all|tests]; "tests" rebuilds the
#          candidate library from the current tree and runs only step 4 (the
#          check outputs of an earlier "all" run are kept).

#SBATCH --job-name=truncation_alignment
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=10
#SBATCH --partition=msismall
#SBATCH --output=diagnosis/logs/%x-%j.out
#SBATCH --error=diagnosis/logs/%x-%j.err
#SBATCH --requeue
#SBATCH --signal=B:USR1@60

set -euo pipefail

# sbatch copies the script to a spool directory, so locate the project by the
# submission directory (submit from the project root).
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"
test -f DESCRIPTION || { echo "submit from the project root (DESCRIPTION not found in ${PROJECT_ROOT})" >&2; exit 1; }

module load R/4.2.2-gcc-8.2.0-vp7tyde
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROCE_NUISANCE_CV_THREADS=5
export _R_CHECK_FORCE_SUGGESTS_=false

PHASE="${1:-all}"
case "${PHASE}" in all|tests) ;; *) echo "usage: $0 [all|tests]" >&2; exit 1 ;; esac
TASK="truncation_alignment"
OUT="diagnosis/out/${TASK}"
BASELINE_LIB="${PROJECT_ROOT}/results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/lib"
CANDIDATE_SRC="${OUT}/RoCE"
CANDIDATE_LIB="${OUT}/lib"
mkdir -p "diagnosis/logs" "${OUT}"
if [[ "${PHASE}" == "tests" ]]; then
  # Build artifacts only; the check outputs (baseline/, candidate/) are kept.
  rm -rf "${CANDIDATE_SRC}" "${CANDIDATE_LIB}" "${OUT}/package_check"
fi
test ! -e "${CANDIDATE_SRC}" || { echo "candidate source already exists: ${CANDIDATE_SRC}" >&2; exit 1; }

# 1. Isolated copy of the package source with the patch applied.
mkdir -p "${CANDIDATE_SRC}"
rsync -a --exclude='*.o' --exclude='*.so' \
  DESCRIPTION NAMESPACE LICENSE README.md .Rbuildignore R src inst man tests \
  "${CANDIDATE_SRC}/"
patch -p1 -d "${CANDIDATE_SRC}" < "diagnosis/${TASK}/${TASK}.patch"
Rscript -e "Rcpp::compileAttributes('${CANDIDATE_SRC}')"

# 2. Install into an isolated library.
mkdir -p "${CANDIDATE_LIB}"
R CMD INSTALL --preclean --library="${CANDIDATE_LIB}" "${CANDIDATE_SRC}"

# 3. Score check with both libraries, C1 identity refit with the candidate.
if [[ "${PHASE}" == "all" ]]; then
  R_LIBS="${BASELINE_LIB}" Rscript "diagnosis/${TASK}/${TASK}.R" "${OUT}" baseline
  R_LIBS="${CANDIDATE_LIB}" Rscript "diagnosis/${TASK}/${TASK}.R" "${OUT}" candidate
fi

# 4. Full installed test suite and R CMD check on the candidate.
(
  cd "${CANDIDATE_SRC}"
  R_LIBS="${CANDIDATE_LIB}" ROCE_TEST_INSTALLED=1 Rscript -e '
    options(testthat.summary.max_reports = 1000L)
    results <- testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = FALSE)
    d <- as.data.frame(results)
    cat(sprintf("testthat: passed=%d failed=%d errors=%d skipped=%d\n",
                sum(d$passed), sum(d$failed), sum(d$error), sum(d$skipped)))
    stopifnot(sum(d$failed) == 0L, !any(d$error))
  '
)
mkdir -p "${OUT}/package_check"
(
  cd "${OUT}/package_check"
  R_LIBS="${CANDIDATE_LIB}" R CMD build --no-build-vignettes --no-manual "${PROJECT_ROOT}/${CANDIDATE_SRC}"
  R_LIBS="${CANDIDATE_LIB}" ROCE_TEST_INSTALLED=1 R CMD check --no-manual --no-build-vignettes RoCE_0.1.0.tar.gz
  grep -x 'Status: OK' RoCE.Rcheck/00check.log
)
echo "truncation_alignment: all steps completed"
