#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=log/test_%j.out
#SBATCH --error=log/test_%j.err
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --job-name=FACEC_test
#SBATCH -p msismall,amdsmall,agsmall
#SBATCH --nice=5

# ============================================================================
# FACE-C Test / Build SLURM Batch Script
# ============================================================================
# Paired with test.sh; invoked via sbatch with --export carrying:
#   COMPILE_ONLY   TRUE|FALSE (default FALSE) — skip the test suite
#   TEST_FILTER    regex passed to testthat::test_local(filter=) (may be empty)
#   TEST_REPORTER  testthat reporter name (default: summary)
#   USE_INSTALLED  TRUE|FALSE (default FALSE) — test against installed FACEC
# ============================================================================

# Change directory to project directory (where sbatch was invoked)
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Load R module
module load R

export R_LIBS_USER="${R_LIBS_USER:-${HOME}/Rlibs}"

# Prevent BLAS/OpenMP oversubscription during compilation and tests
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# ============================================================================
# Parameters (exported by test.sh)
# ============================================================================
COMPILE_ONLY="${COMPILE_ONLY:-false}"
TEST_FILTER="${TEST_FILTER:-}"
TEST_REPORTER="${TEST_REPORTER:-summary}"
USE_INSTALLED="${USE_INSTALLED:-false}"

# Honour the test helper's switch between load_all and installed package
if [[ "${USE_INSTALLED}" == "true" || "${USE_INSTALLED}" == "TRUE" ]]; then
    export FACEC_TEST_INSTALLED=1
else
    unset FACEC_TEST_INSTALLED
fi

mkdir -p log

echo "=============================================="
echo "FACE-C Test / Build Job"
echo "=============================================="
echo "Job ID:        ${SLURM_JOB_ID}"
echo "Mode:          $([[ "${COMPILE_ONLY}" == "true" ]] && echo "compile only" || echo "build + test")"
echo "Test filter:   ${TEST_FILTER:-<none>}"
echo "Test reporter: ${TEST_REPORTER}"
echo "Use installed: ${USE_INSTALLED}"
echo "Start time:    $(date)"
echo "=============================================="

# ============================================================================
# 1. Build the package (triggers C++ compilation via Rcpp)
# ============================================================================
echo ""
echo "[$(date)] Building FACEC package (R CMD INSTALL --preclean)"
R CMD INSTALL --preclean --no-docs --no-help --no-demo .
INSTALL_EXIT=$?

if [[ ${INSTALL_EXIT} -ne 0 ]]; then
    echo "[$(date)] R CMD INSTALL failed with exit code ${INSTALL_EXIT}"
    exit ${INSTALL_EXIT}
fi
echo "[$(date)] Build succeeded."

if [[ "${COMPILE_ONLY}" == "true" ]]; then
    echo ""
    echo "=============================================="
    echo "Compile-only run complete at: $(date)"
    echo "=============================================="
    exit 0
fi

# ============================================================================
# 2. Run the testthat suite
# ============================================================================
echo ""
echo "[$(date)] Running testthat test suite"

# Rscript invocation constructs the filter/reporter arguments from env vars.
# Using testthat::test_local() is consistent with tests/testthat/helper-load.R,
# which prefers devtools::load_all() against the source tree unless
# FACEC_TEST_INSTALLED is set.
Rscript - <<'RSCRIPT'
suppressPackageStartupMessages({
    if (!requireNamespace("testthat", quietly = TRUE)) {
        stop("testthat is required to run the FACE-C test suite")
    }
})

test_filter   <- Sys.getenv("TEST_FILTER", "")
test_reporter <- Sys.getenv("TEST_REPORTER", "summary")

# stop_on_failure=TRUE: testthat runs the full suite, prints the summary, and
# raises at the end if any expectation failed or any test errored out. With
# Rscript, that raised error propagates as a non-zero exit code, which SLURM
# reports as job failure; silent passes leave exit 0.
args <- list(reporter = test_reporter, stop_on_failure = TRUE)
if (nzchar(test_filter)) args$filter <- test_filter

cat(sprintf("[R] testthat::test_local(filter = %s, reporter = %s)\n",
            if (nzchar(test_filter)) shQuote(test_filter) else "NULL",
            shQuote(test_reporter)))

do.call(testthat::test_local, args)
cat("[R] All tests passed.\n")
RSCRIPT

TEST_EXIT=$?

echo ""
echo "=============================================="
echo "Test job completed at: $(date)"
echo "Exit code: ${TEST_EXIT}"
echo "=============================================="

exit ${TEST_EXIT}
