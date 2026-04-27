#!/bin/bash
# ============================================================================
# FACE-C Test / Build Job Submission Script
# ============================================================================
#
# Submits a single SLURM job that (1) rebuilds the FACEC R package (which
# triggers C++ recompilation via Rcpp), and (2) runs the testthat suite.
# Use this script whenever the C++ sources in src/ or R code in R/ change
# and you want to verify the package still builds and tests still pass on
# MSI compute nodes (never on a login node).
#
# Usage:
#   ./test.sh                          # Build + run full test suite
#   ./test.sh --compile-only           # Build the package only, skip tests
#   ./test.sh --filter weight          # Run only test files matching 'weight'
#   ./test.sh --installed              # Test against installed FACEC (not load_all)
#   ./test.sh --reporter progress      # Override testthat reporter
#   ./test.sh --dry-run                # Show what would be submitted
#
# Directory layout:
#   log/            SLURM stdout/stderr (test_*.out / test_*.err)
#   tests/testthat/ Unit tests
#
# ============================================================================

set -euo pipefail

# ============================================================================
# 1. Parse Command-Line Arguments
# ============================================================================

DRY_RUN=false
COMPILE_ONLY=false
USE_INSTALLED=false
TEST_FILTER=""
TEST_REPORTER="summary"

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --compile-only             Build the package only; skip the test suite"
    echo "  --filter PATTERN           Regex filter passed to testthat::test_local(filter=)"
    echo "  --installed                Run tests against installed FACEC (sets FACEC_TEST_INSTALLED=1)"
    echo "  --reporter REPORTER        testthat reporter (default: summary; e.g. progress, minimal)"
    echo "  --dry-run                  Show the sbatch command without submitting"
    echo "  -h, --help                 Show this message"
    echo ""
    echo "Examples:"
    echo "  $0                         # default: build + run all tests"
    echo "  $0 --filter cross_fitting  # run only test-cross_fitting.R"
    echo "  $0 --compile-only          # just rebuild package (sanity-check C++)"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)       DRY_RUN=true; shift ;;
        --compile-only)  COMPILE_ONLY=true; shift ;;
        --filter)        TEST_FILTER="${2:?Missing filter regex}"; shift 2 ;;
        --installed)     USE_INSTALLED=true; shift ;;
        --reporter)      TEST_REPORTER="${2:?Missing reporter name}"; shift 2 ;;
        -h|--help)       usage ;;
        *)               echo "Error: unknown option '$1'"; usage ;;
    esac
done

# ============================================================================
# 2. Print Summary
# ============================================================================

mkdir -p log

echo "======================================================"
echo " FACE-C Test / Build Job Submission"
echo "======================================================"
echo ""
echo " Mode:        $([[ "$COMPILE_ONLY" == true ]] && echo "compile only" || echo "build + test")"
echo " Filter:      ${TEST_FILTER:-<none>}"
echo " Reporter:    ${TEST_REPORTER}"
echo " Use installed package: ${USE_INSTALLED}"
echo " Dry run:     ${DRY_RUN}"
echo ""
echo "======================================================"

# ============================================================================
# 3. Build --export list and Submit
# ============================================================================

EXPORT_VARS="COMPILE_ONLY=${COMPILE_ONLY},TEST_FILTER=${TEST_FILTER},TEST_REPORTER=${TEST_REPORTER},USE_INSTALLED=${USE_INSTALLED}"
JOB_NAME="FACEC_test$([[ "$COMPILE_ONLY" == true ]] && echo "_build" || echo "")"

if [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "[DRY RUN] Would submit:"
    echo "  sbatch --export=${EXPORT_VARS} --job-name=${JOB_NAME} test.cmd"
    echo ""
else
    submit_output=$(sbatch \
        --export="${EXPORT_VARS}" \
        --job-name="${JOB_NAME}" \
        test.cmd)
    echo ""
    echo "Submitted: ${submit_output}"
fi

# ============================================================================
# 4. Final Summary
# ============================================================================

echo ""
echo "======================================================"
if [[ "$DRY_RUN" == true ]]; then
    echo " DRY RUN COMPLETE — nothing submitted"
else
    echo " DONE — test/build job submitted"
fi
echo "======================================================"
echo ""
echo " Useful commands:"
echo "   squeue -u \$USER                # monitor queue"
echo "   tail -f log/test_*.out          # watch live output"
echo "   scancel -u \$USER --name=${JOB_NAME}"
echo "======================================================"
