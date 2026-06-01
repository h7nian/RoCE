#!/bin/bash
# ============================================================================
# FACE-HD Test / Build Job Submission Script
# ============================================================================
#
# Submits a single SLURM job that (1) rebuilds the FACEHD R package (which
# triggers C++ recompilation via Rcpp), and (2) runs the testthat suite.
# Use this script whenever the C++ sources in src/ or R code in R/ change
# and you want to verify the package still builds and tests still pass on
# MSI compute nodes (never on a login node).
#
# Usage:
#   ./test.sh                          # Build + run full test suite
#   ./test.sh --compile-only           # Build the package only, skip tests
#   ./test.sh --filter weight          # Run only test files matching 'weight'
#   ./test.sh --installed              # Test existing installed FACEHD; skip reinstall
#   ./test.sh --source                 # Test source tree via load_all; skip reinstall
#   ./test.sh --env FACEHD_FOO=bar      # Forward one FACEHD_* diagnostic variable
#   ./test.sh --reporter progress      # Override testthat reporter
#   ./test.sh --dry-run                # Show what would be submitted
#   FACEHD_* env vars are forwarded to the Slurm job for gated diagnostics.
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
USE_SOURCE=false
TEST_FILTER=""
TEST_REPORTER="summary"
EXTRA_FACEHD_ENVS=()

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --compile-only             Build the package only; skip the test suite"
    echo "  --filter PATTERN           Regex filter passed to testthat::test_local(filter=)"
    echo "  --installed                Run tests against installed FACEHD and skip R CMD INSTALL"
    echo "  --source                   Run tests against source tree via devtools::load_all and skip install"
    echo "  --env FACEHD_NAME=VALUE     Forward one FACEHD_* diagnostic variable (repeatable)"
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
        --source)        USE_SOURCE=true; shift ;;
        --env)
            env_pair="${2:?Missing FACEHD_NAME=VALUE for --env}"
            if [[ "${env_pair}" != *=* ]]; then
                echo "Error: --env expects FACEHD_NAME=VALUE, got '${env_pair}'"
                exit 1
            fi
            env_name="${env_pair%%=*}"
            if [[ ! "${env_name}" =~ ^FACEHD_[A-Za-z0-9_]*$ ]]; then
                echo "Error: --env only accepts FACEHD_* variables, got '${env_name}'"
                exit 1
            fi
            EXTRA_FACEHD_ENVS+=("${env_pair}")
            shift 2
            ;;
        --reporter)      TEST_REPORTER="${2:?Missing reporter name}"; shift 2 ;;
        -h|--help)       usage ;;
        *)               echo "Error: unknown option '$1'"; usage ;;
    esac
done

if [[ "${USE_INSTALLED}" == true && "${USE_SOURCE}" == true ]]; then
    echo "Error: --installed and --source are mutually exclusive"
    exit 1
fi

# ============================================================================
# 2. Print Summary
# ============================================================================

mkdir -p log

echo "======================================================"
echo " FACE-HD Test / Build Job Submission"
echo "======================================================"
echo ""
echo " Mode:        $([[ "$COMPILE_ONLY" == true ]] && echo "compile only" || echo "build + test")"
echo " Filter:      ${TEST_FILTER:-<none>}"
echo " Reporter:    ${TEST_REPORTER}"
echo " Use installed package: ${USE_INSTALLED}"
echo " Use source tree: ${USE_SOURCE}"
if [[ "${#EXTRA_FACEHD_ENVS[@]}" -gt 0 ]]; then
    echo " Extra FACEHD env:"
    for env_pair in "${EXTRA_FACEHD_ENVS[@]}"; do
        echo "   ${env_pair}"
    done
else
    echo " Extra FACEHD env: <none>"
fi
echo " Dry run:     ${DRY_RUN}"
echo ""
echo "======================================================"

# ============================================================================
# 3. Build job environment and Submit
# ============================================================================

write_job_env_file() {
    local env_file
    env_file=$(mktemp "log/test_env_XXXXXX.sh")
    chmod 600 "${env_file}"
    {
        printf 'export COMPILE_ONLY=%q\n' "${COMPILE_ONLY}"
        printf 'export TEST_FILTER=%q\n' "${TEST_FILTER}"
        printf 'export TEST_REPORTER=%q\n' "${TEST_REPORTER}"
        printf 'export USE_INSTALLED=%q\n' "${USE_INSTALLED}"
        printf 'export USE_SOURCE=%q\n' "${USE_SOURCE}"
        while IFS='=' read -r name value; do
            if [[ "${name}" =~ ^FACEHD_[A-Za-z0-9_]*$ ]]; then
                printf 'export %s=%q\n' "${name}" "${value}"
            fi
        done < <(env)
        for env_pair in "${EXTRA_FACEHD_ENVS[@]}"; do
            name="${env_pair%%=*}"
            value="${env_pair#*=}"
            printf 'export %s=%q\n' "${name}" "${value}"
        done
    } > "${env_file}"
    printf '%s\n' "${env_file}"
}

JOB_NAME="FACEHD_test$([[ "$COMPILE_ONLY" == true ]] && echo "_build" || echo "")"

if [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "[DRY RUN] Would submit:"
    echo "  sbatch --export=FACEHD_TEST_ENV_FILE=<generated-env-file> --job-name=${JOB_NAME} test.cmd"
    echo ""
    echo "[DRY RUN] Job env would include COMPILE_ONLY/TEST_FILTER/TEST_REPORTER/USE_INSTALLED/USE_SOURCE and current FACEHD_* variables."
    echo ""
else
    ENV_FILE=$(write_job_env_file)
    submit_output=$(sbatch \
        --export="FACEHD_TEST_ENV_FILE=${ENV_FILE}" \
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
