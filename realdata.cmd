#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=log/realdata_%j.out
#SBATCH --error=log/realdata_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# ============================================================================
# CPU ALLOCATION:
# The RHC real-data experiment is a single run (no MC loop), but the inner
# two-level cross-fitting loop over K source sites benefits from source-
# level parallelism. realdata.R calls run_rhc_experiment(n_cores = -1L),
# which consumes cpus-per-task-1 workers. --time is set conservatively to
# 8 h to absorb rare node-level slowdowns (a typical run completes well
# under 1 h with this CPU allocation).
# ============================================================================
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --job-name=FACE_realdata
#SBATCH -p preempt,saffo-2tb,msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb

# ============================================================================
# FACE-C Real-Data SLURM Batch Script
# ============================================================================
# Paired with realdata.sh; invoked via sbatch with --export carrying:
#   K_arg               number of sites (target + sources)
#   OUTCOME             death30 | death180 | los
#   N_FOLDS             outer cross-fitting folds
#   A_VAL               target treatment level (0 or 1)
#   SITE_VAR            raw RHC column used for site partitioning
#                         (default "cat1"; alt: ninsclas, income, race, ca)
#   M_TAU_INFERENCE     inference-time truncation cap on phi^T gamma
#                         (default "3"; accepts "Inf")
#   SITE_PRESET         optional site-level recode preset (empty for none).
#                         Known values: "ninsclas4" — drop "No insurance",
#                         merge "Medicaid" + "Medicare & Medicaid" into
#                         "Medicaid_any"; use with SITE_VAR=ninsclas.
#   PHI                 covariate basis map: "identity" (default) or "bspline"
#                         (cubic B-spline with knots=5 on each continuous
#                         variable; pushes p from ~75 to ~200).
#   TARGET_SITE         optional override (empty for default = largest level)
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Pin to the spack/centos7-ivybridge build that FACEC.so is compiled against;
# unqualified `module load R` can resolve to a different build on Rocky 8 nodes.
module load R/4.2.2-gcc-8.2.0-vp7tyde

# Override unconditionally: the spack `R/4.2.2-gcc-8.2.0-vp7tyde` module sets
# its own R_LIBS_USER (~/R/library), which is empty on this account. All 270
# installed packages — including FACEC and its deps — live in ~/Rlibs.
export R_LIBS_USER="${HOME}/Rlibs"

# Prevent BLAS/OpenMP oversubscription when using R-level parallelism.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# Performance knobs matching main.cmd; harmless when ignored.
export NLAMBDA_INIT=100
export NESTED_PARALLEL=1

# ----------------------------------------------------------------------------
# Parameters (with defaults if the exports are absent — enables running
# realdata.cmd directly via `sbatch realdata.cmd` during debugging).
# ----------------------------------------------------------------------------
K_arg="${K_arg:-5}"
OUTCOME="${OUTCOME:-death30}"
N_FOLDS="${N_FOLDS:-10}"
A_VAL="${A_VAL:-1}"
SITE_VAR="${SITE_VAR:-cat1}"
M_TAU_INFERENCE="${M_TAU_INFERENCE:-3}"
SITE_PRESET="${SITE_PRESET:-}"
PHI="${PHI:-identity}"
TARGET_SITE="${TARGET_SITE:-}"

# Export SITE_PRESET and PHI for Rscript / R CMD BATCH to read via Sys.getenv().
export SITE_PRESET
export PHI

case "${M_TAU_INFERENCE}" in
    Inf|inf|Infinity|infinity) M_TAU_TAG="Mtinf" ;;
    *)                          M_TAU_TAG="Mt${M_TAU_INFERENCE/./p}" ;;
esac

SETTING_ID="rhc_K${K_arg}_${OUTCOME}_kf${N_FOLDS}_A${A_VAL}_${SITE_VAR}_${M_TAU_TAG}"
if [[ -n "${SITE_PRESET}" ]]; then
    SETTING_ID="${SETTING_ID}_${SITE_PRESET}"
fi
if [[ "${PHI}" != "identity" ]]; then
    SETTING_ID="${SETTING_ID}_${PHI}"
fi
LOG_FILE="log/${SETTING_ID}_${SLURM_JOB_ID}.out"

mkdir -p log results/real_data

echo "=============================================="
echo "FACE-C Real-Data (RHC) Job"
echo "=============================================="
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Setting ID:  ${SETTING_ID}"
echo "K:           ${K_arg}"
echo "Outcome:     ${OUTCOME}"
echo "n_folds:     ${N_FOLDS}"
echo "A_val:       ${A_VAL}"
echo "Site var:    ${SITE_VAR}"
echo "M_tau_inf:   ${M_TAU_INFERENCE}"
echo "Site preset: ${SITE_PRESET:-<none>}"
echo "Phi basis:   ${PHI}"
echo "Target site: ${TARGET_SITE:-<auto>}"
echo "Start time:  $(date)"
echo "=============================================="

# ----------------------------------------------------------------------------
# Verify the RHC CSV is vendored. If absent, bail out early with a helpful
# message rather than letting realdata.R discover the missing file deep in
# the cohort builder.
# ----------------------------------------------------------------------------
if [[ ! -s inst/extdata/rhc.csv ]]; then
    echo "ERROR: inst/extdata/rhc.csv not found or empty."
    echo "       Run ./data/download_rhc.sh first to fetch the CSV."
    exit 2
fi

# ----------------------------------------------------------------------------
# FACEC is expected to be installed already via `./test.sh --compile-only`
# (or an equivalent out-of-band `R CMD INSTALL`). Skipping the in-job install
# lets multiple realdata experiments run concurrently without racing on the
# $R_LIBS_USER/00LOCK-FACE-C directory.
#
# Sanity-check that FACEC is importable before the experiment begins, so a
# stale / missing install is surfaced immediately instead of at the first
# `library(FACEC)` call deep inside realdata.R.
# ----------------------------------------------------------------------------
echo ""
echo "[$(date)] Sanity check: FACEC package loadable?"
Rscript -e 'if (!requireNamespace("FACEC", quietly = FALSE)) {
  stop("FACEC is not installed; run ./test.sh --compile-only first.")
} else {
  cat(sprintf("FACEC %s loaded.\n", as.character(utils::packageVersion("FACEC"))))
}'
LOAD_EXIT=$?
if [[ ${LOAD_EXIT} -ne 0 ]]; then
    echo "[$(date)] FACEC load check failed with exit code ${LOAD_EXIT}"
    exit ${LOAD_EXIT}
fi

# ----------------------------------------------------------------------------
# Run realdata.R — positional args match the parser in realdata.R.
# Order: K outcome n_folds A_val job_id site_var M_tau_inference target_site.
# target_site is the single truly optional trailing argument, so empty shell
# expansion cannot shift earlier (always-present) arguments.
# --no-restore prevents stale .RData from masking updated package functions.
# ----------------------------------------------------------------------------
R CMD BATCH --no-restore \
    "--args ${K_arg} ${OUTCOME} ${N_FOLDS} ${A_VAL} ${SLURM_JOB_ID} ${SITE_VAR} ${M_TAU_INFERENCE} ${TARGET_SITE}" \
    realdata.R "${LOG_FILE}"

R_EXIT_CODE=$?

echo "=============================================="
echo "Real-data job completed at: $(date)"
echo "Exit code: ${R_EXIT_CODE}"
echo "=============================================="

exit ${R_EXIT_CODE}
