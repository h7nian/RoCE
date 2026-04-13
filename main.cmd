#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=log/%A_%a.out
#SBATCH --error=log/%A_%a.err
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# ============================================================================
# CPU ALLOCATION GUIDE:
# With NESTED_PARALLEL=1, cores are split: outer (simulation) × inner (source)
# Example: 64 cores → 32 outer × 2 inner for K=3-5 source sites
# Memory: ~4GB per core for large datasets (p=200)
# ============================================================================
#SBATCH --cpus-per-task=64
#SBATCH --mem=96g
#SBATCH --job-name=FACE-C
#SBATCH -p ag2tb,amd2tb,agsmall,amdsmall,msismall,msibigmem,msilong,saffo-2tb
#SBATCH --nice=5
#SBATCH --requeue
#SBATCH --signal=B:USR1@120

# ============================================================================
# Checkpoint/Restart Support for Preempt Partition
# ============================================================================
# This script supports automatic requeue on preemption with checkpoint saving.
# When preempted, SLURM sends SIGUSR1 120 seconds before termination.
# The signal handler saves the current state and requests a requeue.
# ============================================================================

# Trap SIGUSR1 for graceful shutdown with checkpoint
# When preempted, SLURM sends this signal before killing the job
handle_preemption() {
    echo "$(date): Received preemption signal (SIGUSR1). Saving checkpoint..."
    local checkpoint_dir="${CHECKPOINT_DIR:-${SLURM_SUBMIT_DIR}/checkpoints}"
    mkdir -p "${checkpoint_dir}"
    # Touch a flag file to signal R to save checkpoint immediately
    touch "${checkpoint_dir}/.preempt_signal_${SLURM_JOB_ID}"
    # Wait for R to save checkpoint (max 100 seconds)
    for i in {1..100}; do
        if [ -f "${checkpoint_dir}/.checkpoint_saved_${SLURM_JOB_ID}" ]; then
            echo "$(date): Checkpoint saved successfully."
            rm -f "${checkpoint_dir}/.checkpoint_saved_${SLURM_JOB_ID}"
            break
        fi
        sleep 1
    done
    # Remove preempt signal to avoid immediate exit on requeue
    rm -f "${checkpoint_dir}/.preempt_signal_${SLURM_JOB_ID}"
    echo "$(date): Requesting requeue..."
    scontrol requeue ${SLURM_JOB_ID}
    exit 0
}

# Set trap for preemption signal
trap 'handle_preemption' USR1

# Change directory to project directory (where sbatch was invoked)
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Load required modules
module load R

export R_LIBS_USER="${R_LIBS_USER:-${HOME}/Rlibs}"

# Prevent BLAS/OpenMP oversubscription when using R-level parallelism
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# ============================================================================
# PERFORMANCE OPTIMIZATION
# ============================================================================
# NLAMBDA_INIT: Number of lambda values for initial outcome model CV (default 100)
# Higher values give finer lambda selection, reducing nuisance estimation bias.
# NESTED_PARALLEL: Enable both simulation-level and source-site parallelism
# ============================================================================
export NLAMBDA_INIT=100
export NESTED_PARALLEL=1

# ============================================================================
# PARAMETERS (passed from main.sh via --export)
# arg1: n_total (sample size)
# arg2: K (number of source sites)
# arg3: p (number of covariates)
# arg4: config (C1, C2, C3, C4)
# arg5: (reserved for job_id, auto-set to SLURM_JOB_ID)
# arg6: estimand_type (superpopulation, sample) - default: superpopulation
# arg7: site_allocation (model, uniform, balanced) - default: model
# arg8: transform_type (strong, mild, none) - default: mild
# arg9: outcome_type (binary, continuous) - default: binary
# arg10: heterogeneity_type (none, mild, strong, partial) - default: none
# arg11: shift_strength (numeric) - default: 1.0
# arg12: n_folds (integer) - default: 10
# arg13: n_sims (integer) - default: 500
# arg14: dgp_type (facec, face) - default: facec
# arg15: ate_deviation (numeric, FACE paper only) - default: 0.0
# arg16: n_deviated_sites (integer, FACE paper only) - default: 0
# arg17: use_lambda_cache (TRUE/FALSE) - default: TRUE
# arg18: verbose_every (integer, sequential log interval) - default: 10
# arg19: parallel_strategy (outer_priority|balanced|outer_only) - default: outer_priority
# ============================================================================

# Optional job-array combo file mode
if [[ -n "${COMBO_FILE:-}" ]]; then
    if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        echo "ERROR: COMBO_FILE is set but SLURM_ARRAY_TASK_ID is missing."
        exit 2
    fi
    if [[ ! -f "${COMBO_FILE}" ]]; then
        echo "ERROR: COMBO_FILE not found: ${COMBO_FILE}"
        exit 2
    fi

    combo_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COMBO_FILE}")
    if [[ -z "${combo_line}" ]]; then
        echo "ERROR: No combo entry for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} in ${COMBO_FILE}"
        exit 2
    fi

    IFS=$'\t' read -r arg1 arg2 arg3 arg4 arg9 arg11 arg15 arg16 SETTING_ID_FROM_FILE <<< "${combo_line}"
fi

# Set defaults for optional data generation parameters
ESTIMAND_TYPE="${arg6:-superpopulation}"
SITE_ALLOCATION="${arg7:-model}"
TRANSFORM_TYPE="${arg8:-mild}"
OUTCOME_TYPE="${arg9:-binary}"
HETEROGENEITY_TYPE="${arg10:-none}"
SHIFT_STRENGTH="${arg11:-1.0}"
N_FOLDS="${arg12:-10}"
N_SIMS="${arg13:-500}"
DGP_TYPE="${arg14:-facec}"
ATE_DEVIATION="${arg15:-0.0}"
N_DEVIATED_SITES="${arg16:-0}"
USE_LAMBDA_CACHE="${arg17:-TRUE}"
VERBOSE_EVERY="${arg18:-10}"
PARALLEL_STRATEGY="${arg19:-outer_priority}"

# Create setting identifier — always includes ALL parameters for full traceability.
# Format (facec):       n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_<alloc>_tf<transform>_ht<het>_ss<shift>_kf<folds>
# Format (face):  n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_face_dev<dev>_nd<nd>_kf<folds>
# Must match main.R and main.sh::build_setting_id()
if [[ -n "${SETTING_ID_FROM_FILE:-}" ]]; then
    SETTING_ID="${SETTING_ID_FROM_FILE}"
elif [[ "$DGP_TYPE" == "face" ]]; then
    SETTING_ID="n${arg1}_K${arg2}_p${arg3}_${arg4}_${ESTIMAND_TYPE}_${OUTCOME_TYPE}_face_dev${ATE_DEVIATION}_nd${N_DEVIATED_SITES}_kf${N_FOLDS}"
else
    SETTING_ID="n${arg1}_K${arg2}_p${arg3}_${arg4}_${ESTIMAND_TYPE}_${OUTCOME_TYPE}_${SITE_ALLOCATION}_tf${TRANSFORM_TYPE}_ht${HETEROGENEITY_TYPE}_ss${SHIFT_STRENGTH}_kf${N_FOLDS}"
fi

# Set directories with setting-specific paths
export CHECKPOINT_DIR="checkpoints/${arg4}"
export RESULTS_DIR="results"
mkdir -p ${CHECKPOINT_DIR} ${RESULTS_DIR} log

# Create unique log file name using setting id
ARRAY_TAG="${SLURM_ARRAY_TASK_ID:-0}"
LOG_FILE="log/${SETTING_ID}_${SLURM_JOB_ID}_${ARRAY_TAG}.out"

echo "=============================================="
echo "FACE-C Simulation Job"
echo "=============================================="
echo "Job ID:        ${SLURM_JOB_ID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Setting ID:    ${SETTING_ID}"
echo "Parameters:"
echo "  n_total:             ${arg1}"
echo "  K:                   ${arg2}"
echo "  p:                   ${arg3}"
echo "  config:              ${arg4}"
echo "Data Generation:"
echo "  site_allocation:     ${SITE_ALLOCATION}"
echo "  estimand_type:       ${ESTIMAND_TYPE}"
echo "  transform_type:      ${TRANSFORM_TYPE}"
echo "Design:"
echo "  outcome_type:        ${OUTCOME_TYPE}"
echo "  heterogeneity_type:  ${HETEROGENEITY_TYPE}"
echo "  shift_strength:      ${SHIFT_STRENGTH}"
echo "  n_folds:             ${N_FOLDS}"
echo "  n_sims:              ${N_SIMS}"
echo "  use_lambda_cache:    ${USE_LAMBDA_CACHE}"
echo "  verbose_every:       ${VERBOSE_EVERY}"
echo "  parallel_strategy:   ${PARALLEL_STRATEGY}"
echo "  dgp_type:            ${DGP_TYPE}"
if [[ "$DGP_TYPE" == "face" ]]; then
echo "  ate_deviation:       ${ATE_DEVIATION}"
echo "  n_deviated_sites:    ${N_DEVIATED_SITES}"
fi
echo "Performance:"
echo "  NLAMBDA_INIT:    ${NLAMBDA_INIT:-100} (lambda grid size for initial outcome CV)"
echo "  NESTED_PARALLEL: ${NESTED_PARALLEL:-0} (sim + source-site parallel)"
echo "Start Time:    $(date)"
echo "Checkpoint:    ${CHECKPOINT_DIR}"
echo "Results:       ${RESULTS_DIR}"
echo "Log File:      ${LOG_FILE}"
echo "=============================================="

# Run the R script in background to allow signal handling
# Pass all 19 arguments
# --no-restore prevents loading stale .RData that masks updated package functions
R CMD BATCH --no-restore "--args $arg1 $arg2 $arg3 $arg4 ${SLURM_JOB_ID} ${ESTIMAND_TYPE} ${SITE_ALLOCATION} ${TRANSFORM_TYPE} ${OUTCOME_TYPE} ${HETEROGENEITY_TYPE} ${SHIFT_STRENGTH} ${N_FOLDS} ${N_SIMS} ${DGP_TYPE} ${ATE_DEVIATION} ${N_DEVIATED_SITES} ${USE_LAMBDA_CACHE} ${VERBOSE_EVERY} ${PARALLEL_STRATEGY}" main.R ${LOG_FILE} &

# Wait for R process to complete
wait $!
R_EXIT_CODE=$?

echo "=============================================="
echo "Job completed at: $(date)"
echo "Exit code: ${R_EXIT_CODE}"
echo "=============================================="

# Clean up signal file if exists
rm -f "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}"

exit ${R_EXIT_CODE}
