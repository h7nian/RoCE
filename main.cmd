#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=/scratch.global/zhan9381/FACE-HD/%A_%a.out
#SBATCH --error=/scratch.global/zhan9381/FACE-HD/%A_%a.err
#SBATCH --time=2-00:00:00
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
#SBATCH --job-name=RoCE
#SBATCH -p preempt,msismall,agsmall,amdsmall,amd512,ag2tb,msibigmem,saffo-2tb
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
    local checkpoint_dir="${CHECKPOINT_DIR:-/scratch.global/${USER}/FACE-HD/runs/checkpoints}"
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
# Pin to the spack/centos7-ivybridge build that the installed RoCE.so was
# compiled against. On Rocky 8 nodes (e.g. saffo-2tb's acl4x), an unqualified
# `module load R` can resolve to a different build, breaking ABI on RoCE.so.
module load R/4.2.2-gcc-8.2.0-vp7tyde

# Override unconditionally: the spack `R/4.2.2-gcc-8.2.0-vp7tyde` module sets
# its own R_LIBS_USER (~/R/library), which is empty on this account. All 270
# installed packages — including RoCE and its deps — live in ~/Rlibs.
export R_LIBS_USER="${HOME}/Rlibs"

# Run FACE core from the installed RoCE package by default. Set
# ROCE_MAIN_USE_SOURCE=TRUE explicitly when a source-tree run is needed.
export ROCE_MAIN_USE_SOURCE="${ROCE_MAIN_USE_SOURCE:-FALSE}"

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
# arg11: shift_strength (numeric) - default: 0.5
# arg12: n_folds (integer) - default: 10
#        NOTE: diagnosis/bias/ recommends n_folds=3 (or 5) when
#        n_site_arm / (n_folds * p) < ~5; see comment in main.R.
#        Verified C1 at n=5000,K=3,p=10: K_f=3 raised one/two-round
#        coverage from 0.886/0.878 to 0.945/0.955.
# arg13: n_sims (integer) - default: 500
# arg14: dgp_type (bounded, roce, face) - default: bounded
# arg15: ate_deviation (numeric, FACE paper only) - default: 0.0
# arg16: n_deviated_sites (integer, FACE paper only) - default: 0
# arg17: use_lambda_cache (TRUE/FALSE) - default: TRUE
# arg18: verbose_every (integer, sequential log interval) - default: 10
# arg19: parallel_strategy (outer_priority|balanced|outer_only) - default: outer_priority
# arg20: estimate_ate (TRUE/FALSE) - default: FALSE
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
SHIFT_STRENGTH="${arg11:-0.5}"
N_FOLDS="${arg12:-10}"
N_SIMS="${arg13:-500}"
DGP_TYPE="${arg14:-bounded}"
ATE_DEVIATION="${arg15:-0.0}"
N_DEVIATED_SITES="${arg16:-0}"
USE_LAMBDA_CACHE="${arg17:-TRUE}"
VERBOSE_EVERY="${arg18:-10}"
PARALLEL_STRATEGY="${arg19:-outer_priority}"
ESTIMATE_ATE="${arg20:-FALSE}"
if [[ "$DGP_TYPE" == "bounded" && -z "${arg20:-}" ]]; then ESTIMATE_ATE="TRUE"; fi

normalize_bool() {
    case "$1" in
        1|true|TRUE|True|t|T|yes|YES|Yes|y|Y) echo "TRUE" ;;
        0|false|FALSE|False|f|F|no|NO|No|n|N) echo "FALSE" ;;
        *) return 1 ;;
    esac
}

if ! USE_LAMBDA_CACHE=$(normalize_bool "$USE_LAMBDA_CACHE"); then
    echo "ERROR: use_lambda_cache must be true/false, yes/no, or 1/0."
    exit 2
fi
if ! ESTIMATE_ATE=$(normalize_bool "$ESTIMATE_ATE"); then
    echo "ERROR: estimate_ate must be true/false, yes/no, or 1/0."
    exit 2
fi
if [[ ! "$PARALLEL_STRATEGY" =~ ^(outer_priority|balanced|outer_only)$ ]]; then
    echo "ERROR: parallel_strategy must be one of: outer_priority, balanced, outer_only."
    exit 2
fi
if [[ ! "$DGP_TYPE" =~ ^(bounded|roce|face)$ ]]; then
    echo "ERROR: dgp_type must be one of: bounded, roce, face."
    exit 2
fi
if [[ ! "$OUTCOME_TYPE" =~ ^(binary|continuous)$ ]]; then
    echo "ERROR: outcome_type must be one of: binary, continuous."
    exit 2
fi
if [[ ! "$ESTIMAND_TYPE" =~ ^(superpopulation|sample)$ ]]; then
    echo "ERROR: estimand_type must be one of: superpopulation, sample."
    exit 2
fi
if ! [[ "$N_FOLDS" =~ ^[0-9]+$ ]] || [[ "$N_FOLDS" -lt 3 ]]; then
    echo "ERROR: n_folds must be an integer >= 3."
    exit 2
fi
if ! [[ "$N_SIMS" =~ ^[0-9]+$ ]] || [[ "$N_SIMS" -lt 1 ]]; then
    echo "ERROR: n_sims must be a positive integer."
    exit 2
fi
if ! [[ "$VERBOSE_EVERY" =~ ^[0-9]+$ ]] || [[ "$VERBOSE_EVERY" -lt 1 ]]; then
    echo "ERROR: verbose_every must be a positive integer."
    exit 2
fi

# Create setting identifier from active DGP/output parameters.
# Format (roce): n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_<alloc>_tf<transform>_ht<het>_ss<shift>_kf<folds>_ate<flag>
# Format (face):  n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_face_dev<dev>_nd<nd>_kf<folds>_ate<flag>
# Must match main.R and main.sh::build_setting_id()
if [[ -n "${SETTING_ID_FROM_FILE:-}" ]]; then
    SETTING_ID="${SETTING_ID_FROM_FILE}"
elif [[ "$DGP_TYPE" == "face" || "$DGP_TYPE" == "bounded" ]]; then
    SETTING_ID="n${arg1}_K${arg2}_p${arg3}_${arg4}_${ESTIMAND_TYPE}_${OUTCOME_TYPE}_${DGP_TYPE}_dev${ATE_DEVIATION}_nd${N_DEVIATED_SITES}_kf${N_FOLDS}_ate${ESTIMATE_ATE}"
else
    SETTING_ID="n${arg1}_K${arg2}_p${arg3}_${arg4}_${ESTIMAND_TYPE}_${OUTCOME_TYPE}_${SITE_ALLOCATION}_tf${TRANSFORM_TYPE}_ht${HETEROGENEITY_TYPE}_ss${SHIFT_STRENGTH}_kf${N_FOLDS}_ate${ESTIMATE_ATE}"
fi

# Set directories with setting-specific paths
RUN_ROOT="${ROCE_RUN_ROOT:-/scratch.global/${USER}/FACE-HD/runs}"
export CHECKPOINT_DIR="${CHECKPOINT_DIR:-${RUN_ROOT}/checkpoints/${arg4}}"
export RESULTS_DIR="${RESULTS_DIR:-${RUN_ROOT}/results}"
LOG_DIR="${RUN_ROOT}/log"
mkdir -p "$CHECKPOINT_DIR" "$RESULTS_DIR" "$LOG_DIR"

# Create unique log file name using setting id
ARRAY_TAG="${SLURM_ARRAY_TASK_ID:-0}"
LOG_FILE="${LOG_DIR}/${SETTING_ID}_${SLURM_JOB_ID}_${ARRAY_TAG}.out"

echo "=============================================="
echo "RoCE Simulation Job"
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
echo "  estimate_ate:        ${ESTIMATE_ATE}"
echo "  dgp_type:            ${DGP_TYPE}"
echo "  ROCE_MAIN_USE_SOURCE: ${ROCE_MAIN_USE_SOURCE}"
if [[ "$DGP_TYPE" == "face" || "$DGP_TYPE" == "bounded" ]]; then
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
# Pass all 20 arguments
# --no-restore prevents loading stale .RData that masks updated package functions
R CMD BATCH --no-restore "--args $arg1 $arg2 $arg3 $arg4 ${SLURM_JOB_ID} ${ESTIMAND_TYPE} ${SITE_ALLOCATION} ${TRANSFORM_TYPE} ${OUTCOME_TYPE} ${HETEROGENEITY_TYPE} ${SHIFT_STRENGTH} ${N_FOLDS} ${N_SIMS} ${DGP_TYPE} ${ATE_DEVIATION} ${N_DEVIATED_SITES} ${USE_LAMBDA_CACHE} ${VERBOSE_EVERY} ${PARALLEL_STRATEGY} ${ESTIMATE_ATE}" main.R ${LOG_FILE} &

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
