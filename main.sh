#!/bin/bash
# ============================================================================
# FACE-C Simulation Job Submission Script
# ============================================================================
#
# Submits SLURM jobs for each parameter combination in the simulation study.
# Each job runs main.R with a unique parameter setting and supports
# checkpoint/restart for preemption-resilient execution.
#
# Usage:
#   ./main.sh                          # Submit core experiments (facec DGP)
#   ./main.sh --subset minimal         # Single test job (facec)
#   ./main.sh --subset full            # Full factorial (many jobs!)
#   ./main.sh --subset face            # FACE DGP (Han et al. JASA 2023)
#   ./main.sh --dry-run                # Preview without submitting
#   ./main.sh --dry-run --subset full  # Count full factorial jobs
#
# CLI overrides (applied AFTER subset defaults):
#   --n-sims 100                       # Override simulation count
#   --outcome-type continuous          # Override outcome type(s)
#   --n-total "5000 10000"             # Override sample sizes
#   --K "3 5"                          # Override site counts
#   --p "10 50"                        # Override covariate counts
#   --use-lambda-cache true            # Toggle lambda cache in cross-fitting
#   --verbose-every 10                 # Sequential-mode detailed log interval
#   --parallel-strategy outer_priority  # nested parallel policy
#   --estimate-ate false               # Also estimate target ATE via A=0 arm
#   --array-concurrency 16             # max concurrent array tasks
#
# Subsets:                           jobs   formula
#   minimal              1           1×1×1×1  (single facec setting)
#   core               108           3n × 3K × 4p × 3cfg
#   continuous         216           core × 2 outcome types
#   shift              324           core × 3 shift strengths
#   full               648           core × 2ot × 3ss (4p variant)
#   face               144           1n × 1K × 4p × 3cfg × 4dev × 3nd
#
# Directory layout:
#   log/            SLURM stdout/stderr per job
#   checkpoints/    per-config checkpoint .rds files
#   results/        final simulation result .csv files
#
# ============================================================================

set -euo pipefail

# ============================================================================
# 1. Parse Command-Line Arguments
# ============================================================================

DRY_RUN=false
SUBSET="core"

# CLI overrides (applied after subset defaults)
OVERRIDE_N_SIMS=""
OVERRIDE_OUTCOME_TYPE=""
OVERRIDE_N_TOTAL=""
OVERRIDE_K=""
OVERRIDE_P=""
OVERRIDE_USE_LAMBDA_CACHE=""
OVERRIDE_VERBOSE_EVERY=""
OVERRIDE_PARALLEL_STRATEGY=""
OVERRIDE_ESTIMATE_ATE=""
OVERRIDE_ARRAY_CONCURRENCY=""

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --dry-run                  Count jobs without submitting"
    echo "  --subset SUBSET            One of: minimal, core, continuous, shift, full, face"
    echo "  --n-sims N                 Override number of MC simulations (e.g. 100)"
    echo "  --outcome-type TYPE        Override outcome type: binary, continuous, or 'binary continuous'"
    echo "  --n-total 'N1 N2 ...'      Override sample sizes (space-separated, quoted)"
    echo "  --K 'K1 K2 ...'            Override site counts (space-separated, quoted)"
    echo "  --p 'P1 P2 ...'            Override covariate counts (space-separated, quoted)"
    echo "  --use-lambda-cache BOOL    Override lambda cache flag (true/false)"
    echo "  --verbose-every N          Override sequential verbose interval (integer >=1)"
    echo "  --parallel-strategy STR    one of: outer_priority, balanced, outer_only"
    echo "  --estimate-ate BOOL        Also estimate ATE via A=0 arm (true/false)"
    echo "  --array-concurrency N      Max concurrent SLURM array tasks (integer >=1)"
    echo ""
    echo "Examples:"
    echo "  $0 --subset face --outcome-type binary --n-sims 100 --dry-run"
    echo "  $0 --subset minimal --n-total '500 1000' --p '10 20'"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)         DRY_RUN=true; shift ;;
        --subset)          SUBSET="${2:?Missing subset name}"; shift 2 ;;
        --n-sims)          OVERRIDE_N_SIMS="${2:?Missing n-sims value}"; shift 2 ;;
        --outcome-type)    OVERRIDE_OUTCOME_TYPE="${2:?Missing outcome-type}"; shift 2 ;;
        --n-total)         OVERRIDE_N_TOTAL="${2:?Missing n-total values}"; shift 2 ;;
        --K)               OVERRIDE_K="${2:?Missing K values}"; shift 2 ;;
        --p)               OVERRIDE_P="${2:?Missing p values}"; shift 2 ;;
        --use-lambda-cache) OVERRIDE_USE_LAMBDA_CACHE="${2:?Missing use-lambda-cache value}"; shift 2 ;;
        --verbose-every)   OVERRIDE_VERBOSE_EVERY="${2:?Missing verbose-every value}"; shift 2 ;;
        --parallel-strategy) OVERRIDE_PARALLEL_STRATEGY="${2:?Missing parallel-strategy value}"; shift 2 ;;
        --estimate-ate)     OVERRIDE_ESTIMATE_ATE="${2:?Missing estimate-ate value}"; shift 2 ;;
        --array-concurrency) OVERRIDE_ARRAY_CONCURRENCY="${2:?Missing array-concurrency value}"; shift 2 ;;
        -h|--help)         usage ;;
        *)                 echo "Error: unknown option '$1'"; usage ;;
    esac
done

# ============================================================================
# 2. Define Parameter Grid
# ============================================================================
# Core parameters (always varied unless subset=minimal)
N_TOTAL_VALUES=(5000 10000 20000)
K_VALUES=(3 4 5)
P_VALUES=(10 50 100 200)
CONFIG_VALUES=("C1" "C2" "C3")

# New design parameters (defaults match main.R defaults)
OUTCOME_TYPE_VALUES=("binary")
SHIFT_STRENGTH_VALUES=(1.0)
N_SIMS=500
USE_LAMBDA_CACHE="TRUE"
VERBOSE_EVERY="10"
PARALLEL_STRATEGY="outer_priority"
ESTIMATE_ATE="FALSE"
ARRAY_CONCURRENCY="16"

# DGP selection: "facec" (default) or "face" (FACE JASA 2023 Section 5.1)
DGP_TYPE="facec"
ATE_DEVIATION_VALUES=(0.0)
N_DEVIATED_SITES_VALUES=(0)

case $SUBSET in
    minimal)
        N_TOTAL_VALUES=(1000)
        K_VALUES=(3)
        P_VALUES=(10)
        CONFIG_VALUES=("C1")
        ;;
    core)
        ;;  # defaults above
    continuous)
        OUTCOME_TYPE_VALUES=("binary" "continuous")
        ;;
    shift)
        SHIFT_STRENGTH_VALUES=(0.5 1.0 2.0)
        ;;
    full)
        OUTCOME_TYPE_VALUES=("binary" "continuous")
        SHIFT_STRENGTH_VALUES=(0.5 1.0 2.0)
        P_VALUES=(10 20 50 100)
        ;;
    face)
        # FACE DGP (Han et al. JASA 2023, Section 5.1)
        # Supports both continuous and binary outcomes
        DGP_TYPE="face"
        N_TOTAL_VALUES=(20000)
        K_VALUES=(5)
        P_VALUES=(10 50 100 200)
        CONFIG_VALUES=("C1" "C2" "C3")
        OUTCOME_TYPE_VALUES=("binary")
        SHIFT_STRENGTH_VALUES=(1.0)
        ATE_DEVIATION_VALUES=(0.0 1.0 2.0 3.0)
        N_DEVIATED_SITES_VALUES=(0 2 4)
        ;;
    *)
        echo "Error: unknown subset '${SUBSET}'"
        echo "Available: minimal, core, continuous, shift, full, face"
        exit 1
        ;;
esac

# ============================================================================
# 2b. Apply CLI Overrides (after subset defaults)
# ============================================================================
[[ -n "$OVERRIDE_N_SIMS" ]]        && N_SIMS=$OVERRIDE_N_SIMS
[[ -n "$OVERRIDE_OUTCOME_TYPE" ]]  && read -ra OUTCOME_TYPE_VALUES <<< "$OVERRIDE_OUTCOME_TYPE"
[[ -n "$OVERRIDE_N_TOTAL" ]]       && read -ra N_TOTAL_VALUES <<< "$OVERRIDE_N_TOTAL"
[[ -n "$OVERRIDE_K" ]]             && read -ra K_VALUES <<< "$OVERRIDE_K"
[[ -n "$OVERRIDE_P" ]]             && read -ra P_VALUES <<< "$OVERRIDE_P"
[[ -n "$OVERRIDE_USE_LAMBDA_CACHE" ]] && USE_LAMBDA_CACHE="$OVERRIDE_USE_LAMBDA_CACHE"
[[ -n "$OVERRIDE_VERBOSE_EVERY" ]] && VERBOSE_EVERY="$OVERRIDE_VERBOSE_EVERY"
[[ -n "$OVERRIDE_PARALLEL_STRATEGY" ]] && PARALLEL_STRATEGY="$OVERRIDE_PARALLEL_STRATEGY"
[[ -n "$OVERRIDE_ESTIMATE_ATE" ]] && ESTIMATE_ATE="$OVERRIDE_ESTIMATE_ATE"
[[ -n "$OVERRIDE_ARRAY_CONCURRENCY" ]] && ARRAY_CONCURRENCY="$OVERRIDE_ARRAY_CONCURRENCY"

normalize_bool() {
    case "$1" in
        1|true|TRUE|True|t|T|yes|YES|Yes|y|Y) echo "TRUE" ;;
        0|false|FALSE|False|f|F|no|NO|No|n|N) echo "FALSE" ;;
        *) return 1 ;;
    esac
}

if ! USE_LAMBDA_CACHE=$(normalize_bool "$USE_LAMBDA_CACHE"); then
    echo "Error: --use-lambda-cache must be true/false, yes/no, or 1/0"
    exit 1
fi
if ! ESTIMATE_ATE=$(normalize_bool "$ESTIMATE_ATE"); then
    echo "Error: --estimate-ate must be true/false, yes/no, or 1/0"
    exit 1
fi

if [[ ! "$PARALLEL_STRATEGY" =~ ^(outer_priority|balanced|outer_only)$ ]]; then
    echo "Error: --parallel-strategy must be one of: outer_priority, balanced, outer_only"
    exit 1
fi
if ! [[ "$N_SIMS" =~ ^[0-9]+$ ]] || [[ "$N_SIMS" -lt 1 ]]; then
    echo "Error: --n-sims must be a positive integer"
    exit 1
fi
if ! [[ "$VERBOSE_EVERY" =~ ^[0-9]+$ ]] || [[ "$VERBOSE_EVERY" -lt 1 ]]; then
    echo "Error: --verbose-every must be a positive integer"
    exit 1
fi
if ! [[ "$ARRAY_CONCURRENCY" =~ ^[0-9]+$ ]] || [[ "$ARRAY_CONCURRENCY" -lt 1 ]]; then
    echo "Error: --array-concurrency must be a positive integer"
    exit 1
fi

validate_positive_int_values() {
    local name=$1
    shift
    local value
    for value in "$@"; do
        if ! [[ "$value" =~ ^[0-9]+$ ]] || [[ "$value" -lt 1 ]]; then
            echo "Error: ${name} values must be positive integers; got '${value}'"
            exit 1
        fi
    done
}

validate_nonnegative_int_values() {
    local name=$1
    shift
    local value
    for value in "$@"; do
        if ! [[ "$value" =~ ^[0-9]+$ ]]; then
            echo "Error: ${name} values must be non-negative integers; got '${value}'"
            exit 1
        fi
    done
}

validate_numeric_values() {
    local name=$1
    shift
    local value
    for value in "$@"; do
        if ! [[ "$value" =~ ^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$ ]]; then
            echo "Error: ${name} values must be numeric; got '${value}'"
            exit 1
        fi
    done
}

validate_choice_values() {
    local name=$1 choices_regex=$2
    shift 2
    local value
    for value in "$@"; do
        if ! [[ "$value" =~ $choices_regex ]]; then
            echo "Error: ${name} has invalid value '${value}'"
            exit 1
        fi
    done
}

validate_positive_int_values "--n-total" "${N_TOTAL_VALUES[@]}"
validate_positive_int_values "--K" "${K_VALUES[@]}"
validate_positive_int_values "--p" "${P_VALUES[@]}"
validate_choice_values "--outcome-type" "^(binary|continuous)$" "${OUTCOME_TYPE_VALUES[@]}"
validate_choice_values "config" "^(C1|C2|C3|C4)$" "${CONFIG_VALUES[@]}"
validate_numeric_values "shift_strength" "${SHIFT_STRENGTH_VALUES[@]}"
validate_numeric_values "ate_deviation" "${ATE_DEVIATION_VALUES[@]}"
validate_nonnegative_int_values "n_deviated_sites" "${N_DEVIATED_SITES_VALUES[@]}"

# ============================================================================
# 3. Helper: Build Setting ID
# ============================================================================
# Constructs a unique, human-readable identifier for each parameter combination.
# Includes active DGP parameters plus output-changing analysis flags.
# Format: n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_<alloc>_tf<transform>_ht<het>_ss<shift>_kf<folds>_ate<flag>
# Must match main.R and main.cmd

build_setting_id() {
    local n_total=$1 K=$2 p=$3 config=$4
    local outcome_type=$5 shift_strength=$6
    local heterogeneity_type=${7:-none}
    local estimand_type=${8:-superpopulation}
    local site_allocation=${9:-model}
    local transform_type=${10:-mild}
    local n_folds=${11:-10}
    local dgp_type=${12:-facec}
    local ate_deviation=${13:-0.0}
    local n_deviated_sites=${14:-0}
    local estimate_ate=${15:-FALSE}

    if [[ "$dgp_type" == "face" ]]; then
        echo "n${n_total}_K${K}_p${p}_${config}_${estimand_type}_${outcome_type}_face_dev${ate_deviation}_nd${n_deviated_sites}_kf${n_folds}_ate${estimate_ate}"
    else
        echo "n${n_total}_K${K}_p${p}_${config}_${estimand_type}_${outcome_type}_${site_allocation}_tf${transform_type}_ht${heterogeneity_type}_ss${shift_strength}_kf${n_folds}_ate${estimate_ate}"
    fi
}

# ============================================================================
# 4. Generate All Parameter Combinations
# ============================================================================
# Flatten the multi-dimensional grid into a list of argument strings.

COMBOS=()
for n in "${N_TOTAL_VALUES[@]}"; do
for K in "${K_VALUES[@]}"; do
for p in "${P_VALUES[@]}"; do
for cfg in "${CONFIG_VALUES[@]}"; do
for ot in "${OUTCOME_TYPE_VALUES[@]}"; do
for ss in "${SHIFT_STRENGTH_VALUES[@]}"; do
for ad in "${ATE_DEVIATION_VALUES[@]}"; do
for nd in "${N_DEVIATED_SITES_VALUES[@]}"; do
    COMBOS+=("${n} ${K} ${p} ${cfg} ${ot} ${ss} ${ad} ${nd}")
done; done; done; done; done; done; done; done

TOTAL=${#COMBOS[@]}

# ============================================================================
# 5. Print Summary and Confirm
# ============================================================================

mkdir -p log checkpoints results

echo "======================================================"
echo " FACE-C Simulation Job Submission"
echo "======================================================"
echo ""
echo " Subset:     ${SUBSET}"
echo " Total jobs: ${TOTAL}"
echo " Dry run:    ${DRY_RUN}"
echo ""
echo " Core parameters:"
echo "   n_total:           ${N_TOTAL_VALUES[*]}"
echo "   K:                 ${K_VALUES[*]}"
echo "   p:                 ${P_VALUES[*]}"
echo "   config:            ${CONFIG_VALUES[*]}"
echo ""
echo " Design parameters:"
echo "   outcome_type:      ${OUTCOME_TYPE_VALUES[*]}"
echo "   shift_strength:    ${SHIFT_STRENGTH_VALUES[*]}"
echo "   n_folds:           10 (fixed)"
echo "   n_sims:            ${N_SIMS}"
echo ""
echo " DGP:"
echo "   dgp_type:          ${DGP_TYPE}"
if [[ "$DGP_TYPE" == "face" ]]; then
echo "   ate_deviation:     ${ATE_DEVIATION_VALUES[*]}"
echo "   n_deviated_sites:  ${N_DEVIATED_SITES_VALUES[*]}"
else
echo "   allocation:        model (fixed)"
echo "   heterogeneity:     none (fixed)"
echo "   estimand:          superpopulation (fixed)"
echo "   transform:         mild (fixed)"
fi
echo ""
echo " Performance (set in main.cmd):"
echo "   NLAMBDA_INIT=100  lambda grid size for initial outcome CV"
echo "   NESTED_PARALLEL=1  sim-level + source-site parallel"
echo "   use_lambda_cache: ${USE_LAMBDA_CACHE}"
echo "   verbose_every:    ${VERBOSE_EVERY}"
echo "   parallel_strategy:${PARALLEL_STRATEGY}"
echo "   estimate_ate:     ${ESTIMATE_ATE}"
echo "   array_concurrency:${ARRAY_CONCURRENCY}"
echo ""
echo "======================================================"

if [[ "$DRY_RUN" == true ]]; then
    echo ""
    echo "[DRY RUN] Listing all ${TOTAL} jobs (nothing will be submitted):"
    echo ""
fi

# ============================================================================
# 6. Submit (or List) Jobs
# ============================================================================

if [[ "$DRY_RUN" == true ]]; then
    job_count=0
    for combo in "${COMBOS[@]}"; do
        read -r n_total K p config outcome_type shift_strength ate_deviation n_deviated_sites <<< "$combo"
        job_count=$((job_count + 1))
        setting_id=$(build_setting_id "$n_total" "$K" "$p" "$config" \
                                      "$outcome_type" "$shift_strength" \
                                      "none" "superpopulation" "model" "mild" "10" \
                                      "$DGP_TYPE" "$ate_deviation" "$n_deviated_sites" \
                                      "$ESTIMATE_ATE")
        printf "  [%3d/%d] %s\n" "$job_count" "$TOTAL" "$setting_id"
    done
else
    mkdir -p log/job_arrays
    combo_file="log/job_arrays/${SUBSET}_$(date +%Y%m%d_%H%M%S)_$$.tsv"

    for combo in "${COMBOS[@]}"; do
        read -r n_total K p config outcome_type shift_strength ate_deviation n_deviated_sites <<< "$combo"
        setting_id=$(build_setting_id "$n_total" "$K" "$p" "$config" \
                                      "$outcome_type" "$shift_strength" \
                                      "none" "superpopulation" "model" "mild" "10" \
                                      "$DGP_TYPE" "$ate_deviation" "$n_deviated_sites" \
                                      "$ESTIMATE_ATE")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
               "$n_total" "$K" "$p" "$config" "$outcome_type" "$shift_strength" \
               "$ate_deviation" "$n_deviated_sites" "$setting_id" >> "$combo_file"
    done

    submit_output=$(sbatch \
        --array="1-${TOTAL}%${ARRAY_CONCURRENCY}" \
        --export="COMBO_FILE=${combo_file},arg6=superpopulation,arg7=model,arg8=mild,arg10=none,arg12=10,arg13=${N_SIMS},arg14=${DGP_TYPE},arg17=${USE_LAMBDA_CACHE},arg18=${VERBOSE_EVERY},arg19=${PARALLEL_STRATEGY},arg20=${ESTIMATE_ATE}" \
        --job-name="FACE_${SUBSET}" \
        main.cmd)

    echo "Submitted job array: ${submit_output}"
    echo "Combo file: ${combo_file}"
    job_count=$TOTAL
fi

# ============================================================================
# 7. Final Summary
# ============================================================================

echo ""
echo "======================================================"
if [[ "$DRY_RUN" == true ]]; then
    echo " DRY RUN COMPLETE — would submit ${job_count} jobs"
else
    echo " DONE — submitted job array with ${job_count} tasks"
fi
echo "======================================================"
echo ""
echo " Output directories:"
echo "   Logs:        log/<array_job_id>_<array_task_id>.out"
echo "   Checkpoints: checkpoints/<config>/"
echo "   Results:     results/"
echo ""
echo " Useful commands:"
echo "   squeue -u \$USER                # monitor queue"
echo "   tail -f log/*.out               # watch live output"
echo "   scancel -u \$USER               # cancel all jobs"
echo "======================================================"
