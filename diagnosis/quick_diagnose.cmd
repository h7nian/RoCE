#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --job-name=FACE-diag
#SBATCH --output=diagnosis/log/%x_%j.out
#SBATCH --error=diagnosis/log/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --nice=5
#SBATCH --requeue
#SBATCH --signal=B:USR1@90

# ============================================================================
# Checkpoint / Requeue Support
# ----------------------------------------------------------------------------
# Each diagnostic step is identified by its name. After a step succeeds we
# append "step:<name>" to ${CHECKPOINT_FILE}. When SLURM requeues the job
# (same SLURM_JOB_ID), completed steps are skipped on the next attempt.
#
# Preemption flow (matches main.cmd):
#   1. SLURM sends SIGUSR1 90s before termination on the preempt partition.
#   2. handle_preemption() touches the saved-signal file (no R kernel to
#      flush — checkpoint is the text file itself) and calls scontrol
#      requeue. SLURM re-launches this script under the same JOB_ID.
#   3. On restart, run_step() skips any step whose name is in the
#      checkpoint file.
#
# Use DIAG_RESET_CHECKPOINT=TRUE to start fresh (used by --reset-checkpoint
# in quick_diagnose.sh). The reset only triggers on the first attempt; on
# requeue we keep the checkpoint we just wrote.
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"

mkdir -p diagnosis/log diagnosis/results diagnosis/tmp

JOB_ID="${SLURM_JOB_ID:-manual}"
CHECKPOINT_FILE="diagnosis/tmp/checkpoint_diag_${JOB_ID}.txt"
RESET_MARKER="diagnosis/tmp/checkpoint_diag_${JOB_ID}.reset_processed"
SAVED_SIGNAL="diagnosis/tmp/.checkpoint_saved_${JOB_ID}"
PREEMPT_SIGNAL="diagnosis/tmp/.preempt_signal_${JOB_ID}"
REPORT="diagnosis/results/quick_diagnosis_${JOB_ID}.txt"
PREEMPTED=0

if [[ "${DIAG_RESET_CHECKPOINT:-FALSE}" == "TRUE" && ! -f "${RESET_MARKER}" ]]; then
    rm -f "${CHECKPOINT_FILE}" "${SAVED_SIGNAL}" "${PREEMPT_SIGNAL}"
    touch "${RESET_MARKER}"
fi
touch "${CHECKPOINT_FILE}"

handle_preemption() {
    echo "$(date): Received preemption signal (SIGUSR1). Marking state..."
    touch "${PREEMPT_SIGNAL}"
    PREEMPTED=1
    # No R kernel running between steps — the checkpoint file is the state.
    touch "${SAVED_SIGNAL}"
    echo "$(date): Requesting requeue (JOB_ID=${JOB_ID})..."
    if ! scontrol requeue "${JOB_ID}"; then
        echo "$(date): WARNING: scontrol requeue failed; exiting without requeue." >&2
    fi
    exit 0
}
trap 'handle_preemption' USR1

mark_done() {
    local step_id=$1
    grep -Fxq "${step_id}" "${CHECKPOINT_FILE}" 2>/dev/null || \
        echo "${step_id}" >> "${CHECKPOINT_FILE}"
}

is_done() {
    local step_id=$1
    grep -Fxq "${step_id}" "${CHECKPOINT_FILE}" 2>/dev/null
}

# Append-tee so requeued runs continue the same report file.
exec > >(tee -a "${REPORT}") 2>&1

ATTEMPT_TS="$(date +%Y%m%d_%H%M%S)"
echo "======================================================"
echo " FACE-HD Quick Diagnosis"
echo "======================================================"
echo " Attempt start:  $(date)"
echo " Attempt tag:    ${ATTEMPT_TS}"
echo " Host:           $(hostname)"
echo " Job ID:         ${JOB_ID}"
echo " Submit dir:     ${SLURM_SUBMIT_DIR:-$(pwd)}"
echo " Report file:    ${REPORT}"
echo " Checkpoint:     ${CHECKPOINT_FILE}"
echo " Install check:  ${DIAG_INSTALL_CHECK:-FALSE}"
echo " Reset on start: ${DIAG_RESET_CHECKPOINT:-FALSE}"
if [[ -s "${CHECKPOINT_FILE}" ]]; then
    echo " Completed steps from prior attempt(s):"
    sed 's/^/   - /' "${CHECKPOINT_FILE}"
fi
echo "======================================================"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

# ----------------------------------------------------------------------------
# run_step: idempotent, checkpoint-aware step runner.
# - Skips if the step is already in the checkpoint file.
# - Runs the command in the FOREGROUND so heredoc stdin (used by
#   `Rscript --vanilla - <<'RS' ... RS`) is correctly inherited. Backgrounding
#   with `& wait` would have bash detach stdin from the child in some
#   configurations, causing the Rscript to read EOF immediately and exit
#   silently (which is what happened in job 8976393's first run).
# - SIGUSR1 is still serviced between steps; diagnostic steps are short
#   relative to the 90s preemption grace period, so a 1-step latency is fine.
# - Uses `set +e` around the call so we can capture the exit code, then
#   re-raise it ourselves with full step context.
# ----------------------------------------------------------------------------
run_step() {
    local name=$1
    shift
    local step_id="step:${name}"
    if is_done "${step_id}"; then
        echo
        echo "===== ${name} ===== [skipped: completed in prior attempt]"
        return 0
    fi
    if [[ ${PREEMPTED} -eq 1 ]]; then
        echo
        echo "===== ${name} ===== [skipped: preemption signal received]"
        return 0
    fi
    echo
    echo "===== ${name} ====="
    set +e
    "$@"
    local rc=$?
    set -e
    if [[ ${PREEMPTED} -eq 1 ]]; then
        # handle_preemption will have called scontrol requeue and exited.
        return 0
    fi
    if [[ ${rc} -ne 0 ]]; then
        echo "[FAIL] step '${name}' exited with ${rc}" >&2
        exit ${rc}
    fi
    mark_done "${step_id}"
}

run_step "Repository status" bash -lc '
    git rev-parse --show-toplevel
    git status --short
'

run_step "Shell syntax" bash -lc '
    files=(main.sh main.cmd realdata.sh realdata.cmd test.sh test.cmd diagnosis/*.sh diagnosis/*.cmd)
    for f in "${files[@]}"; do
        [[ -f "$f" ]] || continue
        echo "[bash -n] $f"
        bash -n "$f"
    done
'

run_step "R version and library paths" Rscript --vanilla -e '
    cat(R.version.string, "\n")
    cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n", sep = "")
    cat("libPaths:\n")
    cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")
'

run_step "R parse all project R files" Rscript --vanilla - <<'RS'
files <- unique(c(
  "main.R",
  "realdata.R",
  list.files("R", pattern = "\\.R$", full.names = TRUE),
  list.files("scripts", pattern = "\\.R$", full.names = TRUE),
  list.files("tests", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
))
files <- files[file.exists(files)]
bad <- character()
for (f in files) {
  cat("[parse]", f, "\n")
  ok <- tryCatch({
    parse(file = f)
    TRUE
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    FALSE
  })
  if (!ok) bad <- c(bad, f)
}
if (length(bad)) {
  stop("Parse failed for: ", paste(bad, collapse = ", "))
}
cat("Parsed ", length(files), " R files successfully.\n", sep = "")
RS

run_step "Static argument and fail-fast checks" Rscript --vanilla - <<'RS'
read_file <- function(path) paste(readLines(path, warn = FALSE), collapse = "\n")
must_have <- function(path, pattern, label) {
  txt <- read_file(path)
  if (!grepl(pattern, txt, perl = TRUE)) {
    stop(label, " missing in ", path)
  }
  cat("[OK]", label, "\n")
}
must_not_have <- function(path, pattern, label) {
  txt <- read_file(path)
  if (grepl(pattern, txt, perl = TRUE)) {
    stop(label, " still present in ", path)
  }
  cat("[OK]", label, "\n")
}

must_have("main.R", "estimate_ate_arg\\s*<-", "main.R explicit estimate_ate parser")
must_have("main.R", "estimate_ate\\s*=\\s*estimate_ate_arg", "main.R passes explicit estimate_ate")
must_not_have("main.R", "estimate_ate\\s*=\\s*\\(heterogeneity_type_arg\\s*!=\\s*\"none\"\\)", "implicit heterogeneity-to-ATE coupling")
must_have("main.sh", "arg20=\"?\\$\\{ESTIMATE_ATE\\}\"?", "main.sh exports arg20 estimate_ate")
must_have("main.cmd", "ESTIMATE_ATE=\"?\\$\\{arg20:-FALSE\\}\"?", "main.cmd receives arg20 estimate_ate")
must_have("R/validation.R", "warn_ignored", "central ignored-argument warning control")
must_not_have("R/simulation.R", "tryCatch\\s*\\(\\s*\\{\\s*if \\(n_cores_outer", "setting-level tryCatch swallowing")
must_not_have("R/comparison_methods.R", "safe_run\\s*<-", "comparison safe_run swallowing")
must_not_have("R/model_fitting.R", "Falling back to target-only", "optimizer target-only fallback")
must_not_have("R/cross_fitting_algorithms.R", "zero initial params", "zero-initial nuisance fallback")
must_not_have("R/checkpoint.R", "starting fresh", "corrupt checkpoint silent restart")
cat("Static argument and fail-fast checks passed.\n")
RS

run_step "NAMESPACE export existence check" Rscript --vanilla - <<'RS'
# Verify that every NAMESPACE-exported name is actually defined in R/.
# We must accept BOTH function definitions (`name <- function(...)`)
# and top-level constant assignments (`name <- value`), because R/constants.R
# exports symbols like FAMILY_GAUSSIAN, Z_ALPHA_05, etc.
# Anchoring on `(?m)^\s*<name>\s*<-` avoids false positives from substrings
# inside function bodies, comments, or string literals.
ns <- readLines("NAMESPACE", warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
exports <- exports[nzchar(exports)]
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
r_text <- paste(vapply(r_files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1)), collapse = "\n")

is_defined <- function(name) {
  escaped <- gsub("([.|()\\^{}+$*?\\[\\]\\\\])", "\\\\\\1", name)
  # multiline mode (?m): ^ matches start of any line within r_text
  grepl(paste0("(?m)^\\s*", escaped, "\\s*<-"), r_text, perl = TRUE)
}

missing <- exports[!vapply(exports, is_defined, logical(1))]
if (length(missing)) {
  stop("NAMESPACE exports missing top-level definitions: ",
       paste(missing, collapse = ", "))
}
cat("Checked ", length(exports), " exported names (functions + constants).\n", sep = "")
RS

run_step "Formula-alignment static scan" bash -lc '
    if command -v rg >/dev/null 2>&1; then
        GREP_CMD=(rg -n)
        RECURSIVE_FLAG=()
    else
        GREP_CMD=(grep -En)
        RECURSIVE_FLAG=(-r)
    fi
    formula_targets=()
    for f in R/model_fitting.R R/cross_fitting_aggregation.R R/cross_fitting_algorithms.R R/estimators_oracle.R docs/main.tex; do
        [[ -f "$f" ]] && formula_targets+=("$f")
    done
    echo "[eq:final_opt references]"
    if (( ${#formula_targets[@]} > 0 )); then
        "${GREP_CMD[@]}" "eq:final_opt|clip_weights|optimize_weights|C_cross|lambda_rule" "${formula_targets[@]}" | head -n 80 || true
    else
        echo "  (no target files present)"
    fi
    echo
    echo "[remaining tryCatch contexts]"
    "${GREP_CMD[@]}" "${RECURSIVE_FLAG[@]}" "tryCatch" R main.R 2>/dev/null || true
    echo
    echo "[remaining warning contexts]"
    "${GREP_CMD[@]}" "${RECURSIVE_FLAG[@]}" "warning\\(" R main.R 2>/dev/null || true
'

run_step "Package load from source" Rscript --vanilla - <<'RS'
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is not available; cannot load current source tree.")
}
devtools::load_all(".", quiet = FALSE)
required <- c(
  "run_crossfit",
  "calculate_crossfit_aggregation",
  "run_single_simulation",
  "run_simulation_study",
  "validate_simulation_params",
  "generate_simulation_data"
)
missing <- required[!vapply(required, exists, logical(1), mode = "function")]
if (length(missing)) {
  stop("Missing loaded functions: ", paste(missing, collapse = ", "))
}
cat("Source package load succeeded.\n")
RS

run_step "Run testthat suite" Rscript --vanilla - <<'RS'
if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is not available; cannot run testthat suite.")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is not available; cannot run testthat suite.")
}
# devtools::test() re-runs load_all() and then dispatches to testthat.
# stop_on_failure = TRUE makes failures / errors propagate as an R error,
# which run_step turns into a non-zero exit code that fails the Slurm step.
# We deliberately do NOT use the "summary" reporter alone because we want
# every failed test to be visible in the Slurm log, not just totals.
devtools::test(".", reporter = c("progress", "fail"), stop_on_failure = TRUE)
cat("All testthat tests passed.\n")
RS

if [[ "${DIAG_INSTALL_CHECK:-FALSE}" == "TRUE" ]]; then
    run_step "Temporary R CMD INSTALL" bash -lc '
        install_lib="diagnosis/tmp/Rlib_${SLURM_JOB_ID:-manual}"
        rm -rf "${install_lib}"
        mkdir -p "${install_lib}"
        R CMD INSTALL --library="${install_lib}" .
    '
fi

echo
echo "======================================================"
echo " Quick diagnosis completed successfully"
echo " End: $(date)"
echo " Report: ${REPORT}"
echo "======================================================"

# Successful exit — wipe checkpoint artifacts so a follow-up submission starts clean.
rm -f "${CHECKPOINT_FILE}" "${RESET_MARKER}" "${SAVED_SIGNAL}" "${PREEMPT_SIGNAL}"
echo "Checkpoint files removed after successful completion."
