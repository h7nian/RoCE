# main.R - RoCE Monte Carlo Simulation Study
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# Primary estimand for the bounded default: target average treatment effect,
#   tau_t = E_t[Y(1)] - E_t[Y(0)].
# Both arms share the same data and folds; TATE inference retains their covariance.
# Potential-outcome means remain available as secondary outputs.
# =============================================================================

args <- commandArgs(TRUE)

print("Command line arguments:")
print(args)

# ============================================================================
# Command Line Arguments
# ============================================================================
# arg1: n_total (sample size)
# arg2: K (number of source sites)
# arg3: p (number of covariates)
# arg4: config ("C1", "C2", "C3", "C4")
# arg5: job_id (SLURM job ID, optional)
# arg6: estimand_type ("sample" or "superpopulation", optional, default "superpopulation")
# arg7: site_allocation ("model", "uniform", "balanced", etc., default "model")
# arg8: transform_type ("strong", "mild", "none", default "mild")
# arg9: outcome_type ("binary", "continuous", default "binary")
# arg10: heterogeneity_type ("none", "mild", "strong", "partial", default "none")
# arg11: shift_strength (numeric, default 0.5)
# arg12: n_folds (integer, cross-fitting folds, default 10)
#        The bounded default uses three fold roles on these ten original folds.
# arg13: n_sims (integer, number of Monte Carlo sims, default 500)
# arg14: dgp_type ("bounded", "roce", or "face"; default "bounded")
# arg15: ate_deviation (numeric ATE deviation for FACE paper non-informative sites, default 0.0)
# arg16: n_deviated_sites (integer, deviated source sites for FACE paper DGP, default 0)
# arg17: use_lambda_cache (TRUE/FALSE, default TRUE)
# arg18: verbose_every (integer, sequential logging interval, default 10)
# arg19: parallel_strategy (outer_priority|balanced|outer_only, default outer_priority)
# arg20: estimate_ate (TRUE/FALSE, default TRUE for bounded, otherwise FALSE)
# arg21: nuisance_cv_certificate (TRUE/FALSE, default FALSE)
# ============================================================================

.arg_display <- function(value) {
  if (is.null(value) || length(value) == 0L) return("<missing>")
  paste(value, collapse = ",")
}

.abort_arg <- function(name, value, reason) {
  stop(sprintf("Invalid %s='%s': %s", name, .arg_display(value), reason),
       call. = FALSE)
}

.raw_arg <- function(args, idx, name, default = NULL, required = FALSE) {
  if (length(args) >= idx) return(args[[idx]])
  if (required) {
    .abort_arg(name, NULL, sprintf("argument %d is required.", idx))
  }
  default
}

.parse_numeric_arg <- function(args, idx, name, default = NULL,
                               required = FALSE, integer = FALSE,
                               min_value = NULL) {
  raw <- .raw_arg(args, idx, name, default = default, required = required)
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    .abort_arg(name, raw, "must be a finite numeric scalar.")
  }
  if (integer) {
    if (abs(value - round(value)) > sqrt(.Machine$double.eps)) {
      .abort_arg(name, raw, "must be an integer scalar.")
    }
    value <- as.integer(round(value))
  }
  if (!is.null(min_value) && value < min_value) {
    .abort_arg(name, raw, sprintf("must be >= %s.", min_value))
  }
  value
}

.parse_bool_value <- function(raw, name) {
  if (is.logical(raw) && length(raw) == 1L && !is.na(raw)) {
    return(isTRUE(raw))
  }
  raw_chr <- tolower(trimws(as.character(raw)))
  if (raw_chr %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (raw_chr %in% c("0", "false", "f", "no", "n")) return(FALSE)
  .abort_arg(name, raw, "must be TRUE/FALSE, yes/no, or 1/0.")
}

.parse_bool_arg <- function(args, idx, name, default = FALSE) {
  raw <- .raw_arg(args, idx, name, default = default)
  .parse_bool_value(raw, name)
}

.parse_choice_arg <- function(args, idx, name, choices, default) {
  raw <- .raw_arg(args, idx, name, default = default)
  raw_chr <- tolower(trimws(as.character(raw)))
  choices_chr <- tolower(choices)
  matched <- match(raw_chr, choices_chr)
  if (is.na(matched)) {
    .abort_arg(name, raw, sprintf("must be one of: %s.",
                                  paste(choices, collapse = ", ")))
  }
  choices[[matched]]
}

.parse_env_integer <- function(name, default, min_value = NULL) {
  raw <- Sys.getenv(name, unset = NA_character_)
  if (is.na(raw) || identical(raw, "")) raw <- default
  value <- suppressWarnings(as.numeric(raw))
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      abs(value - round(value)) > sqrt(.Machine$double.eps)) {
    .abort_arg(name, raw, "must be an integer scalar.")
  }
  value <- as.integer(round(value))
  if (!is.null(min_value) && value < min_value) {
    .abort_arg(name, raw, sprintf("must be >= %s.", min_value))
  }
  value
}

n_total_vec <- .parse_numeric_arg(args, 1, "n_total", required = TRUE,
                                  integer = TRUE, min_value = 1L)
K_vec <- .parse_numeric_arg(args, 2, "K", required = TRUE,
                            integer = TRUE, min_value = 1L)
p_vec <- .parse_numeric_arg(args, 3, "p", required = TRUE,
                            integer = TRUE, min_value = 1L)
config_arg <- if (length(args) >= 4) args[4] else "C1"
job_id <- if (length(args) >= 5) args[5] else format(Sys.time(), "%Y%m%d_%H%M%S")
estimand_type_arg <- if (length(args) >= 6) args[6] else "superpopulation"
site_allocation_arg <- if (length(args) >= 7) args[7] else "model"
transform_type_arg <- if (length(args) >= 8) args[8] else "mild"
outcome_type_arg <- if (length(args) >= 9) args[9] else "binary"
heterogeneity_type_arg <- if (length(args) >= 10) args[10] else "none"
shift_strength_arg <- .parse_numeric_arg(args, 11, "shift_strength", default = 0.5)
n_folds_arg <- .parse_numeric_arg(args, 12, "n_folds", default = 10L,
                                  integer = TRUE, min_value = 1L)
n_sims_arg <- .parse_numeric_arg(args, 13, "n_sims", default = 500L,
                                 integer = TRUE, min_value = 1L)
dgp_type_arg <- .parse_choice_arg(args, 14, "dgp_type",
                                  choices = c("bounded", "roce", "face"),
                                  default = "bounded")
ate_deviation_arg <- .parse_numeric_arg(args, 15, "ate_deviation", default = 0.0)
n_deviated_sites_arg <- .parse_numeric_arg(args, 16, "n_deviated_sites",
                                           default = 0L, integer = TRUE)
use_lambda_cache_arg <- .parse_bool_arg(args, 17, "use_lambda_cache", default = TRUE)
verbose_every_arg <- .parse_numeric_arg(args, 18, "verbose_every", default = 10L,
                                        integer = TRUE, min_value = 1L)
parallel_strategy_arg <- .parse_choice_arg(
  args, 19, "parallel_strategy",
  choices = c("outer_priority", "balanced", "outer_only"),
  default = "outer_priority"
)
estimate_ate_arg <- .parse_bool_arg(args, 20, "estimate_ate", default = dgp_type_arg == "bounded")
nuisance_cv_certificate_arg <- .parse_bool_arg(args, 21, "nuisance_cv_certificate", default = FALSE)

# ============================================================================
# Performance Configuration
# ============================================================================
# NLAMBDA_INIT: Number of lambda values for initial outcome model CV (glmnet).
# Lower values (e.g. 20) speed up fitting; higher values (default 100) give finer tuning.
# Set via environment variable NLAMBDA_INIT=20 (integer)
NLAMBDA_INIT <- .parse_env_integer("NLAMBDA_INIT", default = "100",
                                   min_value = 1L)
if (NLAMBDA_INIT != 100L) {
  cat(sprintf("NLAMBDA_INIT = %d (default 100)\n", NLAMBDA_INIT))
}

# NESTED_PARALLEL: When TRUE, enables both simulation-level and source-site parallelism
# - Outer layer: parLapply (PSOCK cluster) for simulations
# - Inner layer: mclapply (fork) for source sites within each simulation
# Total cores = outer_cores × inner_cores (auto-balanced to not exceed available)
# Set via environment variable NESTED_PARALLEL=1 or NESTED_PARALLEL=TRUE
NESTED_PARALLEL <- .parse_bool_value(Sys.getenv("NESTED_PARALLEL", "0"),
                                     "NESTED_PARALLEL")
if (NESTED_PARALLEL) {
  cat("[PERF] NESTED PARALLEL ENABLED: simulation + source-site parallelism\n")
}

# RoCE implementation status:
# - Canonical entry point: run_crossfit(..., communication_mode = "two_round"|"one_round")
# - Legacy non-crossfitting algorithm labels (one_round, two_round) are not used
#   in this script; cross-fitting variants are the supported workflow.

# Prefer an explicitly selected checked library, or the installed project default.
project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (!nzchar(project_library)) {
  candidate_library <- file.path("/scratch.global", Sys.getenv("USER"), "FACE-HD", "Rlib_default")
  if (dir.exists(file.path(candidate_library, "RoCE"))) project_library <- candidate_library
}
if (nzchar(project_library)) {
  if (!dir.exists(file.path(project_library, "RoCE"))) stop("ROCE_PROJECT_LIB does not contain RoCE.")
  .libPaths(c(project_library, .libPaths()))
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
}

# Load the RoCE package (all R/ code + compiled C++)
# Use the checked installed library by default. Source-tree loading is an
# explicit development option and can compile native code in the checkout.
use_source_tree <- .parse_bool_value(Sys.getenv("ROCE_MAIN_USE_SOURCE", "FALSE"),
                                     "ROCE_MAIN_USE_SOURCE")
if (use_source_tree) {
  if (!requireNamespace("devtools", quietly = FALSE)) {
    stop("ROCE_MAIN_USE_SOURCE=TRUE but devtools is not available.",
         call. = FALSE)
  }
  cat("ROCE_MAIN_USE_SOURCE=TRUE; loading source tree via devtools::load_all()...\n")
  devtools::load_all(".")
} else if (requireNamespace("RoCE", quietly = FALSE)) {
  library(RoCE)
} else if (requireNamespace("devtools", quietly = FALSE)) {
  cat("RoCE not installed; loading via devtools::load_all()...\n")
  devtools::load_all(".")
} else {
  stop("Neither installed RoCE package nor devtools found. ",
       "Install the package with: R CMD INSTALL .")
}

# Validate all simulation parameters via centralized function (R/validation.R)
validate_simulation_params(
  estimand_type    = estimand_type_arg,
  site_allocation  = site_allocation_arg,
  transform_type   = transform_type_arg,
  outcome_type     = outcome_type_arg,
  heterogeneity_type = heterogeneity_type_arg,
  shift_strength   = shift_strength_arg,
  n_folds          = n_folds_arg,
  n_sims           = n_sims_arg,
  n_total          = n_total_vec,
  K                = K_vec,
  p                = p_vec,
  config           = config_arg,
  dgp_type         = dgp_type_arg,
  ate_deviation    = ate_deviation_arg,
  n_deviated_sites = n_deviated_sites_arg,
  warn_ignored     = FALSE
)

# Build setting identifier from active DGP/output parameters.
# For the FACE paper DGP, the identifier also includes ate_deviation/n_deviated_sites.
# Format (roce): n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_<alloc>_tf<transform>_ht<het>_ss<shift>_kf<folds>_ate<flag>
# Format (face):  n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_face_dev<dev>_nd<nd>_kf<folds>_ate<flag>
# Must match main.sh::build_setting_id() and main.cmd
estimate_ate_label <- if (estimate_ate_arg) "TRUE" else "FALSE"
if (dgp_type_arg %in% c("face", "bounded")) {
  setting_id <- sprintf(
    "n%d_K%d_p%d_%s_%s_%s_%s_dev%.1f_nd%d_kf%d_ate%s",
    n_total_vec, K_vec, p_vec, config_arg, estimand_type_arg,
    outcome_type_arg, dgp_type_arg, ate_deviation_arg, n_deviated_sites_arg, n_folds_arg,
    estimate_ate_label
  )
} else {
  setting_id <- sprintf(
    "n%d_K%d_p%d_%s_%s_%s_%s_tf%s_ht%s_ss%.1f_kf%d_ate%s",
    n_total_vec, K_vec, p_vec, config_arg,
    estimand_type_arg, outcome_type_arg, site_allocation_arg,
    transform_type_arg, heterogeneity_type_arg,
    shift_strength_arg, n_folds_arg, estimate_ate_label
  )
}

cat(sprintf("\n============================================\n"))
cat(sprintf("RoCE Simulation: %s\n", setting_id))
cat(sprintf("  n_total: %d\n", n_total_vec))
cat(sprintf("  K: %d\n", K_vec))
cat(sprintf("  p: %d\n", p_vec))
cat(sprintf("  config: %s\n", config_arg))
cat(sprintf("  estimand_type: %s\n", estimand_type_arg))
cat(sprintf("  dgp_type: %s\n", dgp_type_arg))
if (dgp_type_arg %in% c("face", "bounded")) {
  cat(sprintf("  outcome_type: %s\n", outcome_type_arg))
  cat(sprintf("  ate_deviation: %.1f\n", ate_deviation_arg))
  cat(sprintf("  n_deviated_sites: %d\n", n_deviated_sites_arg))
} else {
  cat(sprintf("  site_allocation: %s\n", site_allocation_arg))
  cat(sprintf("  transform_type: %s\n", transform_type_arg))
  cat(sprintf("  outcome_type: %s\n", outcome_type_arg))
  cat(sprintf("  heterogeneity_type: %s\n", heterogeneity_type_arg))
  cat(sprintf("  shift_strength: %.1f\n", shift_strength_arg))
}
cat(sprintf("  n_folds: %d\n", n_folds_arg))
cat(sprintf("  n_sims: %d\n", n_sims_arg))
cat(sprintf("  use_lambda_cache: %s\n", if (use_lambda_cache_arg) "TRUE" else "FALSE"))
cat(sprintf("  verbose_every: %d\n", verbose_every_arg))
cat(sprintf("  parallel_strategy: %s\n", parallel_strategy_arg))
cat(sprintf("  estimate_ate: %s\n", estimate_ate_label))
cat(sprintf("  job_id: %s\n", job_id))
cat(sprintf("============================================\n\n"))

# ============================================================================
# Checkpoint/Restart Configuration
# ============================================================================
# checkpoint functions (init_checkpoint_config, save_checkpoint, load_checkpoint,
# check_preempt_signal, signal_checkpoint_saved, cleanup_checkpoint) are
# provided by the RoCE package (R/checkpoint.R).
run_root <- Sys.getenv("ROCE_RUN_ROOT", file.path("/scratch.global", Sys.getenv("USER"), "FACE-HD", "runs"))
CHECKPOINT_DIR <- Sys.getenv("CHECKPOINT_DIR", file.path(run_root, "checkpoints"))
ckpt_config <- init_checkpoint_config(CHECKPOINT_DIR, setting_id, job_id)

# Set up parallel processing (already imported by RoCE, but needed for main.R scope)
library(parallel)
library(doParallel)

# ============================================================================
# SIMULATION FUNCTIONS
# ============================================================================
# run_single_simulation(), summarize_results(), and run_simulation_study() are
# defined in R/simulation.R and are part of the RoCE package.
# ============================================================================

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Output directory with date stamp
OUTPUT_DIR <- Sys.getenv("RESULTS_DIR", file.path(run_root, "results"))
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
}
if (!dir.exists(OUTPUT_DIR)) {
  stop(sprintf("Failed to create output directory '%s'.", OUTPUT_DIR),
       call. = FALSE)
}

# Run simulation for this specific setting
# Note: Each job runs ONE config + heterogeneity_type combination
cat(sprintf("\nRunning simulation for: %s\n", setting_id))

# Auto-detect available cores from SLURM or system
# SLURM sets SLURM_CPUS_PER_TASK, otherwise use detectCores() - 1
available_cores <- .parse_env_integer("SLURM_CPUS_PER_TASK", default = "0",
                                      min_value = 0L)
if (available_cores == 0) {
  available_cores <- max(1, parallel::detectCores() - 1)
}

# Guard against oversubscription in non-SLURM environments.
# By default, cap cores to 8 unless ROCE_MAX_CORES is explicitly set.
# Also never use more cores than n_sims to avoid idle workers.
roce_max_cores <- .parse_env_integer("ROCE_MAX_CORES", default = "8",
                                      min_value = 1L)

available_cores <- min(available_cores, n_sims_arg, roce_max_cores)
available_cores <- max(1, as.integer(available_cores))
cat(sprintf("Available cores (after cap): %d\n", available_cores))

sim_args <- list(
  n_sims = n_sims_arg,
  n_total_vec = n_total_vec,
  K_vec = K_vec,
  p_vec = p_vec,
  configs = config_arg,
  n_cores = available_cores,
  checkpoint_config = ckpt_config,
  nlambda_init = NLAMBDA_INIT,
  nested_parallel = NESTED_PARALLEL,
  estimand_type = estimand_type_arg,
  site_allocation = site_allocation_arg,
  transform_type = transform_type_arg,
  outcome_type = outcome_type_arg,
  heterogeneity_type = heterogeneity_type_arg,
  shift_strength = shift_strength_arg,
  n_folds = n_folds_arg,
  use_lambda_cache = use_lambda_cache_arg,
  verbose_every = verbose_every_arg,
  parallel_strategy = parallel_strategy_arg,
  estimate_ate = estimate_ate_arg,
  dgp_type = dgp_type_arg,
  ate_deviation = ate_deviation_arg,
  n_deviated_sites = n_deviated_sites_arg
)

if (dgp_type_arg == "bounded") {
  sim_args$methods <- c("one_round_crossfit", "target_only")
  sim_args$target_nuisance_method <- "hou_calibrated"
  sim_args$source_validation_method <- "calibrated"
  sim_args$calibration_layout <- "compact"
  sim_args$nuisance_solver <- "proximal_newton"
  sim_args$nuisance_tol <- 1e-10
  sim_args$M_tau <- sim_args$M_tau_inference <- 12
  sim_args$calibration_control <- list(recipe = "score_derivative",
    target_propensity_initialization = "calibrated", target_radius = 12)
  if (estimate_ate_arg) sim_args$additional_aggregation_modes <- c("separate_arms", "joint_tate")
}
if ("nuisance_cv_certificate" %in% names(formals(run_simulation_study))) {
  sim_args$nuisance_cv_certificate <- nuisance_cv_certificate_arg
} else if (nuisance_cv_certificate_arg) {
  stop("The selected RoCE library does not support nuisance_cv_certificate; select a checked updated library.")
}
results <- do.call(run_simulation_study, sim_args)

# Summarize results
summary_df <- summarize_results(results)

if (nrow(summary_df) > 0) {
  print(summary_df)
} else {
  cat("No summary statistics available (empty results)\n")
}

# Save results with setting-specific filename
save_path <- file.path(OUTPUT_DIR, sprintf("%s.RData", setting_id))
save(results, summary_df, file = save_path)
cat(sprintf("Results saved to: %s\n", save_path))

# Save results to CSV
csv_path <- file.path(OUTPUT_DIR, sprintf("%s.csv", setting_id))
write.csv(results, file = csv_path, row.names = FALSE)
cat(sprintf("CSV saved to: %s\n", csv_path))

# Save summary to CSV
summary_path <- file.path(OUTPUT_DIR, sprintf("%s_summary.csv", setting_id))
write.csv(summary_df, file = summary_path, row.names = FALSE)
cat(sprintf("Summary saved to: %s\n", summary_path))

# Clean up checkpoint after successful completion
cleanup_checkpoint(ckpt_config$file, ckpt_config$preempt_signal, ckpt_config$saved_signal)
  
cat(sprintf("\n============================================\n"))
cat(sprintf("Simulation completed: %s\n", setting_id))
cat(sprintf("============================================\n"))
