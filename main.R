# main.R - FACE-C Monte Carlo Simulation Study
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# Primary estimand: POTENTIAL OUTCOME MEAN at the target site:
#   μ¹_t = E_t[Y(1)]  (expected outcome under treatment)
#
# When estimate_ate=TRUE, the code also computes ATE = E_t[Y(1)] - E_t[Y(0)]
# by rerunning each method with A_val=0 for the mu0 arm.
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
# arg11: shift_strength (numeric, default 1.0)
# arg12: n_folds (integer, cross-fitting folds, default 10)
# arg13: n_sims (integer, number of Monte Carlo sims, default 500)
# arg14: dgp_type ("facec" or "face", default "facec")
# arg15: ate_deviation (numeric ATE deviation for FACE paper non-informative sites, default 0.0)
# arg16: n_deviated_sites (integer, deviated source sites for FACE paper DGP, default 0)
# arg17: use_lambda_cache (TRUE/FALSE, default TRUE)
# arg18: verbose_every (integer, sequential logging interval, default 10)
# arg19: parallel_strategy (outer_priority|balanced|outer_only, default outer_priority)
# ============================================================================
n_total_vec <- as.numeric(args[1])
K_vec <- as.numeric(args[2])
p_vec <- as.numeric(args[3])
config_arg <- if (length(args) >= 4) args[4] else "C1"
job_id <- if (length(args) >= 5) args[5] else format(Sys.time(), "%Y%m%d_%H%M%S")
estimand_type_arg <- if (length(args) >= 6) args[6] else "superpopulation"
site_allocation_arg <- if (length(args) >= 7) args[7] else "model"
transform_type_arg <- if (length(args) >= 8) args[8] else "mild"
outcome_type_arg <- if (length(args) >= 9) args[9] else "binary"
heterogeneity_type_arg <- if (length(args) >= 10) args[10] else "none"
shift_strength_arg <- if (length(args) >= 11) as.numeric(args[11]) else 1.0
n_folds_arg <- if (length(args) >= 12) as.integer(args[12]) else 10L
n_sims_arg <- if (length(args) >= 13) as.integer(args[13]) else 500L
dgp_type_arg <- if (length(args) >= 14) args[14] else "facec"
ate_deviation_arg <- if (length(args) >= 15) as.numeric(args[15]) else 0.0
n_deviated_sites_arg <- if (length(args) >= 16) as.integer(args[16]) else 0L
use_lambda_cache_arg <- if (length(args) >= 17) {
  tolower(args[17]) %in% c("1", "true", "t", "yes", "y")
} else {
  TRUE
}
verbose_every_arg <- if (length(args) >= 18) as.integer(args[18]) else 10L
if (is.na(verbose_every_arg) || verbose_every_arg < 1L) {
  warning(sprintf("Invalid verbose_every='%s'; falling back to 10.",
                  if (length(args) >= 18) args[18] else "NA"))
  verbose_every_arg <- 10L
}
parallel_strategy_arg <- if (length(args) >= 19) tolower(args[19]) else "outer_priority"
if (!parallel_strategy_arg %in% c("outer_priority", "balanced", "outer_only")) {
  warning(sprintf("Invalid parallel_strategy='%s'; falling back to 'outer_priority'.", parallel_strategy_arg))
  parallel_strategy_arg <- "outer_priority"
}

# ============================================================================
# Performance Configuration
# ============================================================================
# NLAMBDA_INIT: Number of lambda values for initial outcome model CV (glmnet).
# Lower values (e.g. 20) speed up fitting; higher values (default 100) give finer tuning.
# Set via environment variable NLAMBDA_INIT=20 (integer)
NLAMBDA_INIT <- as.integer(Sys.getenv("NLAMBDA_INIT", "100"))
if (NLAMBDA_INIT != 100L) {
  cat(sprintf("NLAMBDA_INIT = %d (default 100)\n", NLAMBDA_INIT))
}

# NESTED_PARALLEL: When TRUE, enables both simulation-level and source-site parallelism
# - Outer layer: parLapply (PSOCK cluster) for simulations
# - Inner layer: mclapply (fork) for source sites within each simulation
# Total cores = outer_cores × inner_cores (auto-balanced to not exceed available)
# Set via environment variable NESTED_PARALLEL=1 or NESTED_PARALLEL=TRUE
NESTED_PARALLEL <- Sys.getenv("NESTED_PARALLEL", "0") %in% c("1", "TRUE", "true", "True")
if (NESTED_PARALLEL) {
  cat("[PERF] NESTED PARALLEL ENABLED: simulation + source-site parallelism\n")
}

# FACE-C implementation status:
# - Canonical entry point: run_crossfit(..., communication_mode = "two_round"|"one_round")
# - Legacy non-crossfitting algorithm labels (one_round, two_round) are not used
#   in this script; cross-fitting variants are the supported workflow.

# Load the FACEC package (all R/ code + compiled C++)
# Prefer installed package; fall back to devtools::load_all for development
if (requireNamespace("FACEC", quietly = TRUE)) {
  library(FACEC)
} else if (requireNamespace("devtools", quietly = TRUE)) {
  cat("FACEC not installed; loading via devtools::load_all()...\n")
  devtools::load_all(".")
} else {
  stop("Neither installed FACEC package nor devtools found. ",
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
  config           = config_arg,
  dgp_type         = dgp_type_arg,
  ate_deviation    = ate_deviation_arg,
  n_deviated_sites = n_deviated_sites_arg
)

# Build setting identifier — always includes ALL parameters for full traceability.
# For the FACE paper DGP, the identifier also includes ate_deviation/n_deviated_sites.
# Format (facec):       n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_<alloc>_tf<transform>_ht<het>_ss<shift>_kf<folds>
# Format (face):  n<N>_K<K>_p<P>_<config>_<estimand>_<outcome>_face_dev<dev>_nd<nd>_kf<folds>
# Must match main.sh::build_setting_id() and main.cmd
if (dgp_type_arg == "face") {
  setting_id <- sprintf(
    "n%d_K%d_p%d_%s_%s_%s_face_dev%.1f_nd%d_kf%d",
    n_total_vec, K_vec, p_vec, config_arg, estimand_type_arg,
    outcome_type_arg, ate_deviation_arg, n_deviated_sites_arg, n_folds_arg
  )
} else {
  setting_id <- sprintf(
    "n%d_K%d_p%d_%s_%s_%s_%s_tf%s_ht%s_ss%.1f_kf%d",
    n_total_vec, K_vec, p_vec, config_arg,
    estimand_type_arg, outcome_type_arg, site_allocation_arg,
    transform_type_arg, heterogeneity_type_arg,
    shift_strength_arg, n_folds_arg
  )
}

cat(sprintf("\n============================================\n"))
cat(sprintf("FACE-C Simulation: %s\n", setting_id))
cat(sprintf("  n_total: %d\n", n_total_vec))
cat(sprintf("  K: %d\n", K_vec))
cat(sprintf("  p: %d\n", p_vec))
cat(sprintf("  config: %s\n", config_arg))
cat(sprintf("  estimand_type: %s\n", estimand_type_arg))
cat(sprintf("  dgp_type: %s\n", dgp_type_arg))
if (dgp_type_arg == "face") {
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
cat(sprintf("  job_id: %s\n", job_id))
cat(sprintf("============================================\n\n"))

# ============================================================================
# Checkpoint/Restart Configuration
# ============================================================================
# checkpoint functions (init_checkpoint_config, save_checkpoint, load_checkpoint,
# check_preempt_signal, signal_checkpoint_saved, cleanup_checkpoint) are
# provided by the FACEC package (R/checkpoint.R).
CHECKPOINT_DIR <- Sys.getenv("CHECKPOINT_DIR", "checkpoints")
ckpt_config <- init_checkpoint_config(CHECKPOINT_DIR, setting_id, job_id)

# Set up parallel processing (already imported by FACEC, but needed for main.R scope)
library(parallel)
library(doParallel)

# ============================================================================
# SIMULATION FUNCTIONS
# ============================================================================
# run_single_simulation(), summarize_results(), and run_simulation_study() are
# defined in R/simulation.R and are part of the FACEC package.
# ============================================================================

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Output directory with date stamp
OUTPUT_DIR <- "results"
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
}

# Run simulation for this specific setting
# Note: Each job runs ONE config + heterogeneity_type combination
cat(sprintf("\nRunning simulation for: %s\n", setting_id))

# Auto-detect available cores from SLURM or system
# SLURM sets SLURM_CPUS_PER_TASK, otherwise use detectCores() - 1
available_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "0"))
if (available_cores == 0) {
  available_cores <- max(1, parallel::detectCores() - 1)
}

# Guard against oversubscription in non-SLURM environments.
# By default, cap cores to 8 unless FACEC_MAX_CORES is explicitly set.
# Also never use more cores than n_sims to avoid idle workers.
facec_max_cores <- suppressWarnings(as.integer(Sys.getenv("FACEC_MAX_CORES", "8")))
if (is.na(facec_max_cores) || facec_max_cores < 1) {
  facec_max_cores <- Inf
}

available_cores <- min(available_cores, n_sims_arg, facec_max_cores)
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
  estimate_ate = (heterogeneity_type_arg != "none"),
  dgp_type = dgp_type_arg,
  ate_deviation = ate_deviation_arg,
  n_deviated_sites = n_deviated_sites_arg
)

results <- do.call(run_simulation_study, sim_args)

# Summarize results (with error handling)
summary_df <- tryCatch({
  summarize_results(results)
}, error = function(e) {
  cat(sprintf("[WARN] Could not summarize results: %s\n", e$message))
  data.frame()
})

if (nrow(summary_df) > 0) {
  print(summary_df)
} else {
  cat("No summary statistics available (empty results)\n")
}

# Save results with setting-specific filename
save_path <- file.path(OUTPUT_DIR, sprintf("%s.RData", setting_id))
tryCatch({
  save(results, summary_df, file = save_path)
  cat(sprintf("Results saved to: %s\n", save_path))
}, error = function(e) {
  cat(sprintf("[WARN] Failed to save RData: %s\n", e$message))
})

# Save results to CSV
csv_path <- file.path(OUTPUT_DIR, sprintf("%s.csv", setting_id))
tryCatch({
  write.csv(results, file = csv_path, row.names = FALSE)
  cat(sprintf("CSV saved to: %s\n", csv_path))
}, error = function(e) {
  cat(sprintf("[WARN] Failed to save CSV: %s\n", e$message))
})

# Save summary to CSV
summary_path <- file.path(OUTPUT_DIR, sprintf("%s_summary.csv", setting_id))
tryCatch({
  write.csv(summary_df, file = summary_path, row.names = FALSE)
  cat(sprintf("Summary saved to: %s\n", summary_path))
}, error = function(e) {
  cat(sprintf("[WARN] Failed to save summary CSV: %s\n", e$message))
})

# Clean up checkpoint after successful completion
cleanup_checkpoint(ckpt_config$file, ckpt_config$preempt_signal, ckpt_config$saved_signal)
  
cat(sprintf("\n============================================\n"))
cat(sprintf("Simulation completed: %s\n", setting_id))
cat(sprintf("============================================\n"))
