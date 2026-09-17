#!/usr/bin/env Rscript
# One rho=0 replicate of the production TATE design with an explicit fold count,
# for the n_folds coverage pilot. Mirrors run_direct_tate_rho_group_task.R's
# simulation_args except n_folds, rho (fixed at 0) and the pilot library.
#
# usage: pilot_nfolds_task.R SEED CONFIG K N_FOLDS OUT_DIR
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: pilot_nfolds_task.R SEED CONFIG K N_FOLDS OUT_DIR", call. = FALSE)
}
seed <- as.integer(args[[1L]])
config <- args[[2L]]
K <- as.integer(args[[3L]])
n_folds <- as.integer(args[[4L]])
out_dir <- args[[5L]]

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "resource_topology.R"))

allocated_cores <- roce_read_positive_integer_env("SLURM_CPUS_PER_TASK", 1L)
nuisance_cv_threads <- roce_read_positive_integer_env(
  "ROCE_NUISANCE_CV_THREADS", 1L, n_folds
)
resource_plan <- roce_slurm_resource_plan(
  source_count = K, n_folds = n_folds,
  allocated_cores = allocated_cores, nuisance_cv_threads = nuisance_cv_threads
)

out_path <- file.path(out_dir, sprintf("%s_K%d_f%d_seed%04d.csv", config, K, n_folds, seed))
if (file.exists(out_path)) {
  message("[skip] exists: ", out_path)
  quit(save = "no", status = 0L)
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

started <- proc.time()[["elapsed"]]
result <- run_single_simulation(
  sim_id = seed,
  n_total = 1000L * (K + 1L),
  K = K,
  p = 100L,
  config = config,
  methods = c("one_round_crossfit", "target_only"),
  verbose = FALSE,
  n_cores_internal = resource_plan$source_workers,
  nlambda_init = 100L,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  n_folds = n_folds,
  aggregation_lambda = 1,
  n_bootstrap = 5000L,
  M_tau = 5,
  M_tau_inference = 5,
  estimate_ate = TRUE,
  parallel_treatment_arms = resource_plan$parallel_treatment_arms,
  include_hard_threshold_diagnostic = FALSE,
  include_quadratic_bias_rule = TRUE,
  dgp_type = "face",
  ate_deviation = 0,
  n_deviated_sites = 0L,
  deviation_mechanism = "treated_arm"
)
result$pilot_n_folds <- n_folds
result$pilot_elapsed_seconds <- proc.time()[["elapsed"]] - started
result$pilot_allocated_cores <- allocated_cores
result$pilot_cv_threads <- nuisance_cv_threads
result$nuisance_solver <- RoCE:::nuisance_solver_cpp()
tmp <- paste0(out_path, ".tmp")
write.csv(result, tmp, row.names = FALSE)
file.rename(tmp, out_path)
message(sprintf("[done] seed=%d config=%s K=%d n_folds=%d elapsed=%.0fs -> %s",
                seed, config, K, n_folds, result$pilot_elapsed_seconds[[1L]], out_path))
