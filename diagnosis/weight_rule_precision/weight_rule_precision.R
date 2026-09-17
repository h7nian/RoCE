#!/usr/bin/env Rscript
# HISTORY: 2026-09-10 #0011 Precision of the weight rules at production scale
# Task:    Run the production simulation path for N replications at one
#          (config, K, rho) cell and report bias, empirical SD, mean SE, the
#          SE/SD ratio, the rejection rate and coverage for the primary
#          truncated-Wald rule A, the smooth quadratic-bias rule B, the
#          arm-wise diagnostic and the target-only anchor.
#
#          It calls run_single_simulation() - the same entry point production
#          uses - with only the cross-fit and target-only methods requested, so
#          the comparison estimators' multiplier bootstrap is skipped. The
#          rows are therefore identical to production's for these methods.
#
#          Usage: weight_rule_precision.R OUTPUT_DIR CONFIG K RHO N_REPLICATIONS
#                 [SEED_OFFSET]

source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

P <- 100L
N_SITE <- 1000L
N_FOLDS <- 5L
RULE_LABELS <- c(
  one_round_crossfit_ate = "A_truncated_wald",
  one_round_crossfit_ate_quadratic_bias = "B_quadratic_bias",
  one_round_crossfit_ate_armwise = "armwise_diagnostic",
  target_only_ate = "target_only"
)

.one_replication <- function(sim_id, config, K, rho) {
  run_single_simulation(
    sim_id = sim_id,
    n_total = N_SITE * (K + 1L),
    K = K,
    p = P,
    config = config,
    methods = c("one_round_crossfit", "target_only"),
    verbose = FALSE,
    n_cores_internal = 1L,
    nlambda_init = 100L,
    estimand_type = "superpopulation",
    outcome_type = "binary",
    n_folds = N_FOLDS,
    aggregation_lambda = 1,
    M_tau = 5,
    M_tau_inference = 5,
    estimate_ate = TRUE,
    dgp_type = "face",
    ate_deviation = rho,
    n_deviated_sites = if (rho > 0) 1L else 0L,
    nuisance_lambda_rule = "min",
    deviation_mechanism = "treated_arm",
    include_quadratic_bias_rule = TRUE
  )
}

.summarize <- function(rows, config, K, rho) {
  z <- stats::qnorm(0.975)
  do.call(rbind, lapply(names(RULE_LABELS), function(method) {
    g <- rows[rows$method == method, , drop = FALSE]
    if (nrow(g) < 2L) return(NULL)
    error <- g$estimate - g$truth
    sd_hat <- stats::sd(g$estimate)
    data.frame(
      config = config, K = K, rho = rho, rule = RULE_LABELS[[method]],
      method = method, n = nrow(g),
      bias = mean(error), bias_mcse = sd_hat / sqrt(nrow(g)),
      empirical_sd = sd_hat, mean_se = mean(g$se),
      se_to_empirical_sd = mean(g$se) / sd_hat,
      # Report interval rejection directly alongside the SE/SD ratio.
      rejection_rate = mean(abs(error) > z * g$se),
      coverage = mean(g$coverage),
      rmse = sqrt(mean(error^2)),
      excess_kurtosis = mean((error - mean(error))^4) / stats::var(error)^2 - 3,
      stringsAsFactors = FALSE
    )
  }))
}

.validate_rule_rows <- function(rows, sim_id, settings) {
  required <- c("sim_id", "method", "estimate", "se", "truth", "coverage",
                names(settings))
  if (!is.data.frame(rows) || !all(required %in% names(rows)) ||
      nrow(rows) != length(RULE_LABELS) || anyDuplicated(rows$method) ||
      !setequal(rows$method, names(RULE_LABELS)) ||
      anyNA(rows$sim_id) || any(rows$sim_id != sim_id)) {
    stop("incomplete or mismatched replicate ", sim_id, call. = FALSE)
  }
  for (field in names(settings)) {
    if (anyNA(rows[[field]]) ||
        any(as.character(rows[[field]]) != as.character(settings[[field]]))) {
      stop("replicate ", sim_id, " has mismatched ", field, call. = FALSE)
    }
  }
  if (!all(vapply(rows[c("estimate", "se", "truth")], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(rows[c("estimate", "se", "truth")])) ) ||
      any(rows$se <= 0) || !is.logical(rows$coverage) || anyNA(rows$coverage) ||
      any(rows$coverage != (abs(rows$estimate - rows$truth) <= stats::qnorm(0.975) * rows$se))) {
    stop("invalid numerical result in replicate ", sim_id, call. = FALSE)
  }
  rows
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% c(5L, 6L)) {
    stop("usage: weight_rule_precision.R OUTPUT_DIR CONFIG K RHO N_REPLICATIONS [SEED_OFFSET]",
         call. = FALSE)
  }
  config <- args[[2L]]
  K <- as.integer(args[[3L]])
  rho <- as.numeric(args[[4L]])
  n_replications <- as.integer(args[[5L]])
  seed_offset <- if (length(args) == 6L) as.integer(args[[6L]]) else 0L
  if (!config %in% c("C1", "C2", "C3", "C4") || !K %in% c(2L, 4L, 8L) ||
      !is.finite(rho) || rho < 0 || is.na(n_replications) || n_replications < 2L ||
      is.na(seed_offset) || seed_offset < 0L) {
    stop("invalid CONFIG, K, RHO, N_REPLICATIONS or SEED_OFFSET.", call. = FALSE)
  }
  output <- file.path(
    args[[1L]], sprintf("%s_K%d_rho%s_offset%05d", config, K, format(rho), seed_offset)
  )
  if (file.exists(output)) stop("output already exists: ", output, call. = FALSE)

  project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
  if (!nzchar(project_library)) stop("ROCE_PROJECT_LIB is required.", call. = FALSE)
  .libPaths(c(project_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  provenance <- roce_runtime_package_provenance(project_library)
  workflow <- roce_sha256_environment("ROCE_RULE_WORKFLOW_FINGERPRINT")
  settings <- list(config = config, K = K, p = P, n_total = N_SITE * (K + 1L),
                   n_folds = N_FOLDS, nlambda_init = 100L, M_tau = 5,
                   M_tau_inference = 5, nuisance_lambda_rule = "min",
                   aggregation_lambda = 1, deviation_mechanism = "treated_arm",
                   package_library = provenance$library,
                   package_fingerprint = provenance$fingerprint,
                   workflow_fingerprint = workflow, rho = rho)
  workers <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  if (is.na(workers)) stop("invalid SLURM_CPUS_PER_TASK.", call. = FALSE)
  ids <- seed_offset + seq_len(n_replications)
  # Each finished replicate is written immediately to its own file, so a job
  # that reaches its time limit keeps every replicate it completed and a rerun
  # skips them. (The first version wrote only at the end and lost all 2400
  # core-hours of a timed-out run.)
  replicate_directory <- file.path(args[[1L]], "replicates",
                                   sprintf("%s_K%d_rho%s", config, K, format(rho)))
  dir.create(replicate_directory, recursive = TRUE, showWarnings = FALSE)
  replicate_path <- function(sim_id) {
    file.path(replicate_directory, sprintf("replicate_%06d.csv", sim_id))
  }
  read_replicate <- function(sim_id) {
    .validate_rule_rows(utils::read.csv(replicate_path(sim_id), stringsAsFactors = FALSE),
                        sim_id, settings)
  }
  pending <- ids[!file.exists(vapply(ids, replicate_path, ""))]
  # Validate every cached row before launching any replacement work.
  invisible(lapply(setdiff(ids, pending), read_replicate))
  message(sprintf(
    "[weight rules] %s K=%d rho=%s seeds %d-%d: %d done, %d pending, workers=%d",
    config, K, format(rho), min(ids), max(ids),
    length(ids) - length(pending), length(pending), workers
  ))
  failures <- parallel::mclapply(pending, function(sim_id) {
    tryCatch({
      rows <- .one_replication(sim_id, config, K, rho)
      rows <- rows[rows$method %in% names(RULE_LABELS), , drop = FALSE]
      for (field in c("package_library", "package_fingerprint", "workflow_fingerprint", "rho", "n_folds")) {
        rows[[field]] <- settings[[field]]
      }
      rows <- .validate_rule_rows(rows, sim_id, settings)
      path <- replicate_path(sim_id)
      temporary <- tempfile(paste0(basename(path), "."), tmpdir = dirname(path))
      on.exit(unlink(temporary), add = TRUE)
      write.csv(rows, temporary, row.names = FALSE)
      if (!file.rename(temporary, path)) stop("could not commit ", path)
      NULL
    }, error = function(e) sprintf("%d: %s", sim_id, conditionMessage(e)))
  }, mc.cores = max(1L, min(workers, length(pending))))
  failures <- unlist(failures)
  if (length(failures) > 0L) {
    stop("replication failure(s): ", paste(failures, collapse = "; "), call. = FALSE)
  }
  rows <- do.call(rbind, lapply(ids, read_replicate))
  summary <- .summarize(rows, config, K, rho)
  print(summary[, c("rule", "n", "bias", "empirical_sd", "mean_se",
                    "se_to_empirical_sd", "rejection_rate", "coverage", "rmse")],
        row.names = FALSE, digits = 3)

  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "replication_rows.csv"), row.names = FALSE)
    write.csv(summary, file.path(stage, "rule_summary.csv"), row.names = FALSE)
    writeLines(c("task=weight_rule_precision", "history_entry=0011",
                 paste0("config=", config), paste0("K=", K), paste0("rho=", rho),
                 paste0("n_replications=", n_replications),
                 paste0("seed_offset=", seed_offset),
                 paste0(names(settings), "=", unlist(settings))),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "weight rule precision diagnostic")
}

if (sys.nframe() == 0L) main()
