#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: source_file plan_csv scratch_output")
source_file <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(source_file)
plan <- read.csv(arguments[[2L]], stringsAsFactors = FALSE)
output <- arguments[[3L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "cases"), showWarnings = FALSE)
records <- list()
candidate_records <- list()
availability <- list()
for (index in seq_len(nrow(plan))) {
  task <- plan[index, ]
  directory <- file.path(task$root, "tasks", task$task_id)
  input <- file.path(directory, "score_derivative.rds")
  available <- file.exists(file.path(directory, "COMPLETE")) && file.exists(input)
  availability[[index]] <- data.frame(case_id = task$case_id, available = available)
  if (!available) next
  saved <- readRDS(input)
  artifacts <- attr(saved, "roce_simulation_artifacts")
  if (is.null(artifacts)) stop("Planned first-seed result lacks a complete fit: ", input)
  actual_hash <- digest::digest(artifacts$data_split, algo = "sha256")
  expected_hash <- readLines(file.path(directory, "data_sha256.txt"), warn = FALSE)
  stopifnot(identical(actual_hash, expected_hash))
  fit <- artifacts$direct_tate_results[[paste0(task$protocol, "_crossfit")]]
  if (is.null(fit)) stop("The requested communication protocol is absent from the saved artifact")
  sources <- setdiff(names(artifacts$data_split), "t")
  shifted <- if (task$rho == 0 || task$n_deviated_sites == 0) character() else
    sources[seq_len(task$n_deviated_sites)]
  compatible <- setdiff(sources, shifted)
  summary <- extract_candidate_score_summary(fit, artifacts$data_split, compatible)
  population_scope <- task$config %in% c("C1", "C2", "C3")
  factor <- summary$common_target
  within <- summary$common_target_within_folds
  covariance_error <- max_se_relative_error <- NA_real_
  if (length(compatible) >= 2L) {
    covariance <- summary$covariance[compatible, compatible, drop = FALSE]
    approximation <- matrix(factor$common_variance, length(compatible), length(compatible)) +
      diag(summary$source_private_variances[compatible], nrow = length(compatible))
    covariance_error <- sqrt(sum((covariance - approximation)^2) / sum(covariance^2))
    max_se_relative_error <- max(abs(sqrt(diag(approximation) / diag(covariance)) - 1))
  }
  metadata <- data.frame(case_id = task$case_id, campaign = task$campaign, config = task$config,
    K = task$K, p = task$p, rho = task$rho, seed = task$sim_id, protocol = task$protocol,
    n_deviated_sites = task$n_deviated_sites, nuisance_correctness_supported = population_scope,
    comparison_count = length(compatible), data_sha256 = actual_hash)
  records[[length(records) + 1L]] <- cbind(metadata,
    data.frame(estimate_identity_error = summary$estimate_identity_error,
      anchor_variance_identity_error = summary$anchor_variance_identity_error,
      source_message_identity_error = summary$source_message_identity_error,
      target_covariance_relative_error = factor$relative_covariance_error,
      within_fold_target_covariance_relative_error = within$relative_covariance_error,
      leading_eigenvalue_fraction = factor$leading_eigenvalue_fraction,
      minimum_target_correlation = factor$minimum_correlation,
      prediction_disagreement_rmse = factor$prediction_disagreement_rmse,
      total_covariance_relative_error = covariance_error, max_se_relative_error = max_se_relative_error))
  candidate_records[[length(candidate_records) + 1L]] <- data.frame(
    case_id = task$case_id, candidate = names(summary$estimates),
    estimate = as.numeric(summary$estimates), truth = artifacts$tate_truth,
    empirical_score_se = sqrt(diag(summary$covariance)),
    target_variance = diag(summary$target_covariance),
    private_variance = c(0, summary$source_private_variances),
    outcome_compatible = c(TRUE, sources %in% compatible))
  saveRDS(list(task = task, summary = summary, data_sha256 = actual_hash,
    source_file_md5 = unname(tools::md5sum(source_file)), input_md5 = unname(tools::md5sum(input))),
    file.path(output, "cases", paste0(task$case_id, ".rds")))
  rm(saved, artifacts, fit)
  gc(FALSE)
  cat("Checked candidate-score case", task$case_id, "of", nrow(plan), "\n")
}
write.csv(do.call(rbind, availability), file.path(output, "availability.csv"), row.names = FALSE)
if (length(records)) {
  write.csv(do.call(rbind, records), file.path(output, "factor_diagnostics.csv"), row.names = FALSE)
  write.csv(do.call(rbind, candidate_records), file.path(output, "candidate_estimates.csv"), row.names = FALSE)
}
writeLines(paste("Checked", length(records), "of", nrow(plan),
  "prespecified first-seed fits. Structural diagnostics only; not an MC coverage estimate."),
  file.path(output, "COMPLETE"))
cat("CANDIDATE_SCORE_DIAGNOSTICS_PASSED\n")
