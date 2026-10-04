#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: source_directory plan_csv scratch_output")
source_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
source_names <- c("source_confidence_sets.R", "candidate_score_summary.R",
                  "arm_candidate_score_summary.R", "arm_pair_confidence_sets.R")
for (name in source_names) source(file.path(source_directory, name))
plan <- read.csv(arguments[[2L]], stringsAsFactors = FALSE)
output <- arguments[[3L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
if (dir.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
dir.create(file.path(output, "cases"))
rows <- availability <- list()
for (index in seq_len(nrow(plan))) {
  task <- plan[index, ]
  directory <- file.path(task$root, "tasks", task$task_id)
  input <- file.path(directory, "score_derivative.rds")
  ready <- file.exists(file.path(directory, "COMPLETE")) && file.exists(input)
  availability[[index]] <- data.frame(case_id = task$case_id, available = ready)
  if (!ready) next
  saved <- readRDS(input)
  artifacts <- attr(saved, "roce_simulation_artifacts")
  if (is.null(artifacts)) stop("Planned saved fit has no artifacts: ", input)
  actual_hash <- digest::digest(artifacts$data_split, algo = "sha256")
  stopifnot(identical(actual_hash, trimws(readLines(file.path(directory, "data_sha256.txt")))))
  fit <- artifacts$direct_tate_results[[paste0(task$protocol, "_crossfit")]]
  summary <- extract_joint_arm_score_summary(fit, artifacts$data_split)
  # The count is prespecified from the design, not inferred from agreement.
  # These plug-in intervals are structural diagnostics, not proved coverage.
  valid_mu1 <- ceiling(task$K / 2)
  valid_mu0 <- if (task$deviation_mechanism == "treated_arm") task$K else ceiling(task$K / 2)
  interval <- arm_pair_search_interval(summary$means, summary$covariance, valid_mu1, valid_mu0)
  row <- data.frame(case_id = task$case_id, config = task$config, K = task$K,
    p = task$p, rho = task$rho, protocol = task$protocol,
    estimate_identity_error = summary$estimate_identity_error,
    covariance_identity_error = summary$covariance_identity_error,
    source_message_identity_error = summary$source_message_identity_error,
    minimum_covariance_eigenvalue = min(eigen(summary$covariance, symmetric = TRUE, only.values = TRUE)$values),
    assumed_valid_mu1 = valid_mu1, assumed_valid_mu0 = valid_mu0,
    pair_count = nrow(interval$pairs), valid_pairs = interval$valid_pairs,
    votes_required = interval$votes_required, lower = interval$lower, upper = interval$upper,
    used_anchor_fallback = interval$used_anchor_fallback)
  rows[[length(rows) + 1L]] <- row
  saveRDS(list(task = task, summary = summary, interval = interval,
    data_sha256 = actual_hash, source_md5 = tools::md5sum(file.path(source_directory, source_names)),
    input_md5 = tools::md5sum(input)), file.path(output, "cases", paste0(task$case_id, ".rds")))
  rm(saved, artifacts, fit)
  gc(FALSE)
  cat("Checked joint-arm case", task$case_id, "(", index, "of", nrow(plan), ")\n")
}
write.csv(do.call(rbind, availability), file.path(output, "availability.csv"), row.names = FALSE)
if (length(rows)) write.csv(do.call(rbind, rows), file.path(output, "diagnostics.csv"), row.names = FALSE)
writeLines(paste("Reconstructed", length(rows), "of", nrow(plan),
  "saved joint-arm score problems. Empirical-score diagnostics, not validated fitted-method inference."),
  file.path(output, "COMPLETE"))
