#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (!(length(arguments) %in% c(3L, 4L))) stop("Usage: source_directory plan_csv scratch_output [profile]")
source_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(file.path(source_directory, "candidate_score_summary.R"))
source(file.path(source_directory, "target_common_projection.R"))
plan <- read.csv(arguments[[2L]], stringsAsFactors = FALSE)
profile <- if (length(arguments) == 4L) arguments[[4L]] else "main"
if (profile == "main") {
  selected <- plan[(plan$campaign == "lowdim_cutoff2_control_mc50_v1" & plan$p == 10) |
                   (plan$campaign == "validation_highdim_mc200_cutoff2_v1" & plan$p == 100), ]
  stopifnot(all(selected$K == 2L))
} else if (profile == "source_count") {
  selected <- plan[plan$campaign == "source_count_validation_v1/source_count_p100" & plan$rho == 0, ]
  stopifnot(all(selected$K %in% c(4L, 6L)))
} else stop("Unknown projection diagnostic profile")
stopifnot(nrow(selected) == 6L, all(selected$rho == 0))
output <- arguments[[3L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(selected, file.path(output, "selected_cases.csv"), row.names = FALSE)
records <- list()
for (index in seq_len(nrow(selected))) {
  task <- selected[index, ]
  directory <- file.path(task$root, "tasks", task$task_id)
  stopifnot(file.exists(file.path(directory, "COMPLETE")))
  saved <- readRDS(file.path(directory, "score_derivative.rds"))
  artifacts <- attr(saved, "roce_simulation_artifacts")
  stopifnot(identical(digest::digest(artifacts$data_split, algo = "sha256"),
                      readLines(file.path(directory, "data_sha256.txt"), warn = FALSE)))
  fit <- artifacts$direct_tate_results[[paste0(task$protocol, "_crossfit")]]
  candidate <- extract_candidate_score_summary(fit, artifacts$data_split)
  sources <- names(fit$source_estimates)
  covariance <- candidate$covariance[sources, sources, drop = FALSE]
  for (weighting in c("ipw", "unweighted")) {
    started <- proc.time()[["elapsed"]]
    projection <- crossfit_target_common_projection(fit, artifacts$data_split$t, weighting)
    elapsed <- proc.time()[["elapsed"]] - started
    approximation <- matrix(projection$common_variance, length(sources), length(sources)) +
      diag(candidate$source_private_variances[sources], nrow = length(sources))
    maximum_clipping <- max(vapply(unlist(projection$folds, recursive = FALSE),
      function(fold) fold$propensity_clipped_fraction, numeric(1L)))
    records[[length(records) + 1L]] <- data.frame(case_id = task$case_id, config = task$config,
      K = task$K, p = task$p, weighting = weighting, elapsed_seconds = elapsed,
      target_projection_variance = projection$common_variance,
      target_projection_variance_within_folds = projection$common_variance_within_folds,
      compatible_source_mean_variance = candidate$common_target$common_variance,
      total_covariance_relative_error = sqrt(sum((covariance - approximation)^2) / sum(covariance^2)),
      max_se_relative_error = max(abs(sqrt(diag(approximation) / diag(covariance)) - 1)),
      maximum_propensity_clipped_fraction = maximum_clipping)
    saveRDS(list(task = task, projection = projection, candidate = candidate),
      file.path(output, paste0("case_", task$case_id, "_", weighting, ".rds")))
    cat("Completed target projection", task$case_id, weighting, "in", elapsed, "seconds\n")
  }
  rm(saved, artifacts, fit, projection)
  gc(FALSE)
}
write.csv(do.call(rbind, records), file.path(output, "projection_diagnostics.csv"), row.names = FALSE)
writeLines("Six prespecified data settings and both projection conventions completed; not an inference validation.",
           file.path(output, "COMPLETE"))
