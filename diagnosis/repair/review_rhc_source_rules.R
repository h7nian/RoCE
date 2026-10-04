#!/usr/bin/env Rscript
# Compare predefined source-rule variants within each common outer partition.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
study <- normalizePath(arguments[1L], mustWork = TRUE)
output <- arguments[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
manifest <- read.csv(file.path(study, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(manifest) > 0L, !anyDuplicated(manifest$config), all(manifest$role == "method"),
          all(file.exists(file.path(study, "tasks", manifest$task_id, "COMPLETE"))))
manifest$variant <- sub("_fold[0-9]+$", "", sub("^RHC_", "", manifest$config))
variants <- unique(manifest$variant)
seeds <- unique(manifest$fold_seed)
stopifnot("min" %in% variants,
  all(variants %in% c("min", "final_weight_1se", "final_outcome_1se", "initial_weight_1se", "source_all_1se", "global_1se")),
  all(manifest$config == paste0("RHC_", manifest$variant, "_fold", manifest$fold_seed)),
  all(table(manifest$fold_seed, manifest$variant) == 1L))
fits <- signatures <- list()
estimates <- list()
for (index in seq_len(nrow(manifest))) {
  task <- manifest[index, ]
  directory <- file.path(study, "tasks", task$task_id)
  stopifnot(identical(readLines(file.path(directory, "status.txt")), "COMPLETE"))
  fit <- readRDS(file.path(directory, "checkpoints/analysis.rds"))$value
  if (index == 1L) expected_hash <- fit$metadata$data_sha256
  stopifnot(identical(fit$metadata$data_sha256, expected_hash))
  fits[[task$config]] <- fit
  signatures[[task$config]] <- readRDS(file.path(directory, "stage_signatures.rds"))
  table <- fit$methods
  table$profile <- task$config
  table$fold_seed <- task$fold_seed
  estimates[[task$config]] <- table
}

maximum_difference <- function(first, second) {
  if (is.list(first) || is.list(second)) {
    stopifnot(is.list(first), is.list(second), identical(names(first), names(second)), length(first) == length(second))
    return(if (length(first)) max(vapply(seq_along(first), function(i) maximum_difference(first[[i]], second[[i]]), numeric(1L))) else 0)
  }
  if (is.numeric(first) && is.numeric(second)) {
    stopifnot(length(first) == length(second), all(is.finite(first)), all(is.finite(second)))
    return(if (length(first)) max(abs(first-second)) else 0)
  }
  stopifnot(identical(first, second))
  0
}

checks <- comparisons <- list()
for (seed in seeds) {
  baseline <- paste0("RHC_min_fold", seed)
  for (variant in variants) {
    name <- paste0("RHC_", variant, "_fold", seed)
    fit <- fits[[name]]
    signature <- signatures[[name]]
    stopifnot(identical(fit$metadata$fold_sha256, fits[[baseline]]$metadata$fold_sha256))
    fields <- if (variant == "global_1se") character() else c("target", "initial_outcome_messages")
    if (variant == "final_weight_1se") fields <- c(fields, "initial_weight", "outcome")
    if (variant == "final_outcome_1se") fields <- c(fields, "initial_weight", "weight")
    for (field in fields) {
      difference <- maximum_difference(signature[[field]], signatures[[baseline]][[field]])
      if (difference > 1e-10) stop("Fixed-stage invariant failed: ", name, "/", field)
      checks[[length(checks)+1L]] <- data.frame(profile = name, stage = field, max_difference = difference)
    }
    expected <- list(initial_weight = "min", weight = "min", outcome = "min")
    if (variant == "initial_weight_1se") expected$initial_weight <- "1se"
    if (variant == "final_weight_1se") expected$weight <- "1se"
    if (variant == "final_outcome_1se") expected$outcome <- "1se"
    if (variant %in% c("source_all_1se", "global_1se")) expected[] <- "1se"
    for (rule in signature$selected_rules) for (stage in names(expected)) {
      stopifnot(length(rule[[stage]]) > 0L, all(rule[[stage]] == expected[[stage]]))
    }
    if (variant != "global_1se") {
      for (method in c("Target-only", "Calibrated target-only")) {
        first <- fit$methods[as.character(fit$methods$method) == method, c("estimate", "se")]
        second <- fits[[baseline]]$methods[as.character(fits[[baseline]]$methods$method) == method, c("estimate", "se")]
        stopifnot(max(abs(as.numeric(first[1L, ]) - as.numeric(second[1L, ]))) < 1e-10)
      }
    }
    comparisons[[length(comparisons)+1L]] <- data.frame(variant = variant, fold_seed = seed,
      estimate = fit$tate_fit$estimate, se = fit$tate_fit$se,
      estimate_difference = fit$tate_fit$estimate - fits[[baseline]]$tate_fit$estimate,
      se_ratio = fit$tate_fit$se/fits[[baseline]]$tate_fit$se)
  }
}
dir.create(output, recursive = TRUE)
write.csv(do.call(rbind, estimates), file.path(output, "estimates.csv"), row.names = FALSE)
write.csv(do.call(rbind, comparisons), file.path(output, "comparisons.csv"), row.names = FALSE)
write.csv(do.call(rbind, checks), file.path(output, "fixed_stage_checks.csv"), row.names = FALSE)
writeLines(c("RHC_SOURCE_RULE_COMPARISON_PASSED",
  sprintf("All %d predefined fits completed. Source-only variants preserve target fits and initial OR messages within each partition.", nrow(manifest)),
  "Final-weight-only changes preserve source OR; final-OR-only changes preserve final source weights.",
  sprintf("%d outer partitions reuse the same cohort; no Monte Carlo coverage, independent-replicate CI or lower-bias claim is implied.", length(seeds))),
  file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, comparisons))
