#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
output <- arguments[2L]
stopifnot(startsWith(output, paste0(root, "/")), !dir.exists(output))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
stopifnot(identical(Sys.setlocale("LC_COLLATE", configuration$collation), configuration$collation))
.libPaths(c(configuration$library, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = configuration$library))
source(file.path(root, "workflow/rhc_stage_signatures.R"))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
task <- manifest[manifest$role == "method", ]
stopifnot(nrow(task) == 2L, all(file.exists(file.path(root, "tasks", task$task_id, "COMPLETE"))))
fits <- setNames(lapply(task$task_id, function(id)
  readRDS(file.path(root, "tasks", id, "checkpoints/analysis.rds"))$value), task$preprocessing)
historical <- readRDS(file.path(dirname(root), "source_truncation_v1/tasks/2/checkpoints/analysis.rds"))$value$tate_fit
stopifnot(isTRUE(all.equal(rhc_stage_signatures(fits$cohort$tate_fit),
  rhc_stage_signatures(historical), tolerance = 1e-10)),
  identical(fits$cohort$metadata$data_sha256, fits$outer_fold$metadata$data_sha256))
indices <- function(fit, arm) lapply(fit$arm_results[[arm]]$intermediates$fold_info,
  function(fold) fold[c("target_idx", "source_idx")])
for (arm in c("mu1", "mu0")) {
  stopifnot(identical(indices(fits$cohort$tate_fit, arm), indices(fits$outer_fold$tate_fit, arm)),
    identical(indices(fits$cohort$tate_fit, arm), indices(historical, arm)))
}
cohort <- build_rhc_cohort(site_var = "ninsclas",
  site_recode = c("No insurance" = NA_character_, "Medicare & Medicaid" = NA_character_))
data <- build_rhc_data_split(cohort, K = 4L, target_site = "Private")
rows <- RoCE:::.rhc_cohort_rows_by_site(cohort, data)
parameters <- fits$outer_fold$tate_fit$preprocessing$parameters
outer <- indices(fits$outer_fold$tate_fit, "mu1")
for (fold in seq_along(parameters)) {
  held_out <- c(rows$t[outer[[fold]]$target_idx], unlist(Map(function(site, index)
    rows[[site]][index], names(data)[-1L], outer[[fold]]$source_idx), use.names = FALSE))
  stopifnot(identical(parameters[[fold]]$training_rows, setdiff(seq_len(nrow(cohort)), held_out)))
}
weights <- do.call(rbind, lapply(names(fits), function(mode) {
  fit <- fits[[mode]]$tate_fit
  data.frame(preprocessing = mode, outer_fold = seq_len(fit$n_folds), fit$fold_weights,
    check.names = FALSE)
}))
dir.create(output, recursive = TRUE)
write.csv(weights, file.path(output, "paired_fold_weights.csv"), row.names = FALSE)
jsonlite::write_json(list(historical_model_parity = TRUE, unchanged_patient_data = TRUE,
  unchanged_fold_memberships = TRUE, outer_preprocessing_rows_verified = TRUE,
  site_sizes = fits$cohort$metadata$site_n, initial_and_calibration_cv_preprocessing_refitted = FALSE),
  file.path(output, "pairing_checks.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines("HISTORICAL_MODEL_PARITY_AND_OUTER_BOUNDARIES_PASSED", file.path(output, "CHECKS_PASSED"))
