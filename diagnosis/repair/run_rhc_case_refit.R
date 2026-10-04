#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
stopifnot(startsWith(root, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
.libPaths(c(configuration$library, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = configuration$library))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(configuration$library, "RoCE"))))
source(file.path(root, "workflow/rhc_validation_helpers.R"))
source(file.path(root, "workflow/rhc_stage_signatures.R"))
source(file.path(root, "workflow/layer_comparison_summary.R"))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
task <- manifest[manifest$task_id == as.integer(arguments[2L]), ]
stopifnot(nrow(task) == 1L)
output <- file.path(root, "tasks", arguments[2L])
main <- function() {
  data <- rhc_validation_data()
  reference <- readRDS(file.path(task$reference, "checkpoints/analysis.rds"))$value$tate_fit
  folds <- rhc_saved_folds(data, reference)
  if (task$case_site != "none") {
    removed <- rhc_remove_case(data, folds, task$case_site, task$case_row)
    data <- removed$data; folds <- removed$folds
    saveRDS(removed$retained_original_rows, file.path(output, "retained_original_rows.rds"))
  }
  input_hash <- digest::digest(list(configuration, task, data, folds), algo = "sha256")
  checkpoint <- file.path(output, "checkpoints")
  dir.create(checkpoint, showWarnings = FALSE)
  path <- file.path(checkpoint, "analysis.rds")
  writeLines("RUNNING", file.path(output, "status.txt"))
  if (file.exists(path)) {
    saved <- readRDS(path); stopifnot(identical(saved$input_hash, input_hash)); fit <- saved$value
  } else {
    cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", configuration$cpus)) %/% 2L
    fit <- rhc_validation_fit(data, folds, task$radius, file.path(checkpoint, "nuisance"), cores)
    saveRDS(list(input_hash = input_hash, value = fit), paste0(path, ".tmp"))
    stopifnot(file.rename(paste0(path, ".tmp"), path))
  }
  result <- rhc_fit_result(fit)
  stopifnot(all(is.finite(result$estimate)), all(is.finite(result$se)), all(result$se > 0))
  if (task$case_site == "none") {
    stopifnot(isTRUE(all.equal(result, rhc_fit_result(reference), tolerance = 1e-10)),
      isTRUE(all.equal(rhc_stage_signatures(fit), rhc_stage_signatures(reference), tolerance = 1e-10)),
      isTRUE(all.equal(fit$fold_weights, reference$fold_weights, tolerance = 1e-10)))
    writeLines("FULL_COHORT_REFERENCE_PARITY_PASSED", file.path(output, "REFERENCE_PARITY_PASSED"))
  } else if (task$case_site != "t") {
    current_signature <- rhc_stage_signatures(fit)
    reference_signature <- rhc_stage_signatures(reference)
    for (field in c("target", "initial_outcome_messages")) stopifnot(isTRUE(all.equal(
      current_signature[[field]], reference_signature[[field]], tolerance = 1e-10)))
  }
  write.csv(result, file.path(output, "methods.csv"), row.names = FALSE)
  write.csv(fit$nuisance_fit_diagnostics, file.path(output, "nuisance_fits.csv"), row.names = FALSE)
  saveRDS(summarize_layer_arms(fit, data, keep_influence = TRUE), file.path(output, "arm_influence.rds"))
  saveRDS(rhc_stage_signatures(fit), file.path(output, "stage_signatures.rds"))
  jsonlite::write_json(list(task = task, input_hash = input_hash,
    data_sha256 = digest::digest(data, algo = "sha256"), site_sizes = vapply(data, `[[`, integer(1L), "n"),
    preprocessing = "Original full-cohort preprocessing fixed; all nuisances and eta refitted",
    scope = "Case sensitivity only; full-cohort primary results retained"),
    file.path(output, "metadata.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines("COMPLETE", file.path(output, "status.txt"))
  writeLines("RHC complete-refit case diagnostic finished", file.path(output, "COMPLETE"))
}
tryCatch(main(), error = function(error) {
  writeLines(conditionMessage(error), file.path(output, "analysis_FAILED.txt"))
  writeLines("FAILED", file.path(output, "status.txt")); stop(error)
})
