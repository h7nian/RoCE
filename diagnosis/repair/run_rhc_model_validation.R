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
source(file.path(root, "workflow/rhc_model_validation.R"))
source(file.path(root, "workflow/rhc_stage_signatures.R"))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
task <- manifest[manifest$task_id == as.integer(arguments[2L]), ]
stopifnot(nrow(task) == 1L)
output <- file.path(root, "tasks", arguments[2L])
main <- function() {
  template <- readRDS(configuration$template)
  generated <- rhc_generate_population_sample(template, task$config, task$population_seed)
  folds <- RoCE:::build_crossfit_folds(generated$data, n_folds = 10L)
  input_hash <- digest::digest(list(configuration, task, generated, folds), algo = "sha256")
  data_hash <- digest::digest(generated$data, algo = "sha256")
  writeLines("RUNNING", file.path(output, "status.txt"))
  records <- list(); reference_signature <- NULL; target_reference <- NULL
  cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", configuration$cpus)) %/% 2L)
  for (radius in configuration$source_radii) {
    directory <- file.path(output, paste0("radius", radius))
    dir.create(directory, showWarnings = FALSE)
    path <- file.path(directory, "fit.rds")
    if (file.exists(path)) {
      saved <- readRDS(path)
      stopifnot(identical(saved$input_hash, input_hash), identical(saved$radius, radius))
      fit <- saved$value
    } else {
      set.seed(task$fit_seed)
      # Persistent entries include all fitting inputs in their keys. Reuse only
      # unchanged nuisances across radii through one checkpoint directory.
      fit <- rhc_validation_fit(generated$data, folds, radius, file.path(output, "nuisance_checkpoints"),
        cores = cores)
      saveRDS(list(input_hash = input_hash, radius = radius, value = fit), paste0(path, ".tmp"))
      stopifnot(file.rename(paste0(path, ".tmp"), path))
    }
    result <- rhc_fit_result(fit)
    stopifnot(all(is.finite(result$estimate)), all(is.finite(result$se)), all(result$se > 0))
    signature <- rhc_stage_signatures(fit)[c("target", "initial_outcome_messages")]
    if (is.null(reference_signature)) {
      reference_signature <- signature
      target_reference <- result[result$method == "Target-only", ]
    } else {
      stopifnot(isTRUE(all.equal(signature, reference_signature, tolerance = 1e-10)),
        isTRUE(all.equal(result[result$method == "Target-only", ], target_reference, tolerance = 1e-10)))
    }
    result <- result[result$method == "RoCE", ]
    result$source_radius <- radius
    records[[as.character(radius)]] <- result
    write.csv(fit$nuisance_fit_diagnostics, file.path(directory, "nuisance_fits.csv"), row.names = FALSE)
    retries <- list()
    for (arm in c("mu1", "mu0")) for (fold in seq_len(fit$n_folds)) {
      sources <- fit$arm_results[[arm]]$fold_results[[fold]]$source_results
      for (site in names(sources)) {
        models <- sources[[site]]$per_k2_gamma
        models$fold_sum <- sources[[site]]$gamma_s
        for (block in names(models)) {
          coefficients <- models[[block]]
          retry <- attr(coefficients, "cv_grid_retry")
          if (!is.null(retry)) retries[[length(retries) + 1L]] <- data.frame(arm = arm, site = site,
            stage = if (block == "fold_sum") "calibrated_weight" else "initial_weight",
            outer_fold = fold, calibration_block = block,
            training_folds = paste(attr(coefficients, "training_folds"), collapse = ","),
            original_max = retry$original_max, retry_max = retry$retry_max,
            added_lambdas = retry$added_lambdas, selected_lambda = attr(coefficients, "lambda_used"))
        }
      }
    }
    saveRDS(retries, file.path(directory, "cv_grid_retries.rds"))
    stopifnot(identical(digest::digest(generated$data, algo = "sha256"), data_hash))
  }
  target_reference$source_radius <- NA_real_
  oracle <- generated$oracle
  oracle_result <- data.frame(method = "Oracle-target", estimate = unname(oracle["estimate"]),
    se = unname(oracle["se"]), ci_lower = unname(oracle["estimate"] - qnorm(.975) * oracle["se"]),
    ci_upper = unname(oracle["estimate"] + qnorm(.975) * oracle["se"]), source_radius = NA_real_)
  result <- do.call(rbind, c(records, list(target_reference, oracle_result)))
  result$truth <- generated$truth; result$scenario <- task$config; result$repeat_id <- task$sim_id
  result$covered <- result$ci_lower <= result$truth & result$ci_upper >= result$truth
  write.csv(result, file.path(output, "methods.csv"), row.names = FALSE)
  jsonlite::write_json(list(input_hash = input_hash, data_sha256 = data_hash, task = task,
    truth = generated$truth, oracle_population_se = unname(oracle["population_se"])),
    file.path(output, "metadata.json"), auto_unbox = TRUE, pretty = TRUE)
  writeLines("PAIRED_TARGET_FITS_PASSED", file.path(output, "PAIRED_TARGET_FITS_PASSED"))
  writeLines("COMPLETE", file.path(output, "status.txt"))
  writeLines("All three radii and target/oracle checks complete", file.path(output, "COMPLETE"))
}
tryCatch(main(), error = function(error) {
  writeLines(conditionMessage(error), file.path(output, "analysis_FAILED.txt"))
  writeLines("FAILED", file.path(output, "status.txt")); stop(error)
})
