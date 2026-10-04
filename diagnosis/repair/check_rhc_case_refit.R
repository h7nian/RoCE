#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
.libPaths(c(arguments[1L], .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = arguments[1L]))
repository <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !dir.exists(output))
source(file.path(repository, "diagnosis/repair/rhc_validation_helpers.R"))
data <- rhc_validation_data()
reference <- readRDS("/scratch.global/zhan9381/FACE-HD/real_data/rhc/current_two_layer_v1/tasks/1/checkpoints/analysis.rds")$value$tate_fit
folds <- rhc_saved_folds(data, reference)
stopifnot(identical(folds, RoCE:::build_crossfit_folds(data, 10L)))
for (case in list(list(site = "s3", row = 423L), list(site = "t", row = 303L))) {
  changed <- rhc_remove_case(data, folds, case$site, case$row)
  stopifnot(sum(vapply(changed$data, `[[`, integer(1L), "n")) == 5038L)
  for (site in names(data)) {
    if (site != case$site) stopifnot(identical(changed$data[[site]], data[[site]]))
  }
  retained <- changed$retained_original_rows
  stopifnot(identical(changed$data[[case$site]]$W_outcome, data[[case$site]]$W_outcome[retained, , drop = FALSE]),
    identical(changed$data[[case$site]]$Y, data[[case$site]]$Y[retained]),
    identical(attr(changed$data, "feature_center"), attr(data, "feature_center")),
    identical(attr(changed$data, "feature_scale"), attr(data, "feature_scale")))
}
stopifnot(inherits(try(rhc_remove_case(data, folds, "s9", 1L), silent = TRUE), "try-error"),
          inherits(try(rhc_remove_case(data, folds, "t", 1.5), silent = TRUE), "try-error"))
dir.create(output, recursive = TRUE)
paths <- file.path(repository, "diagnosis/repair", c("rhc_validation_helpers.R", "run_rhc_case_refit.R", "check_rhc_case_refit.R"))
jsonlite::write_json(setNames(lapply(paths, function(path) digest::digest(path, algo = "sha256", file = TRUE)), paths),
                    file.path(output, "checked_files.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines("RHC_CASE_REFIT_INPUT_AND_OUTER_BOUNDARIES_PASSED", file.path(output, "CHECKS_PASSED"))
