#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
stopifnot(startsWith(root, "/scratch.global/zhan9381/FACE-HD/implementation/"))
library_path <- file.path(root, "Rlib")
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
source(file.path(root, "parity_workflow/rhc_validation_helpers.R"))
source(file.path(root, "parity_workflow/rhc_stage_signatures.R"))
study <- "/scratch.global/zhan9381/FACE-HD/real_data/rhc/source_truncation_v1"
manifest <- read.csv(file.path(study, "manifest.csv"), stringsAsFactors = FALSE)
data <- rhc_validation_data()
for (radius in c(3, 5)) {
  row <- manifest[manifest$config == paste0("RHC_source_radius", radius, "_fold0"), ]
  stopifnot(nrow(row) == 1L)
  reference <- readRDS(file.path(study, "tasks", row$task_id, "checkpoints/analysis.rds"))$value$tate_fit
  directory <- file.path(root, "rhc_parity", paste0("radius", radius))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, "fit.rds")
  if (file.exists(path)) fit <- readRDS(path) else {
    fit <- rhc_validation_fit(data, rhc_saved_folds(data, reference), radius,
      file.path(root, "rhc_parity/nuisance_checkpoints"), cores = 2L)
    saveRDS(fit, paste0(path, ".tmp")); stopifnot(file.rename(paste0(path, ".tmp"), path))
  }
  stopifnot(isTRUE(all.equal(rhc_fit_result(fit), rhc_fit_result(reference), tolerance = 1e-10)),
    isTRUE(all.equal(rhc_stage_signatures(fit), rhc_stage_signatures(reference), tolerance = 1e-10)),
    isTRUE(all.equal(fit$fold_weights, reference$fold_weights, tolerance = 1e-10)))
  writeLines("Saved real-data point/SE, all target/source models and eta reproduced.",
    file.path(directory, "PARITY_PASSED"))
}
writeLines("Both candidate radius3 and original radius5 reproduce saved Private/min fits.",
  file.path(root, "RHC_REFERENCE_PARITY_PASSED"))
