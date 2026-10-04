#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
output <- arguments[1L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
library_path <- "/scratch.global/zhan9381/FACE-HD/implementation/r12/density_cv_retry_v2/Rlib"
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
source("diagnosis/repair/rhc_validation_helpers.R")
source("diagnosis/repair/tate_mass_stabilization.R")
reference <- "/scratch.global/zhan9381/FACE-HD/real_data/rhc/source_truncation_v1/tasks/4/checkpoints/analysis.rds"
fit <- readRDS(reference)$value$tate_fit
data <- rhc_validation_data()
unscaled <- rescale_joint_tate(fit, data, rep(1, fit$n_folds))
anchor <- rescale_joint_tate(fit, data, rep(0, fit$n_folds))
stopifnot(abs(unscaled$estimate - fit$estimate) < 1e-12,
  abs(unscaled$se_first_order - fit$se) < 1e-12,
  abs(unscaled$se_fixed_weights - fit$se_fixed_weights) < 1e-12,
  abs(anchor$estimate - fit$target_only$estimate) < 1e-12,
  abs(anchor$se_first_order - fit$target_only$se) < 1e-12)
for (factor in c(.2, .5, .8)) {
  scaled <- rescale_joint_tate(fit, data, rep(factor, fit$n_folds))
  expected <- factor * unscaled$influence + (1 - factor) * anchor$influence
  stopifnot(max(abs(scaled$influence - expected)) < 1e-12,
    abs(scaled$estimate - (factor * fit$estimate + (1 - factor) * fit$target_only$estimate)) < 1e-12,
    scaled$se_first_order <= max(unscaled$se_first_order, anchor$se_first_order) + 1e-12)
}
candidate <- stabilize_joint_tate(fit, data, "weight_mass")
stopifnot(all(candidate$fold_scale > 0 & candidate$fold_scale <= 1))
dir.create(output, recursive = TRUE)
saveRDS(list(reference = unscaled, anchor = anchor, candidate = candidate), file.path(output, "diagnostic.rds"))
summary <- data.frame(method = c("Original joint TATE", "Calibrated target anchor", "Exploratory mass stabilization"),
  estimate = c(unscaled$estimate, anchor$estimate, candidate$estimate),
  se_first_order = c(unscaled$se_first_order, anchor$se_first_order, candidate$se_first_order))
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
write.csv(data.frame(fold = seq_len(fit$n_folds), scale = candidate$fold_scale,
  candidate$weight_mass, check.names = FALSE), file.path(output, "fold_scales.csv"), row.names = FALSE)
writeLines("Identity/anchor/influence-linearity and deployed residual reconstruction checks passed; this is not an empirical coverage validation.",
  file.path(output, "CHECKS_PASSED"))
print(summary)
