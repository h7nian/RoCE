#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
root <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(file.exists(file.path(root, "reviews/complete_v1/CHECKS_PASSED")),
  file.exists(file.path(root, "reviews/population_projection_v1/CHECKS_PASSED")),
  startsWith(output, paste0(root, "/reviews/")), !dir.exists(output))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
model_checks <- read.csv(file.path(dirname(configuration$template), "population_model_checks.csv"))
stopifnot(nrow(model_checks) == 6L, all(model_checks$radius3_log_margin > 0),
  all(model_checks$max_log_weight_linear_residual > 1e-3))
limits <- read.csv(file.path(root, "reviews/population_projection_v1/population_limits.csv"))
stopifnot(nrow(limits) == 12L, all(limits$gradient_max < 2e-7),
  all(limits$hessian_minimum_eigenvalue > 1e-9), all(limits$minimum_clipping_margin > .9))
original <- "/scratch.global/zhan9381/FACE-HD/real_data/rhc/model_validation_mc200_extension_v1/combined_mc200"
old_manifest <- read.csv(file.path(original, "manifest.csv"), stringsAsFactors = FALSE)
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 200L, identical(sort(manifest$sim_id), 1:200),
  all(manifest$config == "O_moderate_mix"))
parity <- list()
for (i in seq_len(nrow(manifest))) {
  task <- manifest[i, ]; directory <- file.path(root, "tasks", task$task_id)
  stopifnot(file.exists(file.path(directory, "COMPLETE")))
  saved <- read.csv(file.path(directory, "methods.csv"))
  old_task <- old_manifest[old_manifest$config == "O_case_mix" & old_manifest$sim_id == task$sim_id, ]
  stopifnot(nrow(old_task) == 1L)
  old <- read.csv(file.path(original, "tasks", old_task$task_id, "methods.csv"))
  for (method in c("Target-only", "Oracle-target")) {
    current <- saved[saved$method == method, c("estimate", "se", "truth")]
    previous <- old[old$method == method, c("estimate", "se", "truth")]
    delta <- max(abs(as.matrix(current) - as.matrix(previous)))
    stopifnot(delta < 1e-12)
    parity[[length(parity) + 1L]] <- data.frame(repeat_id = task$sim_id, method = method,
      max_difference = delta)
  }
}
dir.create(output, recursive = TRUE)
write.csv(do.call(rbind, parity), file.path(output, "target_law_parity.csv"), row.names = FALSE)
file.copy(file.path(root, "reviews/complete_v1/summary.csv"), file.path(output, "summary.csv"))
writeLines(c("All200 repeats complete; target draws/anchor/oracle andtruth reproduce the original OR law.",
  "Joint weights are misspecified in the working linear basis; the true OR is correct.",
  "At radius3,population calibration limits match the required moments,have positive Hessians andat least.9 log-predictor clipping margin.",
  "These are model/regularity diagnostics,not a clinical bias guarantee."), file.path(output, "CHECKS_PASSED"))
print(read.csv(file.path(output, "summary.csv")))
