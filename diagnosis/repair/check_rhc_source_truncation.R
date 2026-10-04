#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
library_path <- normalizePath(arguments[1L], mustWork = TRUE)
repository <- normalizePath(arguments[2L], mustWork = TRUE)
output <- normalizePath(arguments[3L], mustWork = TRUE)
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"),
          file.exists(file.path(output, "INTERFACE_CHECKS_PASSED")))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
source(file.path(repository, "diagnosis/repair/rhc_stage_signatures.R"))
set.seed(92801)
generated <- generate_simulation_data(n_total = 600L, K = 1L, p = 4L,
  config = "C1", dgp_type = "bounded", outcome_type = "binary",
  estimand_type = "superpopulation", n_target = 300L, n_source_sizes = 300L)
data <- split_data_by_site(generated)
fits <- lapply(c(5, 2), function(radius) run_tate_crossfit(data, n_folds = 3L,
  communication_mode = "one_round", crossfit_layers = 2L,
  target_nuisance_method = "hou_calibrated", aggregation_mode = "joint_tate",
  nlambda_init = 6L, nuisance_lambda_rule = "min", nuisance_solver = "proximal_newton",
  nuisance_tol = 1e-10, nuisance_cv_certificate = TRUE,
  M_tau = radius, M_tau_inference = radius, lambda_selection = .5,
  calibration_layout = "compact", n_cores = 1L, verbose = FALSE,
  calibration_control = list(recipe = "score_derivative", target_radius = 5)))
signatures <- lapply(fits, rhc_stage_signatures)
for (field in c("target", "initial_outcome_messages")) {
  stopifnot(identical(signatures[[1L]][[field]], signatures[[2L]][[field]]))
}
stopifnot(identical(fits[[1L]]$target_only, fits[[2L]]$target_only),
          fits[[2L]]$M_tau == 2, fits[[2L]]$M_tau_inference == 2,
          fits[[2L]]$calibration_control$target_radius == 5)
saveRDS(lapply(fits, function(fit) list(estimate = fit$estimate, se = fit$se,
  source_radius = fit$M_tau, target = fit$target_only)), file.path(output, "source_truncation_check.rds"))
writeLines(c("SOURCE_TRUNCATION_TARGET_ISOLATION_PASSED",
  "Source fitting/evaluation radius5 versus2 preserves target scores and initial OR messages exactly."),
  file.path(output, "SOURCE_TRUNCATION_CHECKS_PASSED"))
