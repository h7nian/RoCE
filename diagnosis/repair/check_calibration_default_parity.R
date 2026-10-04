#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: installed_library scratch_output_rds")
library_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
if (!startsWith(arguments[[2L]], "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
set.seed(4604)
data <- split_data_by_site(generate_bounded_data(n_target = 1000, n_source_sizes = c(1000, 1000),
  K = 2, p = 4, config = "C2", n_deviated_sites = 0))
fit <- run_tate_crossfit(data, n_folds = 4, communication_mode = "one_round", nlambda_init = 6,
  target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
  nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, calibration_layout = "compact",
  calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "calibrated",
                             target_radius = 12),
  M_tau = 12, M_tau_inference = 12, lambda_selection = .5, aggregation_mode = "joint_tate", verbose = FALSE)
parameters <- lapply(fit$arm_results, function(arm) lapply(arm$fold_results, function(fold)
  lapply(fold$source_results, function(source) list(weight = as.numeric(source$gamma_s),
                                                  outcome = as.numeric(source$alpha_ts)))))
saveRDS(list(data_hash = digest::digest(data, algo = "sha256"), estimate = fit$estimate,
  se = fit$se, weights = fit$fold_weights, target = fit$target_only,
  source_estimates = fit$source_estimates, parameters = parameters), arguments[[2L]])
cat("DEFAULT_CALIBRATION_REFERENCE_COMPLETED\n")
