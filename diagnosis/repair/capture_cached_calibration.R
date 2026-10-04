#!/usr/bin/env Rscript
# Reconstruct a completed calibration task using only exact cache hits. This
# diagnostic never refits a model or modifies the live worker cache.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: capture_cached_calibration.R LIBRARY PILOT_ROOT TASK SITE ARM CALIBRATION_FOLDS NEW_OUTPUT")
}
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
capture_task <- function() {
  root <- args[2L]
  configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
  manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
  task <- manifest[manifest$task_id == as.integer(args[3L]), , drop = FALSE]
  stopifnot(nrow(task) == 1L, task$protocol == "one_round", isTRUE(configuration$shared_cache),
            identical(configuration$recipes, "score_derivative"))
  site <- args[4L]
  arm <- as.integer(args[5L])
  calibration_folds <- sort(as.integer(strsplit(args[6L], ",", fixed = TRUE)[[1L]]))
  output <- args[7L]
  stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
  dir.create(output, recursive = TRUE)
  dir.create(file.path(output, "cache_entries"))
  cache_root <- file.path(root, "tasks", args[3L], "nuisance_cache")
  directories <- list.dirs(cache_root, recursive = FALSE, full.names = TRUE)
  stopifnot(length(directories) > 0L)
  RoCE:::set_nuisance_solver_cpp("proximal_newton")
  set.seed(task$sim_id)
  data <- split_data_by_site(generate_simulation_data(
    n_total = configuration$n_per_site * (task$K + 1L), K = task$K, p = configuration$p,
    config = task$config, estimand_type = "superpopulation", dgp_type = "face",
    outcome_type = "binary", ate_deviation = task$rho, n_deviated_sites = 1L,
    deviation_mechanism = task$deviation_mechanism, n_target = configuration$n_per_site,
    n_source_sizes = rep(configuration$n_per_site, task$K), warn_ignored = FALSE))
  folds <- RoCE:::build_crossfit_folds(data, configuration$n_folds)
  stopifnot(site %in% names(data))
  original <- RoCE:::.fit_nuisance_training_subset
  on.exit(assignInNamespace(".fit_nuisance_training_subset", original, ns = "RoCE"), add = TRUE)
  requests <- list()
  lookup <- function(fitter, arguments, site, training_folds, fit_cache = NULL) {
    key <- RoCE:::.nuisance_training_key(fitter, site, training_folds, arguments$A_val)
    input_hash <- digest::digest(list(policy = RoCE:::NUISANCE_TRAINING_POLICY, key = key,
      inputs = list(arguments = arguments, solver = RoCE:::nuisance_solver_cpp())), algo = "sha256")
    paths <- file.path(directories, paste0(input_hash, ".rds"))
    paths <- paths[file.exists(paths)]
    if (length(paths) != 1L) stop("Exact cached task unavailable; no refit attempted: ", key)
    entry <- readRDS(paths)
    stopifnot(identical(entry$input_hash, input_hash),
              identical(entry$fit_hash, digest::digest(entry$fit, algo = "sha256")))
    saveRDS(entry, file.path(output, "cache_entries", basename(paths)))
    requests[[length(requests) + 1L]] <<- data.frame(fitter = fitter, site = site,
      arm = arguments$A_val, training_folds = paste(training_folds, collapse = ","), input_hash = input_hash)
    if (fitter == "fit_unified_outcome") {
      saveRDS(list(arguments = arguments, reference_fit = entry$fit,
        input_hash = input_hash, cv_seed = attr(entry$fit, "cv_seed"),
        site = site, training_folds = training_folds), file.path(output, "outcome_inputs.rds"))
    }
    entry$fit
  }
  assignInNamespace(".fit_nuisance_training_subset", lookup, ns = "RoCE")
  if (site == "t") {
    fit <- RoCE:::.fit_target_calibration(folds$target_folds, calibration_folds, arm,
      "binomial", configuration$source_radius, configuration$nlambda,
      lambda_rule = "min", layout = "compact", tol = configuration$nuisance_tol,
      calibration_control = list(recipe = "score_derivative",
        target_propensity_initialization = "calibrated", target_radius = configuration$target_score_radius))
  } else {
    messages <- RoCE:::.prepare_source_calibration_messages(folds$target_folds,
      folds$source_folds[[site]], site, calibration_folds, task$protocol, arm,
      "binomial", configuration$source_radius, configuration$nlambda, "min",
      calibration_recipe = "score_derivative")
    fit <- RoCE:::.fit_source_calibration(folds$source_folds[[site]], site,
      calibration_folds, function(site, fold) messages[[paste0("k2_", fold)]],
      arm, 1L, 1L, configuration$source_radius, configuration$nlambda, 10000L,
      "min", layout = "compact", tol = configuration$nuisance_tol,
      calibration_recipe = "score_derivative")
  }
  saveRDS(fit, file.path(output, "calibration.rds"))
  write.csv(do.call(rbind, requests), file.path(output, "exact_cache_hits.csv"), row.names = FALSE)
  writeLines(c("PASSED", "Every initial and final fitting input matched its live cache hash exactly.",
    "No fitting was executed; original caches were read only."), file.path(output, "status.txt"))
}
capture_task()
