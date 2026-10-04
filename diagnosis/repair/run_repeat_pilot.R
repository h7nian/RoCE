#!/usr/bin/env Rscript
# Compare calibration recipes on identical simulated data and compare all
# three aggregation objectives on the same fitted nuisance models per recipe.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1L], mustWork = TRUE)
stopifnot(startsWith(root, "/scratch.global/zhan9381/FACE-HD/"))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
.libPaths(c(configuration$library, .libPaths()))
library(RoCE)
stopifnot(identical(normalizePath(find.package("RoCE")),
                    normalizePath(file.path(configuration$library, "RoCE"))))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
task <- manifest[manifest$task_id == as.integer(args[2L]), , drop = FALSE]
stopifnot(nrow(task) == 1L)
output <- file.path(root, "tasks", args[2L])
options(error = function() {
  interrupted <- file.exists(file.path(output, "INTERRUPTED"))
  writeLines(if (interrupted) "INTERRUPTED" else "FAILED", file.path(output, "status.txt"))
  traceback(20L)
  quit(status = 1L)
})
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = as.character(configuration$cpus)))
stopifnot(is.finite(cores), cores >= 1L)
parallel_arms <- cores >= 2L
cores_per_arm <- if (parallel_arms) max(1L, cores %/% 2L) else cores
cache_parent <- NULL
if (isTRUE(configuration$shared_cache)) {
  if (!"nuisance_cache_dir" %in% names(formals(run_single_simulation))) {
    stop("The selected installed library does not support shared nuisance caches.")
  }
  cache_parent <- file.path(output, "nuisance_cache")
  stopifnot(dir.create(cache_parent))
}
writeLines("RUNNING", file.path(output, "status.txt"))
saveRDS(list(configuration = configuration, task = task, cores = cores,
             cores_per_arm = cores_per_arm,
             job_id = Sys.getenv("SLURM_JOB_ID"), host = Sys.info()[["nodename"]]),
        file.path(output, paste0("execution_", Sys.getenv("SLURM_RESTART_COUNT", "0"), ".rds")))
warnings <- list()
rows <- list()
failures <- character(0)
data_digest <- NULL
layer_rows <- layer_weights <- list()
if (isTRUE(configuration$layer_diagnostics)) source(file.path(root, "workflow", "layer_comparison_summary.R"))
dgp_type <- if (is.null(configuration$dgp_type)) "face" else configuration$dgp_type
dimension <- if (is.null(task$p)) configuration$p else task$p
save_fits <- if (is.null(configuration$save_fits)) "all" else configuration$save_fits
keep_fits <- save_fits == "all" || (save_fits == "first" && task$sim_id == configuration$first_seed)
target_method <- if (is.null(configuration$target_nuisance_method)) "hou_calibrated" else configuration$target_nuisance_method
source_validation <- if (is.null(configuration$source_validation_method)) "calibrated" else configuration$source_validation_method
source_program <- if (is.null(configuration$source_nuisance_method)) "calibrated" else configuration$source_nuisance_method
n_deviated_sites <- if (is.null(configuration$n_deviated_sites)) 1L else configuration$n_deviated_sites
stopifnot(n_deviated_sites >= 0L, n_deviated_sites <= task$K,
          n_deviated_sites == as.integer(n_deviated_sites))
for (recipe in configuration$recipes) {
  cat("Starting recipe", recipe, "for task", task$task_id, "at", format(Sys.time()), "\n")
  control <- if (recipe == "score_derivative") {
    if (target_method == "hou_calibrated") {
      list(recipe = recipe, target_propensity_initialization = "calibrated",
           target_radius = configuration$target_score_radius)
    } else list(recipe = recipe, target_propensity_initialization = "logistic")
  } else list(recipe = "legacy")
  if (source_program != "calibrated") control$source_nuisance_method <- source_program
  simulation_args <- list(
    sim_id = task$sim_id, n_total = configuration$n_per_site * (task$K + 1L),
    n_target = configuration$n_per_site, n_source_sizes = rep(configuration$n_per_site, task$K),
    K = task$K, p = dimension, config = task$config,
    n_folds = configuration$n_folds, nlambda_init = configuration$nlambda,
    dgp_type = dgp_type, outcome_type = "binary", estimand_type = "superpopulation",
    ate_deviation = task$rho, n_deviated_sites = as.integer(n_deviated_sites),
    deviation_mechanism = task$deviation_mechanism,
    methods = c(paste0(task$protocol, "_crossfit"), "target_only"),
    estimate_ate = TRUE, return_fitted_tate = TRUE,
    # The shared API requires >= 2 draws even when these requested estimators
    # do not use the comparison-estimator bootstrap. Keep its study setting.
    include_quadratic_bias_rule = FALSE, n_bootstrap = 5000L,
    n_cores_internal = cores_per_arm, parallel_treatment_arms = parallel_arms,
    nuisance_solver = "proximal_newton", nuisance_tol = configuration$nuisance_tol,
    target_nuisance_method = target_method, source_validation_method = source_validation,
    calibration_control = control, calibration_layout = "compact",
    additional_aggregation_modes = c("separate_arms", "joint_tate"),
    M_tau = configuration$source_radius, M_tau_inference = configuration$source_radius,
    aggregation_lambda = configuration$aggregation_lambda, verbose = TRUE)
  if (isTRUE(configuration$checkpoint)) {
    if (!"checkpoint_dir" %in% names(formals(run_single_simulation))) {
      stop("The frozen library does not support nuisance checkpoints.")
    }
    simulation_args$checkpoint_dir <- file.path(output, "checkpoints", recipe)
  }
  if (!is.null(cache_parent)) simulation_args$nuisance_cache_dir <- cache_parent
  if (!is.null(configuration$crossfit_layers)) simulation_args$crossfit_layers <- configuration$crossfit_layers
  if (isTRUE(configuration$nuisance_cv_certificate)) {
    if (!"nuisance_cv_certificate" %in% names(formals(run_single_simulation))) {
      stop("The frozen library does not support the CV certificate control.")
    }
    simulation_args$nuisance_cv_certificate <- TRUE
  }
  if (dgp_type == "bounded") simulation_args$dgp_control <- configuration$dgp_control
  result <- tryCatch(withCallingHandlers(do.call(run_single_simulation, simulation_args),
    warning = function(condition) {
      warnings[[length(warnings) + 1L]] <<- data.frame(recipe = recipe,
        message = conditionMessage(condition), stringsAsFactors = FALSE)
    }), error = identity)
  if (inherits(result, "error")) {
    saveRDS(result, file.path(output, paste0(recipe, ".rds")))
    failures <- c(failures, recipe)
    writeLines(conditionMessage(result), file.path(output, paste0(recipe, "_FAILED.txt")))
    next
  }
  if (!is.data.frame(result) || !nrow(result) ||
      any(!is.finite(result$estimate)) || any(!is.finite(result$se)) || any(result$se <= 0)) {
    saveRDS(result, file.path(output, paste0(recipe, ".rds")))
    writeLines("Nonfinite estimates or nonpositive standard errors.",
      file.path(output, paste0(recipe, "_FAILED.txt")))
    stop("Invalid numerical results; repeat retained for diagnosis.")
  }
  expected_methods <- c(paste0(task$protocol, "_crossfit_ate",
    c("", "_separate_arms", "_joint_tate")), "target_only_ate")
  if (target_method != "lasso") expected_methods <- c(expected_methods, "target_anchor_ate")
  if (!all(expected_methods %in% result$method)) {
    saveRDS(result, file.path(output, paste0(recipe, ".rds")))
    stop("A requested TATE method is missing; repeat retained for diagnosis.")
  }
  artifacts <- attr(result, "roce_simulation_artifacts")
  stopifnot(all(vapply(artifacts$data_split, function(site) site$n == 1000L, logical(1L))))
  current_digest <- digest::digest(artifacts$data_split, algo = "sha256")
  if (is.null(data_digest)) data_digest <- current_digest
  stopifnot(identical(data_digest, current_digest))
  if (isTRUE(configuration$layer_diagnostics)) {
    base_fit <- artifacts$direct_tate_results[[paste0(task$protocol, "_crossfit")]]
    if (is.null(artifacts$arm_truth)) stop("The layer-diagnostic library must retain both arm truths.")
    writeLines(layer_outer_fit_fingerprint(base_fit),
               file.path(output, paste0(recipe, "_outer_fit_sha256.txt")))
    saveRDS(layer_inference_packet(base_fit, artifacts$data_split, artifacts$arm_truth),
            file.path(output, paste0(recipe, "_layer_records.rds")))
    for (mode in c("common_tate", "separate_arms", "joint_tate")) {
      selected <- if (mode == "common_tate") base_fit else
        reaggregate_tate_crossfit(artifacts$data_split, base_fit, M_tau_inference = base_fit$M_tau_inference,
                                 aggregation_mode = mode, verbose = FALSE)
      detail <- summarize_layer_arms(selected, artifacts$data_split)
      values <- c(detail$estimates, tate = selected$estimate)
      variances <- c(diag(detail$covariance), tate = selected$variance)
      fixed <- c(diag(detail$fixed_covariance), tate = selected$variance_fixed_weights)
      truths <- c(artifacts$arm_truth, tate = artifacts$tate_truth)
      layer_rows[[length(layer_rows) + 1L]] <- data.frame(recipe = recipe, aggregation_mode = mode,
        crossfit_layers = selected$crossfit_levels, quantity = names(values), estimate = as.numeric(values),
        truth = as.numeric(truths[names(values)]), variance = as.numeric(variances),
        variance_fixed_weights = as.numeric(fixed), covariance_mu1_mu0 = detail$covariance[1L, 2L],
        covariance_mu1_mu0_fixed_weights = detail$fixed_covariance[1L, 2L])
      weight_matrix <- selected$fold_weights
      layer_weights[[length(layer_weights) + 1L]] <- data.frame(recipe = recipe, aggregation_mode = mode,
        crossfit_layers = selected$crossfit_levels, outer_fold = rep(seq_len(nrow(weight_matrix)), times = ncol(weight_matrix)),
        coordinate = rep(colnames(weight_matrix), each = nrow(weight_matrix)), weight = as.vector(weight_matrix))
    }
  }
  if (keep_fits) saveRDS(result, file.path(output, paste0(recipe, ".rds")))
  attr(result, "roce_simulation_artifacts") <- NULL
  if (!keep_fits) saveRDS(result, file.path(output, paste0(recipe, ".rds")))
  result$pilot_recipe <- recipe
  result$pilot_task_id <- task$task_id
  result$n_deviated_sites <- n_deviated_sites
  rows[[recipe]] <- result
}
if (length(rows)) {
  results <- do.call(rbind, rows)
  write.csv(results, file.path(output, "results.csv"), row.names = FALSE)
  references <- results[results$method == "target_only_ate", c("estimate", "se")]
  if (nrow(references) > 1L) {
    stopifnot(all(references$estimate == references$estimate[1L]),
              all(references$se == references$se[1L]))
  }
}
if (length(warnings)) write.csv(do.call(rbind, warnings), file.path(output, "warnings.csv"), row.names = FALSE)
if (length(layer_rows)) write.csv(do.call(rbind, layer_rows), file.path(output, "layer_estimates.csv"), row.names = FALSE)
if (length(layer_weights)) write.csv(do.call(rbind, layer_weights), file.path(output, "layer_weights.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines(if (length(failures)) "FAILED" else "COMPLETE", file.path(output, "status.txt"))
if (length(failures)) stop("Failed recipes: ", paste(failures, collapse = ", "), "; diagnostic outputs retained.")
writeLines(data_digest, file.path(output, "data_sha256.txt"))
completion <- file.path(output, "COMPLETE.tmp")
writeLines("All requested recipes completed; development pilot only.", completion)
stopifnot(file.rename(completion, file.path(output, "COMPLETE")))
# Committed results are skipped on restart; intermediate fits are no longer needed.
if (isTRUE(configuration$checkpoint)) unlink(file.path(output, "checkpoints"), recursive = TRUE)
