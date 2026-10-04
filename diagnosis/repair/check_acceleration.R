#!/usr/bin/env Rscript
# Paired three-level acceleration check. This is a numerical regression case,
# not a Monte Carlo coverage experiment. All variants use identical data/folds.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:4) {
  stop("usage: check_acceleration.R INSTALLED_LIBRARY NEW_OUTPUT [NUISANCE_TOL] [legacy|score_derivative]")
}
recipe <- if (length(args) == 4L) match.arg(args[4L], c("legacy", "score_derivative")) else "legacy"
control <- if (recipe == "score_derivative") {
  list(recipe = recipe, target_propensity_initialization = "calibrated", target_radius = log(9))
} else list(recipe = "legacy")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
cache_parent <- file.path(output, "nuisance_cache")
dir.create(cache_parent)
configuration <- list(seed = 1081L, n_per_site = 1000L, K = 2L, p = 4L,
  n_folds = 4L, nlambda = 8L, outcome_family = "binomial",
  nuisance_tol = if (length(args) >= 3L) as.numeric(args[3L]) else 1e-10,
  calibration_control = control,
  estimate_tolerance = 1e-6, se_tolerance = 1e-6, weight_tolerance = 1e-5,
  coefficient_tolerance = 1e-4, lambda_tolerance = 1e-10,
  scope = "small-basis numerical check of three-level pipeline; not p100 performance evidence")
saveRDS(configuration, file.path(output, "configuration.rds"))
writeLines(capture.output(dput(configuration)), file.path(output, "configuration.txt"))
set.seed(configuration$seed)
data <- split_data_by_site(generate_simulation_data(
  n_total = 3000L, K = 2L, p = 4L, config = "C1", dgp_type = "face",
  outcome_type = "binary", n_target = 1000L, n_source_sizes = c(1000L, 1000L),
  warn_ignored = FALSE))
folds <- RoCE:::build_crossfit_folds(data, configuration$n_folds)
saveRDS(list(data = data, folds = folds), file.path(output, "inputs.rds"))
variants <- list(
  reference = list(nuisance_solver = "coordinate_descent", use_lambda_cache = FALSE,
                   calibration_layout = "block", n_cores = 1L, parallel_arms = FALSE),
  newton = list(nuisance_solver = "proximal_newton", use_lambda_cache = FALSE,
                calibration_layout = "block", n_cores = 1L, parallel_arms = FALSE),
  cached = list(nuisance_solver = "proximal_newton", use_lambda_cache = TRUE,
                calibration_layout = "block", n_cores = 1L, parallel_arms = FALSE),
  compact = list(nuisance_solver = "proximal_newton", use_lambda_cache = TRUE,
                 calibration_layout = "compact", n_cores = 1L, parallel_arms = FALSE),
  parallel = list(nuisance_solver = "proximal_newton", use_lambda_cache = TRUE,
                  calibration_layout = "compact", n_cores = 2L, parallel_arms = TRUE),
  shared_cache = list(nuisance_solver = "proximal_newton", use_lambda_cache = TRUE,
                      calibration_layout = "compact", n_cores = 2L, parallel_arms = TRUE,
                      nuisance_cache_dir = cache_parent))
driver_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(driver_file)), "acceleration_helpers.R"))

reference <- NULL
comparisons <- model_comparisons <- timings <- list()
for (name in names(variants)) {
  cat("Starting", name, "\n")
  started <- proc.time()[["elapsed"]]
  fitted <- tryCatch(do.call(run_tate_crossfit, c(list(
    data_split = data, precomputed_folds = folds, n_folds = configuration$n_folds,
    communication_mode = "one_round", nlambda_init = configuration$nlambda,
    family = "binomial", target_nuisance_method = "hou_calibrated",
    source_validation_method = "calibrated", nuisance_tol = configuration$nuisance_tol,
    calibration_control = configuration$calibration_control, verbose = FALSE), variants[[name]])), error = identity)
  if (inherits(fitted, "error")) {
    writeLines(conditionMessage(fitted), file.path(output, paste0(name, "_FAILED.txt")))
    stop(conditionMessage(fitted))
  }
  stopifnot(identical(fitted$crossfit_levels, 3L))
  timings[[name]] <- data.frame(variant = name, elapsed_seconds = proc.time()[["elapsed"]] - started)
  modes <- setNames(lapply(c("common_tate", "separate_arms", "joint_tate"), function(mode) {
    calculate_tate_crossfit_aggregation(data, fitted$arm_results$mu1, fitted$arm_results$mu0,
      aggregation_mode = mode, verbose = FALSE)
  }), c("common_tate", "separate_arms", "joint_tate"))
  models <- collect_calibration_models(fitted)
  saveRDS(list(fit = fitted, modes = modes), file.path(output, paste0(name, ".rds")))
  if (is.null(reference)) reference <- list(models = models, modes = modes)
  stopifnot(identical(names(models), names(reference$models)))
  model_comparisons[[name]] <- do.call(rbind, lapply(names(models), function(key) {
    a <- models[[key]]
    b <- reference$models[[key]]
    lambda_a <- attr(a, "lambda_used")
    lambda_b <- attr(b, "lambda_used")
    lambda_equal <- (is.na(lambda_a) && is.na(lambda_b)) ||
      abs(lambda_a - lambda_b) <= configuration$lambda_tolerance * max(1, abs(lambda_b))
    data.frame(variant = name, model = key, lambda = lambda_a, reference_lambda = lambda_b,
      lambda_equal = lambda_equal,
      max_coefficient_difference = max(abs(as.numeric(a) - as.numeric(b))),
      kkt = if (is.null(attr(a, "kkt_residual"))) NA_real_ else attr(a, "kkt_residual"))
  }))
  comparisons[[name]] <- do.call(rbind, lapply(names(modes), function(mode) {
    a <- modes[[mode]]
    b <- reference$modes[[mode]]
    data.frame(variant = name, aggregation_mode = mode,
      estimate_difference = abs(a$estimate - b$estimate), se_difference = abs(a$se - b$se),
      max_weight_difference = max(abs(a$fold_weights - b$fold_weights)))
  }))
  write.csv(do.call(rbind, timings), file.path(output, "timing.csv"), row.names = FALSE)
  write.csv(do.call(rbind, comparisons), file.path(output, "result_comparisons.csv"), row.names = FALSE)
  write.csv(do.call(rbind, model_comparisons), file.path(output, "model_comparisons.csv"), row.names = FALSE)
}
results <- do.call(rbind, comparisons)
models <- do.call(rbind, model_comparisons)
passed <- all(models$lambda_equal) &&
  all(models$max_coefficient_difference <= configuration$coefficient_tolerance) &&
  all(results$estimate_difference <= configuration$estimate_tolerance) &&
  all(results$se_difference <= configuration$se_tolerance) &&
  all(results$max_weight_difference <= configuration$weight_tolerance)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
stopifnot(length(list.files(cache_parent, all.files = TRUE, no.. = TRUE)) == 0L)
writeLines(if (passed) "PASSED" else "FAILED", file.path(output, "status.txt"))
if (!passed) stop("Acceleration equivalence criteria failed; all comparisons retained.")
