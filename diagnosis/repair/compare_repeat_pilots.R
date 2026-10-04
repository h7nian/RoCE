#!/usr/bin/env Rscript
# Compare completed recipe outputs from two frozen pilot versions. Missing or
# failed runs remain explicit; this gate never substitutes another recipe.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
  stop("usage: compare_repeat_pilots.R LIBRARY REFERENCE_ROOT REFERENCE_TASK CANDIDATE_ROOT CANDIDATE_TASK RECIPE NEW_OUTPUT")
}
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
driver_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(driver_file)), "acceleration_helpers.R"))
output <- args[7L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
options(error = function() {
  writeLines("FAILED", file.path(output, "status.txt"))
  traceback(15L)
  quit(status = 1L)
})
recipe <- match.arg(args[6L], c("legacy", "score_derivative"))
inputs <- list(reference = list(root = args[2L], task = args[3L]),
               candidate = list(root = args[4L], task = args[5L]))
thresholds <- c(lambda = 1e-10, coefficient = 1e-4, estimate = 1e-6,
                se = 1e-6, weight = 1e-5)
saveRDS(list(inputs = inputs, recipe = recipe, thresholds = thresholds),
        file.path(output, "configuration.rds"))
configuration <- tasks <- results <- list()
states <- lapply(names(inputs), function(name) {
  input <- inputs[[name]]
  configuration[[name]] <<- jsonlite::fromJSON(file.path(input$root, "configuration.json"))
  manifest <- read.csv(file.path(input$root, "manifest.csv"), stringsAsFactors = FALSE)
  tasks[[name]] <<- manifest[manifest$task_id == as.integer(input$task), , drop = FALSE]
  stopifnot(nrow(tasks[[name]]) == 1L, recipe %in% configuration[[name]]$recipes)
  rownames(tasks[[name]]) <<- NULL
  directory <- file.path(input$root, "tasks", input$task)
  path <- file.path(directory, paste0(recipe, ".rds"))
  state_path <- file.path(directory, "status.txt")
  state <- if (file.exists(state_path)) paste(readLines(state_path), collapse = "; ") else "NOT_STARTED"
  if (file.exists(path)) results[[name]] <<- readRDS(path)
  error <- if (inherits(results[[name]], "error")) conditionMessage(results[[name]]) else ""
  data.frame(version = name, task_state = state, result_exists = file.exists(path),
             recipe_error = error, path = path)
})
states <- do.call(rbind, states)
write.csv(states, file.path(output, "input_status.csv"), row.names = FALSE)
if (!all(states$result_exists) || any(nzchar(states$recipe_error))) {
  writeLines(c("# Paired acceleration check: incomplete", "",
    "A requested recipe is missing or failed. No numerical equivalence claim is made.",
    "See input_status.csv; inspect the parent job logs and Slurm state."),
    file.path(output, "report.md"))
  writeLines("INCOMPLETE", file.path(output, "status.txt"))
  quit(status = 2L)
}
scientific_fields <- c("n_per_site", "p", "n_folds", "nlambda", "nuisance_tol",
                       "source_radius", "target_score_radius", "aggregation_lambda")
stopifnot(all(scientific_fields %in% names(configuration$reference)),
          all(scientific_fields %in% names(configuration$candidate)),
          identical(configuration$reference[scientific_fields], configuration$candidate[scientific_fields]),
          identical(tasks$reference[names(tasks$reference) != "task_id"],
                    tasks$candidate[names(tasks$candidate) != "task_id"]))
artifacts <- lapply(results, attr, which = "roce_simulation_artifacts")
stopifnot(!is.null(artifacts$reference$data_split),
          identical(artifacts$reference$data_split, artifacts$candidate$data_split),
          identical(artifacts$reference$tate_truth, artifacts$candidate$tate_truth))
protocol <- tasks$reference$protocol
fits <- lapply(artifacts, function(x) x$direct_tate_results[[paste0(protocol, "_crossfit")]])
stopifnot(!is.null(fits$reference), !is.null(fits$candidate))
fit_fields <- c("crossfit_levels", "communication_mode", "nuisance_tol", "calibration_control",
                "nuisance_training_policy", "nuisance_lambda_rule", "family", "M_tau", "M_tau_inference")
stopifnot(identical(fits$reference[fit_fields], fits$candidate[fit_fields]))
models <- compare_calibration_models(fits$reference, fits$candidate, thresholds[["lambda"]])
write.csv(models, file.path(output, "models.csv"), row.names = FALSE)

rows <- lapply(results, function(x) {
  x <- x[x$estimand_scope == "tate", c("method", "estimate", "se", "truth")]
  stopifnot(!anyDuplicated(x$method), nrow(x) > 0L)
  x[order(x$method), , drop = FALSE]
})
stopifnot(identical(rows$reference$method, rows$candidate$method),
          identical(rows$reference$truth, rows$candidate$truth))
estimates <- data.frame(method = rows$reference$method,
  estimate_difference = abs(rows$candidate$estimate - rows$reference$estimate),
  se_difference = abs(rows$candidate$se - rows$reference$se))
write.csv(estimates, file.path(output, "estimates.csv"), row.names = FALSE)

# Recompute only aggregation with the same checked implementation. Also compare
# each reconstruction to its independently saved estimate/SE row.
comparisons <- list()
for (mode in c("common_tate", "separate_arms", "joint_tate")) {
  aggregated <- lapply(fits, function(fit) calculate_tate_crossfit_aggregation(
    artifacts$reference$data_split, fit$arm_results$mu1, fit$arm_results$mu0,
    lambda_selection = configuration$reference$aggregation_lambda,
    lambda_rule = fit$aggregation_lambda_rule, screening_rule = fit$aggregation_screening_rule,
    aggregation_mode = mode, verbose = FALSE))
  method <- paste0(protocol, "_crossfit_ate", if (mode == "common_tate") "" else paste0("_", mode))
  reconstruction_error <- vapply(names(aggregated), function(name) {
    saved <- rows[[name]][rows[[name]]$method == method, , drop = FALSE]
    stopifnot(nrow(saved) == 1L)
    max(abs(c(aggregated[[name]]$estimate - saved$estimate, aggregated[[name]]$se - saved$se)))
  }, numeric(1L))
  a <- aggregated$candidate
  b <- aggregated$reference
  stopifnot(identical(dim(a$fold_weights), dim(b$fold_weights)),
            identical(dimnames(a$fold_weights), dimnames(b$fold_weights)))
  comparisons[[mode]] <- data.frame(aggregation_mode = mode,
    estimate_difference = abs(a$estimate - b$estimate), se_difference = abs(a$se - b$se),
    max_weight_difference = max(abs(a$fold_weights - b$fold_weights)),
    max_reconstruction_error = max(reconstruction_error))
}
comparisons <- do.call(rbind, comparisons)
write.csv(comparisons, file.path(output, "aggregation.csv"), row.names = FALSE)
passed <- all(models$lambda_equal) &&
  all(models$max_coefficient_difference <= thresholds[["coefficient"]]) &&
  all(is.finite(as.matrix(estimates[, -1L]))) &&
  all(estimates$estimate_difference <= thresholds[["estimate"]]) &&
  all(estimates$se_difference <= thresholds[["se"]]) &&
  all(is.finite(as.matrix(comparisons[, -1L]))) &&
  all(comparisons$estimate_difference <= thresholds[["estimate"]]) &&
  all(comparisons$se_difference <= thresholds[["se"]]) &&
  all(comparisons$max_weight_difference <= thresholds[["weight"]]) &&
  all(comparisons$max_reconstruction_error <= min(thresholds[c("estimate", "se")]))
writeLines(c("# Paired acceleration check", "",
  paste("Result:", if (passed) "PASSED" else "FAILED"),
  paste("Recipe:", recipe), paste("Compared nuisance models:", nrow(models)), "",
  "Identical data, scenario/seed and statistical configuration were required.",
  "All three aggregation modes and saved TATE rows were compared.",
  "This is one numerical pairing, not coverage or RMSE evidence."), file.path(output, "report.md"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines(if (passed) "PASSED" else "FAILED", file.path(output, "status.txt"))
if (!passed) stop("Paired numerical gates failed; detailed differences retained.")
