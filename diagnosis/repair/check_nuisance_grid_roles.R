#!/usr/bin/env Rscript
# Component-wise grid comparison. Final calibration keeps its initial models
# fixed, so changing a grid is not confounded with changed plug-in predictions.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: check_nuisance_grid_roles.R LIBRARY INPUT_RDS CALIBRATION_REFERENCE NEW_OUTPUT")
}
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[4L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
saved <- readRDS(args[2L])
folds <- saved$folds
reference_source <- readRDS(file.path(args[3L], "source_fit.rds"))
reference_target <- readRDS(file.path(args[3L], "target_fit.rds"))
messages <- readRDS(file.path(args[3L], "target_messages.rds"))
target_train <- RoCE:::combine_folds(folds$target_folds, 4:10)
source_train <- RoCE:::combine_folds(folds$source_folds$s1, 4:10)
invisible(RoCE:::set_nuisance_solver_cpp("proximal_newton"))
writeLines(c("1000/site, raw p100, 200-column basis, ten original folds; C3 fixture seed1081",
  "One treated-arm source/target task; same data and per-task CV seed for all grids",
  "Final calibration comparisons hold initial coefficients/predictions fixed at the 100-point reference",
  "Not an end-to-end TATE or coverage study"), file.path(output, "scope.txt"))

profile_state <- new.env(parent = emptyenv())
profile_state$events <- list()
profile_state$role <- ""
profile_state$points <- NA_integer_
assign(".roce_grid_profile", profile_state, globalenv())
trace("glmnet", where = asNamespace("glmnet"), print = FALSE,
  tracer = quote(.roce_grid_started <- proc.time()[["elapsed"]]),
  exit = quote({
    .roce_grid_value <- returnValue()
    .roce_grid_state <- get(".roce_grid_profile", globalenv())
    if (is.list(.roce_grid_value) && !is.null(.roce_grid_value$lambda)) {
      .roce_grid_state$events[[length(.roce_grid_state$events) + 1L]] <- data.frame(
        role = .roce_grid_state$role, requested_grid = .roce_grid_state$points,
        n_train = nrow(x), p = ncol(x), actual_lambda_count = length(.roce_grid_value$lambda),
        path_max = max(.roce_grid_value$lambda), path_min = min(.roce_grid_value$lambda),
        passes = .roce_grid_value$npasses, error_code = .roce_grid_value$jerr,
        seconds = proc.time()[["elapsed"]] - .roce_grid_started)
    }
  }))
fits <- rows <- list()
add_result <- function(role, points, fit, seconds, estimate = NA_real_) {
  key <- paste(role, points, sep = "_")
  fits[[key]] <<- fit
  baseline <- fits[[paste(role, 100L, sep = "_")]]
  scalar_attr <- function(name) {
    value <- attr(fit, name)
    if (is.null(value)) NA_real_ else as.numeric(value)
  }
  rows[[key]] <<- data.frame(role = role, points = points, seconds = seconds,
    selected_lambda = scalar_attr("lambda_used"), cv_seconds = scalar_attr("cv_seconds"),
    final_fit_seconds = scalar_attr("final_fit_seconds"),
    invalid_lambdas = scalar_attr("cv_invalid_lambdas"),
    skipped_fold_fits = scalar_attr("cv_path_tail_skipped_fold_fits"),
    max_coefficient_difference = max(abs(as.numeric(fit) - as.numeric(baseline))),
    validation_mu1 = estimate)
  saveRDS(fit, file.path(output, paste0(key, ".rds")))
  write.csv(do.call(rbind, rows), file.path(output, "comparison.csv"), row.names = FALSE)
}

for (points in c(100L, 50L, 25L)) {
  initial_tasks <- list(
    target_initial_or = list(fitter = "fit_initial_outcome", site = "t", arguments = list(
      W_outcome = target_train$W_outcome, Y = target_train$Y, A = target_train$A, A_val = 1L,
      nlambda = points, family = "binomial", lambda_rule = "min")),
    target_initial_ps_logistic = list(fitter = ".fit_initial_target_propensity", site = "t", arguments = list(
      Z_site = target_train$Z_site, A = target_train$A, A_val = 1L, nlambda = points,
      lambda_rule = "min", initialization = "logistic", M_tau = 5, tol = 1e-10)),
    target_initial_ps_calibrated = list(fitter = ".fit_initial_target_propensity", site = "t", arguments = list(
      Z_site = target_train$Z_site, A = target_train$A, A_val = 1L, nlambda = points,
      lambda_rule = "min", initialization = "calibrated", M_tau = log(9), tol = 1e-10)),
    source_initial_weight = list(fitter = "fit_initial_density_ratio", site = "s1", arguments = list(
      Z_site = source_train$Z_site, A = source_train$A, A_val = 1L,
      mean_phi = c(1, colMeans(target_train$Z_site)), nlambda = points, M_tau = 5, tol = 1e-10)))
  for (role in names(initial_tasks)) {
    task <- initial_tasks[[role]]
    profile_state$role <- role
    profile_state$points <- points
    started <- proc.time()[["elapsed"]]
    fit <- RoCE:::.fit_nuisance_training_subset(task$fitter, task$arguments, task$site, 4:10)
    add_result(role, points, fit, proc.time()[["elapsed"]] - started)
    cat(role, points, "finished\n")
  }
}
untrace("glmnet", where = asNamespace("glmnet"))
write.csv(do.call(rbind, profile_state$events), file.path(output, "glmnet_calls.csv"), row.names = FALSE)
# Independently check that tracing did not alter the initial-model fit.
untraced <- RoCE:::.fit_nuisance_training_subset("fit_initial_outcome", list(
  W_outcome = target_train$W_outcome, Y = target_train$Y, A = target_train$A, A_val = 1L,
  nlambda = 100L, family = "binomial", lambda_rule = "min"), "t", 4:10)
stopifnot(identical(as.numeric(untraced), as.numeric(fits$target_initial_or_100)))

for (site in c("s1", "t")) {
  site_folds <- if (site == "t") folds$target_folds else folds$source_folds$s1
  reference <- if (site == "t") reference_target else reference_source
  blocks <- setNames(lapply(3:10, function(fold) RoCE:::materialize_fold(site_folds, fold)),
                     paste0("k2_", 3:10))
  radius <- if (site == "t") log(9) else 5
  if (site == "t") {
    moments <- Map(function(block, beta) {
      opposite <- block$A == 0L
      sum(opposite) * RoCE:::.mean_glm_gradient_site_basis(
        block$W_outcome[opposite, , drop = FALSE], block$Z_site[opposite, , drop = FALSE],
        beta, 1L, 1L, Inf)
    }, blocks, reference$initial_outcome)
    linear <- Reduce(`+`, moments) / sum(vapply(blocks, `[[`, numeric(1L), "n"))
  } else {
    linear <- RoCE:::.average_numeric_list(lapply(messages, `[[`, "mean_grad_psi_init"))
  }
  for (points in c(100L, 50L, 25L)) {
    started <- proc.time()[["elapsed"]]
    fit <- RoCE:::.fit_fold_summed_calibration(blocks, reference$initial_outcome,
      reference$initial_weight, linear, site, 3:10, 1L, 1L, 1L, radius, points,
      10000L, "min", layout = "compact", tol = 1e-10, calibration_recipe = "score_derivative")
    elapsed <- proc.time()[["elapsed"]] - started
    target_fold <- RoCE:::materialize_fold(folds$target_folds, 2L)
    evaluated <- if (site == "t") {
      RoCE:::.evaluate_target_calibration(fit, target_fold, 1L, "binomial", radius)
    } else RoCE:::.compute_inner_source_arm_info(RoCE:::materialize_fold(site_folds, 2L),
      target_fold, fit$weight, fit$outcome, 5, 1L, 1L, 1L)
    for (type in c("weight", "outcome")) {
      model <- fit[[type]]
      add_result(paste(site, "calibrated", type, sep = "_"), points, model,
        attr(model, "cv_seconds") + attr(model, "final_fit_seconds"), evaluated$estimate)
    }
    cat(site, "calibration", points, "finished in", elapsed, "seconds\n")
  }
}
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("COMPLETED; component grid study only", file.path(output, "status.txt"))
