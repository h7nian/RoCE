# Shared loss assembly for source transport and target propensity calibration.
# Adapters supply the linear moment and each block's initial models; this layer
# only assembles and solves their common, full-sample normalized objectives.

.calibration_plugin <- function(designs, coefficients, layout) {
  layout <- match.arg(layout, c("compact", "block"))
  if (!length(designs) || length(designs) != length(coefficients)) {
    stop(".calibration_plugin: each design requires one initial model.", call. = FALSE)
  }
  designs <- lapply(designs, as.matrix)
  for (index in seq_along(designs)) {
    design <- as.matrix(designs[[index]])
    coefficient <- as.numeric(coefficients[[index]])
    if (length(coefficient) != ncol(design) + 1L ||
        any(!is.finite(design)) || any(!is.finite(coefficient))) {
      stop(".calibration_plugin: initial model and design must be finite and conformable.",
           call. = FALSE)
    }
  }
  if (layout == "compact") {
    predictors <- Map(function(design, coefficient) {
      drop(cbind(1, design) %*% as.numeric(coefficient))
    }, designs, coefficients)
    list(design = matrix(unlist(predictors, use.names = FALSE), ncol = 1L),
         coefficients = c(0, 1))
  } else {
    list(design = .make_plugin_block_design(designs),
         coefficients = c(0, unlist(lapply(coefficients, as.numeric), use.names = FALSE)))
  }
}

.fit_fold_summed_calibration <- function(
    calibration_data, initial_outcome, initial_weight, linear_moment,
    site, training_folds, A_val, family_int, link_int, M_tau,
    nlambda, max_iter, lambda_rule, layout = c("compact", "block"),
    fit_cache = NULL, weight_lambda = NULL, outcome_lambda = NULL,
    tol = TOL_DEFAULT, calibration_recipe = c("legacy", "score_derivative"),
    weight_lambda_rule = lambda_rule, outcome_lambda_rule = lambda_rule) {
  layout <- match.arg(layout)
  weight_lambda_rule <- .match_nuisance_lambda_rule(weight_lambda_rule, ".fit_fold_summed_calibration")
  outcome_lambda_rule <- .match_nuisance_lambda_rule(outcome_lambda_rule, ".fit_fold_summed_calibration")
  calibration_recipe <- match.arg(calibration_recipe)
  if (calibration_recipe == "score_derivative") .validate_score_calibration_radius(M_tau)
  count <- length(calibration_data)
  if (count < 2L || length(initial_outcome) != count || length(initial_weight) != count ||
      length(training_folds) != count) {
    stop(".fit_fold_summed_calibration: require at least two aligned calibration blocks.",
         call. = FALSE)
  }
  for (block in calibration_data) {
    if (calibration_recipe == "score_derivative" &&
        !identical(unname(as.matrix(block$W_outcome)), unname(as.matrix(block$Z_site)))) {
      stop("Score-derivative calibration currently requires the same outcome and weight basis.", call. = FALSE)
    }
    if (!is.numeric(block$n) || length(block$n) != 1L || !is.finite(block$n) ||
        block$n < 1L || block$n != floor(block$n)) {
      stop(".fit_fold_summed_calibration: block n must be a positive integer.", call. = FALSE)
    }
    sizes <- c(nrow(block$W_outcome), nrow(block$Z_site), length(block$A), length(block$Y))
    if (length(sizes) != 4L || any(sizes != block$n) || block$n < 1L) {
      stop(".fit_fold_summed_calibration: calibration block row counts disagree.", call. = FALSE)
    }
  }
  if (!identical(names(calibration_data), names(initial_outcome)) ||
      !identical(names(calibration_data), names(initial_weight))) {
    stop(".fit_fold_summed_calibration: initial models must follow the named calibration blocks.",
         call. = FALSE)
  }
  caller <- ".fit_fold_summed_calibration"
  outcome_design <- .stack_fold_field(calibration_data, "W_outcome", caller)
  weight_design <- .stack_fold_field(calibration_data, "Z_site", caller)
  treatment <- .stack_fold_field(calibration_data, "A", caller)
  outcome <- .stack_fold_field(calibration_data, "Y", caller)
  grouped <- vapply(calibration_data, function(block) !is.null(block$cv_group_id), logical(1L))
  if (any(grouped) && !all(grouped)) {
    stop(caller, ": group metadata must be present in every calibration block or none.", call. = FALSE)
  }
  cv_group_id <- if (all(grouped)) .stack_fold_field(calibration_data, "cv_group_id", caller) else NULL
  outcome_plugin <- .calibration_plugin(lapply(calibration_data, `[[`, "W_outcome"),
                                         initial_outcome, layout)
  weight_plugin <- .calibration_plugin(lapply(calibration_data, `[[`, "Z_site"),
                                        initial_weight, layout)
  weight_started_at <- proc.time()[["elapsed"]]
  weight_coefficients <- .fit_nuisance_training_subset("fit_unified_density_ratio", list(
    Z_site = weight_design, A = treatment, mean_grad_psi = linear_moment,
    alpha_init = outcome_plugin$coefficients, W_outcome = outcome_plugin$design,
    lambda = weight_lambda, calibrated = TRUE, M_tau = M_tau, A_val = A_val,
    truncate_initial_outcome = calibration_recipe == "legacy",
    family_int = family_int, link_int = link_int,
    warm_start = as.numeric(.average_numeric_list(initial_weight)),
    nlambda = nlambda, max_iter = max_iter, tol = tol, lambda_rule = weight_lambda_rule,
    cv_group_id = cv_group_id
  ), site = site, training_folds = training_folds, fit_cache = fit_cache)
  weight_seconds <- .elapsed_process_seconds(weight_started_at)
  outcome_started_at <- proc.time()[["elapsed"]]
  outcome_coefficients <- .fit_nuisance_training_subset("fit_unified_outcome", list(
    W_outcome = outcome_design, Y = outcome, A = treatment, A_val = A_val,
    gamma_s = weight_plugin$coefficients, Z_site = weight_plugin$design,
    lambda = outcome_lambda, calibrated = TRUE, M_tau = M_tau,
    use_weight_derivative = calibration_recipe == "score_derivative",
    family_int = family_int, link_int = link_int,
    warm_start = as.numeric(.average_numeric_list(initial_outcome)),
    nlambda = nlambda, max_iter = max_iter, tol = tol, lambda_rule = outcome_lambda_rule,
    cv_group_id = cv_group_id
  ), site = site, training_folds = training_folds, fit_cache = fit_cache)
  list(
    outcome = outcome_coefficients, weight = weight_coefficients,
    calibration_recipe = calibration_recipe,
    layout = layout, calibration_folds = training_folds,
    calibration_sizes = vapply(calibration_data, `[[`, numeric(1L), "n"),
    timing = c(
      calibrated_density_ratio = weight_seconds,
      calibrated_density_ratio_cv = attr(weight_coefficients, "cv_seconds"),
      calibrated_density_ratio_final_fit = attr(weight_coefficients, "final_fit_seconds"),
      calibrated_outcome = .elapsed_process_seconds(outcome_started_at),
      calibrated_outcome_cv = attr(outcome_coefficients, "cv_seconds"),
      calibrated_outcome_final_fit = attr(outcome_coefficients, "final_fit_seconds")
    )
  )
}
