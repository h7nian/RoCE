# Complete source training programs on explicitly allowed original folds.
# The standard ablation retains the initial merged-weight balancing loss and
# replaces the final score-calibration updates with ordinary outcome fitting.

.complete_source_validation <- function(method) {
  method %in% c("calibrated", "complete")
}

.match_source_validation <- function(method, calibration_control, communication_mode = NULL) {
  if (identical(method, c("initial", "calibrated")) ||
      identical(method, c("initial", "calibrated", "complete"))) method <- method[1L]
  method <- match.arg(method, c("initial", "calibrated", "complete", "outer_fit"))
  if (calibration_control$source_nuisance_method == "standard") {
    if (method == "initial") {
      stop("Standard source fitting requires complete-program weight learning (complete or outer_fit).", call. = FALSE)
    }
    if (!is.null(communication_mode) && communication_mode != "one_round") {
      stop("The standard source ablation currently uses communication_mode='one_round'.", call. = FALSE)
    }
    if (method != "outer_fit") method <- "complete"
  }
  method
}

.source_inner_fits <- function(source_result) {
  source_result$inner_fits %||% source_result$inner_calibrated
}

.prepare_source_training_message <- function(target_folds, training_folds) {
  training_folds <- .validate_calibration_folds(training_folds, length(target_folds),
                                              ".prepare_source_training_message")
  target <- combine_folds(target_folds, training_folds)
  list(mean_phi = c(1, colMeans(target$Z_site)), training_folds = training_folds,
       excluded_folds = setdiff(seq_along(target_folds), training_folds), target_size = target$n)
}

.fit_standard_source_nuisances <- function(site_folds, site, training_folds, target_message,
                                            A_val, family_int, link_int, M_tau, nlambda,
                                            max_iter, lambda_rule, fit_cache = NULL,
                                            tol = TOL_DEFAULT, weight_lambda = NULL,
                                            outcome_lambda = NULL) {
  training_folds <- .validate_calibration_folds(training_folds, length(site_folds),
                                              ".fit_standard_source_nuisances")
  if (!identical(as.integer(target_message$training_folds), training_folds)) {
    stop("Standard source training and target message use different folds.", call. = FALSE)
  }
  train <- combine_folds(site_folds, training_folds)
  weight <- .fit_nuisance_training_subset("fit_initial_density_ratio", list(
    Z_site = train$Z_site, A = train$A, mean_phi = target_message$mean_phi,
    lambda = weight_lambda, A_val = A_val, M_tau = M_tau, nlambda = nlambda,
    max_iter = max_iter, tol = tol, lambda_rule = lambda_rule, cv_group_id = train$cv_group_id
  ), site, training_folds, fit_cache, cv_seed_fitter = "fit_unified_density_ratio")
  # A zero log tilt gives exactly unit outcome-loss weights. This uses the
  # same solver, full-site normalization and CV grid machinery as calibration.
  outcome <- .fit_nuisance_training_subset("fit_unified_outcome", list(
    W_outcome = train$W_outcome, Y = train$Y, A = train$A, A_val = A_val,
    gamma_s = rep(0, ncol(train$Z_site) + 1L), Z_site = train$Z_site,
    lambda = outcome_lambda, calibrated = FALSE, M_tau = M_tau,
    family_int = family_int, link_int = link_int, nlambda = nlambda,
    max_iter = max_iter, tol = tol, lambda_rule = lambda_rule, cv_group_id = train$cv_group_id
  ), site, training_folds, fit_cache)
  timing <- c(initial_density_ratio = attr(weight, "cv_seconds") + attr(weight, "final_fit_seconds"),
    initial_density_ratio_cv = attr(weight, "cv_seconds"),
    initial_density_ratio_final_fit = attr(weight, "final_fit_seconds"),
    calibrated_density_ratio = 0, calibrated_density_ratio_cv = 0,
    calibrated_density_ratio_final_fit = 0, calibrated_outcome = 0,
    calibrated_outcome_cv = 0, calibrated_outcome_final_fit = 0,
    standard_outcome = attr(outcome, "cv_seconds") + attr(outcome, "final_fit_seconds"))
  list(weight = weight, outcome = outcome,
    initial_weight = list(training = weight), initial_outcome = list(training = outcome),
    training_folds = training_folds, calibration_folds = integer(),
    layout = "not_applicable", calibration_recipe = "none", timing = timing,
    source_nuisance_method = "standard")
}

.fit_complete_source_program <- function(target_folds, source_folds, site, training_folds,
                                           communication_mode, A_val, family, M_tau, nlambda,
                                           max_iter, lambda_rule, fit_cache = NULL,
                                           layout = "compact", tol = TOL_DEFAULT,
                                           calibration_control = NULL) {
  control <- .validate_calibration_control(calibration_control)
  spec <- resolve_glm_family(family)
  if (control$source_nuisance_method == "standard") {
    return(.fit_standard_source_nuisances(source_folds, site, training_folds,
      .prepare_source_training_message(target_folds, training_folds), A_val,
      spec$family_int, spec$link_int, M_tau, nlambda, max_iter, lambda_rule, fit_cache, tol))
  }
  messages <- .prepare_source_calibration_messages(target_folds, source_folds, site, training_folds,
    communication_mode, A_val, family, M_tau, nlambda, lambda_rule, fit_cache,
    calibration_recipe = control$recipe)
  fit <- .fit_source_calibration(source_folds, site, training_folds,
    function(.site, fold) messages[[paste0("k2_", fold)]], A_val,
    spec$family_int, spec$link_int, M_tau, nlambda, max_iter, lambda_rule, fit_cache,
    layout = layout, tol = tol, calibration_recipe = control$recipe,
    source_lambda_rules = control$source_lambda_rules)
  fit$training_folds <- training_folds
  fit$source_nuisance_method <- "calibrated"
  fit
}
