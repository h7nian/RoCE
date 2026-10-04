# Site-specific inputs to the shared fold-summed calibration solver. Fold IDs
# always refer to the original partition, including in nested validation.

.validate_calibration_folds <- function(calibration_folds, n_folds, caller) {
  if (!is.numeric(calibration_folds) || length(calibration_folds) < 2L ||
      anyNA(calibration_folds) || any(!is.finite(calibration_folds)) ||
      any(calibration_folds != floor(calibration_folds)) ||
      any(!calibration_folds %in% seq_len(n_folds)) || anyDuplicated(calibration_folds)) {
    stop(caller, ": require at least two distinct original calibration fold IDs.", call. = FALSE)
  }
  sort(as.integer(calibration_folds))
}

.prepare_source_calibration_messages <- function(
    target_folds, source_folds, site, calibration_folds, communication_mode,
    A_val, family, M_tau, nlambda, lambda_rule, fit_cache = NULL,
    calibration_recipe = c("legacy", "score_derivative")) {
  calibration_recipe <- match.arg(calibration_recipe)
  calibration_folds <- .validate_calibration_folds(
    calibration_folds, length(target_folds), ".prepare_source_calibration_messages"
  )
  communication_mode <- match.arg(communication_mode, c("one_round", "two_round"))
  spec <- resolve_glm_family(family)
  messages <- vector("list", length(calibration_folds))
  names(messages) <- paste0("k2_", calibration_folds)
  for (fold in calibration_folds) {
    initial_training_folds <- setdiff(calibration_folds, fold)
    target_train <- combine_folds(target_folds, initial_training_folds)
    target_calibration <- materialize_fold(target_folds, fold)
    initial_site <- if (communication_mode == "one_round") "t" else site
    outcome_train <- if (initial_site == "t") target_train else
      combine_folds(source_folds, initial_training_folds)
    initial_outcome <- .fit_nuisance_training_subset("fit_initial_outcome", list(
      W_outcome = outcome_train$W_outcome, Y = outcome_train$Y, A = outcome_train$A,
      A_val = A_val, lambda = NULL, nlambda = nlambda, family = family,
      lambda_rule = lambda_rule, cv_group_id = outcome_train$cv_group_id
    ), initial_site, initial_training_folds, fit_cache)
    messages[[paste0("k2_", fold)]] <- list(
      mean_phi = c(1, colMeans(target_train$Z_site)),
      mean_grad_psi_init = .mean_glm_gradient_site_basis(
        target_calibration$W_outcome, target_calibration$Z_site, initial_outcome,
        spec$family_int, spec$link_int,
        if (calibration_recipe == "legacy") M_tau else Inf
      ),
      alpha_init = initial_outcome,
      initial_site = initial_site,
      calibration_recipe = calibration_recipe,
      initial_training_folds = initial_training_folds,
      calibration_fold = fold,
      excluded_folds = setdiff(seq_along(target_folds), calibration_folds),
      target_calibration_size = target_calibration$n
    )
  }
  messages
}

.fit_source_calibration <- function(
    site_folds, site, calibration_folds, get_fold_inputs, A_val,
    family_int, link_int, M_tau, nlambda, max_iter, lambda_rule,
    fit_cache = NULL, combine_cache = NULL, layout = c("block", "compact"),
    tol = TOL_DEFAULT, calibration_recipe = c("legacy", "score_derivative"),
    source_lambda_rules = NULL) {
  layout <- match.arg(layout)
  calibration_recipe <- match.arg(calibration_recipe)
  rules <- .resolve_source_lambda_rules(source_lambda_rules, lambda_rule)
  calibration_folds <- .validate_calibration_folds(
    calibration_folds, length(site_folds), ".fit_source_calibration"
  )
  keys <- paste0("k2_", calibration_folds)
  blocks <- initial_outcome <- initial_weight <- linear_moments <-
    setNames(vector("list", length(keys)), keys)
  timing <- c(initial_density_ratio = 0, initial_density_ratio_cv = 0,
              initial_density_ratio_final_fit = 0)
  for (fold in calibration_folds) {
    key <- paste0("k2_", fold)
    training_folds <- setdiff(calibration_folds, fold)
    inputs <- get_fold_inputs(site, fold)
    if (!is.null(inputs$calibration_recipe) && !identical(inputs$calibration_recipe, calibration_recipe)) {
      stop(".fit_source_calibration: message calibration recipe disagrees with the training task.", call. = FALSE)
    }
    if (!is.null(inputs$initial_training_folds) &&
        (!identical(as.integer(inputs$initial_training_folds), training_folds) ||
         !identical(as.integer(inputs$calibration_fold), fold))) {
      stop(".fit_source_calibration: message fold provenance disagrees with training task.",
           call. = FALSE)
    }
    cache_key <- paste(site, paste(training_folds, collapse = "_"), sep = ":")
    train <- if (is.null(combine_cache)) NULL else combine_cache[[cache_key]]
    if (is.null(train)) {
      train <- combine_folds(site_folds, training_folds)
      if (!is.null(combine_cache)) combine_cache[[cache_key]] <- train
    }
    started_at <- proc.time()[["elapsed"]]
    initial_weight[[key]] <- .fit_nuisance_training_subset("fit_initial_density_ratio", list(
      Z_site = train$Z_site, A = train$A, mean_phi = inputs$mean_phi,
      lambda = NULL, A_val = A_val, M_tau = M_tau, nlambda = nlambda,
      max_iter = max_iter, tol = tol, lambda_rule = rules$initial_weight, cv_group_id = train$cv_group_id
    ), site, training_folds, fit_cache)
    timing <- timing + c(
      initial_density_ratio = .elapsed_process_seconds(started_at),
      initial_density_ratio_cv = attr(initial_weight[[key]], "cv_seconds"),
      initial_density_ratio_final_fit = attr(initial_weight[[key]], "final_fit_seconds")
    )
    initial_outcome[[key]] <- inputs$alpha_init
    linear_moments[[key]] <- inputs$mean_grad_psi_init
    blocks[[key]] <- materialize_fold(site_folds, fold)
  }
  # Retain the manuscript's mean of per-fold target moments. This is a distinct
  # normalization from weighting observations in the stacked source loss.
  calibrated <- .fit_fold_summed_calibration(
    blocks, initial_outcome, initial_weight, .average_numeric_list(linear_moments),
    site, calibration_folds, A_val, family_int, link_int, M_tau,
    nlambda, max_iter, lambda_rule, layout = layout, tol = tol, fit_cache = fit_cache,
    calibration_recipe = calibration_recipe, weight_lambda_rule = rules$weight,
    outcome_lambda_rule = rules$outcome
  )
  calibrated$initial_outcome <- initial_outcome
  calibrated$initial_weight <- initial_weight
  calibrated$timing <- c(timing, calibrated$timing)
  calibrated
}

# Target PS initialization is either a full-target logistic likelihood or the
# unweighted clipped-odds calibration loss. The logistic path rejects the
# outcome-only degenerate fallback: that constant fallback is not a PS fit.
.fit_initial_target_propensity <- function(
    Z_site, A, A_val, nlambda, lambda_rule, cv_group_id = NULL,
    initialization = c("logistic", "calibrated"), M_tau = M_TAU_DEFAULT,
    tol = TOL_DEFAULT) {
  initialization <- match.arg(initialization)
  A_val <- .validate_A_val(A_val, ".fit_initial_target_propensity")
  if (length(unique(A)) != 2L || any(!A %in% c(0, 1))) {
    stop(".fit_initial_target_propensity: both treatment classes are required.", call. = FALSE)
  }
  if (initialization == "calibrated") {
    opposite_moment <- colMeans(cbind(1, Z_site) * as.numeric(A != A_val))
    return(fit_initial_density_ratio(Z_site, A, opposite_moment,
      A_val = A_val, nlambda = nlambda, lambda_rule = lambda_rule,
      cv_group_id = cv_group_id, M_tau = M_tau, tol = tol))
  }
  fit <- fit_initial_outcome(Z_site, as.numeric(A == A_val), rep(1L, length(A)),
    A_val = 1L, nlambda = nlambda, family = "binomial",
    lambda_rule = lambda_rule, cv_group_id = cv_group_id)
  if (identical(attr(fit, "lambda_rule"), "degenerate_constant")) {
    stop(".fit_initial_target_propensity: insufficient treatment variation for PS tuning.",
         call. = FALSE)
  }
  fit
}

.fit_target_calibration <- function(
    target_folds, calibration_folds, A_val, family, M_tau, nlambda,
    max_iter = MAX_ITER_DEFAULT, lambda_rule = "min", fit_cache = NULL,
    layout = c("block", "compact"), tol = TOL_DEFAULT, calibration_control = NULL) {
  layout <- match.arg(layout)
  calibration_control <- .validate_calibration_control(calibration_control)
  M_tau <- calibration_control$target_radius %||% M_tau
  calibration_folds <- .validate_calibration_folds(
    calibration_folds, length(target_folds), ".fit_target_calibration"
  )
  spec <- resolve_glm_family(family)
  keys <- paste0("k2_", calibration_folds)
  blocks <- initial_outcome <- initial_weight <-
    setNames(vector("list", length(keys)), keys)
  linear_sum <- NULL
  total_size <- 0
  for (fold in calibration_folds) {
    key <- paste0("k2_", fold)
    training_folds <- setdiff(calibration_folds, fold)
    train <- combine_folds(target_folds, training_folds)
    block <- materialize_fold(target_folds, fold)
    if (!identical(unname(as.matrix(train$W_outcome)), unname(as.matrix(train$Z_site))) ||
        !identical(unname(as.matrix(block$W_outcome)), unname(as.matrix(block$Z_site)))) {
      stop(".fit_target_calibration: Hou calibration currently requires the same OR and PS basis.",
           call. = FALSE)
    }
    initial_outcome[[key]] <- .fit_nuisance_training_subset("fit_initial_outcome", list(
      W_outcome = train$W_outcome, Y = train$Y, A = train$A, A_val = A_val,
      lambda = NULL, nlambda = nlambda, family = family,
      lambda_rule = lambda_rule, cv_group_id = train$cv_group_id
    ), "t", training_folds, fit_cache)
    initial_weight[[key]] <- .fit_nuisance_training_subset(".fit_initial_target_propensity", list(
      Z_site = train$Z_site, A = train$A, A_val = A_val, nlambda = nlambda,
      lambda_rule = lambda_rule, cv_group_id = train$cv_group_id,
      initialization = calibration_control$target_propensity_initialization,
      M_tau = M_tau, tol = tol
    ), "t", training_folds, fit_cache)
    opposite_arm <- which(block$A != A_val)
    if (!length(opposite_arm)) {
      stop(".fit_target_calibration: each calibration fold must contain both arms.", call. = FALSE)
    }
    # Full-target empirical E[(1-I_a) m'_initial(X) Z]. The source side
    # E[I_a m'_initial(X) exp(-Z theta)] is assembled by the shared solver.
    moment <- length(opposite_arm) * .mean_glm_gradient_site_basis(
      block$W_outcome[opposite_arm, , drop = FALSE],
      block$Z_site[opposite_arm, , drop = FALSE], initial_outcome[[key]],
      spec$family_int, spec$link_int,
      if (calibration_control$recipe == "legacy") M_tau else Inf
    )
    linear_sum <- if (is.null(linear_sum)) moment else linear_sum + moment
    total_size <- total_size + block$n
    blocks[[key]] <- block
  }
  calibrated <- .fit_fold_summed_calibration(
    blocks, initial_outcome, initial_weight, linear_sum / total_size,
    "t", calibration_folds, A_val, spec$family_int, spec$link_int, M_tau,
    nlambda, max_iter, lambda_rule, layout = layout, tol = tol, fit_cache = fit_cache,
    calibration_recipe = calibration_control$recipe
  )
  calibrated$initial_outcome <- initial_outcome
  calibrated$initial_weight <- initial_weight
  calibrated$training_radius <- M_tau
  calibrated$propensity_initialization <- calibration_control$target_propensity_initialization
  calibrated
}

.evaluate_target_calibration <- function(
    calibrated, evaluation_data, A_val, family, M_tau_inference) {
  spec <- resolve_glm_family(family)
  predictor <- drop(cbind(1, evaluation_data$Z_site) %*% calibrated$weight)
  arm_propensity <- stats::plogis(pmax(-M_tau_inference, pmin(M_tau_inference, predictor)))
  outcome_prediction <- predict_glm_cpp(
    evaluation_data$W_outcome, calibrated$outcome, spec$family_int, spec$link_int
  )
  score <- outcome_prediction + as.numeric(evaluation_data$A == A_val) *
    (evaluation_data$Y - outcome_prediction) / arm_propensity
  if (any(!is.finite(score))) {
    stop(".evaluate_target_calibration: nonfinite target score.", call. = FALSE)
  }
  estimate <- mean(score)
  centered <- as.numeric(score - estimate)
  list(
    estimate = estimate, V_ot = mean(centered^2),
    variance = mean(centered^2) / evaluation_data$n, varphi_ot = centered,
    # Legacy prop_scores is the arm-implied P(A=1), not a shared fit across arms.
    prop_scores = if (A_val == 1L) arm_propensity else 1 - arm_propensity,
    arm_propensity = arm_propensity, m_pred = outcome_prediction,
    outcome_degenerate = sum(vapply(calibrated$initial_outcome, function(fit)
      identical(attr(fit, "lambda_rule"), "degenerate_constant"), logical(1L))),
    method = "hou_calibrated", calibration = calibrated,
    inference_scope = "empirical_score; clipping and nuisance remainder require validation"
  )
}
