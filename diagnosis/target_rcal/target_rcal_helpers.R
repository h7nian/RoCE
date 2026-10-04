# Diagnostic-only target nuisance alternatives. The installed RoCE namespace and
# production fits are never modified. All training excludes the evaluation fold.

target_truth_predictions <- function(data) {
  stopifnot(identical(data$dgp_type, "face"),
            identical(data$outcome_type, "binary"))
  target_rows <- data$R == "t"
  strengths <- RoCE:::.face_misspecification_strengths(
    data$config, data$misspecification_strength
  )
  calibration <- RoCE:::get_face_binary_calibration(
    p = data$p, config = data$config, kappa = RoCE:::FACE_KAPPA,
    misspecification_strength = strengths$outcome
  )
  predictor <- RoCE:::.face_mixed_predictor(
    data$X, data$X_dagger, data$out_params$beta_linear,
    data$out_params$beta_squared, RoCE:::FACE_KAPPA, strengths$outcome
  )
  logit <- RoCE:::face_binary_logit(predictor, calibration)[target_rows]
  list(mu1 = plogis(logit + unname(data$ate_map[["t"]])),
       mu0 = plogis(logit), propensity = data$p_treat_true[target_rows])
}

target_prediction_fit <- function(eval_data, A_val, outcome, arm_probability) {
  if (length(A_val) != 1L || !(A_val %in% c(0L, 1L)) ||
      length(outcome) != eval_data$n ||
      length(arm_probability) != eval_data$n ||
      any(!is.finite(outcome)) || any(!is.finite(arm_probability)) ||
      any(arm_probability <= 0 | arm_probability >= 1)) {
    stop("target_prediction_fit: invalid arm or prediction vectors.")
  }
  pseudo_outcome <- outcome + (eval_data$A == A_val) *
    (eval_data$Y - outcome) / arm_probability
  estimate <- mean(pseudo_outcome)
  influence <- pseudo_outcome - estimate
  list(estimate = estimate, V_ot = mean(influence^2),
       variance = mean(influence^2) / eval_data$n,
       varphi_ot = influence,
       prop_scores = if (A_val == 1L) arm_probability else 1 - arm_probability,
       m_pred = outcome, outcome_degenerate = 0L, method = "complement",
       nuisance_lambda_rule = "min")
}

fit_target_rcal_predictions <- function(train_data, eval_data, A_val,
                                        cv_seed, nlambda = 100L,
                                        rcal_backend = c("reference", "native")) {
  rcal_backend <- match.arg(rcal_backend)
  if (!identical(train_data$Z_site, train_data$W_outcome) ||
      !identical(eval_data$Z_site, eval_data$W_outcome)) {
    stop("This RCAL diagnostic requires the common FACE working basis.")
  }
  if (length(A_val) != 1L || !(A_val %in% c(0L, 1L)) ||
      length(nlambda) != 1L || nlambda < 2L || nlambda != as.integer(nlambda)) {
    stop("fit_target_rcal_predictions: invalid arm or penalty-grid size.")
  }
  # RCAL does not standardize x internally. Use only complement-training rows;
  # the same transformation applies to both models and held-out predictions.
  feature_center <- colMeans(train_data$Z_site)
  centered <- sweep(train_data$Z_site, 2L, feature_center, "-")
  feature_scale <- sqrt(colMeans(centered^2))
  if (any(!is.finite(feature_scale)) || any(feature_scale <= 0)) {
    stop("fit_target_rcal_predictions: constant or non-finite training column.")
  }
  x_train <- sweep(centered, 2L, feature_scale, "/")
  x_eval <- sweep(sweep(eval_data$Z_site, 2L, feature_center, "-"),
                  2L, feature_scale, "/")
  indicator <- as.numeric(train_data$A == A_val)
  arm_rows <- which(indicator == 1L)
  if (min(table(factor(indicator, levels = c(0, 1)))) < 10L ||
      min(table(factor(train_data$Y[arm_rows], levels = c(0, 1)))) < 8L) {
    stop("fit_target_rcal_predictions: insufficient support for binary CV.")
  }
  # Match the declared 100-candidate, 1e-4 relative range. RCAL's default
  # tune.fac=0.5 with nrho=100 would instead extend to approximately 1e-30.
  tune_factor <- 1e-4^(1 / (nlambda - 1L))
  fit_model <- function(x, y, loss, seed, weights = NULL) {
    if (rcal_backend == "native") {
      if (loss == "cal") return(fit_rcal_calibration_cv(x, y, nlambda, seed))
      return(fit_rcal_weighted_outcome_cv(x, y, weights, nlambda, seed))
    }
    RoCE:::with_seed(seed, RCAL::glm.regu.cv(
      fold = RoCE:::get_cv_fold_count(nrow(x), min_per_fold = 10L),
      nrho = nlambda, y = as.numeric(y), x = x, iw = weights,
      loss = loss, tune.fac = tune_factor, n.iter = 100L, eps = 1e-6
    ))
  }
  ps_fit <- fit_model(x_train, indicator, "cal", cv_seed)
  ps_coefficients <- ps_fit$sel.bet[, 1L]
  ps_train <- plogis(drop(cbind(1, x_train) %*% ps_coefficients))
  if (any(!is.finite(ps_coefficients)) || any(ps_train <= 0 | ps_train >= 1) ||
      max(abs(ps_train - ps_fit$sel.fit[, 1L])) > 1e-9) {
    stop("RCAL propensity coefficient/prediction contract failed.")
  }
  outcome_weights <- (1 - ps_train[arm_rows]) / ps_train[arm_rows]
  outcome_fit <- fit_model(
    x_train[arm_rows, , drop = FALSE], train_data$Y[arm_rows], "ml",
    cv_seed + 1L, outcome_weights
  )
  outcome_coefficients <- outcome_fit$sel.bet[, 1L]
  outcome_train <- plogis(drop(cbind(1, x_train[arm_rows, , drop = FALSE]) %*%
                                outcome_coefficients))
  if (any(!is.finite(outcome_coefficients)) ||
      max(abs(outcome_train - outcome_fit$sel.fit[, 1L])) > 1e-9) {
    stop("RCAL outcome coefficient/prediction contract failed.")
  }
  raw_probability <- plogis(drop(cbind(1, x_eval) %*% ps_coefficients))
  outcome <- plogis(drop(cbind(1, x_eval) %*% outcome_coefficients))
  diagnostics <- data.frame(
    A_val = A_val, n_train = train_data$n, n_arm = length(arm_rows),
    ps_lambda = ps_fit$sel.rho[1L], outcome_lambda = outcome_fit$sel.rho[1L],
    ps_invalid_candidates = sum(ps_fit$non.conv != 0),
    outcome_invalid_candidates = sum(outcome_fit$non.conv != 0),
    probability_clip_fraction = mean(raw_probability < RoCE:::PROP_SCORE_LOWER |
                                      raw_probability > RoCE:::PROP_SCORE_UPPER),
    calibration_intercept = mean(indicator / ps_train - 1),
    calibration_max = max(abs(colMeans(x_train * (indicator / ps_train - 1))))
  )
  list(arm_probability = RoCE:::clip_propensity(raw_probability),
       outcome = RoCE:::clip_outcome_pred(outcome, "binomial"),
       diagnostics = diagnostics)
}

target_fold_variants <- function(target_folds, truth, k1, A_val, k2 = NULL,
                                 baseline = NULL, nlambda = 100L,
                                 rcal_backend = c("reference", "native")) {
  rcal_backend <- match.arg(rcal_backend)
  n_folds <- length(target_folds)
  excluded <- if (is.null(k2)) k1 else c(k1, k2)
  eval_fold <- if (is.null(k2)) k1 else k2
  train_data <- RoCE:::combine_folds(target_folds, setdiff(seq_len(n_folds), excluded))
  eval_data <- RoCE:::materialize_fold(target_folds, eval_fold)
  stopifnot(!any(train_data$original_idx %in% eval_data$original_idx))
  if (is.null(baseline)) {
    baseline <- RoCE::estimate_target_only_from_complement(
      target_folds, k1, n_folds, k2 = k2, family = "binomial", A_val = A_val,
      nuisance_lambda_rule = "min"
    )
  }
  p_fitted <- if (A_val == 1L) baseline$prop_scores else 1 - baseline$prop_scores
  rows <- eval_data$original_idx
  m_true <- truth[[if (A_val == 1L) "mu1" else "mu0"]][rows]
  p_true <- if (A_val == 1L) truth$propensity[rows] else 1 - truth$propensity[rows]
  rcal <- fit_target_rcal_predictions(
    train_data, eval_data, A_val,
    cv_seed = RoCE:::.target_only_cv_seed(k1, k2, "outcome", A_val) + 2000000L,
    nlambda = nlambda, rcal_backend = rcal_backend
  )
  predictions <- list(
    lasso = list(baseline$m_pred, p_fitted),
    oracle_outcome = list(m_true, p_fitted),
    oracle_propensity = list(baseline$m_pred, p_true),
    oracle_both = list(m_true, p_true),
    rcal_propensity = list(baseline$m_pred, rcal$arm_probability),
    rcal_weighted = list(rcal$outcome, rcal$arm_probability)
  )
  fits <- lapply(predictions, function(prediction) {
    fit <- target_prediction_fit(eval_data, A_val, prediction[[1L]], prediction[[2L]])
    fit$expected_remainder <- (1 - p_true / prediction[[2L]]) *
      (prediction[[1L]] - m_true)
    fit
  })
  # This is a reconstruction check against the actual production helper.
  stopifnot(isTRUE(all.equal(fits$lasso$estimate, baseline$estimate, tolerance = 1e-12)),
            isTRUE(all.equal(fits$lasso$varphi_ot, baseline$varphi_ot, tolerance = 1e-12)))
  attr(fits, "diagnostics") <- transform(rcal$diagnostics, k1 = k1,
                                        k2 = if (is.null(k2)) NA_integer_ else k2)
  fits
}

reaggregate_target_variant <- function(data_split, folds, reference, variants,
                                      variant) {
  glm_spec <- RoCE:::resolve_glm_family("binomial")
  arms <- lapply(c(mu1 = 1L, mu0 = 0L), function(A_val) {
    arm_name <- if (A_val == 1L) "mu1" else "mu0"
    fitted_arm <- reference$arm_results[[arm_name]]
    fold_results <- fitted_arm$fold_results
    for (k1 in seq_along(fold_results)) {
      fold_results[[k1]]$target_only <- variants[[arm_name]][[k1]]$outer[[variant]]
      for (k2_name in names(fold_results[[k1]]$target_only_inner)) {
        fold_results[[k1]]$target_only_inner[[k2_name]] <-
          variants[[arm_name]][[k1]]$inner[[k2_name]][[variant]]
      }
    }
    RoCE:::.aggregate_crossfit_from_fitted_folds(
      data_split = data_split, target_folds = folds$target_folds,
      source_folds = folds$source_folds, fold_results = fold_results,
      n_folds = fitted_arm$n_folds,
      M_tau = fitted_arm$M_tau, M_tau_inference = fitted_arm$M_tau_inference,
      lambda_selection = fitted_arm$aggregation_lambda_selection,
      lambda_rule = "min", aggregation_lambda_grid = NULL,
      communication_mode = fitted_arm$communication_mode,
      family_int = glm_spec$family_int, link_int = glm_spec$link_int,
      A_val = A_val, verbose = FALSE
    )
  })
  RoCE:::calculate_tate_crossfit_aggregation(
    data_split = data_split, mu1_result = arms$mu1, mu0_result = arms$mu0,
    lambda_selection = 1, lambda_rule = "min", screening_rule = "soft_penalty",
    verbose = FALSE
  )
}
