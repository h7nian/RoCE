# Fold-specific eta learning: reuse outer final fits (two layers), or evaluate
# separately calibrated inner fits (three layers). Outer evaluation is excluded
# from both procedures. Legacy initial-model validation remains explicit.

.resolve_crossfit_validation <- function(method, crossfit_layers, calibration_control,
                                          target_nuisance_method, communication_mode = NULL) {
  if (!is.null(crossfit_layers)) {
    if (!is.numeric(crossfit_layers) || length(crossfit_layers) != 1L ||
        !is.finite(crossfit_layers) || !crossfit_layers %in% c(2, 3)) {
      stop("crossfit_layers must be 2, 3, or NULL.", call. = FALSE)
    }
    if (calibration_control$source_nuisance_method == "standard" && target_nuisance_method == "lasso") {
      stop("crossfit_layers selects a calibrated procedure; at least one nuisance program must use final calibration.", call. = FALSE)
    }
    legacy_default <- identical(method, c("initial", "calibrated")) ||
      identical(method, c("initial", "calibrated", "complete"))
    if (!legacy_default) method <- match.arg(method, c("initial", "calibrated", "complete", "outer_fit"))
    if (!legacy_default && identical(method, "initial")) {
      stop("crossfit_layers cannot be combined with legacy source_validation_method='initial'.", call. = FALSE)
    }
    if (crossfit_layers == 2L) {
      method <- "outer_fit"
    } else {
      if (identical(method, "outer_fit")) stop("outer_fit weight learning uses two layers, not three.", call. = FALSE)
      if (legacy_default) method <- "calibrated"
    }
  }
  method <- .match_source_validation(method, calibration_control, communication_mode)
  if (method == "outer_fit" && calibration_control$source_nuisance_method == "standard" &&
      target_nuisance_method == "lasso") {
    stop("outer_fit requires a calibrated source or target nuisance program.", call. = FALSE)
  }
  method
}

.reuse_outer_target_scores <- function(outer_fit, target_folds, k1, k2, A_val, family, radius) {
  allowed <- setdiff(seq_along(target_folds), k1)
  if (!k2 %in% allowed) stop("Target weight-learning scores must exclude the outer evaluation fold.", call. = FALSE)
  evaluation <- materialize_fold(target_folds, k2)
  if (!is.null(outer_fit$calibration)) {
    if (!identical(as.integer(outer_fit$calibration$calibration_folds), as.integer(allowed))) {
      stop("The reused target calibration must use exactly the outer-training folds.", call. = FALSE)
    }
    result <- .evaluate_target_calibration(outer_fit$calibration, evaluation, A_val, family, radius)
  } else {
    records <- outer_fit$training_scores
    if (is.null(records) || !identical(as.integer(records$training_folds), as.integer(allowed))) {
      stop("The ordinary outer target fit lacks its matching training scores.", call. = FALSE)
    }
    expected <- unlist(lapply(target_folds[allowed], `[[`, "original_idx"), use.names = FALSE)
    if (!identical(records$original_idx, expected) || anyDuplicated(expected)) {
      stop("Target training-score observation indices do not match the allowed folds.", call. = FALSE)
    }
    rows <- match(evaluation$original_idx, records$original_idx)
    if (anyNA(rows)) stop("Target weight-learning fold is missing from the stored training scores.", call. = FALSE)
    phi <- records$phi[rows]
    estimate <- mean(phi)
    centered <- phi - estimate
    result <- list(estimate = estimate, V_ot = mean(centered^2), variance = mean(centered^2)/length(phi),
      varphi_ot = centered, prop_scores = records$prop_scores[rows], m_pred = records$m_pred[rows],
      method = "outer_fit", nuisance_lambda_rule = outer_fit$nuisance_lambda_rule)
  }
  # Prediction reuse does not constitute another fitted degenerate model.
  result$outcome_degenerate <- 0L
  result$reused_outer_fit <- TRUE
  result$training_folds <- as.integer(allowed)
  result
}
