# Target-local IPW outcome projections for the next-paper covariance diagnostic.
# These do not replace the current calibrated target anchor or source fits.
# Source candidate_score_summary.R before using the cross-fitted wrapper.

fit_target_projection_fold <- function(target_data, fold_labels, evaluation_fold,
                                        propensity_coefficients, arm,
                                        weighting = c("ipw", "unweighted"),
                                        nlambda = 100L, tolerance = 1e-10,
                                        propensity_radius = 12) {
  weighting <- match.arg(weighting)
  if (length(arm) != 1L || !is.numeric(arm) || !is.finite(arm) || !(arm %in% 0:1) ||
      length(evaluation_fold) != 1L || !is.finite(evaluation_fold) ||
      evaluation_fold != floor(evaluation_fold) ||
      length(fold_labels) != target_data$n || anyNA(fold_labels) ||
      any(fold_labels < 1 | fold_labels != floor(fold_labels)) ||
      !(evaluation_fold %in% fold_labels)) stop("Invalid arm or fold labels")
  if (length(nlambda) != 1L || !is.finite(nlambda) || nlambda < 2 || nlambda != floor(nlambda) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      length(propensity_radius) != 1L || !is.finite(propensity_radius) || propensity_radius <= 0) {
    stop("Invalid projection solver controls")
  }
  outcome_design <- as.matrix(target_data$W_outcome)
  weight_design <- cbind(1, target_data$Z_site)
  if (length(propensity_coefficients) != ncol(weight_design) ||
      any(!is.finite(propensity_coefficients))) stop("Invalid propensity coefficients")
  evaluation <- which(fold_labels == evaluation_fold)
  training <- which(fold_labels != evaluation_fold & target_data$A == arm)
  groups <- match(fold_labels[training], sort(unique(fold_labels[training])))
  if (length(unique(groups)) < 3L || length(unique(target_data$Y[training])) != 2L) {
    stop("Projection fitting needs three training folds and a nondegenerate binary response")
  }
  predictor <- drop(weight_design[training, , drop = FALSE] %*% propensity_coefficients)
  propensity <- plogis(pmax(-propensity_radius, pmin(propensity_radius, predictor)))
  weights <- if (weighting == "ipw") 1 / propensity else rep(1, length(training))
  fitted <- glmnet::cv.glmnet(outcome_design[training, , drop = FALSE], target_data$Y[training],
    family = "binomial", alpha = 1, weights = weights, foldid = groups,
    nlambda = nlambda, lambda.min.ratio = if (length(training) > ncol(outcome_design)) 1e-4 else .01,
    standardize = TRUE, thresh = tolerance, maxit = 1000000L)
  if (fitted$glmnet.fit$jerr != 0L) stop("Target projection path did not converge")
  prediction <- drop(predict(fitted, newx = outcome_design[evaluation, , drop = FALSE],
                             s = "lambda.min", type = "response"))
  if (any(!is.finite(prediction))) stop("Nonfinite target projection prediction")
  list(prediction = prediction, coefficients = as.numeric(coef(fitted, s = "lambda.min")),
    selected_lambda = fitted$lambda.min, evaluation_indices = evaluation,
    training_indices = training, cv_fold_labels = groups, weighting = weighting,
    propensity_clipped_fraction = mean(abs(predictor) > propensity_radius),
    weight_range = range(weights))
}

crossfit_target_common_projection <- function(fitted_tate, target_data,
                                               weighting = c("ipw", "unweighted"),
                                               nlambda = 100L, tolerance = 1e-10) {
  weighting <- match.arg(weighting)
  if (!identical(fitted_tate$family, "binomial") ||
      !identical(fitted_tate$target_nuisance_method, "hou_calibrated")) {
    stop("This diagnostic requires the retained binomial Hou target fits")
  }
  labels <- integer(target_data$n)
  assignments <- integer(target_data$n)
  for (fold in seq_len(fitted_tate$n_folds)) {
    index <- check_score_indices(fitted_tate$intermediates$fold_info[[fold]]$target_idx,
                                 target_data$n, "Target projection")
    labels[index] <- fold
    assignments <- assignments + tabulate(index, nbins = target_data$n)
  }
  if (any(assignments != 1L)) stop("Target folds must cover observations exactly once")
  predictions <- matrix(NA_real_, target_data$n, 2L, dimnames = list(NULL, c("mu0", "mu1")))
  fits <- vector("list", fitted_tate$n_folds)
  for (fold in seq_len(fitted_tate$n_folds)) {
    fits[[fold]] <- vector("list", 2L)
    names(fits[[fold]]) <- colnames(predictions)
    for (arm in 0:1) {
      name <- paste0("mu", arm)
      calibration <- fitted_tate$arm_results[[name]]$fold_results[[fold]]$target_only$calibration
      if (is.null(calibration$weight) || fold %in% calibration$calibration_folds ||
          !setequal(calibration$calibration_folds, setdiff(seq_len(fitted_tate$n_folds), fold))) {
        stop("Retained propensity model has inconsistent training-fold provenance")
      }
      fit <- fit_target_projection_fold(target_data, labels, fold, calibration$weight, arm,
        weighting, nlambda, tolerance, fitted_tate$M_tau_inference)
      predictions[fit$evaluation_indices, name] <- fit$prediction
      fits[[fold]][[name]] <- fit
    }
  }
  contrast <- predictions[, "mu1"] - predictions[, "mu0"]
  centered_by_fold <- contrast - ave(contrast, labels)
  list(predictions = predictions, contrast = contrast, folds = fits, weighting = weighting,
    fold_labels = labels,
    common_variance = as.numeric(score_mean_covariance(contrast)),
    common_variance_within_folds = as.numeric(score_mean_covariance(centered_by_fold)),
    scope = "Target-local projection diagnostic; finite-sample variance inference is not established")
}
