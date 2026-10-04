# The RCAL propensity objective is an existing RoCE exponential-tilting
# objective with mean_phi = mean((1 - I[A=a]) * (1, X)). Reusing the compiled
# optimizer avoids adding another solver. Fixed-penalty equivalence is tested
# against RCAL::glm.regu before this backend is used in simulations.

fit_rcal_calibration <- function(x, indicator, lambda, warm_start = numeric(),
                                 max_iter = 1000L) {
  design <- cbind(1, x)
  fit <- RoCE:::fit_initial_density_ratio_cpp(
    Z_site = x, A_source = indicator,
    mean_phi = colMeans(design * (1 - indicator)),
    lambda = lambda, max_iter = max_iter, tol = 1e-8,
    A_val = 1L, M_tau = Inf, warm_start = warm_start
  )
  coefficients <- as.numeric(fit$gamma)
  predictor <- drop(design %*% coefficients)
  # Reject numerical guards/boundaries; they are not solutions to the
  # untruncated calibrated propensity objective being compared here.
  odds_weights <- exp(-predictor)
  valid <- isTRUE(fit$converged) && all(is.finite(odds_weights)) &&
    all(odds_weights > RoCE:::WEIGHT_MIN & odds_weights < RoCE:::WEIGHT_MAX)
  kkt_error <- Inf
  if (valid) {
    gradient <- colMeans(design * (1 - indicator - indicator * odds_weights))
    penalized <- coefficients[-1L]
    residual <- ifelse(abs(penalized) > 1e-7,
                       abs(gradient[-1L] + lambda * sign(penalized)),
                       pmax(abs(gradient[-1L]) - lambda, 0))
    kkt_error <- max(abs(gradient[1L]), residual)
    valid <- kkt_error < 1e-5
  }
  list(coefficients = coefficients, converged = valid,
       kkt_error = kkt_error, iterations = fit$iterations)
}

fit_rcal_calibration_cv <- function(x, indicator, nlambda = 100L,
                                    cv_seed, n_folds = 5L) {
  n <- nrow(x)
  lambda_max <- max(abs(colMeans(x * (indicator / mean(indicator) - 1))))
  lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * 1e-4), length.out = nlambda))
  permutation <- RoCE:::with_seed(cv_seed, sample.int(n))
  fold_id <- integer(n)
  fold_id[permutation] <- rep(seq_len(n_folds), each = ceiling(n / n_folds))[seq_len(n)]
  validation_loss <- matrix(Inf, nlambda, n_folds)
  for (fold in seq_len(n_folds)) {
    training <- fold_id != fold
    evaluation <- !training
    x_training <- x[training, , drop = FALSE]
    x_evaluation <- x[evaluation, , drop = FALSE]
    indicator_training <- indicator[training]
    indicator_evaluation <- indicator[evaluation]
    warm_start <- numeric()
    for (index in seq_along(lambda_grid)) {
      fit <- fit_rcal_calibration(x_training, indicator_training,
                                 lambda_grid[index], warm_start, max_iter = 100L)
      if (!fit$converged) break
      warm_start <- fit$coefficients
      predictor <- drop(cbind(1, x_evaluation) %*% fit$coefficients)
      loss <- mean(indicator_evaluation * exp(-predictor) +
                     (1 - indicator_evaluation) * predictor)
      if (is.finite(loss)) validation_loss[index, fold] <- loss
    }
  }
  loss <- rowMeans(validation_loss)
  if (!any(is.finite(loss))) stop("RCAL calibration CV: no valid penalty.")
  selected <- which.min(loss)
  fit <- fit_rcal_calibration(x, indicator, lambda_grid[selected])
  if (!fit$converged) stop("RCAL calibration final fit failed its KKT check.")
  fitted <- plogis(drop(cbind(1, x) %*% fit$coefficients))
  list(sel.bet = matrix(fit$coefficients, ncol = 1L),
       sel.fit = matrix(fitted, ncol = 1L), sel.rho = lambda_grid[selected],
       non.conv = as.integer(!is.finite(loss)), kkt_error = fit$kkt_error)
}

fit_rcal_weighted_outcome_cv <- function(x, y, weights, nlambda, cv_seed) {
  fit <- RoCE:::with_seed(cv_seed, glmnet::cv.glmnet(
    x = x, y = y, weights = weights, family = "binomial", alpha = 1,
    nfolds = RoCE:::get_cv_fold_count(nrow(x), min_per_fold = 10L),
    nlambda = nlambda, lambda.min.ratio = 1e-4,
    standardize = FALSE, thresh = 1e-9, maxit = RoCE:::GLMNET_MAX_ITER
  ))
  coefficients <- drop(as.matrix(coef(fit, s = "lambda.min")))
  list(sel.bet = matrix(coefficients, ncol = 1L),
       sel.fit = matrix(plogis(drop(cbind(1, x) %*% coefficients)), ncol = 1L),
       # glmnet does not expose RCAL's per-candidate convergence vector. Do not
       # equate a returned prediction with zero failures over the entire path.
       sel.rho = fit$lambda.min, non.conv = rep(NA_integer_, length(fit$lambda)))
}
