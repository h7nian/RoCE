# Training-only coefficient policy, fixed by TARGET_PROJECTION_VALIDATION_PROTOCOL.md.
.target_score_derivative_rows <- function(block, theta, limits, system) {
  index <- system$index; W <- block$W; Z <- block$Z; A <- block$A; Y <- block$Y
  raw_p <- plogis(drop(Z %*% theta[index$propensity]))
  raw1 <- plogis(drop(W %*% theta[index$outcome_treated]))
  raw0 <- plogis(drop(W %*% theta[index$outcome_control]))
  p <- system$predictions$p; m1 <- system$predictions$m1; m0 <- system$predictions$m0
  dp <- raw_p*(1-raw_p)*(raw_p > limits$ps[1] & raw_p < limits$ps[2])
  d1 <- raw1*(1-raw1)*(raw1 > limits$outcome[1] & raw1 < limits$outcome[2])
  d0 <- raw0*(1-raw0)*(raw0 > limits$outcome[1] & raw0 < limits$outcome[2])
  rows <- cbind(Z*((-A*(Y-m1)/p^2-(1-A)*(Y-m0)/(1-p)^2)*dp),
                W*((1-A/p)*d1), -W*((1-(1-A)/(1-p))*d0))
  stopifnot(all(is.finite(rows)), max(abs(colMeans(rows)-system$gradient)) < 1e-10)
  rows
}

.fit_target_projection_candidate <- function(system, derivative_rows, penalty_scale) {
  if (!is.list(system)) stop("invalid target nuisance system")
  dimension <- length(system$gradient)
  expected_blocks <- c("propensity", "outcome_treated", "outcome_control")
  if (!is.numeric(system$gradient) || dimension < 1L ||
      any(!is.finite(system$gradient)) || !is.null(dim(system$gradient)) ||
      !is.matrix(system$jacobian) || !is.numeric(system$jacobian) ||
      !identical(dim(system$jacobian), c(dimension, dimension)) ||
      any(!is.finite(system$jacobian)) || !is.list(system$index) ||
      !identical(names(system$index), expected_blocks)) stop("invalid target nuisance system")
  indices <- unlist(system$index, use.names = FALSE)
  if (any(lengths(system$index) == 0L) || !is.numeric(indices) ||
      any(!is.finite(indices)) || any(indices != floor(indices)) ||
      !identical(sort(as.integer(indices)), seq_len(dimension))) stop("invalid nuisance block partition")
  if (!is.numeric(penalty_scale) || length(penalty_scale) != 1L ||
      !is.null(dim(penalty_scale)) || !is.finite(penalty_scale) || penalty_scale <= 0) {
    stop("projection penalty_scale must be positive")
  }
  if (!is.matrix(derivative_rows) || !is.numeric(derivative_rows) ||
      ncol(derivative_rows) != dimension || nrow(derivative_rows) < 2L ||
      any(!is.finite(derivative_rows))) stop("invalid derivative rows")
  penalties <- penalty_scale*apply(derivative_rows, 2L, sd)*sqrt(log(2*dimension)/nrow(derivative_rows))
  coefficients <- numeric(dimension)
  fits <- list()
  failure <- tryCatch({
    for (block in names(system$index)) {
      index <- system$index[[block]]
      fit <- .solve_sparse_moment_projection(-system$jacobian[index, index, drop = FALSE],
        -system$gradient[index], penalties[index], tolerance = 1e-8, max_iterations = 2000L)
      coefficients[index] <- fit$coefficients
      fits[[block]] <- fit
    }
    NULL
  }, error = function(error) {
    if (!inherits(error, "roce_projection_error")) stop(error)
    error
  })
  list(coefficients = coefficients, penalties = penalties, block_fits = fits, failure = failure)
}

.select_target_projection_candidate <- function(losses, validation_folds = 2:5) {
  if (!is.numeric(validation_folds) || !length(validation_folds) %in% c(3L, 4L) ||
      anyNA(validation_folds) || anyDuplicated(validation_folds) ||
      any(!validation_folds %in% 1:5)) stop("invalid validation folds")
  candidates <- c("null", "c0.5", "c1", "c2")
  required <- c("candidate", "fold", "n_validation", "loss", "failed")
  if (!is.data.frame(losses) || !all(required %in% names(losses)) ||
      nrow(losses) != length(candidates)*length(validation_folds) ||
      anyDuplicated(paste(losses$candidate, losses$fold)) ||
      !setequal(losses$candidate, candidates) || !setequal(losses$fold, validation_folds) ||
      any(!is.finite(losses$n_validation)) || any(losses$n_validation <= 0) ||
      any(losses$n_validation != floor(losses$n_validation)) ||
      !is.logical(losses$failed) || anyNA(losses$failed)) stop("invalid complete CV loss table")
  for (fold in validation_folds) {
    x <- losses[losses$fold == fold, ]
    if (nrow(x) != 4L || length(unique(x$n_validation)) != 1L) stop("inconsistent CV fold sizes")
  }
  if (any(!is.finite(losses$loss[!losses$failed]))) stop("nonfinite successful CV loss")
  null <- losses[losses$candidate == "null", ]
  if (any(null$failed) || any(null$loss != 0)) stop("null projection must have zero loss in every fold")
  scores <- do.call(rbind, lapply(candidates, function(candidate) {
    x <- losses[losses$candidate == candidate, ]
    eligible <- nrow(x) == length(validation_folds) && !any(x$failed)
    data.frame(candidate, eligible, failed_folds = sum(x$failed),
      validation_risk = if (eligible) weighted.mean(x$loss, x$n_validation) else NA_real_)
  }))
  valid <- scores[scores$eligible, ]
  best <- min(valid$validation_risk)
  ties <- valid[valid$validation_risk <= best+1e-10*max(1, max(abs(valid$validation_risk))), ]
  priority <- c("null", "c2", "c1", "c0.5")
  selected <- ties$candidate[which.min(match(ties$candidate, priority))]
  list(selected = selected, scores = scores,
       all_nonnull_ineligible = !any(scores$eligible[scores$candidate != "null"]))
}
