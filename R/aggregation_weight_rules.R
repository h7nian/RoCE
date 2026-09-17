# Smooth quadratic-bias source weights: the pre-specified sensitivity rule
# computed alongside the truncated-Wald rule (main.tex
# rem:quadratic_bias_rule). With n_t the inner training size, the weights
# minimize n_t Var(eta) + n_t^power sum_j delta_j^2 eta_j^2, whose stationarity
# condition is (Q + n_t^(power - 1) diag(delta^2)) eta = -l in the notation of
# .variance_quadratic_form(). No source weight is set exactly to zero and no
# ridge is added: like the Stage-1 candidate, the rule requires the variance
# quadratic itself to be positive definite.
.quadratic_bias_weights <- function(moments, power = AGG_QUADRATIC_BIAS_POWER) {
  form <- .variance_quadratic_form(moments, psd_ridge = 0)
  discrepancy <- moments$avg_source_est - moments$avg_target_est
  K <- length(discrepancy)
  curvature <- form$Q + diag(moments$n_t^(power - 1) * discrepancy^2, K)
  if (any(!is.finite(form$l)) || any(!is.finite(curvature))) {
    stop(".quadratic_bias_weights: non-finite weight moments.", call. = FALSE)
  }
  .require_positive_definite(form$Q, "the variance quadratic")
  scale <- .require_positive_definite(curvature, "the penalized weight objective")
  normalized <- curvature / outer(scale, scale)
  cholesky <- chol(normalized)
  weights <- drop(backsolve(cholesky, forwardsolve(t(cholesky), -form$l / scale)) / scale)
  residual <- max(abs(drop(curvature %*% weights) + form$l))
  residual_scale <- max(abs(form$l), abs(curvature) * max(abs(weights)), .Machine$double.xmin)
  if (any(!is.finite(weights)) || residual > WEIGHT_LAYER_KKT_TOLERANCE * residual_scale) {
    stop(sprintf(".quadratic_bias_weights: normal equations not solved (residual %.3e).", residual),
         call. = FALSE)
  }
  weights
}

# Checks that a symmetric matrix is positive definite after diagonal scaling
# and returns the scaling (square roots of the diagonal).
.require_positive_definite <- function(matrix, label) {
  curvature <- diag(matrix)
  if (any(!is.finite(curvature)) || any(curvature <= 0)) {
    stop(sprintf(".quadratic_bias_weights: %s lacks positive diagonal curvature.", label),
         call. = FALSE)
  }
  scale <- sqrt(curvature)
  positive_definite <- !inherits(
    try(chol(matrix / outer(scale, scale)), silent = TRUE), "try-error"
  )
  if (!positive_definite) {
    stop(sprintf(".quadratic_bias_weights: %s is not positive definite.", label),
         call. = FALSE)
  }
  scale
}
