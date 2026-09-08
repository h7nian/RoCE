# Experimental structural-score coefficient solver. Not a deployed estimator.
# Minimize 0.5*a' H a - d' a + sum_j penalty_j*abs(a_j).
# No matrix inverse, statistical ridge, or confidence-interval adjustment.

.projection_failure <- function(message, diagnostics = list()) {
  stop(structure(list(message = message, call = NULL, diagnostics = diagnostics),
                 class = c("roce_projection_error", "error", "condition")))
}

.solve_sparse_moment_projection <- function(curvature, gradient, penalty,
                                             tolerance = 1e-9,
                                             max_iterations = 10000L) {
  if (!is.matrix(curvature) || !is.numeric(curvature) ||
      nrow(curvature) < 1L || nrow(curvature) != ncol(curvature) ||
      any(!is.finite(curvature))) stop("curvature must be a finite nonempty square matrix")
  dimension <- nrow(curvature)
  if (!is.numeric(gradient) || !is.null(dim(gradient)) ||
      length(gradient) != dimension || any(!is.finite(gradient))) {
    stop("gradient must be a finite vector matching curvature")
  }
  if (!is.numeric(penalty) || !is.null(dim(penalty)) ||
      !length(penalty) %in% c(1L, dimension) ||
      any(!is.finite(penalty)) || any(penalty < 0)) {
    stop("penalty must be a nonnegative scalar or a matching finite vector")
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance <= 0 || tolerance >= 1 ||
      !is.numeric(max_iterations) || length(max_iterations) != 1L ||
      !is.finite(max_iterations) || max_iterations < 1 ||
      max_iterations != floor(max_iterations) || max_iterations > .Machine$integer.max) {
    stop("invalid solver tolerance or iteration limit")
  }
  penalty <- rep(penalty, length.out = dimension)
  scale <- max(abs(curvature))
  if (scale == 0) {
    if (any(abs(gradient) > penalty)) .projection_failure(
      "zero curvature with linear drift beyond the L1 penalty")
    return(list(coefficients = rep(0, dimension), iterations = 0L,
      max_kkt_residual = 0, relative_kkt_residual = 0,
      max_moment_residual = max(abs(gradient)), max_original_kkt_residual = 0,
      objective = 0,
      numerical_psd_adjustment = 0, input_minimum_eigenvalue = 0,
      numerical_rank = 0L, penalty = penalty))
  }
  H <- curvature/scale
  d <- gradient/scale
  lambda <- penalty/scale
  numerical_tolerance <- 100*.Machine$double.eps*dimension
  if (max(abs(H-t(H))) > numerical_tolerance) stop("curvature must be symmetric")
  H <- (H+t(H))/2
  spectral <- eigen(H, symmetric = TRUE)
  minimum_eigenvalue <- min(spectral$values)
  eigen_tolerance <- numerical_tolerance*max(1, max(abs(spectral$values)))
  if (minimum_eigenvalue < -eigen_tolerance) stop("curvature must be positive semidefinite")
  adjustment <- 0
  if (minimum_eigenvalue < 0) {
    # Remove negative eigenvalues only within floating-point tolerance; never
    # add a positive ridge to turn an unidentified direction into a fitted one.
    original_H <- H
    H <- tcrossprod(sweep(spectral$vectors, 2L,
                          sqrt(pmax(spectral$values, 0)), "*"))
    adjustment <- max(abs(H-original_H))*scale
  }
  near_null <- which(spectral$values <= eigen_tolerance)
  if (length(near_null)) {
    for (j in near_null) {
      direction <- spectral$vectors[, j]
      if (abs(sum(d*direction)) > sum(lambda*abs(direction)) +
          numerical_tolerance*max(1, max(abs(d)))) {
        .projection_failure("unsupported near-null direction: linear drift exceeds the L1 penalty",
                            list(input_minimum_eigenvalue = minimum_eigenvalue*scale))
      }
    }
  }
  coefficients <- numeric(dimension)
  residual <- -d
  kkt_scale <- max(1, max(abs(d)), max(lambda))
  evaluate_kkt <- function() {
    residual <- drop(H %*% coefficients)-d
    active <- coefficients != 0
    violations <- pmax(abs(residual)-lambda, 0)
    violations[active] <- abs(residual[active]+lambda[active]*sign(coefficients[active]))
    list(residual = residual, maximum = max(violations))
  }
  finish <- function(iterations, kkt) {
    # Also report residuals against the original supplied moment system.
    original_residual <- drop(curvature %*% coefficients)-gradient
    original_violations <- pmax(abs(original_residual)-penalty, 0)
    active <- coefficients != 0
    original_violations[active] <- abs(original_residual[active]+
                                       penalty[active]*sign(coefficients[active]))
    if (kkt$maximum <= tolerance*kkt_scale &&
        max(original_violations) > 10*tolerance*kkt_scale*scale) {
      .projection_failure("numerical PSD adjustment failed the original-system KKT check",
                          list(max_original_kkt_residual = max(original_violations)))
    }
    list(coefficients = coefficients, iterations = iterations,
      max_kkt_residual = kkt$maximum*scale,
      relative_kkt_residual = kkt$maximum/kkt_scale,
      max_moment_residual = max(abs(original_residual)),
      max_original_kkt_residual = max(original_violations),
      objective = .5*sum(coefficients*drop(curvature %*% coefficients))-
        sum(gradient*coefficients)+sum(penalty*abs(coefficients)),
      numerical_psd_adjustment = adjustment,
      input_minimum_eigenvalue = minimum_eigenvalue*scale,
      numerical_rank = sum(spectral$values > eigen_tolerance), penalty = penalty)
  }
  kkt <- evaluate_kkt()
  if (kkt$maximum <= tolerance*kkt_scale) return(finish(0L, kkt))
  for (iteration in seq_len(as.integer(max_iterations))) {
    for (j in seq_len(dimension)) {
      curvature_j <- H[j, j]
      if (curvature_j <= numerical_tolerance) {
        if (abs(residual[j]) > lambda[j]+tolerance*kkt_scale) {
          .projection_failure("insufficient coordinate curvature",
                              list(coordinate = j, iterations = iteration))
        }
        next
      }
      partial <- curvature_j*coefficients[j]-residual[j]
      updated <- sign(partial)*max(abs(partial)-lambda[j], 0)/curvature_j
      if (!is.finite(updated)) .projection_failure("nonfinite projection coefficient")
      change <- updated-coefficients[j]
      if (change != 0) {
        coefficients[j] <- updated
        residual <- residual+H[, j]*change
      }
    }
    # Recalculate from scratch each sweep to prevent cached-gradient drift.
    kkt <- evaluate_kkt()
    residual <- kkt$residual
    if (kkt$maximum <= tolerance*kkt_scale) return(finish(iteration, kkt))
  }
  .projection_failure("projection solver reached its iteration limit without satisfying KKT",
                      finish(as.integer(max_iterations), kkt))
}

.solve_joint_nuisance_projection <- function(curvatures, couplings, gradients,
                                              penalties, tolerance = 1e-9,
                                              max_iterations = 10000L) {
  blocks <- c("alpha_initial", "gamma_initial", "alpha_final", "gamma_final")
  if (!is.list(curvatures) || !identical(names(curvatures), blocks) ||
      !is.list(penalties) || !identical(names(penalties), blocks) ||
      !is.list(couplings) || !identical(names(couplings),
        c("alpha_final_gamma_initial", "gamma_final_alpha_initial")) ||
      !is.list(gradients) || !identical(names(gradients), c("alpha_final", "gamma_final"))) {
    stop("joint projection inputs must contain the exact named nuisance blocks")
  }
  dimensions <- vapply(curvatures, function(x) {
    if (!is.matrix(x) || nrow(x) != ncol(x) || nrow(x) < 1L) stop("invalid curvature block")
    nrow(x)
  }, integer(1L))
  check_coupling <- function(x, rows, columns) {
    if (!is.matrix(x) || !is.numeric(x) || any(!is.finite(x)) ||
        !identical(dim(x), as.integer(c(rows, columns)))) stop("invalid coupling dimensions/values")
  }
  check_coupling(couplings$alpha_final_gamma_initial, dimensions[["alpha_final"]],
                 dimensions[["gamma_initial"]])
  check_coupling(couplings$gamma_final_alpha_initial, dimensions[["gamma_final"]],
                 dimensions[["alpha_initial"]])
  fit <- function(block, gradient) .solve_sparse_moment_projection(
    curvatures[[block]], gradient, penalties[[block]], tolerance, max_iterations)
  fits <- setNames(vector("list", length(blocks)), blocks)
  fits$alpha_final <- fit("alpha_final", -gradients$alpha_final)
  fits$gamma_final <- fit("gamma_final", gradients$gamma_final)
  fits$alpha_initial <- fit("alpha_initial", drop(crossprod(
    couplings$gamma_final_alpha_initial, fits$gamma_final$coefficients)))
  fits$gamma_initial <- fit("gamma_initial", -drop(crossprod(
    couplings$alpha_final_gamma_initial, fits$alpha_final$coefficients)))
  list(coefficients = unlist(lapply(fits, `[[`, "coefficients"), use.names = FALSE),
       block_fits = fits)
}
