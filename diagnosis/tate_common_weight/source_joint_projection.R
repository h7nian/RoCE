# Backward source adjoint solve using training-only candidate penalties.
# Requires source_final_projection.R and sparse_moment_projection.R.
.source_initial_coupling_rows <- function(state, key, final_fits) {
  index <- state$indices; theta <- state$theta; piece <- state$pieces[[key]]
  tc <- piece$target_calibration; sc <- piece$source_calibration
  a0 <- theta[index[[paste0("alpha_initial_", key)]]]
  g0 <- theta[index[[paste0("gamma_initial_", key)]]]
  at <- drop(tc$W %*% a0); as <- drop(sc$W %*% a0); gs <- drop(sc$Z %*% g0)
  mt <- plogis(at); ms <- plogis(pmax(-state$M_fit, pmin(state$M_fit, as)))
  mf <- plogis(drop(sc$W %*% theta[index$alpha_final]))
  w <- exp(-drop(sc$Z %*% theta[index$gamma_final]))
  w0 <- exp(-pmax(-state$M_fit, pmin(state$M_fit, gs)))
  fraction <- nrow(sc$W)/nrow(state$source$W)
  alpha_final <- final_fits$alpha_final$coefficients
  gamma_final <- final_fits$gamma_final$coefficients
  list(alpha_initial = list(
    target = tc$W*(mt*(1-mt)*(1-2*mt)*drop(tc$Z %*% gamma_final)/length(state$keys)),
    source = -sc$W*(sc$A*w*ms*(1-ms)*(1-2*ms)*(abs(as) < state$M_fit)*
      drop(sc$Z %*% gamma_final)*fraction)),
    gamma_initial = list(target = matrix(0, nrow(tc$Z), ncol(tc$Z)),
      source = sc$Z*(sc$A*w0*(sc$Y-mf)*(abs(gs) < state$M_fit)*
        drop(sc$W %*% alpha_final)*fraction)))
}

.fit_source_joint_projection <- function(state, evaluated, penalty_scale,
                                         initial_penalty_scale = penalty_scale) {
  if (!is.null(state$source_calibration_denominators)) stop("fit coefficients on training states, not validation views")
  index <- state$indices; J <- evaluated$jacobian; D <- evaluated$score_gradient
  final_fits <- .fit_source_final_projection(state, evaluated, penalty_scale)
  fits <- final_fits
  coefficients <- penalties <- numeric(length(D))
  for (name in names(final_fits)) {
    coefficients[index[[name]]] <- final_fits[[name]]$coefficients
    penalties[index[[name]]] <- final_fits[[name]]$penalty
  }
  maximum_coupling_error <- 0
  for (key in state$keys) {
    rows <- .source_initial_coupling_rows(state, key, final_fits)
    for (type in names(rows)) {
      name <- paste0(type, "_", key); pos <- index[[name]]
      sign <- if (type == "alpha_initial") -1 else 1
      downstream <- if (type == "alpha_initial") "gamma_final" else "alpha_final"
      rhs <- -sign*drop(crossprod(J[index[[downstream]], pos, drop = FALSE], coefficients[index[[downstream]]]))
      reconstructed <- colMeans(rows[[type]]$target)+colMeans(rows[[type]]$source)
      maximum_coupling_error <- max(maximum_coupling_error, abs(rhs-reconstructed))
      stopifnot(maximum_coupling_error < 1e-10, all(D[pos] == 0))
      if (is.null(initial_penalty_scale)) {
        # Explicit no-initial-correction candidate, not a solved adjoint block.
        penalty <- rep(Inf, length(pos))
        answer <- list(coefficients = numeric(length(pos)), penalty = penalty,
                       null_projection = TRUE)
      } else {
        penalty <- .source_projection_noise_penalty(rows[[type]]$target, rows[[type]]$source, initial_penalty_scale)
        answer <- .solve_sparse_moment_projection(sign*J[pos, pos, drop = FALSE], rhs,
          penalty, tolerance = 1e-8, max_iterations = 2000L)
      }
      fits[[name]] <- answer
      coefficients[pos] <- answer$coefficients; penalties[pos] <- penalty
    }
  }
  stopifnot(setequal(names(fits), names(index)))
  residual <- drop(crossprod(J, coefficients))-D
  # Unpenalized adjoint residual is allowed only within the selected L1 bands.
  maximum_excess <- max(pmax(abs(residual)-penalties, 0))
  stopifnot(maximum_excess < 1e-6)
  list(coefficients = coefficients, penalties = penalties, block_fits = fits,
    maximum_coupling_error = maximum_coupling_error, maximum_adjoint_excess = maximum_excess)
}
