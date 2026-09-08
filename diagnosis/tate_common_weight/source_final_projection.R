# Final calibration blocks only; initial-fit coupling corrections remain required.
.source_final_derivative_rows <- function(state, evaluated) {
  index <- state$indices; theta <- state$theta
  target <- state$target; source <- state$source
  mt <- plogis(drop(target$W %*% theta[index$alpha_final]))
  ms <- plogis(drop(source$W %*% theta[index$alpha_final]))
  log_weight <- drop(source$Z %*% theta[index$gamma_final])
  weight <- exp(-pmax(-state$M_inference, pmin(state$M_inference, log_weight)))
  rows <- list(alpha_final = list(target = target$W*(mt*(1-mt)),
    source = -source$W*(source$A*weight*ms*(1-ms))),
    gamma_final = list(target = matrix(0, nrow(target$Z), ncol(target$Z)),
      source = -source$Z*(source$A*weight*(source$Y-ms)*(abs(log_weight) < state$M_inference))))
  for (name in names(rows)) {
    x <- rows[[name]]
    stopifnot(all(is.finite(x$target)), all(is.finite(x$source)),
      max(abs(colMeans(x$target)+colMeans(x$source)-evaluated$score_gradient[index[[name]]])) < 1e-10)
  }
  rows
}

.source_projection_noise_penalty <- function(target_rows, source_rows, penalty_scale) {
  if (!is.matrix(target_rows) || !is.matrix(source_rows) || !is.numeric(target_rows) ||
      !is.numeric(source_rows) || ncol(target_rows) < 1L || ncol(target_rows) != ncol(source_rows) ||
      min(nrow(target_rows), nrow(source_rows)) < 2L ||
      any(!is.finite(target_rows)) || any(!is.finite(source_rows))) stop("invalid two-site derivative rows")
  if (!is.numeric(penalty_scale) || length(penalty_scale) != 1L ||
      !is.null(dim(penalty_scale)) || !is.finite(penalty_scale) || penalty_scale <= 0)
    stop("penalty_scale must be positive")
  penalty_scale*sqrt(log(2*ncol(target_rows))*(
    apply(target_rows, 2L, var)/nrow(target_rows)+apply(source_rows, 2L, var)/nrow(source_rows)))
}

.fit_source_final_projection <- function(state, evaluated, penalty_scale) {
  rows <- .source_final_derivative_rows(state, evaluated)
  fits <- list()
  for (name in names(rows)) {
    index <- state$indices[[name]]
    sign <- if (name == "alpha_final") -1 else 1
    penalty <- .source_projection_noise_penalty(rows[[name]]$target, rows[[name]]$source, penalty_scale)
    fits[[name]] <- .solve_sparse_moment_projection(
      sign*evaluated$jacobian[index, index, drop = FALSE],
      sign*evaluated$score_gradient[index], penalty, tolerance = 1e-8, max_iterations = 2000L)
  }
  fits
}
