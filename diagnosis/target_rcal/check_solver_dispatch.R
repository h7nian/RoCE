#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L || !(args[1L] %in% c("coordinate_descent", "proximal_newton"))) {
    stop("usage: check_solver_dispatch.R SOLVER OUTPUT_CSV (one fresh R process per solver)")
  }
  solver <- args[1L]
  Sys.setenv(ROCE_NUISANCE_SOLVER = solver)
  .libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  # The C++ switch is a static local initialized on first use. Changing the
  # environment inside an already initialized process does not switch solvers.
  stopifnot(identical(RoCE:::nuisance_solver_cpp(), solver))
  set.seed(507)
  x <- matrix(rnorm(1500), 300, 5)
  a <- rbinom(300, 1, plogis(.2 + .3 * x[, 1L]))
  y <- rbinom(300, 1, plogis(.1 + .4 * x[, 2L]))
  weights <- exp(.1 * x[, 3L])
  rows <- list()
  fits <- list(
    final_outcome = RoCE:::fit_general_glm_cpp(
      X = x, Y = y, weights = weights, family_int = 1L, link_int = 1L,
      lambda = .03, max_iter = 500L, tol = 1e-8, warm_start = numeric()
    ),
    initial_tilting = RoCE:::fit_initial_density_ratio_cpp(
      Z_site = x, A_source = a, mean_phi = colMeans(cbind(1, x)),
      lambda = .1, max_iter = 500L, tol = 1e-8, A_val = 1L,
      M_tau = Inf, warm_start = numeric()
    )
  )
  for (model in names(fits)) {
    fit <- fits[[model]]
    coefficients <- if (model == "final_outcome") fit$alpha else fit$gamma
    rows[[length(rows) + 1L]] <- data.frame(
      model, requested_solver = solver, recorded_solver = RoCE:::nuisance_solver_cpp(),
      converged = fit$converged, iterations = fit$iterations,
      max_update = fit$max_update,
      coefficients = paste(sprintf("%.17g", coefficients), collapse = ";")
    )
  }
  rows <- do.call(rbind, rows)
  write.csv(rows, args[[2L]], row.names = FALSE)
  print(rows[, setdiff(names(rows), "coefficients")])
}
if (sys.nframe() == 0L) main()
