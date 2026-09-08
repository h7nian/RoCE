#!/usr/bin/env Rscript

# Deterministic one-source Gaussian counterexample, not a DGP simulation or an
# inference correction. Numerical safeguards are inactive at the chosen n=1000.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: check_wald_weight_limit.R V19_LIBRARY")
  library_root <- normalizePath(args[[1L]], mustWork = TRUE)
  .libPaths(unique(c(library_root, .libPaths())))
  suppressPackageStartupMessages(library("RoCE", lib.loc = library_root))
  stopifnot(identical(normalizePath(find.package("RoCE")), file.path(library_root, "RoCE")))

  # N_all * Var_hat(eta) = 1 - eta + eta^2 for these compatible score moments.
  # Adding (|u|-1)_+ |eta| yields this actual soft-threshold special case.
  gain <- function(u) pmax(0, 0.5 - 0.5 * pmax(abs(u) - 1, 0))
  variance_discrepancy <- (0.5 + 0.25 + 0.25 - 2 * 0.25) / 1000
  stopifnot(variance_discrepancy > RoCE:::VARIANCE_MIN)
  u <- seq(-3, 3, by = 0.25)
  observed <- vapply(u, function(z) as.numeric(RoCE::optimize_weights(
    estimates = c(s1 = 0.2 + z * sqrt(variance_discrepancy)),
    variances = list(V_ot = 0.5, V_t = 0.25, V_s = 0.25),
    C_ot = 0.25, n_samples = list(n_t = 1000, n_s = 1000),
    lambda = 1, mu_ot = 0.2, C_cross = matrix(0, 1L, 1L), clip_weights = FALSE
  )), numeric(1L))
  optimizer_error <- max(abs(observed - gain(u)))
  stopifnot(optimizer_error < 1e-12)

  integral <- function(fun, lower, upper) stats::integrate(
    fun, lower, upper, rel.tol = 1e-11, abs.tol = 1e-12,
    subdivisions = 1000L, stop.on.error = TRUE
  )$value
  mean_gain <- 2 * (integral(function(z) gain(z) * dnorm(z), 0, 1) +
                     integral(function(z) gain(z) * dnorm(z), 1, 2))
  mean_squared_gain <- 2 * (integral(function(z) gain(z)^2 * dnorm(z), 0, 1) +
                             integral(function(z) gain(z)^2 * dnorm(z), 1, 2))
  mean_proxy <- 1 - mean_gain + mean_squared_gain
  q_closed <- function(r, m) {
    a <- sqrt(m - 1)
    -(pnorm(2 * a - r) - pnorm(a - r) - pnorm(-a - r) +
        pnorm(-2 * a - r)) / (2 * a)
  }
  q_direct <- function(r, m) {
    a <- sqrt(m - 1)
    endpoints <- c(-Inf, -2 * a - r, -a - r, a - r, 2 * a - r, Inf)
    sum(vapply(seq_len(length(endpoints) - 1L), function(index) integral(
      function(z) z * gain((z + r) / a) * dnorm(z),
      endpoints[[index]], endpoints[[index + 1L]]
    ), numeric(1L)))
  }
  q_error <- max(vapply(c(2L, 3L, 5L), function(m) {
    max(vapply(c(-2, -0.7, 0, 0.7, 2), function(r)
      abs(q_closed(r, m) - q_direct(r, m)), numeric(1L)))
  }, numeric(1L)))
  stopifnot(q_error < 1e-9)
  result <- do.call(rbind, lapply(c(2L, 3L, 5L), function(m) {
    cross_covariance <- if (m == 2L) 0 else integral(
      function(r) q_closed(r, m)^2 * dnorm(r, sd = sqrt(m - 2)), -Inf, Inf
    )
    extra <- (m - 1) * cross_covariance
    data.frame(folds = m, mean_fixed_proxy = mean_proxy,
               pair_cross_covariance = cross_covariance,
               unconditional_variance = mean_proxy + extra,
               extra_variance_relative_to_proxy = extra / mean_proxy)
  }))
  stopifnot(result$pair_cross_covariance[[1L]] == 0,
            all(result$pair_cross_covariance[-1L] > 0))
  cat("optimizer_identity_max_error=", optimizer_error, "\n", sep = "")
  cat("q_formula_max_error=", q_error, "\n", sep = "")
  print(result, row.names = FALSE, digits = 12L)
  cat("No Monte Carlo draws or nuisance fits; this is not a RoCE CLT or an SE inflation rule.\n")
  invisible(result)
}

if (sys.nframe() == 0L) main()
