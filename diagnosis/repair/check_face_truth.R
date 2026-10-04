#!/usr/bin/env Rscript
# Audit the current binary DGP truth without changing its fixed calibration.
# Quadrature and fresh draws integrate only the four signal coordinates.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: check_face_truth.R INSTALLED_LIBRARY NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
kappa <- RoCE:::FACE_KAPPA
beta <- seq(.4, 1.2, length.out = 4L)
delta <- RoCE:::FACE_BINARY_ATE_TARGET
stopifnot(identical(RoCE:::get_face_outcome_parameters(100L)$beta_linear[1:4], beta),
          identical(RoCE:::get_face_outcome_parameters(100L)$beta_squared[1:4], beta),
          delta == 1, RoCE:::FACE_BINARY_SIGNAL_SD == 1)
references <- setNames(lapply(c("C1", "C2", "C3"), function(config) {
  RoCE:::.face_reference_population(100L, kappa, config,
    RoCE:::FACE_MISSPECIFICATION_STRENGTH, "binary")
}), c("C1", "C2", "C3"))
saveRDS(references, file.path(output, "fixed_reference_parameters.rds"))
reference_rows <- do.call(rbind, lapply(names(references), function(config) {
  r <- references[[config]]
  data.frame(config = config, mu0 = r$mu0_superpop, mu1 = r$mu1_superpop,
    tate = r$ate_superpop, eta_mean = r$calibration$eta_mean,
    eta_sd = r$calibration$eta_sd, n_reference = r$n_ref, reference_seed = r$ref_seed)
}))
write.csv(reference_rows, file.path(output, "reported_truth.csv"), row.names = FALSE)
stopifnot(identical(references$C1$calibration, references$C3$calibration),
          identical(references$C1$ate_superpop, references$C3$ate_superpop))

normal_quadrature <- function(order) {
  jacobi <- matrix(0, order, order)
  for (j in seq_len(order - 1L)) jacobi[j, j + 1L] <- jacobi[j + 1L, j] <- sqrt(j)
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(nodes = decomposition$values, weights = decomposition$vectors[1L, ]^2)
}
integrate_quadratic <- function(order, calibration) {
  rule <- normal_quadrature(order)
  z <- rule$nodes
  w <- rule$weights
  stopifnot(abs(sum(w) - 1) < 1e-12, abs(sum(w * z)) < 1e-12,
            abs(sum(w * z^2) - 1) < 1e-12)
  x <- kappa + z
  terms <- outer(x - kappa + x^2, beta)
  index <- expand.grid(first = seq_len(order), second = seq_len(order), third = seq_len(order))
  partial <- terms[index$first, 1L] + terms[index$second, 2L] + terms[index$third, 3L]
  partial_weight <- w[index$first] * w[index$second] * w[index$third]
  mu0 <- mu1 <- 0
  for (j in seq_len(order)) {
    g <- (partial + terms[j, 4L] - calibration$eta_mean) / calibration$eta_sd
    mu0 <- mu0 + w[j] * sum(partial_weight * plogis(g))
    mu1 <- mu1 + w[j] * sum(partial_weight * plogis(g + delta))
  }
  data.frame(order = order, mu0 = mu0, mu1 = mu1, tate = mu1 - mu0)
}
quadrature <- do.call(rbind, lapply(c(16L, 24L, 32L, 48L, 64L),
  integrate_quadratic, calibration = references$C3$calibration))
write.csv(quadrature, file.path(output, "c1_c3_quadrature.csv"), row.names = FALSE)

# Fresh independent integration, with the original transform/logit calibration
# held fixed. Recalibrating on these draws would define a different DGP.
set.seed(314159L)
sample_size <- 1000000L
batch_size <- 50000L
totals <- setNames(lapply(c("C1", "C2"), function(x) rep(0, 5L)), c("C1", "C2"))
for (batch in seq_len(sample_size %/% batch_size)) {
  x <- matrix(rnorm(batch_size * 4L), ncol = 4L) + kappa
  raw <- drop((sweep(x, 2L, kappa, "-") + x^2) %*% beta)
  transformed <- cbind(exp(x[, 1L] / 2), x[, 2L] / (1 + exp(x[, 1L])) + 10,
    (x[, 1L] * x[, 3L] / 25 + .6)^3, (x[, 2L] + x[, 4L] + 20)^2)
  standardization <- references$C2$standardization
  dagger <- sweep(sweep(transformed, 2L, standardization$mean, "-"),
                  2L, standardization$sd, "/") + kappa
  transformed_raw <- drop((sweep(dagger, 2L, kappa, "-") + dagger^2) %*% beta)
  for (config in names(totals)) {
    r <- references[[config]]
    eta <- if (config == "C1") raw else .25 * raw + .75 * transformed_raw
    g <- (eta - r$calibration$eta_mean) / r$calibration$eta_sd
    m0 <- plogis(g)
    m1 <- plogis(g + delta)
    difference <- m1 - m0
    totals[[config]] <- totals[[config]] + c(sum(m0), sum(m1), sum(difference),
      sum(difference^2), length(difference))
  }
}
independent <- do.call(rbind, lapply(names(totals), function(config) {
  x <- totals[[config]]
  n <- x[5L]
  average <- x[3L] / n
  mc_se <- sqrt((x[4L] - n * average^2) / (n - 1) / n)
  data.frame(config = config, n_integration = n, mu0 = x[1L] / n, mu1 = x[2L] / n,
    tate = average, integration_mc_se = mc_se,
    reported_minus_independent = references[[config]]$ate_superpop - average)
}))
write.csv(independent, file.path(output, "independent_integration.csv"), row.names = FALSE)
print(reference_rows, row.names = FALSE)
print(quadrature, row.names = FALSE)
print(independent, row.names = FALSE)
error <- references$C3$ate_superpop - tail(quadrature$tate, 1L)
convergence <- max(abs(as.numeric(quadrature[nrow(quadrature), -1L]) -
                      as.numeric(quadrature[nrow(quadrature) - 1L, -1L])))
writeLines(c("# Binary FACE truth audit", "",
  "Calibration and transformation constants are frozen at their original reference values.",
  "No study data, DGP parameters, sample sizes or historical result labels were changed.", "",
  sprintf("C1/C3 reported TATE: %.12f", references$C3$ate_superpop),
  sprintf("C1/C3 64-node-per-coordinate quadrature TATE: %.12f", tail(quadrature$tate, 1L)),
  sprintf("Reported minus quadrature: %.12g", error),
  sprintf("Largest last-two-order change across mu0/mu1/TATE: %.12g", convergence),
  "",
  "The quadrature convergence comparison is a numerical check, not a rigorous error bound.",
  "The independent million-draw integration uses only four outcome signal coordinates",
  "and has Monte Carlo uncertainty recorded in independent_integration.csv.",
  "The log-odds shift is 1; the marginal target risk difference is an integral, not 1.",
  "C1 and C3 share the same outcome law, so changing propensity alone does not change TATE."),
  file.path(output, "report.md"))
stopifnot(convergence < 1e-8)
writeLines("PASSED: numerical truth audit; study configuration unchanged", file.path(output, "status.txt"))
