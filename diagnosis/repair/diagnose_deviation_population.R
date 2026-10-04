#!/usr/bin/env Rscript
# Population source-outcome discrepancy under the frozen bounded DGP.
# This is a quadrature diagnostic, not a simulation of the fitted estimator.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: library_directory output_directory")
.libPaths(c(arguments[[1L]], .libPaths()))
suppressPackageStartupMessages(library(RoCE))
output <- arguments[[2L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Output must be on FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

evaluate_population <- function(order) {
  grid <- RoCE:::.bounded_quadrature(order)
  stopifnot(abs(sum(grid$weights) - 1) < 1e-12)
  rows <- list()
  for (scenario in c("C1", "C2", "C3")) {
    strength <- RoCE:::.face_misspecification_strengths(scenario, .75)$outcome
    predictor <- drop(RoCE:::.bounded_active_features(grid$X, strength) %*% c(.25, .20, .10, .05))
    target_mu0 <- sum(grid$weights * plogis(-.5 + predictor))
    target_mu1 <- sum(grid$weights * plogis(.5 + predictor))
    for (deviation in c(0, .5, 1, 2)) {
      source_mu1_on_target <- sum(grid$weights * plogis(.5 + predictor + deviation))
      discrepancy <- source_mu1_on_target - target_mu1
      rows[[length(rows) + 1L]] <- data.frame(
        config = scenario, rho = deviation, target_mu0 = target_mu0, target_mu1 = target_mu1,
        target_tate = target_mu1 - target_mu0, source_mu1_on_target = source_mu1_on_target,
        source_control_discrepancy = 0, source_treated_discrepancy = discrepancy,
        fixed_weight_point24_bias = .24 * discrepancy)
    }
  }
  do.call(rbind, rows)
}

coarse <- evaluate_population(24L)
fine <- evaluate_population(32L)
maximum_error <- max(abs(as.matrix(coarse[, -(1:2)]) - as.matrix(fine[, -(1:2)])))
stopifnot(maximum_error < 1e-10)
write.csv(fine, file.path(output, "population_discrepancies.csv"), row.names = FALSE)
writeLines(c(
  sprintf("Maximum quadrature difference (orders 24 and 32): %.3g", maximum_error),
  "The calculation uses the four active target coordinates, so it does not depend on p >= 4.",
  "An assisted mean with oracle merged weights targets E_target[m_source,a(X)].",
  "It still differs from the target mean when that source's outcome model is nontransportable.",
  "The .24-weight column is a fixed-weight illustration, not the bias of the learned-weight estimator."
), file.path(output, "checks.txt"))
print(fine)
cat("Maximum quadrature difference:", maximum_error, "\n")
