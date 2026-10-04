#!/usr/bin/env Rscript
# Numerical case-weight derivatives provide a check independent of the coded IF.
.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
set.seed(281)
n <- 300L
x <- matrix(rnorm(n * 3L), n, 3L)
a <- rbinom(n, 1L, plogis(.1 + .6 * x[, 1L]))
outcomes <- list(
  binomial = rbinom(n, 1L, plogis(-.2 + .4 * a + .5 * x[, 2L])),
  gaussian = -.2 + .4 * a + .5 * x[, 2L] + rnorm(n)
)
design <- cbind(1, x)
weighted_estimate <- function(case_weights, A_val, y, family) {
  propensity_fit <- glm.fit(design, a, weights = case_weights, family = binomial(),
                            control = glm.control(epsilon = 1e-12))
  arm <- a == A_val
  outcome_fit <- glm.fit(design[arm, ], y[arm], weights = case_weights[arm],
                         family = if (family == "binomial") binomial() else gaussian(),
                         control = glm.control(epsilon = 1e-12))
  propensity <- plogis(drop(design %*% propensity_fit$coefficients))
  predictor <- drop(design %*% outcome_fit$coefficients)
  outcome <- if (family == "binomial") plogis(predictor) else predictor
  probability <- if (A_val == 1L) propensity else 1 - propensity
  pseudo_outcome <- outcome + arm * (y - outcome) / probability
  list(estimate = weighted.mean(pseudo_outcome, case_weights),
       outcome = outcome, propensity = propensity)
}
settings <- expand.grid(A_val = 0:1, family = names(outcomes), stringsAsFactors = FALSE)
rows <- lapply(seq_len(nrow(settings)), function(setting) {
  A_val <- settings$A_val[setting]
  family <- settings$family[setting]
  y <- outcomes[[family]]
  fit <- weighted_estimate(rep(1, n), A_val, y, family)
  actual <- RoCE:::calculate_aipw_influence(y, a, x, fit$outcome, fit$propensity,
                                          A_val = A_val, family = family)
  stopifnot(abs(actual$estimate - fit$estimate) < 1e-12)
  base_influence <- actual$phi - actual$estimate
  plus_correction <- 2 * base_influence - actual$influence
  indices <- seq(1L, n, by = 13L)
  derivative <- vapply(indices, function(index) {
    step <- 1e-4
    upper <- lower <- rep(1, n)
    upper[index] <- 1 + step
    lower[index] <- 1 - step
    n * (weighted_estimate(upper, A_val, y, family)$estimate -
           weighted_estimate(lower, A_val, y, family)$estimate) / (2 * step)
  }, numeric(1L))
  data.frame(family, A_val, index = indices, numerical_derivative = derivative,
             implemented = actual$influence[indices],
             plus_correction = plus_correction[indices])
})
rows <- do.call(rbind, rows)
write.csv(rows, args[[1L]], row.names = FALSE)
print(aggregate(cbind(implemented_error = abs(rows$implemented - rows$numerical_derivative),
                      plus_error = abs(rows$plus_correction - rows$numerical_derivative)),
                by = list(family = rows$family, A_val = rows$A_val), FUN = max))
stopifnot(max(abs(rows$plus_correction - rows$numerical_derivative)) < 1e-5)
