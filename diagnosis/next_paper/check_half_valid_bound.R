#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !startsWith(arguments[[1L]], "/scratch.global/zhan9381/FACE-HD/")) {
  stop("Supply an output directory on FACE-HD scratch")
}
output <- arguments[[1L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
settings <- expand.grid(num_sources = c(4L, 16L, 64L, 256L), valid_fraction = c(.25, .5),
                        shared_correlation = c(0, .5, .9))
alpha <- .05
sample_size <- 1000L
distance <- 1
records <- list()
for (index in seq_len(nrow(settings))) {
  count <- settings$num_sources[index]
  valid_minimum <- as.integer(count * settings$valid_fraction[index])
  correlation <- settings$shared_correlation[index]
  covariance <- (correlation * matrix(1, count + 1L, count + 1L) +
                  (1 - correlation) * diag(count + 1L)) / sample_size
  anchor_variance <- covariance[1L, 1L]
  delta <- distance * sqrt(anchor_variance)
  loading <- covariance[-1L, 1L] / anchor_variance
  first_valid <- seq_len(valid_minimum)
  second_valid <- seq.int(valid_minimum + 1L, 2L * valid_minimum)
  mean_zero <- numeric(count + 1L)
  mean_zero[second_valid + 1L] <- (1 - loading[second_valid]) * delta
  difference <- delta * covariance[, 1L] / anchor_variance
  mean_one <- mean_zero + difference
  mahalanobis <- as.numeric(crossprod(difference, solve(covariance, difference)))
  stopifnot(abs(mean_zero[1L]) < 1e-14, abs(mean_one[1L] - delta) < 1e-14,
            all(abs(mean_zero[first_valid + 1L]) < 1e-14),
            all(abs(mean_one[second_valid + 1L] - delta) < 1e-14),
            abs(mahalanobis - distance^2) < 1e-10)
  total_variation <- 2 * pnorm(sqrt(mahalanobis) / 2) - 1
  lower_bound <- delta * max(0, 1 - 2 * alpha - total_variation)
  oracle_indices <- c(1L, second_valid + 1L)
  oracle_covariance <- covariance[oracle_indices, oracle_indices]
  oracle_se <- sqrt(1 / sum(solve(oracle_covariance)))
  records[[index]] <- data.frame(settings[index, ], valid_minimum = valid_minimum,
    delta = delta, mahalanobis_squared = mahalanobis, total_variation = total_variation,
    expected_length_lower_bound = lower_bound,
    oracle_normal_length = 2 * qnorm(1 - alpha / 2) * oracle_se)
}
records <- do.call(rbind, records)
write.csv(records, file.path(output, "bound_checks.csv"), row.names = FALSE)
# A covariance with unequal variances and signed target correlations checks
# that the construction is not restricted to the common-factor example.
set.seed(92501)
factor <- matrix(rnorm(49), 7, 7)
covariance <- (tcrossprod(factor) + diag(7)) / sample_size
difference <- sqrt(covariance[1, 1]) * covariance[, 1] / covariance[1, 1]
stopifnot(abs(as.numeric(crossprod(difference, solve(covariance, difference))) - 1) < 1e-12)
saveRDS(list(settings = settings, records = records, arbitrary_covariance = covariance),
        file.path(output, "bound_checks.rds"))
writeLines("HALF_VALID_BOUND_CHECKS_PASSED: 24 factor cases and one arbitrary covariance.",
           file.path(output, "COMPLETE"))
cat("HALF_VALID_BOUND_CHECKS_PASSED\n")
