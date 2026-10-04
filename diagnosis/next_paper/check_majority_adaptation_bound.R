#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !startsWith(arguments[[1L]], "/scratch.global/zhan9381/FACE-HD/")) {
  stop("Supply an output directory on FACE-HD scratch")
}
output <- arguments[[1L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
mixture_chisquare <- function(shift, probability) {
  ratio <- probability / (1 - probability)
  probability^2 * expm1(shift^2) + 2 * probability * (1 - probability) * expm1(-ratio * shift^2) +
    (1 - probability)^2 * expm1(ratio^2 * shift^2)
}
quadrature_checks <- list()
for (fraction in c(.25, .5, .75)) {
  probability <- (1 + fraction) / 2
  ratio <- probability / (1 - probability)
  for (shift in c(.02, .05, .1)) {
    integrand <- function(x) {
      likelihood_difference <- probability * expm1(shift * x - shift^2 / 2) +
        (1 - probability) * expm1(-ratio * shift * x - ratio^2 * shift^2 / 2)
      likelihood_difference^2 * dnorm(x)
    }
    numerical <- integrate(integrand, -12, 12, rel.tol = 1e-9, abs.tol = 1e-13)$value
    exact <- mixture_chisquare(shift, probability)
    stopifnot(exact >= 0, abs(numerical - exact) <= 1e-10 * max(1, exact))
    quadrature_checks[[length(quadrature_checks) + 1L]] <- data.frame(valid_fraction = fraction,
      shift = shift, exact = exact, numerical = numerical)
  }
}
settings <- expand.grid(num_sources = c(4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 1000000L),
                        valid_fraction = c(.25, .5, .75))
alpha <- .05
sample_size <- 1000L
records <- list()
for (index in seq_len(nrow(settings))) {
  count <- settings$num_sources[index]
  fraction <- settings$valid_fraction[index]
  probability <- (1 + fraction) / 2
  ratio <- probability / (1 - probability)
  shift <- .5 / sqrt(ratio) * count^(-1/4)
  source_divergence <- mixture_chisquare(shift, probability)
  joint_divergence <- expm1(shift^2 + count * log1p(source_divergence))
  invalid_count_probability <- pbinom(ceiling(fraction * count) - 1L, count, probability)
  total_variation_bound <- min(1, .5 * sqrt(joint_divergence) + invalid_count_probability)
  length_bound <- shift / sqrt(sample_size) * max(0, 1 - 2 * alpha - total_variation_bound)
  stopifnot(is.finite(length_bound), source_divergence >= 0)
  records[[index]] <- data.frame(settings[index, ], total_variation_bound = total_variation_bound,
    prior_conditioning_error = invalid_count_probability, length_lower_bound = length_bound,
    all_valid_oracle_length = 2 * qnorm(.975) / sqrt(sample_size * (count + 1)),
    scaled_length_bound = length_bound * sqrt(sample_size) * count^(1/4))
}
records <- do.call(rbind, records)
stopifnot(all(records$length_lower_bound > 0),
          all(subset(records, num_sources == 1000000)$length_lower_bound >
                subset(records, num_sources == 1000000)$all_valid_oracle_length))
write.csv(do.call(rbind, quadrature_checks), file.path(output, "quadrature_checks.csv"), row.names = FALSE)
write.csv(records, file.path(output, "length_bound_checks.csv"), row.names = FALSE)
writeLines("MAJORITY_ADAPTATION_BOUND_CHECKS_PASSED: 9 density-integral checks and 24 finite-K bounds.",
           file.path(output, "COMPLETE"))
cat("MAJORITY_ADAPTATION_BOUND_CHECKS_PASSED\n")
