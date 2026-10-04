#!/usr/bin/env Rscript
# Algebra, exact finite-support checks and a fixed-observation radius diagnostic.
# No nuisance fits or Monte Carlo coverage experiment are performed here.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L,
          startsWith(arguments[1L], "/scratch.global/zhan9381/FACE-HD/"),
          !dir.exists(arguments[1L]))
output_directory <- arguments[1L]
dir.create(output_directory, recursive = TRUE)
assertion_count <- 0L
verify <- function(condition) {
  stopifnot(isTRUE(condition))
  assertion_count <<- assertion_count + 1L
}
assert_close <- function(actual, expected, tolerance = 1e-9) {
  verify(max(abs(actual - expected)) <= tolerance * max(1, max(abs(expected))))
}
dispersion_upper <- function(observed, linear_bound, quadratic_bound, alpha) {
  multiplier <- (1 - alpha) / alpha
  pmax(0, observed + multiplier * linear_bound / 2 +
    sqrt(pmax(0, multiplier * linear_bound * observed +
      multiplier^2 * linear_bound^2 / 4 + multiplier * quadratic_bound)))
}
variance_bounds <- function(upper_variances, sample_sizes) {
  source_count <- length(upper_variances)
  if (source_count == 1L) return(c(linear = 0, quadratic = 0))
  correction <- sample_sizes / (sample_sizes - 1)
  c(linear = 4 * max(upper_variances) / source_count,
    quadratic = 2 / source_count^4 * (sum(upper_variances)^2 +
      sum((correction * (source_count - 1)^2 - 1) * upper_variances^2)))
}

set.seed(20261002)
for (source_count in c(1L, 2L, 4L, 8L, 20L, 64L)) {
  for (valid_count in unique(c(1L, ceiling(source_count / 2), source_count))) {
    for (repeat_id in seq_len(20L)) {
      bias_budget <- runif(1, 0, .05)
      biases <- c(runif(valid_count, -bias_budget, bias_budget),
                  rnorm(source_count - valid_count))
      dispersion <- mean((biases - mean(biases))^2)
      bound <- bias_budget + sqrt((source_count - valid_count) / valid_count * dispersion)
      verify(abs(mean(biases)) <= bound + 1e-12)
      if (valid_count < source_count) {
        sharp_biases <- c(rep(bias_budget, valid_count), rep(bias_budget + .4, source_count - valid_count))
        assert_close(mean(sharp_biases), bias_budget + sqrt((source_count - valid_count) /
          valid_count * mean((sharp_biases - mean(sharp_biases))^2)))
      }
    }
  }
}
# The moment condition is a strict relaxation of a count condition.
population_centers <- c(0, 0, 1, 1)
candidate_mean <- .5
verify(abs(mean(population_centers) - candidate_mean) <= sqrt(mean((population_centers - .5)^2)))
verify(sum(population_centers == candidate_mean) < 2L)

for (linear_bound in c(0, 1e-4, 1, 20)) for (quadratic_bound in c(0, 1e-6, .5, 10)) {
  for (energy in c(0, .01, 1, 100)) for (alpha in c(.01, .05, .2)) {
    boundary <- energy - sqrt((1 - alpha) / alpha * (linear_bound * energy + quadratic_bound))
    verify(all(dispersion_upper(boundary + c(0, .1, 1), linear_bound,
                                quadratic_bound, alpha) >= energy - 1e-8))
  }
}

# Three independent sites, two iid patient scores per site, four states per
# patient. Both arm scores have within-patient dependence: 4^6 = 4096 cases.
source_count <- 3L
sample_sizes <- rep(2L, source_count)
supports <- list(rbind(c(-.8, 0), c(.6, 0), c(0, -.5), c(0, 1.2)),
                 rbind(c(-.6, 0), c(.9, 0), c(0, -.7), c(0, 1.1)),
                 rbind(c(-.7, 0), c(.8, 0), c(0, -.9), c(0, .6)))
probabilities <- list(c(.3, .2, .1, .4), c(.2, .3, .3, .2), c(.1, .4, .2, .3))
means <- Map(function(x, probability) colSums(x * probability), supports, probabilities)
covariances <- Map(function(x, probability, expectation) {
  centered <- sweep(x, 2L, expectation, "-")
  crossprod(centered, centered * probability)
}, supports, probabilities, means)
verify(any(abs(vapply(covariances, function(x) x[1L, 2L], numeric(1L))) > 1e-3))
site_cases <- lapply(seq_len(source_count), function(site) {
  pairs <- expand.grid(first = 1:4, second = 1:4)
  first <- supports[[site]][pairs$first, , drop = FALSE]
  second <- supports[[site]][pairs$second, , drop = FALSE]
  list(mean = (first + second) / 2, variance = (first - second)^2 / 4,
       product = first * second,
       probability = probabilities[[site]][pairs$first] * probabilities[[site]][pairs$second])
})
configurations <- expand.grid(rep(list(1:16), source_count))
probability <- Reduce(`*`, lapply(seq_len(source_count), function(site)
  site_cases[[site]]$probability[configurations[[site]]]))
assert_close(sum(probability), 1)
centering <- diag(source_count) / source_count - matrix(1 / source_count^2, source_count, source_count)
enumeration_summary <- vector("list", 2L)
for (arm in 1:2) {
  sample_means <- sapply(seq_len(source_count), function(site)
    site_cases[[site]]$mean[configurations[[site]], arm])
  unbiased_variances <- sapply(seq_len(source_count), function(site)
    site_cases[[site]]$variance[configurations[[site]], arm])
  statistic <- rowMeans((sample_means - rowMeans(sample_means))^2) -
    (source_count - 1) / source_count^2 * rowSums(unbiased_variances)
  direct <- Reduce(`+`, lapply(seq_len(source_count), function(site)
    centering[site, site] * site_cases[[site]]$product[configurations[[site]], arm]))
  for (site in seq_len(source_count)) for (other in seq_len(source_count)) if (site != other) {
    direct <- direct + centering[site, other] * sample_means[, site] * sample_means[, other]
  }
  assert_close(statistic, direct)
  expectation <- vapply(means, `[`, numeric(1L), arm)
  variances <- vapply(covariances, function(x) x[arm, arm], numeric(1L)) / sample_sizes
  dispersion <- mean((expectation - mean(expectation))^2)
  assert_close(sum(probability * statistic), dispersion)
  exact_variance <- 4 * sum((centering %*% expectation)^2 * variances) +
    2 * sum((centering^2 * outer(variances, variances))[row(centering) != col(centering)]) +
    2 * sum(sample_sizes / (sample_sizes - 1) * diag(centering)^2 * variances^2)
  assert_close(sum(probability * (statistic - dispersion)^2), exact_variance)
  bounds <- variance_bounds(variances * 1.2, sample_sizes)
  verify(exact_variance <= bounds["linear"] * dispersion + bounds["quadratic"] + 1e-12)
  for (alpha in c(.01, .05, .2)) {
    upper <- dispersion_upper(statistic, bounds["linear"], bounds["quadratic"], alpha)
    verify(sum(probability[upper >= dispersion - 1e-12]) >= 1 - alpha)
  }
  enumeration_summary[[arm]] <- data.frame(arm = arm, dispersion = dispersion,
    exact_variance = exact_variance, variance_bound = bounds["linear"] * dispersion + bounds["quadratic"])
}

# Joint-ellipsoid witnesses imply simultaneous protection over the entire
# borrowing family. A dense grid checks the implementation of that identity.
contrast_map <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
for (repeat_id in seq_len(100L)) {
  matrix_draw <- matrix(rnorm(16), 4L)
  upper_covariance <- crossprod(matrix_draw) + diag(.2, 4L)
  direction <- rnorm(4L)
  critical <- 3
  noise <- drop(t(chol(upper_covariance)) %*% direction) / sqrt(sum(direction^2)) * critical * runif(1)
  budgets <- runif(4L, 0, .2)
  drift <- budgets * runif(4L, -1, 1)
  treatment_means <- c(.6, .3)
  summaries <- rep(treatment_means, 2L) + drift + noise
  verify(drop(crossprod(noise, solve(upper_covariance, noise))) <= critical^2 + 1e-10)
  contrast_covariance <- contrast_map %*% upper_covariance %*% t(contrast_map)
  contrast_estimates <- drop(contrast_map %*% summaries)
  contrast_budgets <- c(sum(budgets[1:2]), sum(budgets[3:4]))
  fractions <- seq(0, 1, length.out = 401L)
  coefficients <- cbind(1 - fractions, fractions)
  radii <- critical * sqrt(rowSums((coefficients %*% contrast_covariance) * coefficients)) +
    drop(coefficients %*% contrast_budgets)
  true_effect <- treatment_means[1L] - treatment_means[2L]
  verify(all(abs(drop(coefficients %*% contrast_estimates) - true_effect) <= radii + 1e-10))
  verify(min(radii) <= radii[1L] + 1e-12)
}

# Oracle Gaussian formula diagnostic at observed dispersion zero. These are
# radii at one possible observation, not average lengths or coverage rates.
gaussian_dispersion_upper <- function(observed, variance, source_count, alpha) {
  degrees <- source_count - 1L
  statistic <- max(0, source_count * observed / variance + degrees)
  if (pchisq(statistic, degrees) <= alpha) return(0)
  objective <- function(noncentrality) pchisq(statistic, degrees, ncp = noncentrality) - alpha
  upper <- 1
  while (objective(upper) > 0) upper <- 2 * upper
  variance / source_count * uniroot(objective, c(0, upper), tol = 1e-10)$root
}
radius_table <- do.call(rbind, lapply(c(250L, 1000L), function(evaluation_size) {
  do.call(rbind, lapply(c(4L, 8L, 20L, 64L, 200L, 1024L), function(source_count) {
    valid_count <- ceiling(source_count / 2)
    variance <- .5 / evaluation_size
    bounds <- variance_bounds(rep(variance, source_count), rep(evaluation_size, source_count))
    upper_cantelli <- dispersion_upper(0, bounds["linear"], bounds["quadratic"], .0125)
    upper_gaussian <- gaussian_dispersion_upper(0, variance, source_count, .0125)
    noise <- sqrt(qchisq(.975, 4L) * 2 * variance / source_count)
    data.frame(evaluation_per_site = evaluation_size, source_count = source_count,
      valid_count = valid_count, source_noise_radius = noise,
      cantelli_source_radius = noise + 2 * sqrt((source_count - valid_count) / valid_count * upper_cantelli),
      gaussian_source_radius = noise + 2 * sqrt((source_count - valid_count) / valid_count * upper_gaussian),
      full_target_normal_95_radius = qnorm(.975) * sqrt(1 / 1000))
  }))
}))
write.csv(radius_table, file.path(output_directory, "ideal_radius_diagnostic.csv"), row.names = FALSE)
write.csv(do.call(rbind, enumeration_summary), file.path(output_directory, "finite_support_variance.csv"), row.names = FALSE)
summary <- data.frame(assertions = assertion_count, enumerated_score_configurations = nrow(configurations),
                      joint_ellipsoid_witnesses = 100L)
write.csv(summary, file.path(output_directory, "summary.csv"), row.names = FALSE)
writeLines("CENTER_DISPERSION_CHECKS_PASSED", file.path(output_directory, "CHECKS_PASSED"))
print(summary)
print(radius_table[radius_table$source_count %in% c(20L, 200L), ], row.names = FALSE)
