#!/usr/bin/env Rscript
# Deterministic algebra and exact finite-support checks for the design note.
# This is not a fitted RoCE-K implementation or a Monte Carlo coverage study.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L, startsWith(arguments[1L], "/scratch.global/zhan9381/FACE-HD/"),
          !dir.exists(arguments[1L]))
output <- arguments[1L]
dir.create(output, recursive = TRUE)
checks <- 0L
verify <- function(condition) {
  stopifnot(isTRUE(condition))
  checks <<- checks + 1L
}
close <- function(actual, expected, tolerance = 1e-9) {
  verify(max(abs(actual - expected)) <= tolerance * max(1, max(abs(expected))))
}
block_diagonal <- function(blocks) {
  sizes <- vapply(blocks, nrow, integer(1L))
  result <- matrix(0, sum(sizes), sum(sizes))
  offset <- 0L
  for (block in blocks) {
    index <- offset + seq_len(nrow(block))
    result[index, index] <- block
    offset <- offset + nrow(block)
  }
  result
}
assembly <- function(K) {
  dimension <- 2L * (K + 1L)
  L <- matrix(0, dimension, 4L + 2L * K)
  L[1L, 1L] <- L[K + 2L, 2L] <- 1
  for (j in seq_len(K)) {
    L[j + 1L, 3L] <- L[K + 2L + j, 4L] <- 1
    L[j + 1L, 4L + 2L * j - 1L] <- 1
    L[K + 2L + j, 4L + 2L * j] <- 1
  }
  list(L = L, H = cbind(c(rep(1, K + 1L), rep(0, K + 1L)),
                        c(rep(0, K + 1L), rep(1, K + 1L))))
}
metric <- function(W, H, known) {
  inverse <- solve(W)
  coefficients <- drop(inverse %*% H %*% solve(crossprod(H, inverse %*% H), c(1, -1)))
  precision <- inverse - inverse %*% H %*% solve(crossprod(H, inverse %*% H), t(H) %*% inverse)
  W_known <- W[known, known, drop = FALSE]
  H_known <- H[known, , drop = FALSE]
  known_inverse <- solve(W_known)
  reference <- matrix(0, nrow(W), ncol(W))
  reference[known, known] <- known_inverse - known_inverse %*% H_known %*%
    solve(crossprod(H_known, known_inverse %*% H_known), t(H_known) %*% known_inverse)
  list(a = coefficients, D = precision - reference, R = reference,
       variance = drop(crossprod(coefficients, W %*% coefficients)),
       df = nrow(W) - length(known))
}
energy_upper <- function(value, df, kappa, alpha) {
  if (df == 0L) return(rep(0, length(value)))
  cD <- (1 - alpha) / alpha
  pmax(0, value + 2 * cD + sqrt(pmax(0, 4 * cD * value + 4 * cD^2 + 2 * cD * kappa * df)))
}

# Verify inversion on the whole concentration event, including negative
# observations, near-zero energy and large noncentrality.
for (energy in c(0, 1e-8, .1, 1, 50, 10000)) for (df in c(1, 4, 32)) {
  for (kappa in c(1.002, 2)) for (alpha in c(.01, .05, .2)) {
    boundary <- energy - sqrt((1 - alpha) / alpha * (4 * energy + 2 * kappa * df))
    values <- boundary + c(0, .01, 1, energy)
    verify(all(energy_upper(values, df, kappa, alpha) + 1e-8 * max(1, energy) >= energy))
  }
}
close(energy_upper(c(-1, 0, 1), 0, 2, .05), c(0, 0, 0))

# The Bernstein radius solves its two-sided tail inequality exactly.
for (variance_bound in c(1e-6, .2, 10)) for (maximum in c(0, .01, 4)) {
  for (alpha in c(.01, .05, .2)) {
    threshold <- log(2 / alpha)
    radius <- maximum * threshold / 3 +
      sqrt(2 * variance_bound * threshold + (maximum * threshold / 3)^2)
    close(radius^2 / (2 * (variance_bound + maximum * radius / 3)), threshold)
  }
}
exact_variance <- function(L, D, means, covariances, sizes) {
  Sigma <- L %*% block_diagonal(Map(`/`, covariances, sizes)) %*% t(L)
  mean <- drop(L %*% unlist(means, use.names = FALSE))
  A <- t(L) %*% D %*% L
  indices <- split(seq_len(ncol(L)), rep(seq_along(sizes), vapply(means, length, integer(1L))))
  quadratic <- 0
  for (b in seq_along(sizes)) for (c in seq_along(sizes)) {
    Abc <- A[indices[[b]], indices[[c]], drop = FALSE]
    denominator <- if (b == c) sizes[b] * (sizes[b] - 1) else sizes[b] * sizes[c]
    quadratic <- quadratic + 2 * sum(diag(Abc %*% covariances[[c]] %*%
      t(Abc) %*% covariances[[b]])) / denominator
  }
  as.numeric(4 * crossprod(mean, D %*% Sigma %*% D %*% mean) + quadratic)
}

set.seed(20261001)
for (K in c(1L, 2L, 4L, 8L)) for (repeat_id in seq_len(15L)) {
  map <- assembly(K)
  dimensions <- c(4L, rep(2L, K))
  sizes <- sample(2:9, K + 1L, replace = TRUE)
  covariances <- lapply(dimensions, function(d) {
    x <- matrix(rnorm(d * d), d, d)
    crossprod(x) / d + diag(.1, d)
  })
  upper <- lapply(covariances, function(x) x + diag(.15, nrow(x)))
  Sigma <- map$L %*% block_diagonal(Map(`/`, covariances, sizes)) %*% t(map$L)
  W <- map$L %*% block_diagonal(Map(`/`, upper, sizes)) %*% t(map$L)
  known <- if (repeat_id %% 2L) c(1L, K + 2L) else c(1L, seq.int(K + 2L, 2L * (K + 1L)))
  fit <- metric(W, map$H, known)
  close(crossprod(map$H, fit$a), c(1, -1))
  close(fit$D %*% map$H, matrix(0, nrow(W), 2L))
  close(fit$D %*% W %*% fit$D, fit$D)
  close(sum(diag(fit$D %*% W)), fit$df)
  verify(min(eigen(W - Sigma, symmetric = TRUE, only.values = TRUE)$values) > -1e-10)
  means <- lapply(dimensions, function(d) rnorm(d))
  mean <- drop(map$L %*% unlist(means, use.names = FALSE))
  energy <- drop(crossprod(mean, fit$D %*% mean))
  variance <- exact_variance(map$L, fit$D, means, covariances, sizes)
  verify(variance <= 4 * energy + 2 * max(sizes / (sizes - 1)) * fit$df + 1e-8)
  valid <- sort(unique(c(known, 2L, 2L * (K + 1L))))
  restricted <- metric(W[valid, valid], map$H[valid, , drop = FALSE], seq_along(valid))
  d <- numeric(nrow(W)); d[valid] <- restricted$a
  projected <- drop((diag(nrow(W)) - fit$R %*% W) %*% d)
  gap <- drop(crossprod(projected, W %*% projected)) - fit$variance
  close(projected - fit$a, fit$D %*% W %*% d)
  close(gap, crossprod(d, W %*% fit$D %*% W %*% d))
  budget <- abs(mean[valid] - drop(map$H[valid, , drop = FALSE] %*% c(.5, .3)))
  drift <- mean - drop(map$H %*% c(.5, .3))
  bound <- sum(abs(projected[valid]) * budget) + sqrt(max(0, gap * energy))
  verify(abs(sum(fit$a * drift)) <= bound + 1e-9)
}

# An actual randomized binary-observation example. Two sources have different
# valid arm sets. Patient arms are dependent and treatment denominators are fixed.
K <- 2L; sizes <- c(2L, 2L, 2L); map <- assembly(K)
states <- expand.grid(X = 0:1, A = 0:1, Y = 0:1)
treatment_probability <- c(.45, .55, .35)
departures <- list(c(0, 0), c(0, .10), c(.10, 0))
common <- cbind(.42 + .06 * states$X, .20 + .10 * states$X)
supports <- probabilities <- vector("list", 3L)
for (b in seq_len(3L)) {
  e <- treatment_probability[b]
  pA <- ifelse(states$A == 1, e, 1 - e)
  pY <- ifelse(states$A == 1, .50 + .10 * states$X + departures[[b]][1],
                .25 + .10 * states$X + departures[[b]][2])
  probabilities[[b]] <- .5 * pA * ifelse(states$Y == 1, pY, 1 - pY)
  residual <- cbind((states$A == 1) / e * (states$Y - common[, 1]),
                    (states$A == 0) / (1 - e) * (states$Y - common[, 2]))
  anchor <- cbind(.40 + .12 * states$X, .28 + .04 * states$X)
  anchor[, 1] <- anchor[, 1] + (states$A == 1) / e * (states$Y - anchor[, 1])
  anchor[, 2] <- anchor[, 2] + (states$A == 0) / (1 - e) * (states$Y - anchor[, 2])
  supports[[b]] <- if (b == 1L) cbind(anchor, common) else residual
}
means <- Map(function(x, p) colSums(x * p), supports, probabilities)
covariances <- Map(function(x, p, mu) {
  centered <- sweep(x, 2, mu, "-")
  crossprod(centered, centered * p)
}, supports, probabilities, means)

# Check range domination and row-wise moment-error inflation, including an
# added known-constant coordinate. These are matrix certificates, not tests
# of the probability of the pilot concentration event.
for (b in seq_len(3L)) {
  support <- cbind(supports[[b]], constant = 1)
  center <- (apply(support, 2, max) + apply(support, 2, min)) / 2
  centered_support <- sweep(support, 2, center, "-")
  ranges <- apply(abs(centered_support), 2, max)
  dimension <- ncol(support)
  population_covariance <- matrix(0, dimension, dimension)
  population_covariance[-dimension, -dimension] <- covariances[[b]]
  range_upper <- dimension * diag(ranges^2)
  verify(min(eigen(range_upper - population_covariance, symmetric = TRUE,
                   only.values = TRUE)$values) > -1e-10)
  pilot_mean <- colMeans(centered_support)
  pilot_covariance <- crossprod(sweep(centered_support, 2, pilot_mean, "-")) / nrow(support)
  error <- abs(population_covariance - pilot_covariance)
  moment_upper <- pilot_covariance + diag(rowSums(error))
  verify(min(eigen(moment_upper - population_covariance, symmetric = TRUE,
                   only.values = TRUE)$values) > -1e-10)
  close(moment_upper[dimension, ], numeric(dimension))
  close(range_upper[dimension, ], numeric(dimension))
}
Sigma <- map$L %*% block_diagonal(Map(`/`, covariances, sizes)) %*% t(map$L)
upper <- lapply(covariances, function(x) x + diag(.02, nrow(x)))
W <- map$L %*% block_diagonal(Map(`/`, upper, sizes)) %*% t(map$L)
fit <- metric(W, map$H, c(1L, 4L))
population_mean <- drop(map$L %*% unlist(means, use.names = FALSE))
energy <- drop(crossprod(population_mean, fit$D %*% population_mean))
A <- t(map$L) %*% fit$D %*% map$L
indices <- list(1:4, 5:6, 7:8)
site_cases <- lapply(seq_len(3L), function(b) {
  pairs <- expand.grid(first = seq_len(nrow(states)), second = seq_len(nrow(states)))
  first <- supports[[b]][pairs$first, , drop = FALSE]
  second <- supports[[b]][pairs$second, , drop = FALSE]
  difference <- first - second
  block <- A[indices[[b]], indices[[b]], drop = FALSE]
  list(mean = (first + second) / 2,
       probability = probabilities[[b]][pairs$first] * probabilities[[b]][pairs$second],
       correction = rowSums((difference %*% block) * difference) / 4,
       within_pair = rowSums((first %*% block) * second))
})
grid <- expand.grid(rep(list(seq_len(nrow(states)^2)), 3L))
sample_means <- lapply(seq_len(3L), function(b) site_cases[[b]]$mean[grid[[b]], , drop = FALSE])
Z <- do.call(cbind, sample_means) %*% t(map$L)
values <- rowSums((Z %*% fit$D) * Z) - Reduce(`+`, lapply(seq_len(3L),
  function(b) site_cases[[b]]$correction[grid[[b]]]))
direct <- Reduce(`+`, lapply(seq_len(3L), function(b) site_cases[[b]]$within_pair[grid[[b]]]))
for (b in seq_len(3L)) for (c in seq_len(3L)) if (b != c) {
  direct <- direct + rowSums((sample_means[[b]] %*% A[indices[[b]], indices[[c]], drop = FALSE]) * sample_means[[c]])
}
close(values, direct)
probability <- Reduce(`*`, lapply(seq_len(3L), function(b) site_cases[[b]]$probability[grid[[b]]]))
close(sum(probability), 1)
close(sum(values * probability), energy)
variance <- sum((values - energy)^2 * probability)
close(variance, exact_variance(map$L, fit$D, means, covariances, sizes))
verify(variance <= 4 * energy + 2 * max(sizes / (sizes - 1)) * fit$df)
for (alpha in c(.01, .05, .2)) {
  probability_below <- sum(probability[energy <= energy_upper(values, fit$df, max(sizes / (sizes - 1)), alpha)])
  verify(probability_below >= 1 - alpha - 1e-12)
}
summary <- data.frame(checks = checks, matrix_cases = 60L, enumerated_observed_data_configurations = nrow(grid),
  energy = energy, exact_variance = variance,
  variance_upper = 4 * energy + 2 * max(sizes / (sizes - 1)) * fit$df)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
saveRDS(list(statistic = values, probabilities = probability, summary = summary),
        file.path(output, "enumerated_reference.rds"))
writeLines(c("JOINT_QUADRATIC_IDENTITIES_PASSED", "Algebra and finite-support enumeration only; no fitted-method coverage claim."),
  file.path(output, "CHECKS_PASSED"))
print(summary)
