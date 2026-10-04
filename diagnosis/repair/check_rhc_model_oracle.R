#!/usr/bin/env Rscript
# Independent oracle-only sampling check; never replaces fitted pilot repeats.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
root <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(startsWith(output, paste0(root, "/reviews/")), !dir.exists(output))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
template <- readRDS(configuration$template)
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
n <- template$observed_structure$t$n
records <- list(); summaries <- list()
for (scenario in c("O_case_mix", "W_tail")) {
  population <- template$populations[[scenario]]
  target <- population$sites$t
  delta <- target$outcome[, 2L] - target$outcome[, 1L]
  variance <- sum(target$probability * ((delta - population$truth)^2 +
    target$outcome[, 2L] * (1 - target$outcome[, 2L]) / target$propensity +
    target$outcome[, 1L] * (1 - target$outcome[, 1L]) / (1 - target$propensity))) / n
  draw <- function(seed) {
    set.seed(seed)
    index <- sample.int(length(target$probability), n, replace = TRUE, prob = target$probability)
    a <- rbinom(n, 1L, target$propensity[index])
    y <- rbinom(n, 1L, target$outcome[cbind(index, a + 1L)])
    m0 <- target$outcome[index, 1L]; m1 <- target$outcome[index, 2L]
    e <- target$propensity[index]
    score <- m1 - m0 + a * (y - m1) / e - (1 - a) * (y - m0) / (1 - e)
    c(estimate = mean(score), se = sqrt(var(score) / n))
  }
  pilot <- manifest[manifest$config == scenario, ]
  pilot_draws <- t(vapply(pilot$population_seed, draw, numeric(2L)))
  for (i in seq_len(nrow(pilot))) {
    saved <- read.csv(file.path(root, "tasks", pilot$task_id[i], "methods.csv"))
    oracle <- saved[saved$method == "Oracle-target", ]
    stopifnot(max(abs(pilot_draws[i, ] - c(oracle$estimate, oracle$se))) < 1e-12)
  }
  seeds <- 948000L + seq_len(2000L)
  sampled <- t(vapply(seeds, draw, numeric(2L)))
  error <- sampled[, "estimate"] - population$truth
  ratio <- var(sampled[, "estimate"]) / variance
  coverage <- mean(abs(error) <= qnorm(.975) * sampled[, "se"])
  stopifnot(abs(mean(error)) < 4 * sqrt(variance / length(seeds)),
    abs(ratio - 1) < .1, coverage > .925, coverage < .975)
  summaries[[scenario]] <- data.frame(scenario = scenario, truth = population$truth,
    pilot_repeats = nrow(pilot), pilot_empirical_sd = sd(pilot_draws[, "estimate"]),
    exact_population_sd = sqrt(variance), independent_oracle_repeats = length(seeds),
    independent_bias = mean(error), independent_empirical_sd = sd(sampled[, "estimate"]),
    independent_variance_ratio = ratio, independent_mean_se = mean(sampled[, "se"]),
    independent_coverage = coverage)
  records[[scenario]] <- data.frame(scenario = scenario, seed = seeds,
    estimate = sampled[, "estimate"], se = sampled[, "se"], truth = population$truth)
}
dir.create(output, recursive = TRUE)
write.csv(do.call(rbind, summaries), file.path(output, "summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, records), file.path(output, "oracle_replicates.csv"), row.names = FALSE)
writeLines("Oracle-only generator/variance check. The original fitted 20-repeat pilot is retained; these draws cannot be pooled with it to claim fitted-method performance.",
  file.path(output, "scope.txt"))
writeLines("Original oracle rows reproduced; independent oracle sampling agrees with exact population variance.",
  file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, summaries))
