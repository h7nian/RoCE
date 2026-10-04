#!/usr/bin/env Rscript
# Known-Gaussian references only; this is not a fitted RoCE coverage study.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: source_directory scratch_output repeats")
source_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
source_files <- file.path(source_directory, c("source_confidence_sets.R", "arm_pair_confidence_sets.R"))
for (path in source_files) source(path)
output <- arguments[[2L]]
repeats <- as.integer(arguments[[3L]])
check_integer(repeats, "repeats", 100, 10000)
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "cells"), showWarnings = FALSE)
settings <- expand.grid(num_sources = c(4L, 16L),
  profile = c("treated_half", "disjoint_half", "disjoint_quarter"),
  local_bias = c(0, .5, 2, 4), stringsAsFactors = FALSE)
configuration <- list(settings = settings, repeats = repeats, n_per_site = 1000L,
  mu1 = .6, mu0 = .4, alpha = .05,
  code_md5 = tools::md5sum(c(source_files,
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))))
configuration_file <- file.path(output, "configuration.rds")
if (file.exists(configuration_file)) stopifnot(identical(readRDS(configuration_file), configuration)) else
  saveRDS(configuration, configuration_file)

linear_gls <- function(covariance, arm_design, included) {
  precision_design <- solve(covariance[included, included, drop = FALSE], arm_design[included, , drop = FALSE])
  information <- crossprod(arm_design[included, , drop = FALSE], precision_design)
  contrast <- c(1, -1)
  weights <- numeric(nrow(covariance))
  weights[included] <- drop(precision_design %*% solve(information, contrast))
  list(weights = weights, se = sqrt(drop(crossprod(contrast, solve(information, contrast)))))
}

all_metrics <- list()
truth <- configuration$mu1 - configuration$mu0
for (cell in seq_len(nrow(settings))) {
  path <- file.path(output, "cells", paste0(cell, ".rds"))
  if (file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(identical(saved$configuration, configuration), saved$cell == cell)
    all_metrics[[cell]] <- saved$metrics
    next
  }
  setting <- settings[cell, ]
  count <- setting$num_sources
  width <- count + 1L
  candidates <- c("target_anchor", paste0("s", seq_len(count)))
  axis <- seq(-1, 1, length.out = width)
  # Unequal target loadings and correlated local arm noise give a full joint
  # covariance. A common scalar target factor is not imposed.
  loading <- rbind(cbind(.7, .3 * cos(axis), .25 * axis),
                   cbind(.45, .7 + .2 * sin(axis), -.15 * axis))
  covariance <- tcrossprod(loading)
  for (candidate in seq_len(width)) {
    positions <- c(candidate, width + candidate)
    variance_mu1 <- .7 + candidate / width
    variance_mu0 <- 1.4 - .5 * candidate / width
    correlation <- if (candidate %% 2L) -.3 else .4
    cross <- correlation * sqrt(variance_mu1 * variance_mu0)
    covariance[positions, positions] <- covariance[positions, positions] +
      matrix(c(variance_mu1, cross, cross, variance_mu0), 2L)
  }
  covariance <- covariance / configuration$n_per_site
  labels <- c(paste0("mu1:", candidates), paste0("mu0:", candidates))
  dimnames(covariance) <- list(labels, labels)
  q1 <- if (setting$profile == "disjoint_quarter") count / 4 else count / 2
  q0 <- if (setting$profile == "treated_half") count else count - q1
  valid_mu1 <- seq_len(q1)
  valid_mu0 <- if (setting$profile == "treated_half") seq_len(count) else seq.int(q1 + 1L, count)
  biases <- numeric(2L * width)
  biases[1L + setdiff(seq_len(count), valid_mu1)] <- setting$local_bias / sqrt(configuration$n_per_site)
  biases[width + 1L + setdiff(seq_len(count), valid_mu0)] <- -setting$local_bias / sqrt(configuration$n_per_site)
  population <- c(rep(configuration$mu1, width), rep(configuration$mu0, width)) + biases
  design <- cbind(c(rep(1, width), rep(0, width)), c(rep(0, width), rep(1, width)))
  valid_positions <- if (setting$local_bias == 0) seq_len(2L * width) else
    c(1L, 1L + valid_mu1, width + 1L, width + 1L + valid_mu0)
  oracle <- linear_gls(covariance, design, valid_positions)
  pooled <- linear_gls(covariance, design, seq_len(2L * width))
  differences <- cbind(diag(width), -diag(width))
  tate_covariance <- differences %*% covariance %*% t(differences)
  tate_se <- sqrt(diag(tate_covariance))
  common_valid <- max(q1 + q0 - count, 0)
  scalar_calibration <- make_source_calibration(count, common_valid,
    if (common_valid) max(1L, floor(common_valid / 2)) else 0L)
  set.seed(68000L + match(count, c(4L, 16L)) * 100L +
             match(setting$profile, c("treated_half", "disjoint_half", "disjoint_quarter")))
  draws <- sweep(matrix(rnorm(repeats * 2L * width), repeats, 2L * width) %*% chol(covariance),
                 2L, population, "+")
  records <- vector("list", repeats)
  critical <- qnorm(.975)
  for (iteration in seq_len(repeats)) {
    observed <- draws[iteration, ]
    means <- cbind(mu1 = observed[seq_len(width)], mu0 = observed[width + seq_len(width)])
    rownames(means) <- candidates
    pair_half <- arm_pair_search_interval(means, covariance, q1, q0)
    pair_all <- arm_pair_search_interval(means, covariance, q1, q0,
                                         votes_required = (q1 + 1) * (q0 + 1))
    tate <- drop(differences %*% observed)
    scalar <- source_search_interval(tate[1L], tate_se[1L], tate[-1L], tate_se[-1L], scalar_calibration)
    intervals <- list(pair_half = c(pair_half$lower, pair_half$upper),
      pair_all = c(pair_all$lower, pair_all$upper), scalar_tate = c(scalar$lower, scalar$upper),
      target = tate[1L] + c(-1, 1) * critical * tate_se[1L],
      oracle_arm_gls = sum(oracle$weights * observed) + c(-1, 1) * critical * oracle$se,
      pooled_arm_gls = sum(pooled$weights * observed) + c(-1, 1) * critical * pooled$se)
    records[[iteration]] <- data.frame(iteration = iteration, method = names(intervals),
      lower = vapply(intervals, `[`, numeric(1L), 1L), upper = vapply(intervals, `[`, numeric(1L), 2L),
      anchor_fallback = c(pair_half$used_anchor_fallback, pair_all$used_anchor_fallback,
                          scalar$used_anchor_fallback, FALSE, FALSE, FALSE))
  }
  records <- do.call(rbind, records)
  metrics <- do.call(rbind, lapply(split(records, records$method), function(rows) {
    covered <- rows$lower <= truth & rows$upper >= truth
    ci <- binom.test(sum(covered), nrow(rows))$conf.int
    data.frame(cell = cell, setting, method = rows$method[1L], repeats = nrow(rows),
      assumed_valid_mu1 = q1, assumed_valid_mu0 = q0, coverage = mean(covered),
      coverage_lower = ci[1L], coverage_upper = ci[2L], mean_length = mean(rows$upper - rows$lower),
      fallback_rate = mean(rows$anchor_fallback), row.names = NULL)
  }))
  saved <- list(configuration = configuration, cell = cell, metrics = metrics, records = records,
                covariance = covariance, population_means = population)
  temporary <- paste0(path, ".tmp")
  saveRDS(saved, temporary)
  stopifnot(file.rename(temporary, path))
  all_metrics[[cell]] <- metrics
  write.csv(do.call(rbind, all_metrics), file.path(output, "metrics.csv"), row.names = FALSE)
  cat("Completed arm-pair Gaussian cell", cell, "of", nrow(settings), "\n")
}
write.csv(do.call(rbind, all_metrics), file.path(output, "metrics.csv"), row.names = FALSE)
writeLines("All prespecified known-Gaussian arm-pair cells completed; fitted RoCE coverage remains unproved.",
           file.path(output, "COMPLETE"))
