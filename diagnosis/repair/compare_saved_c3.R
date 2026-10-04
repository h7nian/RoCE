#!/usr/bin/env Rscript
# Quick, read-only comparison of the saved C3/K2/rho0/seed1 results. This compares
# statistical versions; it is not an acceleration-equivalence or coverage test.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: compare_saved_c3.R INSTALLED_LIBRARY NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
base <- "/scratch.global/zhan9381/FACE-HD"
paths <- c(
  production = "/users/0/zhan9381/FACE-HD/results/direct_tate_mc500_b5000/production_20260917_v4/raw/task_006001.csv",
  production_fit = file.path(base, "diagnosis/fixed_n_audit_20260919/full_native/C3/seed_0001/fits.rds"),
  legacy = file.path(base, "implementation/r2/pilot_v2/tasks/3/legacy.rds"),
  score = file.path(base, "implementation/r4/pilot_c3_shared_v1/tasks/1/score_derivative.rds"))
write.csv(data.frame(input = names(paths), path = paths,
  sha256 = vapply(paths, digest::digest, character(1L), algo = "sha256", file = TRUE)),
  file.path(output, "inputs.csv"), row.names = FALSE)
production <- read.csv(paths[["production"]], stringsAsFactors = FALSE)
saved_production <- readRDS(paths[["production_fit"]])
legacy <- readRDS(paths[["legacy"]])
score <- readRDS(paths[["score"]])
stopifnot(!inherits(legacy, "error"), !inherits(score, "error"))
artifacts <- lapply(list(legacy = legacy, score = score), attr, which = "roce_simulation_artifacts")
data <- list(production = saved_production$data_split,
             legacy = artifacts$legacy$data_split, score = artifacts$score$data_split)
fits <- list(production = saved_production$reference,
             legacy = artifacts$legacy$direct_tate_results$one_round_crossfit,
             score = artifacts$score$direct_tate_results$one_round_crossfit)
checks <- list()
for (name in names(data)) for (site in names(data$score)) {
  checks[[paste(name, site)]] <- data.frame(version = name, site = site,
    all_fields_identical = identical(data[[name]][[site]], data$score[[site]]),
    analysis_data_identical = all(vapply(c("W_outcome", "Z_site", "A", "Y", "n"), function(field)
      !is.null(data[[name]][[site]][[field]]) &&
        identical(data[[name]][[site]][[field]], data$score[[site]][[field]]), logical(1L))))
}
checks <- do.call(rbind, checks)
write.csv(checks, file.path(output, "data_identity.csv"), row.names = FALSE)
stopifnot(all(checks$analysis_data_identical))
fold_checks <- lapply(names(fits), function(name) {
  matches <- vapply(c("mu1", "mu0"), function(arm) {
    a <- fits[[name]]$arm_results[[arm]]$intermediates$fold_info
    b <- fits$score$arm_results[[arm]]$intermediates$fold_info
    stopifnot(length(a) == 10L, length(b) == 10L)
    all(vapply(seq_along(a), function(k) {
      stopifnot(length(a[[k]]$target_idx) > 0L, length(a[[k]]$source_idx) > 0L)
      identical(a[[k]]$target_idx, b[[k]]$target_idx) &&
        identical(a[[k]]$source_idx, b[[k]]$source_idx)
    }, logical(1L)))
  }, logical(1L))
  data.frame(version = name, arm = names(matches), indices_identical = unname(matches))
})
fold_checks <- do.call(rbind, fold_checks)
write.csv(fold_checks, file.path(output, "fold_identity.csv"), row.names = FALSE)
stopifnot(all(fold_checks$indices_identical))

rows <- lapply(list(production = production, legacy = legacy, score = score), function(x) {
  x[x$estimand_scope == "tate", c("method", "estimate", "se", "truth")]
})
stopifnot(all(production$config == "C3"), all(production$sim_id == 1L),
          all(production$rho == 0), all(production$K == 2L), all(production$p == 100L),
          all(production$n_site == 1000L), all(production$n_folds == 10L),
          all(production$nlambda_init == 100L), all(production$cutoff == 1))
for (name in names(rows)) {
  direct <- rows[[name]][rows[[name]]$method == "one_round_crossfit_ate", ]
  stopifnot(nrow(direct) == 1L,
            abs(direct$estimate - fits[[name]]$estimate) < 1e-12,
            abs(direct$se - fits[[name]]$se) < 1e-12)
}
all_rows <- do.call(rbind, lapply(names(rows), function(name) {
  x <- rows[[name]]
  x$version <- name
  x$error <- x$estimate - x$truth
  x$ci_lower <- x$estimate - qnorm(.975) * x$se
  x$ci_upper <- x$estimate + qnorm(.975) * x$se
  x$contains_truth <- x$ci_lower <= x$truth & x$ci_upper >= x$truth
  x
}))
write.csv(all_rows, file.path(output, "estimates.csv"), row.names = FALSE)
paired <- do.call(rbind, lapply(c("production", "legacy"), function(name) {
  methods <- intersect(rows[[name]]$method, rows$score$method)
  do.call(rbind, lapply(methods, function(method) {
    a <- rows$score[rows$score$method == method, ]
    b <- rows[[name]][rows[[name]]$method == method, ]
    stopifnot(nrow(a) == 1L, nrow(b) == 1L, abs(a$truth - b$truth) < 1e-14)
    data.frame(reference = name, method = method, estimate_change = a$estimate - b$estimate,
      se_change = a$se - b$se, se_percent_change = 100 * (a$se / b$se - 1))
  }))
}))
write.csv(paired, file.path(output, "paired_changes.csv"), row.names = FALSE)
benchmark <- paired[paired$method == "target_only_ate", ]
stopifnot(all(abs(benchmark$estimate_change) < 1e-14), all(abs(benchmark$se_change) < 1e-14))
weights <- do.call(rbind, lapply(names(fits), function(name) {
  fit <- fits[[name]]
  w <- as.data.frame(fit$fold_weights)
  w$fold <- seq_len(nrow(w))
  w$version <- name
  w
}))
write.csv(weights, file.path(output, "common_weights.csv"), row.names = FALSE)
selected <- all_rows[all_rows$method %in% c("one_round_crossfit_ate",
  "one_round_crossfit_ate_separate_arms", "one_round_crossfit_ate_joint_tate", "target_only_ate"), ]
print(selected, row.names = FALSE)
print(paired, row.names = FALSE)
writeLines(c("# Saved C3 comparison", "",
  "C3/K2/rho0/seed1; 1000/site, p100, ten folds, 100 lambdas and cutoff1.",
  "Analysis data and evaluation-fold indices are identical across all three versions.",
  "The archived original fit reproduces production raw estimates and SEs within 1e-12.",
  "The ordinary target-only benchmark is unchanged within 1e-14.", "",
  "production: frozen v4 two-level method; legacy: three-level R2 legacy calibration;",
  "score: three-level R4 score-derivative calibration with calibrated target PS and target radius log(9).", "",
  "These versions differ statistically. This comparison cannot establish acceleration equivalence.",
  "There is one independent data replicate. Errors and interval inclusion are descriptive;",
  "they are not Monte Carlo bias, RMSE or coverage-rate estimates.", "",
  "See estimates.csv, paired_changes.csv and common_weights.csv for numerical results."),
  file.path(output, "report.md"))
writeLines("PASSED: saved-result identity/arithmetic checks; no performance claim", file.path(output, "status.txt"))
