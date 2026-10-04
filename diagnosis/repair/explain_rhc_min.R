#!/usr/bin/env Rscript
# Decompose the saved min fit; no nuisance or aggregation fit is changed.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
library_path <- normalizePath(arguments[1L], mustWork = TRUE)
primary <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !dir.exists(output))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
saved <- readRDS(file.path(primary, "checkpoints/analysis.rds"))$value
fit <- saved$tate_fit
raw <- load_rhc_raw()
cohort <- build_rhc_cohort(raw, outcome = "death30", site_var = "ninsclas",
  site_recode = c("No insurance" = NA_character_, "Medicare & Medicaid" = NA_character_))
data <- build_rhc_data_split(cohort, K = 4L, target_site = "Private", seed = 42L)
stopifnot(identical(digest::digest(data, algo = "sha256"), saved$metadata$data_sha256))
sizes <- vapply(data, `[[`, integer(1L), "n")
increments <- lapply(fit$arm_results[c("mu1", "mu0")], function(arm)
  RoCE:::.outer_fold_increments(arm$intermediates$fold_info, sizes[1L], sizes[-1L]))
contributions <- Map(function(d, eta) colSums(d * eta), increments, fit$fold_weights_by_arm)
decomposition <- data.frame(source = names(data)[-1L],
  mu1_increment = contributions$mu1, mu0_increment = contributions$mu0,
  tate_increment = contributions$mu1 - contributions$mu0)
stopifnot(abs(saved$target_anchor$estimate + sum(decomposition$tate_increment) - fit$estimate) < 1e-12)
training <- do.call(rbind, lapply(seq_len(fit$n_folds), function(fold) {
  moments <- fit$intermediates$joint_moments[[fold]]
  selected <- fit$intermediates$joint_fits[[fold]]
  eta <- selected$weights
  value <- moments$constant + 2 * sum(moments$l * eta) + drop(crossprod(eta, moments$Q %*% eta))
  stopifnot(selected$objective <= moments$N_all * moments$constant + 1e-8)
  data.frame(fold = fold, training_anchor_variance = moments$constant,
    training_aggregate_variance = value, training_variance_ratio = value/moments$constant,
    kkt_residual = selected$kkt_residual)
}))
site <- "s3"
source_index <- match(site, names(data)[-1L])
record <- 423L
arm <- fit$arm_results$mu1
fold <- which(vapply(arm$intermediates$fold_info, function(info) record %in% info$source_idx[[source_index]], logical(1L)))
stopifnot(length(fold) == 1L, data[[site]]$A[record] == 1L)
ids <- arm$intermediates$fold_info[[fold]]$source_idx[[source_index]]
model <- arm$fold_results[[fold]]$source_results[[site]]
spec <- RoCE:::resolve_glm_family("binomial")
prediction <- as.numeric(RoCE:::predict_glm_cpp(data[[site]]$W_outcome[record, , drop = FALSE],
  model$alpha_ts, spec$family_int, spec$link_int))
correction <- model$correction_components[match(record, ids)]
eta <- fit$fold_weights_by_arm$mu1[fold, site]
point_contribution <- eta * correction/sizes[site]
patient <- data.frame(site = site, record_index = record, fold = fold,
  observed_y = data[[site]]$Y[record], outcome_prediction = prediction,
  merged_weight = correction/(data[[site]]$Y[record] - prediction),
  aggregation_weight = eta, residual_contribution = point_contribution,
  score_only_neutralized_estimate = fit$estimate-point_contribution)
summary <- data.frame(target_anchor = saved$target_anchor$estimate, estimate = fit$estimate,
  reported_se = fit$se, fixed_weight_se = fit$se_fixed_weights,
  indirect_variance = fit$weight_layer$indirect_variance, cross_term = fit$weight_layer$cross_term)
dir.create(output, recursive = TRUE)
for (name in c("decomposition", "training", "patient", "summary")) {
  write.csv(get(name), file.path(output, paste0(name, ".csv")), row.names = FALSE)
}
writeLines(c("SAVED_MIN_POINT_AND_TRAINING_OBJECTIVE_RECONSTRUCTION_PASSED",
  "Neutralization holds every fitted nuisance and eta fixed; it is not a deletion/refit estimate or a bias estimate."),
  file.path(output, "CHECKS_PASSED"))
print(decomposition); print(training); print(patient); print(summary)
