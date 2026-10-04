#!/usr/bin/env Rscript
# Diagnose saved RHC fits without refitting or changing their specification.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
library_path <- normalizePath(arguments[1L], mustWork = TRUE)
study <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
configuration <- jsonlite::fromJSON(file.path(study, "configuration.json"))
interface <- new.env(parent = asNamespace("RoCE"))
sys.source(file.path(study, "workflow/real_data_rhc.R"), envir = interface)
raw <- load_rhc_raw()
excluded <- c("Medicare & Medicaid", "No insurance")
recode <- stats::setNames(rep(NA_character_, 2L), excluded)
cohort_arguments <- list(raw = raw, outcome = "death30", site_var = "ninsclas", site_recode = recode)
if (!is.null(configuration$covariate_profile)) cohort_arguments$covariate_profile <- configuration$covariate_profile
cohort <- do.call(interface$build_rhc_cohort, cohort_arguments)
data <- interface$build_rhc_data_split(cohort, K = 4L, target_site = configuration$target_site, seed = 42L)
raw_kept <- raw[!raw$ninsclas %in% excluded, , drop = FALSE]
stopifnot(nrow(raw_kept) == nrow(cohort))
mapping <- attr(data, "site_mapping")
sizes <- vapply(data, `[[`, integer(1L), "n")
site_rows <- split(seq_len(sum(sizes)), rep(names(data), sizes))
raw_by_site <- lapply(names(data), function(site) raw_kept[raw_kept$ninsclas ==
  mapping[if (site == "t") "target" else site], , drop = FALSE])
names(raw_by_site) <- names(data)
manifest <- read.csv(file.path(study, "manifest.csv"), stringsAsFactors = FALSE)
site_summary <- weight_summary <- top_records <- covariate_balance <- list()
spec <- RoCE:::resolve_glm_family("binomial")
dir.create(output, recursive = TRUE)
site_labels <- unname(mapping[c("target", names(data)[-1L])])
write.csv(data.frame(site = names(data), insurance = site_labels, n = sizes,
  treated = vapply(data, function(sample) sum(sample$A == 1L), numeric(1L))),
  file.path(output, "site_mapping.csv"), row.names = FALSE)
support <- list()
target_mean <- colMeans(data$t$X)
target_sd <- apply(data$t$X, 2L, sd)
variable_target_features <- target_sd > 0
for (site in names(data)[-1L]) for (arm in 0:1) {
  sample <- data[[site]]
  observed <- sample$X[sample$A == arm, , drop = FALSE]
  minimum <- apply(observed, 2L, min)
  maximum <- apply(observed, 2L, max)
  support[[length(support)+1L]] <- data.frame(site = site, arm = arm, feature = names(target_mean),
    target_mean = target_mean, source_arm_min = minimum, source_arm_max = maximum,
    target_mean_outside_range = target_mean < minimum - 1e-10 | target_mean > maximum + 1e-10)
}
write.csv(do.call(rbind, support), file.path(output, "source_arm_support.csv"), row.names = FALSE)

selected_tasks <- manifest$task_id[manifest$role == "method" &
  file.exists(file.path(study, "tasks", manifest$task_id, "COMPLETE"))]
stopifnot(length(selected_tasks) > 0L)
for (task_id in selected_tasks) {
  directory <- file.path(study, "tasks", task_id)
  stopifnot(file.exists(file.path(directory, "COMPLETE")),
    identical(readLines(file.path(directory, "data_sha256.txt")), digest::digest(data, algo = "sha256")))
  profile <- manifest$config[manifest$task_id == task_id]
  fit <- readRDS(file.path(directory, "checkpoints/analysis.rds"))$value$tate_fit
  influence <- readRDS(file.path(directory, "arm_influence.rds"))
  tau_if <- influence$influence[, "mu1"] - influence$influence[, "mu0"]
  fixed_if <- influence$influence_fixed[, "mu1"] - influence$influence_fixed[, "mu0"]
  stopifnot(abs(sum(tau_if^2) - fit$variance) < 1e-12,
            abs(sum(fixed_if^2) - fit$variance_fixed_weights) < 1e-12)
  weights <- lapply(names(data)[-1L], function(site) matrix(NA_real_, data[[site]]$n, 2L,
    dimnames = list(NULL, c("mu1", "mu0"))))
  names(weights) <- names(data)[-1L]
  predictions <- weights
  raw_logits <- weights
  for (arm in c("mu1", "mu0")) {
    a <- if (arm == "mu1") 1L else 0L
    fitted_arm <- fit$arm_results[[arm]]
    for (fold in seq_len(fit$n_folds)) for (site in names(weights)) {
      ids <- fitted_arm$intermediates$fold_info[[fold]]$source_idx[[site]]
      if (is.null(ids)) ids <- fitted_arm$intermediates$fold_info[[fold]]$source_idx[[match(site, names(weights))]]
      selected <- fitted_arm$fold_results[[fold]]$source_results[[site]]
      sample <- data[[site]]
      W <- sample$W_outcome[ids, , drop = FALSE]
      Z <- sample$Z_site[ids, , drop = FALSE]
      m <- RoCE:::predict_glm_cpp(W, selected$alpha_ts, spec$family_int, spec$link_int)
      # A unit residual isolates the exact deployed weight, including all caps.
      unit_residual <- RoCE:::calculate_correction_term_cpp(Z, sample$A[ids], as.numeric(m)+1,
        selected$gamma_s, selected$alpha_ts, W, fit$M_tau_inference,
        spec$family_int, spec$link_int, a)$correction_components
      reconstructed <- as.numeric(unit_residual) * (sample$Y[ids] - as.numeric(m))
      stopifnot(max(abs(reconstructed - selected$correction_components)) < 1e-10)
      weights[[site]][ids, arm] <- as.numeric(unit_residual)
      predictions[[site]][ids, arm] <- as.numeric(m)
      raw_logits[[site]][ids, arm] <- drop(cbind(1, Z) %*% selected$gamma_s)
    }
    for (site in names(weights)) {
      sample <- data[[site]]
      active <- sample$A == a
      w <- weights[[site]][active, arm]
      stopifnot(all(is.finite(w)), all(w > 0), !anyNA(predictions[[site]][, arm]))
      quantiles <- quantile(w, c(.5, .9, .95, .99, 1), names = FALSE)
      weight_summary[[length(weight_summary)+1L]] <- data.frame(profile = profile, site = site,
        arm = arm, n = length(w), effective_n = sum(w)^2/sum(w^2), mean_weight = mean(w),
        median = quantiles[1L], p90 = quantiles[2L], p95 = quantiles[3L], p99 = quantiles[4L],
        max_weight = quantiles[5L], largest_weight_share = max(w)/sum(w),
        clipped_count = sum(abs(raw_logits[[site]][active, arm]) > fit$M_tau_inference),
        at_upper_cap_count = sum(w >= exp(fit$M_tau_inference) * (1 - 1e-12)),
        at_lower_cap_count = sum(w <= exp(-fit$M_tau_inference) * (1 + 1e-12)))
      weighted_mean <- colSums(sample$X[active, , drop = FALSE] * w)/sum(w)
      mean_gap <- weighted_mean - target_mean
      gap <- setNames(rep(NA_real_, length(mean_gap)), names(mean_gap))
      gap[variable_target_features] <- mean_gap[variable_target_features]/target_sd[variable_target_features]
      covariate_balance[[length(covariate_balance)+1L]] <- data.frame(profile = profile, site = site, arm = arm,
        feature = names(gap), weighted_minus_target_smd = as.numeric(gap),
        weighted_minus_target_mean = as.numeric(mean_gap), target_sd = as.numeric(target_sd))
    }
  }
  for (site in names(data)) {
    ids <- site_rows[[site]]
    values <- tau_if[ids]
    variance <- sum(values^2)
    largest <- order(values^2, decreasing = TRUE)
    site_summary[[length(site_summary)+1L]] <- data.frame(profile = profile, site = site,
      n = length(ids), variance = variance, variance_share = variance/fit$variance,
      fixed_variance = sum(fixed_if[ids]^2), top1_total_variance_share = max(values^2)/fit$variance,
      top5_total_variance_share = sum(sort(values^2, decreasing = TRUE)[seq_len(5L)])/fit$variance)
    for (row in head(largest, 5L)) {
      record <- raw_by_site[[site]][row, , drop = FALSE]
      top_records[[length(top_records)+1L]] <- data.frame(profile = profile, site = site,
        record_index = row, A = data[[site]]$A[row], Y = data[[site]]$Y[row],
        influence = values[row], total_variance_share = values[row]^2/fit$variance,
        weight = if (site == "t") NA_real_ else weights[[site]][row, if (data[[site]]$A[row] == 1) "mu1" else "mu0"],
        age = record$age, primary_disease = record$cat1, urine = record$urin1,
        creatinine = record$crea1, bilirubin = record$bili1, white_cells = record$wblc1,
        mean_bp = record$meanbp1, aps = record$aps1)
    }
  }
}
write.csv(do.call(rbind, site_summary), file.path(output, "site_variance.csv"), row.names = FALSE)
write.csv(do.call(rbind, weight_summary), file.path(output, "source_weight_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, top_records), file.path(output, "influential_records.csv"), row.names = FALSE)
write.csv(do.call(rbind, covariate_balance), file.path(output, "weighted_balance.csv"), row.names = FALSE)
writeLines("SAVED_RHC_WEIGHT_AND_VARIANCE_RECONSTRUCTION_PASSED", file.path(output, "CHECKS_PASSED"))
first_profile <- manifest$config[match(selected_tasks[1L], manifest$task_id)]
print(subset(do.call(rbind, site_summary), profile == first_profile))
print(subset(do.call(rbind, weight_summary), profile == first_profile))
