#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1L) args[[1L]] else {
  "results/direct_tate_mc500_b5000/rhc_supported_primary"
}
result_path <- file.path(output_root, "rhc_direct_tate.rds")
if (!file.exists(result_path)) stop("RHC result not found: ", result_path)

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))

result <- readRDS(result_path)
fit <- result$tate_fit
excluded_levels <- result$metadata$excluded_site_levels
excluded_levels <- if (
    is.character(excluded_levels) && length(excluded_levels) == 1L &&
    nzchar(excluded_levels) && !identical(excluded_levels, "none")) {
  strsplit(excluded_levels, ";", fixed = TRUE)[[1L]]
} else {
  character(0)
}
site_recode <- if (length(excluded_levels) > 0L) {
  stats::setNames(rep(NA_character_, length(excluded_levels)), excluded_levels)
} else {
  NULL
}
cohort <- build_rhc_cohort(
  outcome = result$metadata$outcome,
  site_var = result$metadata$site_var,
  site_recode = site_recode
)
data_split <- build_rhc_data_split(
  cohort = cohort,
  K = as.integer(result$metadata$n_sites),
  target_site = unname(result$metadata$target_site),
  seed = 42L
)
expected_sizes <- as.integer(result$metadata$site_n)
actual_sizes <- vapply(data_split, `[[`, integer(1L), "n")
if (!identical(unname(expected_sizes), unname(actual_sizes))) {
  stop("Reconstructed RHC site sizes do not match the fitted result.")
}

n_folds <- as.integer(result$metadata$n_folds)
folds <- RoCE:::build_crossfit_folds(data_split, n_folds)
summary_rows <- list()
flag_rows <- list()
tolerance <- 1e-10

for (arm in intersect(c("mu1", "mu0"), names(fit$arm_results))) {
  treatment <- if (identical(arm, "mu1")) 1L else 0L
  arm_result <- fit$arm_results[[arm]]
  for (outer_fold in seq_len(n_folds)) {
    source_results <- arm_result$fold_results[[outer_fold]]$source_results
    for (source in names(source_results)) {
      gamma_list <- source_results[[source]]$per_k2_gamma
      for (inner_name in names(gamma_list)) {
        inner_fold <- as.integer(sub("^k2_", "", inner_name))
        training_folds <- setdiff(
          seq_len(n_folds), c(outer_fold, inner_fold)
        )
        target_train <- RoCE:::combine_folds(
          folds$target_folds, training_folds
        )
        source_train <- RoCE:::combine_folds(
          folds$source_folds[[source]], training_folds
        )
        arm_rows <- source_train$A == treatment
        design <- cbind(1, source_train$Z_site[arm_rows, , drop = FALSE])
        target_mean <- c(1, colMeans(target_train$Z_site))
        feature_names <- c("(Intercept)", colnames(source_train$Z_site))
        if (is.null(feature_names) || length(feature_names) != ncol(design)) {
          feature_names <- c(
            "(Intercept)", paste0("Z", seq_len(ncol(design) - 1L))
          )
        }
        gamma <- gamma_list[[inner_name]]
        lambda <- as.numeric(attr(gamma, "lambda_used"))
        source_min <- apply(design, 2L, min)
        source_max <- apply(design, 2L, max)
        source_constant <- (source_min + source_max) / 2
        zero_variance <- (source_max - source_min) <= tolerance
        # If a source-arm feature is constant, an intercept/feature shift can
        # keep every source linear predictor unchanged.  The penalized initial
        # DR objective is unbounded along that direction whenever the target-
        # source moment gap exceeds the L1 slope penalty.
        moment_gap <- target_mean - source_constant
        unsupported <- zero_variance &
          seq_along(zero_variance) > 1L &
          abs(moment_gap) > lambda + tolerance
        boundary <- abs(as.numeric(gamma)) >= 0.9 * 100
        max_index <- which.max(abs(as.numeric(gamma)))

        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          arm = arm,
          treatment = treatment,
          outer_fold = outer_fold,
          source = source,
          inner_fold = inner_name,
          n_source_train = source_train$n,
          n_source_arm = sum(arm_rows),
          lambda = lambda,
          converged = isTRUE(attr(gamma, "converged")),
          iterations = as.integer(attr(gamma, "iterations")),
          line_search_failures = as.integer(attr(gamma, "line_search_failures")),
          max_abs_coefficient = abs(as.numeric(gamma[[max_index]])),
          max_coefficient_name = feature_names[[max_index]],
          n_zero_variance_slopes = sum(zero_variance[-1L]),
          n_unbounded_zero_variance_directions = sum(unsupported),
          max_gap_minus_lambda = if (any(unsupported)) {
            max(abs(moment_gap[unsupported]) - lambda)
          } else {
            0
          },
          n_boundary_coefficients = sum(boundary),
          n_boundary_unsupported = sum(boundary & unsupported),
          stringsAsFactors = FALSE
        )

        flagged <- which(unsupported | boundary)
        if (length(flagged) > 0L) {
          flag_rows[[length(flag_rows) + 1L]] <- data.frame(
            arm = arm,
            treatment = treatment,
            outer_fold = outer_fold,
            source = source,
            inner_fold = inner_name,
            coefficient_index = flagged,
            coefficient_name = feature_names[flagged],
            gamma = as.numeric(gamma)[flagged],
            lambda = lambda,
            target_mean = target_mean[flagged],
            source_constant = source_constant[flagged],
            source_range = (source_max - source_min)[flagged],
            moment_gap = moment_gap[flagged],
            zero_variance = zero_variance[flagged],
            unbounded_direction = unsupported[flagged],
            boundary_coefficient = boundary[flagged],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

summary_result <- do.call(rbind, summary_rows)
flag_result <- if (length(flag_rows) > 0L) do.call(rbind, flag_rows) else {
  data.frame()
}
write.csv(
  summary_result,
  file.path(output_root, "rhc_direct_tate_initial_dr_support_summary.csv"),
  row.names = FALSE
)
write.csv(
  flag_result,
  file.path(output_root, "rhc_direct_tate_initial_dr_support_flags.csv"),
  row.names = FALSE
)

cat("Initial DR fits audited:", nrow(summary_result), "\n")
cat("Nonconverged fits:", sum(!summary_result$converged), "\n")
cat(
  "Fits with an empirically unbounded zero-variance direction:",
  sum(summary_result$n_unbounded_zero_variance_directions > 0L), "\n"
)
cat(
  "Fits with a coefficient at >=90% of the hard numerical bound:",
  sum(summary_result$n_boundary_coefficients > 0L), "\n"
)
