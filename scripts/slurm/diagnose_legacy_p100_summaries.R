#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_directory <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "diagnosis/face_probe/validation"
}
output_directory <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  input_directory
}
expected_replications <- if (length(args) >= 3L) {
  as.integer(args[[3L]])
} else {
  200L
}
if (is.na(expected_replications) || expected_replications < 1L) {
  stop("expected_replications must be one positive integer.")
}

files <- list.files(
  input_directory,
  pattern = "^negtransfer_C[123]_p100_K(2|4|8)\\.csv$",
  full.names = TRUE
)
if (length(files) != 9L) {
  stop("expected exactly nine legacy p=100 summary files; found ",
       length(files))
}

summary <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
summary$estimand_scope <- ifelse(
  grepl("_ate$", summary$method), "TATE", "treated_mean"
)
summary$expected_replications <- expected_replications
summary$bias_mcse <- summary$bias_sd / sqrt(summary$n_success)
summary$coverage_mcse <- sqrt(
  summary$coverage * (1 - summary$coverage) / summary$n_success
)
nominal_mcse <- sqrt(0.95 * 0.05 / summary$n_success)
summary$coverage_mc_lower <- pmax(
  0, 0.95 - stats::qnorm(0.975) * nominal_mcse
)
summary$coverage_mc_upper <- pmin(
  1, 0.95 + stats::qnorm(0.975) * nominal_mcse
)
summary$se_to_empirical_sd <- summary$se_mean / summary$bias_sd
summary$normal_reference_coverage <- stats::pnorm(
  (stats::qnorm(0.975) * summary$se_mean - summary$bias_mean) /
    summary$bias_sd
) - stats::pnorm(
  (-stats::qnorm(0.975) * summary$se_mean - summary$bias_mean) /
    summary$bias_sd
)
summary$rmse_from_moments <- sqrt(
  summary$bias_mean^2 +
    ((summary$n_success - 1) / summary$n_success) * summary$bias_sd^2
)
summary$rmse_identity_error <- abs(
  summary$rmse - summary$rmse_from_moments
)

setting_key <- interaction(
  summary[c("config", "p", "K", "rho", "estimand_scope")],
  drop = TRUE, lex.order = TRUE
)
summary$rmse_rank <- ave(
  summary$rmse, setting_key,
  FUN = function(x) rank(x, ties.method = "min", na.last = "keep")
)
target_rmse <- ave(
  seq_len(nrow(summary)), setting_key,
  FUN = function(index) {
    setting <- summary[index, , drop = FALSE]
    target_name <- if (setting$estimand_scope[[1L]] == "TATE") {
      "target_only_ate"
    } else {
      "target_only"
    }
    value <- setting$rmse[setting$method == target_name]
    rep(if (length(value) == 1L) value else NA_real_, length(index))
  }
)
summary$rmse_relative_to_target <- summary$rmse / target_rmse

summary$diagnostic_status <- vapply(seq_len(nrow(summary)), function(index) {
  flags <- c(
    incomplete_replications =
      summary$n_success[[index]] < expected_replications,
    excess_replications =
      summary$n_success[[index]] > expected_replications,
    rmse_identity_failed =
      summary$rmse_identity_error[[index]] > 1e-10,
    coverage_below_mc_band =
      summary$coverage[[index]] < summary$coverage_mc_lower[[index]],
    coverage_above_mc_band =
      summary$coverage[[index]] > summary$coverage_mc_upper[[index]],
    se_empirical_ratio_low =
      summary$se_to_empirical_sd[[index]] < 0.85,
    se_empirical_ratio_high =
      summary$se_to_empirical_sd[[index]] > 1.15,
    bias_exceeds_2_mcse =
      abs(summary$bias_mean[[index]]) > 2 * summary$bias_mcse[[index]]
  )
  active <- names(flags)[flags]
  if (length(active) == 0L) "ok" else paste(active, collapse = ";")
}, character(1L))

summary <- summary[order(
  summary$config, summary$K, summary$rho, summary$estimand_scope,
  summary$rmse_rank, summary$method
), ]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  summary,
  file.path(output_directory, "legacy_p100_setting_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  summary[summary$diagnostic_status != "ok", , drop = FALSE],
  file.path(output_directory, "legacy_p100_flagged_settings.csv"),
  row.names = FALSE
)
write.csv(
  summary[summary$rho == 0, , drop = FALSE],
  file.path(output_directory, "legacy_p100_rho0_rmse_ranking.csv"),
  row.names = FALSE
)
message("diagnosed ", nrow(summary), " legacy p=100 setting-method rows")
