test_that("diagnose completed n5000 C2 coverage summaries", {
  requested_results_dir <- Sys.getenv("ROCE_COVERAGE_DIAG_RESULTS_DIR", unset = NA_character_)
  candidate_dirs <- if (!is.na(requested_results_dir) && nzchar(requested_results_dir)) {
    requested_results_dir
  } else {
    c("results", file.path("..", "..", "results"))
  }
  existing_dirs <- candidate_dirs[dir.exists(candidate_dirs)]
  if (length(existing_dirs) == 0L) {
    skip(sprintf(
      "coverage diagnostics require one of these results directories: %s",
      paste(candidate_dirs, collapse = ", ")
    ))
  }
  results_dir <- existing_dirs[1L]

  summary_files <- list.files(
    results_dir,
    pattern = "^n5000_K[0-9]+_p[0-9]+_C[123]_.*_summary[.]csv$",
    full.names = TRUE
  )
  if (length(summary_files) == 0L) {
    skip(sprintf("no completed n5000 C1/C2/C3 summary CSVs found in %s", results_dir))
  }

  parse_setting <- function(path) {
    name <- basename(path)
    match <- regexec("^n([0-9]+)_K([0-9]+)_p([0-9]+)_(C[0-9])_", name)
    parts <- regmatches(name, match)[[1L]]
    if (length(parts) != 5L) {
      stop(sprintf("Could not parse setting from result filename: %s", name),
           call. = FALSE)
    }
    data.frame(
      result_file = name,
      n_total_file = as.integer(parts[2L]),
      K_file = as.integer(parts[3L]),
      p = as.integer(parts[4L]),
      config_file = parts[5L],
      stringsAsFactors = FALSE
    )
  }

  read_summary <- function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    meta <- parse_setting(path)
    df$result_file <- meta$result_file
    df$p <- meta$p
    df$config_file <- meta$config_file
    df
  }

  summaries <- do.call(rbind, lapply(summary_files, read_summary))
  required <- c("method", "config", "n_total", "K", "p", "bias_mean",
                "bias_sd", "rmse", "se_mean", "coverage", "ci_width_mean",
                "n_success")
  missing_required <- setdiff(required, names(summaries))
  expect_equal(missing_required, character(0))
  expect_true(all(is.finite(summaries$coverage)))
  expect_true(all(summaries$n_success > 0L))

  summaries$abs_bias <- abs(summaries$bias_mean)
  summaries$bias_over_se <- summaries$abs_bias / summaries$se_mean
  summaries$se_over_emp_sd <- summaries$se_mean / summaries$bias_sd
  summaries$half_width_over_abs_bias <- (summaries$ci_width_mean / 2) /
    pmax(summaries$abs_bias, .Machine$double.eps)

  c2 <- summaries[summaries$config == "C2", , drop = FALSE]
  expect_gt(nrow(c2), 0L)

  diagnostic_cols <- c(
    "K", "p", "method", "coverage", "bias_mean", "bias_sd", "se_mean",
    "bias_over_se", "se_over_emp_sd", "ci_width_mean", "rmse", "n_success"
  )

  c2_sorted <- c2[order(c2$coverage, c2$K, c2$p, c2$method),
                  diagnostic_cols, drop = FALSE]
  rownames(c2_sorted) <- NULL
  cat("\n[coverage-diagnostics] n=5000, C2 methods sorted by coverage:\n")
  print(c2_sorted, row.names = FALSE, digits = 4)

  focus_methods <- c("one_round_crossfit", "two_round_crossfit",
                     "federated_dr", "target_only", "pooled_dr",
                     "tilted_aipw", "oracle_dr")
  focus <- c2[c2$method %in% focus_methods,
              diagnostic_cols, drop = FALSE]
  focus <- focus[order(focus$K, focus$p, focus$method), , drop = FALSE]
  rownames(focus) <- NULL
  cat("\n[coverage-diagnostics] focused C2 rows:\n")
  print(focus, row.names = FALSE, digits = 4)

  low <- c2[c2$coverage < 0.93,
            c("K", "p", "method", "coverage", "bias_mean", "se_mean",
              "bias_over_se", "se_over_emp_sd", "ci_width_mean"),
            drop = FALSE]
  low <- low[order(low$coverage, low$K, low$p, low$method), , drop = FALSE]
  rownames(low) <- NULL
  cat("\n[coverage-diagnostics] C2 rows with coverage < 0.93:\n")
  print(low, row.names = FALSE, digits = 4)

  if (nrow(low) > 0L) {
    cat("\n[coverage-diagnostics] interpretation:\n")
    cat("  - bias_over_se near or above 0.5 means first-order bias is large relative to the reported SE.\n")
    cat("  - se_over_emp_sd below 1 means reported SE is smaller than empirical Monte Carlo SD.\n")
    cat("  - Rows with low coverage and se_over_emp_sd near 1 are primarily bias-driven, not SE-underestimation-driven.\n")
  }

  compare <- summaries[summaries$method %in% focus_methods,
                       c("config", "K", "p", "method", "coverage",
                         "bias_mean", "se_mean", "bias_over_se",
                         "se_over_emp_sd"),
                       drop = FALSE]
  compare <- compare[order(compare$K, compare$p, compare$method,
                           compare$config), , drop = FALSE]
  rownames(compare) <- NULL
  cat("\n[coverage-diagnostics] C1/C2/C3 comparison for focused methods:\n")
  print(compare, row.names = FALSE, digits = 4)
})
