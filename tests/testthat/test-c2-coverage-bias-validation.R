find_coverage_results_dir <- function() {
  requested <- Sys.getenv("FACEHD_COVERAGE_DIAG_RESULTS_DIR", unset = NA_character_)
  candidates <- if (!is.na(requested) && nzchar(requested)) {
    requested
  } else {
    c("results", file.path("..", "..", "results"))
  }
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) == 0L) {
    skip(sprintf("coverage diagnostics require one of these results directories: %s",
                 paste(candidates, collapse = ", ")))
  }
  existing[1L]
}

read_n5000_coverage_summaries <- function() {
  results_dir <- find_coverage_results_dir()
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
    list(result_file = name, p = as.integer(parts[4L]), config_file = parts[5L])
  }

  rows <- lapply(summary_files, function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    meta <- parse_setting(path)
    df$result_file <- meta$result_file
    df$p <- meta$p
    df$config_file <- meta$config_file
    df
  })
  summaries <- do.call(rbind, rows)

  required <- c("method", "config", "n_total", "K", "p", "bias_mean",
                "bias_sd", "se_mean", "coverage", "ci_width_mean",
                "n_success")
  missing <- setdiff(required, names(summaries))
  expect_equal(missing, character(0))
  expect_true(all(summaries$n_success > 0L))
  expect_true(all(is.finite(summaries$coverage)))

  summaries$abs_bias <- abs(summaries$bias_mean)
  summaries$bias_over_se <- summaries$abs_bias / summaries$se_mean
  summaries$se_over_emp_sd <- summaries$se_mean / summaries$bias_sd
  summaries$normal_coverage_with_bias <- with(
    summaries,
    pnorm((1.96 * se_mean - bias_mean) / bias_sd) -
      pnorm((-1.96 * se_mean - bias_mean) / bias_sd)
  )
  summaries$normal_coverage_centered <- with(
    summaries,
    pnorm((1.96 * se_mean) / bias_sd) -
      pnorm((-1.96 * se_mean) / bias_sd)
  )
  summaries
}

read_n5000_c2_raw_results <- function() {
  results_dir <- find_coverage_results_dir()
  raw_files <- list.files(
    results_dir,
    pattern = "^n5000_K[0-9]+_p[0-9]+_C2_.*[.]csv$",
    full.names = TRUE
  )
  raw_files <- raw_files[!grepl("_summary[.]csv$", raw_files)]
  if (length(raw_files) == 0L) {
    skip(sprintf("no completed n5000 C2 raw CSVs found in %s", results_dir))
  }

  parse_setting <- function(path) {
    name <- basename(path)
    match <- regexec("^n([0-9]+)_K([0-9]+)_p([0-9]+)_(C[0-9])_", name)
    parts <- regmatches(name, match)[[1L]]
    if (length(parts) != 5L) {
      stop(sprintf("Could not parse setting from result filename: %s", name),
           call. = FALSE)
    }
    list(result_file = name, p = as.integer(parts[4L]), config_file = parts[5L])
  }

  rows <- lapply(raw_files, function(path) {
    df <- read.csv(path, stringsAsFactors = FALSE)
    meta <- parse_setting(path)
    df$result_file <- meta$result_file
    df$p <- meta$p
    df$config_file <- meta$config_file
    df
  })
  raw <- do.call(rbind, rows)

  required <- c("sim_id", "method", "K", "p", "bias", "se", "coverage")
  missing <- setdiff(required, names(raw))
  expect_equal(missing, character(0))
  raw
}

summarise_c2_residual_bias_vs_oracle <- function(raw) {
  oracle <- raw[raw$method == "oracle_dr",
                c("result_file", "K", "p", "sim_id", "bias"),
                drop = FALSE]
  names(oracle)[names(oracle) == "bias"] <- "oracle_bias"
  merged <- merge(
    raw,
    oracle,
    by = c("result_file", "K", "p", "sim_id"),
    all.x = FALSE,
    all.y = FALSE
  )
  merged$residual_vs_oracle <- merged$bias - merged$oracle_bias

  key <- with(merged, paste(K, p, method, sep = "|"))
  rows <- by(merged, key, function(sub) {
    data.frame(
      K = sub$K[1L],
      p = sub$p[1L],
      method = sub$method[1L],
      n_success = nrow(sub),
      bias_mean = mean(sub$bias),
      oracle_bias_mean = mean(sub$oracle_bias),
      residual_mean = mean(sub$residual_vs_oracle),
      residual_sd = stats::sd(sub$residual_vs_oracle),
      coverage = mean(sub$coverage),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, lapply(rows, identity))
  out[order(out$K, out$p, out$method), , drop = FALSE]
}

test_that("C2 low coverage diagnostics separate center bias from SE scale", {
  summaries <- read_n5000_coverage_summaries()
  c2 <- summaries[summaries$config == "C2", , drop = FALSE]
  expect_gt(nrow(c2), 0L)

  focus_methods <- c("federated_dr", "one_round_crossfit", "two_round_crossfit")
  focus <- c2[c2$method %in% focus_methods, , drop = FALSE]
  expect_gt(nrow(focus), 0L)

  diagnostic_cols <- c(
    "K", "p", "method", "coverage", "bias_mean", "bias_sd", "se_mean",
    "bias_over_se", "se_over_emp_sd", "normal_coverage_with_bias",
    "normal_coverage_centered", "n_success"
  )
  focus <- focus[order(focus$K, focus$p, focus$method), diagnostic_cols, drop = FALSE]
  rownames(focus) <- NULL
  cat("\n[c2-bias-validation] Focused C2 rows:\n")
  print(focus, row.names = FALSE, digits = 4)

  low <- focus[focus$coverage < 0.93, , drop = FALSE]
  if (nrow(low) == 0L) {
    cat("\n[c2-bias-validation] No focused C2 rows below 0.93 coverage.\n")
    succeed()
    return(invisible(NULL))
  }

  cat("\n[c2-bias-validation] Focused C2 rows below 0.93 coverage:\n")
  print(low, row.names = FALSE, digits = 4)

  expect_true(
    all(low$bias_over_se > 0.4),
    info = "Low-coverage C2 rows should have material center bias relative to mean SE."
  )
  expect_true(
    all(low$se_over_emp_sd > 0.9 & low$se_over_emp_sd < 1.25),
    info = "Low-coverage C2 rows do not look like pure SE-collapse rows."
  )
  expect_true(
    all(low$normal_coverage_centered > 0.93),
    info = "With the empirical SD and reported SE but zero bias, normal-approx coverage should be near nominal."
  )
  expect_true(
    mean(abs(low$normal_coverage_with_bias - low$coverage)) < 0.06,
    info = "Normal approximation with the observed bias should explain the coverage drop."
  )
  expect_true(
    all(low$normal_coverage_with_bias < low$normal_coverage_centered - 0.03),
    info = "The observed center bias should materially reduce coverage relative to a centered estimator."
  )
})

test_that("C2 residual bias is visible after subtracting oracle simulation noise", {
  summaries <- read_n5000_coverage_summaries()
  raw <- read_n5000_c2_raw_results()
  residuals <- summarise_c2_residual_bias_vs_oracle(raw)

  diagnostic_cols <- c(
    "K", "p", "method", "coverage", "bias_mean", "oracle_bias_mean",
    "residual_mean", "residual_sd", "n_success"
  )
  residual_print <- residuals[residuals$method %in% c(
    "target_only", "federated_dr", "pooled_dr", "tilted_aipw",
    "oracle_dr", "one_round_crossfit", "two_round_crossfit"
  ), diagnostic_cols, drop = FALSE]
  rownames(residual_print) <- NULL
  cat("\n[c2-bias-validation] C2 residual bias after subtracting oracle bias:\n")
  print(residual_print, row.names = FALSE, digits = 4)

  c2 <- summaries[summaries$config == "C2", , drop = FALSE]
  low <- c2[c2$coverage < 0.93 & c2$method != "oracle_dr",
            c("K", "p", "method", "coverage"), drop = FALSE]
  if (nrow(low) == 0L) {
    cat("\n[c2-bias-validation] No non-oracle C2 rows below 0.93 coverage.\n")
    succeed()
    return(invisible(NULL))
  }

  low_with_residual <- merge(
    low,
    residuals[, c("K", "p", "method", "residual_mean", "residual_sd")],
    by = c("K", "p", "method"),
    all.x = TRUE,
    all.y = FALSE
  )
  low_with_residual <- low_with_residual[
    order(low_with_residual$K, low_with_residual$p, low_with_residual$method),
    , drop = FALSE
  ]
  rownames(low_with_residual) <- NULL
  cat("\n[c2-bias-validation] Low-coverage C2 rows with oracle residuals:\n")
  print(low_with_residual, row.names = FALSE, digits = 4)

  expect_true(
    all(abs(low_with_residual$residual_mean) > 0.004),
    info = "Low-coverage non-oracle C2 rows should retain material method-specific bias after subtracting oracle simulation noise."
  )
  expect_true(
    all(is.finite(low_with_residual$residual_mean)) &&
      all(is.finite(low_with_residual$residual_sd)) &&
      all(low_with_residual$residual_sd >= 0),
    info = "Residual summaries should be finite and non-negative where applicable."
  )
})

test_that("C2 oracle benchmark validates that the summary pipeline is not globally broken", {
  summaries <- read_n5000_coverage_summaries()
  c2 <- summaries[summaries$config == "C2", , drop = FALSE]
  comparator_methods <- c("pooled_dr", "tilted_aipw", "oracle_dr")
  comparators <- c2[c2$method %in% comparator_methods, , drop = FALSE]
  expect_gt(nrow(comparators), 0L)

  diagnostic_cols <- c("K", "p", "method", "coverage", "bias_mean",
                       "bias_sd", "se_mean", "se_over_emp_sd", "n_success")
  comparators <- comparators[order(comparators$K, comparators$p, comparators$method),
                             diagnostic_cols, drop = FALSE]
  rownames(comparators) <- NULL
  cat("\n[c2-bias-validation] C2 comparator rows:\n")
  print(comparators, row.names = FALSE, digits = 4)

  oracle <- comparators[comparators$method == "oracle_dr", , drop = FALSE]
  expect_true(all(abs(oracle$bias_mean) < 0.002),
              info = "oracle_dr should be effectively centered in completed C2 summaries.")
  expect_true(all(oracle$coverage >= 0.95),
              info = "oracle_dr should remain near or above nominal in completed C2 summaries.")

  non_oracle <- comparators[comparators$method != "oracle_dr", , drop = FALSE]
  expect_true(all(is.finite(non_oracle$coverage)))
  expect_true(all(is.finite(non_oracle$bias_mean)))
})
