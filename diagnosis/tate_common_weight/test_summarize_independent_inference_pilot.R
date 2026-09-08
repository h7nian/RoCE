#!/usr/bin/env Rscript

library(testthat)
source("diagnosis/tate_common_weight/run_independent_inference_pilot.R")
source("diagnosis/tate_common_weight/summarize_independent_inference_pilot.R")

make_summary_fixture <- function(n = 2L, constant = FALSE) {
  rhos <- .inference_pilot_rhos()
  raw <- inference <- list(); index <- 0L
  for (rho in rhos) {
    truth <- 10 + rho
    soft_error <- if (constant) rep(1, n) else c(1, 3)[seq_len(n)]
    target_error <- rep(2, n)
    for (i in seq_len(n)) {
      sim <- 10000L + i
      index <- index + 1L
      raw[[2L * index - 1L]] <- data.frame(
        sim_id = sim, rho = rho, method = "one_round_crossfit_ate",
        estimate = truth + soft_error[[i]], truth = truth, coverage = i == 1L
      )
      raw[[2L * index]] <- data.frame(
        sim_id = sim, rho = rho, method = "target_only_ate",
        estimate = truth + target_error[[i]], truth = truth, coverage = TRUE
      )
      inference[[index]] <- data.frame(
        sim_id = sim, rho = rho, method = "one_round_crossfit_ate",
        estimate = truth + soft_error[[i]], truth = truth,
        analytic_se = 2, analytic_ci_lower = truth - 1,
        analytic_ci_upper = truth + 3, analytic_coverage = i == 1L,
        fixed_weight_bootstrap_se = 3,
        weight_relearn_bootstrap_se = 4,
        weight_relearn_ci_lower = truth - 2,
        weight_relearn_ci_upper = truth + 6,
        weight_relearn_coverage = TRUE,
        weight_uncertainty_ratio = 4 / 3, n_weight_bootstrap = 1000L,
        inference_status = "diagnostic_only"
      )
    }
  }
  list(raw = do.call(rbind, raw), inference = do.call(rbind, inference))
}

test_that("paired summary algebra and exact binomial intervals are explicit", {
  fixture <- make_summary_fixture()
  result <- roce_summarize_independent_inference(fixture$raw, fixture$inference)
  expect_equal(nrow(result), 6L)
  row <- result[result$rho == 0, ]
  expect_equal(row$mean_error, 2)
  expect_equal(row$mean_error_mcse, 1)
  expect_equal(row$empirical_sd, sqrt(2))
  expect_equal(row$rmse, sqrt(5))
  expect_equal(row$target_rmse, 2)
  expect_equal(row$paired_mse_difference_soft_minus_target, 1)
  expect_equal(row$paired_mse_difference_mcse, 4)
  expect_equal(row$mean_analytic_se, 2)
  expect_equal(row$mean_fixed_weight_bootstrap_se, 3)
  expect_equal(row$mean_weight_relearn_bootstrap_se, 4)
  expect_equal(row$analytic_se_to_empirical_sd, sqrt(2))
  expect_equal(row$mean_analytic_ci_length, 4)
  expect_equal(row$mean_weight_relearn_ci_length, 8)
  interval <- binom.test(1, 2)$conf.int
  expect_equal(c(row$analytic_coverage_exact_lower,
                 row$analytic_coverage_exact_upper), as.numeric(interval))
  expect_equal(row$weight_relearn_coverage_minus_analytic, 0.5)
  expect_equal(row$weight_relearn_coverage_minus_analytic_mcse, 0.5)
  expect_true(is.na(row$warning_count))
  expect_false(row$warning_capture_complete)
  expect_false(row$inference_validated)
  expect_true(row$statistical_review_required)
})

test_that("paired relearn-minus-analytic coverage handles concordance", {
  fixture <- make_summary_fixture()
  fixture$inference$weight_relearn_coverage <-
    fixture$inference$analytic_coverage
  result <- roce_summarize_independent_inference(fixture$raw, fixture$inference)
  expect_true(all(result$weight_relearn_coverage_minus_analytic == 0))
  expect_true(all(result$weight_relearn_coverage_minus_analytic_mcse == 0))
})

test_that("N=1 and zero empirical SD yield NA rather than Inf", {
  one <- make_summary_fixture(1L)
  one_result <- roce_summarize_independent_inference(one$raw, one$inference)
  expect_true(all(is.na(one_result$empirical_sd)))
  expect_true(all(is.na(one_result$mean_error_mcse)))
  expect_true(all(is.na(one_result$analytic_se_to_empirical_sd)))
  expect_true(all(is.na(one_result$paired_mse_difference_mcse)))
  expect_true(all(is.na(one_result$analytic_coverage_minus_target_mcse)))
  expect_true(all(is.na(
    one_result$weight_relearn_coverage_minus_analytic_mcse
  )))

  constant <- make_summary_fixture(2L, constant = TRUE)
  constant_result <- roce_summarize_independent_inference(
    constant$raw, constant$inference
  )
  expect_true(all(constant_result$empirical_sd == 0))
  ratio_fields <- grep("_se_to_empirical_sd$", names(constant_result), value = TRUE)
  expect_true(all(vapply(constant_result[ratio_fields], function(x) all(is.na(x)), logical(1L))))
})

test_that("pairing, truth, and coverage schema fail closed", {
  fixture <- make_summary_fixture()
  missing <- fixture$raw[!(fixture$raw$rho == 0 & fixture$raw$sim_id == 10002), ]
  expect_error(roce_summarize_independent_inference(missing, fixture$inference),
               "incomplete or invalid")
  bad_truth <- fixture$raw
  index <- bad_truth$rho == 0 & bad_truth$sim_id == 10001 &
    bad_truth$method == "target_only_ate"
  bad_truth$truth[index] <- bad_truth$truth[index] + 1
  expect_error(roce_summarize_independent_inference(bad_truth, fixture$inference),
               "disagree on truth")
  bad_coverage <- fixture$raw
  bad_coverage$coverage <- as.integer(bad_coverage$coverage)
  expect_error(roce_summarize_independent_inference(bad_coverage, fixture$inference),
               "incomplete or invalid")
})

test_that("missing prefix fails without publishing a complete summary", {
  root <- tempfile("independent-summary-root-")
  dir.create(root)
  write.csv(roce_inference_pilot_manifest(), file.path(root, "manifest.csv"),
            row.names = FALSE)
  output <- file.path(root, "summary_n001")
  expect_error(
    summarize_independent_inference_pilot_main(c(root, "1", output)),
    "checkpoint incomplete"
  )
  expect_false(file.exists(output) || dir.exists(output))
})
