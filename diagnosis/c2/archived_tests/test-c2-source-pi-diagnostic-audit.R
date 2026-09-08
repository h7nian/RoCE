.c2_source_pi_audit_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_C2_SOURCE_PI_AUDIT", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-source-pi-diagnostic-audit", filter, fixed = TRUE)
}

.c2_source_pi_weight <- function(Z, gamma) {
  exp(-as.numeric(cbind(1, as.matrix(Z)) %*% as.numeric(gamma)))
}

test_that("c2 source diagnostics use source-conditional treatment propensity", {
  skip_if_not(.c2_source_pi_audit_enabled(),
              message = "set ROCE_C2_SOURCE_PI_AUDIT=1 or run with --filter c2-source-pi-diagnostic-audit")

  set.seed(71041)
  data <- generate_simulation_data(
    n_total = 8000L, K = 3L, p = 10L,
    config = "C2",
    estimand_type = "superpopulation",
    site_allocation = "model",
    transform_type = "mild",
    outcome_type = "binary",
    heterogeneity_type = "none",
    shift_strength = 0.5,
    dgp_type = "roce",
    warn_ignored = FALSE
  )

  calculate_site_probabilities <- c2_roce_function("calculate_site_probabilities")
  Z_true <- as.matrix(data$Z_site_true)
  probs <- calculate_site_probabilities(Z_true, data$gamma_params, data$K)

  for (source_name in c("s1", "s2", "s3")) {
    source_idx <- c2_source_index(source_name)
    source_prob <- probs[[paste0("p_s", source_idx, "_0")]] +
      probs[[paste0("p_s", source_idx, "_1")]]
    expected_pi <- probs[[paste0("p_s", source_idx, "_1")]] /
      pmax(source_prob, .Machine$double.eps)
    source_pi <- c2_true_source_treatment_propensity(data, source_name, A_val = 1L)

    expect_equal(source_pi, expected_pi, tolerance = 1e-12)
    expect_gt(mean(abs(source_pi - data$p_treat_true)), 1e-4)

    gamma <- c2_true_source_calibration_gamma(data, source_name, A_val = 1L)
    calibrated_intercept <- mean(
      source_prob * source_pi * .c2_source_pi_weight(Z_true, gamma)
    ) / mean(source_prob)
    global_pi_intercept <- mean(
      source_prob * data$p_treat_true * .c2_source_pi_weight(Z_true, gamma)
    ) / mean(source_prob)

    expect_equal(calibrated_intercept, 1, tolerance = 1e-10)
    expect_gt(abs(global_pi_intercept - 1), 1e-4)
  }
})
