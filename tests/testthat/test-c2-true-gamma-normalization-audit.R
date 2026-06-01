.c2_true_gamma_audit_enabled <- function() {
  env_enabled <- Sys.getenv("FACEHD_C2_TRUE_GAMMA_AUDIT", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-true-gamma-normalization-audit", filter, fixed = TRUE)
}

source({
  helper_path <- file.path("diagnosis", "c2", "c2_true_gamma_utils.R")
  if (file.exists(helper_path)) helper_path else file.path("..", "..", helper_path)
})

.c2_true_gamma_weight <- function(Z, gamma) {
  exp(-as.numeric(cbind(1, as.matrix(Z)) %*% as.numeric(gamma)))
}

.c2_true_gamma_model_intercept <- function(data, source_name, gamma,
                                           A_val = 1L) {
  source_idx <- c2_source_index(source_name)
  calculate_site_probabilities <- c2_facehd_function("calculate_site_probabilities")
  Z_true <- as.matrix(if (is.null(data$Z_site_true)) data$Z_site else data$Z_site_true)
  probs <- calculate_site_probabilities(Z_true, data$gamma_params, data$K)
  source_prob <- probs[[paste0("p_s", source_idx, "_0")]] +
    probs[[paste0("p_s", source_idx, "_1")]]
  source_arm_prob <- probs[[paste0("p_s", source_idx, "_", A_val)]]
  mean(source_arm_prob * .c2_true_gamma_weight(Z_true, gamma)) / mean(source_prob)
}

test_that("c2 true gamma helper converts joint DGP gamma to source-conditional calibration gamma", {
  skip_if_not(.c2_true_gamma_audit_enabled(),
              message = "set FACEHD_C2_TRUE_GAMMA_AUDIT=1 or run with --filter c2-true-gamma-normalization-audit")

  set.seed(71031)
  data <- generate_simulation_data(
    n_total = 6000L, K = 3L, p = 10L,
    config = "C2",
    estimand_type = "superpopulation",
    site_allocation = "model",
    transform_type = "mild",
    outcome_type = "binary",
    heterogeneity_type = "none",
    shift_strength = 0.5,
    dgp_type = "facehd",
    warn_ignored = FALSE
  )

  for (source_name in c("s1", "s2", "s3")) {
    raw_gamma <- c2_true_source_joint_gamma(data, source_name, A_val = 1L)
    calibrated_gamma <- c2_true_source_calibration_gamma(
      data, source_name, A_val = 1L
    )
    expected_shift <- c2_true_source_log_target_over_source(
      data, source_name, ratio = "population"
    )

    expect_equal(calibrated_gamma[1L] - raw_gamma[1L], expected_shift,
                 tolerance = 1e-12)

    raw_intercept <- .c2_true_gamma_model_intercept(
      data, source_name, raw_gamma, A_val = 1L
    )
    calibrated_intercept <- .c2_true_gamma_model_intercept(
      data, source_name, calibrated_gamma, A_val = 1L
    )

    expect_equal(calibrated_intercept, 1, tolerance = 1e-10)
    expect_gt(abs(raw_intercept - 1), 1e-3)
    expect_lt(abs(calibrated_intercept - 1),
              abs(raw_intercept - 1) * 1e-6)
  }
})
