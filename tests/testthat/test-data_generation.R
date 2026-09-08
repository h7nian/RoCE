# test-data_generation.R - Tests for data generation functions

library(testthat)

# Package loaded by helper-load.R (all functions available via RoCE namespace)

test_that("generate_base_covariates creates correct dimensions", {
  set.seed(42)
  n <- 100
  p <- 5
  X <- generate_base_covariates(n, p)

  expect_equal(nrow(X), n)
  expect_equal(ncol(X), p)
})

test_that("generate_base_covariates produces standardized output", {
  set.seed(42)
  n <- 1000
  p <- 4
  X <- generate_base_covariates(n, p)

  # Check approximately zero mean (within tolerance)
  col_means <- colMeans(X)
  expect_true(all(abs(col_means) < 0.1), 
              info = "Column means should be approximately zero")

  # Check approximately unit variance
  col_sds <- apply(X, 2, sd)
  expect_true(all(abs(col_sds - 1) < 0.1),
              info = "Column SDs should be approximately 1")
})

test_that("generate_base_covariates handles different dimensions", {
  set.seed(42)
  n <- 1000
  p <- 4
  
  X <- generate_base_covariates(n, p)
  
  # Output should be a numeric matrix with correct dimensions
  expect_true(is.matrix(X))
  expect_equal(nrow(X), n)
  expect_equal(ncol(X), p)
})

test_that("transform_covariates applies nonlinear transformations", {
  set.seed(42)
  n <- 100
  p <- 4
  X <- generate_base_covariates(n, p)
  X_dagger <- transform_covariates(X)
  
  expect_equal(nrow(X_dagger), n)
  expect_equal(ncol(X_dagger), p)
  
  # Transformed covariates should be different from original
  expect_false(all(X == X_dagger))
})

test_that("generate_simulation_data creates valid structure", {
  set.seed(42)
  n_total <- 200
  K <- 2  # 2 source sites
  p <- 4
  
  data <- generate_simulation_data(n_total, K, p, config = "C1",
                                   estimand_type = "sample", dgp_type = "roce")
  
  # Check that data has required components
  expect_true("X" %in% names(data))
  expect_true("X_dagger" %in% names(data))
  expect_true("R" %in% names(data))
  expect_true("A" %in% names(data))
  expect_true("Y" %in% names(data))
  expect_true("Z_site_true" %in% names(data))
  expect_true("W_outcome_true" %in% names(data))
  expect_true("Z_site" %in% names(data))
  expect_true("W_outcome" %in% names(data))
  
  # Check dimensions
  expect_equal(nrow(data$X), n_total)
  expect_equal(nrow(data$Z_site_true), n_total)
  expect_equal(nrow(data$W_outcome_true), n_total)
  expect_equal(nrow(data$Z_site), n_total)
  expect_equal(nrow(data$W_outcome), n_total)
  expect_equal(length(data$R), n_total)
  expect_equal(length(data$A), n_total)
  expect_equal(length(data$Y), n_total)
  expect_equal(data$Y, ifelse(data$A == 1, data$Y_1, data$Y_0))
  
  # Check treatment is binary
  expect_true(all(data$A %in% c(0, 1)))
  
  # Check outcome is binary (for logit link)
  expect_true(all(data$Y %in% c(0, 1)))
})

test_that("default dgp_type is 'face' at the simulation entry points", {
  # The package default DGP is the FACE negative-transfer DGP.
  expect_identical(formals(generate_simulation_data)$dgp_type, "face")
  expect_identical(formals(run_single_simulation)$dgp_type,    "face")
  expect_identical(formals(run_simulation_study)$dgp_type,     "face")
  expect_identical(formals(validate_simulation_params)$dgp_type, "face")

  # Calling the dispatcher without dgp_type produces FACE-DGP data.
  set.seed(1)
  data <- generate_simulation_data(n_total = 300, K = 2, p = 10,
                                   config = "C1", estimand_type = "sample",
                                   warn_ignored = FALSE)
  expect_equal(data$dgp_type, "face")
})

test_that("split_data_by_site correctly separates sites", {
  set.seed(42)
  n_total <- 300
  K <- 2
  p <- 4
  
  data <- generate_simulation_data(n_total, K, p, config = "C1",
                                   estimand_type = "sample", dgp_type = "roce")
  data_split <- split_data_by_site(data)
  
  # Should have target + K source sites
  expect_equal(length(data_split), K + 1)
  expect_true("t" %in% names(data_split))
  
  # Check source site names
  for (j in 1:K) {
    expect_true(paste0("s", j) %in% names(data_split))
  }
  
  # Check that total observations sum to original n
  total_n <- sum(sapply(data_split, function(x) x$n))
  expect_equal(total_n, n_total)

  for (site in names(data_split)) {
    expect_true("Z_site_true" %in% names(data_split[[site]]))
    expect_true("W_outcome_true" %in% names(data_split[[site]]))
  }
})

test_that("RoCE configs change fitted bases without changing the DGP", {
  set.seed(42)
  n_total <- 200
  K <- 2
  p <- 4
  
  data_c1 <- generate_simulation_data(n_total, K, p, config = "C1",
                                      estimand_type = "sample", dgp_type = "roce")
  set.seed(42)
  data_c2 <- generate_simulation_data(n_total, K, p, config = "C2",
                                      estimand_type = "sample", dgp_type = "roce")
  set.seed(42)
  data_c3 <- generate_simulation_data(n_total, K, p, config = "C3",
                                      estimand_type = "sample", dgp_type = "roce")
  set.seed(42)
  data_c4 <- generate_simulation_data(n_total, K, p, config = "C4",
                                      estimand_type = "sample", dgp_type = "roce")

  # The true data-generating process is fixed across configurations.
  expect_identical(data_c1$R, data_c2$R)
  expect_identical(data_c1$R, data_c3$R)
  expect_identical(data_c1$R, data_c4$R)
  expect_identical(data_c1$A, data_c2$A)
  expect_identical(data_c1$A, data_c3$A)
  expect_identical(data_c1$A, data_c4$A)
  expect_identical(data_c1$Y, data_c2$Y)
  expect_identical(data_c1$Y, data_c3$Y)
  expect_identical(data_c1$Y, data_c4$Y)
  expect_equal(data_c1$mu1_true, data_c2$mu1_true)
  expect_equal(data_c1$mu1_true, data_c3$mu1_true)
  expect_equal(data_c1$mu1_true, data_c4$mu1_true)
  expect_equal(data_c1$Z_site_true, data_c2$Z_site_true)
  expect_equal(data_c1$Z_site_true, data_c3$Z_site_true)
  expect_equal(data_c1$Z_site_true, data_c4$Z_site_true)
  expect_equal(data_c1$W_outcome_true, data_c2$W_outcome_true)
  expect_equal(data_c1$W_outcome_true, data_c3$W_outcome_true)
  expect_equal(data_c1$W_outcome_true, data_c4$W_outcome_true)

  # Configs only control the fitted working bases exposed to estimators.
  expect_equal(data_c1$Z_site, data_c2$Z_site)
  expect_false(isTRUE(all.equal(data_c1$W_outcome, data_c2$W_outcome)))
  expect_equal(data_c1$W_outcome, data_c3$W_outcome)
  expect_false(isTRUE(all.equal(data_c1$Z_site, data_c3$Z_site)))
  expect_equal(data_c2$W_outcome, data_c4$W_outcome)
  expect_equal(data_c3$Z_site, data_c4$Z_site)
  expect_false(isTRUE(all.equal(data_c1$Z_site, data_c4$Z_site)))
  expect_false(isTRUE(all.equal(data_c1$W_outcome, data_c4$W_outcome)))
})

test_that("shift_strength affects RoCE site model under model allocation", {
  K <- 3
  p <- 10
  n_ref <- 2000

  # Fixed covariates so differences come only from gamma parameters
  set.seed(123)
  X_ref <- generate_base_covariates(n_ref, p)
  X_dag <- transform_covariates(X_ref, transform_type = "mild")
  Z_site <- X_dag

  # With the same RNG seed, changing shift_strength should change gamma
  # and therefore the site-probability mapping.
  set.seed(999)
  gamma_1 <- generate_site_model_parameters(K, p, shift_strength = 1.0)
  set.seed(999)
  gamma_2 <- generate_site_model_parameters(K, p, shift_strength = 2.0)

  # At least one gamma vector should differ materially.
  diffs <- vapply(names(gamma_1), function(nm) {
    max(abs(gamma_1[[nm]] - gamma_2[[nm]]))
  }, numeric(1))
  expect_true(any(diffs > 1e-6), info = "gamma should change with shift_strength")

  probs_1 <- calculate_site_probabilities(Z_site, gamma_1, K)
  probs_2 <- calculate_site_probabilities(Z_site, gamma_2, K)
  expect_true(mean(abs(probs_1$p_target - probs_2$p_target)) > 1e-4,
              info = "P(R=t|X) should change with shift_strength")
})

test_that("superpopulation truth respects independent target covariate distribution", {
  # Under independent allocation, each site has its own covariate distribution.
  # The target-specific estimand E[Y(1) | R=t] should therefore integrate over
  # the target site's covariate distribution, which differs from the pooled
  # (generate_base_covariates) reference used for uniform/balanced allocations.
  res_uniform <- calculate_superpopulation_truth(
    p = 4,
    K = 2,
    config = "C4",
    n_ref = 20000,
    ref_seed = 12345,
    transform_type = "none",
    site_allocation = "uniform",
    outcome_type = "binary",
    shift_strength = 1.0
  )

  res_indep <- calculate_superpopulation_truth(
    p = 4,
    K = 2,
    config = "C4",
    n_ref = 20000,
    ref_seed = 12345,
    transform_type = "none",
    site_allocation = "independent",
    outcome_type = "binary",
    shift_strength = 1.0
  )

  expect_true(abs(res_uniform$mu1_superpop - res_indep$mu1_superpop) > 5e-4,
              info = "Independent allocation should use target covariate distribution, not pooled")
})

# =============================================================================
# FACE Paper DGP Tests (Han et al. JASA 2023, Section 5.1)
# =============================================================================

test_that("generate_skewed_normal produces correct dimensions and moments", {
  set.seed(42)
  n <- 5000
  X_sym <- generate_skewed_normal(n, kappa = 0, phi = 1, nu = 0)
  X_skew <- generate_skewed_normal(n, kappa = 0, phi = 1, nu = 5)
  
  expect_length(X_sym, n)
  expect_length(X_skew, n)
  
  # Symmetric (nu=0) should be approximately normal: near-zero skewness
  expect_true(abs(mean(X_sym)) < 0.1)
  
  # Skewed (nu=5) should have positive skewness
  skew <- mean((X_skew - mean(X_skew))^3) / sd(X_skew)^3
  expect_true(skew > 0.1, info = "nu=5 should produce positive skewness")
})

test_that("generate_face_covariates returns correct structure", {
  set.seed(42)
  K <- 3
  p <- 10
  n_target <- 200
  n_source <- 150
  
  cov_data <- generate_face_covariates(n_target, rep(n_source, K), p,
                                              kappa = 0.125, nu_source_max = 0.2)
  
  expect_true(is.matrix(cov_data$X))
  expect_equal(ncol(cov_data$X), p)
  expect_equal(length(cov_data$R), nrow(cov_data$X))
  
  # Should have target (R="t") and K source sites (R="s1".."sK")
  expect_true("t" %in% cov_data$R)
  expect_true(all(paste0("s", 1:K) %in% cov_data$R))
  
  # Target site should have n_target observations
  expect_equal(sum(cov_data$R == "t"), n_target)
  # Each source site should have n_source observations
  for (j in 1:K) {
    expect_equal(sum(cov_data$R == paste0("s", j)), n_source)
  }
})

test_that("calculate_face_propensity returns valid probabilities", {
  set.seed(42)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  params <- get_face_ps_parameters(p)
  
  pi_vec <- calculate_face_propensity(X, params$alpha1, params$alpha2)
  
  expect_length(pi_vec, n)
  expect_true(all(pi_vec > 0 & pi_vec < 1))
  # Should be clipped to positivity bounds
  expect_true(all(pi_vec >= POSITIVITY_LOWER))
  expect_true(all(pi_vec <= POSITIVITY_UPPER))
})

test_that("build_face_ate_map builds correct ATE structure", {
  K <- 5
  ate_map_0 <- build_face_ate_map(K, ate_deviation = 0.0, n_deviated = 0)
  ate_map_d <- build_face_ate_map(K, ate_deviation = 2.0, n_deviated = 3)
  
  # All informative (no deviated)
  expect_equal(length(ate_map_0), K + 1)  # target + K sources
  expect_true(all(ate_map_0 == 3.0))
  expect_equal(names(ate_map_0)[1], "t")
  
  # 3 deviated sites
  expect_equal(ate_map_d[["t"]], 3.0)
  deviated_vals <- ate_map_d[paste0("s", 1:3)]
  expect_true(all(deviated_vals == 5.0))  # 3.0 + 2.0
  informative_vals <- ate_map_d[paste0("s", 4:5)]
  expect_true(all(informative_vals == 3.0))

  # Binary base ATE is a first-class parameter (continuous default unchanged)
  ate_map_bin <- build_face_ate_map(K, base_ate = FACE_BINARY_ATE_TARGET)
  expect_true(all(ate_map_bin == FACE_BINARY_ATE_TARGET))
})

test_that("generate_face_data returns valid structure", {
  set.seed(42)
  n_total <- 500
  K <- 2
  p <- 10
  
  data <- generate_face_data(n_total, K, p, config = "C1",
                                   estimand_type = "superpopulation")
  
  # Basic structure checks
  expect_true("X" %in% names(data))
  expect_true("R" %in% names(data))
  expect_true("A" %in% names(data))
  expect_true("Y" %in% names(data))
  expect_true("W_outcome" %in% names(data))
  expect_true("Z_site" %in% names(data))
  expect_true("dgp_type" %in% names(data))
  expect_equal(data$dgp_type, "face")
  expect_equal(data$outcome_type, "continuous")
  
  # Dimensions
  expect_equal(nrow(data$X), n_total)
  expect_equal(ncol(data$X), p)
  expect_equal(length(data$R), n_total)
  expect_equal(length(data$A), n_total)
  expect_equal(length(data$Y), n_total)
  expect_equal(data$Y, ifelse(data$A == 1, data$Y_1, data$Y_0))
  
  # Treatment binary
  expect_true(all(data$A %in% c(0, 1)))
  
  # Continuous outcomes should NOT be binary
  expect_false(all(data$Y %in% c(0, 1)))
  
  # RoCE specific params should be NULL
  expect_null(data$gamma_params)
  expect_null(data$beta1_true)
  expect_null(data$beta0_true)
})

test_that("generate_face_data supports binary outcomes", {
  set.seed(42)
  data_bin <- generate_face_data(500, K = 2, p = 10, config = "C1",
                                 outcome_type = "binary")
  
  expect_equal(data_bin$outcome_type, "binary")
  expect_equal(data_bin$dgp_type, "face")
  
  # Binary outcomes must be 0/1
  expect_true(all(data_bin$Y %in% c(0, 1)))
  expect_true(all(data_bin$Y_1 %in% c(0, 1)))
  expect_true(all(data_bin$Y_0 %in% c(0, 1)))
  
  # True values should be between 0 and 1 (probabilities)
  expect_true(data_bin$mu1_true > 0 && data_bin$mu1_true < 1)
  expect_true(data_bin$mu0_true > 0 && data_bin$mu0_true < 1)
})

test_that("FACE binary outcome is calibrated (no expit saturation), data matches estimand", {
  # Standardizing the logit signal keeps both arms away from 0/1, so the outcome
  # regression's treated-arm minority class is healthy (the saturation that made
  # the raw Δ_T = 3 log-odds binary DGP degenerate). p = 50 is the dimension that
  # previously saturated. estimand_type = "sample": truth uses the realized
  # target X, so the generated potential outcomes must match it up to Bernoulli
  # sampling noise.
  set.seed(7)
  d <- generate_face_data(n_total = 8000, K = 2, p = 50, config = "C1",
                          outcome_type = "binary", estimand_type = "sample")
  t_idx <- d$R == "t"
  expect_equal(mean(d$Y_1[t_idx]), d$mu1_true, tolerance = 0.05)
  expect_equal(mean(d$Y_0[t_idx]), d$mu0_true, tolerance = 0.05)
  # Both arms bounded away from saturation (vs. the old ~0.98 / ~0.83)
  expect_true(d$mu1_true < 0.9 && d$mu0_true > 0.1)
})

test_that("FACE binary superpopulation truth matches generated data", {
  set.seed(11)
  d <- generate_face_data(n_total = 20000, K = 2, p = 50, config = "C1",
                          outcome_type = "binary",
                          estimand_type = "superpopulation")
  t_idx <- d$R == "t"
  expect_equal(mean(d$Y_1[t_idx]), d$mu1_true, tolerance = 0.05)
  expect_equal(mean(d$Y_0[t_idx]), d$mu0_true, tolerance = 0.05)
})

test_that("FACE binary calibration and truth come from one reference population", {
  calibration <- RoCE:::get_face_binary_calibration(
    p = 4, n_ref = 2000L, ref_seed = 731L, config = "C2"
  )
  truth <- RoCE:::calculate_face_truth(
    p = 4, outcome_type = "binary", n_ref = 2000L, ref_seed = 731L, config = "C2"
  )
  expect_equal(calibration$mu1_superpop, truth$mu1_superpop, tolerance = 0)
  expect_equal(calibration$mu0_superpop, truth$mu0_superpop, tolerance = 0)
  expect_equal(truth$ate_superpop, truth$mu1_superpop - truth$mu0_superpop)
  correct <- RoCE:::calculate_face_truth(
    p = 4, outcome_type = "binary", n_ref = 2000L, ref_seed = 731L, config = "C1"
  )
  expect_false(isTRUE(all.equal(correct$mu1_superpop, truth$mu1_superpop)))
})

test_that("FACE continuous estimand is exactly FACE_ATE_TARGET and the binary patch is inert", {
  # The binary calibration must not perturb the continuous path: its estimand is
  # still exactly Δ_T, and the continuous RNG stream is unchanged (deterministic).
  set.seed(3)
  d_sup <- generate_face_data(2000, K = 2, p = 10, config = "C1",
                              outcome_type = "continuous",
                              estimand_type = "superpopulation")
  expect_equal(d_sup$mu1_true - d_sup$mu0_true, FACE_ATE_TARGET, tolerance = 1e-8)

  set.seed(3)
  d_samp <- generate_face_data(2000, K = 2, p = 10, config = "C1",
                               outcome_type = "continuous",
                               estimand_type = "sample")
  expect_equal(d_samp$mu1_true - d_samp$mu0_true, FACE_ATE_TARGET, tolerance = 1e-8)

  set.seed(3)
  d_rep <- generate_face_data(2000, K = 2, p = 10, config = "C1",
                              outcome_type = "continuous",
                              estimand_type = "sample")
  expect_identical(d_samp$Y, d_rep$Y)
})

test_that("generate_face_data uses one working basis and misspecifies the true mechanism", {
  generate <- function(config, strength = FACE_MISSPECIFICATION_STRENGTH) {
    set.seed(42)
    generate_face_data(300, K = 2, p = 10, config = config, outcome_type = "binary",
                       misspecification_strength = strength)
  }
  data_c1 <- generate("C1")
  data_c2 <- generate("C2")
  data_c3 <- generate("C3")
  data_c4 <- generate("C4")

  # Both nuisance models see the same quadratic basis in every configuration.
  for (data in list(data_c1, data_c2, data_c3, data_c4)) {
    expect_identical(data$Z_site, data$W_outcome)
    expect_equal(ncol(data$Z_site), 2L * ncol(data$X))
    expect_equal(data$Z_site, cbind(sweep(data$X, 2L, FACE_KAPPA, "-"), data$X^2))
  }
  expect_identical(data_c1$X, data_c2$X)
  expect_identical(data_c1$X_dagger, data_c3$X_dagger)

  # C2 changes the outcome mechanism (and its truth), C3 the treatment mechanism.
  expect_identical(data_c1$A, data_c2$A)
  expect_false(identical(data_c1$Y_1, data_c2$Y_1))
  expect_false(isTRUE(all.equal(data_c1$mu1_true, data_c2$mu1_true)))
  expect_false(identical(data_c1$A, data_c3$A))
  expect_equal(data_c1$mu1_true, data_c3$mu1_true)
  expect_false(identical(data_c1$A, data_c4$A))
  expect_equal(data_c2$mu1_true, data_c4$mu1_true)

  # Zero strength reproduces C1 exactly for every configuration.
  for (config in c("C2", "C3", "C4")) {
    data_zero <- generate(config, strength = 0)
    expect_identical(data_zero$Y, data_c1$Y)
    expect_identical(data_zero$A, data_c1$A)
    expect_equal(data_zero$mu1_true, data_c1$mu1_true)
  }

  # Transformed signal coordinates are standardized on the reference population.
  reference <- RoCE:::.face_reference_population(
    p = 10, kappa = FACE_KAPPA, config = "C2",
    misspecification_strength = FACE_MISSPECIFICATION_STRENGTH,
    outcome_type = "binary", n_ref = 20000L, ref_seed = 99999L
  )
  expect_length(reference$standardization$mean, FACE_SIGNAL_COORDINATES)
  expect_true(all(reference$standardization$sd > 0))
  expect_error(generate("C2", strength = 1.5), "misspecification_strength")
  expect_identical(generate("C1", strength = 0.5)$misspecification_strength, 0)
  expect_identical(generate("C2", strength = 0.5)$misspecification_strength, 0.5)
  expect_error(
    generate_simulation_data(
      n_total = 60, K = 1, p = 3, config = "C2", dgp_type = "roce",
      estimand_type = "sample", misspecification_strength = 0.5
    ),
    "applies only to dgp_type = 'face'"
  )
})

test_that("generate_simulation_data dispatches face DGP correctly", {
  set.seed(42)
  data <- generate_simulation_data(500, K = 2, p = 10, config = "C1",
                                   dgp_type = "face",
                                   outcome_type = "continuous")
  
  expect_equal(data$dgp_type, "face")
  expect_equal(data$outcome_type, "continuous")
  expect_true(nrow(data$X) > 0)
})

test_that("generate_simulation_data dispatches binary face DGP correctly", {
  set.seed(42)
  data <- generate_simulation_data(500, K = 2, p = 10, config = "C1",
                                   dgp_type = "face",
                                   outcome_type = "binary")
  
  expect_equal(data$dgp_type, "face")
  expect_equal(data$outcome_type, "binary")
  expect_true(all(data$Y %in% c(0, 1)))
  expect_true(nrow(data$X) > 0)
})

test_that("source site-treatment categories map to the intended treatment labels", {
  site_probs <- list(
    p_target = rep(0, 4),
    p_s1_0 = c(1, 0, 0, 0),
    p_s1_1 = c(0, 1, 0, 0),
    p_s2_0 = c(0, 0, 1, 0),
    p_s2_1 = c(0, 0, 0, 1)
  )

  assigned <- generate_site_treatment_assignments(
    n = 4L,
    site_probs = site_probs,
    K = 2L
  )

  expect_equal(assigned$R, c("s1", "s1", "s2", "s2"))
  expect_equal(assigned$A, c(0L, 1L, 0L, 1L))
})

test_that("get_face_outcome_parameters handles varying p", {
  # p < 4: uses min(4, p) non-zero coefficients

  params3 <- get_face_outcome_parameters(3)
  expect_length(params3$beta_linear, 3)
  expect_length(params3$beta_squared, 3)
  expect_true(all(params3$beta_linear != 0))

  # p >= 4: first 4 non-zero, rest zero
  params10 <- get_face_outcome_parameters(10)
  expect_length(params10$beta_linear, 10)
  expect_length(params10$beta_squared, 10)
  expect_true(all(params10$beta_linear[5:10] == 0))
  expect_true(all(params10$beta_squared[5:10] == 0))
  expect_true(all(params10$beta_linear[1:4] != 0))
})

test_that("get_face_ps_parameters handles varying p", {
  # p < 4: uses min(4, p) non-zero coefficients
  params3 <- get_face_ps_parameters(3)
  expect_length(params3$alpha1, 3)
  expect_length(params3$alpha2, 3)

  # p >= 4: first 4 non-zero for alpha1, alpha2 has only first non-zero
  params10 <- get_face_ps_parameters(10)
  expect_length(params10$alpha1, 10)
  expect_length(params10$alpha2, 10)
  expect_true(params10$alpha2[1] != 0)
  expect_true(all(params10$alpha2[2:10] == 0))
})

test_that("the both-arms deviation shifts only the deviated source, in both arms", {
  generate <- function(...) {
    set.seed(5)
    generate_face_data(300, K = 2, p = 10, config = "C1", outcome_type = "binary",
                       n_deviated_sites = 1L, ...)
  }
  null_data <- generate(ate_deviation = 0)
  treated_arm <- generate(ate_deviation = 1.5)
  both_arms <- generate(ate_deviation = 1.5, deviation_mechanism = "both_arms")

  expect_identical(both_arms$deviation_mechanism, "both_arms")
  expect_identical(treated_arm$deviation_mechanism, "treated_arm")
  expect_null(treated_arm$outcome_shift_map)
  expect_identical(both_arms$outcome_shift_map, c(t = 0, s1 = 1.5, s2 = 0))
  expect_identical(both_arms$ate_map, null_data$ate_map)
  expect_identical(treated_arm$ate_map[["s1"]], null_data$ate_map[["s1"]] + 1.5)
  # The target-population truth never depends on a source deviation.
  for (field in c("X", "A", "mu1_true", "mu0_true")) {
    expect_identical(both_arms[[field]], null_data[[field]])
    expect_identical(treated_arm[[field]], null_data[[field]])
  }
  unchanged <- null_data$R != "s1"
  deviated <- !unchanged
  for (arm in c("Y_1", "Y_0")) {
    expect_identical(both_arms[[arm]][unchanged], null_data[[arm]][unchanged])
    expect_identical(treated_arm[[arm]][unchanged], null_data[[arm]][unchanged])
    expect_false(identical(both_arms[[arm]][deviated], null_data[[arm]][deviated]))
  }
  expect_identical(treated_arm$Y_0[deviated], null_data$Y_0[deviated])
  expect_false(identical(treated_arm$Y_1[deviated], null_data$Y_1[deviated]))
  # A shared shift raises both arms' risks at the deviated source.
  expect_gt(mean(both_arms$Y_0[deviated]), mean(null_data$Y_0[deviated]))
  expect_gt(mean(both_arms$Y_1[deviated]), mean(null_data$Y_1[deviated]))

  expect_error(
    generate_simulation_data(
      n_total = 60, K = 1, p = 3, config = "C1", dgp_type = "roce",
      estimand_type = "sample", deviation_mechanism = "both_arms"
    ),
    "applies only to dgp_type = 'face'"
  )
  expect_error(generate(ate_deviation = 1, deviation_mechanism = "shared"), "arg")
  # Forwarded defaults must not trip the roce guard (drivers forward them).
  roce_data <- generate_simulation_data(
    n_total = 60, K = 1, p = 3, config = "C1", dgp_type = "roce",
    estimand_type = "sample", deviation_mechanism = "treated_arm",
    misspecification_strength = FACE_MISSPECIFICATION_STRENGTH
  )
  expect_identical(roce_data$dgp_type, "roce")
})
