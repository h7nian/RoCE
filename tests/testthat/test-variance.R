# test-variance.R - Tests for variance utilities

library(testthat)

# Package loaded by helper-load.R (all functions available via FACEC namespace)

test_that("GLM family support is restricted to gaussian/binomial", {
	expect_equal(VALID_GLM_FAMILIES, c("gaussian", "binomial"))
	expect_error(resolve_glm_family("poisson"), "Unsupported GLM family")
	expect_error(resolve_glm_family("gamma"), "Unsupported GLM family")
})

test_that("calculate_aggregated_variance matches target-only baseline at eta=0", {
	skip_if_not(exists("calculate_aggregated_variance_cpp"), message = "C++ not compiled")

	eta <- c(0, 0)
	variances <- list(
		V_t = c(0.10, 0.12),
		V_s = c(0.05, 0.04),
		V_ot = 0.20
	)
	C_ot <- c(0.01, 0.02)
	n_samples <- list(n_t = 100, n_s = c(200, 250))

	v <- calculate_aggregated_variance(eta, variances, C_ot, n_samples)
	expect_equal(v, variances$V_ot / n_samples$n_t, tolerance = 1e-10)
})

test_that("calculate_aggregated_variance is monotone in V_ot when eta is fixed", {
	skip_if_not(exists("calculate_aggregated_variance_cpp"), message = "C++ not compiled")

	eta <- c(0.3, 0.2)
	variances_low <- list(V_t = c(0.10, 0.12), V_s = c(0.05, 0.04), V_ot = 0.10)
	variances_high <- list(V_t = c(0.10, 0.12), V_s = c(0.05, 0.04), V_ot = 0.50)
	C_ot <- c(0.01, 0.02)
	n_samples <- list(n_t = 100, n_s = c(200, 250))

	v_low <- calculate_aggregated_variance(eta, variances_low, C_ot, n_samples)
	v_high <- calculate_aggregated_variance(eta, variances_high, C_ot, n_samples)

	expect_true(is.finite(v_low) && is.finite(v_high))
	expect_gte(v_high, v_low)
})

test_that("public defaults keep lambda cache on", {
	expect_true(isTRUE(formals(run_crossfit)$use_lambda_cache))
	expect_true(isTRUE(formals(run_single_simulation)$use_lambda_cache))
	expect_true(isTRUE(formals(run_simulation_study)$use_lambda_cache))
})
