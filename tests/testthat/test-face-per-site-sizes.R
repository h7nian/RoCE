# test-face-per-site-sizes.R - Per-site sample sizes for the FACE paper DGP
#
# Covers resolve_face_site_sizes() (the allocation resolver) and its use by
# generate_face_data() / generate_simulation_data(): the default equal split
# must reproduce the previous behaviour, while explicit per-site sizes
# (n_target + n_source_sizes) must be honoured exactly.

test_that("resolve_face_site_sizes equal-split matches the previous behaviour", {
  # Divisible: 1200 / (2 + 1) = 400 per site, target gets the remainder (0).
  res <- resolve_face_site_sizes(n_total = 1200L, K = 2L)
  expect_equal(res$K, 2L)
  expect_equal(res$n_total, 1200L)
  expect_equal(res$n_target, 400L)
  expect_equal(res$n_source_sizes, c(400L, 400L))

  # Non-divisible: floor(1000 / 3) = 333 per source, remainder (334) to target.
  res2 <- resolve_face_site_sizes(n_total = 1000L, K = 2L)
  expect_equal(res2$n_source_sizes, c(333L, 333L))
  expect_equal(res2$n_target, 334L)
  expect_equal(res2$n_target + sum(res2$n_source_sizes), 1000L)
})

test_that("resolve_face_site_sizes honours explicit per-site sizes", {
  res <- resolve_face_site_sizes(n_target = 500L, n_source_sizes = c(100L, 200L, 300L))
  expect_equal(res$K, 3L)
  expect_equal(res$n_target, 500L)
  expect_equal(res$n_source_sizes, c(100L, 200L, 300L))
  expect_equal(res$n_total, 1100L)
})

test_that("resolve_face_site_sizes validates its inputs", {
  # Per-site mode requires BOTH n_target and n_source_sizes.
  expect_error(resolve_face_site_sizes(n_source_sizes = c(100L, 200L)), "BOTH")
  expect_error(resolve_face_site_sizes(n_target = 100L), "BOTH")
  # Non-positive source sizes are rejected.
  expect_error(resolve_face_site_sizes(n_target = 100L, n_source_sizes = c(100L, -5L)),
               "positive integers")
  # Equal-split mode needs n_total large enough for K + 1 sites.
  expect_error(resolve_face_site_sizes(n_total = 2L, K = 5L), "too small")
  # Equal-split mode needs both n_total and K.
  expect_error(resolve_face_site_sizes(), "requires both")
})

test_that("generate_face_data honours explicit per-site sizes", {
  set.seed(42)
  data <- generate_face_data(config = "C1", p = 10, estimand_type = "sample",
                             n_target = 300L, n_source_sizes = c(100L, 200L))
  expect_equal(data$K, 2L)
  expect_equal(data$n, 600L)
  expect_equal(nrow(data$X), 600L)
  expect_equal(sum(data$R == "t"),  300L)
  expect_equal(sum(data$R == "s1"), 100L)
  expect_equal(sum(data$R == "s2"), 200L)
})

test_that("generate_face_data default equal split is unchanged (backward compatible)", {
  set.seed(42)
  data <- generate_face_data(n_total = 900L, K = 2L, p = 10, estimand_type = "sample")
  expect_equal(data$n, 900L)
  expect_equal(sum(data$R == "t"),  300L)
  expect_equal(sum(data$R == "s1"), 300L)
  expect_equal(sum(data$R == "s2"), 300L)
})

test_that("generate_simulation_data routes per-site sizes to the FACE DGP", {
  set.seed(42)
  data <- generate_simulation_data(
    dgp_type = "face", config = "C1", p = 10, outcome_type = "binary",
    estimand_type = "sample", n_target = 250L, n_source_sizes = c(150L, 400L),
    warn_ignored = FALSE
  )
  expect_equal(data$dgp_type, "face")
  expect_equal(sum(data$R == "t"),  250L)
  expect_equal(sum(data$R == "s1"), 150L)
  expect_equal(sum(data$R == "s2"), 400L)
})

test_that("generate_simulation_data rejects per-site sizes for the roce DGP", {
  expect_error(
    generate_simulation_data(n_total = 900L, K = 2L, p = 10, dgp_type = "roce",
                             n_source_sizes = c(300L, 300L), warn_ignored = FALSE),
    "only.*face"
  )
})
