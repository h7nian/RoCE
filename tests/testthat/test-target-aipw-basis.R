test_that("target-only AIPW uses Z_site for propensity and W_outcome for outcome", {
  set.seed(74201)
  n <- 500L
  Z <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
  W <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
  colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  colnames(W) <- paste0("W", seq_len(ncol(W)))

  ps_true <- plogis(2.8 * Z[, 1L] - 1.8 * Z[, 2L])
  A <- rbinom(n, 1L, ps_true)
  y_prob <- plogis(-0.3 + 1.4 * W[, 1L] - 1.1 * W[, 2L] + 0.4 * A)
  Y <- rbinom(n, 1L, y_prob)

  target_data <- list(
    n = n,
    W_outcome = W,
    Z_site = Z,
    A = A,
    Y = Y
  )

  res <- estimate_target_only_crossfit(
    target_data,
    n_folds = 5L,
    family = "binomial",
    A_val = 1L
  )

  expect_gt(abs(stats::cor(res$prop_scores, Z[, 1L])), 0.45)
  expect_lt(abs(stats::cor(res$prop_scores, W[, 1L])), 0.25)
  expect_gt(abs(stats::cor(res$m_pred, W[, 1L])), 0.35)
})
