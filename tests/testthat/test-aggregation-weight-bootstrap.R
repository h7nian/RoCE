library(testthat)

# Contract tests for the planned adaptive aggregation-weight bootstrap.
#
# The public entry point is intentionally small:
#   estimate_tate_weight_bootstrap(
#     tate_result, B, seed, relearn_weights, screening_rule
#   )
#
# The two internal helpers named below make the statistically important parts
# independently testable.  Tests skip while the API is not yet implemented;
# once the functions exist, a partial implementation must fail closed rather
# than silently falling back to the conditional plug-in variance.

.weight_bootstrap_api_available <- function() {
  exists("estimate_tate_weight_bootstrap", mode = "function", inherits = TRUE)
}

.weight_bootstrap_helpers_available <- function() {
  all(vapply(
    c(
      ".reweight_tate_inner_info",
      ".weight_bootstrap_outer_estimate"
    ),
    exists,
    logical(1L),
    mode = "function",
    inherits = TRUE
  ))
}

.make_weight_bootstrap_component <- function(
    target_id, source_ids, target_value, source_values) {
  K <- length(source_ids)
  varphi <- target_value - mean(target_value)
  zeta <- lapply(source_values, function(value) {
    rep(mean(value), length(target_id)) - mean(value)
  })
  xi <- lapply(source_values, function(value) value - mean(value))
  target_estimate <- mean(target_value)
  source_estimate <- vapply(source_values, mean, numeric(1L))

  list(
    target_idx = as.integer(target_id),
    source_idx = lapply(source_ids, as.integer),
    V_ot = mean(varphi^2),
    V_t = vapply(zeta, function(x) mean(x^2), numeric(1L)),
    V_s = vapply(xi, function(x) mean(x^2), numeric(1L)),
    C_ot = vapply(zeta, function(x) mean(varphi * x), numeric(1L)),
    C_cross = matrix(0, K, K),
    avg_target_est = target_estimate,
    avg_source_est = source_estimate,
    fold_target_estimate = target_estimate,
    fold_source_estimates = source_estimate,
    mu_pred_ts = source_estimate,
    delta_ts = rep(0, K),
    n_t = length(target_id),
    n_s = vapply(source_ids, length, integer(1L)),
    varphi_ot = varphi,
    zeta_components = zeta,
    xi_components = xi
  )
}

.make_weight_bootstrap_toy <- function(hard = FALSE,
                                       heterogeneous_fold_means = FALSE) {
  source_names <- c("s1", "s2")
  target_folds <- list(1:2, 3:4, 5:6)
  source_folds <- list(
    s1 = list(1:2, 3:4, 5:6),
    s2 = list(1:2, 3:4, 5:6)
  )
  if (isTRUE(heterogeneous_fold_means)) {
    target_values <- c(0, 2, 8, 12, -7, -1)
    source_values <- list(
      s1 = c(1, 2, 5, 9, -4, 0),
      s2 = c(10, 12, 18, 22, 3, 9)
    )
  } else {
    target_values <- c(-1, 1, -2, 2, -3, 3)
    source_values <- list(
      s1 = c(-0.5, 0.5, -1, 1, -1.5, 1.5),
      # A large location shift makes s2 a deterministic hard-screening case.
      s2 = c(9, 11, 8, 12, 7, 13)
    )
  }

  fold_info <- lapply(seq_along(target_folds), function(k) {
    component <- .make_weight_bootstrap_component(
      target_folds[[k]],
      lapply(source_folds, `[[`, k),
      target_values[target_folds[[k]]],
      lapply(seq_along(source_names), function(j) {
        source_values[[j]][source_folds[[j]][[k]]]
      })
    )
    component
  })

  # For outer fold k, Phase 1b contains only global observations outside k.
  inner_fold_info <- lapply(seq_along(target_folds), function(k) {
    inner_ids <- setdiff(seq_along(target_folds), k)
    inner <- lapply(inner_ids, function(k2) {
      component <- fold_info[[k2]]
      component$outer_fold <- k
      component$inner_fold <- k2
      component
    })
    names(inner) <- paste0("k2_", inner_ids)
    pooled <- RoCE:::.average_aggregation_components(inner)
    pooled$components <- inner
    pooled
  })

  fold_weights <- matrix(
    rep(if (hard) c(0.25, 0) else c(0.25, 0.10), 3L),
    nrow = 3L, byrow = TRUE,
    dimnames = list(NULL, source_names)
  )
  fold_source_included <- matrix(
    rep(if (hard) c(TRUE, FALSE) else c(TRUE, TRUE), 3L),
    nrow = 3L, byrow = TRUE,
    dimnames = list(NULL, source_names)
  )
  n_t <- length(target_values)
  n_s <- vapply(source_values, length, integer(1L))
  N_all <- n_t + sum(n_s)
  all_phi <- RoCE:::.compute_phase3_all_phi(
    n_folds = 3L,
    fold_weights = fold_weights,
    fold_info = fold_info,
    n_t = n_t,
    n_source_full = n_s,
    N_all = N_all,
    source_sites = source_names,
    K = 2L,
    verbose = FALSE
  )

  list(
    estimate = mean(all_phi),
    variance = RoCE:::.multisite_pseudovalue_variance(
      all_phi, c(n_t, n_s)
    ),
    se = sqrt(RoCE:::.multisite_pseudovalue_variance(
      all_phi, c(n_t, n_s)
    )),
    weights = colMeans(fold_weights),
    fold_weights = fold_weights,
    fold_source_included = fold_source_included,
    fold_lambdas = rep(1, 3L),
    aggregation_lambda_selection = 1,
    aggregation_lambda_rule = "min",
    aggregation_screening_rule = if (hard) "hard_threshold" else "soft_penalty",
    n_folds = 3L,
    N_all = N_all,
    all_phi_tau = all_phi,
    source_sites = source_names,
    intermediates = list(
      fold_info = fold_info,
      inner_fold_info = inner_fold_info,
      fold_weights = fold_weights,
      fold_source_included = fold_source_included,
      sample_sizes = list(n_t = n_t, n_source = n_s, N_all = N_all)
    )
  )
}

test_that("weight-bootstrap public API has the reviewable five-argument contract", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  expect_identical(
    names(formals(estimate_tate_weight_bootstrap))[1:5],
    c("tate_result", "B", "seed", "relearn_weights", "screening_rule")
  )
})

test_that("all-one multipliers reproduce the fitted estimator exactly", {
  skip_if_not(.weight_bootstrap_helpers_available(), "weight-bootstrap helpers not implemented")
  toy <- .make_weight_bootstrap_toy()
  ones <- list(t = rep(1, 6L), s1 = rep(1, 6L), s2 = rep(1, 6L))

  fixed_estimate <- .weight_bootstrap_outer_estimate(
    toy$intermediates$fold_info,
    toy$fold_weights,
    ones$t,
    ones[c("s1", "s2")],
    c("s1", "s2"),
    6L,
    c(s1 = 6L, s2 = 6L)
  )
  expect_equal(fixed_estimate, toy$estimate, tolerance = 1e-14)

  for (k in seq_len(toy$n_folds)) {
    reweighted <- .reweight_tate_inner_info(
      toy$intermediates$inner_fold_info[[k]],
      ones$t,
      ones[c("s1", "s2")],
      c("s1", "s2")
    )
    original <- toy$intermediates$inner_fold_info[[k]]
    for (field in c(
      "V_ot", "V_t", "V_s", "C_ot", "C_cross",
      "avg_target_est", "avg_source_est", "mu_pred_ts", "delta_ts"
    )) {
      expect_equal(reweighted[[field]], original[[field]], tolerance = 1e-14,
                   info = paste("outer fold", k, "field", field))
    }
  }
})

test_that("all-one relearning preserves heterogeneous-fold moments and weights", {
  skip_if_not(.weight_bootstrap_helpers_available(), "weight-bootstrap helpers not implemented")
  toy <- .make_weight_bootstrap_toy(heterogeneous_fold_means = TRUE)
  ones <- list(t = rep(1, 6L), s1 = rep(1, 6L), s2 = rep(1, 6L))
  reweighted <- lapply(seq_len(toy$n_folds), function(k) {
    .reweight_tate_inner_info(
      toy$intermediates$inner_fold_info[[k]], ones$t,
      ones[c("s1", "s2")], c("s1", "s2")
    )
  })

  # Fold means deliberately differ, so an incorrect global recentering of
  # already fold-centered IFs would add a between-fold variance term here.
  expect_gt(stats::var(vapply(
    toy$intermediates$fold_info,
    `[[`, numeric(1L), "fold_target_estimate"
  )), 0)
  for (k in seq_len(toy$n_folds)) {
    original <- toy$intermediates$inner_fold_info[[k]]
    for (field in c(
      "V_ot", "V_t", "V_s", "C_ot", "C_cross",
      "avg_target_est", "avg_source_est", "mu_pred_ts", "delta_ts"
    )) {
      expect_equal(reweighted[[k]][[field]], original[[field]],
                   tolerance = 1e-14,
                   info = paste("outer fold", k, "field", field))
    }
  }

  original_phase2 <- RoCE:::.compute_phase2_weights(
    n_folds = toy$n_folds,
    inner_fold_info = toy$intermediates$inner_fold_info,
    fold_info = toy$intermediates$fold_info,
    lambda_selection = 1,
    K = 2L,
    verbose = FALSE,
    screening_rule = "soft_penalty"
  )
  reweighted_phase2 <- RoCE:::.compute_phase2_weights(
    n_folds = toy$n_folds,
    inner_fold_info = reweighted,
    fold_info = toy$intermediates$fold_info,
    lambda_selection = 1,
    K = 2L,
    verbose = FALSE,
    screening_rule = "soft_penalty"
  )
  expect_equal(reweighted_phase2$fold_weights,
               original_phase2$fold_weights, tolerance = 1e-12)
  expect_equal(reweighted_phase2$fold_wald_statistics,
               original_phase2$fold_wald_statistics, tolerance = 1e-12)

  original_estimate <- .weight_bootstrap_outer_estimate(
    toy$intermediates$fold_info, original_phase2$fold_weights,
    ones$t, ones[c("s1", "s2")], c("s1", "s2"), 6L,
    c(s1 = 6L, s2 = 6L)
  )
  reweighted_estimate <- .weight_bootstrap_outer_estimate(
    toy$intermediates$fold_info, reweighted_phase2$fold_weights,
    ones$t, ones[c("s1", "s2")], c("s1", "s2"), 6L,
    c(s1 = 6L, s2 = 6L)
  )
  expect_equal(reweighted_estimate, original_estimate, tolerance = 1e-14)
})

test_that("multiplier draws use one shared global observation ID map", {
  skip_if_not(.weight_bootstrap_helpers_available(), "weight-bootstrap helpers not implemented")
  toy <- .make_weight_bootstrap_toy()
  ones <- list(t = rep(1, 6L), s1 = rep(1, 6L), s2 = rep(1, 6L))
  perturbed <- ones
  perturbed$t[[1L]] <- 2

  original <- lapply(seq_len(3L), function(k) {
    .reweight_tate_inner_info(
      toy$intermediates$inner_fold_info[[k]], ones$t,
      ones[c("s1", "s2")], c("s1", "s2")
    )
  })
  changed <- lapply(seq_len(3L), function(k) {
    .reweight_tate_inner_info(
      toy$intermediates$inner_fold_info[[k]], perturbed$t,
      perturbed[c("s1", "s2")], c("s1", "s2")
    )
  })

  # Global target ID 1 is the evaluation observation of outer fold 1.  It is
  # absent from fold 1's training set but present in the training sets used to
  # learn weights for folds 2 and 3.  Reusing one multiplier by global ID must
  # therefore change exactly the latter two Phase 1b summaries.
  expect_equal(changed[[1L]]$avg_target_est,
               original[[1L]]$avg_target_est, tolerance = 0)
  expect_false(isTRUE(all.equal(changed[[2L]]$avg_target_est,
                                original[[2L]]$avg_target_est, tolerance = 0)))
  expect_false(isTRUE(all.equal(changed[[3L]]$avg_target_est,
                                original[[3L]]$avg_target_est, tolerance = 0)))
})

test_that("disabling weight relearning gives the fixed-weight result", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy()

  fixed_a <- estimate_tate_weight_bootstrap(
    toy, B = 20L, seed = 901L, relearn_weights = FALSE,
    screening_rule = "soft_penalty", save_draws = TRUE
  )
  fixed_b <- estimate_tate_weight_bootstrap(
    toy, B = 20L, seed = 901L, relearn_weights = FALSE,
    screening_rule = toy$aggregation_screening_rule, save_draws = TRUE
  )
  expect_equal(fixed_a$draws$fixed_weight_estimate,
               fixed_b$draws$fixed_weight_estimate,
               tolerance = 0)
  expect_equal(fixed_a$se_fixed_weight_bootstrap,
               fixed_b$se_fixed_weight_bootstrap, tolerance = 0)
  expect_false(fixed_a$weights_relearned)
  expect_null(fixed_a$draws$weight_relearn_estimate)
  expect_null(fixed_a$draws$fold_weights)
})

test_that("source permutation is equivariant", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy()
  permuted <- toy
  permutation <- c("s2", "s1")
  permuted$source_sites <- permutation
  permuted$weights <- toy$weights[permutation]
  permuted$fold_weights <- toy$fold_weights[, permutation, drop = FALSE]
  permuted$fold_source_included <-
    toy$fold_source_included[, permutation, drop = FALSE]
  permuted$intermediates$fold_weights <- permuted$fold_weights
  permuted$intermediates$fold_source_included <- permuted$fold_source_included
  permuted$intermediates$sample_sizes$n_source <-
    toy$intermediates$sample_sizes$n_source[permutation]
  permuted$intermediates$fold_info <- lapply(
    toy$intermediates$fold_info,
    function(x) {
      permutation_index <- match(permutation, c("s1", "s2"))
      x$source_idx <- x$source_idx[permutation]
      x$fold_source_estimates <- x$fold_source_estimates[permutation_index]
      x$avg_source_est <- x$avg_source_est[permutation_index]
      x$mu_pred_ts <- x$mu_pred_ts[permutation_index]
      x$delta_ts <- x$delta_ts[permutation_index]
      x$V_t <- x$V_t[permutation_index]
      x$V_s <- x$V_s[permutation_index]
      x$C_ot <- x$C_ot[permutation_index]
      x$C_cross <- x$C_cross[
        permutation_index, permutation_index, drop = FALSE
      ]
      x$zeta_components <- x$zeta_components[permutation_index]
      x$xi_components <- x$xi_components[permutation_index]
      x
    }
  )
  # Rebuild Phase 1b from the permuted outer raw components so no stale source
  # ordering can accidentally make the test pass.
  permuted$intermediates$inner_fold_info <- lapply(1:3, function(k) {
    inner_ids <- setdiff(1:3, k)
    inner <- lapply(inner_ids, function(k2) {
      component <- permuted$intermediates$fold_info[[k2]]
      component$outer_fold <- k
      component$inner_fold <- k2
      component
    })
    names(inner) <- paste0("k2_", inner_ids)
    pooled <- RoCE:::.average_aggregation_components(inner)
    pooled$components <- inner
    pooled
  })

  original_result <- estimate_tate_weight_bootstrap(
    toy, B = 20L, seed = 123L, relearn_weights = TRUE,
    screening_rule = "soft_penalty", save_draws = TRUE
  )
  permuted_result <- estimate_tate_weight_bootstrap(
    permuted, B = 20L, seed = 123L, relearn_weights = TRUE,
    screening_rule = "soft_penalty", save_draws = TRUE
  )
  # The coordinate-descent optimizer can take a different, but numerically
  # equivalent, path after coordinates are permuted.  The resulting differences
  # are below 1e-7; use a tolerance above optimizer arithmetic noise while still
  # detecting any material source-order dependence.
  expect_equal(permuted_result$draws$weight_relearn_estimate,
               original_result$draws$weight_relearn_estimate,
               tolerance = 1e-6)
  expect_equal(permuted_result$se_weight_relearn_bootstrap,
               original_result$se_weight_relearn_bootstrap,
               tolerance = 1e-6)
})

test_that("hard-screened sources have exactly zero bootstrap weights", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy(hard = TRUE)
  result <- estimate_tate_weight_bootstrap(
    toy, B = 20L, seed = 88L, relearn_weights = TRUE,
    screening_rule = "hard_threshold", save_draws = TRUE
  )

  expect_true(all(
    result$draws$fold_weights[!result$draws$fold_source_included] == 0
  ))
})

test_that("weight bootstrap is reproducible and RNG-neutral", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy()
  set.seed(314L)
  expected_next <- runif(1L)

  set.seed(314L)
  first <- estimate_tate_weight_bootstrap(
    toy, B = 20L, seed = 29L, relearn_weights = TRUE,
    screening_rule = "soft_penalty", save_draws = TRUE
  )
  actual_next <- runif(1L)
  second <- estimate_tate_weight_bootstrap(
    toy, B = 20L, seed = 29L, relearn_weights = TRUE,
    screening_rule = "soft_penalty", save_draws = TRUE
  )

  expect_equal(actual_next, expected_next, tolerance = 0)
  expect_identical(first$draws$weight_relearn_estimate,
                   second$draws$weight_relearn_estimate)
  expect_identical(first$draws$fold_weights,
                   second$draws$fold_weights)
})

test_that("weight bootstrap fails closed on malformed or incomplete inputs", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy()

  for (bad_B in list(1L, 2.5, NA_real_, Inf)) {
    expect_error(
      estimate_tate_weight_bootstrap(
        toy, B = bad_B, seed = 1L, relearn_weights = TRUE,
        screening_rule = "soft_penalty"
      ),
      "B|replicate"
    )
  }
  expect_error(
    estimate_tate_weight_bootstrap(
      toy, B = 5L, seed = NA_integer_, relearn_weights = TRUE,
      screening_rule = "soft_penalty"
    ),
    "seed"
  )
  expect_error(
    estimate_tate_weight_bootstrap(
      toy, B = 5L, seed = 1L, relearn_weights = NA,
      screening_rule = "soft_penalty"
    ),
    "relearn_weights|flags"
  )
  expect_error(
    estimate_tate_weight_bootstrap(
      toy, B = 5L, seed = 1L, relearn_weights = TRUE,
      screening_rule = "unknown"
    ),
    "screening_rule|soft_penalty|hard_threshold"
  )

  missing_ids <- toy
  missing_ids$intermediates$fold_info[[1L]]$target_idx <- NULL
  expect_error(
    estimate_tate_weight_bootstrap(
      missing_ids, B = 5L, seed = 1L, relearn_weights = TRUE,
      screening_rule = "soft_penalty"
    ),
    "target_idx|observation|ID|fold metadata"
  )

  duplicate_id <- toy
  duplicate_id$intermediates$fold_info[[2L]]$target_idx[1L] <- 1L
  expect_error(
    estimate_tate_weight_bootstrap(
      duplicate_id, B = 5L, seed = 1L, relearn_weights = TRUE,
      screening_rule = "soft_penalty"
    ),
    "exactly once|each observation once|duplicate|partition"
  )

  nonfinite <- toy
  nonfinite$intermediates$inner_fold_info[[1L]]$components[[1L]]$
    varphi_ot[1L] <- NA_real_
  expect_error(
    estimate_tate_weight_bootstrap(
      nonfinite, B = 5L, seed = 1L, relearn_weights = TRUE,
      screening_rule = "soft_penalty"
    ),
    "finite|NA"
  )
})

test_that("weight bootstrap validates the complete outer-inner ID topology", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy()
  call_bootstrap <- function(object) {
    estimate_tate_weight_bootstrap(
      object, B = 2L, seed = 1L, relearn_weights = TRUE,
      screening_rule = "soft_penalty"
    )
  }

  # Regression for the original bug: the first component of outer fold 1 is
  # inner fold 2, so IDs 1:2 are valid globally but belong to the wrong fold.
  wrong_inner_target <- toy
  wrong_inner_target$intermediates$inner_fold_info[[1L]]$
    components[[1L]]$target_idx <- 1:2
  expect_error(call_bootstrap(wrong_inner_target), "exactly match.*target")

  wrong_inner_source <- toy
  wrong_inner_source$intermediates$inner_fold_info[[1L]]$
    components[[1L]]$source_idx$s1 <- rev(3:4)
  expect_error(call_bootstrap(wrong_inner_source), "exactly match.*source")

  duplicate_inner_fold <- toy
  duplicate_inner_fold$intermediates$inner_fold_info[[1L]]$
    components[[2L]]$inner_fold <- 2L
  expect_error(call_bootstrap(duplicate_inner_fold),
               "cover every fold|exactly once|exactly match|fold labels")

  fractional_inner_fold <- toy
  fractional_inner_fold$intermediates$inner_fold_info[[1L]]$
    components[[1L]]$inner_fold <- 2.5
  expect_error(call_bootstrap(fractional_inner_fold), "inner_fold.*integer")

  duplicate_outer_id <- toy
  duplicate_outer_id$intermediates$fold_info[[1L]]$target_idx <- c(1L, 1L)
  expect_error(call_bootstrap(duplicate_outer_id), "unique integer observation IDs")

  out_of_range_source_id <- toy
  out_of_range_source_id$intermediates$fold_info[[1L]]$source_idx$s1 <-
    c(1L, 7L)
  expect_error(call_bootstrap(out_of_range_source_id), "source_idx.*s1")

  wrong_varphi_length <- toy
  wrong_varphi_length$intermediates$fold_info[[1L]]$varphi_ot <- 0
  expect_error(call_bootstrap(wrong_varphi_length), "varphi_ot.*exactly 2")

  wrong_zeta_length <- toy
  wrong_zeta_length$intermediates$inner_fold_info[[1L]]$
    components[[1L]]$zeta_components[[1L]] <- 0
  expect_error(call_bootstrap(wrong_zeta_length), "zeta.*exactly 2")

  wrong_xi_length <- toy
  wrong_xi_length$intermediates$inner_fold_info[[1L]]$
    components[[1L]]$xi_components[[1L]] <- 0
  expect_error(call_bootstrap(wrong_xi_length), "xi.*exactly 2")
})

test_that("weight bootstrap validates sizes, source columns, lambda mode, and capacity", {
  skip_if_not(.weight_bootstrap_api_available(), "weight-bootstrap API not implemented")
  toy <- .make_weight_bootstrap_toy()
  call_bootstrap <- function(object, B = 2L) {
    estimate_tate_weight_bootstrap(
      object, B = B, seed = 1L, relearn_weights = TRUE,
      screening_rule = "soft_penalty"
    )
  }

  fractional_size <- toy
  fractional_size$intermediates$sample_sizes$n_t <- 6.5
  expect_error(call_bootstrap(fractional_size), "sample_sizes.*integer")

  overflow_size <- toy
  overflow_size$intermediates$sample_sizes$n_source[[1L]] <-
    .Machine$integer.max + 1
  expect_error(call_bootstrap(overflow_size), "sample_sizes.*integer")

  wrong_size_names <- toy
  names(wrong_size_names$intermediates$sample_sizes$n_source) <- c("s2", "s1")
  expect_error(call_bootstrap(wrong_size_names), "names and order")

  wrong_weight_columns <- toy
  colnames(wrong_weight_columns$fold_weights) <- c("s2", "s1")
  expect_error(call_bootstrap(wrong_weight_columns), "columns exactly match")

  missing_weight_columns <- toy
  colnames(missing_weight_columns$fold_weights) <- NULL
  expect_error(call_bootstrap(missing_weight_columns), "columns exactly match")

  cv_fit <- toy
  cv_fit$aggregation_lambda_selection <- "cv"
  expect_error(call_bootstrap(cv_fit), "fixed finite numeric|weighted inner CV")

  expect_error(call_bootstrap(toy, B = .Machine$integer.max + 1),
               "integer.max")
  capacity_B <- floor(.Machine$integer.max / (toy$n_folds * 2L)) + 1
  expect_error(call_bootstrap(toy, B = capacity_B),
               "B.*n_folds.*K|array dimension")
})
