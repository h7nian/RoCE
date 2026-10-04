test_that("arm variance diagnostics match the TATE and empirical-mass derivatives", {
  skip_on_cran()
  source(file.path(repo_root, "diagnosis", "repair", "layer_comparison_summary.R"), local = TRUE)
  set.seed(4821)
  data <- split_data_by_site(generate_bounded_data(n_target = 500, n_source_sizes = 500,
    K = 1, p = 4, config = "C1"))
  sizes <- vapply(data, function(site) as.integer(site$n), integer(1L))
  for (layers in c(2L, 3L)) {
    fit <- run_tate_crossfit(data, n_folds = 4L, communication_mode = "one_round",
      target_nuisance_method = "hou_calibrated", crossfit_layers = layers, nlambda_init = 4L,
      nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, calibration_layout = "compact",
      calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "calibrated", target_radius = 12),
      M_tau = 12, M_tau_inference = 12, lambda_selection = .5, aggregation_mode = "joint_tate", verbose = FALSE)
    for (mode in c("common_tate", "separate_arms", "joint_tate")) {
      selected <- if (mode == "joint_tate") fit else reaggregate_tate_crossfit(data, fit,
        M_tau_inference = fit$M_tau_inference, aggregation_mode = mode, verbose = FALSE)
      summary <- summarize_layer_arms(selected, data)
      expect_equal(summary$estimates["mu1"] - summary$estimates["mu0"],
                   unname(selected$estimate), tolerance = 1e-12, ignore_attr = TRUE)
      expect_gte(min(eigen(summary$covariance, symmetric = TRUE, only.values = TRUE)$values), -1e-12)
    }
    summary <- summarize_layer_arms(fit, data, keep_influence = TRUE)
    records <- lapply(seq_len(fit$n_folds), function(fold)
      .joint_inner_records(fit$arm_results$mu1$intermediates$inner_fold_info[[fold]],
                           fit$arm_results$mu0$intermediates$inner_fold_info[[fold]]))
    # Relearn eta from perturbed empirical masses, then recompute each arm's
    # weighted estimating equation. This does not use the analytic derivatives.
    point <- function(mass) {
      numerators <- matrix(0, length(sizes), 2L)
      for (fold in seq_len(fit$n_folds)) {
        moments <- .joint_inner_moments(records[[fold]], mass)
        eta <- .solve_joint_weights(moments, .5, "soft_penalty")$weights
        for (arm_index in 1:2) {
          info <- fit$arm_results[[arm_index]]$intermediates$fold_info[[fold]]
          weight <- eta[arm_index]
          target <- (1 - weight) * (info$varphi_ot + info$fold_target_estimate) +
            weight * (info$zeta_components[[1L]] + info$mu_pred_ts[1L])
          source <- weight * (info$xi_components[[1L]] + info$delta_ts[1L])
          numerators[1L, arm_index] <- numerators[1L, arm_index] + sum(target * mass[[1L]][info$target_idx])
          numerators[2L, arm_index] <- numerators[2L, arm_index] + sum(source * mass[[2L]][info$source_idx[[1L]]])
        }
      }
      colSums(numerators / vapply(mass, sum, numeric(1L)))
    }
    mass <- lapply(sizes, function(size) rep(1, size))
    expect_equal(point(mass), as.numeric(summary$estimates), tolerance = 1e-12)
    for (site in seq_along(sizes)) {
      step <- 1e-4
      upper <- lower <- mass
      upper[[site]][1L] <- 1 + step
      lower[[site]][1L] <- 1 - step
      difference <- (point(upper) - point(lower)) / (2 * step)
      index <- if (site == 1L) 1L else sizes[1L] + 1L
      expect_equal(difference, as.numeric(summary$influence[index, ]), tolerance = 1e-7)
    }
  }
})
