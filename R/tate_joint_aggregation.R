# Joint two-arm aggregation. Columns are treated-arm increments followed by
# negative control-arm increments. A source remains one site with a full
# covariance matrix, even though it contributes to two weight coordinates.

.joint_inner_records <- function(inner1, inner0) {
  arm1 <- .inner_fold_records(inner1)
  arm0 <- .inner_fold_records(inner0)
  if (!identical(arm1$target_ids, arm0$target_ids) ||
      !identical(arm1$source_ids, arm0$source_ids)) {
    stop(".joint_inner_records: the two arms must use identical observation order.", call. = FALSE)
  }
  K <- length(arm1$source)
  target <- Map(function(x1, x0) {
    cbind(x1[, 1L] - x0[, 1L],
      sweep(x1[, -1L, drop = FALSE], 1L, x1[, 1L], "-"),
      -sweep(x0[, -1L, drop = FALSE], 1L, x0[, 1L], "-"))
  }, arm1$target, arm0$target)
  source <- lapply(seq_len(K), function(j) {
    Map(function(x1, x0) {
      x <- matrix(0, length(x1), 1L + 2L * K)
      x[, 1L + j] <- x1
      x[, 1L + K + j] <- -x0
      x
    }, arm1$source[[j]], arm0$source[[j]])
  })
  list(blocks = c(list(target), source),
       ids = c(list(arm1$target_ids), arm1$source_ids))
}

.joint_inner_moments <- function(records, mass = NULL) {
  if (is.null(mass)) mass <- lapply(records$ids, function(ids) rep(1, max(unlist(ids))))
  sites <- lapply(seq_along(records$blocks), function(site) {
    blocks <- records$blocks[[site]]
    weights <- lapply(records$ids[[site]], function(ids) mass[[site]][ids])
    if (any(!is.finite(unlist(weights))) || any(unlist(weights) <= 0)) {
      stop(".joint_inner_moments: observation masses must be finite and positive.", call. = FALSE)
    }
    centered <- do.call(rbind, Map(function(x, w) {
      sweep(x, 2L, colSums(x * w) / sum(w), "-")
    }, blocks, weights))
    raw <- do.call(rbind, blocks)
    weight <- unlist(weights, use.names = FALSE)
    list(raw = raw, centered = centered, n = nrow(raw),
      mean = colSums(raw * weight) / sum(weight),
      covariance = crossprod(centered, centered * weight) / sum(weight))
  })
  covariance <- Reduce(`+`, lapply(sites, function(site) site$covariance / site$n))
  mean <- Reduce(`+`, lapply(sites, `[[`, "mean"))
  list(Q = covariance[-1L, -1L, drop = FALSE], l = covariance[-1L, 1L],
    constant = covariance[1L, 1L], discrepancy = mean[-1L],
    N_all = sum(vapply(sites, `[[`, integer(1L), "n")), n_t = sites[[1L]]$n,
    sites = sites)
}

.solve_joint_weights <- function(moments, lambda, screening_rule,
                                 tol = 1e-10, max_iter = 100000L) {
  lambda <- .validate_lambda_scalar(lambda, ".solve_joint_weights")
  if (lambda > LAMBDA_MAX) stop("Joint aggregation lambda exceeds LAMBDA_MAX.", call. = FALSE)
  screening_rule <- match.arg(screening_rule, c("soft_penalty", "hard_threshold", "quadratic_bias"))
  Q <- moments$Q
  q <- ncol(Q)
  wald <- abs(moments$discrepancy) / sqrt(pmax(diag(Q), VARIANCE_MIN))
  penalty <- pmax(lambda * wald - 1, 0)
  included <- rep(TRUE, q)
  if (screening_rule == "hard_threshold") {
    if (lambda <= 0) stop("Hard thresholding requires a positive lambda.", call. = FALSE)
    included <- wald <= 1 / lambda
    penalty[] <- 0
  } else if (screening_rule == "quadratic_bias") {
    diag(Q) <- diag(Q) + moments$n_t^(AGG_QUADRATIC_BIAS_POWER - 1) * moments$discrepancy^2
    penalty[] <- 0
  }
  eta <- numeric(q)
  active <- which(included)
  iterations <- 1L
  residual <- 0
  if (length(active)) {
    H <- 2 * Q[active, active, drop = FALSE]
    linear <- 2 * moments$l[active]
    lasso <- penalty[active] / moments$N_all
    eigenvalues <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(H)) || min(eigenvalues) <= 0) {
      stop(".solve_joint_weights: selected weight coordinates have no positive definite curvature.",
           call. = FALSE)
    }
    beta <- numeric(length(active))
    if (all(lasso == 0)) {
      beta <- drop(solve(H, -linear))
    } else {
      for (iterations in seq_len(max_iter)) {
        for (j in seq_along(beta)) {
          partial <- -linear[j] - sum(H[j, ] * beta) + H[j, j] * beta[j]
          beta[j] <- sign(partial) * max(abs(partial) - lasso[j], 0) / H[j, j]
        }
        gradient <- drop(H %*% beta) + linear
        kkt <- ifelse(beta != 0, gradient + lasso * sign(beta),
                      pmax(abs(gradient) - lasso, 0))
        residual <- max(abs(kkt))
        if (residual <= tol * max(1, max(abs(linear)), max(abs(H)))) break
      }
      if (residual > tol * max(1, max(abs(linear)), max(abs(H)))) {
        stop(".solve_joint_weights: coordinate descent did not satisfy KKT tolerance.", call. = FALSE)
      }
    }
    gradient <- drop(H %*% beta) + linear
    residual <- max(abs(ifelse(beta != 0, gradient + lasso * sign(beta),
                               pmax(abs(gradient) - lasso, 0))))
    eta[active] <- beta
  }
  list(weights = eta, wald = wald, penalty = penalty, included = included,
    iterations = iterations, kkt_residual = residual, lambda = lambda,
    objective = moments$N_all * (moments$constant + 2 * sum(moments$l * eta) +
      drop(crossprod(eta, Q %*% eta))) + sum(penalty * abs(eta)))
}

.select_joint_lambda <- function(records, lambda_selection, lambda_grid,
                                 lambda_rule, screening_rule) {
  if (!identical(lambda_selection, "cv")) {
    return(.validate_lambda_scalar(lambda_selection, ".select_joint_lambda"))
  }
  grid <- .aggregation_lambda_grid(NULL, lambda_grid)
  folds <- seq_along(records$blocks[[1L]])
  subset_records <- function(indices) list(
    blocks = lapply(records$blocks, `[`, indices), ids = lapply(records$ids, `[`, indices))
  scores <- vapply(grid, function(lambda) {
    vapply(folds, function(fold) {
      train <- .joint_inner_moments(subset_records(setdiff(folds, fold)))
      validation <- .joint_inner_moments(subset_records(fold))
      eta <- .solve_joint_weights(train, lambda, screening_rule)$weights
      validation$N_all * (validation$constant + 2 * sum(validation$l * eta) +
        drop(crossprod(eta, validation$Q %*% eta)))
    }, numeric(1L))
  }, numeric(length(folds)))
  means <- colMeans(scores)
  standard_errors <- apply(scores, 2L, stats::sd) / sqrt(length(folds))
  best <- which.min(means)
  within <- which(means <= means[best] + standard_errors[best])
  one_se <- within[which.max(grid[within])]
  selected <- grid[if (lambda_rule == "1se") one_se else best]
  attr(selected, "lambda_min") <- grid[best]
  attr(selected, "lambda_1se") <- grid[one_se]
  attr(selected, "cv_scores") <- means
  attr(selected, "cv_se") <- standard_errors
  selected
}

# Exact empirical-mass derivative within a fixed active set and selected
# lambda. This includes the off-diagonal covariance of both arms at each
# source. It does not assert validity at selection kinks or remove nuisance
# estimation remainder terms.
.joint_weight_derivative <- function(moments, fit, screening_rule, increment) {
  eta <- fit$weights
  active <- which(abs(eta) > WEIGHT_LAYER_ACTIVE_TOLERANCE)
  gradients <- lapply(moments$sites, function(site) numeric(site$n))
  quadratic <- screening_rule == "quadratic_bias"
  penalized <- screening_rule == "soft_penalty"
  kinks <- if (quadratic) 0L else sum(abs(fit$lambda * fit$wald - 1) < WEIGHT_LAYER_KINK_TOLERANCE)
  if (!length(active)) return(list(sites = gradients, kinks = kinks))
  rate <- if (quadratic) moments$n_t^(AGG_QUADRATIC_BIAS_POWER - 1) else 0
  curvature <- moments$Q
  diag(curvature) <- diag(curvature) + rate * moments$discrepancy^2
  h <- drop(solve(curvature[active, active, drop = FALSE], increment[active]))
  s <- sqrt(pmax(diag(moments$Q), VARIANCE_MIN))
  for (site_index in seq_along(moments$sites)) {
    site <- moments$sites[[site_index]]
    d <- site$centered[, -1L, drop = FALSE]
    b <- site$centered[, 1L]
    covariance <- site$covariance[-1L, -1L, drop = FALSE]
    dQ_eta <- sweep(d * drop(d %*% eta), 2L, drop(covariance %*% eta), "-") / site$n^2
    dl <- sweep(d * b, 2L, site$covariance[-1L, 1L], "-") / site$n^2
    dm <- sweep(site$raw[, -1L, drop = FALSE], 2L, site$mean[-1L], "-") / site$n
    penalty_derivative <- matrix(0, site$n, length(eta))
    if (quadratic) {
      penalty_derivative <- sweep(dm, 2L, 2 * rate * moments$discrepancy * eta, "*")
    } else if (penalized) {
      dQ_diag <- sweep(d^2, 2L, diag(covariance), "-") / site$n^2
      dQ_diag[, diag(moments$Q) < VARIANCE_MIN] <- 0
      dt <- sweep(dm, 2L, sign(moments$discrepancy) / s, "*") -
        sweep(dQ_diag, 2L, abs(moments$discrepancy) / (2 * s^3), "*")
      penalty_derivative <- sweep(dt, 2L,
        sign(eta) * fit$lambda * (fit$lambda * fit$wald > 1) / (2 * moments$N_all), "*")
    }
    gradients[[site_index]] <- -drop((dQ_eta + dl + penalty_derivative)[, active, drop = FALSE] %*% h)
  }
  list(sites = gradients, kinks = kinks)
}

.aggregate_tate_arm_weights <- function(
    data_split, mu1_result, mu0_result, fold_info, inner_fold_info,
    aggregation_mode, lambda_selection, lambda_rule,
    aggregation_lambda_grid, screening_rule, verbose) {
  arms <- list(mu1 = mu1_result, mu0 = mu0_result)
  source_sites <- setdiff(names(data_split), "t")
  K <- length(source_sites)
  n_folds <- length(fold_info)
  n_t <- data_split$t$n
  n_source <- vapply(source_sites, function(site) data_split[[site]]$n, numeric(1L))
  sizes <- c(n_t, n_source)
  N_all <- sum(sizes)
  phase2 <- vector("list", 2L)
  names(phase2) <- names(arms)
  joint_fits <- joint_records <- joint_moments <- vector("list", n_folds)
  if (aggregation_mode == "separate_arms") {
    for (arm in names(arms)) {
      phase2[[arm]] <- .compute_phase2_weights(
        n_folds, arms[[arm]]$intermediates$inner_fold_info,
        arms[[arm]]$intermediates$fold_info, lambda_selection, K, verbose,
        lambda_rule, aggregation_lambda_grid, screening_rule
      )
    }
  } else {
    weights <- matrix(0, n_folds, 2L * K)
    for (fold in seq_len(n_folds)) {
      records <- .joint_inner_records(
        mu1_result$intermediates$inner_fold_info[[fold]],
        mu0_result$intermediates$inner_fold_info[[fold]])
      moments <- .joint_inner_moments(records)
      lambda <- .select_joint_lambda(records, lambda_selection,
        aggregation_lambda_grid, lambda_rule, screening_rule)
      fit <- .solve_joint_weights(moments, lambda, screening_rule)
      weights[fold, ] <- fit$weights
      joint_fits[[fold]] <- fit
      joint_records[[fold]] <- records
      joint_moments[[fold]] <- moments
    }
    for (arm in names(arms)) {
      columns <- if (arm == "mu1") seq_len(K) else K + seq_len(K)
      phase2[[arm]] <- list(
        fold_weights = weights[, columns, drop = FALSE],
        fold_lambdas = vapply(joint_fits, `[[`, numeric(1L), "lambda"),
        fold_wald_statistics = do.call(rbind, lapply(joint_fits, `[[`, "wald"))[, columns, drop = FALSE],
        fold_penalty_coefficients = do.call(rbind, lapply(joint_fits, `[[`, "penalty"))[, columns, drop = FALSE],
        fold_source_included = do.call(rbind, lapply(joint_fits, `[[`, "included"))[, columns, drop = FALSE],
        fold_weight_optimizer_iterations = vapply(joint_fits, `[[`, integer(1L), "iterations"),
        fold_weight_psd_ridge = rep(0, n_folds)
      )
    }
  }
  all_phi <- lapply(names(arms), function(arm) {
    .compute_phase3_all_phi(n_folds, phase2[[arm]]$fold_weights,
      arms[[arm]]$intermediates$fold_info, n_t, n_source, N_all,
      source_sites, K, verbose)
  })
  all_phi_tau <- all_phi[[1L]] - all_phi[[2L]]
  if (aggregation_mode == "separate_arms") {
    gradients <- lapply(names(arms), function(arm) {
      .weight_layer_gradient(arms[[arm]]$intermediates$fold_info,
        arms[[arm]]$intermediates$inner_fold_info, phase2[[arm]]$fold_weights,
        phase2[[arm]]$fold_lambdas, screening_rule,
        phase2[[arm]]$fold_weight_psd_ridge, n_t, n_source)
    })
    gradient <- list(target = gradients[[1L]]$target - gradients[[2L]]$target,
      source = Map(`-`, gradients[[1L]]$source, gradients[[2L]]$source),
      kink_cells = gradients[[1L]]$kink_cells + gradients[[2L]]$kink_cells,
      fold_source_cells = 2L * n_folds * K)
  } else {
    increments <- cbind(
      .outer_fold_increments(mu1_result$intermediates$fold_info, n_t, n_source),
      -.outer_fold_increments(mu0_result$intermediates$fold_info, n_t, n_source))
    gradient_sites <- lapply(sizes, function(n) numeric(n))
    kink_cells <- 0L
    for (fold in seq_len(n_folds)) {
      derivative <- .joint_weight_derivative(joint_moments[[fold]], joint_fits[[fold]],
        screening_rule, increments[fold, ])
      kink_cells <- kink_cells + derivative$kinks
      outer_ids <- c(list(fold_info[[fold]]$target_idx), fold_info[[fold]]$source_idx)
      for (site in seq_along(sizes)) {
        ids <- unlist(joint_records[[fold]]$ids[[site]], use.names = FALSE)
        if (anyDuplicated(ids) || any(ids %in% outer_ids[[site]])) {
          stop("Joint weight learning reused an outer evaluation observation.", call. = FALSE)
        }
        gradient_sites[[site]][ids] <- gradient_sites[[site]][ids] + derivative$sites[[site]]
      }
    }
    gradient <- list(target = gradient_sites[[1L]], source = gradient_sites[-1L],
      kink_cells = kink_cells, fold_source_cells = 2L * n_folds * K)
  }
  variance <- .weight_layer_variance(all_phi_tau, n_t, n_source, gradient)
  estimate <- mean(all_phi_tau)
  se <- sqrt(variance$variance)
  target_raw <- .assemble_target_pseudovalues(fold_info, n_t, ".aggregate_tate_arm_weights")
  target_variance <- .multisite_pseudovalue_variance(target_raw, as.integer(n_t))
  for (arm in names(phase2)) colnames(phase2[[arm]]$fold_weights) <- source_sites
  fold_weights <- do.call(cbind, lapply(phase2, `[[`, "fold_weights"))
  colnames(fold_weights) <- c(paste0("mu1:", source_sites), paste0("mu0:", source_sites))
  field_matrix <- function(field) {
    result <- do.call(cbind, lapply(phase2, `[[`, field))
    colnames(result) <- colnames(fold_weights)
    result
  }
  outer_estimates <- lapply(names(arms), function(arm) {
    vapply(seq_len(n_folds), function(fold) {
      info <- arms[[arm]]$intermediates$fold_info[[fold]]
      info$fold_target_estimate + sum(phase2[[arm]]$fold_weights[fold, ] *
        (info$fold_source_estimates - info$fold_target_estimate))
    }, numeric(1L))
  })
  source_estimates <- mu1_result$source_estimates - mu0_result$source_estimates
  list(estimate = estimate, se = se, variance = variance$variance,
    se_fixed_weights = sqrt(variance$fixed_variance), variance_fixed_weights = variance$fixed_variance,
    weight_layer = variance[c("indirect_variance", "cross_term", "kink_cells", "fold_source_cells")],
    ci_lower = estimate - Z_ALPHA_05 * se, ci_upper = estimate + Z_ALPHA_05 * se,
    target_only = list(estimate = mean(target_raw), variance = target_variance, se = sqrt(target_variance)),
    source_estimates = source_estimates, weights = colMeans(fold_weights), fold_weights = fold_weights,
    weights_by_arm = lapply(phase2, function(x) colMeans(x$fold_weights)),
    fold_weights_by_arm = lapply(phase2, `[[`, "fold_weights"),
    fold_lambdas = do.call(cbind, lapply(phase2, `[[`, "fold_lambdas")),
    fold_wald_statistics = field_matrix("fold_wald_statistics"),
    fold_penalty_coefficients = field_matrix("fold_penalty_coefficients"),
    fold_source_included = field_matrix("fold_source_included"),
    fold_weight_optimizer_iterations = lapply(phase2, `[[`, "fold_weight_optimizer_iterations"),
    fold_weight_psd_ridge = lapply(phase2, `[[`, "fold_weight_psd_ridge"),
    fold_aggregated_estimates = outer_estimates[[1L]] - outer_estimates[[2L]],
    aggregation_mode = aggregation_mode, aggregation_lambda_selection = lambda_selection,
    aggregation_lambda_rule = lambda_rule, aggregation_lambda_grid = aggregation_lambda_grid,
    aggregation_screening_rule = screening_rule,
    inference_scope = "fixed nuisance fits and selected active sets; full validity remains under review",
    clip_diagnostics = .combine_tate_clip_diagnostics(mu1_result, mu0_result),
    n_sites = K + 1L, n_folds = n_folds, N_all = N_all,
    crossfit_levels = max(mu1_result$crossfit_levels %||% 2L, mu0_result$crossfit_levels %||% 2L),
    all_phi_agg = all_phi_tau, all_phi_tau = all_phi_tau,
    estimand = "TATE", method = paste0(aggregation_mode, "_crossfit"), arm_results = arms,
    intermediates = list(fold_info = fold_info, inner_fold_info = inner_fold_info,
      target_only_phi = target_raw, phase2_by_arm = phase2, joint_fits = joint_fits,
      joint_moments = joint_moments,
      sample_sizes = list(n_t = n_t, n_source = n_source, N_all = N_all)))
}
