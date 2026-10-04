# Arm-level variance diagnostics at the actual selected aggregation weights.
# Uses the checked RoCE kernels; nuisance fits and selected active sets are held
# fixed in the eta-sensitivity calculation, as in the estimator's reported SE.

layer_outer_fit_fingerprint <- function(fit) {
  # Exclude eta, inner fits, timing, and training-score retention metadata.
  # Comparing paired repeats checks the actual final source coefficients and
  # target predictions used on the outer evaluation observations.
  values <- lapply(fit$arm_results[c("mu1", "mu0")], function(arm)
    lapply(arm$fold_results, function(fold) list(
      target = lapply(fold$target_only[c("estimate", "m_pred", "prop_scores", "varphi_ot")], as.numeric),
      sources = lapply(fold$source_results, function(source)
        list(outcome = as.numeric(source$alpha_ts), weight = as.numeric(source$gamma_s))))))
  digest::digest(values, algo = "sha256")
}

summarize_layer_arms <- function(fit, data_split, keep_influence = FALSE) {
  arms <- fit$arm_results[c("mu1", "mu0")]
  sources <- setdiff(names(data_split), "t")
  sizes <- vapply(data_split[c("t", sources)], function(site) as.integer(site$n), integer(1L))
  total <- sum(sizes)
  count <- length(sources)
  folds <- fit$n_folds
  mode <- fit$aggregation_mode
  weights <- if (is.null(fit$fold_weights_by_arm)) {
    list(mu1 = fit$fold_weights, mu0 = fit$fold_weights)
  } else fit$fold_weights_by_arm
  raw <- lapply(names(arms), function(arm) RoCE:::.compute_phase3_all_phi(folds, weights[[arm]],
    arms[[arm]]$intermediates$fold_info, sizes[1L], sizes[-1L], total, sources, count, FALSE))
  names(raw) <- names(arms)
  direct <- lapply(raw, function(values) {
    groups <- split(values, rep(seq_along(sizes), sizes))
    unlist(lapply(groups, function(group) (group - mean(group))/total), use.names = FALSE)
  })
  rule <- fit$aggregation_screening_rule
  if (is.null(rule)) rule <- "soft_penalty"
  if (mode == "joint_tate") {
    increments <- lapply(arms, function(arm)
      RoCE:::.outer_fold_increments(arm$intermediates$fold_info, sizes[1L], sizes[-1L]))
    gradients <- lapply(arms, function(arm) lapply(sizes, function(size) numeric(size)))
    for (fold in seq_len(folds)) {
      records <- RoCE:::.joint_inner_records(arms$mu1$intermediates$inner_fold_info[[fold]],
                                           arms$mu0$intermediates$inner_fold_info[[fold]])
      for (arm in names(arms)) {
        increment <- numeric(2L * count)
        positions <- if (arm == "mu1") seq_len(count) else count + seq_len(count)
        increment[positions] <- increments[[arm]][fold, ]
        derivative <- RoCE:::.joint_weight_derivative(fit$intermediates$joint_moments[[fold]],
          fit$intermediates$joint_fits[[fold]], rule, increment)
        for (site in seq_along(sizes)) {
          ids <- unlist(records$ids[[site]], use.names = FALSE)
          outer <- arms$mu1$intermediates$fold_info[[fold]]
          held_out <- if (site == 1L) outer$target_idx else outer$source_idx[[site - 1L]]
          if (anyDuplicated(ids) || any(ids %in% held_out)) stop("Eta diagnostics reused an outer evaluation row")
          gradients[[arm]][[site]][ids] <- gradients[[arm]][[site]][ids] + derivative$sites[[site]]
        }
      }
    }
    indirect <- lapply(gradients, unlist, use.names = FALSE)
  } else {
    indirect <- lapply(names(arms), function(arm) {
      phase <- if (mode == "separate_arms") fit$intermediates$phase2_by_arm[[arm]] else fit
      training <- if (mode == "separate_arms") arms[[arm]]$intermediates$inner_fold_info else
        fit$intermediates$inner_fold_info
      gradient <- RoCE:::.weight_layer_gradient(arms[[arm]]$intermediates$fold_info, training,
        weights[[arm]], phase$fold_lambdas, rule, phase$fold_weight_psd_ridge, sizes[1L], sizes[-1L])
      unlist(c(list(gradient$target), gradient$source), use.names = FALSE)
    })
    names(indirect) <- names(arms)
  }
  corrected <- Map(`+`, direct, indirect)
  for (arm in names(corrected)) {
    groups <- split(corrected[[arm]], rep(seq_along(sizes), sizes))
    if (any(vapply(groups, function(group) abs(sum(group)) > 1e-10 * max(1, sum(abs(group))), logical(1L)))) {
      stop("An arm's estimator contributions do not center within site")
    }
  }
  covariance <- crossprod(do.call(cbind, corrected))
  fixed_covariance <- crossprod(do.call(cbind, direct))
  estimates <- vapply(raw, mean, numeric(1L))
  contrast <- c(1, -1)
  variance <- drop(crossprod(contrast, covariance %*% contrast))
  fixed_variance <- drop(crossprod(contrast, fixed_covariance %*% contrast))
  if (abs(estimates[["mu1"]] - estimates[["mu0"]] - fit$estimate) > 1e-12 ||
      abs(variance - fit$variance) > 1e-12 || abs(fixed_variance - fit$variance_fixed_weights) > 1e-12) {
    stop("Arm summaries do not reproduce the fitted TATE estimate and both variance formulas")
  }
  result <- list(estimates = estimates, covariance = covariance, fixed_covariance = fixed_covariance,
       scope = "Eta sensitivity holds nuisance fits and active sets fixed; compare with empirical repeated-simulation variance")
  if (keep_influence) {
    result$influence <- do.call(cbind, corrected)
    result$influence_fixed <- do.call(cbind, direct)
  }
  result
}

layer_inference_packet <- function(fit, data_split, arm_truth) {
  # Retain score records for later variance audits without duplicating the
  # full high-dimensional nuisance fits and raw feature matrices per repeat.
  list(n_folds = fit$n_folds, crossfit_layers = fit$crossfit_levels,
    source_validation_method = fit$source_validation_method,
    site_sizes = vapply(data_split, function(site) as.integer(site$n), integer(1L)),
    arm_truth = arm_truth,
    arms = lapply(fit$arm_results[c("mu1", "mu0")], function(arm)
      arm$intermediates[c("fold_info", "inner_fold_info")]),
    calibration_control = fit$calibration_control,
    aggregation_lambda = fit$aggregation_lambda_selection,
    aggregation_screening_rule = fit$aggregation_screening_rule)
}
