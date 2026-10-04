# Research-only postprocessing. The nuisance fits and original joint optimizer
# stay fixed. A common fold multiplier scales both arms and all sources.
rescale_joint_tate <- function(fit, data, fold_scale) {
  stopifnot(fit$aggregation_mode == "joint_tate", length(fold_scale) == fit$n_folds,
    all(is.finite(fold_scale)), all(fold_scale >= 0 & fold_scale <= 1))
  arms <- fit$arm_results[c("mu1", "mu0")]
  sources <- setdiff(names(data), "t")
  sizes <- vapply(data[c("t", sources)], function(site) as.integer(site$n), integer(1L))
  count <- length(sources); folds <- fit$n_folds; total <- sum(sizes)
  stopifnot(count > 0L)
  weights <- lapply(fit$fold_weights_by_arm, function(value) value * fold_scale)
  raw <- lapply(names(arms), function(arm) RoCE:::.compute_phase3_all_phi(
    folds, weights[[arm]], arms[[arm]]$intermediates$fold_info,
    sizes[1L], sizes[-1L], total, sources, count, FALSE))
  contrast <- raw[[1L]] - raw[[2L]]
  increments <- cbind(RoCE:::.outer_fold_increments(arms$mu1$intermediates$fold_info, sizes[1L], sizes[-1L]),
    -RoCE:::.outer_fold_increments(arms$mu0$intermediates$fold_info, sizes[1L], sizes[-1L]))
  gradients <- lapply(sizes, function(size) numeric(size)); kinks <- 0L
  for (fold in seq_len(folds)) {
    records <- RoCE:::.joint_inner_records(arms$mu1$intermediates$inner_fold_info[[fold]],
      arms$mu0$intermediates$inner_fold_info[[fold]])
    derivative <- RoCE:::.joint_weight_derivative(fit$intermediates$joint_moments[[fold]],
      fit$intermediates$joint_fits[[fold]], fit$aggregation_screening_rule,
      fold_scale[fold] * increments[fold, ])
    kinks <- kinks + derivative$kinks
    outer <- arms$mu1$intermediates$fold_info[[fold]]
    excluded <- c(list(outer$target_idx), outer$source_idx)
    for (site in seq_along(sizes)) {
      ids <- unlist(records$ids[[site]], use.names = FALSE)
      stopifnot(!anyDuplicated(ids), !any(ids %in% excluded[[site]]))
      gradients[[site]][ids] <- gradients[[site]][ids] + derivative$sites[[site]]
    }
  }
  variance <- RoCE:::.weight_layer_variance(contrast, sizes[1L], sizes[-1L],
    list(target = gradients[[1L]], source = gradients[-1L], kink_cells = kinks,
      fold_source_cells = 2L * folds * count))
  direct <- unlist(lapply(split(contrast, rep(seq_along(sizes), sizes)),
    function(values) (values - mean(values)) / total), use.names = FALSE)
  influence <- direct + unlist(gradients, use.names = FALSE)
  stopifnot(abs(sum(influence^2) - variance$variance) < 1e-12)
  list(estimate = mean(contrast), se_first_order = sqrt(variance$variance),
    se_fixed_weights = sqrt(variance$fixed_variance), influence = influence,
    fold_scale = fold_scale, weights_by_arm = weights,
    scope = "Fixed nuisances/active sets; mass-factor uncertainty is first-order zero only under the stated consistent valid-candidate conditions; requires empirical validation")
}

source_weight_mass <- function(fit, data) {
  stopifnot(fit$aggregation_mode == "joint_tate", fit$family == "binomial")
  sources <- setdiff(names(data), "t")
  masses <- matrix(NA_real_, fit$n_folds, 2L * length(sources),
    dimnames = list(NULL, c(paste0("mu1:", sources), paste0("mu0:", sources))))
  for (arm in c("mu1", "mu0")) {
    a <- if (arm == "mu1") 1L else 0L
    result <- fit$arm_results[[arm]]
    for (fold in seq_len(fit$n_folds)) for (j in seq_along(sources)) {
      name <- sources[j]; selected <- result$fold_results[[fold]]$source_results[[name]]
      ids <- result$intermediates$fold_info[[fold]]$source_idx[[j]]
      site <- data[[name]]; outcome_design <- site$W_outcome[ids, , drop = FALSE]
      prediction <- RoCE:::predict_glm_cpp(outcome_design, selected$alpha_ts, 1L, 1L)
      weight <- RoCE:::calculate_correction_term_cpp(site$Z_site[ids, , drop = FALSE],
        site$A[ids], as.numeric(prediction) + 1, selected$gamma_s, selected$alpha_ts,
        outcome_design, fit$M_tau_inference, 1L, 1L, a)$correction_components
      reconstructed <- as.numeric(weight) * (site$Y[ids] - as.numeric(prediction))
      stopifnot(max(abs(reconstructed - selected$correction_components)) < 1e-10)
      masses[fold, paste0(arm, ":", name)] <- mean(weight)
    }
  }
  stopifnot(all(is.finite(masses)), all(masses >= 0))
  masses
}

stabilize_joint_tate <- function(fit, data, method = c("none", "weight_mass")) {
  method <- match.arg(method)
  mass <- source_weight_mass(fit, data)
  scale <- if (method == "none") rep(1, fit$n_folds) else 1 / pmax(1, apply(mass, 1L, max))
  result <- rescale_joint_tate(fit, data, scale)
  result$weight_mass <- mass; result$method <- method
  result
}
