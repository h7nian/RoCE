# Shared inputs and fitting controls for RHC validation diagnostics.
rhc_validation_data <- function() {
  cohort <- RoCE::build_rhc_cohort(outcome = "death30", site_var = "ninsclas",
    site_recode = c("No insurance" = NA_character_, "Medicare & Medicaid" = NA_character_))
  data <- RoCE::build_rhc_data_split(cohort, K = 4L, target_site = "Private", seed = 42L)
  stopifnot(identical(digest::digest(data, algo = "sha256"),
    "a6d8135c01cd4c1113e134b024b2b440a2e31c3e21e2e0a41b9b4d14d6c5aa1a"))
  data
}

rhc_fold_views <- function(data, indices) {
  stopifnot(length(indices) >= 3L,
    identical(sort(as.integer(unlist(indices, use.names = FALSE))), seq_len(data$n)))
  views <- lapply(indices, function(index) {
    index <- as.integer(index)
    stopifnot(length(index) > 0L, length(unique(data$A[index])) == 2L)
    list(original_idx = index, n = length(index))
  })
  attr(views, ".data_ref") <- data
  views
}

rhc_saved_folds <- function(data, fit) {
  arms <- fit$arm_results[c("mu1", "mu0")]
  first <- arms$mu1$intermediates$fold_info
  second <- arms$mu0$intermediates$fold_info
  fields <- c("target_idx", "source_idx")
  stopifnot(length(first) == fit$n_folds, length(second) == fit$n_folds,
    identical(lapply(first, `[`, fields), lapply(second, `[`, fields)))
  sources <- setdiff(names(data), "t")
  list(target_folds = rhc_fold_views(data$t, lapply(first, `[[`, "target_idx")),
    source_folds = setNames(lapply(seq_along(sources), function(j) {
      rhc_fold_views(data[[sources[j]]], lapply(first, function(fold) fold$source_idx[[j]]))
    }), sources))
}

rhc_remove_case <- function(data, folds, site, record) {
  stopifnot(site %in% names(data), length(record) == 1L, is.numeric(record),
    is.finite(record), record == as.integer(record), record >= 1L, record <= data[[site]]$n)
  retained <- setdiff(seq_len(data[[site]]$n), as.integer(record))
  updated <- data
  for (field in c("X", "X_dagger", "W_outcome", "Z_site")) {
    updated[[site]][[field]] <- data[[site]][[field]][retained, , drop = FALSE]
  }
  for (field in c("A", "Y")) updated[[site]][[field]] <- data[[site]][[field]][retained]
  if (!is.null(data[[site]]$cv_group_id)) updated[[site]]$cv_group_id <- data[[site]]$cv_group_id[retained]
  updated[[site]]$n <- length(retained)
  old_views <- if (site == "t") folds$target_folds else folds$source_folds[[site]]
  indices <- lapply(old_views, function(fold) {
    matched <- match(fold$original_idx, retained)
    as.integer(matched[!is.na(matched)])
  })
  updated_folds <- folds
  replacement <- rhc_fold_views(updated[[site]], indices)
  if (site == "t") updated_folds$target_folds <- replacement else updated_folds$source_folds[[site]] <- replacement
  # The original preprocessing is deliberately fixed; only the designated record
  # is removed for this diagnostic. Every other patient's outer fold is retained.
  for (k in seq_along(indices)) stopifnot(identical(retained[indices[[k]]], setdiff(old_views[[k]]$original_idx, record)))
  list(data = updated, folds = updated_folds, retained_original_rows = retained)
}

rhc_validation_fit <- function(data, folds, radius, checkpoint, cores = 2L, nlambda = 100L,
                               shared_cache = NULL) {
  stopifnot(radius %in% c(2, 3, 5), cores >= 1L)
  RoCE::run_tate_crossfit(data_split = data, n_folds = length(folds$target_folds),
    communication_mode = "one_round", family = "binomial", n_cores = as.integer(cores),
    parallel_arms = cores > 1L, M_tau = radius, M_tau_inference = radius,
    lambda_selection = .5, nlambda_init = as.integer(nlambda), nuisance_lambda_rule = "min",
    precomputed_folds = folds, aggregation_mode = "joint_tate",
    target_nuisance_method = "hou_calibrated", source_validation_method = "outer_fit",
    crossfit_layers = 2L, calibration_control = list(recipe = "score_derivative", target_radius = 5),
    calibration_layout = "compact", nuisance_solver = "proximal_newton",
    nuisance_tol = 1e-10, nuisance_cv_certificate = TRUE,
    nuisance_cache_dir = shared_cache, checkpoint_dir = checkpoint, verbose = FALSE)
}

rhc_fit_result <- function(fit) {
  data.frame(method = c("RoCE", "Target-only"),
    estimate = c(fit$estimate, fit$target_only$estimate),
    se = c(fit$se, fit$target_only$se),
    ci_lower = c(fit$ci_lower, fit$target_only$estimate - qnorm(.975) * fit$target_only$se),
    ci_upper = c(fit$ci_upper, fit$target_only$estimate + qnorm(.975) * fit$target_only$se))
}
