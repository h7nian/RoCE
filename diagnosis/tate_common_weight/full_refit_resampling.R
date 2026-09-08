# Diagnostic full nuisance refits conditional on the original outer partition.
# This is NOT a validated confidence-interval procedure. See FULL_REFIT_REVIEW.md.
# Requires the explicitly selected, installed RoCE package; never load_all().

roce_refit_resample <- function(data_split, fold_indices, seed,
                              identity = FALSE) {
  fail <- function(...) stop("full-refit resampling: ", ..., call. = FALSE)
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 1 || seed != floor(seed) || seed > .Machine$integer.max) {
    fail("seed must be one positive integer.")
  }
  if (!is.logical(identity) || length(identity) != 1L || is.na(identity)) {
    fail("identity must be TRUE or FALSE.")
  }
  sites <- names(data_split)
  if (!is.list(data_split) || length(sites) < 2L || anyNA(sites) ||
      any(!nzchar(sites)) || anyDuplicated(sites) || !"t" %in% sites ||
      !is.list(fold_indices) || !identical(names(fold_indices), sites)) {
    fail("data and fold lists must have identical unique site names including t.")
  }
  matrix_fields <- c("X", "X_dagger", "Z_site_true", "W_outcome_true",
                     "Z_site", "W_outcome")
  vector_fields <- c("A", "Y")
  expected_fields <- c("n", matrix_fields, vector_fields)
  fold_counts <- lengths(fold_indices)
  if (any(fold_counts < 3L) || length(unique(fold_counts)) != 1L) {
    fail("every site must have the same number of folds, at least three.")
  }
  for (site in sites) {
    data <- data_split[[site]]
    n <- data$n
    if (!is.list(data) || anyDuplicated(names(data)) ||
        !setequal(names(data), expected_fields) ||
        !is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n < 1 || n != floor(n) || n > .Machine$integer.max) {
      fail("invalid site schema for ", site, "; no fields may be silently dropped.")
    }
    for (field in matrix_fields) {
      value <- data[[field]]
      if (!is.matrix(value) || !is.numeric(value) || nrow(value) != n ||
          any(!is.finite(value))) fail("invalid matrix ", site, "$", field, ".")
    }
    for (field in vector_fields) {
      value <- data[[field]]
      if (!is.numeric(value) || length(value) != n || any(!is.finite(value))) {
        fail("invalid vector ", site, "$", field, ".")
      }
    }
    if (any(!data$A %in% c(0, 1))) fail("treatment must be binary.")
    indices <- fold_indices[[site]]
    if (any(vapply(indices, function(idx) {
      !is.numeric(idx) || length(idx) == 0L || any(!is.finite(idx)) ||
        any(idx != floor(idx)) || any(idx < 1L) || any(idx > n)
    }, logical(1L))) ||
        !identical(sort(as.integer(unlist(indices))), seq_len(n))) {
      fail("folds must partition each observation exactly once at ", site, ".")
    }
  }

  # Sample whole rows once per site/fold. Both arms and every overlapping
  # training fold see the same copies. Duplicates never cross outer folds.
  maps <- RoCE:::with_seed(as.integer(seed), {
    out <- setNames(vector("list", length(sites)), sites)
    for (site in sort(sites)) {
      rows <- seq_len(data_split[[site]]$n)
      for (idx in fold_indices[[site]]) {
        if (!identity) rows[idx] <- idx[sample.int(length(idx), replace = TRUE)]
      }
      out[[site]] <- rows
    }
    out
  })
  sampled <- data_split
  folds <- setNames(vector("list", length(sites)), sites)
  for (site in sites) {
    rows <- maps[[site]]
    for (field in matrix_fields) {
      sampled[[site]][[field]] <- data_split[[site]][[field]][rows, , drop = FALSE]
    }
    for (field in vector_fields) {
      sampled[[site]][[field]] <- data_split[[site]][[field]][rows]
    }
    folds[[site]] <- lapply(fold_indices[[site]], function(idx) {
      list(original_idx = as.integer(idx), n = length(idx))
    })
    attr(folds[[site]], ".data_ref") <- sampled[[site]]
  }
  list(
    data_split = sampled,
    precomputed_folds = list(target_folds = folds$t,
                            source_folds = folds[setdiff(sites, "t")]),
    original_row_ids = maps,
    multiplicities = stats::setNames(lapply(sites, function(site) {
      tabulate(maps[[site]], nbins = data_split[[site]]$n)
    }), sites),
    seed = as.integer(seed), identity = identity,
    resampling_scheme = "site_by_original_outer_fold_multinomial"
  )
}

roce_refit_fixed_weight_estimate <- function(fit, reference_weights) {
  RoCE:::.validate_weight_bootstrap_result(fit)
  if (!is.matrix(reference_weights) || !is.numeric(reference_weights) ||
      any(!is.finite(reference_weights)) ||
      !identical(dim(reference_weights), dim(fit$fold_weights)) ||
      !identical(colnames(reference_weights), names(fit$weights))) {
    stop("full-refit: reference weights must match fitted folds and source order.",
         call. = FALSE)
  }
  sizes <- fit$intermediates$sample_sizes
  phi <- RoCE:::.compute_phase3_all_phi(
    n_folds = fit$n_folds, fold_weights = reference_weights,
    fold_info = fit$intermediates$fold_info, n_t = sizes$n_t,
    n_source_full = sizes$n_source, N_all = sizes$N_all,
    source_sites = names(fit$weights), K = length(fit$weights), verbose = FALSE
  )
  mean(phi)
}

roce_refit_matched_fixed_nuisance <- function(reference_fit, multiplicities) {
  RoCE:::.validate_weight_bootstrap_result(reference_fit)
  sources <- names(reference_fit$weights)
  sizes <- reference_fit$intermediates$sample_sizes
  sites <- c("t", sources)
  expected_sizes <- c(sizes$n_t, sizes$n_source)
  if (!is.list(multiplicities) || !identical(names(multiplicities), sites)) {
    stop("matched resampling: multiplicity names must match reference sites.")
  }
  for (j in seq_along(sites)) {
    counts <- multiplicities[[j]]
    if (!is.numeric(counts) || length(counts) != expected_sizes[[j]] ||
        any(!is.finite(counts)) || any(counts < 0) ||
        any(counts != floor(counts)) || sum(counts) != expected_sizes[[j]]) {
      stop("matched resampling: invalid counts for ", sites[[j]], ".")
    }
    for (info in reference_fit$intermediates$fold_info) {
      idx <- if (j == 1L) info$target_idx else info$source_idx[[j - 1L]]
      if (sum(counts[idx]) != length(idx)) {
        stop("matched resampling: counts must preserve each outer-fold size.")
      }
    }
  }
  target_counts <- multiplicities$t
  source_counts <- multiplicities[sources]
  inner <- lapply(reference_fit$intermediates$inner_fold_info, function(info) {
    RoCE:::.reweight_tate_inner_info(info, target_counts, source_counts, sources)
  })
  phase2 <- RoCE:::.compute_phase2_weights(
    n_folds = reference_fit$n_folds, inner_fold_info = inner,
    fold_info = reference_fit$intermediates$fold_info,
    lambda_selection = reference_fit$aggregation_lambda_selection,
    K = length(sources), verbose = FALSE,
    lambda_rule = reference_fit$aggregation_lambda_rule,
    screening_rule = reference_fit$aggregation_screening_rule
  )
  estimate <- function(weights) RoCE:::.weight_bootstrap_outer_estimate(
    reference_fit$intermediates$fold_info, weights, target_counts,
    source_counts, sources, sizes$n_t, sizes$n_source
  )
  colnames(phase2$fold_weights) <- sources
  list(estimate_original_weights = estimate(reference_fit$fold_weights),
       estimate_relearned_weights = estimate(phase2$fold_weights),
       fold_weights = phase2$fold_weights)
}

roce_refit_once <- function(data_split, reference_fit, seed, identity = FALSE,
                           n_cores = 1L, parallel_arms = FALSE,
                           nuisance_cv_grouping = c("row", "origin")) {
  nuisance_cv_grouping <- match.arg(nuisance_cv_grouping)
  if (nuisance_cv_grouping == "origin" &&
      (!identical(as.numeric(n_cores), 1) || !identical(parallel_arms, FALSE))) {
    stop("Origin-grouped audit requires sequential R workers so every CV call is recorded.")
  }
  RoCE:::.validate_weight_bootstrap_result(reference_fit)
  if (!identical(reference_fit$communication_mode, "one_round") ||
      !is.numeric(reference_fit$aggregation_lambda_selection) ||
      !identical(reference_fit$aggregation_screening_rule, "soft_penalty")) {
    stop("full-refit diagnostic requires one-round soft weights and fixed lambda.",
         call. = FALSE)
  }
  sites <- c("t", names(reference_fit$weights))
  if (!identical(names(data_split), sites)) {
    stop("full-refit: data site order must match the reference fit.", call. = FALSE)
  }
  info <- reference_fit$intermediates$fold_info
  fold_indices <- c(list(t = lapply(info, `[[`, "target_idx")),
    setNames(lapply(seq_along(reference_fit$weights), function(j) {
      lapply(info, function(fold) fold$source_idx[[j]])
    }), names(reference_fit$weights)))
  sampled <- roce_refit_resample(data_split, fold_indices, seed, identity)
  if (nuisance_cv_grouping == "origin") {
    for (site in sites) sampled$data_split[[site]]$cv_group_id <- sampled$original_row_ids[[site]]
    attr(sampled$precomputed_folds$target_folds, ".data_ref") <- sampled$data_split$t
    for (site in setdiff(sites, "t")) {
      attr(sampled$precomputed_folds$source_folds[[site]], ".data_ref") <- sampled$data_split[[site]]
    }
  }
  matched_fixed_nuisance <- roce_refit_matched_fixed_nuisance(
    reference_fit, sampled$multiplicities
  )
  nlambda <- unique(unlist(lapply(reference_fit$arm_results, function(arm) {
    lapply(arm$fold_results, function(fold) {
      vapply(fold$source_results, function(source) source$nuisance_nlambda,
             numeric(1L))
    })
  })))
  if (length(nlambda) != 1L || !is.finite(nlambda) || nlambda < 2L ||
      nlambda != floor(nlambda)) {
    stop("full-refit: reference nuisance grid is missing or inconsistent.",
         call. = FALSE)
  }
  fit_function <- function() RoCE::run_tate_crossfit(
    data_split = sampled$data_split, n_folds = reference_fit$n_folds,
    communication_mode = "one_round", verbose = FALSE,
    lambda_selection = reference_fit$aggregation_lambda_selection,
    lambda_rule = reference_fit$aggregation_lambda_rule,
    M_tau = reference_fit$M_tau,
    M_tau_inference = reference_fit$M_tau_inference,
    n_cores = n_cores, nlambda_init = as.integer(nlambda),
    family = reference_fit$family, use_lambda_cache = TRUE,
    precomputed_folds = sampled$precomputed_folds,
    nuisance_lambda_rule = reference_fit$nuisance_lambda_rule,
    parallel_arms = parallel_arms, screening_rule = "soft_penalty"
    # Deliberately omit target-only caches: all nuisance fits must be fresh.
  )
  cv_partitions <- NULL
  failure_stage <- "nuisance_refit"
  if (nuisance_cv_grouping == "origin") {
    audited <- roce_with_nuisance_cv_audit(fit_function)
    fitted <- audited$fitted
    cv_partitions <- audited$partitions
    if (!inherits(fitted, "error") &&
        (!length(cv_partitions) || !all(vapply(cv_partitions, roce_cv_partition_record_valid, logical(1L))))) {
      fitted <- simpleError("Nuisance CV partition audit failed or recorded no calls.")
      failure_stage <- "nuisance_cv_partition_audit"
    }
  } else {
    fitted <- tryCatch(fit_function(), error = function(e) e)
  }
  if (inherits(fitted, "error")) {
    return(list(failure_message = conditionMessage(fitted),
                failure_stage = failure_stage,
                matched_fixed_nuisance = matched_fixed_nuisance,
                original_row_ids = sampled$original_row_ids,
                multiplicities = sampled$multiplicities,
                nuisance_cv_grouping = nuisance_cv_grouping,
                nuisance_cv_partitions = cv_partitions))
  }
  fixed <- roce_refit_fixed_weight_estimate(fitted, reference_fit$fold_weights)
  if (identity &&
      (abs(fitted$estimate - reference_fit$estimate) > 1e-10 ||
       abs(fitted$se - reference_fit$se) > 1e-10 ||
       max(abs(fitted$fold_weights - reference_fit$fold_weights)) > 1e-10 ||
       abs(fixed - reference_fit$estimate) > 1e-10 ||
       abs(matched_fixed_nuisance$estimate_original_weights - reference_fit$estimate) > 1e-10 ||
       abs(matched_fixed_nuisance$estimate_relearned_weights - reference_fit$estimate) > 1e-10)) {
    stop("full-refit: identity draw did not reproduce the reference fit.",
         call. = FALSE)
  }
  list(
    fitted = fitted, estimate_refit_relearned_weights = fitted$estimate,
    estimate_refit_original_weights = fixed,
    matched_fixed_nuisance = matched_fixed_nuisance,
    original_row_ids = sampled$original_row_ids,
    multiplicities = sampled$multiplicities, seed = sampled$seed,
    identity = identity, resampling_scheme = sampled$resampling_scheme,
    nuisance_refit = TRUE, inference_validated = FALSE,
    nuisance_cv_grouping = nuisance_cv_grouping,
    nuisance_cv_partitions = cv_partitions
  )
}
