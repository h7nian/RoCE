# cross_fitting_algorithms.R - Two-Layer Cross-fitting variants for FACE algorithm
# Following docs/main.tex (Two-Level Cross-fitting appendix section)
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# This file implements estimators for the POTENTIAL OUTCOME MEAN at the target site:
#
#   mu^1_t = E_t[Y(1)]  (expected outcome under treatment for target population)
#
# This is NOT the Average Treatment Effect (ATE). To compute ATE:
#   ATE = mu^1_t - mu^0_t = E_t[Y(1)] - E_t[Y(0)]
#
# Run the algorithm twice with A_val = 1 and A_val = 0 to get both potential
# outcome means, then compute the difference.
#
# =============================================================================
# This file implements the main FACE-HD algorithms with optional parallelization
# for processing multiple source sites simultaneously.

# Parallelization utilities are in R/parallel_utils.R (setup_parallel, parallel_lapply)

.make_plugin_block_design <- function(feature_list, caller = ".make_plugin_block_design") {
  if (length(feature_list) == 0L) {
    stop(sprintf("%s: feature_list must contain at least one matrix.", caller))
  }

  feature_list <- lapply(feature_list, as.matrix)
  p <- ncol(feature_list[[1L]])
  if (is.null(p)) {
    stop(sprintf("%s: feature entries must be matrices.", caller))
  }
  bad <- which(vapply(feature_list, ncol, integer(1L)) != p)
  if (length(bad) > 0L) {
    stop(sprintf("%s: all feature matrices must have the same number of columns.", caller))
  }

  block_width <- p + 1L
  total_n <- sum(vapply(feature_list, nrow, integer(1L)))
  out <- matrix(0, nrow = total_n, ncol = length(feature_list) * block_width)

  row_start <- 1L
  for (i in seq_along(feature_list)) {
    mat <- feature_list[[i]]
    n_i <- nrow(mat)
    if (n_i == 0L) next
    row_end <- row_start + n_i - 1L
    col_start <- (i - 1L) * block_width + 1L
    col_end <- i * block_width
    out[row_start:row_end, col_start:col_end] <- cbind(1, mat)
    row_start <- row_end + 1L
  }

  out
}

.stack_fold_field <- function(fold_data_list, field, caller = ".stack_fold_field") {
  values <- lapply(fold_data_list, `[[`, field)
  if (length(values) == 0L) {
    stop(sprintf("%s: no values to stack for field '%s'.", caller, field))
  }
  if (is.matrix(values[[1L]])) {
    do.call(rbind, values)
  } else {
    unlist(values, use.names = FALSE)
  }
}

.average_numeric_list <- function(values, caller = ".average_numeric_list") {
  if (length(values) == 0L) {
    stop(sprintf("%s: values must be non-empty.", caller))
  }
  Reduce(`+`, values) / length(values)
}

#' Process a single source site for cross-fitting (shared by both algorithms)
#'
#' Unified helper that computes calibrated parameters and final estimates for
#' a source site. Used by both two-round and one-round algorithms; the only
#' difference is how initial params and target summaries are accessed, which
#' is controlled by the \code{get_fold_inputs} callback.
#'
#' @param s Source site name
#' @param source_folds List of source site folds
#' @param target_folds List of target folds
#' @param k1 Main fold index
#' @param n_folds Total number of folds
#' @param A_val Treatment value
#' @param M_tau Truncation parameter
#' @param data_split Original data split used for source metadata
#' @param combine_cache Environment for caching combine_folds() results (or NULL)
#' @param get_fold_inputs Function(s, k2) returning list(mean_phi, mean_grad_psi_init, alpha_init)
#'   for the given source site and secondary fold. Encapsulates the algorithm-specific
#'   data access pattern (two-round vs one-round).
#' @param family_int Integer code for GLM family (0=gaussian, 1=binomial)
#' @param link_int Integer code for link function (0=identity, 1=logit)
#' @param use_lambda_cache Logical. If TRUE, reuse selected nuisance lambdas
#'   within each source-site k2 loop to reduce repeated CV.
#' @param nuisance_lambda_rule CV selection rule for nuisance fits:
#'   \code{"min"} (default) uses \code{lambda.min}; \code{"1se"} uses
#'   \code{lambda.1se}.
#' @return List with source site results (mu_ts, gamma_s, alpha_ts, delta_ts, mu_pred_ts, n_s)
process_source_site <- function(s, source_folds, target_folds, k1, n_folds,
                                A_val, M_tau, M_tau_inference = M_TAU_INFERENCE_DEFAULT,
                                data_split,
                                combine_cache = NULL, get_fold_inputs,
                                family_int = 1L, link_int = 1L,
                                use_lambda_cache = TRUE,
                                nuisance_lambda_rule = c("min", "1se")) {
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(nuisance_lambda_rule, "process_source_site")

  secondary_folds <- setdiff(1:n_folds, k1)
  label <- "crossfit"
  
  # Initialize fold-specific storage. The calibrated nuisance itself is obtained
  # by one fold-summed optimization after all secondary-fold plug-ins are ready.
  gamma_init_list <- list()
  alpha_init_list <- list()
  source_calib_list <- list()
  mean_grad_psi_list <- list()
  
  # Lambda caching (optional): reuse selected lambda values within the same
  # source site's k2 loop. This remains valid under parallel execution because
  # it is local state (no cross-worker shared mutable cache).
  cached_lambda_init_dr <- NULL
  
  # Warm-start: pass previous k2's solution to accelerate convergence.
  # Across k2 iterations the training/calibration data changes by only 1 fold,
  # so optimal parameters are similar. Warm-starting reduces iterations by 3-10x.
  prev_gamma_init <- NULL
  
  # For each secondary fold k2, compute calibrated parameters
  for (k2 in secondary_folds) {
    source_calib <- materialize_fold(source_folds[[s]], k2)
    
    if (source_calib$n == 0) {
      stop(sprintf("%s: source site '%s' fold k2=%d has 0 calibration observations.",
                   label, s, k2))
    }
    
    training_folds <- setdiff(1:n_folds, c(k1, k2))
    if (length(training_folds) == 0) {
      stop(sprintf("%s: source site '%s' fold k1=%d, k2=%d leaves 0 training folds.",
                   label, s, k1, k2))
    }
    
    # Use cache if available, otherwise compute directly
    cache_key <- paste(s, paste(sort(training_folds), collapse = "_"), sep = ":")
    if (!is.null(combine_cache) && !is.null(combine_cache[[cache_key]])) {
      source_train <- combine_cache[[cache_key]]
    } else {
      source_train <- combine_folds(source_folds[[s]], training_folds)
      if (!is.null(combine_cache)) combine_cache[[cache_key]] <- source_train
    }
    
    if (source_train$n == 0) {
      stop(sprintf("%s: source site '%s' training folds {%s} have 0 observations.",
                   label, s, paste(training_folds, collapse = ",")))
    }
    
    # Get algorithm-specific inputs via callback
    inputs <- get_fold_inputs(s, k2)
    mean_phi_k2 <- inputs$mean_phi
    mean_grad_psi_init_k2 <- inputs$mean_grad_psi_init
    alpha_init_k1_k2 <- inputs$alpha_init
    
    # Compute initial density ratio  — eq:gamma_init in main.tex
    gamma_init_k1_k2 <- fit_initial_density_ratio(
      source_train$Z_site, source_train$A, mean_phi_k2,
      lambda = if (isTRUE(use_lambda_cache)) cached_lambda_init_dr else NULL,
      A_val = A_val, warm_start = prev_gamma_init,
      lambda_rule = nuisance_lambda_rule
    )
    if (isTRUE(use_lambda_cache) && is.null(cached_lambda_init_dr)) {
      cached_lambda_init_dr <- attr(gamma_init_k1_k2, "lambda_used")
    }
    prev_gamma_init <- as.numeric(gamma_init_k1_k2)
    
    k2_key <- paste0("k2_", k2)
    gamma_init_list[[k2_key]] <- gamma_init_k1_k2
    alpha_init_list[[k2_key]] <- alpha_init_k1_k2
    source_calib_list[[k2_key]] <- source_calib
    mean_grad_psi_list[[k2_key]] <- mean_grad_psi_init_k2
  }
  
  # Fold-summed SMMAL-style calibrated optimization over all secondary folds.
  # The existing C++ optimizers accept one plug-in vector, so fold-specific
  # plug-ins are represented through block designs. C++ prepends a global
  # intercept; we set its coefficient to zero and carry fold-specific intercepts
  # inside the block design.
  if (length(source_calib_list) == 0 || length(gamma_init_list) == 0 || length(alpha_init_list) == 0) {
    stop(sprintf("process_source_site: all calibration folds failed for source site '%s' (outer fold k1=%d). Cannot form method-aligned nuisance estimates.",
                 s, k1))
  } else {
    Z_cal_stack <- .stack_fold_field(source_calib_list, "Z_site", "process_source_site")
    W_cal_stack <- .stack_fold_field(source_calib_list, "W_outcome", "process_source_site")
    A_cal_stack <- .stack_fold_field(source_calib_list, "A", "process_source_site")
    Y_cal_stack <- .stack_fold_field(source_calib_list, "Y", "process_source_site")

    mean_grad_psi_avg <- .average_numeric_list(mean_grad_psi_list, "process_source_site")

    W_plugin_block <- .make_plugin_block_design(
      lapply(source_calib_list, `[[`, "W_outcome"),
      caller = "process_source_site(gamma fold-summed calibration)"
    )
    alpha_plugin_block <- c(0, unlist(alpha_init_list, use.names = FALSE))

    gamma_final_k1 <- fit_unified_density_ratio(
      Z_site = Z_cal_stack, A = A_cal_stack,
      mean_grad_psi = mean_grad_psi_avg, alpha_init = alpha_plugin_block,
      lambda = NULL,
      calibrated = TRUE, M_tau = M_tau,
      W_outcome = W_plugin_block, A_val = A_val,
      family_int = family_int, link_int = link_int,
      lambda_rule = nuisance_lambda_rule
    )

    Z_plugin_block <- .make_plugin_block_design(
      lapply(source_calib_list, `[[`, "Z_site"),
      caller = "process_source_site(alpha fold-summed calibration)"
    )
    gamma_plugin_block <- c(0, unlist(gamma_init_list, use.names = FALSE))

    alpha_final_k1 <- fit_unified_outcome(
      W_outcome = W_cal_stack, Y = Y_cal_stack, A = A_cal_stack,
      A_val = A_val, gamma_s = gamma_plugin_block,
      lambda = NULL,
      family_int = family_int, link_int = link_int,
      calibrated = TRUE, M_tau = M_tau, Z_site = Z_plugin_block,
      lambda_rule = nuisance_lambda_rule
    )
  }
  
  # Compute correction term on main fold  — eq:site_membership in main.tex
  source_main_fold <- materialize_fold(source_folds[[s]], k1)
  
  if (source_main_fold$n == 0) {
    stop(sprintf("%s: source site '%s' main fold k1=%d has 0 observations; cannot compute source correction term.",
                 label, s, k1))
  }
  
  correction_result <- calculate_correction_term_cpp(
    source_main_fold$Z_site, source_main_fold$A, source_main_fold$Y,
    gamma_final_k1, alpha_final_k1, source_main_fold$W_outcome, M_tau_inference,
    family_int, link_int, A_val
  )
  delta_ts_k1 <- correction_result$delta_ts
  
  # Calculate outcome model component using target main fold
  target_main_fold <- materialize_fold(target_folds, k1)
  mu_pred_ts_k1 <- mean(predict_glm_cpp(target_main_fold$W_outcome, alpha_final_k1,
                                    family_int, link_int))
  
  # Final source estimate: μ̂_{t,s_j} = mu_pred_ts + δ_{t,s_j}
  mu_ts_k1 <- mu_pred_ts_k1 + delta_ts_k1
  
  return(list(
    mu_ts = mu_ts_k1,
    gamma_s = gamma_final_k1,
    alpha_ts = alpha_final_k1,
    delta_ts = delta_ts_k1,
    mu_pred_ts = mu_pred_ts_k1,
    n_s = source_main_fold$n,
    correction_components = correction_result$correction_components,
    correction_clip_diagnostics = correction_result$clip_diagnostics,
    n_calibrated_folds = length(source_calib_list),
    calibrated_fold_keys = names(source_calib_list),
    # Per-k2 out-of-two-fold plug-ins for aggregation inner validation.
    # These are trained without folds k1 and k2, so validation summaries on k2
    # do not reuse the validation observations.
    per_k2_gamma = gamma_init_list,
    per_k2_alpha = alpha_init_list
  ))
}

#' Partition data into K_f folds with stratified sampling
#'
#' Uses stratified sampling by treatment indicator A to ensure each fold
#' contains both treated and control units. This prevents degenerate folds
#' which would cause propensity score estimation to fail.
#'
#' Returns lightweight fold views storing only indices and a reference to
#' the original data. Data is materialized on demand via \code{materialize_fold},
#' avoiding upfront copies of W_outcome, Z_site, A, Y matrices.
#'
#' @param data Site data (list with W_outcome, Z_site, A, Y, n)
#' @param n_folds Number of folds
#' @param seed Optional RNG seed for reproducible fold assignment. When non-NULL,
#'   the RNG state is saved, \code{set.seed(seed)} is called, fold IDs are
#'   assigned, and the original RNG state is restored. This ensures fold
#'   creation is deterministic regardless of prior random draws.
#' @return List of fold views (each with original_idx, n), with .data_ref attribute
partition_into_folds <- function(data, n_folds, seed = NULL) {
  n <- data$n
  A <- data$A
  
  # Fold assignment logic (may be wrapped in with_seed for reproducibility)
  assign_folds <- function() {
    fold_ids <- integer(n)
    treated_idx <- which(A == 1)
    control_idx <- which(A == 0)
    
    if (length(treated_idx) >= n_folds && length(control_idx) >= n_folds) {
      fold_ids[treated_idx] <- sample(rep(1:n_folds, length.out = length(treated_idx)))
      fold_ids[control_idx] <- sample(rep(1:n_folds, length.out = length(control_idx)))
    } else {
      stop(sprintf(
        "partition_into_folds: too few treated or control units for stratified fold assignment (n_treated=%d, n_control=%d, n_folds=%d). Reduce n_folds or provide data with at least one observation from each treatment arm in every fold.",
        length(treated_idx), length(control_idx), n_folds
      ), call. = FALSE)
    }
    fold_ids
  }
  
  # Isolate RNG state if a seed is requested (reproducible fold assignment)
  fold_ids <- if (!is.null(seed)) with_seed(seed, assign_folds()) else assign_folds()
  
  # Store only indices per fold; data is subset on demand via materialize_fold()
  # The .data_ref attribute on the returned list holds the parent data (zero-copy
  # under R's copy-on-modify semantics).
  folds <- lapply(1:n_folds, function(k) {
    idx <- which(fold_ids == k)
    list(
      original_idx = idx,
      n = length(idx)
    )
  })
  attr(folds, ".data_ref") <- data
  folds
}

#' Materialize a fold view into a concrete dataset
#'
#' Subsets the parent data by the fold's index vector. Call this when
#' W_outcome, Z_site, A, Y fields are actually needed (e.g., before passing
#' to model fitting or C++ functions). Fields \code{n} and \code{original_idx}
#' are available directly on the fold view without materialization.
#'
#' @param fold_list The fold list returned by \code{partition_into_folds}
#' @param k Index of the fold to materialize
#' @return List with W_outcome, Z_site, A, Y, n, original_idx
materialize_fold <- function(fold_list, k) {
  data <- attr(fold_list, ".data_ref")
  if (is.null(data)) {
    stop("materialize_fold() requires fold lists created by partition_into_folds() with .data_ref.")
  }
  fold <- fold_list[[k]]
  idx <- fold$original_idx
  list(
    W_outcome = data$W_outcome[idx, , drop = FALSE],
    Z_site = data$Z_site[idx, , drop = FALSE],
    A = data$A[idx],
    Y = data$Y[idx],
    n = fold$n,
    original_idx = idx
  )
}

#' Extract source estimates matrix from fold results (vectorized)
#'
#' @param fold_results List of fold results
#' @param source_sites Vector of source site names
#' @param n_folds Number of folds
#' @return Matrix of source estimates (n_folds x K)
extract_source_estimates_matrix <- function(fold_results, source_sites, n_folds) {
  K <- length(source_sites)
  
  # Build n_folds x K matrix using vapply for type safety and consistent dimensions
  estimates <- vapply(seq_len(n_folds), function(k) {
    vapply(source_sites, function(s) {
      fold_results[[k]]$source_results[[s]]$mu_ts
    }, FUN.VALUE = numeric(1))
  }, FUN.VALUE = numeric(K))
  
  # vapply always returns K x n_folds matrix (even when K=1, it's 1 x n_folds)
  # Transpose to n_folds x K
  t(estimates)
}

# NOTE: process_source_site_one_round has been merged into the unified
# process_source_site() function above. The algorithm-specific data access
# pattern is now handled via the get_fold_inputs callback parameter.

#' Build deterministic fold partitions for cross-fitting reuse
#'
#' Creates target/source fold objects once so multiple calls (e.g., one-round,
#' two-round, and ATE reruns with A_val=0) can reuse identical partitions.
#'
#' @param data_split split data by site
#' @param n_folds number of cross-fitting folds
#' @return list(target_folds=..., source_folds=...)
#' @keywords internal
build_crossfit_folds <- function(data_split, n_folds) {
  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  fold_seed_base <- target_data$n + sum(sapply(source_sites, function(s) data_split[[s]]$n))

  target_folds <- partition_into_folds(target_data, n_folds, seed = fold_seed_base)
  source_folds <- setNames(
    lapply(seq_along(source_sites), function(i) {
      partition_into_folds(data_split[[source_sites[i]]], n_folds,
                           seed = fold_seed_base + i)
    }),
    source_sites
  )

  list(target_folds = target_folds, source_folds = source_folds)
}


# =============================================================================
# UNIFIED CROSS-FITTING ALGORITHM
# =============================================================================
# Implements both one-round and two-round communication protocols.
# The key algorithmic difference is WHO trains the initial outcome model:
#   - Two-round: each SOURCE site trains on its own data (Algorithm 1)
#   - One-round: the TARGET site trains on its local data (Algorithm 2)
# All other steps (fold partitioning, calibration, aggregation) are identical.
# =============================================================================

#' Unified Cross-fitting Algorithm for FACE-HD
#'
#' Implements Algorithms 1 and 2 from main.tex (Section A.2).
#' The \code{communication_mode} parameter controls who trains the initial
#' outcome model:
#' \itemize{
#'   \item \code{"two_round"}: Source sites train initial outcome models (Round 1),
#'     target computes summaries from source-specific alphas (Round 2).
#'   \item \code{"one_round"}: Target site trains initial outcome models locally
#'     and sends both alphas and summaries to source sites (single round).
#' }
#'
#' @param data_split split data by site (from split_data_by_site)
#' @param n_folds number of cross-fitting folds (default 10, K_f in main.tex).
#'        Minimum 3 required for proper two-level cross-fitting calibration.
#' @param communication_mode character, "two_round" or "one_round"
#' @param lambda_selection method for lambda selection ("cv" or numeric value)
#' @param lambda_rule rule for the aggregation penalty when
#'        \code{lambda_selection = "cv"}: \code{"min"} (default) selects the
#'        bare validation minimizer. \code{"1se"} selects the largest lambda
#'        whose average inner-validation criterion is within one standard error
#'        of the minimum.
#' @param verbose print progress messages
#' @param M_tau truncation parameter for calibrated losses during training (default 10.0)
#' @param M_tau_inference truncation parameter for the inference step (correction term
#'        and variance calculation). Default Inf (no truncation) per the paper:
#'        T(.) appears only in calibrated losses, not in the final DR estimator.
#' @param n_cores number of cores for parallel processing of source sites.
#'        NULL or 1 for sequential, -1 for all cores minus 1.
#' @param nlambda_init Integer. Number of lambda values for cv.glmnet in the
#'        initial outcome model (fit_initial_outcome). Default 100.
#' @param nuisance_lambda_rule CV selection rule for nuisance fits:
#'        \code{"min"} (default) uses \code{lambda.min}; \code{"1se"} uses
#'        \code{lambda.1se}.
#' @param family GLM family specification: "gaussian" or "binomial".
#' @param A_val Integer (0 or 1). Treatment value for potential outcome estimation.
#' @param use_lambda_cache Logical. If TRUE, reuse selected nuisance
#'        lambdas within fold-specific loops to reduce repeated CV.
#' @param precomputed_folds Optional list with \\code{target_folds} and
#'        \\code{source_folds} created in advance for this \\code{data_split}.
#' @param target_only_ps_cache Optional environment used to cache target-only
#'        complement-fold propensity predictions.
#' @return list with estimates and variance components
#' @export
run_crossfit <- function(data_split, n_folds = N_FOLDS_DEFAULT,
                         communication_mode = c("two_round", "one_round"),
                         lambda_selection = "cv",
                         lambda_rule = c("min", "1se"),
                         verbose = TRUE, M_tau = M_TAU_DEFAULT,
                         M_tau_inference = M_TAU_INFERENCE_DEFAULT,
                         n_cores = NULL,
                         nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
                         family = "binomial", A_val = 1L,
                         use_lambda_cache = TRUE,
                         precomputed_folds = NULL,
                         target_only_ps_cache = NULL,
                         nuisance_lambda_rule = c("min", "1se")) {

  communication_mode <- match.arg(communication_mode)
  lambda_rule <- match.arg(lambda_rule)
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(nuisance_lambda_rule, "run_crossfit")

  # ---- SHARED SETUP ----
  glm_spec <- resolve_glm_family(family)
  link <- glm_spec$link
  family_int <- glm_spec$family_int
  link_int <- glm_spec$link_int

  validate_algorithm_inputs(data_split, lambda_selection, family = family, A_val = A_val)
  validate_truncation_parameters(M_tau, M_tau_inference)

  if (n_folds < 3) {
    stop(sprintf("n_folds=%d is too small for two-level cross-fitting (minimum: 3, recommended: >= 5).", n_folds))
  }
  if (n_folds < 5 && verbose) {
    warning("n_folds < 5 may lead to high variance. Consider using n_folds >= 5 for better stability.")
  }

  validate_crossfit_sample_sizes(data_split, n_folds, A_val = A_val)

  A_val <- as.integer(A_val)
  stopifnot("A_val must be 0 or 1" = A_val %in% c(0L, 1L))

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  K <- length(source_sites)

  actual_cores <- setup_parallel(n_cores)
  use_parallel <- actual_cores > 1

  mode_label <- toupper(gsub("_", "-", communication_mode))

  if (verbose) {
    cat(paste0("\n", paste(rep("=", 70), collapse = ""), "\n"))
    cat(sprintf("==> STARTING %s TWO-LAYER CROSS-FITTING ALGORITHM\n", mode_label))
    cat(paste0("  Total sites: ", K + 1, " (1 target + ", K, " source)\n"))
    cat(paste0("  Target site sample size: ", target_data$n, "\n"))
    cat(paste0("  GLM family: ", family, " (link: ", link, ")\n"))
    cat(paste0("  Cross-fitting folds (K_f): ", n_folds, "\n"))
    cat(paste0("  Truncation parameter (M_tau): ", M_tau, "\n"))
    cat(paste0("  Treatment value (A_val): ", A_val, " (estimating mu^", A_val, "_t)\n"))
    if (use_parallel) {
      cat(paste0("  Parallel processing: ", actual_cores, " cores\n"))
    }
    cat(paste0(paste(rep("=", 70), collapse = ""), "\n\n"))
  }

  # ---- PARTITION DATA INTO FOLDS ----
  if (verbose) {
    cat("--> INITIALIZATION: PARTITIONING DATA INTO FOLDS\n")
  }

  if (!is.null(precomputed_folds)) {
    if (is.null(precomputed_folds$target_folds) || is.null(precomputed_folds$source_folds)) {
      stop("run_crossfit: precomputed_folds must contain 'target_folds' and 'source_folds'.")
    }
    target_folds <- precomputed_folds$target_folds
    source_folds <- precomputed_folds$source_folds
  } else {
    fold_seed_base <- target_data$n + sum(sapply(source_sites, function(s) data_split[[s]]$n))
    target_folds <- partition_into_folds(target_data, n_folds, seed = fold_seed_base)
    source_folds <- setNames(
      lapply(seq_along(source_sites), function(i) {
        partition_into_folds(data_split[[source_sites[i]]], n_folds,
                             seed = fold_seed_base + i)
      }),
      source_sites
    )
  }

  fold_results <- list()

  # ---- MAIN LOOP ----
  for (k1 in 1:n_folds) {
    if (verbose) {
      cat(paste0("\n--> PROCESSING MAIN FOLD k1 = ", k1, "/", n_folds, "\n"))
    }

    # Reset lambda caches PER k1 fold to avoid systematic bias.
    # Each fold has a different train/test split, so the optimal lambda
    # can differ. Caching across k1 folds introduces a non-vanishing bias
    # because the first fold's lambda may be suboptimal for later folds,
    # and this error doesn't average out (same direction of bias).
    # Within a k1 fold, caching is restricted to local loops:
    # - two_round: per-source-site across k2
    # - one_round: target-side across k2
    cached_lambda <- NULL
    site_cached_lambdas <- if (isTRUE(use_lambda_cache) && communication_mode == "two_round") {
      setNames(vector("list", length(source_sites)), source_sites)
    } else {
      NULL
    }

    secondary_folds <- setdiff(1:n_folds, k1)

    # ================================================================
    # COMMUNICATION-MODE SPECIFIC: prepare initial models & summaries
    # ================================================================
    if (communication_mode == "two_round") {
      # ROUND 1: Source sites compute initial outcome models
      log_info(verbose, "   --> ROUND 1: INITIAL NUISANCE MODELS (SOURCE)")

      source_initial_params <- list()
      for (s in source_sites) {
        source_initial_params[[s]] <- list()
        cached_lambda_site <- if (isTRUE(use_lambda_cache)) site_cached_lambdas[[s]] else NULL
        for (k2 in secondary_folds) {
          training_folds <- setdiff(1:n_folds, c(k1, k2))
          source_train <- combine_folds(source_folds[[s]], training_folds)

          if (source_train$n == 0) {
            stop(sprintf(
              "Two-round: source site '%s' training folds {%s} have 0 observations for k2=%d; cannot fit the initial outcome nuisance required by main.tex eq:alpha_init.",
              s, paste(training_folds, collapse = ","), k2
            ))
          }

          alpha_init_k1_k2 <- fit_initial_outcome(
            source_train$W_outcome, source_train$Y, source_train$A, A_val,
            lambda = if (isTRUE(use_lambda_cache)) cached_lambda_site else NULL,
            nlambda = nlambda_init, family = family,
            lambda_rule = nuisance_lambda_rule
          )
          if (isTRUE(use_lambda_cache) && is.null(cached_lambda_site)) {
            cached_lambda_site <- attr(alpha_init_k1_k2, "lambda_used")
            if (verbose && !is.null(cached_lambda_site)) {
              cat(sprintf("      Cached lambda (site %s) = %.6f\n", s, cached_lambda_site))
            }
          }
          source_initial_params[[s]][[paste0("k2_", k2)]] <- alpha_init_k1_k2
        }
        if (isTRUE(use_lambda_cache)) {
          site_cached_lambdas[[s]] <- cached_lambda_site
        }
      }

      # Target precomputes mean_phi (shared across sources) and
      # mean_grad_psi_init (per source, depends on source alpha_init)
      target_fold_cache <- list()
      target_train_cache <- list()
      mean_phi_cache <- list()
      for (k2 in secondary_folds) {
        k2_key <- paste0("k2_", k2)
        training_folds <- setdiff(1:n_folds, c(k1, k2))
        target_fold_cache[[k2_key]] <- materialize_fold(target_folds, k2)
        target_train_cache[[k2_key]] <- combine_folds(target_folds, training_folds)
        if (target_train_cache[[k2_key]]$n > 0) {
          mean_phi_cache[[k2_key]] <- c(1, colMeans(target_train_cache[[k2_key]]$Z_site))
        } else {
          mean_phi_cache[[k2_key]] <- rep(0, ncol(target_fold_cache[[k2_key]]$Z_site) + 1)
        }
      }

      target_summaries <- list()
      for (s in source_sites) {
        target_summaries[[s]] <- list()
        for (k2 in secondary_folds) {
          k2_key <- paste0("k2_", k2)
          alpha_init_k1_k2 <- source_initial_params[[s]][[k2_key]]
          target_fold_k2 <- target_fold_cache[[k2_key]]

          if (target_fold_k2$n == 0) {
            stop(sprintf("Two-round: target fold k2=%d has 0 observations for source site '%s'.",
                         k2, s))
          }

          mean_grad_psi_init <- .mean_glm_gradient_site_basis(
            target_fold_k2$W_outcome, target_fold_k2$Z_site,
            alpha_init_k1_k2, family_int, link_int
          )
          target_summaries[[s]][[k2_key]] <- list(
            mean_grad_psi_init = mean_grad_psi_init,
            mean_phi = mean_phi_cache[[k2_key]]
          )
        }
      }

      get_fold_inputs <- function(site, k2) {
        k2_key <- paste0("k2_", k2)
        list(
          mean_phi = target_summaries[[site]][[k2_key]]$mean_phi,
          mean_grad_psi_init = target_summaries[[site]][[k2_key]]$mean_grad_psi_init,
          alpha_init = source_initial_params[[site]][[k2_key]]
        )
      }

      log_info(verbose, "   --> ROUND 2: CALIBRATED NUISANCE MODELS")

    } else {
      # SINGLE ROUND: Target computes initial models & summaries
      log_info(verbose, "   --> SINGLE ROUND: TARGET COMPUTES AND SENDS ALL INFO")

      target_initial_models <- list()
      target_summaries <- list()

      for (k2 in secondary_folds) {
        k2_key <- paste0("k2_", k2)
        training_folds <- setdiff(1:n_folds, c(k1, k2))
        target_train <- combine_folds(target_folds, training_folds)
        target_calib_k2 <- materialize_fold(target_folds, k2)

        if (target_train$n == 0 || target_calib_k2$n == 0) {
          stop(sprintf("One-round: target train (n=%d) or calibration fold k2=%d (n=%d) is empty.",
                       target_train$n, k2, target_calib_k2$n))
        }

        alpha_init_k1_k2 <- fit_initial_outcome(
          target_train$W_outcome, target_train$Y, target_train$A, A_val,
          lambda = if (isTRUE(use_lambda_cache)) cached_lambda else NULL,
          nlambda = nlambda_init, family = family,
          lambda_rule = nuisance_lambda_rule
        )
        if (isTRUE(use_lambda_cache) && is.null(cached_lambda)) {
          cached_lambda <- attr(alpha_init_k1_k2, "lambda_used")
          if (verbose && !is.null(cached_lambda)) {
            cat(sprintf("      Cached lambda = %.6f (reused for remaining calls)\n", cached_lambda))
          }
        }

        target_initial_models[[k2_key]] <- alpha_init_k1_k2

        mean_grad_psi_init <- .mean_glm_gradient_site_basis(
          target_calib_k2$W_outcome, target_calib_k2$Z_site,
          alpha_init_k1_k2, family_int, link_int
        )
        mean_phi <- c(1, colMeans(target_train$Z_site))

        target_summaries[[k2_key]] <- list(
          mean_grad_psi_init = mean_grad_psi_init,
          mean_phi = mean_phi
        )
      }

      get_fold_inputs <- function(site, k2) {
        k2_key <- paste0("k2_", k2)
        list(
          mean_phi = target_summaries[[k2_key]]$mean_phi,
          mean_grad_psi_init = target_summaries[[k2_key]]$mean_grad_psi_init,
          alpha_init = target_initial_models[[k2_key]]
        )
      }
    }

    # ================================================================
    # SHARED: Process source sites via unified helper
    # ================================================================
    if (verbose) {
      if (use_parallel) cat(sprintf("      Processing source sites (parallel: %d cores)\n", actual_cores))
    }

    # combine_cache stores combine_folds() results keyed by (site, training-fold-set).
    # Parallelism note: in sequential mode every k2-iteration for a site can read
    # entries written by earlier k2-iterations of the SAME site.
    # In mclapply (fork) mode, each worker child inherits the environment at fork
    # time (copy-on-write), so entries created by one worker are NOT visible to
    # sibling workers or the parent.  The cache still prevents redundant computation
    # within the sequential k2-loop *inside* each worker, but cross-worker sharing
    # is not possible under fork-based parallelism without explicit IPC.
    combine_cache <- new.env(hash = TRUE, parent = emptyenv())
    source_results_list <- parallel_lapply(
      source_sites,
      function(s) {
        process_source_site(
          s = s,
          source_folds = source_folds,
          target_folds = target_folds,
          k1 = k1,
          n_folds = n_folds,
          A_val = A_val,
          M_tau = M_tau,
          M_tau_inference = M_tau_inference,
          data_split = data_split,
          combine_cache = combine_cache,
          get_fold_inputs = get_fold_inputs,
          family_int = family_int,
          link_int = link_int,
          use_lambda_cache = use_lambda_cache,
          nuisance_lambda_rule = nuisance_lambda_rule
        )
      },
      n_cores = actual_cores
    )

    source_results_k1 <- setNames(source_results_list, source_sites)

    if (verbose) {
      cat(sprintf("      [OK] Processed %d source sites\n", length(source_sites)))
    }

    # Target-only estimate (complement-fold training, avoids nested cross-fitting)
    target_only_k1 <- estimate_target_only_from_complement(target_folds, k1, n_folds,
                                                           family = family, A_val = A_val,
                                                           propensity_cache = target_only_ps_cache)

    # Inner-fold target-only estimates (Version A Step 2)
    # For each k2 in secondary_folds, train target-only model on
    # {folds}\{k1,k2} and evaluate on fold k2. These are used to compute
    # varphi_ot on inner folds for variance component estimation.
    log_info(verbose, "   --> Inner-fold target-only models for Version A Step 2")
    target_only_inner <- list()
    for (k2 in secondary_folds) {
      target_only_inner[[paste0("k2_", k2)]] <-
        estimate_target_only_from_complement(target_folds, k1, n_folds,
                                              k2 = k2, family = family, A_val = A_val,
                                              propensity_cache = target_only_ps_cache)
    }

    fold_results[[k1]] <- list(
      target_only = target_only_k1,
      source_results = source_results_k1,
      target_only_inner = target_only_inner
    )
  }

  # ---- FINAL AGGREGATION ----
  if (verbose) {
    cat(sprintf("\n--> FINAL AGGREGATION (%s)\n", mode_label))
  }

  target_estimates <- sapply(fold_results, function(x) x$target_only$estimate)
  final_target_estimate <- mean(target_estimates)

  source_estimates_matrix <- extract_source_estimates_matrix(fold_results, source_sites, n_folds)
  source_estimates <- colMeans(source_estimates_matrix)
  names(source_estimates) <- source_sites

  result <- calculate_crossfit_aggregation(
    data_split = data_split, target_data = target_data,
    source_sites = source_sites, K = K,
    target_folds = target_folds, source_folds = source_folds,
    fold_results = fold_results,
    n_folds = n_folds, M_tau = M_tau,
    M_tau_inference = M_tau_inference,
    lambda_selection = lambda_selection, verbose = verbose,
    lambda_rule = lambda_rule,
    final_target_estimate = final_target_estimate,
    target_estimates = target_estimates,
    source_estimates = source_estimates,
    source_estimates_matrix = source_estimates_matrix,
    crossfit_type = communication_mode,
    algorithm_label = paste0(communication_mode, "_two_layer_crossfit"),
    family_int = family_int,
    link_int = link_int,
    A_val = A_val
  )
  result$nuisance_lambda_rule <- nuisance_lambda_rule
  return(result)
}

#' Combine multiple folds into a single dataset
#'
#' Merges the data from the specified fold indices into a unified list
#' with the same structure as a single fold (W_outcome, Z_site, A, Y, n).
#' Uses the \code{.data_ref} attribute for efficient index-based subsetting
#' rather than repeated \code{rbind} operations.
#'
#' @param fold_list List of fold views returned by \code{partition_into_folds}.
#' @param fold_indices Integer vector of fold indices to combine.
#' @return List with combined data (W_outcome, Z_site, A, Y, n, original_idx).
combine_folds <- function(fold_list, fold_indices) {
  data <- attr(fold_list, ".data_ref")
  if (is.null(data)) {
    stop("combine_folds() requires fold lists created by partition_into_folds() with .data_ref.")
  }
  
  if (length(fold_indices) == 0) {
    # Return empty dataset with correct structure
    return(list(
      W_outcome = matrix(0, nrow = 0, ncol = ncol(data$W_outcome)),
      Z_site = matrix(0, nrow = 0, ncol = ncol(data$Z_site)),
      A = numeric(0), Y = numeric(0), n = 0, original_idx = integer(0)
    ))
  }
  
  if (length(fold_indices) == 1) {
    return(materialize_fold(fold_list, fold_indices[1]))
  }
  
  # Combine indices from all requested folds
  combined_idx <- do.call(c, lapply(fold_indices, function(i) fold_list[[i]]$original_idx))
  
  # Efficient path: single index-based subset of original data
  # (avoids multiple rbind operations on pre-copied fold subsets)
  return(list(
    W_outcome = data$W_outcome[combined_idx, , drop = FALSE],
    Z_site = data$Z_site[combined_idx, , drop = FALSE],
    A = data$A[combined_idx],
    Y = data$Y[combined_idx],
    n = length(combined_idx),
    original_idx = combined_idx
  ))
}

#' Legacy aggregation-penalty (lambda) variance-curve selector
#'
#' Selects the aggregation penalty lambda for the weight optimization
#' (eq:final_opt in main.tex) over a glmnet-style KKT lambda path. This is NOT a
#' cross-validation procedure: lambda affects only the low-dimensional weight
#' optimization (not nuisance estimation), so there are no held-out fold losses.
#' The current cross-fitting path uses \code{select_aggregation_lambda_inner_cv()}
#' in \code{R/cross_fitting_aggregation.R}. This helper remains for diagnostic
#' variance-curve evaluation and only supports \code{lambda_rule = "min"}.
#'
#' \enumerate{
#'   \item \strong{KKT lambda path.}
#'     In eq:final_opt, the penalty is weighted L1:
#'     \eqn{\lambda |\eta_j| (\mu_{ot} - \mu_{ts,j})^2}.
#'     At \eqn{\eta = 0}, the smooth variance gradient is
#'     \eqn{2(C_{ot,j} - V_{ot})/n_t}, so the zero-solution endpoint is
#'     \eqn{\lambda_{\max} = \max_j |2(C_{ot,j} - V_{ot})/n_t| /
#'     (\mu_{ot} - \mu_{ts,j})^2} over coordinates with positive penalty
#'     weights. The grid descends geometrically from this endpoint.
#'
#'   \item \strong{Selection rule.}
#'     This legacy helper evaluates one variance curve, not fold-level validation
#'     scores, so it cannot compute a valid one-SE rule. Use
#'     \code{select_aggregation_lambda_inner_cv()} through
#'     \code{calculate_crossfit_aggregation()} for glmnet-style
#'     \code{"min"}/\code{"1se"} aggregation tuning.
#' }
#'
#' @param n_folds number of cross-fitting folds used in main algorithm
#' @param fold_target_estimate target-only estimate for current fold
#' @param fold_source_estimates source estimates for current fold
#' @param variances_k1 list of variance components (V_ot, V_t, V_s)
#' @param C_ot_k1 target-source covariance vector
#' @param n_samples_k1 list of sample sizes (n_t, n_s)
#' @param C_cross_k1 cross-site covariance matrix
#' @param verbose print progress
#' @param crossfit_type type of cross-fitting algorithm ("one_round" or "two_round")
#' @param lambda_grid lambda values to test (default NULL for KKT lambda path)
#' @param lambda_rule selection rule. Only \code{"min"} is supported here because
#'   this legacy helper has no fold-level standard errors.
#' @return selected lambda value
select_aggregation_lambda <- function(n_folds, fold_target_estimate,
                                     fold_source_estimates, variances_k1, C_ot_k1,
                                     n_samples_k1, C_cross_k1, verbose = FALSE,
                                     crossfit_type = "one_round",
                                     lambda_grid = NULL,
                                     lambda_rule = c("min", "1se")) {

  lambda_rule <- match.arg(lambda_rule)
  if (identical(lambda_rule, "1se")) {
    stop("select_aggregation_lambda: lambda_rule='1se' requires fold-level validation standard errors; use select_aggregation_lambda_inner_cv() via calculate_crossfit_aggregation(), or set lambda_rule='min' for this legacy variance-curve selector.",
         call. = FALSE)
  }
  K <- length(fold_source_estimates)

  # ==========================================================================
  # STEP 1: glmnet-style lambda grid from KKT lambda_max
  # ==========================================================================
  # In eq:final_opt, the penalty for site j is weighted L1:
  #   λ * |η_j| * (μ_ot − μ_{ts,j})²
  # The variance curvature contribution scales as:
  #   η_j² * (V_{t,j}/n_t + V_{s,j}/n_{s,j})
  #
  # At η = 0, the smooth variance gradient is
  #   ∂V/∂η_j = 2 * (C_ot,j - V_ot) / n_t.
  # The exact zero-solution KKT endpoint is therefore
  #   λ_max = max_j |∂V/∂η_j| / d_j²
  # over coordinates with d_j² > 0. The helper below implements this and then
  # builds a descending geometric path, matching the nuisance CV convention.
  # ==========================================================================
  adaptive_grid <- is.null(lambda_grid)
  lambda_max <- NA_real_

  if (adaptive_grid) {
    grid_component <- list(
      V_ot = variances_k1$V_ot,
      C_ot = C_ot_k1,
      avg_target_est = fold_target_estimate,
      avg_source_est = fold_source_estimates,
      n_t = n_samples_k1$n_t
    )
    lambda_grid <- .aggregation_lambda_grid(grid_component)
    lambda_max <- attr(lambda_grid, "lambda_max")
  } else {
    if (!is.numeric(lambda_grid) || length(lambda_grid) == 0L ||
        any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
      stop("select_aggregation_lambda: lambda_grid must contain positive finite numeric values.",
           call. = FALSE)
    }
    clipped <- pmax(pmin(lambda_grid, LAMBDA_MAX), LAMBDA_MIN)
    if (!isTRUE(all.equal(clipped, lambda_grid, check.attributes = FALSE))) {
      warning(sprintf("select_aggregation_lambda: lambda_grid values were clipped to [%g, %g].",
                      LAMBDA_MIN, LAMBDA_MAX), call. = FALSE)
    }
    lambda_grid <- clipped
  }

  if (verbose) {
    cat("Weight-optimization lambda selection:\n")
    if (adaptive_grid) {
      cat(sprintf("  KKT lambda_max: %.4f  |  Grid: [%.2e, %.2e]  |  %d pts\n",
                  lambda_max, min(lambda_grid), max(lambda_grid),
                  length(lambda_grid)))
    } else {
      cat(sprintf("  User-supplied grid: [%.2e, %.2e]  |  %d pts\n",
                  min(lambda_grid), max(lambda_grid), length(lambda_grid)))
    }
  }

  # ==========================================================================
  # STEP 2: Evaluate aggregated variance for each lambda
  # ==========================================================================
  # Batch C++ call eliminates per-lambda R→C++ round-trips (previously ~2×100
  # transitions per fold for the default 100-point grid).  The same input
  # clipping and matrix preparation from optimize_weights / model_fitting.R is
  # applied here so results are numerically identical to the sequential loop.
  # ==========================================================================
  cross_matrix <- if (is.matrix(C_cross_k1) && nrow(C_cross_k1) == K &&
                        ncol(C_cross_k1) == K) C_cross_k1 else matrix(0, 0, 0)

  cv_variances <- tryCatch(
    evaluate_lambda_grid_weights_cpp(
      estimates   = pmax(pmin(fold_source_estimates, ESTIMATE_MAX), -ESTIMATE_MAX),
      V_t         = pmax(variances_k1$V_t, VARIANCE_MIN),
      V_s         = pmax(variances_k1$V_s, VARIANCE_MIN),
      n_s_vec     = pmax(n_samples_k1$n_s, 1),
      V_ot        = max(variances_k1$V_ot, VARIANCE_MIN),
      n_t         = n_samples_k1$n_t,
      C_ot        = pmax(pmin(C_ot_k1, ESTIMATE_MAX), -ESTIMATE_MAX),
      mu_ot       = max(min(fold_target_estimate, ESTIMATE_MAX), -ESTIMATE_MAX),
      C_cross     = cross_matrix,
      lambda_grid = lambda_grid,
      max_iter    = WEIGHT_OPT_MAX_ITER,
      tol         = WEIGHT_OPT_TOL
    ),
    error = function(e) {
      stop(sprintf(
        "select_aggregation_lambda: batch C++ evaluation failed: %s",
        conditionMessage(e)
      ))
    }
  )

  valid <- which(is.finite(cv_variances) & cv_variances > 0)

  if (length(valid) == 0) {
    stop(sprintf("select_aggregation_lambda: all %d lambdas produced invalid variance; cannot select aggregation lambda.",
                 length(lambda_grid)))
  }

  min_var <- min(cv_variances[valid])
  min_var_idx <- valid[which.min(cv_variances[valid])]

  best_idx <- min_var_idx
  best_lambda <- lambda_grid[best_idx]

  if (verbose) {
    cat(sprintf("  Min variance: %.6f at lambda=%.4f\n",
                min_var, lambda_grid[min_var_idx]))
    cat(sprintf("  Selected:     lambda=%.4f (var=%.6f, +%.1f%% vs min, rule=%s)\n",
                best_lambda, cv_variances[best_idx],
                100 * (cv_variances[best_idx] / min_var - 1), lambda_rule))
  }

  return(best_lambda)
}
