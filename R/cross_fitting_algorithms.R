# cross_fitting_algorithms.R - Two-Layer Cross-fitting variants for FACE algorithm
# Following docs/main.tex (Two-Level Cross-fitting appendix section)
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# run_crossfit() provides the arm-specific building block
# mu^a_t = E_t[Y(a)]. run_tate_crossfit() fits both arms on a common partition,
# contrasts their influence components, and estimates one TATE-level source-
# weight vector for tau_t = mu^1_t - mu^0_t.
#
# =============================================================================
# This file implements the main RoCE algorithms with optional parallelization
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

.elapsed_process_seconds <- function(started_at) {
  as.numeric(proc.time()[["elapsed"]] - started_at)
}

.max_iterations_or_na <- function(iterations) {
  iterations <- suppressWarnings(as.integer(iterations))
  iterations <- iterations[!is.na(iterations) & iterations >= 0L]
  if (length(iterations) == 0L) NA_integer_ else max(iterations)
}

.max_numeric_or_na <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values <- values[is.finite(values)]
  if (length(values) == 0L) NA_real_ else max(values)
}

.safe_density_ratio_warm_start <- function(fit) {
  fit_numeric <- suppressWarnings(as.numeric(fit))
  line_search_failures <- suppressWarnings(
    as.integer(attr(fit, "line_search_failures"))
  )
  if (!isTRUE(attr(fit, "converged")) ||
      length(fit_numeric) == 0L || any(!is.finite(fit_numeric)) ||
      max(abs(fit_numeric)) >= 0.9 * PARAM_MAX ||
      length(line_search_failures) != 1L || is.na(line_search_failures) ||
      line_search_failures > 0L) {
    return(NULL)
  }
  fit_numeric
}

.summarize_source_nuisance_fit_diagnostics <- function(source_results) {
  diagnostics <- lapply(
    source_results,
    function(site_result) site_result$nuisance_fit_diagnostics
  )
  if (length(diagnostics) == 0L || any(vapply(diagnostics, is.null, logical(1L)))) {
    stop(
      ".summarize_source_nuisance_fit_diagnostics: every source result must contain diagnostics.",
      call. = FALSE
    )
  }

  initial_iterations <- unlist(lapply(
    diagnostics,
    `[[`, "initial_density_ratio_iterations"
  ), use.names = FALSE)
  calibrated_dr_iterations <- vapply(
    diagnostics,
    function(diagnostic) diagnostic$calibrated_density_ratio_iterations,
    integer(1L)
  )
  calibrated_outcome_iterations <- vapply(
    diagnostics,
    function(diagnostic) diagnostic$calibrated_outcome_iterations,
    integer(1L)
  )
  initial_update_ratios <- unlist(lapply(
    diagnostics,
    `[[`, "initial_density_ratio_update_ratios"
  ), use.names = FALSE)
  initial_abs_coefficients <- unlist(lapply(
    diagnostics,
    `[[`, "initial_density_ratio_max_abs_coefficients"
  ), use.names = FALSE)
  calibrated_dr_update_ratios <- vapply(
    diagnostics,
    function(diagnostic) as.numeric(
      diagnostic$calibrated_density_ratio_update_ratio %||% NA_real_
    ),
    numeric(1L)
  )
  calibrated_dr_abs_coefficients <- vapply(
    diagnostics,
    function(diagnostic) as.numeric(
      diagnostic$calibrated_density_ratio_max_abs_coefficient %||% NA_real_
    ),
    numeric(1L)
  )
  calibrated_outcome_update_ratios <- vapply(
    diagnostics,
    function(diagnostic) as.numeric(
      diagnostic$calibrated_outcome_update_ratio %||% NA_real_
    ),
    numeric(1L)
  )
  calibrated_outcome_abs_coefficients <- vapply(
    diagnostics,
    function(diagnostic) as.numeric(
      diagnostic$calibrated_outcome_max_abs_coefficient %||% NA_real_
    ),
    numeric(1L)
  )
  initial_support_floors <- unlist(lapply(
    diagnostics,
    `[[`, "initial_density_ratio_support_floors"
  ), use.names = FALSE)
  calibrated_support_floors <- vapply(
    diagnostics,
    function(diagnostic) as.numeric(
      diagnostic$calibrated_density_ratio_support_floor %||% NA_real_
    ),
    numeric(1L)
  )
  sum_site_diagnostic <- function(field) {
    sum(vapply(
      diagnostics,
      function(diagnostic) {
        values <- suppressWarnings(as.numeric(diagnostic[[field]]))
        if (length(values) == 0L || any(!is.finite(values)) ||
            any(values < 0)) {
          stop(sprintf(
            ".summarize_source_nuisance_fit_diagnostics: invalid '%s'.",
            field
          ), call. = FALSE)
        }
        sum(values)
      },
      numeric(1L)
    ))
  }

  c(
    initial_dr_nonconverged = sum(vapply(
      diagnostics,
      function(diagnostic) {
        sum(!diagnostic$initial_density_ratio_converged)
      },
      integer(1L)
    )),
    calibrated_dr_nonconverged = sum(!vapply(
      diagnostics,
      function(diagnostic) diagnostic$calibrated_density_ratio_converged,
      logical(1L)
    )),
    calibrated_outcome_nonconverged = sum(!vapply(
      diagnostics,
      function(diagnostic) diagnostic$calibrated_outcome_converged,
      logical(1L)
    )),
    initial_dr_cv_invalid_fold_fits = sum_site_diagnostic(
      "initial_density_ratio_cv_invalid_fold_fits"
    ),
    initial_dr_cv_invalid_lambdas = sum_site_diagnostic(
      "initial_density_ratio_cv_invalid_lambdas"
    ),
    initial_dr_cv_path_tail_skipped_fold_fits = sum_site_diagnostic(
      "initial_density_ratio_cv_path_tail_skipped_fold_fits"
    ),
    calibrated_dr_cv_invalid_fold_fits = sum_site_diagnostic(
      "calibrated_density_ratio_cv_invalid_fold_fits"
    ),
    calibrated_dr_cv_invalid_lambdas = sum_site_diagnostic(
      "calibrated_density_ratio_cv_invalid_lambdas"
    ),
    calibrated_dr_cv_path_tail_skipped_fold_fits = sum_site_diagnostic(
      "calibrated_density_ratio_cv_path_tail_skipped_fold_fits"
    ),
    calibrated_outcome_cv_invalid_fold_fits = sum_site_diagnostic(
      "calibrated_outcome_cv_invalid_fold_fits"
    ),
    calibrated_outcome_cv_invalid_lambdas = sum_site_diagnostic(
      "calibrated_outcome_cv_invalid_lambdas"
    ),
    calibrated_outcome_cv_path_tail_skipped_fold_fits = sum_site_diagnostic(
      "calibrated_outcome_cv_path_tail_skipped_fold_fits"
    ),
    initial_dr_line_search_failures = sum(vapply(
      diagnostics,
      function(diagnostic) {
        as.integer(sum(
          diagnostic$initial_density_ratio_line_search_failures,
          na.rm = TRUE
        ))
      },
      integer(1L)
    )),
    calibrated_dr_line_search_failures = sum(vapply(
      diagnostics,
      function(diagnostic) {
        as.integer(
          diagnostic$calibrated_density_ratio_line_search_failures %||% 0L
        )
      },
      integer(1L)
    )),
    calibrated_outcome_line_search_failures = sum(vapply(
      diagnostics,
      function(diagnostic) {
        as.integer(
          diagnostic$calibrated_outcome_line_search_failures %||% 0L
        )
      },
      integer(1L)
    )),
    initial_dr_support_floor_applied = sum(vapply(
      diagnostics,
      function(diagnostic) {
        sum(diagnostic$initial_density_ratio_support_floor_applied)
      },
      integer(1L)
    )),
    calibrated_dr_support_floor_applied = sum(vapply(
      diagnostics,
      function(diagnostic) {
        isTRUE(diagnostic$calibrated_density_ratio_support_floor_applied)
      },
      logical(1L)
    )),
    max_initial_dr_iterations = .max_iterations_or_na(initial_iterations),
    max_calibrated_dr_iterations =
      .max_iterations_or_na(calibrated_dr_iterations),
    max_calibrated_outcome_iterations =
      .max_iterations_or_na(calibrated_outcome_iterations),
    max_initial_dr_update_ratio =
      .max_numeric_or_na(initial_update_ratios),
    max_initial_dr_abs_coefficient =
      .max_numeric_or_na(initial_abs_coefficients),
    max_calibrated_dr_update_ratio =
      .max_numeric_or_na(calibrated_dr_update_ratios),
    max_calibrated_dr_abs_coefficient =
      .max_numeric_or_na(calibrated_dr_abs_coefficients),
    max_calibrated_outcome_update_ratio =
      .max_numeric_or_na(calibrated_outcome_update_ratios),
    max_calibrated_outcome_abs_coefficient =
      .max_numeric_or_na(calibrated_outcome_abs_coefficients),
    max_initial_dr_support_floor =
      .max_numeric_or_na(initial_support_floors),
    max_calibrated_dr_support_floor =
      .max_numeric_or_na(calibrated_support_floors)
  )
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
#' @param M_tau_inference Truncation radius used in inference-stage source
#'   corrections and influence-function moments.
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
#' @param nuisance_nlambda Number of candidates in each source nuisance-model
#'   CV lambda path.
#' @param nuisance_max_iter Maximum coordinate-descent iterations for each
#'   nuisance fit. Cross-validation kernels retain their separate internal cap;
#'   this argument primarily controls the final refit at the selected penalty.
#' @return List with source site results (mu_ts, gamma_s, alpha_ts, delta_ts, mu_pred_ts, n_s)
process_source_site <- function(s, source_folds, target_folds, k1, n_folds,
                                A_val, M_tau, M_tau_inference = M_TAU_INFERENCE_DEFAULT,
                                data_split,
                                combine_cache = NULL, get_fold_inputs,
                                family_int = 1L, link_int = 1L,
                                use_lambda_cache = TRUE,
                                nuisance_nlambda = LAMBDA_GRID_SIZE_STANDARD,
                                nuisance_max_iter = MAX_ITER_DEFAULT,
                                nuisance_lambda_rule = c("min", "1se")) {
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(nuisance_lambda_rule, "process_source_site")
  nuisance_max_iter <- suppressWarnings(as.integer(nuisance_max_iter))
  if (length(nuisance_max_iter) != 1L || is.na(nuisance_max_iter) ||
      nuisance_max_iter < 1L) {
    stop("process_source_site: nuisance_max_iter must be one positive integer.",
         call. = FALSE)
  }
  nuisance_nlambda <- suppressWarnings(as.integer(nuisance_nlambda))
  if (length(nuisance_nlambda) != 1L || is.na(nuisance_nlambda) ||
      nuisance_nlambda < 2L) {
    stop("process_source_site: nuisance_nlambda must be one integer >= 2.",
         call. = FALSE)
  }
  source_started_at <- proc.time()[["elapsed"]]
  timing <- c(
    initial_density_ratio = 0,
    initial_density_ratio_cv = 0,
    initial_density_ratio_final_fit = 0,
    calibrated_density_ratio = 0,
    calibrated_density_ratio_cv = 0,
    calibrated_density_ratio_final_fit = 0,
    calibrated_outcome = 0,
    calibrated_outcome_cv = 0,
    calibrated_outcome_final_fit = 0,
    correction_and_prediction = 0
  )

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
    fit_started_at <- proc.time()[["elapsed"]]
    gamma_init_k1_k2 <- fit_initial_density_ratio(
      source_train$Z_site, source_train$A, mean_phi_k2,
      lambda = if (isTRUE(use_lambda_cache)) cached_lambda_init_dr else NULL,
      A_val = A_val, warm_start = prev_gamma_init,
      nlambda = nuisance_nlambda,
      max_iter = nuisance_max_iter,
      lambda_rule = nuisance_lambda_rule,
      cv_group_id = source_train$cv_group_id
    )
    timing[["initial_density_ratio"]] <-
      timing[["initial_density_ratio"]] +
      .elapsed_process_seconds(fit_started_at)
    timing[["initial_density_ratio_cv"]] <-
      timing[["initial_density_ratio_cv"]] +
      as.numeric(attr(gamma_init_k1_k2, "cv_seconds"))
    timing[["initial_density_ratio_final_fit"]] <-
      timing[["initial_density_ratio_final_fit"]] +
      as.numeric(attr(gamma_init_k1_k2, "final_fit_seconds"))
    if (isTRUE(use_lambda_cache) && is.null(cached_lambda_init_dr)) {
      cached_lambda_init_dr <- attr(gamma_init_k1_k2, "lambda_used")
    }
    # Never carry a failed/boundary density-ratio solution into the next
    # secondary fold.  A warm start is only a computational device; resetting
    # it preserves the fitted objective while preventing a sparse support cell
    # in one fold from contaminating the next fold's optimization path.
    prev_gamma_init <- .safe_density_ratio_warm_start(gamma_init_k1_k2)
    
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
    cv_group_cal_stack <- if (is.null(source_calib_list[[1L]]$cv_group_id)) NULL else
      .stack_fold_field(source_calib_list, "cv_group_id", "process_source_site")

    mean_grad_psi_avg <- .average_numeric_list(mean_grad_psi_list, "process_source_site")

    W_plugin_block <- .make_plugin_block_design(
      lapply(source_calib_list, `[[`, "W_outcome"),
      caller = "process_source_site(gamma fold-summed calibration)"
    )
    alpha_plugin_block <- c(0, unlist(alpha_init_list, use.names = FALSE))
    gamma_warm_start <- .average_numeric_list(
      gamma_init_list,
      "process_source_site(gamma final warm start)"
    )

    fit_started_at <- proc.time()[["elapsed"]]
    gamma_final_k1 <- fit_unified_density_ratio(
      Z_site = Z_cal_stack, A = A_cal_stack,
      mean_grad_psi = mean_grad_psi_avg, alpha_init = alpha_plugin_block,
      lambda = NULL,
      calibrated = TRUE, M_tau = M_tau,
      W_outcome = W_plugin_block, A_val = A_val,
      family_int = family_int, link_int = link_int,
      warm_start = gamma_warm_start,
      nlambda = nuisance_nlambda,
      max_iter = nuisance_max_iter,
      lambda_rule = nuisance_lambda_rule,
      cv_group_id = cv_group_cal_stack
    )
    timing[["calibrated_density_ratio"]] <-
      .elapsed_process_seconds(fit_started_at)
    timing[["calibrated_density_ratio_cv"]] <-
      as.numeric(attr(gamma_final_k1, "cv_seconds"))
    timing[["calibrated_density_ratio_final_fit"]] <-
      as.numeric(attr(gamma_final_k1, "final_fit_seconds"))

    Z_plugin_block <- .make_plugin_block_design(
      lapply(source_calib_list, `[[`, "Z_site"),
      caller = "process_source_site(alpha fold-summed calibration)"
    )
    gamma_plugin_block <- c(0, unlist(gamma_init_list, use.names = FALSE))
    alpha_warm_start <- .average_numeric_list(
      alpha_init_list,
      "process_source_site(alpha final warm start)"
    )

    fit_started_at <- proc.time()[["elapsed"]]
    alpha_final_k1 <- fit_unified_outcome(
      W_outcome = W_cal_stack, Y = Y_cal_stack, A = A_cal_stack,
      A_val = A_val, gamma_s = gamma_plugin_block,
      lambda = NULL,
      family_int = family_int, link_int = link_int,
      calibrated = TRUE, M_tau = M_tau, Z_site = Z_plugin_block,
      warm_start = alpha_warm_start,
      nlambda = nuisance_nlambda,
      max_iter = nuisance_max_iter,
      lambda_rule = nuisance_lambda_rule,
      cv_group_id = cv_group_cal_stack
    )
    timing[["calibrated_outcome"]] <-
      .elapsed_process_seconds(fit_started_at)
    timing[["calibrated_outcome_cv"]] <-
      as.numeric(attr(alpha_final_k1, "cv_seconds"))
    timing[["calibrated_outcome_final_fit"]] <-
      as.numeric(attr(alpha_final_k1, "final_fit_seconds"))
  }
  
  # Compute correction term on main fold  — eq:site_membership in main.tex
  source_main_fold <- materialize_fold(source_folds[[s]], k1)
  
  if (source_main_fold$n == 0) {
    stop(sprintf("%s: source site '%s' main fold k1=%d has 0 observations; cannot compute source correction term.",
                 label, s, k1))
  }
  
  fit_started_at <- proc.time()[["elapsed"]]
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
  timing[["correction_and_prediction"]] <-
    .elapsed_process_seconds(fit_started_at)
  
  # Final source estimate: μ̂_{t,s_j} = mu_pred_ts + δ_{t,s_j}
  mu_ts_k1 <- mu_pred_ts_k1 + delta_ts_k1
  
  return(list(
    mu_ts = mu_ts_k1,
    gamma_s = gamma_final_k1,
    alpha_ts = alpha_final_k1,
    delta_ts = delta_ts_k1,
    mu_pred_ts = mu_pred_ts_k1,
    n_s = source_main_fold$n,
    nuisance_nlambda = nuisance_nlambda,
    correction_components = correction_result$correction_components,
    correction_clip_diagnostics = correction_result$clip_diagnostics,
    n_calibrated_folds = length(source_calib_list),
    calibrated_fold_keys = names(source_calib_list),
    # Per-k2 out-of-two-fold plug-ins for aggregation inner validation.
    # These are trained without folds k1 and k2, so validation summaries on k2
    # do not reuse the validation observations.
    per_k2_gamma = gamma_init_list,
    per_k2_alpha = alpha_init_list,
    nuisance_fit_diagnostics = list(
      initial_density_ratio_converged = vapply(
        gamma_init_list,
        function(fit) isTRUE(attr(fit, "converged")),
        logical(1L)
      ),
      initial_density_ratio_iterations = vapply(
        gamma_init_list,
        function(fit) as.integer(attr(fit, "iterations") %||% NA_integer_),
        integer(1L)
      ),
      initial_density_ratio_lambdas = vapply(
        gamma_init_list,
        function(fit) as.numeric(attr(fit, "lambda_used") %||% NA_real_),
        numeric(1L)
      ),
      initial_density_ratio_cv_invalid_fold_fits = vapply(
        gamma_init_list,
        function(fit) {
          as.integer(attr(fit, "cv_invalid_fold_fits") %||% 0L)
        },
        integer(1L)
      ),
      initial_density_ratio_cv_invalid_lambdas = vapply(
        gamma_init_list,
        function(fit) {
          as.integer(attr(fit, "cv_invalid_lambdas") %||% 0L)
        },
        integer(1L)
      ),
      initial_density_ratio_cv_path_tail_skipped_fold_fits = vapply(
        gamma_init_list,
        function(fit) {
          as.integer(attr(fit, "cv_path_tail_skipped_fold_fits") %||% 0L)
        },
        integer(1L)
      ),
      initial_density_ratio_update_ratios = vapply(
        gamma_init_list,
        function(fit) {
          as.numeric(attr(fit, "update_to_threshold_ratio") %||% NA_real_)
        },
        numeric(1L)
      ),
      initial_density_ratio_max_abs_coefficients = vapply(
        gamma_init_list,
        function(fit) {
          as.numeric(attr(fit, "max_abs_coefficient") %||% NA_real_)
        },
        numeric(1L)
      ),
      initial_density_ratio_line_search_failures = vapply(
        gamma_init_list,
        function(fit) {
          as.integer(attr(fit, "line_search_failures") %||% NA_integer_)
        },
        integer(1L)
      ),
      initial_density_ratio_support_floors = vapply(
        gamma_init_list,
        function(fit) {
          as.numeric(attr(fit, "support_penalty_floor") %||% NA_real_)
        },
        numeric(1L)
      ),
      initial_density_ratio_support_floor_applied = vapply(
        gamma_init_list,
        function(fit) isTRUE(attr(fit, "support_penalty_floor_applied")),
        logical(1L)
      ),
      calibrated_density_ratio_converged =
        isTRUE(attr(gamma_final_k1, "converged")),
      calibrated_density_ratio_iterations =
        as.integer(attr(gamma_final_k1, "iterations") %||% NA_integer_),
      calibrated_density_ratio_lambda =
        as.numeric(attr(gamma_final_k1, "lambda_used") %||% NA_real_),
      calibrated_density_ratio_cv_invalid_fold_fits = as.integer(
        attr(gamma_final_k1, "cv_invalid_fold_fits") %||% 0L
      ),
      calibrated_density_ratio_cv_invalid_lambdas = as.integer(
        attr(gamma_final_k1, "cv_invalid_lambdas") %||% 0L
      ),
      calibrated_density_ratio_cv_path_tail_skipped_fold_fits = as.integer(
        attr(gamma_final_k1, "cv_path_tail_skipped_fold_fits") %||% 0L
      ),
      calibrated_density_ratio_update_ratio = as.numeric(
        attr(gamma_final_k1, "update_to_threshold_ratio") %||% NA_real_
      ),
      calibrated_density_ratio_max_abs_coefficient = as.numeric(
        attr(gamma_final_k1, "max_abs_coefficient") %||% NA_real_
      ),
      calibrated_density_ratio_line_search_failures = as.integer(
        attr(gamma_final_k1, "line_search_failures") %||% NA_integer_
      ),
      calibrated_density_ratio_support_floor = as.numeric(
        attr(gamma_final_k1, "support_penalty_floor") %||% NA_real_
      ),
      calibrated_density_ratio_support_floor_applied =
        isTRUE(attr(gamma_final_k1, "support_penalty_floor_applied")),
      calibrated_outcome_converged =
        isTRUE(attr(alpha_final_k1, "converged")),
      calibrated_outcome_iterations =
        as.integer(attr(alpha_final_k1, "iterations") %||% NA_integer_),
      calibrated_outcome_lambda =
        as.numeric(attr(alpha_final_k1, "lambda_used") %||% NA_real_),
      calibrated_outcome_cv_invalid_fold_fits = as.integer(
        attr(alpha_final_k1, "cv_invalid_fold_fits") %||% 0L
      ),
      calibrated_outcome_cv_invalid_lambdas = as.integer(
        attr(alpha_final_k1, "cv_invalid_lambdas") %||% 0L
      ),
      calibrated_outcome_cv_path_tail_skipped_fold_fits = as.integer(
        attr(alpha_final_k1, "cv_path_tail_skipped_fold_fits") %||% 0L
      ),
      calibrated_outcome_update_ratio = as.numeric(
        attr(alpha_final_k1, "update_to_threshold_ratio") %||% NA_real_
      ),
      calibrated_outcome_max_abs_coefficient = as.numeric(
        attr(alpha_final_k1, "max_abs_coefficient") %||% NA_real_
      ),
      calibrated_outcome_line_search_failures = as.integer(
        attr(alpha_final_k1, "line_search_failures") %||% NA_integer_
      )
    ),
    timing = c(
      timing,
      total = .elapsed_process_seconds(source_started_at)
    )
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
  result <- list(
    W_outcome = data$W_outcome[idx, , drop = FALSE],
    Z_site = data$Z_site[idx, , drop = FALSE],
    A = data$A[idx],
    Y = data$Y[idx],
    n = fold$n,
    original_idx = idx
  )
  if (!is.null(data$cv_group_id)) {
    groups <- .validate_nuisance_cv_group_id(data$cv_group_id, data$n, "materialize_fold")
    result$cv_group_id <- groups[idx]
  }
  result
}

#' Extract a scalar source-result field across folds and sites
#'
#' @param fold_results List of fold results
#' @param source_sites Vector of source site names
#' @param n_folds Number of folds
#' @param field Name of the finite scalar field to extract.
#' @return Numeric matrix with folds in rows and source sites in columns.
extract_source_result_matrix <- function(
    fold_results, source_sites, n_folds, field) {
  K <- length(source_sites)
  values <- vapply(seq_len(n_folds), function(k) {
    vapply(source_sites, function(s) {
      value <- fold_results[[k]]$source_results[[s]][[field]]
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop(sprintf(
          "extract_source_result_matrix: field '%s' must be a finite scalar for fold %d, source '%s'.",
          field, k, s
        ))
      }
      value
    }, FUN.VALUE = numeric(1))
  }, FUN.VALUE = numeric(K))
  # vapply drops to a vector when K = 1. Reconstruct the documented
  # K-by-fold layout explicitly before transposing so every K has the same
  # fold-by-source return shape.
  t(matrix(values, nrow = K, ncol = n_folds))
}

extract_source_estimates_matrix <- function(fold_results, source_sites, n_folds) {
  extract_source_result_matrix(
    fold_results, source_sites, n_folds, field = "mu_ts"
  )
}

.aggregate_crossfit_from_fitted_folds <- function(
    data_split, target_folds, source_folds, fold_results, n_folds,
    M_tau, M_tau_inference, lambda_selection, lambda_rule,
    aggregation_lambda_grid,
    communication_mode, family_int, link_int, A_val, verbose = FALSE) {
  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  K <- length(source_sites)

  target_estimates <- vapply(
    fold_results, function(fit) fit$target_only$estimate, numeric(1L)
  )
  target_fold_sizes <- vapply(
    fold_results,
    function(fit) length(fit$target_only$varphi_ot),
    numeric(1L)
  )
  final_target_estimate <- mean(target_estimates)

  source_estimates_matrix <- extract_source_estimates_matrix(
    fold_results, source_sites, n_folds
  )
  source_mu_pred_matrix <- extract_source_result_matrix(
    fold_results, source_sites, n_folds, "mu_pred_ts"
  )
  source_delta_matrix <- extract_source_result_matrix(
    fold_results, source_sites, n_folds, "delta_ts"
  )
  source_fold_sizes <- extract_source_result_matrix(
    fold_results, source_sites, n_folds, "n_s"
  )
  source_estimates <- .pool_fold_pairwise_estimates(
    source_mu_pred_matrix, source_delta_matrix,
    target_fold_sizes, source_fold_sizes
  )
  names(source_estimates) <- source_sites

  calculate_crossfit_aggregation(
    data_split = data_split,
    target_data = target_data,
    source_sites = source_sites,
    K = K,
    target_folds = target_folds,
    source_folds = source_folds,
    fold_results = fold_results,
    n_folds = n_folds,
    M_tau = M_tau,
    M_tau_inference = M_tau_inference,
    lambda_selection = lambda_selection,
    verbose = verbose,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
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
}

.get_target_only_fold_fit <- function(
    target_folds, k1, n_folds, family, A_val, propensity_cache,
    fit_cache = NULL, k2 = NULL,
    nuisance_lambda_rule = c("min", "1se")) {
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, ".get_target_only_fold_fit",
    arg = "nuisance_lambda_rule"
  )
  if (!is.null(fit_cache) && !is.environment(fit_cache)) {
    stop(".get_target_only_fold_fit: fit_cache must be NULL or an environment.")
  }
  cache_key <- sprintf(
    "family=%s|A=%d|k1=%d|k2=%s|nuisance_rule=%s",
    family, as.integer(A_val), as.integer(k1),
    if (is.null(k2)) "outer" else as.character(as.integer(k2)),
    nuisance_lambda_rule
  )
  cv_groups <- attr(target_folds, ".data_ref")$cv_group_id
  if (!is.null(cv_groups)) {
    cv_groups <- .validate_nuisance_cv_group_id(
      cv_groups, attr(target_folds, ".data_ref")$n, ".get_target_only_fold_fit"
    )
    cache_key <- paste0(cache_key, "|cv_groups=", paste(cv_groups, collapse = ","))
  }
  if (!is.null(fit_cache) &&
      exists(cache_key, envir = fit_cache, inherits = FALSE)) {
    return(get(cache_key, envir = fit_cache, inherits = FALSE))
  }

  fit <- estimate_target_only_from_complement(
    target_folds = target_folds,
    k1 = k1,
    n_folds = n_folds,
    k2 = k2,
    family = family,
    A_val = A_val,
    propensity_cache = propensity_cache,
    nuisance_lambda_rule = nuisance_lambda_rule
  )
  if (!is.null(fit_cache)) {
    assign(cache_key, fit, envir = fit_cache)
  }
  fit
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

#' Deterministic random seed for a cross-fitting work unit
#'
#' Nuisance-model cross-validation uses random fold assignments.  Giving each
#' arm/site work unit its own isolated stream makes estimates invariant to
#' whether those units are evaluated sequentially or by fork workers.  The
#' simulation ID need not enter this seed: the fitted data already change by
#' replication, while a common CV partitioning rule improves comparability.
#'
#' @param communication_mode One- or two-round protocol.
#' @param A_val Treatment arm, 0 or 1.
#' @param outer_fold Outer cross-fitting fold; zero denotes an arm-level unit.
#' @param source_index Source-site position; zero denotes an arm-level unit.
#' @return A positive integer seed.
#' @keywords internal
.crossfit_work_seed <- function(
    communication_mode, A_val, outer_fold = 0L, source_index = 0L) {
  mode_index <- match(communication_mode, c("one_round", "two_round"))
  values <- c(A_val, outer_fold, source_index)
  if (is.na(mode_index) || any(!is.finite(values)) ||
      A_val < 0L || A_val > 1L || outer_fold < 0L || source_index < 0L) {
    stop(".crossfit_work_seed: invalid work-unit identifiers.", call. = FALSE)
  }
  as.integer(
    500003L + 100000L * mode_index + 10000L * as.integer(A_val) +
      100L * as.integer(outer_fold) + as.integer(source_index)
  )
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

#' Unified Cross-fitting Algorithm for RoCE
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
#' @param data_split Split data by site (from split_data_by_site). A site may
#'   optionally include positive integer cv_group_id values aligned with its
#'   rows. Repeated observation origins are kept together in nuisance CV.
#'   NULL or all-unique IDs preserve the ordinary nuisance CV path. Repeated
#'   origins require group-consistent precomputed outer folds; stale views or
#'   origins crossing outer folds are rejected. This is resampling metadata,
#'   not a change to the common TATE aggregation objective.
#' @param n_folds number of cross-fitting folds (default 10, K_f in main.tex).
#'        Minimum 3 required for proper two-level cross-fitting calibration.
#' @param communication_mode character, "two_round" or "one_round"
#' @param lambda_selection Aggregation Wald-penalty factor. The default
#'   \code{AGG_WALD_LAMBDA = 1} corresponds to the pilot-selected
#'   penalty-activation cutoff \eqn{c = 1/\lambda = 1} used in
#'   \code{main.tex}. The legacy value
#'   \code{"cv"} remains available for a user-supplied sensitivity grid.
#' @param lambda_rule rule for the aggregation penalty when
#'        \code{lambda_selection = "cv"}: \code{"min"} (default) selects the
#'        bare validation minimizer. \code{"1se"} selects the largest lambda
#'        whose average inner-validation criterion is within one standard error
#'        of the minimum.
#' @param aggregation_lambda_grid Optional positive numeric vector of
#'   truncated-Wald multipliers used when \code{lambda_selection = "cv"}.
#'   Selection is performed using inner folds contained entirely in each
#'   outer-training sample. It must be \code{NULL} for a fixed numeric
#'   \code{lambda_selection}.
#' @param verbose print progress messages
#' @param M_tau truncation radius for the calibrated losses during training
#'        (default \code{M_TAU_DEFAULT} = 5; a single SMMAL-style radius shared
#'        with \code{M_tau_inference}).
#' @param M_tau_inference truncation radius for the inference step (correction term
#'        and source-variance calculation). Default \code{M_TAU_INFERENCE_DEFAULT}
#'        (= \code{M_TAU_DEFAULT} = 5): a single SMMAL-style radius governs both the
#'        calibrated losses and the final DR estimator (the truncated weight
#'        w_T = exp(-T(phi^T gamma)) feeds the estimate and its variance), chosen
#'        large enough to be asymptotically inactive under bounded tilt yet finite
#'        enough to cap rare finite-sample weight excursions.
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
#' @param precomputed_folds Optional list with \code{target_folds} and
#'        \code{source_folds} created in advance for this \code{data_split}.
#' @param target_only_ps_cache Optional environment used to cache target-only
#'        complement-fold propensity predictions.
#' @param target_only_fit_cache Optional environment used to reuse complete
#'        target-only outer/inner fold fits across communication modes. Use the
#'        same cache only with the same data and fold partition.
#' @return list with estimates and variance components
#' @export
run_crossfit <- function(data_split, n_folds = N_FOLDS_DEFAULT,
                         communication_mode = c("two_round", "one_round"),
                         lambda_selection = AGG_WALD_LAMBDA,
                         lambda_rule = c("min", "1se"),
                         verbose = TRUE, M_tau = M_TAU_DEFAULT,
                         M_tau_inference = M_TAU_INFERENCE_DEFAULT,
                         n_cores = NULL,
                         nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
                         family = "binomial", A_val = 1L,
                         use_lambda_cache = TRUE,
                         precomputed_folds = NULL,
                         target_only_ps_cache = NULL,
                         target_only_fit_cache = NULL,
                         nuisance_lambda_rule = c("min", "1se"),
                         aggregation_lambda_grid = NULL) {

  communication_mode <- match.arg(communication_mode)
  lambda_rule <- match.arg(lambda_rule)
  aggregation_lambda_grid <- .validate_aggregation_lambda_grid(
    lambda_selection, aggregation_lambda_grid, "run_crossfit"
  )
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(nuisance_lambda_rule, "run_crossfit")
  algorithm_started_at <- proc.time()[["elapsed"]]

  # ---- SHARED SETUP ----
  glm_spec <- resolve_glm_family(family)
  link <- glm_spec$link
  family_int <- glm_spec$family_int
  link_int <- glm_spec$link_int

  validate_algorithm_inputs(data_split, lambda_selection, family = family, A_val = A_val)
  for (site in names(data_split)) {
    .validate_nuisance_cv_group_id(data_split[[site]]$cv_group_id,
                                  data_split[[site]]$n, "run_crossfit")
  }
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
  partition_started_at <- proc.time()[["elapsed"]]
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
  .validate_nuisance_cv_fold_views(target_data, target_folds, "run_crossfit target")
  for (site in source_sites) {
    .validate_nuisance_cv_fold_views(data_split[[site]], source_folds[[site]],
                                    paste("run_crossfit source", site))
  }
  partition_seconds <- .elapsed_process_seconds(partition_started_at)

  fold_results <- list()
  fold_timing <- vector("list", n_folds)
  fold_fit_diagnostics <- vector("list", n_folds)

  # ---- MAIN LOOP ----
  for (k1 in 1:n_folds) {
    fold_started_at <- proc.time()[["elapsed"]]
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
    initial_outcome_degenerate <- 0L

    # ================================================================
    # COMMUNICATION-MODE SPECIFIC: prepare initial models & summaries
    # ================================================================
    initial_nuisance_started_at <- proc.time()[["elapsed"]]
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
              "Two-round: source site '%s' training folds {%s} have 0 observations for k2=%d; cannot fit the initial outcome nuisance required by main.tex eq:nuisance_initial_losses.",
              s, paste(training_folds, collapse = ","), k2
            ))
          }

          alpha_init_k1_k2 <- fit_initial_outcome(
            source_train$W_outcome, source_train$Y, source_train$A, A_val,
            lambda = if (isTRUE(use_lambda_cache)) cached_lambda_site else NULL,
            nlambda = nlambda_init, family = family,
            lambda_rule = nuisance_lambda_rule,
            cv_group_id = source_train$cv_group_id
          )
          initial_outcome_degenerate <- initial_outcome_degenerate +
            as.integer(identical(
              attr(alpha_init_k1_k2, "lambda_rule"),
              "degenerate_constant"
            ))
          if (isTRUE(use_lambda_cache) && is.null(cached_lambda_site)) {
            lu <- attr(alpha_init_k1_k2, "lambda_used")
            # Skip degenerate folds (non-finite lambda) so they do not poison the
            # cache; a later non-degenerate fold supplies the cached value.
            if (!is.null(lu) && is.finite(lu)) {
              cached_lambda_site <- lu
              if (verbose) {
                cat(sprintf("      Cached lambda (site %s) = %.6f\n", s, cached_lambda_site))
              }
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
          lambda_rule = nuisance_lambda_rule,
          cv_group_id = target_train$cv_group_id
        )
        initial_outcome_degenerate <- initial_outcome_degenerate +
          as.integer(identical(
            attr(alpha_init_k1_k2, "lambda_rule"),
            "degenerate_constant"
          ))
        if (isTRUE(use_lambda_cache) && is.null(cached_lambda)) {
          lu <- attr(alpha_init_k1_k2, "lambda_used")
          # Skip degenerate folds (non-finite lambda) so they do not poison the
          # cache; a later non-degenerate fold supplies the cached value.
          if (!is.null(lu) && is.finite(lu)) {
            cached_lambda <- lu
            if (verbose) {
              cat(sprintf("      Cached lambda = %.6f (reused for remaining calls)\n", cached_lambda))
            }
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
    initial_nuisance_seconds <-
      .elapsed_process_seconds(initial_nuisance_started_at)

    # ================================================================
    # SHARED: Process source sites via unified helper
    # ================================================================
    source_processing_started_at <- proc.time()[["elapsed"]]
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
      seq_along(source_sites),
      function(source_index) {
        s <- source_sites[[source_index]]
        with_seed(
          .crossfit_work_seed(
            communication_mode, A_val, k1, source_index
          ),
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
            nuisance_nlambda = nlambda_init,
            nuisance_lambda_rule = nuisance_lambda_rule
          )
        )
      },
      n_cores = actual_cores
    )

    source_results_k1 <- setNames(source_results_list, source_sites)
    observed_nlambda <- vapply(
      source_results_k1,
      function(result) as.integer(result$nuisance_nlambda),
      integer(1L)
    )
    if (anyNA(observed_nlambda) || any(observed_nlambda != nlambda_init)) {
      stop(
        sprintf(
          paste0(
            "run_crossfit: source nuisance grid mismatch in outer fold %d; ",
            "requested %d but observed {%s}."
          ),
          k1, nlambda_init, paste(unique(observed_nlambda), collapse = ",")
        ),
        call. = FALSE
      )
    }
    source_processing_wall_seconds <-
      .elapsed_process_seconds(source_processing_started_at)

    if (verbose) {
      cat(sprintf("      [OK] Processed %d source sites\n", length(source_sites)))
    }

    # Target-only estimate (complement-fold training, avoids nested cross-fitting)
    target_only_started_at <- proc.time()[["elapsed"]]
    target_only_k1 <- .get_target_only_fold_fit(
      target_folds = target_folds,
      k1 = k1,
      n_folds = n_folds,
      family = family,
      A_val = A_val,
      propensity_cache = target_only_ps_cache,
      fit_cache = target_only_fit_cache,
      nuisance_lambda_rule = nuisance_lambda_rule
    )

    # Inner-fold target-only estimates (Version A Step 2)
    # For each k2 in secondary_folds, train target-only model on
    # {folds}\{k1,k2} and evaluate on fold k2. These are used to compute
    # varphi_ot on inner folds for variance component estimation.
    log_info(verbose, "   --> Inner-fold target-only models for Version A Step 2")
    target_only_inner <- list()
    for (k2 in secondary_folds) {
      target_only_inner[[paste0("k2_", k2)]] <- .get_target_only_fold_fit(
        target_folds = target_folds,
        k1 = k1,
        n_folds = n_folds,
        family = family,
        A_val = A_val,
        propensity_cache = target_only_ps_cache,
        fit_cache = target_only_fit_cache,
        k2 = k2,
        nuisance_lambda_rule = nuisance_lambda_rule
      )
    }
    target_only_outcome_degenerate <-
      as.integer(target_only_k1$outcome_degenerate %||% 0L) +
      sum(vapply(
        target_only_inner,
        function(fit) as.integer(fit$outcome_degenerate %||% 0L),
        integer(1L)
      ))
    target_only_seconds <- .elapsed_process_seconds(target_only_started_at)

    source_timing_matrix <- do.call(
      rbind,
      lapply(source_results_k1, function(site_result) site_result$timing)
    )
    source_timing_sums <- colSums(source_timing_matrix)
    fold_timing[[k1]] <- c(
      fold = k1,
      initial_nuisance_seconds = initial_nuisance_seconds,
      source_processing_wall_seconds = source_processing_wall_seconds,
      target_only_seconds = target_only_seconds,
      source_initial_dr_sum_seconds =
        source_timing_sums[["initial_density_ratio"]],
      source_initial_dr_cv_sum_seconds =
        source_timing_sums[["initial_density_ratio_cv"]],
      source_initial_dr_final_fit_sum_seconds =
        source_timing_sums[["initial_density_ratio_final_fit"]],
      source_calibrated_dr_sum_seconds =
        source_timing_sums[["calibrated_density_ratio"]],
      source_calibrated_dr_cv_sum_seconds =
        source_timing_sums[["calibrated_density_ratio_cv"]],
      source_calibrated_dr_final_fit_sum_seconds =
        source_timing_sums[["calibrated_density_ratio_final_fit"]],
      source_calibrated_outcome_sum_seconds =
        source_timing_sums[["calibrated_outcome"]],
      source_calibrated_outcome_cv_sum_seconds =
        source_timing_sums[["calibrated_outcome_cv"]],
      source_calibrated_outcome_final_fit_sum_seconds =
        source_timing_sums[["calibrated_outcome_final_fit"]],
      source_correction_sum_seconds =
        source_timing_sums[["correction_and_prediction"]],
      fold_total_seconds = .elapsed_process_seconds(fold_started_at)
    )
    fold_fit_diagnostics[[k1]] <- c(
      fold = k1,
      initial_outcome_degenerate = initial_outcome_degenerate,
      target_only_outcome_degenerate = target_only_outcome_degenerate,
      .summarize_source_nuisance_fit_diagnostics(source_results_k1)
    )

    fold_results[[k1]] <- list(
      target_only = target_only_k1,
      source_results = source_results_k1,
      target_only_inner = target_only_inner
    )
  }

  # ---- FINAL AGGREGATION ----
  aggregation_started_at <- proc.time()[["elapsed"]]
  if (verbose) {
    cat(sprintf("\n--> FINAL AGGREGATION (%s)\n", mode_label))
  }

  result <- .aggregate_crossfit_from_fitted_folds(
    data_split = data_split,
    target_folds = target_folds,
    source_folds = source_folds,
    fold_results = fold_results,
    n_folds = n_folds,
    M_tau = M_tau,
    M_tau_inference = M_tau_inference,
    lambda_selection = lambda_selection,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
    communication_mode = communication_mode,
    family_int = family_int,
    link_int = link_int,
    A_val = A_val,
    verbose = verbose
  )
  result$nuisance_lambda_rule <- nuisance_lambda_rule
  result$communication_mode <- communication_mode
  result$family <- family
  result$A_val <- A_val
  result$M_tau <- M_tau
  result$M_tau_inference <- M_tau_inference
  result$aggregation_lambda_selection <- lambda_selection
  result$aggregation_lambda_grid <- aggregation_lambda_grid
  result$timing <- list(
    total_seconds = .elapsed_process_seconds(algorithm_started_at),
    partition_seconds = partition_seconds,
    aggregation_seconds = .elapsed_process_seconds(aggregation_started_at),
    folds = as.data.frame(do.call(rbind, fold_timing), row.names = NULL)
  )
  result$nuisance_fit_diagnostics <- as.data.frame(
    do.call(rbind, fold_fit_diagnostics),
    row.names = NULL
  )
  return(result)
}

#' Direct cross-fitted RoCE estimator of the TATE
#'
#' Fits arm-specific high-dimensional nuisances using a common fold partition,
#' then aggregates the target-only and source-assisted TATE estimators with one
#' shared source-weight vector. The aggregation objective and final variance are
#' computed from treated-minus-control influence-function contrasts.
#'
#' @inheritParams run_crossfit
#' @return A RoCE TATE result. Arm-specific fits are combined before
#'   weight estimation; the returned object contains the TATE-level estimates,
#'   common weights, fold diagnostics, and target-only reference.
#' @param parallel_arms Logical. If \code{TRUE}, fit the treated and control
#'   nuisance pipelines concurrently. Callers must allocate approximately
#'   \code{2 * n_cores} CPUs because each arm can still parallelize over source
#'   sites. Defaults to \code{FALSE} for backward-compatible resource use.
#' @param screening_rule Source-screening rule passed to
#'   \code{\link{calculate_tate_crossfit_aggregation}}. The default retains the
#'   manuscript soft penalty; hard thresholding is a diagnostic and
#'   \code{"quadratic_bias"} the pre-specified smooth sensitivity rule.
#' @export
run_tate_crossfit <- function(
    data_split, n_folds = N_FOLDS_DEFAULT,
    communication_mode = c("two_round", "one_round"),
    lambda_selection = AGG_WALD_LAMBDA,
    lambda_rule = c("min", "1se"),
    verbose = TRUE,
    M_tau = M_TAU_DEFAULT,
    M_tau_inference = M_TAU_INFERENCE_DEFAULT,
    n_cores = NULL,
    nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
    family = "binomial",
    use_lambda_cache = TRUE,
    precomputed_folds = NULL,
    target_only_ps_cache = NULL,
    target_only_fit_cache = NULL,
    nuisance_lambda_rule = c("min", "1se"),
    parallel_arms = FALSE,
    aggregation_lambda_grid = NULL,
    screening_rule = c("soft_penalty", "hard_threshold", "quadratic_bias")) {
  tate_started_at <- proc.time()[["elapsed"]]
  communication_mode <- match.arg(communication_mode)
  lambda_rule <- match.arg(lambda_rule)
  aggregation_lambda_grid <- .validate_aggregation_lambda_grid(
    lambda_selection, aggregation_lambda_grid, "run_tate_crossfit"
  )
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, "run_tate_crossfit"
  )
  screening_rule <- match.arg(screening_rule)
  if (length(parallel_arms) != 1L || is.na(parallel_arms) ||
      !is.logical(parallel_arms)) {
    stop("run_tate_crossfit: parallel_arms must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (isTRUE(parallel_arms) && !is.null(n_cores) && n_cores == -1) {
    stop(
      paste0(
        "run_tate_crossfit: parallel_arms=TRUE requires an explicit positive ",
        "n_cores per arm; n_cores=-1 would oversubscribe the node."
      ),
      call. = FALSE
    )
  }

  if (is.null(precomputed_folds)) {
    precomputed_folds <- build_crossfit_folds(data_split, n_folds)
  }
  if (is.null(target_only_ps_cache)) {
    target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
  }
  if (is.null(target_only_fit_cache)) {
    target_only_fit_cache <- new.env(hash = TRUE, parent = emptyenv())
  }

  common_args <- list(
    data_split = data_split,
    n_folds = n_folds,
    communication_mode = communication_mode,
    lambda_selection = lambda_selection,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
    verbose = verbose,
    M_tau = M_tau,
    M_tau_inference = M_tau_inference,
    n_cores = n_cores,
    nlambda_init = nlambda_init,
    family = family,
    use_lambda_cache = use_lambda_cache,
    precomputed_folds = precomputed_folds,
    target_only_ps_cache = target_only_ps_cache,
    target_only_fit_cache = target_only_fit_cache,
    nuisance_lambda_rule = nuisance_lambda_rule
  )
  arm_values <- c(mu1 = 1L, mu0 = 0L)
  arm_results <- parallel_lapply(
    arm_values,
    function(A_val) {
      with_seed(
        .crossfit_work_seed(communication_mode, A_val),
        do.call(run_crossfit, c(common_args, list(A_val = A_val)))
      )
    },
    n_cores = if (isTRUE(parallel_arms)) 2L else 1L
  )
  mu1_result <- arm_results$mu1
  mu0_result <- arm_results$mu0

  result <- calculate_tate_crossfit_aggregation(
    data_split = data_split,
    mu1_result = mu1_result,
    mu0_result = mu0_result,
    lambda_selection = lambda_selection,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
    screening_rule = screening_rule,
    verbose = verbose
  )
  result$communication_mode <- communication_mode
  result$parallel_arms <- parallel_arms
  result$nuisance_lambda_rule <- nuisance_lambda_rule
  result$family <- family
  result$M_tau <- M_tau
  result$M_tau_inference <- M_tau_inference
  result$aggregation_lambda_selection <- lambda_selection
  result$aggregation_lambda_grid <- aggregation_lambda_grid
  result$aggregation_screening_rule <- screening_rule
  result$method <- paste0(communication_mode, "_direct_tate_crossfit")
  result$timing <- list(
    total_seconds = .elapsed_process_seconds(tate_started_at),
    mu1 = mu1_result$timing,
    mu0 = mu0_result$timing
  )
  result$nuisance_fit_diagnostics <- rbind(
    transform(mu1_result$nuisance_fit_diagnostics, A_val = 1L),
    transform(mu0_result$nuisance_fit_diagnostics, A_val = 0L)
  )
  rownames(result$nuisance_fit_diagnostics) <- NULL
  result
}

.validate_one_round_rho_reuse_data <- function(
    reference_data_split, data_split, changed_sources, A_val,
    caller = ".refit_one_round_crossfit_sources") {
  if (!is.list(reference_data_split) || !is.list(data_split) ||
      !identical(names(reference_data_split), names(data_split))) {
    stop(caller, ": reference and current site lists must have identical names.",
         call. = FALSE)
  }
  source_sites <- setdiff(names(data_split), "t")
  changed_sources <- unique(as.character(changed_sources))
  if (length(changed_sources) == 0L ||
      any(!changed_sources %in% source_sites)) {
    stop(caller, ": changed_sources must identify at least one source site.",
         call. = FALSE)
  }

  stable_fields <- c(
    "n", "X", "X_dagger", "A", "Z_site_true", "W_outcome_true",
    "Z_site", "W_outcome"
  )
  for (site in names(data_split)) {
    reference_site <- reference_data_split[[site]]
    current_site <- data_split[[site]]
    if (!identical(reference_site$cv_group_id, current_site$cv_group_id)) {
      stop(caller, ": nuisance CV groups changed at site '", site, "'.", call. = FALSE)
    }
    missing_fields <- stable_fields[
      !stable_fields %in% names(reference_site) |
        !stable_fields %in% names(current_site)
    ]
    if (length(missing_fields) > 0L) {
      stop(
        caller, ": site '", site, "' lacks invariant field(s): ",
        paste(missing_fields, collapse = ", "), ".", call. = FALSE
      )
    }
    invariant <- vapply(
      stable_fields,
      function(field) identical(reference_site[[field]], current_site[[field]]),
      logical(1L)
    )
    if (!all(invariant)) {
      stop(
        caller, ": site '", site,
        "' differs in field(s) that must be rho-invariant: ",
        paste(names(invariant)[!invariant], collapse = ", "), ".",
        call. = FALSE
      )
    }

    if (identical(site, "t") || !site %in% changed_sources) {
      if (!identical(reference_site$Y, current_site$Y)) {
        stop(caller, ": outcome changed at non-refitted site '", site, "'.",
             call. = FALSE)
      }
    } else {
      unchanged_arm <- current_site$A != as.integer(A_val)
      if (!identical(
        reference_site$Y[unchanged_arm], current_site$Y[unchanged_arm]
      )) {
        stop(
          caller, ": source '", site,
          "' changed outcomes outside the refitted treatment arm.",
          call. = FALSE
        )
      }
    }
  }
  changed_sources
}

.refit_one_round_crossfit_sources <- function(
    reference_data_split, data_split, fitted_arm, changed_sources,
    lambda_selection = fitted_arm$aggregation_lambda_selection,
    lambda_rule = fitted_arm$aggregation_lambda_rule %||% "min",
    aggregation_lambda_grid = fitted_arm$aggregation_lambda_grid %||% NULL,
    verbose = FALSE, n_cores = 1L) {
  started_at <- proc.time()[["elapsed"]]
  caller <- ".refit_one_round_crossfit_sources"
  required <- c(
    "fold_results", "n_folds", "communication_mode", "family", "A_val",
    "M_tau", "M_tau_inference", "nuisance_lambda_rule",
    "nuisance_fit_diagnostics"
  )
  missing <- required[!vapply(required, function(field) {
    !is.null(fitted_arm[[field]])
  }, logical(1L))]
  if (length(missing) > 0L) {
    stop(caller, ": fitted arm lacks field(s): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  if (!identical(fitted_arm$communication_mode, "one_round")) {
    stop(caller, ": only one-round fits can be reused across FACE rho values.",
         call. = FALSE)
  }
  A_val <- as.integer(fitted_arm$A_val)
  changed_sources <- .validate_one_round_rho_reuse_data(
    reference_data_split, data_split, changed_sources, A_val, caller
  )
  n_folds <- as.integer(fitted_arm$n_folds)
  current_folds <- build_crossfit_folds(data_split, n_folds)
  reference_folds <- build_crossfit_folds(reference_data_split, n_folds)
  source_sites <- setdiff(names(data_split), "t")
  .validate_reaggregation_folds(
    fitted_arm, reference_folds$target_folds, reference_folds$source_folds,
    source_sites, caller
  )
  .validate_reaggregation_folds(
    fitted_arm, current_folds$target_folds, current_folds$source_folds,
    source_sites, caller
  )

  first_source_result <- fitted_arm$fold_results[[1L]]$source_results[[1L]]
  nuisance_nlambda <- as.integer(first_source_result$nuisance_nlambda)
  if (length(nuisance_nlambda) != 1L || is.na(nuisance_nlambda) ||
      nuisance_nlambda < 2L) {
    stop(caller, ": fitted arm has an invalid nuisance grid size.",
         call. = FALSE)
  }
  glm_spec <- resolve_glm_family(fitted_arm$family)
  fold_results <- fitted_arm$fold_results
  refit_wall_seconds <- numeric(n_folds)

  for (k1 in seq_len(n_folds)) {
    fold_started_at <- proc.time()[["elapsed"]]
    secondary_folds <- setdiff(seq_len(n_folds), k1)
    combine_cache <- new.env(hash = TRUE, parent = emptyenv())
    refitted <- parallel_lapply(
      changed_sources,
      function(site) {
        source_index <- match(site, source_sites)
        baseline_source <- fitted_arm$fold_results[[k1]]$source_results[[site]]
        get_fold_inputs <- function(.site, k2) {
          k2_key <- paste0("k2_", k2)
          alpha_init <- baseline_source$per_k2_alpha[[k2_key]]
          if (is.null(alpha_init)) {
            stop(
              caller, ": baseline source '", site, "' fold ", k1,
              " lacks ", k2_key, " initial outcome coefficients.",
              call. = FALSE
            )
          }
          training_folds <- setdiff(seq_len(n_folds), c(k1, k2))
          target_train <- combine_folds(
            current_folds$target_folds, training_folds
          )
          target_calibration <- materialize_fold(
            current_folds$target_folds, k2
          )
          list(
            mean_phi = c(1, colMeans(target_train$Z_site)),
            mean_grad_psi_init = .mean_glm_gradient_site_basis(
              target_calibration$W_outcome,
              target_calibration$Z_site,
              alpha_init,
              glm_spec$family_int,
              glm_spec$link_int
            ),
            alpha_init = alpha_init
          )
        }
        with_seed(
          .crossfit_work_seed("one_round", A_val, k1, source_index),
          process_source_site(
            s = site,
            source_folds = current_folds$source_folds,
            target_folds = current_folds$target_folds,
            k1 = k1,
            n_folds = n_folds,
            A_val = A_val,
            M_tau = fitted_arm$M_tau,
            M_tau_inference = fitted_arm$M_tau_inference,
            data_split = data_split,
            combine_cache = combine_cache,
            get_fold_inputs = get_fold_inputs,
            family_int = glm_spec$family_int,
            link_int = glm_spec$link_int,
            use_lambda_cache = TRUE,
            nuisance_nlambda = nuisance_nlambda,
            nuisance_lambda_rule = fitted_arm$nuisance_lambda_rule
          )
        )
      },
      n_cores = min(setup_parallel(n_cores), length(changed_sources))
    )
    names(refitted) <- changed_sources
    for (site in changed_sources) {
      fold_results[[k1]]$source_results[[site]] <- refitted[[site]]
    }
    refit_wall_seconds[[k1]] <- .elapsed_process_seconds(fold_started_at)
  }

  result <- .aggregate_crossfit_from_fitted_folds(
    data_split = data_split,
    target_folds = current_folds$target_folds,
    source_folds = current_folds$source_folds,
    fold_results = fold_results,
    n_folds = n_folds,
    M_tau = fitted_arm$M_tau,
    M_tau_inference = fitted_arm$M_tau_inference,
    lambda_selection = lambda_selection,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
    communication_mode = "one_round",
    family_int = glm_spec$family_int,
    link_int = glm_spec$link_int,
    A_val = A_val,
    verbose = verbose
  )
  result$nuisance_lambda_rule <- fitted_arm$nuisance_lambda_rule
  result$communication_mode <- "one_round"
  result$family <- fitted_arm$family
  result$A_val <- A_val
  result$M_tau <- fitted_arm$M_tau
  result$M_tau_inference <- fitted_arm$M_tau_inference
  result$aggregation_lambda_selection <- lambda_selection
  result$aggregation_lambda_grid <- aggregation_lambda_grid

  original_diagnostics <- fitted_arm$nuisance_fit_diagnostics
  fold_diagnostics <- lapply(seq_len(n_folds), function(k1) {
    original_row <- original_diagnostics[original_diagnostics$fold == k1, , drop = FALSE]
    if (nrow(original_row) != 1L) {
      stop(caller, ": fitted arm has inconsistent fold diagnostics.",
           call. = FALSE)
    }
    c(
      fold = k1,
      initial_outcome_degenerate =
        as.numeric(original_row$initial_outcome_degenerate),
      target_only_outcome_degenerate =
        as.numeric(original_row$target_only_outcome_degenerate),
      .summarize_source_nuisance_fit_diagnostics(
        fold_results[[k1]]$source_results
      )
    )
  })
  result$nuisance_fit_diagnostics <- as.data.frame(
    do.call(rbind, fold_diagnostics), row.names = NULL
  )
  result$rho_reuse <- list(
    reused = TRUE,
    changed_sources = changed_sources,
    reused_sources = setdiff(source_sites, changed_sources),
    refit_wall_seconds = refit_wall_seconds
  )
  result$timing <- list(
    total_seconds = .elapsed_process_seconds(started_at),
    partition_seconds = 0,
    aggregation_seconds = NA_real_,
    folds = fitted_arm$timing$folds,
    rho_reuse_refit_wall_seconds = sum(refit_wall_seconds)
  )
  result
}

.reuse_one_round_tate_across_rho <- function(
    reference_data_split, data_split, fitted_tate, changed_sources,
    lambda_selection,
    aggregation_lambda_grid = fitted_tate$aggregation_lambda_grid %||% NULL,
    verbose = FALSE, n_cores = 1L) {
  started_at <- proc.time()[["elapsed"]]
  caller <- ".reuse_one_round_tate_across_rho"
  if (!is.list(fitted_tate) || is.null(fitted_tate$arm_results) ||
      !all(c("mu1", "mu0") %in% names(fitted_tate$arm_results)) ||
      !identical(fitted_tate$communication_mode, "one_round")) {
    stop(caller, ": fitted_tate must be a reusable one-round TATE fit.",
         call. = FALSE)
  }
  mu1_result <- .refit_one_round_crossfit_sources(
    reference_data_split = reference_data_split,
    data_split = data_split,
    fitted_arm = fitted_tate$arm_results$mu1,
    changed_sources = changed_sources,
    lambda_selection = lambda_selection,
    lambda_rule = fitted_tate$aggregation_lambda_rule %||% "min",
    aggregation_lambda_grid = aggregation_lambda_grid,
    verbose = verbose,
    n_cores = n_cores
  )
  mu0_result <- fitted_tate$arm_results$mu0
  result <- calculate_tate_crossfit_aggregation(
    data_split = data_split,
    mu1_result = mu1_result,
    mu0_result = mu0_result,
    lambda_selection = lambda_selection,
    lambda_rule = fitted_tate$aggregation_lambda_rule %||% "min",
    aggregation_lambda_grid = aggregation_lambda_grid,
    verbose = verbose
  )
  result$communication_mode <- "one_round"
  result$parallel_arms <- FALSE
  result$nuisance_lambda_rule <- fitted_tate$nuisance_lambda_rule
  result$family <- fitted_tate$family
  result$M_tau <- fitted_tate$M_tau
  result$M_tau_inference <- fitted_tate$M_tau_inference
  result$aggregation_lambda_selection <- lambda_selection
  result$aggregation_lambda_grid <- aggregation_lambda_grid
  result$method <- "one_round_direct_tate_crossfit"
  result$timing <- list(
    total_seconds = .elapsed_process_seconds(started_at),
    mu1 = mu1_result$timing,
    mu0 = mu0_result$timing,
    reused_control_arm = TRUE
  )
  result$nuisance_fit_diagnostics <- rbind(
    transform(mu1_result$nuisance_fit_diagnostics, A_val = 1L),
    transform(mu0_result$nuisance_fit_diagnostics, A_val = 0L)
  )
  rownames(result$nuisance_fit_diagnostics) <- NULL
  result$rho_reuse <- list(
    reused = TRUE,
    changed_sources = unique(as.character(changed_sources)),
    reused_control_arm = TRUE
  )
  result
}

.validate_reaggregation_folds <- function(
    fitted_arm, target_folds, source_folds, source_sites, caller) {
  stored <- fitted_arm$intermediates$fold_info
  n_folds <- fitted_arm$n_folds
  if (is.null(stored) || length(stored) != n_folds ||
      length(target_folds) != n_folds ||
      any(vapply(source_folds, length, integer(1L)) != n_folds)) {
    stop(caller, ": fitted arm and deterministic fold partition do not match.",
         call. = FALSE)
  }
  for (k in seq_len(n_folds)) {
    if (!identical(stored[[k]]$target_idx,
                   target_folds[[k]]$original_idx)) {
      stop(caller, ": target evaluation-fold indices do not match the fit.",
           call. = FALSE)
    }
    for (j in seq_along(source_sites)) {
      if (!identical(stored[[k]]$source_idx[[j]],
                     source_folds[[source_sites[[j]]]][[k]]$original_idx)) {
        stop(sprintf(
          "%s: source '%s' evaluation-fold indices do not match the fit.",
          caller, source_sites[[j]]
        ), call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

.refresh_crossfit_inference_components <- function(
    fitted_arm, data_split, target_folds, source_folds,
    M_tau_inference, family_int, link_int, A_val) {
  source_sites <- setdiff(names(data_split), "t")
  fold_results <- fitted_arm$fold_results
  if (is.null(fold_results) || length(fold_results) != fitted_arm$n_folds) {
    stop(
      ".refresh_crossfit_inference_components: fitted arm lacks fold_results.",
      call. = FALSE
    )
  }

  for (k in seq_len(fitted_arm$n_folds)) {
    target_fold <- materialize_fold(target_folds, k)
    for (source_name in source_sites) {
      source_fold <- materialize_fold(source_folds[[source_name]], k)
      source_fit <- fold_results[[k]]$source_results[[source_name]]
      required <- c("gamma_s", "alpha_ts")
      missing <- required[!vapply(required, function(field) {
        !is.null(source_fit[[field]])
      }, logical(1L))]
      if (length(missing) > 0L) {
        stop(sprintf(
          paste0(
            ".refresh_crossfit_inference_components: fold %d source '%s' ",
            "lacks fitted nuisance field(s): %s."
          ),
          k, source_name, paste(missing, collapse = ", ")
        ), call. = FALSE)
      }

      correction <- calculate_correction_term_cpp(
        source_fold$Z_site,
        source_fold$A,
        source_fold$Y,
        source_fit$gamma_s,
        source_fit$alpha_ts,
        source_fold$W_outcome,
        M_tau_inference,
        family_int,
        link_int,
        A_val
      )
      mu_pred <- mean(predict_glm_cpp(
        target_fold$W_outcome,
        source_fit$alpha_ts,
        family_int,
        link_int
      ))
      source_fit$delta_ts <- correction$delta_ts
      source_fit$mu_pred_ts <- mu_pred
      source_fit$mu_ts <- mu_pred + correction$delta_ts
      source_fit$n_s <- source_fold$n
      source_fit$correction_components <- correction$correction_components
      source_fit$correction_clip_diagnostics <- correction$clip_diagnostics
      fold_results[[k]]$source_results[[source_name]] <- source_fit
    }
  }
  fold_results
}

.decorate_reaggregated_tate <- function(
    result, fitted_tate, refreshed_arms, M_tau_inference,
    lambda_selection, started_at, aggregation_lambda_grid = NULL) {
  result$communication_mode <- fitted_tate$communication_mode
  result$parallel_arms <- fitted_tate$parallel_arms %||% FALSE
  result$nuisance_lambda_rule <- fitted_tate$nuisance_lambda_rule
  result$family <- fitted_tate$family
  result$M_tau <- fitted_tate$M_tau
  result$M_tau_inference <- M_tau_inference
  result$aggregation_lambda_selection <- lambda_selection
  result$aggregation_lambda_grid <- aggregation_lambda_grid
  result$method <- paste0(
    fitted_tate$communication_mode, "_direct_tate_crossfit"
  )
  result$arm_results <- refreshed_arms
  result$nuisance_fit_diagnostics <- fitted_tate$nuisance_fit_diagnostics
  result$reused_nuisance_fits <- TRUE
  result$timing <- list(
    total_seconds = .elapsed_process_seconds(started_at),
    nuisance_refit_seconds = 0
  )
  result
}

#' Reaggregate a fitted TATE estimator without refitting nuisances
#'
#' Recomputes inference-stage correction terms, influence functions, common
#' source weights, and the final TATE variance for a new inference truncation
#' radius and/or aggregation cutoff. The fitted high-dimensional nuisance
#' coefficients are reused unchanged. This is suitable for sensitivity analyses
#' in which \code{M_tau} (the fitting radius) is fixed.
#'
#' @param data_split Site-stratified data used for the original fit.
#' @param fitted_tate Result returned by \code{run_tate_crossfit()}.
#' @param M_tau_inference Positive or infinite inference-stage truncation radius.
#' @param lambda_selection Numeric truncated-Wald multiplier, or \code{NULL} to
#'   reuse the original fixed multiplier.
#' @param lambda_rule Aggregation rule, or \code{NULL} to reuse the original.
#' @param verbose Print aggregation progress.
#' @param aggregation_lambda_grid Optional positive numeric grid used when
#'   \code{lambda_selection = "cv"}. When omitted for an original CV fit, the
#'   fitted object's grid is reused. It must be \code{NULL} for a fixed numeric
#'   multiplier.
#' @return A TATE result with refreshed inference components and the same
#'   fitted nuisance coefficients as \code{fitted_tate}.
#' @export
reaggregate_tate_crossfit <- function(
    data_split, fitted_tate, M_tau_inference,
    lambda_selection = NULL, lambda_rule = NULL, verbose = FALSE,
    aggregation_lambda_grid = NULL) {
  started_at <- proc.time()[["elapsed"]]
  caller <- "reaggregate_tate_crossfit"
  required <- c("arm_results", "n_folds", "communication_mode", "family",
                "M_tau")
  missing <- required[!vapply(required, function(field) {
    !is.null(fitted_tate[[field]])
  }, logical(1L))]
  missing_arms <- setdiff(c("mu1", "mu0"), names(fitted_tate$arm_results))
  if (length(missing) > 0L || length(missing_arms) > 0L) {
    missing_fields <- c(
      missing,
      if (length(missing_arms) > 0L) {
        paste0("arm_results$", missing_arms)
      }
    )
    stop(
      caller, ": fitted_tate lacks reusable field(s): ",
      paste(unique(missing_fields), collapse = ", "),
      ". Refit with the current run_tate_crossfit().",
      call. = FALSE
    )
  }

  validate_truncation_parameters(fitted_tate$M_tau, M_tau_inference)
  glm_spec <- resolve_glm_family(fitted_tate$family)
  communication_mode <- match.arg(
    fitted_tate$communication_mode, c("two_round", "one_round")
  )
  if (is.null(lambda_rule)) {
    lambda_rule <- fitted_tate$aggregation_lambda_rule %||% "min"
  }
  lambda_rule <- match.arg(lambda_rule, c("min", "1se"))
  if (is.null(lambda_selection)) {
    lambda_selection <- fitted_tate$aggregation_lambda_selection
    if (is.null(lambda_selection)) {
      observed <- unique(as.numeric(fitted_tate$fold_lambdas))
      if (length(observed) != 1L || !is.finite(observed)) {
        stop(
          caller,
          ": lambda_selection must be supplied when the original fit did not use one fixed multiplier.",
          call. = FALSE
        )
      }
      lambda_selection <- observed
    }
  }
  if (identical(lambda_selection, "cv") &&
      is.null(aggregation_lambda_grid)) {
    aggregation_lambda_grid <-
      fitted_tate$aggregation_lambda_grid %||% NULL
  }
  aggregation_lambda_grid <- .validate_aggregation_lambda_grid(
    lambda_selection, aggregation_lambda_grid, caller
  )
  validate_algorithm_inputs(
    data_split, lambda_selection, family = fitted_tate$family, A_val = 1L
  )

  folds <- build_crossfit_folds(data_split, fitted_tate$n_folds)
  source_sites <- setdiff(names(data_split), "t")
  refreshed_arms <- lapply(c(mu1 = 1L, mu0 = 0L), function(A_val) {
    fitted_arm <- fitted_tate$arm_results[[if (A_val == 1L) "mu1" else "mu0"]]
    .validate_reaggregation_folds(
      fitted_arm, folds$target_folds, folds$source_folds,
      source_sites, caller
    )
    refreshed_fold_results <- .refresh_crossfit_inference_components(
      fitted_arm = fitted_arm,
      data_split = data_split,
      target_folds = folds$target_folds,
      source_folds = folds$source_folds,
      M_tau_inference = M_tau_inference,
      family_int = glm_spec$family_int,
      link_int = glm_spec$link_int,
      A_val = A_val
    )
    refreshed <- .aggregate_crossfit_from_fitted_folds(
      data_split = data_split,
      target_folds = folds$target_folds,
      source_folds = folds$source_folds,
      fold_results = refreshed_fold_results,
      n_folds = fitted_tate$n_folds,
      M_tau = fitted_tate$M_tau,
      M_tau_inference = M_tau_inference,
      lambda_selection = lambda_selection,
      lambda_rule = lambda_rule,
      aggregation_lambda_grid = aggregation_lambda_grid,
      communication_mode = communication_mode,
      family_int = glm_spec$family_int,
      link_int = glm_spec$link_int,
      A_val = A_val,
      verbose = verbose
    )
    refreshed$nuisance_fit_diagnostics <- fitted_arm$nuisance_fit_diagnostics
    refreshed$nuisance_lambda_rule <- fitted_arm$nuisance_lambda_rule
    refreshed$communication_mode <- communication_mode
    refreshed$family <- fitted_tate$family
    refreshed$A_val <- A_val
    refreshed$M_tau <- fitted_tate$M_tau
    refreshed$M_tau_inference <- M_tau_inference
    refreshed$reused_nuisance_fits <- TRUE
    refreshed
  })

  result <- calculate_tate_crossfit_aggregation(
    data_split = data_split,
    mu1_result = refreshed_arms$mu1,
    mu0_result = refreshed_arms$mu0,
    lambda_selection = lambda_selection,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
    verbose = verbose
  )
  .decorate_reaggregated_tate(
    result = result,
    fitted_tate = fitted_tate,
    refreshed_arms = refreshed_arms,
    M_tau_inference = M_tau_inference,
    lambda_selection = lambda_selection,
    started_at = started_at,
    aggregation_lambda_grid = aggregation_lambda_grid
  )
}

.same_sensitivity_value <- function(x, y) {
  length(x) == 1L && length(y) == 1L && !is.na(x) && !is.na(y) &&
    ((is.infinite(x) && is.infinite(y) && sign(x) == sign(y)) ||
       (is.finite(x) && is.finite(y) && abs(x - y) <= 1e-12))
}

.sensitivity_value_key <- function(values) {
  vapply(as.numeric(values), function(value) {
    if (is.infinite(value)) {
      if (value > 0) "Inf" else "-Inf"
    } else {
      sprintf("%.17g", value)
    }
  }, character(1L))
}

#' Evaluate a cutoff/inference-radius grid from one TATE fit
#'
#' Reuses one set of fitted nuisance coefficients. Each distinct inference
#' radius refreshes correction and influence components once; all cutoffs at
#' that radius then reuse those refreshed components and rerun only the
#' low-dimensional common-weight optimization.
#'
#' @param data_split Site-stratified data used for the original fit.
#' @param fitted_tate Result returned by \code{run_tate_crossfit()}.
#' @param sensitivity_grid Data frame with numeric \code{cutoff} and
#'   \code{M_tau_inference} columns. Rows must be unique.
#' @param verbose Print reaggregation progress.
#' @return A list with the validated grid, a result list aligned to its rows,
#'   and counts of inference refreshes and nuisance refits.
#' @export
reaggregate_tate_sensitivity_grid <- function(
    data_split, fitted_tate, sensitivity_grid, verbose = FALSE) {
  required <- c("cutoff", "M_tau_inference")
  if (!is.data.frame(sensitivity_grid) || nrow(sensitivity_grid) < 1L ||
      !all(required %in% names(sensitivity_grid))) {
    stop(
      paste0(
        "reaggregate_tate_sensitivity_grid: sensitivity_grid must be a ",
        "non-empty data frame with cutoff and M_tau_inference columns."
      ),
      call. = FALSE
    )
  }
  grid <- sensitivity_grid[, required, drop = FALSE]
  grid$cutoff <- suppressWarnings(as.numeric(grid$cutoff))
  grid$M_tau_inference <- suppressWarnings(
    as.numeric(grid$M_tau_inference)
  )
  if (any(!is.finite(grid$cutoff)) || any(grid$cutoff <= 0) ||
      any(is.na(grid$M_tau_inference)) ||
      any(grid$M_tau_inference <= 0)) {
    stop(
      paste0(
        "reaggregate_tate_sensitivity_grid: cutoffs must be finite and ",
        "positive; inference radii must be positive or Inf."
      ),
      call. = FALSE
    )
  }
  key <- paste(
    .sensitivity_value_key(grid$cutoff),
    .sensitivity_value_key(grid$M_tau_inference),
    sep = "|"
  )
  if (anyDuplicated(key)) {
    stop(
      "reaggregate_tate_sensitivity_grid: grid rows must be unique.",
      call. = FALSE
    )
  }

  base_lambda <- fitted_tate$aggregation_lambda_selection
  if (!is.numeric(base_lambda) || length(base_lambda) != 1L ||
      !is.finite(base_lambda) || base_lambda <= 0) {
    stop(
      paste0(
        "reaggregate_tate_sensitivity_grid: fitted_tate must use one fixed ",
        "positive aggregation multiplier."
      ),
      call. = FALSE
    )
  }
  base_inference_radius <- fitted_tate$M_tau_inference
  validate_truncation_parameters(fitted_tate$M_tau, base_inference_radius)

  grid_radius_keys <- .sensitivity_value_key(grid$M_tau_inference)
  radius_keys <- unique(grid_radius_keys)
  inference_fits <- setNames(vector("list", length(radius_keys)), radius_keys)
  n_inference_refreshes <- 0L
  for (radius_key in radius_keys) {
    radius <- grid$M_tau_inference[
      match(radius_key, grid_radius_keys)
    ]
    if (.same_sensitivity_value(radius, base_inference_radius)) {
      inference_fits[[radius_key]] <- fitted_tate
    } else {
      inference_fits[[radius_key]] <- reaggregate_tate_crossfit(
        data_split = data_split,
        fitted_tate = fitted_tate,
        M_tau_inference = radius,
        lambda_selection = base_lambda,
        verbose = verbose
      )
      n_inference_refreshes <- n_inference_refreshes + 1L
    }
  }

  results <- vector("list", nrow(grid))
  for (index in seq_len(nrow(grid))) {
    radius_key <- grid_radius_keys[[index]]
    inference_fit <- inference_fits[[radius_key]]
    lambda <- 1 / grid$cutoff[[index]]
    if (.same_sensitivity_value(
      lambda, inference_fit$aggregation_lambda_selection
    )) {
      results[[index]] <- inference_fit
    } else {
      started_at <- proc.time()[["elapsed"]]
      reweighted <- calculate_tate_crossfit_aggregation(
        data_split = data_split,
        mu1_result = inference_fit$arm_results$mu1,
        mu0_result = inference_fit$arm_results$mu0,
        lambda_selection = lambda,
        lambda_rule = inference_fit$aggregation_lambda_rule %||% "min",
        verbose = verbose
      )
      results[[index]] <- .decorate_reaggregated_tate(
        result = reweighted,
        fitted_tate = inference_fit,
        refreshed_arms = inference_fit$arm_results,
        M_tau_inference = grid$M_tau_inference[[index]],
        lambda_selection = lambda,
        started_at = started_at
      )
    }
  }

  list(
    grid = grid,
    results = results,
    n_inference_refreshes = n_inference_refreshes,
    n_nuisance_refits = 0L
  )
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
  groups <- .validate_nuisance_cv_group_id(data$cv_group_id, data$n, "combine_folds")
  
  if (length(fold_indices) == 0) {
    # Return empty dataset with correct structure
    result <- list(
      W_outcome = matrix(0, nrow = 0, ncol = ncol(data$W_outcome)),
      Z_site = matrix(0, nrow = 0, ncol = ncol(data$Z_site)),
      A = numeric(0), Y = numeric(0), n = 0, original_idx = integer(0)
    )
    if (!is.null(groups)) result$cv_group_id <- integer(0)
    return(result)
  }
  
  if (length(fold_indices) == 1) {
    return(materialize_fold(fold_list, fold_indices[1]))
  }
  
  # Combine indices from all requested folds
  combined_idx <- do.call(c, lapply(fold_indices, function(i) fold_list[[i]]$original_idx))
  
  # Efficient path: single index-based subset of original data
  # (avoids multiple rbind operations on pre-copied fold subsets)
  result <- list(
    W_outcome = data$W_outcome[combined_idx, , drop = FALSE],
    Z_site = data$Z_site[combined_idx, , drop = FALSE],
    A = data$A[combined_idx],
    Y = data$Y[combined_idx],
    n = length(combined_idx),
    original_idx = combined_idx
  )
  if (!is.null(groups)) {
    result$cv_group_id <- groups[combined_idx]
  }
  result
}

#' Legacy aggregation-penalty variance-curve diagnostic
#'
#' Evaluates the aggregation variance over a supplied vector of truncated-Wald
#' multipliers. This is not a nuisance-model cross-validation procedure:
#' lambda affects only the low-dimensional source-weight optimization. When
#' \code{lambda_grid = NULL}, the function uses the manuscript's single fixed
#' value \code{AGG_WALD_LAMBDA = 1} (cutoff \code{1}). The helper is retained
#' for the optional oracle diagnostic; the primary cross-fitting path fixes its
#' cutoff in advance.
#'
#' \enumerate{
#'   \item \strong{Truncated-Wald multiplier.}
#'     The source-specific L1 coefficient is
#'     \eqn{(\lambda\widehat t_j-1)_+}, where
#'     \eqn{\widehat t_j} is the standardized target--source discrepancy.
#'     Consequently, \eqn{1/\lambda} is the Wald penalty-activation cutoff.
#'
#'   \item \strong{Selection rule.}
#'     This legacy helper evaluates one variance curve, not fold-level validation
#'     scores, so it cannot compute a valid one-SE rule. Use
#'     \code{select_aggregation_lambda_inner_cv()} through
#'     \code{calculate_crossfit_aggregation()} only for an explicitly requested
#'     sensitivity analysis.
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
  # STEP 1: fixed manuscript multiplier or an explicit sensitivity grid
  # ==========================================================================
  # The default helper returns only AGG_WALD_LAMBDA. A user-supplied grid is
  # allowed for diagnostics; each value maps to penalty-activation cutoff 1/lambda.
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
      cat(sprintf("  Fixed Wald multiplier: %.4f  |  cutoff: %.2f\n",
                  lambda_max, 1 / lambda_max))
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
