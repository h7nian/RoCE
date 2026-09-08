# Diagnosis-only monkey patch for density-ratio CV scaling.
#
# The patch replaces only the R wrappers used to choose density-ratio lambdas
# in the current R session. Final fixed-lambda fitting still calls RoCE's
# compiled fit_*_density_ratio_cpp functions. No package source files are
# modified.

c2_patch_ns_get <- function(name) {
  ns <- asNamespace("RoCE")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(sprintf("RoCE namespace does not contain '%s'.", name), call. = FALSE)
  }
  get(name, envir = ns, inherits = FALSE)
}

c2_replace_namespace_binding <- function(name, value) {
  ns <- asNamespace("RoCE")
  was_locked <- bindingIsLocked(name, ns)
  if (was_locked) unlockBinding(name, ns)
  assign(name, value, envir = ns)
  if (was_locked) lockBinding(name, ns)
  invisible(value)
}

install_c2_dr_cv_scale_patch <- function(
    train_scale = c("production", "arm_fraction"),
    cpp_path = file.path("diagnosis", "c2", "c2_dr_cv_scale_patch.cpp")) {
  train_scale <- match.arg(train_scale)

  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("install_c2_dr_cv_scale_patch requires Rcpp.", call. = FALSE)
  }

  cpp_env <- new.env(parent = baseenv())
  Rcpp::sourceCpp(cpp_path, env = cpp_env, rebuild = FALSE, verbose = FALSE)

  ns <- asNamespace("RoCE")
  original_initial <- get("fit_initial_density_ratio", envir = ns, inherits = FALSE)
  original_unified <- get("fit_unified_density_ratio", envir = ns, inherits = FALSE)

  ns_get <- c2_patch_ns_get

  patched_initial <- function(Z_site, A, mean_phi, lambda = NULL,
                              max_iter = ns_get("MAX_ITER_DEFAULT"),
                              tol = ns_get("TOL_DEFAULT"),
                              A_val = 1L, warm_start = NULL,
                              lambda_rule = c("min", "1se")) {
    if (!is.null(lambda)) {
      return(original_initial(
        Z_site = Z_site, A = A, mean_phi = mean_phi, lambda = lambda,
        max_iter = max_iter, tol = tol, A_val = A_val,
        warm_start = warm_start, lambda_rule = lambda_rule
      ))
    }

    A_val <- ns_get(".validate_A_val")(A_val, "fit_initial_density_ratio")
    lambda_rule <- ns_get(".match_nuisance_lambda_rule")(
      lambda_rule, "fit_initial_density_ratio"
    )

    n_cv_folds <- ns_get(".nuisance_cv_fold_count")(
      A, A_val, "fit_initial_density_ratio"
    )
    lmax <- ns_get("compute_lambda_max_initial_dr")(Z_site, A, mean_phi, A_val)
    lambda_min_ratio <- if (nrow(Z_site) > ncol(Z_site)) {
      ns_get("LAMBDA_MIN_RATIO_LOW_DIM")
    } else {
      ns_get("LAMBDA_MIN_RATIO_HIGH_DIM")
    }
    lambda_grid <- ns_get("build_lambda_grid")(
      lambda_max = lmax, lambda_min_ratio = lambda_min_ratio
    )

    cv_result <- cpp_env$c2_select_lambda_cv_initial_density_ratio_scaled_cpp(
      Z_site, A, mean_phi, lambda_grid, n_cv_folds, max_iter, tol, A_val,
      train_scale
    )
    lambda <- ns_get(".select_nuisance_cv_lambda")(
      cv_result, lambda_rule, "fit_initial_density_ratio"
    )

    ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
    cpp_result <- ns_get("fit_initial_density_ratio_cpp")(
      Z_site, A, mean_phi, lambda, max_iter, tol, A_val, ws
    )
    result <- cpp_result$gamma
    attr(lambda, "cv_source_scale") <- cv_result$cv_source_scale
    attr(lambda, "cv_train_scale") <- cv_result$cv_train_scale
    ns_get(".attach_nuisance_lambda_attrs")(result, lambda)
  }

  patched_unified <- function(Z_site, A, mean_grad_psi, alpha_init,
                              lambda = NULL,
                              max_iter = ns_get("MAX_ITER_DEFAULT"),
                              tol = ns_get("TOL_DEFAULT"),
                              calibrated = FALSE,
                              M_tau = ns_get("M_TAU_DEFAULT"),
                              W_outcome = NULL, A_val = 1L,
                              family_int = 1L, link_int = 1L,
                              warm_start = NULL,
                              lambda_rule = c("min", "1se")) {
    if (!is.null(lambda)) {
      return(original_unified(
        Z_site = Z_site, A = A, mean_grad_psi = mean_grad_psi,
        alpha_init = alpha_init, lambda = lambda, max_iter = max_iter,
        tol = tol, calibrated = calibrated, M_tau = M_tau,
        W_outcome = W_outcome, A_val = A_val,
        family_int = family_int, link_int = link_int,
        warm_start = warm_start, lambda_rule = lambda_rule
      ))
    }

    A_val <- ns_get(".validate_A_val")(A_val, "fit_unified_density_ratio")
    lambda_rule <- ns_get(".match_nuisance_lambda_rule")(
      lambda_rule, "fit_unified_density_ratio"
    )
    if (is.null(W_outcome)) {
      stop("fit_unified_density_ratio: W_outcome must be provided explicitly.",
           call. = FALSE)
    }
    W_outcome_use <- as.matrix(W_outcome)
    if (nrow(W_outcome_use) != nrow(Z_site)) {
      stop(sprintf(
        "fit_unified_density_ratio: row mismatch between Z_site (%d) and W_outcome (%d).",
        nrow(Z_site), nrow(W_outcome_use)
      ), call. = FALSE)
    }

    n_cv_folds <- ns_get(".nuisance_cv_fold_count")(
      A, A_val, "fit_unified_density_ratio"
    )
    lmax <- ns_get("compute_lambda_max_refined_dr")(
      Z_site, A, mean_grad_psi, alpha_init,
      A_val = A_val, family_int = family_int, link_int = link_int,
      W_outcome = W_outcome_use, calibrated = calibrated, M_tau = M_tau
    )
    lambda_min_ratio <- if (nrow(Z_site) > ncol(Z_site)) {
      ns_get("LAMBDA_MIN_RATIO_LOW_DIM")
    } else {
      ns_get("LAMBDA_MIN_RATIO_HIGH_DIM")
    }
    lambda_grid <- ns_get("build_lambda_grid")(
      lambda_max = lmax, lambda_min_ratio = lambda_min_ratio
    )

    cv_M_tau <- if (isTRUE(calibrated)) M_tau else Inf
    cv_result <- cpp_env$c2_select_lambda_cv_calibrated_density_ratio_scaled_cpp(
      Z_site, A, mean_grad_psi, alpha_init, lambda_grid, n_cv_folds,
      max_iter, tol, cv_M_tau, W_outcome_use, A_val, family_int, link_int,
      train_scale
    )
    lambda <- ns_get(".select_nuisance_cv_lambda")(
      cv_result, lambda_rule, "fit_unified_density_ratio"
    )

    ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
    cpp_result <- ns_get("fit_unified_density_ratio_cpp")(
      Z_site, A, mean_grad_psi, alpha_init, lambda, max_iter, tol,
      calibrated, M_tau, W_outcome_use, A_val, family_int, link_int, ws
    )
    result <- cpp_result$gamma
    attr(lambda, "cv_source_scale") <- cv_result$cv_source_scale
    attr(lambda, "cv_train_scale") <- cv_result$cv_train_scale
    ns_get(".attach_nuisance_lambda_attrs")(result, lambda)
  }

  c2_replace_namespace_binding("fit_initial_density_ratio", patched_initial)
  c2_replace_namespace_binding("fit_unified_density_ratio", patched_unified)

  list(
    train_scale = train_scale,
    cpp_path = normalizePath(cpp_path, mustWork = FALSE),
    original_initial = original_initial,
    original_unified = original_unified
  )
}
