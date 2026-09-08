roce_make_tate_result_row <- function(
    task, fit, method, truth, cutoff, M_tau_inference,
    primary_cutoff, n_total, nlambda_init, n_bootstrap, outcome_family) {
  outcome_family <- match.arg(outcome_family, c("binomial", "gaussian"))
  estimate <- as.numeric(fit$estimate)
  se <- as.numeric(fit$se)
  z <- stats::qnorm(0.975)
  ci_lower <- estimate - z * se
  ci_upper <- estimate + z * se
  heterogeneity_type <- RoCE:::.face_heterogeneity_type(
    ate_deviation = as.numeric(task$rho),
    n_deviated_sites = as.integer(as.numeric(task$rho) > 0)
  )
  result <- RoCE:::.make_simulation_result_row(
    sim_id = as.integer(task$sim_id),
    method = method,
    estimate = estimate,
    se = se,
    truth = truth,
    n_total = n_total,
    K = as.integer(task$K),
    p = as.integer(task$p),
    config = as.character(task$config),
    heterogeneity_type = heterogeneity_type,
    estimand_type = "superpopulation"
  )
  result$estimand_scope <- "tate"
  result$dgp_type <- "face"
  result$outcome_family <- outcome_family
  result$truth <- truth
  result$ci_lower <- ci_lower
  result$ci_upper <- ci_upper
  result$aggregation_lambda <- 1 / cutoff
  result$aggregation_cutoff <- cutoff
  result$primary_cutoff <- primary_cutoff
  result$nlambda_init <- nlambda_init
  result$n_bootstrap <- n_bootstrap
  result$M_tau <- as.numeric(fit$M_tau)
  result$M_tau_inference <- M_tau_inference
  result$task_id <- as.integer(task$task_id)
  result$experiment <- "c3_reused_sensitivity"
  result$primary_experiment <- as.character(task$experiment)
  result$rho <- as.numeric(task$rho)
  result$cutoff <- cutoff
  result$n_site <- as.integer(task$n_site)
  result$n_folds <- as.integer(task$n_folds)
  result$reused_primary_task <- TRUE
  result$nuisance_refit_count <- 0L
  result
}

roce_build_reused_sensitivity_rows <- function(
    task, artifacts, sensitivity_grid, nlambda_init, n_bootstrap,
    outcome_family = c("binomial", "gaussian")) {
  outcome_family <- match.arg(outcome_family)
  if (is.null(artifacts$data_split) ||
      is.null(artifacts$direct_tate_results$one_round_crossfit) ||
      !is.numeric(artifacts$tate_truth) ||
      length(artifacts$tate_truth) != 1L ||
      !is.finite(artifacts$tate_truth)) {
    stop(
      "reused sensitivity requires complete one-round simulation artifacts.",
      call. = FALSE
    )
  }
  fitted_tate <- artifacts$direct_tate_results$one_round_crossfit
  primary_cutoff <- 1 / as.numeric(fitted_tate$aggregation_lambda_selection)
  if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
      primary_cutoff <= 0) {
    stop("fitted TATE artifacts have no valid primary cutoff.",
         call. = FALSE)
  }
  evaluated <- reaggregate_tate_sensitivity_grid(
    data_split = artifacts$data_split,
    fitted_tate = fitted_tate,
    sensitivity_grid = sensitivity_grid,
    verbose = FALSE
  )
  if (evaluated$n_nuisance_refits != 0L ||
      length(evaluated$results) != nrow(evaluated$grid)) {
    stop("sensitivity reaggregation unexpectedly refitted nuisances.",
         call. = FALSE)
  }

  n_total <- sum(vapply(
    artifacts$data_split, function(site) site$n, integer(1L)
  ))
  crossfit_diagnostics <- c(
    RoCE:::.summarize_tate_crossfit_timing(
      fitted_tate$arm_results$mu1,
      fitted_tate$arm_results$mu0
    ),
    RoCE:::.summarize_tate_nuisance_diagnostics(
      fitted_tate$arm_results$mu1,
      fitted_tate$arm_results$mu0
    )
  )
  cell_diagnostics <- RoCE:::.binary_cell_diagnostics(
    artifacts$data_split, outcome_family
  )

  rows <- vector("list", 2L * nrow(evaluated$grid))
  for (index in seq_len(nrow(evaluated$grid))) {
    fit <- evaluated$results[[index]]
    cutoff <- evaluated$grid$cutoff[[index]]
    M_tau_inference <- evaluated$grid$M_tau_inference[[index]]
    face <- roce_make_tate_result_row(
      task = task,
      fit = fit,
      method = "one_round_crossfit_ate",
      truth = artifacts$tate_truth,
      cutoff = cutoff,
      M_tau_inference = M_tau_inference,
      primary_cutoff = primary_cutoff,
      n_total = n_total,
      nlambda_init = nlambda_init,
      n_bootstrap = n_bootstrap,
      outcome_family = outcome_family
    )
    aggregation_diagnostics <-
      RoCE:::.summarize_tate_aggregation_diagnostics(fit)
    face_diagnostics <- c(aggregation_diagnostics, crossfit_diagnostics)
    if (anyDuplicated(names(face_diagnostics))) {
      stop("sensitivity diagnostic names must be unique.", call. = FALSE)
    }
    for (name in names(face_diagnostics)) {
      face[[name]] <- face_diagnostics[[name]]
    }
    for (name in names(cell_diagnostics)) {
      face[[name]] <- cell_diagnostics[[name]]
    }
    face$nuisance_lambda_rule <- fitted_tate$nuisance_lambda_rule
    face$inference_refresh_count <- evaluated$n_inference_refreshes

    target <- roce_make_tate_result_row(
      task = task,
      fit = list(
        estimate = fit$target_only$estimate,
        se = fit$target_only$se,
        M_tau = fit$M_tau
      ),
      method = "target_only_ate",
      truth = artifacts$tate_truth,
      cutoff = cutoff,
      M_tau_inference = M_tau_inference,
      primary_cutoff = primary_cutoff,
      n_total = n_total,
      nlambda_init = nlambda_init,
      n_bootstrap = n_bootstrap,
      outcome_family = outcome_family
    )
    for (name in names(face_diagnostics)) {
      target[[name]] <- NA_real_
    }
    for (name in names(cell_diagnostics)) {
      target[[name]] <- cell_diagnostics[[name]]
    }
    target$nuisance_lambda_rule <- fitted_tate$nuisance_lambda_rule
    target$target_anchor_weight <- 1
    target$mean_abs_source_weight <- 0
    target$max_abs_source_weight <- 0
    target$inference_refresh_count <- evaluated$n_inference_refreshes

    rows[[2L * index - 1L]] <- face
    rows[[2L * index]] <- target
  }
  result <- RoCE:::.bind_sim_result_list(rows)
  result$sensitivity_kind <- ifelse(
    result$M_tau_inference != fitted_tate$M_tau_inference,
    "inference_radius",
    ifelse(
      abs(result$cutoff - primary_cutoff) >
        1e-12,
      "aggregation_cutoff",
      "reference_identity"
    )
  )
  reference_rows <- result$sensitivity_kind == "reference_identity"
  if (sum(reference_rows) != 2L ||
      !setequal(
        result$method[reference_rows],
        c("one_round_crossfit_ate", "target_only_ate")
      )) {
    stop(
      "sensitivity grid must contain exactly one two-method primary identity.",
      call. = FALSE
    )
  }
  result
}
