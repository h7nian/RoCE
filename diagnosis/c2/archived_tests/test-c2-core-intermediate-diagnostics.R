.c2_core_intermediate_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_RUN_C2_CORE_INTERMEDIATE", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-core-intermediate", filter, fixed = TRUE)
}

.c2_core_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_core_num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || !is.finite(parsed)) default else parsed
}

.c2_core_char_vector_env <- function(name, default, choices = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  parsed <- parsed[nzchar(parsed)]
  if (!is.null(choices)) {
    bad <- setdiff(parsed, choices)
    if (length(bad) > 0L) {
      stop(sprintf("%s contains invalid value(s): %s; choices: %s",
                   name, paste(bad, collapse = ", "),
                   paste(choices, collapse = ", ")), call. = FALSE)
    }
  }
  if (length(parsed) == 0L) default else parsed
}

.c2_core_log <- function(path, state, detail = "") {
  line <- sprintf("[%s] %s%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  state, if (nzchar(detail)) paste0(" ", detail) else "")
  cat(line, "\n", sep = "")
  flush.console()
  cat(line, "\n", file = path, append = TRUE, sep = "")
}

.c2_core_append_csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(FALSE))
  write.table(df, file = path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path),
              qmethod = "double")
  invisible(TRUE)
}

.c2_core_coef_attr <- function(x, name) {
  value <- attr(x, name)
  if (is.null(value)) NA_real_ else as.numeric(value)
}

.c2_core_support <- function(x) {
  x <- as.numeric(x)
  if (length(x) <= 1L) return(0L)
  sum(abs(x[-1L]) > 1e-8)
}

.c2_core_method_row <- function(seed, method, res, truth, n_total, K_sites, p,
                                n_folds, nlambda_init, elapsed_sec) {
  estimate <- as.numeric(res$estimate)
  se <- as.numeric(res$se)
  ci_lower <- as.numeric(res$ci_lower %||% (estimate - Z_ALPHA_05 * se))
  ci_upper <- as.numeric(res$ci_upper %||% (estimate + Z_ALPHA_05 * se))
  data.frame(
    seed = seed,
    method = method,
    truth = truth,
    estimate = estimate,
    bias = estimate - truth,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    covered = truth >= ci_lower && truth <= ci_upper,
    n_total = n_total,
    K = K_sites,
    p = p,
    config = "C2",
    n_folds = n_folds,
    nlambda_init = nlambda_init,
    elapsed_sec = elapsed_sec,
    stringsAsFactors = FALSE
  )
}

.c2_core_crossfit_rows <- function(seed, method, result, truth) {
  weights <- as.numeric(result$weights)
  sources <- names(result$source_estimates)
  source_vals <- as.numeric(result$source_estimates)

  run_row <- data.frame(
    seed = seed,
    method = method,
    truth = truth,
    estimate = as.numeric(result$estimate),
    bias = as.numeric(result$estimate) - truth,
    se = as.numeric(result$se),
    target_estimate = as.numeric(result$target_only$estimate),
    target_bias = as.numeric(result$target_only$estimate) - truth,
    source_estimate_mean = mean(source_vals),
    source_estimate_min = min(source_vals),
    source_estimate_max = max(source_vals),
    weight_sum = sum(weights),
    weight_l1 = sum(abs(weights)),
    weight_min = min(weights),
    weight_max = max(weights),
    aggregation_lambda_mean = mean(as.numeric(result$fold_lambdas)),
    score_mean_minus_truth = mean(as.numeric(result$all_phi_agg)) - truth,
    n_folds = as.integer(result$n_folds),
    N_all = as.integer(result$N_all),
    nuisance_lambda_rule = as.character(result$nuisance_lambda_rule %||% NA_character_),
    aggregation_lambda_rule = as.character(result$aggregation_lambda_rule %||% NA_character_),
    stringsAsFactors = FALSE
  )

  source_row <- data.frame(
    seed = seed,
    method = method,
    source = sources,
    source_estimate = source_vals,
    source_bias = source_vals - truth,
    avg_weight = weights,
    stringsAsFactors = FALSE
  )

  fold_weight_row <- do.call(rbind, lapply(seq_len(nrow(result$fold_weights)), function(k1) {
    data.frame(
      seed = seed,
      method = method,
      k1 = k1,
      source = sources,
      weight = as.numeric(result$fold_weights[k1, ]),
      fold_lambda = as.numeric(result$fold_lambdas[k1]),
      stringsAsFactors = FALSE
    )
  }))

  fold_source <- result$intermediates$phase1$fold_source_estimates
  fold_source_row <- do.call(rbind, lapply(seq_len(nrow(fold_source)), function(k1) {
    data.frame(
      seed = seed,
      method = method,
      k1 = k1,
      source = sources,
      fold_target_estimate = as.numeric(result$intermediates$phase1$fold_target_estimate[k1]),
      fold_source_estimate = as.numeric(fold_source[k1, ]),
      fold_source_bias = as.numeric(fold_source[k1, ]) - truth,
      mu_pred_ts = as.numeric(result$intermediates$phase1$mu_pred_ts[k1, ]),
      delta_ts = as.numeric(result$intermediates$phase1$delta_ts[k1, ]),
      V_t = as.numeric(result$intermediates$phase1$V_t[k1, ]),
      V_s = as.numeric(result$intermediates$phase1$V_s[k1, ]),
      C_ot = as.numeric(result$intermediates$phase1$C_ot[k1, ]),
      stringsAsFactors = FALSE
    )
  }))

  lambda_row <- do.call(rbind, lapply(seq_along(result$fold_lambda_info), function(k1) {
    info <- result$fold_lambda_info[[k1]]
    cv_scores <- as.numeric(info$cv_scores %||% NA_real_)
    cv_se <- as.numeric(info$cv_se %||% NA_real_)
    idx_min <- as.integer(info$idx_min %||% NA_integer_)
    idx_1se <- as.integer(info$idx_1se %||% NA_integer_)
    lambda_rule <- as.character(info$lambda_rule %||% NA_character_)
    idx_used <- if (identical(lambda_rule, "1se")) idx_1se else idx_min
    data.frame(
      seed = seed,
      method = method,
      k1 = k1,
      lambda_used = as.numeric(info$lambda_used %||% NA_real_),
      lambda_min = as.numeric(info$lambda_min %||% NA_real_),
      lambda_1se = as.numeric(info$lambda_1se %||% NA_real_),
      idx_min = idx_min,
      idx_1se = idx_1se,
      idx_used = idx_used,
      selected_rule = lambda_rule,
      n_grid = length(cv_scores),
      cv_score_min = min(cv_scores, na.rm = TRUE),
      cv_score_selected = if (!is.na(idx_used) && length(cv_scores) >= idx_used) {
        cv_scores[idx_used]
      } else {
        NA_real_
      },
      cv_se_at_min = if (!is.na(idx_min) && length(cv_se) >= idx_min) {
        cv_se[idx_min]
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))

  nuisance_rows <- list()
  add_coef <- function(stage, coef, k1, source, k2 = NA_character_) {
    nuisance_rows[[length(nuisance_rows) + 1L]] <<- data.frame(
      seed = seed,
      method = method,
      k1 = k1,
      source = source,
      k2 = k2,
      stage = stage,
      lambda_used = .c2_core_coef_attr(coef, "lambda_used"),
      lambda_min = .c2_core_coef_attr(coef, "lambda_min"),
      lambda_1se = .c2_core_coef_attr(coef, "lambda_1se"),
      lambda_rule = as.character(attr(coef, "lambda_rule") %||% NA_character_),
      support = .c2_core_support(coef),
      coef_l2 = sqrt(mean(as.numeric(coef)^2)),
      stringsAsFactors = FALSE
    )
  }

  for (k1 in seq_along(result$fold_results)) {
    source_results <- result$fold_results[[k1]]$source_results
    for (source in names(source_results)) {
      src <- source_results[[source]]
      add_coef("gamma_final", src$gamma_s, k1, source)
      add_coef("alpha_final", src$alpha_ts, k1, source)
      for (k2 in names(src$per_k2_gamma)) {
        add_coef("gamma_init", src$per_k2_gamma[[k2]], k1, source, k2)
      }
      for (k2 in names(src$per_k2_alpha)) {
        add_coef("alpha_init", src$per_k2_alpha[[k2]], k1, source, k2)
      }
    }
  }

  list(
    run = run_row,
    source = source_row,
    fold_weights = fold_weight_row,
    fold_source = fold_source_row,
    aggregation_lambda = lambda_row,
    nuisance_lambda = do.call(rbind, nuisance_rows)
  )
}

.c2_core_federated_site_rows <- function(seed, data_split, truth, family, A_val = 1L) {
  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  Z_target <- as.matrix(target_data$Z_site)
  n_target <- target_data$n
  mean_phi_target <- c(1, colMeans(Z_target))
  Z_target_int <- cbind(1, Z_target)
  Z_target_centered <- sweep(Z_target_int, 2, mean_phi_target, "-")

  target_res <- calculate_weighted_site_aipw(
    y = target_data$Y, a = target_data$A,
    X = as.matrix(target_data$W_outcome), weights = NULL,
    family = family, A_val = A_val
  )
  rows <- list(data.frame(
    seed = seed,
    method = "federated_dr_site_detail",
    site = "t",
    estimate = target_res$estimate,
    bias = target_res$estimate - truth,
    site_variance = target_res$variance,
    n = target_data$n,
    dr_lambda = NA_real_,
    weight_mean = 1,
    weight_sd = 0,
    weight_min = 1,
    weight_max = 1,
    weight_ess = target_data$n,
    target_if_component_sd = NA_real_,
    source_if_component_sd = stats::sd(target_res$varphi_ot),
    stringsAsFactors = FALSE
  ))

  for (site in source_sites) {
    source_data <- data_split[[site]]
    Z_source <- as.matrix(source_data$Z_site)
    dr_weights <- calculate_dr_weights(Z_source, Z_target, lambda = NULL)
    source_res <- calculate_weighted_site_aipw(
      y = source_data$Y, a = source_data$A,
      X = as.matrix(source_data$W_outcome), weights = dr_weights,
      family = family, A_val = A_val
    )

    n_source <- source_data$n
    Z_int <- cbind(1, Z_source)
    Z_int_centered <- sweep(Z_int, 2, mean_phi_target, "-")
    d_alpha <- max(mean(dr_weights), DIVISION_FLOOR)
    phi_centered <- source_res$phi - source_res$estimate
    A_s <- -1 / d_alpha * colMeans(as.numeric(dr_weights) * Z_int_centered * as.numeric(phi_centered))
    density_score_weights <- as.numeric(dr_weights)
    M_alpha <- t(Z_int_centered) %*%
      (Z_int_centered * density_score_weights) / max(1, n_source)
    adj_alpha <- solve_with_ridge(M_alpha) %*% A_s
    infl_alpha_source <- as.numeric(
      density_score_weights * (Z_int_centered %*% adj_alpha)
    )
    target_if_component <- -as.numeric(Z_target_centered %*% adj_alpha)
    varphi_source <- as.numeric(source_res$varphi_ot) + infl_alpha_source
    varphi_source <- varphi_source - mean(varphi_source)
    source_variance <- mean(varphi_source^2) / max(1, n_source)

    rows[[length(rows) + 1L]] <- data.frame(
      seed = seed,
      method = "federated_dr_site_detail",
      site = site,
      estimate = source_res$estimate,
      bias = source_res$estimate - truth,
      site_variance = source_variance,
      n = source_data$n,
      dr_lambda = as.numeric(attr(dr_weights, "lambda_used") %||% NA_real_),
      weight_mean = mean(dr_weights),
      weight_sd = stats::sd(dr_weights),
      weight_min = min(dr_weights),
      weight_max = max(dr_weights),
      weight_ess = sum(dr_weights)^2 / sum(dr_weights^2),
      target_if_component_sd = stats::sd(target_if_component),
      source_if_component_sd = stats::sd(varphi_source),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.c2_core_component_row <- function(seed, method, res) {
  components <- res$components %||% list()
  scalar_component <- function(name) {
    x <- components[[name]]
    if (is.null(x) || length(x) != 1L) NA_real_ else as.numeric(x)
  }
  data.frame(
    seed = seed,
    method = method,
    estimate = as.numeric(res$estimate),
    se = as.numeric(res$se),
    variance = as.numeric(res$variance),
    var_target_component = scalar_component("var_target_component"),
    var_source_component = scalar_component("var_source_component"),
    target_variance_component = scalar_component("target_variance_component"),
    source_variance_component = scalar_component("source_variance_component"),
    mean_weight = scalar_component("mean_weight"),
    variance_wss = scalar_component("variance_wss"),
    ivw_weight_sum = if (!is.null(components$ivw_weights)) sum(as.numeric(components$ivw_weights)) else NA_real_,
    target_weight = scalar_component("target_weight"),
    source_weight_sum = if (!is.null(components$source_weights)) sum(as.numeric(components$source_weights)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

.c2_core_contract_rows <- function() {
  data.frame(
    item = c(
      "aggregation_objective",
      "aggregation_lambda_cv",
      "aggregation_lambda_grid",
      "two_level_calibration",
      "calibrated_gamma_cv",
      "calibrated_alpha_cv"
    ),
    expected = c(
      "eq:agg_penalized_objective TATE variance + cross-site covariance + truncated-Wald L1 penalty",
      "FACE-style inner validation: train weights on inner-training summaries and score on validation variance",
      "glmnet-style path starts at KKT lambda_max and descends geometrically",
      "SMMAL-style outer k1, secondary k2, plug-ins trained without k1/k2, fold-summed calibrated loss",
      "calibrated gamma loss uses calibrated density-ratio CV when lambda=NULL",
      "calibrated alpha loss uses calibrated outcome CV when lambda=NULL"
    ),
    implementation = c(
      "optimize_weights() includes V_ot, V_t, V_s, C_ot, C_cross and (lambda*t_j-1)_+*|eta_j| on the N_all variance scale",
      "select_aggregation_lambda_inner_cv() optimizes eta on components[-m] and scores .validation_aggregation_objective() on components[[m]]",
      ".aggregation_lambda_max() and .aggregation_lambda_grid() build the path",
      "process_source_site() loops k2 != k1, uses training_folds=setdiff(1:n_folds,c(k1,k2)), then fit_unified_* with calibrated=TRUE",
      "fit_unified_density_ratio(calibrated=TRUE) calls select_lambda_cv_calibrated_density_ratio_cpp()",
      "fit_unified_outcome(calibrated=TRUE) calls select_lambda_cv_calibrated_outcome_cpp()"
    ),
    status = "checked",
    stringsAsFactors = FALSE
  )
}

.c2_core_repo_file <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("..", ...),
    file.path("..", "..", ...)
  )
  hit <- candidates[file.exists(candidates)][1L]
  if (is.na(hit)) {
    stop(sprintf("Could not locate repository file: %s", file.path(...)),
         call. = FALSE)
  }
  hit
}

test_that("c2-core-intermediate writes single-seed intermediate estimators", {
  skip_if_not(.c2_core_intermediate_enabled(),
              message = "set ROCE_RUN_C2_CORE_INTERMEDIATE=1 or run with --filter c2-core-intermediate")

  n_total <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_N", 5000L)
  K_sites <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_K", 5L)
  p <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_P", 50L)
  n_folds <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_FOLDS", 10L)
  nlambda_init <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_NLAMBDA_INIT", 100L)
  seed <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_SEED", 240001L)
  shift_strength <- .c2_core_num_env("ROCE_C2_CORE_INTERMEDIATE_SHIFT", 0.5)
  n_cores <- .c2_core_int_env("ROCE_C2_CORE_INTERMEDIATE_CORES", 1L)
  nuisance_rule <- .c2_core_char_vector_env(
    "ROCE_C2_CORE_INTERMEDIATE_NUISANCE_RULE", "min", c("min", "1se")
  )[1L]
  aggregation_rule <- .c2_core_char_vector_env(
    "ROCE_C2_CORE_INTERMEDIATE_AGG_RULE", "min", c("min", "1se")
  )[1L]
  methods <- .c2_core_char_vector_env(
    "ROCE_C2_CORE_INTERMEDIATE_METHODS",
    c("target_only", "federated_dr", "pooled_dr", "tilted_aipw", "oracle_dr", "one_round_crossfit"),
    c("target_only", "sample_size", "inverse_variance", "federated_dr", "pooled_dr",
      "tilted_aipw", "oracle_dr", "one_round_crossfit", "two_round_crossfit")
  )

  out_dir <- file.path("c2_core_intermediate_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  slurm_id <- Sys.getenv("SLURM_JOB_ID", "local")
  run_id <- paste(
    ts,
    paste0("job", gsub("[^A-Za-z0-9_.-]", "_", slurm_id)),
    paste0("seed", seed),
    paste0("nuis", nuisance_rule),
    paste0("agg", aggregation_rule),
    sep = "_"
  )
  progress_path <- file.path(out_dir, paste0(run_id, "_c2_core_intermediate_progress.log"))
  run_path <- file.path(out_dir, paste0(run_id, "_method_runs.csv"))
  component_path <- file.path(out_dir, paste0(run_id, "_baseline_components.csv"))
  site_path <- file.path(out_dir, paste0(run_id, "_baseline_site_details.csv"))
  cf_run_path <- file.path(out_dir, paste0(run_id, "_crossfit_runs.csv"))
  cf_source_path <- file.path(out_dir, paste0(run_id, "_crossfit_sources.csv"))
  cf_weight_path <- file.path(out_dir, paste0(run_id, "_crossfit_fold_weights.csv"))
  cf_fold_path <- file.path(out_dir, paste0(run_id, "_crossfit_fold_components.csv"))
  cf_agg_path <- file.path(out_dir, paste0(run_id, "_crossfit_aggregation_lambdas.csv"))
  cf_nuisance_path <- file.path(out_dir, paste0(run_id, "_crossfit_nuisance_lambdas.csv"))
  contract_path <- file.path(out_dir, paste0(run_id, "_paper_contracts.csv"))

  .c2_core_log(
    progress_path, "config",
    sprintf("seed=%d n=%d K=%d p=%d C2 folds=%d nlambda_init=%d methods=%s nuisance_rule=%s agg_rule=%s",
            seed, n_total, K_sites, p, n_folds, nlambda_init,
            paste(methods, collapse = ","), nuisance_rule, aggregation_rule)
  )
  .c2_core_append_csv(.c2_core_contract_rows(), contract_path)

  set.seed(seed)
  data <- generate_simulation_data(
    n_total = n_total,
    K = K_sites,
    p = p,
    config = "C2",
    estimand_type = "superpopulation",
    site_allocation = "model",
    transform_type = "mild",
    outcome_type = "binary",
    heterogeneity_type = "none",
    shift_strength = shift_strength,
    dgp_type = "roce",
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  truth <- as.numeric(data$mu1_true)
  family <- switch(data$outcome_type,
    "binary" = "binomial",
    "continuous" = "gaussian",
    tolower(data$outcome_type)
  )
  .c2_core_log(
    progress_path, "data",
    sprintf("truth=%.8f mu1_realized=%.8f n_t=%d n_sources=%s",
            truth, as.numeric(data$mu1_realized), data_split$t$n,
            paste(vapply(setdiff(names(data_split), "t"),
                         function(s) data_split[[s]]$n, integer(1L)),
                  collapse = ","))
  )

  run_method <- function(method, expr) {
    .c2_core_log(progress_path, "method_start", method)
    start <- Sys.time()
    res <- force(expr)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    row <- .c2_core_method_row(seed, method, res, truth, n_total, K_sites, p,
                               n_folds, nlambda_init, elapsed)
    .c2_core_append_csv(row, run_path)
    .c2_core_append_csv(.c2_core_component_row(seed, method, res), component_path)
    .c2_core_log(
      progress_path, "method_done",
      sprintf("%s estimate=%.8f bias=%.8f se=%.8f elapsed_sec=%.1f",
              method, row$estimate, row$bias, row$se, elapsed)
    )
    res
  }

  precomputed_folds <- NULL
  target_only_ps_cache <- NULL
  if (any(methods %in% c("one_round_crossfit", "two_round_crossfit"))) {
    precomputed_folds <- build_crossfit_folds(data_split, n_folds)
    target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
  }

  results <- list()
  if ("target_only" %in% methods) {
    results$target_only <- run_method("target_only", estimate_target_only(
      data_split, family, use_rcal = FALSE, use_crossfit = TRUE,
      n_folds = n_folds, A_val = 1L
    ))
  }

  if (any(methods %in% c("sample_size", "inverse_variance"))) {
    site_fits <- .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = FALSE,
      use_crossfit = TRUE,
      n_folds = n_folds,
      A_val = 1L
    )
    if ("sample_size" %in% methods) {
      results$sample_size <- run_method("sample_size", estimate_sample_size_weighted(
        data_split, family, use_rcal = FALSE, use_crossfit = TRUE,
        n_folds = n_folds, A_val = 1L, site_fits = site_fits
      ))
    }
    if ("inverse_variance" %in% methods) {
      results$inverse_variance <- run_method("inverse_variance", estimate_inverse_variance_weighted(
        data_split, family, use_rcal = FALSE, use_crossfit = TRUE,
        n_folds = n_folds, A_val = 1L, site_fits = site_fits
      ))
    }
  }

  if ("federated_dr" %in% methods) {
    results$federated_dr <- run_method("federated_dr", estimate_federated_dr(
      data_split, dr_lambda = NULL, A_val = 1L, family = family
    ))
    .c2_core_append_csv(
      .c2_core_federated_site_rows(seed, data_split, truth, family, A_val = 1L),
      site_path
    )
  }
  if ("pooled_dr" %in% methods) {
    results$pooled_dr <- run_method("pooled_dr", estimate_pooled_dr(
      data_split, dr_lambda = NULL, A_val = 1L, family = family
    ))
  }
  if ("tilted_aipw" %in% methods) {
    results$tilted_aipw <- run_method("tilted_aipw", estimate_tilted_aipw(
      data_split, family = family, A_val = 1L
    ))
  }
  if ("oracle_dr" %in% methods) {
    results$oracle_dr <- run_method("oracle_dr", estimate_oracle_dr(
      data_split, data$alpha1_true, data$gamma_params,
      outcome_type = data$outcome_type,
      target_propensity_true = data$p_treat_true[which(data$R == "t")]
    ))
  }

  run_crossfit_method <- function(method, communication_mode) {
    .c2_core_log(progress_path, "crossfit_start", method)
    start <- Sys.time()
    res <- run_crossfit(
      data_split,
      n_folds = n_folds,
      communication_mode = communication_mode,
      lambda_selection = "cv",
      lambda_rule = aggregation_rule,
      verbose = FALSE,
      n_cores = n_cores,
      nlambda_init = nlambda_init,
      family = family,
      A_val = 1L,
      use_lambda_cache = TRUE,
      precomputed_folds = precomputed_folds,
      target_only_ps_cache = target_only_ps_cache,
      nuisance_lambda_rule = nuisance_rule
    )
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    rows <- .c2_core_crossfit_rows(seed, method, res, truth)
    rows$run$elapsed_sec <- elapsed
    .c2_core_append_csv(rows$run, cf_run_path)
    .c2_core_append_csv(rows$source, cf_source_path)
    .c2_core_append_csv(rows$fold_weights, cf_weight_path)
    .c2_core_append_csv(rows$fold_source, cf_fold_path)
    .c2_core_append_csv(rows$aggregation_lambda, cf_agg_path)
    .c2_core_append_csv(rows$nuisance_lambda, cf_nuisance_path)
    .c2_core_append_csv(.c2_core_method_row(seed, method, res, truth, n_total,
                                            K_sites, p, n_folds, nlambda_init,
                                            elapsed),
                        run_path)
    .c2_core_log(
      progress_path, "crossfit_done",
      sprintf("%s estimate=%.8f bias=%.8f se=%.8f target_bias=%.8f weight_l1=%.4f elapsed_sec=%.1f",
              method, res$estimate, res$estimate - truth, res$se,
              rows$run$target_bias, rows$run$weight_l1, elapsed)
    )
    res
  }

  if ("one_round_crossfit" %in% methods) {
    results$one_round_crossfit <- run_crossfit_method("one_round_crossfit", "one_round")
  }
  if ("two_round_crossfit" %in% methods) {
    results$two_round_crossfit <- run_crossfit_method("two_round_crossfit", "two_round")
  }

  .c2_core_log(progress_path, "output", sprintf("runs=%s", run_path))
  .c2_core_log(progress_path, "output", sprintf("contracts=%s", contract_path))
  if (file.exists(run_path)) {
    print(read.csv(run_path), row.names = FALSE, digits = 5)
  }
  if (file.exists(cf_run_path)) {
    print(read.csv(cf_run_path), row.names = FALSE, digits = 5)
  }

  expect_true(file.exists(progress_path))
  expect_true(file.exists(contract_path))
  expect_true(file.exists(run_path))
  expect_true(all(vapply(results, function(x) is.finite(as.numeric(x$estimate)), logical(1))))
})

test_that("c2-core-intermediate static contract matches FACE/SMMAL-style implementation", {
  skip_if_not(.c2_core_intermediate_enabled(),
              message = "set ROCE_RUN_C2_CORE_INTERMEDIATE=1 or run with --filter c2-core-intermediate")

  agg <- paste(readLines(.c2_core_repo_file("R", "cross_fitting_aggregation.R"),
                         warn = FALSE), collapse = "\n")
  alg <- paste(readLines(.c2_core_repo_file("R", "cross_fitting_algorithms.R"),
                         warn = FALSE), collapse = "\n")
  mf <- paste(readLines(.c2_core_repo_file("R", "model_fitting.R"),
                        warn = FALSE), collapse = "\n")

  expect_match(agg, "select_aggregation_lambda_inner_cv", fixed = TRUE)
  expect_match(agg, "components\\[-m\\]")
  expect_match(agg, "components\\[\\[m\\]\\]")
  expect_match(agg, ".validation_aggregation_objective", fixed = TRUE)
  expect_match(agg, "n_val \\* var_term")
  expect_match(agg, ".aggregation_lambda_max", fixed = TRUE)

  expect_match(alg, "secondary_folds <- setdiff\\(1:n_folds, k1\\)")
  expect_match(alg, "training_folds <- setdiff\\(1:n_folds, c\\(k1, k2\\)\\)")
  expect_match(alg, "fit_unified_density_ratio", fixed = TRUE)
  expect_match(alg, "calibrated = TRUE", fixed = TRUE)
  expect_match(alg, "fit_unified_outcome", fixed = TRUE)
  expect_match(alg, "alpha_plugin_block", fixed = TRUE)
  expect_match(alg, "gamma_plugin_block", fixed = TRUE)

  expect_match(mf, "select_lambda_cv_calibrated_density_ratio_cpp", fixed = TRUE)
  expect_match(mf, "select_lambda_cv_calibrated_outcome_cpp", fixed = TRUE)
})
