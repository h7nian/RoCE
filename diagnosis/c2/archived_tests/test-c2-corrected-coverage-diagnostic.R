.c2_corrected_coverage_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_RUN_C2_CORRECTED_COVERAGE", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-corrected-coverage", filter, fixed = TRUE)
}

.c2_corrected_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_corrected_num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || !is.finite(parsed)) default else parsed
}

.c2_corrected_char_env <- function(name, default, choices) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  parsed <- parsed[nzchar(parsed)]
  bad <- setdiff(parsed, choices)
  if (length(bad) > 0L) {
    stop(sprintf("%s must only contain: %s; got invalid value(s): %s",
                 name, paste(choices, collapse = ", "), paste(bad, collapse = ", ")),
         call. = FALSE)
  }
  if (length(parsed) == 0L) default else parsed
}

.c2_corrected_summarise <- function(results) {
  split_key <- with(results, paste(n_total, K, p, config, method, sep = "|"))
  rows <- by(results, split_key, function(sub) {
    bias_sd <- stats::sd(sub$bias)
    se_mean <- mean(sub$se)
    bias_mean <- mean(sub$bias)
    normal_with_bias <- NA_real_
    normal_centered <- NA_real_
    if (is.finite(bias_sd) && bias_sd > 0 && is.finite(se_mean)) {
      half_width <- Z_ALPHA_05 * se_mean
      normal_with_bias <- stats::pnorm((half_width - bias_mean) / bias_sd) -
        stats::pnorm((-half_width - bias_mean) / bias_sd)
      normal_centered <- stats::pnorm(half_width / bias_sd) -
        stats::pnorm(-half_width / bias_sd)
    }
    data.frame(
      n_total = sub$n_total[1L],
      K = sub$K[1L],
      p = sub$p[1L],
      config = sub$config[1L],
      method = sub$method[1L],
      n_success = nrow(sub),
      bias_mean = bias_mean,
      bias_sd = bias_sd,
      rmse = sqrt(mean(sub$bias^2)),
      se_mean = se_mean,
      coverage = mean(sub$coverage),
      ci_width_mean = mean(sub$ci_width),
      bias_over_se = abs(bias_mean) / se_mean,
      se_over_emp_sd = se_mean / bias_sd,
      normal_coverage_with_bias = normal_with_bias,
      normal_coverage_centered = normal_centered,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, lapply(rows, identity))
  out[order(out$p, out$coverage, out$method), , drop = FALSE]
}

.c2_corrected_row <- function(sim_id, method, res, truth, job,
                              estimand_type = "superpopulation",
                              heterogeneity_type = "none",
                              dgp_type = "roce") {
  estimate <- as.numeric(res$estimate)
  se <- as.numeric(res$se)
  ci_lower <- as.numeric(res$ci_lower %||% (estimate - Z_ALPHA_05 * se))
  ci_upper <- as.numeric(res$ci_upper %||% (estimate + Z_ALPHA_05 * se))
  data.frame(
    sim_id = sim_id,
    method = method,
    estimate = estimate,
    se = se,
    bias = estimate - truth,
    coverage = (truth >= ci_lower) & (truth <= ci_upper),
    ci_width = ci_upper - ci_lower,
    n_total = job$n_total,
    K = job$K,
    p = job$p,
    config = "C2",
    heterogeneity_type = heterogeneity_type,
    estimand_type = estimand_type,
    dgp_type = dgp_type,
    stringsAsFactors = FALSE
  )
}

.c2_corrected_run_one_direct <- function(job, progress_line) {
  crossfit_verbose <- Sys.getenv("ROCE_C2_CORRECTED_COVERAGE_CROSSFIT_VERBOSE", "0") %in%
    c("1", "TRUE", "true", "True")
  set.seed(job$seed)
  data <- generate_simulation_data(
    n_total = job$n_total,
    K = job$K,
    p = job$p,
    config = "C2",
    estimand_type = "superpopulation",
    site_allocation = "model",
    transform_type = "mild",
    outcome_type = "binary",
    heterogeneity_type = "none",
    shift_strength = job$shift_strength,
    dgp_type = "roce",
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  truth <- as.numeric(data$mu1_true)
  methods <- job$methods[[1L]]
  family <- switch(data$outcome_type,
    "binary" = "binomial",
    "continuous" = "gaussian",
    tolower(data$outcome_type)
  )

  rows <- list()
  add_method <- function(method, expr) {
    progress_line("method_start", sprintf(" method=%s", method))
    method_start <- Sys.time()
    res <- force(expr)
    elapsed <- as.numeric(difftime(Sys.time(), method_start, units = "secs"))
    progress_line(
      "method_done",
      sprintf(" method=%s estimate=%.6f se=%.6f elapsed_sec=%.1f",
              method, as.numeric(res$estimate), as.numeric(res$se), elapsed)
    )
    rows[[length(rows) + 1L]] <<- .c2_corrected_row(
      sim_id = job$seed,
      method = method,
      res = res,
      truth = truth,
      job = job
    )
  }

  precomputed_folds <- NULL
  target_only_ps_cache <- NULL
  if (any(methods %in% c("one_round_crossfit", "two_round_crossfit"))) {
    precomputed_folds <- build_crossfit_folds(data_split, job$n_folds)
    target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
  }

  if ("one_round_crossfit" %in% methods) {
    add_method("one_round_crossfit", run_crossfit(
      data_split,
      n_folds = job$n_folds,
      communication_mode = "one_round",
      lambda_selection = "cv",
      verbose = crossfit_verbose,
      n_cores = 1L,
      nlambda_init = job$nlambda_init,
      family = family,
      A_val = 1L,
      use_lambda_cache = TRUE,
      precomputed_folds = precomputed_folds,
      target_only_ps_cache = target_only_ps_cache
    ))
  }

  if ("two_round_crossfit" %in% methods) {
    add_method("two_round_crossfit", run_crossfit(
      data_split,
      n_folds = job$n_folds,
      communication_mode = "two_round",
      lambda_selection = "cv",
      verbose = crossfit_verbose,
      n_cores = 1L,
      nlambda_init = job$nlambda_init,
      family = family,
      A_val = 1L,
      use_lambda_cache = TRUE,
      precomputed_folds = precomputed_folds,
      target_only_ps_cache = target_only_ps_cache
    ))
  }

  if ("target_only" %in% methods) {
    add_method("target_only", estimate_target_only(
      data_split, family, use_rcal = FALSE, use_crossfit = TRUE,
      n_folds = job$n_folds, A_val = 1L
    ))
  }

  if (any(methods %in% c("sample_size", "inverse_variance"))) {
    site_fits <- .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = FALSE,
      use_crossfit = TRUE,
      n_folds = job$n_folds,
      A_val = 1L
    )
    if ("sample_size" %in% methods) {
      add_method("sample_size", estimate_sample_size_weighted(
        data_split, family, use_rcal = FALSE, use_crossfit = TRUE,
        n_folds = job$n_folds, A_val = 1L, site_fits = site_fits
      ))
    }
    if ("inverse_variance" %in% methods) {
      add_method("inverse_variance", estimate_inverse_variance_weighted(
        data_split, family, use_rcal = FALSE, use_crossfit = TRUE,
        n_folds = job$n_folds, A_val = 1L, site_fits = site_fits
      ))
    }
  }

  if ("federated_dr" %in% methods) {
    add_method("federated_dr", estimate_federated_dr(
      data_split, dr_lambda = NULL, A_val = 1L, family = family
    ))
  }
  if ("pooled_dr" %in% methods) {
    add_method("pooled_dr", estimate_pooled_dr(
      data_split, dr_lambda = NULL, A_val = 1L, family = family
    ))
  }
  if ("tilted_aipw" %in% methods) {
    add_method("tilted_aipw", estimate_tilted_aipw(
      data_split, family = family, A_val = 1L
    ))
  }
  if ("oracle_dr" %in% methods) {
    target_idx <- which(data$R == "t")
    target_propensity_true <- NULL
    if (!is.null(data$p_treat_true)) {
      target_propensity_true <- data$p_treat_true[target_idx]
    }
    add_method("oracle_dr", estimate_oracle_dr(
      data_split,
      data$alpha1_true,
      data$gamma_params,
      outcome_type = data$outcome_type,
      target_propensity_true = target_propensity_true
    ))
  }

  do.call(rbind, rows)
}

.c2_corrected_run_one <- function(job, job_index = NA_integer_,
                                  total_jobs = NA_integer_,
                                  progress_path = NULL) {
  progress_line <- function(state, extra = "") {
    line <- sprintf(
      "[%s] %s job=%s/%s rep=%s seed=%s n=%s K=%s p=%s%s",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      state,
      as.character(job_index),
      as.character(total_jobs),
      as.character(job$rep),
      as.character(job$seed),
      as.character(job$n_total),
      as.character(job$K),
      as.character(job$p),
      extra
    )
    cat(line, "\n", sep = "")
    flush.console()
    if (!is.null(progress_path)) {
      cat(line, "\n", file = progress_path, append = TRUE, sep = "")
    }
  }

  progress_line("start")
  start_time <- Sys.time()
  runner <- Sys.getenv("ROCE_C2_CORRECTED_COVERAGE_RUNNER", "single")
  out <- if (identical(runner, "direct")) {
    .c2_corrected_run_one_direct(job, progress_line = progress_line)
  } else {
    run_single_simulation(
      sim_id = job$seed,
      n_total = job$n_total,
      K = job$K,
      p = job$p,
      config = "C2",
      methods = job$methods[[1L]],
      verbose = FALSE,
      n_cores_internal = 1L,
      nlambda_init = job$nlambda_init,
      estimand_type = "superpopulation",
      site_allocation = "model",
      transform_type = "mild",
      outcome_type = "binary",
      heterogeneity_type = "none",
      shift_strength = job$shift_strength,
      n_folds = job$n_folds,
      use_lambda_cache = TRUE,
      estimate_ate = FALSE,
      dgp_type = "roce"
    )
  }
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  progress_line("done", sprintf(" runner=%s rows=%d elapsed_sec=%.1f",
                                runner, nrow(out), elapsed))
  out
}

test_that("c2-corrected-coverage rechecks coverage after source-label correction", {
  skip_if_not(.c2_corrected_coverage_enabled(),
              message = "set ROCE_RUN_C2_CORRECTED_COVERAGE=1 or run with --filter c2-corrected-coverage")

  n_reps <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_REPS", 3L)
  n_total <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_N", 2000L)
  K_sites <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_K", 3L)
  n_folds <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_FOLDS", 3L)
  nlambda_init <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_NLAMBDA_INIT", 20L)
  n_cores <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_CORES", 1L)
  base_seed <- .c2_corrected_int_env("ROCE_C2_CORRECTED_COVERAGE_SEED", 88000L)
  shift_strength <- .c2_corrected_num_env("ROCE_C2_CORRECTED_COVERAGE_SHIFT", 0.5)

  p_values <- as.integer(strsplit(
    Sys.getenv("ROCE_C2_CORRECTED_COVERAGE_P", "10,50"),
    ",", fixed = TRUE
  )[[1L]])
  p_values <- p_values[is.finite(p_values) & p_values > 0L]
  if (length(p_values) == 0L) p_values <- c(10L, 50L)

  valid_methods <- c(
    "one_round_crossfit", "two_round_crossfit", "target_only",
    "sample_size", "inverse_variance", "federated_dr", "pooled_dr",
    "tilted_aipw", "oracle_dr"
  )
  methods <- .c2_corrected_char_env(
    "ROCE_C2_CORRECTED_COVERAGE_METHODS",
    c("one_round_crossfit", "two_round_crossfit", "federated_dr",
      "pooled_dr", "tilted_aipw", "oracle_dr"),
    valid_methods
  )

  jobs <- do.call(rbind, lapply(seq_along(p_values), function(p_idx) {
    p <- p_values[[p_idx]]
    data.frame(
      rep = seq_len(n_reps),
      seed = base_seed + 100000L * p_idx + seq_len(n_reps),
      n_total = n_total,
      K = K_sites,
      p = p,
      n_folds = n_folds,
      nlambda_init = nlambda_init,
      shift_strength = shift_strength,
      stringsAsFactors = FALSE
    )
  }))
  jobs$methods <- I(rep(list(methods), nrow(jobs)))

  out_dir <- file.path("c2_corrected_coverage_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  progress_path <- file.path(out_dir, paste0(ts, "_c2_corrected_coverage_progress.log"))

  cat(sprintf(
    "\n[c2-corrected-coverage] reps=%d, n_total=%d, K=%d, p={%s}, folds=%d, nlambda_init=%d, cores=%d, methods={%s}\n",
    n_reps, n_total, K_sites, paste(p_values, collapse = ","),
    n_folds, nlambda_init, n_cores, paste(methods, collapse = ",")
  ))

  run_list <- if (n_cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(
      seq_len(nrow(jobs)),
      function(i) .c2_corrected_run_one(
        jobs[i, , drop = FALSE],
        job_index = i,
        total_jobs = nrow(jobs),
        progress_path = progress_path
      ),
      mc.cores = min(n_cores, nrow(jobs)),
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(nrow(jobs)),
           function(i) .c2_corrected_run_one(
             jobs[i, , drop = FALSE],
             job_index = i,
             total_jobs = nrow(jobs),
             progress_path = progress_path
           ))
  }

  results <- do.call(rbind, run_list)
  summary <- .c2_corrected_summarise(results)

  run_path <- file.path(out_dir, paste0(ts, "_c2_corrected_coverage_runs.csv"))
  summary_path <- file.path(out_dir, paste0(ts, "_c2_corrected_coverage_summary.csv"))
  write.csv(results, run_path, row.names = FALSE)
  write.csv(summary, summary_path, row.names = FALSE)

  cat(sprintf("\n[c2-corrected-coverage] Run-level CSV: %s\n", run_path))
  cat(sprintf("[c2-corrected-coverage] Summary CSV: %s\n", summary_path))
  cat(sprintf("[c2-corrected-coverage] Progress log: %s\n", progress_path))
  cat("\n[c2-corrected-coverage] Summary sorted by p and coverage:\n")
  print(summary, row.names = FALSE, digits = 4)

  expect_equal(nrow(results), nrow(jobs) * length(methods))
  expect_true(all(is.finite(results$estimate)))
  expect_true(all(is.finite(results$se)))
  expect_true(all(results$se >= 0))
  expect_true(all(is.finite(results$bias)))
  expect_true(all(results$method %in% methods))
  expect_true(all(summary$n_success == n_reps))

  focus <- summary[summary$method %in%
                     c("one_round_crossfit", "two_round_crossfit", "federated_dr"),
                   , drop = FALSE]
  if (n_reps >= 10L) {
    expect_true(all(focus$bias_over_se < 1.5, na.rm = TRUE),
                info = "Corrected C2 diagnostics should not show bias above 1.5 mean SE in focused methods.")
  }
})
