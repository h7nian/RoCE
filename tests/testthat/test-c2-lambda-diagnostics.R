.c2_lambda_diag_enabled <- function() {
  env_enabled <- Sys.getenv("FACEHD_RUN_C2_LAMBDA_DIAG", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-lambda", filter, fixed = TRUE)
}

.c2_lambda_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_lambda_char_env <- function(name, default, choices) {
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

.coef_attr <- function(x, name) {
  value <- attr(x, name)
  if (is.null(value)) NA_real_ else as.numeric(value)
}

.coef_support <- function(x) {
  x <- as.numeric(x)
  if (length(x) <= 1L) return(0L)
  sum(abs(x[-1L]) > 1e-8)
}

.extract_nuisance_lambdas <- function(result, rep_idx, seed, p, mode,
                                      nuisance_rule, aggregation_rule) {
  rows <- list()
  add_coef <- function(stage, coef, k1, source, k2 = NA_character_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      rep = rep_idx,
      seed = seed,
      p = p,
      mode = mode,
      nuisance_rule = nuisance_rule,
      aggregation_rule = aggregation_rule,
      k1 = k1,
      source = source,
      k2 = k2,
      stage = stage,
      lambda_used = .coef_attr(coef, "lambda_used"),
      lambda_min = .coef_attr(coef, "lambda_min"),
      lambda_1se = .coef_attr(coef, "lambda_1se"),
      lambda_rule = as.character(attr(coef, "lambda_rule") %||% NA_character_),
      support = .coef_support(coef),
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
  do.call(rbind, rows)
}

.extract_aggregation_lambdas <- function(result, rep_idx, seed, p, mode,
                                         nuisance_rule, aggregation_rule) {
  rows <- lapply(seq_along(result$fold_lambda_info), function(k1) {
    info <- result$fold_lambda_info[[k1]]
    data.frame(
      rep = rep_idx,
      seed = seed,
      p = p,
      mode = mode,
      nuisance_rule = nuisance_rule,
      aggregation_rule = aggregation_rule,
      k1 = k1,
      lambda_used = as.numeric(info$lambda_used %||% NA_real_),
      lambda_min = as.numeric(info$lambda_min %||% NA_real_),
      lambda_1se = as.numeric(info$lambda_1se %||% NA_real_),
      idx_min = as.integer(info$idx_min %||% NA_integer_),
      idx_1se = as.integer(info$idx_1se %||% NA_integer_),
      selected_rule = as.character(info$lambda_rule %||% NA_character_),
      mean_cv_score = mean(as.numeric(info$cv_scores %||% NA_real_), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.run_c2_lambda_case <- function(rep_idx, seed, p, mode, nuisance_rule,
                                aggregation_rule, n_total, K_sites,
                                n_folds, nlambda_init) {
  set.seed(seed)
  data <- generate_simulation_data(
    n_total = n_total, K = K_sites, p = p,
    config = "C2",
    estimand_type = "superpopulation",
    site_allocation = "model",
    transform_type = "mild",
    outcome_type = "binary",
    heterogeneity_type = "none",
    shift_strength = 0.5,
    dgp_type = "facehd",
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  truth <- as.numeric(data$mu1_true)

  n_cores <- .c2_lambda_int_env("FACEHD_C2_LAMBDA_CORES", 1L)

  result <- run_crossfit(
    data_split,
    n_folds = n_folds,
    communication_mode = mode,
    lambda_selection = "cv",
    lambda_rule = aggregation_rule,
    verbose = FALSE,
    n_cores = n_cores,
    nlambda_init = nlambda_init,
    family = "binomial",
    A_val = 1L,
    use_lambda_cache = TRUE,
    nuisance_lambda_rule = nuisance_rule
  )

  run_row <- data.frame(
    rep = rep_idx,
    seed = seed,
    p = p,
    mode = mode,
    nuisance_rule = nuisance_rule,
    aggregation_rule = aggregation_rule,
    truth = truth,
    estimate = result$estimate,
    bias = result$estimate - truth,
    abs_bias = abs(result$estimate - truth),
    se = result$se,
    covered = truth >= result$ci_lower && truth <= result$ci_upper,
    target_estimate = result$target_only$estimate,
    target_bias = result$target_only$estimate - truth,
    source_estimate_mean = mean(as.numeric(result$source_estimates)),
    source_estimate_min = min(as.numeric(result$source_estimates)),
    source_estimate_max = max(as.numeric(result$source_estimates)),
    weight_sum = sum(as.numeric(result$weights)),
    weight_l1 = sum(abs(as.numeric(result$weights))),
    weight_min = min(as.numeric(result$weights)),
    weight_max = max(as.numeric(result$weights)),
    aggregation_lambda_mean = mean(as.numeric(result$fold_lambdas)),
    score_mean_minus_truth = mean(as.numeric(result$all_phi_agg)) - truth,
    stringsAsFactors = FALSE
  )

  list(
    run = run_row,
    nuisance = .extract_nuisance_lambdas(
      result, rep_idx, seed, p, mode, nuisance_rule, aggregation_rule
    ),
    aggregation = .extract_aggregation_lambdas(
      result, rep_idx, seed, p, mode, nuisance_rule, aggregation_rule
    )
  )
}

.mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
.median_or_na <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

.summarise_runs <- function(df) {
  key <- with(df, paste(p, mode, nuisance_rule, aggregation_rule, sep = "|"))
  rows <- by(df, key, function(sub) {
    data.frame(
      p = sub$p[1L],
      mode = sub$mode[1L],
      nuisance_rule = sub$nuisance_rule[1L],
      aggregation_rule = sub$aggregation_rule[1L],
      n_reps = nrow(sub),
      bias_mean = mean(sub$bias),
      abs_bias_mean = mean(abs(sub$bias)),
      bias_sd = stats::sd(sub$bias),
      se_mean = mean(sub$se),
      coverage = mean(sub$covered),
      weight_sum_mean = mean(sub$weight_sum),
      weight_l1_mean = mean(sub$weight_l1),
      aggregation_lambda_mean = mean(sub$aggregation_lambda_mean),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, lapply(rows, identity))
}

.summarise_nuisance <- function(df) {
  df$lambda_1se_over_min <- df$lambda_1se / df$lambda_min
  key <- with(df, paste(p, mode, nuisance_rule, stage, sep = "|"))
  rows <- by(df, key, function(sub) {
    data.frame(
      p = sub$p[1L],
      mode = sub$mode[1L],
      nuisance_rule = sub$nuisance_rule[1L],
      stage = sub$stage[1L],
      n_rows = nrow(sub),
      lambda_used_median = .median_or_na(sub$lambda_used),
      lambda_min_median = .median_or_na(sub$lambda_min),
      lambda_1se_median = .median_or_na(sub$lambda_1se),
      lambda_1se_over_min_median = .median_or_na(sub$lambda_1se_over_min),
      support_median = .median_or_na(sub$support),
      coef_l2_median = .median_or_na(sub$coef_l2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, lapply(rows, identity))
}

test_that("c2-lambda diagnostics compare nuisance lambda rules without algorithm changes", {
  skip_if_not(.c2_lambda_diag_enabled(),
              message = "set FACEHD_RUN_C2_LAMBDA_DIAG=1 or run with --filter c2-lambda")

  n_reps <- .c2_lambda_int_env("FACEHD_C2_LAMBDA_REPS", 3L)
  n_total <- .c2_lambda_int_env("FACEHD_C2_LAMBDA_N", 2000L)
  n_folds <- .c2_lambda_int_env("FACEHD_C2_LAMBDA_FOLDS", 5L)
  nlambda_init <- .c2_lambda_int_env("FACEHD_C2_LAMBDA_NLAMBDA_INIT", 30L)
  K_sites <- .c2_lambda_int_env("FACEHD_C2_LAMBDA_K", 3L)
  p_values <- as.integer(strsplit(Sys.getenv("FACEHD_C2_LAMBDA_P", "50"),
                                  ",", fixed = TRUE)[[1L]])
  p_values <- p_values[is.finite(p_values) & p_values > 0L]
  if (length(p_values) == 0L) p_values <- 50L

  modes <- .c2_lambda_char_env(
    "FACEHD_C2_LAMBDA_MODES",
    c("one_round", "two_round"),
    c("one_round", "two_round")
  )
  nuisance_rules <- .c2_lambda_char_env(
    "FACEHD_C2_LAMBDA_NUISANCE_RULES",
    c("min", "1se"),
    c("min", "1se")
  )
  aggregation_rule <- "min"

  cat(sprintf(
    "\n[c2-lambda] Running diagnostic grid: reps=%d, n_total=%d, K=%d, p={%s}, modes={%s}, nuisance={%s}, folds=%d, nlambda_init=%d\n",
    n_reps, n_total, K_sites, paste(p_values, collapse = ","),
    paste(modes, collapse = ","), paste(nuisance_rules, collapse = ","),
    n_folds, nlambda_init
  ))
  cat("[c2-lambda] Grid compares nuisance lambda.min vs lambda.1se; aggregation lambda rule is held at min.\n")

  run_rows <- list()
  nuisance_rows <- list()
  aggregation_rows <- list()

  case_idx <- 0L
  for (p in p_values) {
    for (rep_idx in seq_len(n_reps)) {
      seed <- 73000L + 1009L * p + 37L * rep_idx
      for (mode in modes) {
        for (nuisance_rule in nuisance_rules) {
          case_idx <- case_idx + 1L
          cat(sprintf("[c2-lambda] case %02d: rep=%d p=%d mode=%s nuisance=%s seed=%d\n",
                      case_idx, rep_idx, p, mode, nuisance_rule, seed))
          out <- .run_c2_lambda_case(
            rep_idx = rep_idx,
            seed = seed,
            p = p,
            mode = mode,
            nuisance_rule = nuisance_rule,
            aggregation_rule = aggregation_rule,
            n_total = n_total,
            K_sites = K_sites,
            n_folds = n_folds,
            nlambda_init = nlambda_init
          )
          run_rows[[length(run_rows) + 1L]] <- out$run
          nuisance_rows[[length(nuisance_rows) + 1L]] <- out$nuisance
          aggregation_rows[[length(aggregation_rows) + 1L]] <- out$aggregation
        }
      }
    }
  }

  run_df <- do.call(rbind, run_rows)
  nuisance_df <- do.call(rbind, nuisance_rows)
  aggregation_df <- do.call(rbind, aggregation_rows)

  out_dir <- file.path("c2_lambda_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  run_path <- file.path(out_dir, paste0(ts, "_c2_lambda_runs.csv"))
  nuisance_path <- file.path(out_dir, paste0(ts, "_c2_lambda_nuisance.csv"))
  aggregation_path <- file.path(out_dir, paste0(ts, "_c2_lambda_aggregation.csv"))
  write.csv(run_df, run_path, row.names = FALSE)
  write.csv(nuisance_df, nuisance_path, row.names = FALSE)
  write.csv(aggregation_df, aggregation_path, row.names = FALSE)

  cat(sprintf("\n[c2-lambda] Run-level CSV: %s\n", run_path))
  cat(sprintf("[c2-lambda] Nuisance-lambda CSV: %s\n", nuisance_path))
  cat(sprintf("[c2-lambda] Aggregation-lambda CSV: %s\n", aggregation_path))

  run_summary <- .summarise_runs(run_df)
  run_summary <- run_summary[order(run_summary$p, run_summary$mode,
                                   run_summary$nuisance_rule), , drop = FALSE]
  rownames(run_summary) <- NULL
  cat("\n[c2-lambda] Run summary by nuisance lambda rule:\n")
  print(run_summary, row.names = FALSE, digits = 4)

  nuisance_summary <- .summarise_nuisance(nuisance_df)
  nuisance_summary <- nuisance_summary[order(nuisance_summary$p, nuisance_summary$mode,
                                             nuisance_summary$stage,
                                             nuisance_summary$nuisance_rule),
                                       , drop = FALSE]
  rownames(nuisance_summary) <- NULL
  cat("\n[c2-lambda] Nuisance lambda/support summary:\n")
  print(nuisance_summary, row.names = FALSE, digits = 4)

  contrast_rows <- list()
  for (p in p_values) {
    for (mode in modes) {
      sub_min <- run_df[run_df$p == p & run_df$mode == mode &
                          run_df$nuisance_rule == "min", , drop = FALSE]
      sub_1se <- run_df[run_df$p == p & run_df$mode == mode &
                          run_df$nuisance_rule == "1se", , drop = FALSE]
      matched <- merge(sub_min, sub_1se, by = c("rep", "seed", "p", "mode"),
                       suffixes = c("_min", "_1se"))
      if (nrow(matched) == 0L) next
      contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
        p = p,
        mode = mode,
        n_pairs = nrow(matched),
        mean_bias_min = mean(matched$bias_min),
        mean_bias_1se = mean(matched$bias_1se),
        mean_abs_bias_min = mean(abs(matched$bias_min)),
        mean_abs_bias_1se = mean(abs(matched$bias_1se)),
        mean_abs_bias_delta_1se_minus_min =
          mean(abs(matched$bias_1se) - abs(matched$bias_min)),
        mean_weight_l1_min = mean(matched$weight_l1_min),
        mean_weight_l1_1se = mean(matched$weight_l1_1se),
        stringsAsFactors = FALSE
      )
    }
  }
  contrast <- do.call(rbind, contrast_rows)
  cat("\n[c2-lambda] Matched contrast: nuisance lambda.1se minus lambda.min\n")
  print(contrast, row.names = FALSE, digits = 4)

  expected_runs <- length(p_values) * n_reps * length(modes) * length(nuisance_rules)
  expect_equal(nrow(run_df), expected_runs)
  expect_true(all(is.finite(run_df$estimate)))
  expect_true(all(is.finite(run_df$se)))
  expect_true(all(is.finite(run_df$bias)))
  expect_true(any(is.finite(nuisance_df$lambda_used)))
  expect_true(any(nuisance_df$nuisance_rule == "1se" &
                    nuisance_df$lambda_used >= nuisance_df$lambda_min,
                  na.rm = TRUE))
})
