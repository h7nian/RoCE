.c2_ps_basis_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_RUN_C2_TARGET_PS_BASIS", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-target-ps-basis", filter, fixed = TRUE)
}

.c2_ps_basis_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_ps_basis_num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || !is.finite(parsed)) default else parsed
}

.c2_ps_basis_log <- function(path, state, detail = "") {
  line <- sprintf("[%s] %s%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  state, if (nzchar(detail)) paste0(" ", detail) else "")
  cat(line, "\n", sep = "")
  flush.console()
  cat(line, "\n", file = path, append = TRUE, sep = "")
}

.c2_ps_basis_append_csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(FALSE))
  write.table(df, file = path, sep = ",", row.names = FALSE,
              col.names = !file.exists(path), append = file.exists(path),
              qmethod = "double")
  invisible(TRUE)
}

.target_only_ps_basis_crossfit <- function(target_data,
                                           ps_basis = c("W", "Z", "oracle"),
                                           or_basis = c("W", "oracle"),
                                           n_folds = 10L, family = "binomial",
                                           A_val = 1L,
                                           alpha1_true = NULL,
                                           alpha0_true = NULL) {
  ps_basis <- match.arg(ps_basis)
  or_basis <- match.arg(or_basis)
  glm_spec <- resolve_glm_family(family)

  y <- target_data$Y
  a <- target_data$A
  x_or <- as.matrix(target_data$W_outcome)
  x_ps <- switch(ps_basis,
    W = as.matrix(target_data$W_outcome),
    Z = as.matrix(target_data$Z_site),
    oracle = NULL
  )
  n <- target_data$n

  oracle_m_pred <- NULL
  if (identical(or_basis, "oracle")) {
    alpha_true <- if (A_val == 1L) alpha1_true else alpha0_true
    if (is.null(alpha_true)) {
      stop("target-only PS-basis diagnostic: oracle OR requires true alpha.",
           call. = FALSE)
    }
    w_true <- if ("W_outcome_true" %in% names(target_data)) {
      as.matrix(target_data$W_outcome_true)
    } else {
      as.matrix(target_data$W_outcome)
    }
    if (length(alpha_true) != ncol(w_true) + 1L) {
      stop("target-only PS-basis diagnostic: true alpha has incompatible length.",
           call. = FALSE)
    }
    eta <- as.numeric(cbind(1, w_true) %*% alpha_true)
    oracle_m_pred <- if (identical(glm_spec$glmnet_family, "binomial")) {
      logistic(eta)
    } else {
      eta
    }
    oracle_m_pred <- clip_outcome_pred(oracle_m_pred, family)
  }

  oracle_prop_scores <- NULL
  if (identical(ps_basis, "oracle")) {
    if (!"p_treat_true" %in% names(target_data)) {
      stop("target-only PS-basis diagnostic: oracle PS requires p_treat_true.",
           call. = FALSE)
    }
    oracle_prop_scores <- clip_propensity(as.numeric(target_data$p_treat_true))
    if (length(oracle_prop_scores) != n) {
      stop("target-only PS-basis diagnostic: p_treat_true has incompatible length.",
           call. = FALSE)
    }
  }

  treated_idx <- which(a == 1L)
  control_idx <- which(a == 0L)
  if (length(treated_idx) < n_folds || length(control_idx) < n_folds) {
    stop(sprintf("target-only PS-basis diagnostic: too few observations for %d folds.",
                 n_folds), call. = FALSE)
  }

  fold_ids <- integer(n)
  fold_ids[treated_idx] <- sample(rep(seq_len(n_folds), length.out = length(treated_idx)))
  fold_ids[control_idx] <- sample(rep(seq_len(n_folds), length.out = length(control_idx)))

  prop_scores <- numeric(n)
  m_pred <- numeric(n)

  for (k in seq_len(n_folds)) {
    val_idx <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)

    if (identical(ps_basis, "oracle")) {
      prop_scores[val_idx] <- oracle_prop_scores[val_idx]
    } else {
      prop_scores[val_idx] <- fit_glmnet_cv(
        x_train = x_ps[train_idx, , drop = FALSE],
        y_train = as.numeric(a[train_idx]),
        x_predict = x_ps[val_idx, , drop = FALSE],
        family = "binomial",
        clip_fn = clip_propensity,
        caller_name = "target_only_ps_basis_crossfit",
        model_name = paste0("PS_", ps_basis),
        fold_id = k
      )
    }

    if (identical(or_basis, "oracle")) {
      m_pred[val_idx] <- oracle_m_pred[val_idx]
    } else {
      train_arm_idx <- train_idx[a[train_idx] == A_val]
      if (length(train_arm_idx) < MIN_TREATED_FOR_MODEL) {
        stop(sprintf("target-only PS-basis diagnostic: fold %d has only %d A_val=%d training observations.",
                     k, length(train_arm_idx), A_val), call. = FALSE)
      }

      m_pred[val_idx] <- fit_glmnet_cv(
        x_train = x_or[train_arm_idx, , drop = FALSE],
        y_train = y[train_arm_idx],
        x_predict = x_or[val_idx, , drop = FALSE],
        family = glm_spec$glmnet_family,
        clip_fn = function(pred) clip_outcome_pred(pred, family),
        caller_name = "target_only_ps_basis_crossfit",
        model_name = paste0("OR_", or_basis, "_ps_", ps_basis),
        fold_id = k
      )
    }
  }

  p_a <- if (A_val == 1L) prop_scores else (1 - prop_scores)
  phi <- calculate_aipw_pseudo_outcome(as.numeric(y), a, m_pred, p_a, A_val = A_val)
  if (any(!is.finite(phi))) {
    stop("target-only PS-basis diagnostic: non-finite AIPW pseudo-outcome.",
         call. = FALSE)
  }

  estimate <- mean(phi)
  influence <- phi - estimate
  variance <- mean(influence^2) / max(1L, n)

  list(
    estimate = as.numeric(estimate),
    se = as.numeric(sqrt(variance)),
    variance = as.numeric(variance),
    prop_scores_range = range(prop_scores),
    m_pred_range = range(m_pred),
    n = n,
    ps_basis = ps_basis,
    or_basis = or_basis
  )
}

.c2_ps_basis_row <- function(rep_idx, seed, method, res, truth, n_total, K, p,
                             n_folds) {
  estimate <- as.numeric(res$estimate)
  se <- as.numeric(res$se)
  ci_lower <- estimate - Z_ALPHA_05 * se
  ci_upper <- estimate + Z_ALPHA_05 * se
  data.frame(
    rep = rep_idx,
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
    K = K,
    p = p,
    n_folds = n_folds,
    ps_min = min(res$prop_scores_range),
    ps_max = max(res$prop_scores_range),
    m_min = min(res$m_pred_range),
    m_max = max(res$m_pred_range),
    stringsAsFactors = FALSE
  )
}

.c2_ps_basis_summary <- function(runs) {
  pieces <- split(runs, runs$method)
  do.call(rbind, lapply(names(pieces), function(method) {
    x <- pieces[[method]]
    bias <- x$bias
    se <- x$se
    data.frame(
      method = method,
      n_success = nrow(x),
      bias_mean = mean(bias),
      bias_sd = if (nrow(x) > 1L) stats::sd(bias) else NA_real_,
      rmse = sqrt(mean(bias^2)),
      se_mean = mean(se),
      coverage = mean(x$covered),
      bias_over_se = mean(bias) / mean(se),
      se_over_emp_sd = if (nrow(x) > 1L && stats::sd(bias) > 0) {
        mean(se) / stats::sd(bias)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
}

test_that("c2-target-ps-basis diagnoses target-only propensity basis", {
  skip_if_not(.c2_ps_basis_enabled(),
              message = "set ROCE_RUN_C2_TARGET_PS_BASIS=1 or run with --filter c2-target-ps-basis")

  n_reps <- .c2_ps_basis_int_env("ROCE_C2_TARGET_PS_BASIS_REPS", 5L)
  n_total <- .c2_ps_basis_int_env("ROCE_C2_TARGET_PS_BASIS_N", 5000L)
  K <- .c2_ps_basis_int_env("ROCE_C2_TARGET_PS_BASIS_K", 5L)
  p <- .c2_ps_basis_int_env("ROCE_C2_TARGET_PS_BASIS_P", 50L)
  n_folds <- .c2_ps_basis_int_env("ROCE_C2_TARGET_PS_BASIS_FOLDS", 10L)
  base_seed <- .c2_ps_basis_int_env("ROCE_C2_TARGET_PS_BASIS_SEED", 250000L)
  shift_strength <- .c2_ps_basis_num_env("ROCE_C2_TARGET_PS_BASIS_SHIFT", 0.5)

  out_dir <- file.path("c2_target_ps_basis_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  slurm_id <- Sys.getenv("SLURM_JOB_ID", "local")
  run_id <- paste(
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    paste0("job", gsub("[^A-Za-z0-9_.-]", "_", slurm_id)),
    paste0("n", n_total),
    paste0("K", K),
    paste0("p", p),
    sep = "_"
  )
  progress_path <- file.path(out_dir, paste0(run_id, "_progress.log"))
  run_path <- file.path(out_dir, paste0(run_id, "_runs.csv"))
  summary_path <- file.path(out_dir, paste0(run_id, "_summary.csv"))

  .c2_ps_basis_log(
    progress_path, "config",
    sprintf("reps=%d n=%d K=%d p=%d C2 folds=%d base_seed=%d",
            n_reps, n_total, K, p, n_folds, base_seed)
  )

  all_rows <- list()
  for (rep_idx in seq_len(n_reps)) {
    seed <- base_seed + rep_idx - 1L
    set.seed(seed)
    data <- generate_simulation_data(
      n_total = n_total,
      K = K,
      p = p,
      config = "C2",
      estimand_type = "superpopulation",
      dgp_type = "roce",
      outcome_type = "binary",
      site_allocation = "model",
      transform_type = "mild",
      shift_strength = shift_strength
    )
    target_data <- split_data_by_site(data)$t
    target_idx <- which(data$R == "t")
    target_data$p_treat_true <- data$p_treat_true[target_idx]
    truth <- as.numeric(data$mu1_true)

    current <- .target_only_ps_basis_crossfit(
      target_data, ps_basis = "W", or_basis = "W", n_folds = n_folds,
      family = "binomial", A_val = 1L,
      alpha1_true = data$alpha1_true, alpha0_true = data$alpha0_true
    )
    diagnostic <- .target_only_ps_basis_crossfit(
      target_data, ps_basis = "Z", or_basis = "W", n_folds = n_folds,
      family = "binomial", A_val = 1L,
      alpha1_true = data$alpha1_true, alpha0_true = data$alpha0_true
    )
    oracle_ps <- .target_only_ps_basis_crossfit(
      target_data, ps_basis = "oracle", or_basis = "W", n_folds = n_folds,
      family = "binomial", A_val = 1L,
      alpha1_true = data$alpha1_true, alpha0_true = data$alpha0_true
    )
    oracle_or <- .target_only_ps_basis_crossfit(
      target_data, ps_basis = "W", or_basis = "oracle", n_folds = n_folds,
      family = "binomial", A_val = 1L,
      alpha1_true = data$alpha1_true, alpha0_true = data$alpha0_true
    )
    oracle_both <- .target_only_ps_basis_crossfit(
      target_data, ps_basis = "oracle", or_basis = "oracle", n_folds = n_folds,
      family = "binomial", A_val = 1L,
      alpha1_true = data$alpha1_true, alpha0_true = data$alpha0_true
    )

    rows <- rbind(
      .c2_ps_basis_row(rep_idx, seed, "target_ps_W_current", current,
                       truth, n_total, K, p, n_folds),
      .c2_ps_basis_row(rep_idx, seed, "target_ps_Z_diagnostic", diagnostic,
                       truth, n_total, K, p, n_folds),
      .c2_ps_basis_row(rep_idx, seed, "target_ps_oracle", oracle_ps,
                       truth, n_total, K, p, n_folds),
      .c2_ps_basis_row(rep_idx, seed, "target_oracle_or_ps_W", oracle_or,
                       truth, n_total, K, p, n_folds),
      .c2_ps_basis_row(rep_idx, seed, "target_oracle_both", oracle_both,
                       truth, n_total, K, p, n_folds)
    )
    .c2_ps_basis_append_csv(rows, run_path)
    all_rows[[rep_idx]] <- rows

    .c2_ps_basis_log(
      progress_path, "rep_done",
      sprintf("rep=%d seed=%d W_bias=%.6f Z_bias=%.6f oracle_ps_bias=%.6f oracle_or_bias=%.6f oracle_both_bias=%.6f",
              rep_idx, seed, rows$bias[1], rows$bias[2], rows$bias[3],
              rows$bias[4], rows$bias[5])
    )
  }

  runs <- do.call(rbind, all_rows)
  summary <- .c2_ps_basis_summary(runs)
  .c2_ps_basis_append_csv(summary, summary_path)

  .c2_ps_basis_log(progress_path, "output", paste0("runs=", run_path))
  .c2_ps_basis_log(progress_path, "output", paste0("summary=", summary_path))
  print(summary)

  expect_true(file.exists(run_path))
  expect_true(file.exists(summary_path))
  expect_true(all(is.finite(runs$estimate)))
})
