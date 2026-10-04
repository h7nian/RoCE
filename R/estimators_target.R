# estimators_target.R - Target-only AIPW estimator variants
#
# Estimators that use only target site data (no source sites).
# Includes cross-fitted and complement-fold variants for different
# inference settings.
#
# Contents:
#   1. estimate_target_only (main entry point)
#   2. estimate_target_only_crossfit (standalone cross-fitting)
#   3. estimate_target_only_from_complement (for use inside RoCE outer loop)

.target_only_cv_seed <- function(k1, k2 = NULL, model, A_val = 0L) {
  model_index <- match(model, c("propensity", "outcome"))
  if (is.na(model_index)) {
    stop(".target_only_cv_seed: model must be propensity or outcome.",
         call. = FALSE)
  }
  inner_fold <- if (is.null(k2)) 0L else as.integer(k2)
  as.integer(
    800011L + 10000L * model_index + 1000L * as.integer(A_val) +
      100L * as.integer(k1) + inner_fold
  )
}

#' Target-only estimator using doubly robust AIPW
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances
#'   nuisance fitting. Combining it with \code{use_rcal = TRUE} is unsupported
#'   and raises an error.
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param A_val Treatment arm, either 0 or 1.
#' @return estimate with variance
#' @export
estimate_target_only <- function(data_split, family = "binomial", 
                                 use_rcal = FALSE, use_crossfit = TRUE,
                                 n_folds = NULL, A_val = 1L) {
  
  target_data <- data_split[["t"]]
  res <- fit_site_aipw(target_data, family = family, use_rcal = use_rcal,
                       use_crossfit = use_crossfit, n_folds = n_folds, A_val = A_val)
  res$method <- "target_only"
  return(res)
}

#' Cross-fitted Target-Only AIPW Estimator (standalone)
#' 
#' This function implements cross-fitting for the target-only estimator to avoid
#' overfitting bias in the influence function and variance estimation.
#' Each observation's IF is computed using nuisance parameters trained on
#' out-of-fold data, following standard cross-fitting methodology.
#'
#' @note This function performs its OWN internal K-fold cross-fitting and is
#'   intended for **standalone** use (e.g., comparison methods, ATE estimation
#'   in main.R). When called inside the outer cross-fitting loop of RoCE
#'   algorithms (\code{run_crossfit(..., communication_mode = "one_round"/"two_round")}), use
#'   \code{\link{estimate_target_only_from_complement}} instead to avoid
#'   nested cross-fitting that inflates \code{V_ot}.
#' 
#' @param target_data Target site data (list with W_outcome, Z_site, A, Y, n)
#' @param n_folds Number of cross-fitting folds (default 3)
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param A_val Treatment value to estimate (default 1)
#' @param nuisance_lambda_rule Cross-validation rule used for both target
#'   propensity and outcome nuisance fits.
#' @return List with estimate (= mu_hat_ot), V_ot, variance,
#'   varphi_ot (CENTERED influence function: hat{varphi}_{ot,i} - hat{mu}^1_{ot}),
#'   and other components
#' @export
estimate_target_only_crossfit <- function(target_data, n_folds = 3, 
                                          family = "binomial", A_val = 1,
                                          nuisance_lambda_rule = c("min", "1se")) {
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, "estimate_target_only_crossfit",
    arg = "nuisance_lambda_rule"
  )
  
  # resolve_glm_family provides canonical family/link mapping.
  glm_spec <- resolve_glm_family(family)
  link <- glm_spec$link
  
  y <- target_data$Y
  tr <- target_data$A
  x_or <- as.matrix(target_data$W_outcome)
  x_ps <- if (!is.null(target_data$Z_site)) {
    as.matrix(target_data$Z_site)
  } else {
    x_or
  }
  n <- target_data$n
  
  if (sum(tr == A_val) == 0) {
    stop(sprintf("estimate_target_only_crossfit: no observations with A_val=%d in target data.",
                 A_val), call. = FALSE)
  }
  
  # Assign folds with stratification by treatment to ensure each fold has
  # both treated and control units (mirrors partition_into_folds logic)
  fold_ids <- integer(n)
  treated_idx <- which(tr == 1)
  control_idx <- which(tr == 0)
  
  if (length(treated_idx) < n_folds || length(control_idx) < n_folds) {
    stop(sprintf(
      "estimate_target_only_crossfit: too few treated or control units for %d stratified folds (treated=%d, control=%d). Reduce n_folds or provide data with at least one observation per treatment arm in every fold.",
      n_folds, length(treated_idx), length(control_idx)
    ), call. = FALSE)
  }
  # Stratified assignment: shuffle each arm separately
  fold_ids[treated_idx] <- sample(rep(1:n_folds, length.out = length(treated_idx)))
  fold_ids[control_idx] <- sample(rep(1:n_folds, length.out = length(control_idx)))
  
  # Initialize output vectors
  prop_scores_oof <- numeric(n)
  m_pred_oof <- numeric(n)
  
  # Cross-fitting loop
  for (k in 1:n_folds) {
    val_idx <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)

    training_preprocessor <- attr(target_data, "training_preprocessor", exact = TRUE)
    if (!is.null(training_preprocessor)) {
      if (!is.function(training_preprocessor)) stop("training_preprocessor must be a function.", call. = FALSE)
      transformed <- training_preprocessor(train_idx)
      if (!identical(dim(transformed$W_outcome), dim(x_or)) ||
          !identical(dim(transformed$Z_site), dim(x_ps)) ||
          any(!is.finite(transformed$W_outcome)) || any(!is.finite(transformed$Z_site))) {
        stop("Training preprocessing must return finite aligned feature matrices.", call. = FALSE)
      }
      x_or <- transformed$W_outcome
      x_ps <- transformed$Z_site
    }
    
    if (length(val_idx) == 0 || length(train_idx) == 0) {
      stop(sprintf(
        "estimate_target_only_crossfit: fold %d has an empty validation or training split (n_val=%d, n_train=%d).",
        k, length(val_idx), length(train_idx)
      ), call. = FALSE)
    }
    
    X_or_train <- x_or[train_idx, , drop = FALSE]
    X_or_val <- x_or[val_idx, , drop = FALSE]
    X_ps_train <- x_ps[train_idx, , drop = FALSE]
    X_ps_val <- x_ps[val_idx, , drop = FALSE]
    y_train <- y[train_idx]
    tr_train <- tr[train_idx]
    
    # Train propensity score model on training fold
    if (length(unique(tr_train)) < 2) {
      stop(sprintf(
        "estimate_target_only_crossfit: fold %d training data has only one treatment class; propensity model is not identifiable.",
        k
      ), call. = FALSE)
    }
    ps_fit <- fit_glmnet_cv(
      x_train = X_ps_train, y_train = as.numeric(tr_train), x_predict = X_ps_val,
      family = "binomial",
      caller_name = "estimate_target_only_crossfit", model_name = "PS",
      fold_id = k,
      lambda_rule = nuisance_lambda_rule
    )
    
    prop_scores_oof[val_idx] <- clip_propensity(ps_fit)
    
    # Train outcome model on training fold (only on the requested treatment arm)
    treated_train_idx <- which(tr_train == A_val)
    if (length(treated_train_idx) < MIN_TREATED_FOR_MODEL) {
      stop(sprintf(
        "estimate_target_only_crossfit: fold %d has only %d training observations with A_val=%d; need at least %d for outcome model.",
        k, length(treated_train_idx), A_val, MIN_TREATED_FOR_MODEL
      ), call. = FALSE)
    }
    X_treated_train <- X_or_train[treated_train_idx, , drop = FALSE]
    y_treated_train <- y_train[treated_train_idx]

    or_fit <- fit_glmnet_cv(
      x_train = X_treated_train, y_train = y_treated_train, x_predict = X_or_val,
      family = glm_spec$glmnet_family,
      clip_fn = function(pred) clip_outcome_pred(pred, family),
      caller_name = "estimate_target_only_crossfit", model_name = "OR",
      fold_id = k,
      on_degenerate_response = "constant",
      lambda_rule = nuisance_lambda_rule
    )

    m_pred_oof[val_idx] <- or_fit
  }
  
  # Compute AIPW estimator with cross-fitted nuisance
  p_a_oof <- if (A_val == 1L) prop_scores_oof else (1 - prop_scores_oof)
  phi_i <- calculate_aipw_pseudo_outcome(as.numeric(y), tr, m_pred_oof, p_a_oof, A_val = A_val)
  
  # Handle any NA or Inf values
  bad_idx <- is.na(phi_i) | is.infinite(phi_i)
  if (any(bad_idx)) {
    stop(sprintf("estimate_target_only_crossfit: AIPW pseudo-outcome produced %d non-finite value(s).",
                 sum(bad_idx)), call. = FALSE)
  }
  
  mu_hat_ot <- mean(phi_i)
  varphi_ot <- phi_i - mu_hat_ot
  V_ot <- mean(varphi_ot^2)
  variance <- V_ot / max(1, n)
  
  return(list(
    estimate = as.numeric(mu_hat_ot),
    V_ot = as.numeric(V_ot),
    variance = as.numeric(variance),
    varphi_ot = as.numeric(varphi_ot),
    prop_scores = prop_scores_oof,
    m_pred = m_pred_oof,
    method = "crossfit",
    nuisance_lambda_rule = nuisance_lambda_rule
  ))
}

#' Target-only estimator using complement-fold training (no nested cross-fitting)
#'
#' Trains propensity and outcome models on the COMPLEMENT of fold k1
#' (all other folds combined) and evaluates AIPW pseudo-outcomes on fold k1.
#' This avoids the nested cross-fitting that inflates V_ot when called on
#' small per-fold data.
#'
#' Intended to be used **inside** the main cross-fitting loop of RoCE
#' algorithms (\code{run_crossfit(..., communication_mode = "one_round"/"two_round")}) where the
#' outer loop already provides cross-fitting structure. For standalone use
#' (no outer loop), use \code{\link{estimate_target_only_crossfit}} instead.
#'
#' @param target_folds Fold list from partition_into_folds (target site)
#' @param k1 Main estimation fold index (outer fold to exclude from training)
#' @param n_folds Total number of folds
#' @param k2 Optional inner fold index. When non-NULL, the function operates in
#'   "inner-fold mode": training excludes both k1 and k2, and evaluation is
#'   performed on fold k2 (instead of k1). This is used for Version A Step 2
#'   (inner M-fold variance component estimation). Default NULL = outer-fold mode.
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param A_val Treatment value to estimate (default 1)
#' @param propensity_cache Optional environment used to cache complement-fold
#'   propensity predictions keyed by excluded/evaluation folds.
#' @param return_training_scores If TRUE, retain predictions and AIPW scores
#'   on the same training observations for two-layer aggregation-weight learning.
#' @param nuisance_lambda_rule Cross-validation rule used for both target
#'   propensity and outcome nuisance fits.
#' @return List with estimate, V_ot, variance, varphi_ot (same interface as
#'   \code{\link{estimate_target_only_crossfit}})
#' @export
estimate_target_only_from_complement <- function(target_folds, k1, n_folds,
                                                  k2 = NULL,
                                                  family = "binomial", A_val = 1L,
                                                  propensity_cache = NULL,
                                                  nuisance_lambda_rule = c("min", "1se"),
                                                  return_training_scores = FALSE) {
  return_training_scores <- .validate_logical_control(return_training_scores, "return_training_scores")
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, "estimate_target_only_from_complement",
    arg = "nuisance_lambda_rule"
  )
  glm_spec <- resolve_glm_family(family)
  link <- glm_spec$link

  # Training data = all folds except k1 (and k2 if inner-fold mode)
  exclude_folds <- if (is.null(k2)) k1 else c(k1, k2)
  training_folds <- setdiff(1:n_folds, exclude_folds)
  train_data <- combine_folds(target_folds, training_folds)

  # Evaluation data = k2 if inner-fold mode, otherwise k1
  eval_fold <- if (is.null(k2)) k1 else k2
  eval_data <- materialize_fold(target_folds, eval_fold)

  x_or_train <- as.matrix(train_data$W_outcome)
  x_ps_train <- if (!is.null(train_data$Z_site)) {
    as.matrix(train_data$Z_site)
  } else {
    x_or_train
  }
  y_train <- train_data$Y
  tr_train <- train_data$A
  train_cv_group_id <- train_data$cv_group_id %||% NULL
  x_or_eval <- as.matrix(eval_data$W_outcome)
  x_ps_eval <- if (!is.null(eval_data$Z_site)) {
    as.matrix(eval_data$Z_site)
  } else {
    x_or_eval
  }
  y_eval <- eval_data$Y
  tr_eval <- eval_data$A
  n_eval <- eval_data$n
  x_or_predict <- if (return_training_scores) rbind(x_or_eval, x_or_train) else x_or_eval
  x_ps_predict <- if (return_training_scores) rbind(x_ps_eval, x_ps_train) else x_ps_eval

  treated_train_idx <- which(tr_train == A_val)
  if (length(treated_train_idx) == 0) {
    stop(sprintf("estimate_target_only_from_complement: no training observations with A_val=%d after excluding fold(s) %s.",
                 A_val, paste(exclude_folds, collapse = ",")), call. = FALSE)
  }

  cache_key <- paste0(
    "k1=", k1,
    ";k2=", if (is.null(k2)) "NA" else as.character(k2),
    ";eval=", eval_fold,
    ";exclude=", paste(sort(exclude_folds), collapse = ","),
    ";ps_basis=Z_site",
    ";lambda_rule=", nuisance_lambda_rule
  )
  if (!is.null(train_cv_group_id)) {
    validated_cache_groups <- .validate_nuisance_cv_group_id(
      train_cv_group_id, nrow(x_ps_train),
      "estimate_target_only_from_complement propensity cache"
    )
    # Exact IDs are included only for the opt-in grouped path. The legacy
    # NULL cache key above remains byte-for-byte unchanged.
    cache_key <- paste0(
      cache_key, ";cv_group_id=",
      paste(validated_cache_groups, collapse = ",")
    )
  }

  # Cache keys locate candidates; exact fitting/prediction inputs decide reuse.
  propensity_inputs <- list(
    x_train = x_ps_train, treatment = tr_train, x_predict = x_ps_predict,
    cv_group_id = train_cv_group_id, lambda_rule = nuisance_lambda_rule,
    training_policy = NUISANCE_TRAINING_POLICY
  )
  cached_propensity <- if (is.null(propensity_cache)) NULL else propensity_cache[[cache_key]]
  if (is.list(cached_propensity) && identical(cached_propensity$inputs, propensity_inputs)) {
    propensity_predictions <- cached_propensity$predictions
  } else {
    if (length(unique(tr_train)) < 2) {
      stop(sprintf(
        "estimate_target_only_from_complement: training data after excluding fold(s) %s has only one treatment class; propensity model is not identifiable.",
        paste(exclude_folds, collapse = ",")
      ), call. = FALSE)
    }
    propensity_predictions <- with_seed(
      .target_only_cv_seed(k1, k2, "propensity"),
      fit_glmnet_cv(
        x_train = x_ps_train, y_train = as.numeric(tr_train),
        x_predict = x_ps_predict,
        family = "binomial",
        clip_fn = clip_propensity,
        caller_name = "estimate_complement_fold_aipw", model_name = "PS",
        lambda_rule = nuisance_lambda_rule,
        cv_group_id = train_cv_group_id
      )
    )
    if (!is.null(propensity_cache)) {
      propensity_cache[[cache_key]] <- list(inputs = propensity_inputs, predictions = propensity_predictions)
    }
  }
  propensity_predictions <- clip_propensity(propensity_predictions)

  # --- Outcome model: train on treated in complement, predict on fold k1 ---
  X_treated_train <- x_or_train[treated_train_idx, , drop = FALSE]
  y_treated_train <- y_train[treated_train_idx]

  if (length(treated_train_idx) < MIN_TREATED_FOR_MODEL) {
    stop(sprintf(
      "estimate_target_only_from_complement: only %d training observations with A_val=%d after excluding fold(s) %s; need at least %d.",
      length(treated_train_idx), A_val, paste(exclude_folds, collapse = ","),
      MIN_TREATED_FOR_MODEL
    ), call. = FALSE)
  }
  outcome_predictions <- with_seed(
    .target_only_cv_seed(k1, k2, "outcome", A_val),
    fit_glmnet_cv(
      x_train = X_treated_train, y_train = y_treated_train,
      x_predict = x_or_predict,
      family = glm_spec$glmnet_family,
      clip_fn = function(pred) clip_outcome_pred(pred, family),
      caller_name = "estimate_complement_fold_aipw", model_name = "OR",
      on_degenerate_response = "constant",
      lambda_rule = nuisance_lambda_rule,
      cv_group_id = if (is.null(train_cv_group_id)) {
        NULL
      } else {
        train_cv_group_id[treated_train_idx]
      }
    )
  )
  outcome_degenerate <- as.integer(
    attr(outcome_predictions, "outcome_degenerate") %||% 0L
  )
  outcome_predictions <- clip_outcome_pred(outcome_predictions, family)

  training_scores <- NULL
  if (return_training_scores) {
    rows <- n_eval + seq_len(train_data$n)
    training_propensity <- propensity_predictions[rows]
    training_outcome <- outcome_predictions[rows]
    training_scores <- list(
      training_folds = as.integer(training_folds),
      original_idx = unlist(lapply(target_folds[training_folds], `[[`, "original_idx"), use.names = FALSE),
      prop_scores = training_propensity, m_pred = training_outcome,
      phi = calculate_aipw_pseudo_outcome(as.numeric(y_train), tr_train, training_outcome,
        if (A_val == 1L) training_propensity else 1 - training_propensity, A_val = A_val))
  }

  prop_scores_eval <- if (return_training_scores) propensity_predictions[seq_len(n_eval)] else propensity_predictions
  m_pred_eval <- if (return_training_scores) outcome_predictions[seq_len(n_eval)] else outcome_predictions

  # --- AIPW pseudo-outcomes on fold k1 ---
  p_a <- if (A_val == 1L) prop_scores_eval else (1 - prop_scores_eval)
  phi_i <- calculate_aipw_pseudo_outcome(as.numeric(y_eval), tr_eval, m_pred_eval, p_a, A_val = A_val)

  bad_idx <- is.na(phi_i) | is.infinite(phi_i)
  if (any(bad_idx)) {
    stop(sprintf("estimate_target_only_from_complement: AIPW pseudo-outcome produced %d non-finite value(s).",
                 sum(bad_idx)), call. = FALSE)
  }

  mu_hat_ot <- mean(phi_i)
  varphi_ot <- phi_i - mu_hat_ot
  V_ot <- mean(varphi_ot^2)
  variance <- V_ot / max(1, n_eval)

  result <- list(
    estimate = as.numeric(mu_hat_ot),
    V_ot = as.numeric(V_ot),
    variance = as.numeric(variance),
    varphi_ot = as.numeric(varphi_ot),
    prop_scores = prop_scores_eval,
    m_pred = m_pred_eval,
    outcome_degenerate = outcome_degenerate,
    method = "complement",
    nuisance_lambda_rule = nuisance_lambda_rule
  )
  if (return_training_scores) result$training_scores <- training_scores
  result
}

# Fixed ordinary-AIPW reference for paired method comparisons. Selecting a
# new RoCE target anchor must not change the benchmark used to assess it.
.target_tate_reference <- function(target_folds, family, nuisance_lambda_rule) {
  n_folds <- length(target_folds)
  n <- sum(vapply(target_folds, `[[`, numeric(1L), "n"))
  propensity_cache <- new.env(parent = emptyenv())
  score <- matrix(0, n, 2L)
  counts <- integer(n)
  for (fold in seq_len(n_folds)) {
    rows <- target_folds[[fold]]$original_idx
    counts[rows] <- counts[rows] + 1L
    for (arm in 0:1) {
      fit <- .get_target_only_fold_fit(
        target_folds, fold, n_folds, family, arm, propensity_cache,
        nuisance_lambda_rule = nuisance_lambda_rule)
      score[rows, arm + 1L] <- fit$varphi_ot + fit$estimate
    }
  }
  if (any(counts != 1L)) stop("Target reference folds do not form a partition.", call. = FALSE)
  contrast <- score[, 2L] - score[, 1L]
  variance <- .multisite_pseudovalue_variance(contrast, as.integer(n))
  list(estimate = mean(contrast), se = sqrt(variance), variance = variance)
}
