# estimators_target.R - Target-only AIPW estimator variants
#
# Estimators that use only target site data (no source sites).
# Includes cross-fitted and complement-fold variants for different
# inference settings.
#
# Contents:
#   1. estimate_target_only (main entry point)
#   2. estimate_target_only_crossfit (standalone cross-fitting)
#   3. estimate_target_only_from_complement (for use inside FACE-HD outer loop)

#' Target-only estimator using doubly robust AIPW
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances
#'   for variance-valid inference (use_rcal is ignored in this mode).
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
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
#'   in main.R). When called inside the outer cross-fitting loop of FACE-HD
#'   algorithms (\code{run_crossfit(..., communication_mode = "one_round"/"two_round")}), use
#'   \code{\link{estimate_target_only_from_complement}} instead to avoid
#'   nested cross-fitting that inflates \code{V_ot}.
#' 
#' @param target_data Target site data (list with W_outcome, Z_site, A, Y, n)
#' @param n_folds Number of cross-fitting folds (default 3)
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param A_val Treatment value to estimate (default 1)
#' @return List with estimate (= mu_hat_ot), V_ot, variance,
#'   varphi_ot (CENTERED influence function: hat{varphi}_{ot,i} - hat{mu}^1_{ot}),
#'   and other components
#' @export
estimate_target_only_crossfit <- function(target_data, n_folds = 3, 
                                          family = "binomial", A_val = 1) {
  
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
      fold_id = k
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
      fold_id = k
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
    method = "crossfit"
  ))
}

#' Target-only estimator using complement-fold training (no nested cross-fitting)
#'
#' Trains propensity and outcome models on the COMPLEMENT of fold k1
#' (all other folds combined) and evaluates AIPW pseudo-outcomes on fold k1.
#' This avoids the nested cross-fitting that inflates V_ot when called on
#' small per-fold data.
#'
#' Intended to be used **inside** the main cross-fitting loop of FACE-HD
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
#' @return List with estimate, V_ot, variance, varphi_ot (same interface as
#'   \code{\link{estimate_target_only_crossfit}})
#' @export
estimate_target_only_from_complement <- function(target_folds, k1, n_folds,
                                                  k2 = NULL,
                                                  family = "binomial", A_val = 1L,
                                                  propensity_cache = NULL) {
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
  x_or_eval <- as.matrix(eval_data$W_outcome)
  x_ps_eval <- if (!is.null(eval_data$Z_site)) {
    as.matrix(eval_data$Z_site)
  } else {
    x_or_eval
  }
  y_eval <- eval_data$Y
  tr_eval <- eval_data$A
  n_eval <- eval_data$n

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
    ";ps_basis=Z_site"
  )

  # --- Propensity score model: train on complement, predict on eval fold ---
  if (!is.null(propensity_cache) && exists(cache_key, envir = propensity_cache, inherits = FALSE)) {
    prop_scores_eval <- get(cache_key, envir = propensity_cache, inherits = FALSE)
  } else {
    if (length(unique(tr_train)) < 2) {
      stop(sprintf(
        "estimate_target_only_from_complement: training data after excluding fold(s) %s has only one treatment class; propensity model is not identifiable.",
        paste(exclude_folds, collapse = ",")
      ), call. = FALSE)
    }
    prop_scores_eval <- fit_glmnet_cv(
      x_train = x_ps_train, y_train = as.numeric(tr_train), x_predict = x_ps_eval,
      family = "binomial",
      clip_fn = clip_propensity,
      caller_name = "estimate_complement_fold_aipw", model_name = "PS"
    )
    if (!is.null(propensity_cache)) {
      assign(cache_key, prop_scores_eval, envir = propensity_cache)
    }
  }
  prop_scores_eval <- clip_propensity(prop_scores_eval)

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
  m_pred_eval <- fit_glmnet_cv(
    x_train = X_treated_train, y_train = y_treated_train, x_predict = x_or_eval,
    family = glm_spec$glmnet_family,
    clip_fn = function(pred) clip_outcome_pred(pred, family),
    caller_name = "estimate_complement_fold_aipw", model_name = "OR"
  )
  m_pred_eval <- clip_outcome_pred(m_pred_eval, family)

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

  return(list(
    estimate = as.numeric(mu_hat_ot),
    V_ot = as.numeric(V_ot),
    variance = as.numeric(variance),
    varphi_ot = as.numeric(varphi_ot),
    prop_scores = prop_scores_eval,
    m_pred = m_pred_eval,
    method = "complement"
  ))
}
