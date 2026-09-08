# Fixed-nuisance bootstrap for data-adaptive TATE aggregation weights.

.weight_bootstrap_weighted_mean <- function(x, w, caller) {
  if (!is.numeric(x) || !is.numeric(w) || length(x) != length(w) ||
      length(x) == 0L || any(!is.finite(x)) || any(!is.finite(w)) ||
      any(w < 0) || sum(w) <= 0) {
    stop(caller, ": values and non-negative multipliers must be finite, ",
         "non-empty, equally sized, and have positive total weight.",
         call. = FALSE)
  }
  sum(w * x) / sum(w)
}

.weight_bootstrap_centered_moment <- function(x, w, caller) {
  center <- .weight_bootstrap_weighted_mean(x, w, caller)
  .weight_bootstrap_weighted_mean((x - center)^2, w, caller)
}

.weight_bootstrap_covariance <- function(x, y, w, caller) {
  center_x <- .weight_bootstrap_weighted_mean(x, w, caller)
  center_y <- .weight_bootstrap_weighted_mean(y, w, caller)
  .weight_bootstrap_weighted_mean(
    (x - center_x) * (y - center_y), w, caller
  )
}

.permute_weight_bootstrap_fold_info <- function(info, order_index) {
  vector_fields <- c(
    "V_t", "V_s", "C_ot", "fold_source_estimates", "mu_pred_ts",
    "delta_ts"
  )
  list_fields <- c(
    "zeta_components", "xi_components", "source_idx"
  )
  for (field in vector_fields) {
    if (!is.null(info[[field]])) info[[field]] <- info[[field]][order_index]
  }
  for (field in list_fields) {
    if (!is.null(info[[field]])) info[[field]] <- info[[field]][order_index]
  }
  if (!is.null(info$C_cross)) {
    info$C_cross <- info$C_cross[order_index, order_index, drop = FALSE]
  }
  info
}

.permute_weight_bootstrap_component <- function(component, order_index) {
  vector_fields <- c(
    "V_t", "V_s", "C_ot", "avg_source_est", "mu_pred_ts", "delta_ts",
    "n_s"
  )
  list_fields <- c("zeta_components", "xi_components", "source_idx")
  for (field in vector_fields) {
    if (!is.null(component[[field]])) {
      component[[field]] <- component[[field]][order_index]
    }
  }
  for (field in list_fields) {
    if (!is.null(component[[field]])) {
      component[[field]] <- component[[field]][order_index]
    }
  }
  if (!is.null(component$C_cross)) {
    component$C_cross <- component$C_cross[
      order_index, order_index, drop = FALSE
    ]
  }
  component
}

.validate_weight_bootstrap_result <- function(tate_result) {
  caller <- "estimate_tate_weight_bootstrap"
  fail <- function(...) stop(caller, ": ", ..., call. = FALSE)
  integer_scalar <- function(x, label, minimum = 1L) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < minimum || x != floor(x) || x > .Machine$integer.max) {
      fail(label, " must be one integer in [", minimum, ", ",
           .Machine$integer.max, "].")
    }
    as.integer(x)
  }
  integer_vector <- function(x, label, expected_length) {
    if (!is.numeric(x) || length(x) != expected_length || anyNA(x) ||
        any(!is.finite(x)) || any(x < 1) || any(x != floor(x)) ||
        any(x > .Machine$integer.max)) {
      fail(label, " must contain exactly ", expected_length,
           " positive integer sample sizes.")
    }
    as.integer(x)
  }
  indices <- function(x, size, label) {
    if (!is.numeric(x) || length(x) == 0L || anyNA(x) ||
        any(!is.finite(x)) || any(x != floor(x)) || any(x < 1) ||
        any(x > size) || anyDuplicated(x)) {
      fail(label, " must contain unique integer observation IDs in [1, ",
           size, "].")
    }
    as.integer(x)
  }
  finite_values <- function(x, expected_length, label) {
    if (!is.numeric(x) || length(x) != expected_length ||
        any(!is.finite(x))) {
      fail(label, " must contain exactly ", expected_length,
           " finite numeric value(s).")
    }
  }

  required <- c(
    "estimate", "se", "weights", "fold_weights", "fold_lambdas",
    "aggregation_lambda_selection", "aggregation_lambda_rule",
    "aggregation_screening_rule", "n_folds", "N_all", "intermediates"
  )
  missing <- required[!vapply(required, function(field) {
    !is.null(tate_result[[field]])
  }, logical(1L))]
  if (length(missing) > 0L) {
    fail("tate_result is missing field(s): ",
         paste(missing, collapse = ", "), ".")
  }
  intermediates <- tate_result$intermediates
  required_intermediates <- c("fold_info", "inner_fold_info", "sample_sizes")
  missing_intermediates <- required_intermediates[
    !vapply(required_intermediates, function(field) {
      !is.null(intermediates[[field]])
    }, logical(1L))
  ]
  if (length(missing_intermediates) > 0L) {
    fail("tate_result$intermediates is missing field(s): ",
         paste(missing_intermediates, collapse = ", "), ".")
  }

  source_names <- names(tate_result$weights)
  if (!is.numeric(tate_result$weights) || length(tate_result$weights) == 0L ||
      any(!is.finite(tate_result$weights)) || is.null(source_names) ||
      anyNA(source_names) || any(source_names == "") ||
      anyDuplicated(source_names)) {
    fail("tate_result$weights must be a finite vector with unique non-empty ",
         "source names.")
  }
  K <- length(source_names)
  n_folds <- integer_scalar(tate_result$n_folds, "n_folds", 2L)
  sample_sizes <- intermediates$sample_sizes
  n_t <- integer_scalar(sample_sizes$n_t, "sample_sizes$n_t")
  n_source <- integer_vector(
    sample_sizes$n_source, "sample_sizes$n_source", K
  )
  if (is.null(names(sample_sizes$n_source)) ||
      !identical(names(sample_sizes$n_source), source_names)) {
    fail("sample_sizes$n_source names and order must exactly match weights.")
  }
  names(n_source) <- source_names
  N_all <- integer_scalar(tate_result$N_all, "N_all")
  if (as.double(n_t) + sum(as.double(n_source)) != as.double(N_all)) {
    fail("N_all does not equal n_t + sum(n_source).")
  }
  if (length(intermediates$fold_info) != n_folds ||
      length(intermediates$inner_fold_info) != n_folds) {
    fail("outer and inner fold lists must both have length n_folds.")
  }
  fold_weights <- as.matrix(tate_result$fold_weights)
  if (!is.numeric(fold_weights) ||
      !identical(dim(fold_weights), c(n_folds, K)) ||
      any(!is.finite(fold_weights)) || is.null(colnames(fold_weights)) ||
      !identical(colnames(fold_weights), source_names)) {
    fail("fold_weights must be a finite n_folds-by-K matrix whose source ",
         "columns exactly match weights.")
  }
  finite_values(tate_result$fold_lambdas, n_folds, "fold_lambdas")
  if (!is.numeric(tate_result$aggregation_lambda_selection) ||
      length(tate_result$aggregation_lambda_selection) != 1L ||
      !is.finite(tate_result$aggregation_lambda_selection)) {
    fail("only a fixed finite numeric aggregation_lambda_selection is ",
         "supported; lambda_selection = 'cv' requires weighted inner CV.")
  }

  outer_target_ids <- vector("list", n_folds)
  outer_source_ids <- stats::setNames(
    lapply(seq_len(K), function(j) vector("list", n_folds)), source_names
  )
  # Validate outer folds first, so every inner component can be checked against
  # one authoritative global-ID sequence rather than merely against its length.
  for (k in seq_len(n_folds)) {
    info <- intermediates$fold_info[[k]]
    if (is.null(info$target_idx)) {
      fail("outer fold ", k, " is missing target_idx observation IDs.")
    }
    if (!is.list(info$source_idx) || length(info$source_idx) != K ||
        (!is.null(names(info$source_idx)) &&
         !identical(names(info$source_idx), source_names))) {
      fail("outer fold ", k,
           " source_idx order must match weights; when names are present they ",
           "must match exactly.")
    }
    outer_target_ids[[k]] <- indices(
      info$target_idx, n_t, paste0("outer fold ", k, " target_idx")
    )
    finite_values(info$varphi_ot, length(outer_target_ids[[k]]),
                  paste0("outer fold ", k, " varphi_ot"))
    finite_values(info$fold_target_estimate, 1L,
                  paste0("outer fold ", k, " fold_target_estimate"))
    for (field in c("fold_source_estimates", "mu_pred_ts", "delta_ts")) {
      finite_values(info[[field]], K, paste0("outer fold ", k, " ", field))
    }
    if (!is.list(info$zeta_components) ||
        length(info$zeta_components) != K ||
        !is.list(info$xi_components) || length(info$xi_components) != K) {
      fail("outer fold ", k, " must contain K zeta and xi components.")
    }
    for (j in seq_len(K)) {
      source <- source_names[[j]]
      outer_source_ids[[source]][[k]] <- indices(
        info$source_idx[[j]], n_source[[j]],
        paste0("outer fold ", k, " source_idx$", source)
      )
      finite_values(
        info$zeta_components[[j]], length(outer_target_ids[[k]]),
        paste0("outer fold ", k, " zeta_components$", source)
      )
      finite_values(
        info$xi_components[[j]], length(outer_source_ids[[source]][[k]]),
        paste0("outer fold ", k, " xi_components$", source)
      )
    }
  }
  if (!identical(sort(unlist(outer_target_ids, use.names = FALSE)),
                 seq_len(n_t))) {
    fail("outer target IDs must partition 1:n_t exactly once.")
  }
  for (j in seq_len(K)) {
    source <- source_names[[j]]
    if (!identical(
      sort(unlist(outer_source_ids[[source]], use.names = FALSE)),
      seq_len(n_source[[j]])
    )) {
      fail("outer source IDs for ", source,
           " must partition its observations exactly once.")
    }
  }

  for (k in seq_len(n_folds)) {
    components <- intermediates$inner_fold_info[[k]]$components
    if (is.null(components) || length(components) != n_folds - 1L) {
      fail("incomplete outer/inner fold metadata at fold ", k, ".")
    }
    inner_ids <- integer(length(components))
    for (component_index in seq_along(components)) {
      component <- components[[component_index]]
      if (is.null(component$outer_fold) || is.null(component$inner_fold) ||
          is.null(component$target_idx) || !is.list(component$source_idx) ||
          length(component$source_idx) != K ||
          is.null(names(component$source_idx)) ||
          !identical(names(component$source_idx), source_names)) {
        fail("incomplete explicit inner-fold indices at fold ", k,
             ". Refit with the current RoCE version.")
      }
      outer_id <- integer_scalar(
        component$outer_fold,
        paste0("inner component outer_fold at outer fold ", k)
      )
      inner_id <- integer_scalar(
        component$inner_fold,
        paste0("inner component inner_fold at outer fold ", k)
      )
      if (outer_id != k || inner_id == k || inner_id > n_folds) {
        fail("inner component fold labels are inconsistent at outer fold ", k,
             ".")
      }
      inner_ids[[component_index]] <- inner_id
      target_idx <- indices(
        component$target_idx, n_t,
        paste0("outer fold ", k, " inner fold ", inner_id, " target_idx")
      )
      if (!identical(target_idx, outer_target_ids[[inner_id]])) {
        fail("outer fold ", k, " inner fold ", inner_id,
             " target_idx must exactly match that fold's outer target IDs ",
             "in order.")
      }
      finite_values(component$varphi_ot, length(target_idx),
                    paste0("inner fold ", inner_id, " varphi_ot"))
      finite_values(component$avg_target_est, 1L,
                    paste0("inner fold ", inner_id, " avg_target_est"))
      finite_values(component$n_t, 1L,
                    paste0("inner fold ", inner_id, " n_t"))
      if (component$n_t != length(target_idx)) {
        fail("inner fold ", inner_id,
             " n_t must equal the target_idx length.")
      }
      for (field in c("avg_source_est", "mu_pred_ts", "delta_ts",
                      "V_t", "V_s", "C_ot", "n_s")) {
        finite_values(component[[field]], K,
                      paste0("inner fold ", inner_id, " ", field))
      }
      if (!is.matrix(component$C_cross) ||
          !identical(dim(component$C_cross), c(K, K)) ||
          any(!is.finite(component$C_cross)) ||
          !is.list(component$zeta_components) ||
          length(component$zeta_components) != K ||
          !is.list(component$xi_components) ||
          length(component$xi_components) != K) {
        fail("inner fold ", inner_id,
             " has inconsistent covariance or influence-component dimensions.")
      }
      for (j in seq_len(K)) {
        source <- source_names[[j]]
        source_idx <- indices(
          component$source_idx[[j]], n_source[[j]],
          paste0("outer fold ", k, " inner fold ", inner_id,
                 " source_idx$", source)
        )
        if (!identical(source_idx, outer_source_ids[[source]][[inner_id]])) {
          fail("outer fold ", k, " inner fold ", inner_id,
               " source_idx$", source, " must exactly match that fold's ",
               "outer source IDs in order.")
        }
        if (component$n_s[[j]] != length(source_idx)) {
          fail("inner fold ", inner_id, " n_s$", source,
               " must equal its source_idx length.")
        }
        finite_values(component$zeta_components[[j]], length(target_idx),
                      paste0("inner fold ", inner_id, " zeta$", source))
        finite_values(component$xi_components[[j]], length(source_idx),
                      paste0("inner fold ", inner_id, " xi$", source))
      }
    }
    if (!identical(sort(inner_ids), setdiff(seq_len(n_folds), k))) {
      fail("inner_fold labels at outer fold ", k,
           " must cover every fold except the outer fold exactly once.")
    }
  }

  list(
    source_names = source_names,
    K = K,
    n_folds = n_folds,
    n_t = n_t,
    n_source = n_source,
    N_all = N_all
  )
}

.reweight_tate_inner_component <- function(
    component, target_multiplier, source_multipliers, source_names) {
  caller <- ".reweight_tate_inner_component"
  K <- length(source_names)
  target_idx <- as.integer(component$target_idx)
  source_idx <- component$source_idx
  if (!is.null(names(source_idx))) source_idx <- source_idx[source_names]
  target_w <- target_multiplier[target_idx]

  target_raw <- component$varphi_ot + component$avg_target_est
  target_estimate <- .weight_bootstrap_weighted_mean(
    target_raw, target_w, paste0(caller, " target estimate")
  )
  zeta_raw <- vector("list", K)
  xi_raw <- vector("list", K)
  source_estimate <- mu_pred <- delta <- numeric(K)
  V_t <- V_s <- C_ot <- numeric(K)

  for (j in seq_len(K)) {
    source_w <- source_multipliers[[source_names[[j]]]][
      as.integer(source_idx[[j]])
    ]
    zeta_raw[[j]] <- component$zeta_components[[j]] +
      component$mu_pred_ts[[j]]
    xi_raw[[j]] <- component$xi_components[[j]] + component$delta_ts[[j]]
    mu_pred[[j]] <- .weight_bootstrap_weighted_mean(
      zeta_raw[[j]], target_w, paste0(caller, " target source component")
    )
    delta[[j]] <- .weight_bootstrap_weighted_mean(
      xi_raw[[j]], source_w, paste0(caller, " source residual component")
    )
    source_estimate[[j]] <- mu_pred[[j]] + delta[[j]]
    V_t[[j]] <- .weight_bootstrap_centered_moment(
      zeta_raw[[j]], target_w, paste0(caller, " V_t")
    )
    V_s[[j]] <- .weight_bootstrap_centered_moment(
      xi_raw[[j]], source_w, paste0(caller, " V_s")
    )
    C_ot[[j]] <- .weight_bootstrap_covariance(
      target_raw, zeta_raw[[j]], target_w, paste0(caller, " C_ot")
    )
  }

  C_cross <- matrix(0, nrow = K, ncol = K)
  if (K > 1L) {
    for (j in seq_len(K - 1L)) {
      for (l in (j + 1L):K) {
        value <- .weight_bootstrap_covariance(
          zeta_raw[[j]], zeta_raw[[l]], target_w,
          paste0(caller, " C_cross")
        )
        C_cross[j, l] <- C_cross[l, j] <- value
      }
    }
  }

  list(
    outer_fold = component$outer_fold,
    inner_fold = component$inner_fold,
    target_idx = component$target_idx,
    source_idx = component$source_idx,
    V_ot = .weight_bootstrap_centered_moment(
      target_raw, target_w, paste0(caller, " V_ot")
    ),
    V_t = V_t,
    V_s = V_s,
    C_ot = C_ot,
    C_cross = C_cross,
    avg_target_est = target_estimate,
    avg_source_est = source_estimate,
    mu_pred_ts = mu_pred,
    delta_ts = delta,
    n_t = component$n_t,
    n_s = component$n_s,
    varphi_ot = target_raw - target_estimate,
    zeta_components = lapply(seq_len(K), function(j) {
      zeta_raw[[j]] - mu_pred[[j]]
    }),
    xi_components = lapply(seq_len(K), function(j) {
      xi_raw[[j]] - delta[[j]]
    })
  )
}

.reweight_tate_inner_info <- function(
    inner_info, target_multiplier, source_multipliers, source_names) {
  components <- lapply(inner_info$components, function(component) {
    .reweight_tate_inner_component(
      component, target_multiplier, source_multipliers, source_names
    )
  })
  names(components) <- names(inner_info$components)
  K <- length(source_names)
  target_raw <- unlist(lapply(components, function(component) {
    component$varphi_ot + component$avg_target_est
  }), use.names = FALSE)
  target_w <- unlist(lapply(components, function(component) {
    target_multiplier[as.integer(component$target_idx)]
  }), use.names = FALSE)
  zeta_raw <- lapply(seq_len(K), function(j) {
    unlist(lapply(components, function(component) {
      component$zeta_components[[j]] + component$mu_pred_ts[[j]]
    }), use.names = FALSE)
  })
  xi_raw <- lapply(seq_len(K), function(j) {
    unlist(lapply(components, function(component) {
      component$xi_components[[j]] + component$delta_ts[[j]]
    }), use.names = FALSE)
  })
  source_w <- lapply(seq_len(K), function(j) {
    unlist(lapply(components, function(component) {
      idx <- component$source_idx
      if (!is.null(names(idx))) idx <- idx[source_names]
      source_multipliers[[source_names[[j]]]][as.integer(idx[[j]])]
    }), use.names = FALSE)
  })
  target_estimate <- .weight_bootstrap_weighted_mean(
    target_raw, target_w, ".reweight_tate_inner_info target estimate"
  )
  mu_pred <- vapply(seq_len(K), function(j) {
    .weight_bootstrap_weighted_mean(
      zeta_raw[[j]], target_w, ".reweight_tate_inner_info mu_pred"
    )
  }, numeric(1L))
  delta <- vapply(seq_len(K), function(j) {
    .weight_bootstrap_weighted_mean(
      xi_raw[[j]], source_w[[j]], ".reweight_tate_inner_info delta"
    )
  }, numeric(1L))
  # Production Phase 1b pools influence components centered within each
  # inner fold. Keep that convention after reweighting: pooling raw values
  # here would add between-fold location variation to the weight objective
  # and would change the fitted weights even when every multiplier is one.
  target_centered <- unlist(lapply(components, `[[`, "varphi_ot"),
                           use.names = FALSE)
  zeta_centered <- lapply(seq_len(K), function(j) {
    unlist(lapply(components, function(component) {
      component$zeta_components[[j]]
    }), use.names = FALSE)
  })
  xi_centered <- lapply(seq_len(K), function(j) {
    unlist(lapply(components, function(component) {
      component$xi_components[[j]]
    }), use.names = FALSE)
  })
  C_cross <- matrix(0, nrow = K, ncol = K)
  if (K > 1L) {
    for (j in seq_len(K - 1L)) {
      for (l in (j + 1L):K) {
        value <- .weight_bootstrap_covariance(
          zeta_centered[[j]], zeta_centered[[l]], target_w,
          ".reweight_tate_inner_info C_cross"
        )
        C_cross[j, l] <- C_cross[l, j] <- value
      }
    }
  }
  result <- list(
    V_ot = .weight_bootstrap_centered_moment(
      target_centered, target_w, ".reweight_tate_inner_info V_ot"
    ),
    V_t = vapply(seq_len(K), function(j) {
      .weight_bootstrap_centered_moment(
        zeta_centered[[j]], target_w, ".reweight_tate_inner_info V_t"
      )
    }, numeric(1L)),
    V_s = vapply(seq_len(K), function(j) {
      .weight_bootstrap_centered_moment(
        xi_centered[[j]], source_w[[j]], ".reweight_tate_inner_info V_s"
      )
    }, numeric(1L)),
    C_ot = vapply(seq_len(K), function(j) {
      .weight_bootstrap_covariance(
        target_centered, zeta_centered[[j]], target_w,
        ".reweight_tate_inner_info C_ot"
      )
    }, numeric(1L)),
    C_cross = C_cross,
    avg_target_est = target_estimate,
    avg_source_est = mu_pred + delta,
    mu_pred_ts = mu_pred,
    delta_ts = delta,
    n_t = sum(vapply(components, `[[`, numeric(1L), "n_t")),
    n_s = colSums(do.call(rbind, lapply(components, `[[`, "n_s"))),
    varphi_ot = target_centered,
    zeta_components = zeta_centered,
    xi_components = xi_centered
  )
  result$components <- components
  result
}

.weight_bootstrap_outer_estimate <- function(
    fold_info, fold_weights, target_multiplier, source_multipliers,
    source_names, n_t, n_source) {
  K <- length(source_names)
  target_values <- numeric(n_t)
  source_values <- stats::setNames(
    lapply(n_source, numeric), source_names
  )
  target_count <- integer(n_t)
  source_count <- stats::setNames(lapply(n_source, integer), source_names)

  for (k in seq_along(fold_info)) {
    info <- fold_info[[k]]
    eta <- as.numeric(fold_weights[k, ])
    target_raw <- info$varphi_ot + info$fold_target_estimate
    target_psi <- (1 - sum(eta)) * target_raw
    for (j in seq_len(K)) {
      target_psi <- target_psi + eta[[j]] *
        (info$zeta_components[[j]] + info$mu_pred_ts[[j]])
    }
    target_idx <- as.integer(info$target_idx)
    target_count <- target_count + tabulate(target_idx, nbins = n_t)
    target_values[target_idx] <- target_psi

    source_idx <- info$source_idx
    if (!is.null(names(source_idx))) source_idx <- source_idx[source_names]
    for (j in seq_len(K)) {
      idx <- as.integer(source_idx[[j]])
      values <- eta[[j]] *
        (info$xi_components[[j]] + info$delta_ts[[j]])
      source_count[[j]] <- source_count[[j]] +
        tabulate(idx, nbins = n_source[[j]])
      source_values[[j]][idx] <- values
    }
  }
  if (any(target_count != 1L) || any(vapply(source_count, function(x) {
    any(x != 1L)
  }, logical(1L)))) {
    stop(".weight_bootstrap_outer_estimate: outer folds must cover each observation once.",
         call. = FALSE)
  }

  estimate <- .weight_bootstrap_weighted_mean(
    target_values, target_multiplier, ".weight_bootstrap_outer_estimate target"
  )
  for (j in seq_len(K)) {
    estimate <- estimate + .weight_bootstrap_weighted_mean(
      source_values[[j]], source_multipliers[[source_names[[j]]]],
      ".weight_bootstrap_outer_estimate source"
    )
  }
  estimate
}

#' Fixed-nuisance multiplier bootstrap for TATE aggregation weights
#'
#' Uses one global set of site-stratified exponential multipliers per draw.
#' Multipliers are shared across overlapping inner folds, outer evaluation
#' folds, and treatment arms. High-dimensional nuisance fits remain fixed, but
#' all inner aggregation moments, Wald screens, and source weights are relearned.
#'
#' @param tate_result Result from \code{run_tate_crossfit()}.
#' @param B Number of multiplier-bootstrap draws.
#' @param seed Positive integer seed used without perturbing the global RNG.
#' @param relearn_weights If false, compute only the fixed-weight comparison.
#' @param screening_rule Optional soft/hard rule for relearned weights. The
#'   fitted result's rule is used when omitted.
#'   Only fits using one fixed numeric aggregation multiplier are supported;
#'   aggregation-lambda inner CV is rejected because its weighted-bootstrap
#'   analogue is not implemented.
#' @param save_draws Retain draw-level estimates and weights.
#' @param verbose Print progress every 100 draws.
#' @return Bootstrap variance diagnostics for the fitted TATE estimator.
#' @export
estimate_tate_weight_bootstrap <- function(
    tate_result, B = 500L, seed = 1L, relearn_weights = TRUE,
    screening_rule = NULL, save_draws = FALSE, verbose = FALSE) {
  structure <- .validate_weight_bootstrap_result(tate_result)
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || !is.finite(B) ||
      B < 2 || B != floor(B) || B > .Machine$integer.max) {
    stop("estimate_tate_weight_bootstrap: B must be one integer in [2, integer.max].",
         call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 1 || seed != floor(seed) ||
      seed > .Machine$integer.max) {
    stop("estimate_tate_weight_bootstrap: seed must be one positive integer.",
         call. = FALSE)
  }
  B <- as.integer(B)
  seed <- as.integer(seed)
  cube_cells_per_draw <- as.double(structure$n_folds) * structure$K
  if (isTRUE(relearn_weights) &&
      as.double(B) * cube_cells_per_draw > .Machine$integer.max) {
    stop(paste0(
      "estimate_tate_weight_bootstrap: B * n_folds * K exceeds the safe ",
      "R array dimension limit."
    ), call. = FALSE)
  }
  if (length(relearn_weights) != 1L || is.na(relearn_weights) ||
      !is.logical(relearn_weights)) {
    stop("estimate_tate_weight_bootstrap: relearn_weights must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(save_draws) != 1L || is.na(save_draws) ||
      !is.logical(save_draws)) {
    stop("estimate_tate_weight_bootstrap: save_draws must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (length(verbose) != 1L || is.na(verbose) || !is.logical(verbose)) {
    stop("estimate_tate_weight_bootstrap: verbose must be TRUE or FALSE.",
         call. = FALSE)
  }
  if (is.null(screening_rule)) {
    screening_rule <- tate_result$aggregation_screening_rule
  }
  screening_rule <- match.arg(
    screening_rule, c("soft_penalty", "hard_threshold", "quadratic_bias")
  )

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  fixed_draws <- numeric(B)
  relearned_draws <- if (isTRUE(relearn_weights)) numeric(B) else NULL
  canonical_sources <- sort(structure$source_names)
  canonical_order <- match(canonical_sources, structure$source_names)
  canonical_n_source <- structure$n_source[canonical_sources]
  fold_info <- lapply(
    tate_result$intermediates$fold_info,
    .permute_weight_bootstrap_fold_info,
    order_index = canonical_order
  )
  inner_fold_info <- lapply(tate_result$intermediates$inner_fold_info,
                            function(inner) {
    inner$components <- lapply(
      inner$components, .permute_weight_bootstrap_component,
      order_index = canonical_order
    )
    inner
  })
  fixed_weights <- as.matrix(tate_result$fold_weights)[
    , canonical_order, drop = FALSE
  ]

  weight_draws <- if (isTRUE(relearn_weights)) {
    array(NA_real_, dim = c(B, structure$n_folds, structure$K),
          dimnames = list(NULL, NULL, canonical_sources))
  } else NULL
  inclusion_draws <- if (isTRUE(relearn_weights)) {
    array(NA, dim = c(B, structure$n_folds, structure$K),
          dimnames = list(NULL, NULL, canonical_sources))
  } else NULL
  wald_draws <- if (isTRUE(relearn_weights)) {
    array(NA_real_, dim = c(B, structure$n_folds, structure$K),
          dimnames = list(NULL, NULL, canonical_sources))
  } else NULL
  for (b in seq_len(B)) {
    target_multiplier <- stats::rexp(structure$n_t)
    source_multipliers <- stats::setNames(lapply(
      canonical_n_source, stats::rexp
    ), canonical_sources)

    fixed_draws[[b]] <- .weight_bootstrap_outer_estimate(
      fold_info, fixed_weights, target_multiplier, source_multipliers,
      canonical_sources, structure$n_t, canonical_n_source
    )
    if (isTRUE(relearn_weights)) {
      reweighted_inner <- lapply(inner_fold_info, function(inner) {
        .reweight_tate_inner_info(
          inner, target_multiplier, source_multipliers,
          canonical_sources
        )
      })
      phase2 <- .compute_phase2_weights(
        n_folds = structure$n_folds,
        inner_fold_info = reweighted_inner,
        fold_info = fold_info,
        lambda_selection = tate_result$aggregation_lambda_selection,
        K = structure$K,
        verbose = FALSE,
        lambda_rule = tate_result$aggregation_lambda_rule,
        aggregation_lambda_grid = tate_result$aggregation_lambda_grid %||% NULL,
        screening_rule = screening_rule
      )
      colnames(phase2$fold_weights) <- canonical_sources
      relearned_draws[[b]] <- .weight_bootstrap_outer_estimate(
        fold_info, phase2$fold_weights, target_multiplier,
        source_multipliers, canonical_sources, structure$n_t,
        canonical_n_source
      )
      weight_draws[b, , ] <- phase2$fold_weights
      inclusion_draws[b, , ] <- phase2$fold_source_included
      wald_draws[b, , ] <- phase2$fold_wald_statistics
    }
    if (isTRUE(verbose) && (b %% 100L == 0L || b == B)) {
      message("estimate_tate_weight_bootstrap: completed ", b, "/", B,
              " draws")
    }
  }

  summarize_cube <- function(x, fun) {
    stats::setNames(vapply(seq_len(structure$K), function(j) {
      fun(as.numeric(x[, , j]))
    }, numeric(1L)), canonical_sources)
  }
  fixed_se <- stats::sd(fixed_draws)
  result <- list(
    estimate = as.numeric(tate_result$estimate),
    se_analytic = as.numeric(tate_result$se),
    variance_analytic = as.numeric(tate_result$se)^2,
    se_fixed_weight_bootstrap = fixed_se,
    variance_fixed_weight_bootstrap = fixed_se^2,
    n_bootstrap = B,
    bootstrap_multiplier = "site_stratified_exponential",
    bootstrap_seed = seed,
    nuisance_refit = FALSE,
    weights_relearned = isTRUE(relearn_weights),
    screening_rule = screening_rule,
    bootstrap_failures = 0L
  )
  if (isTRUE(relearn_weights)) {
    relearned_se <- stats::sd(relearned_draws)
    result$se_weight_relearn_bootstrap <- relearned_se
    result$variance_weight_relearn_bootstrap <- relearned_se^2
    result$weight_uncertainty_ratio <- relearned_se / fixed_se
    result$bootstrap_weight_mean <- summarize_cube(weight_draws, mean)
    result$bootstrap_weight_sd <- summarize_cube(weight_draws, stats::sd)
    result$bootstrap_inclusion_probability <- summarize_cube(
      inclusion_draws, mean
    )
    result$bootstrap_wald_mean <- summarize_cube(wald_draws, mean)
  }
  if (isTRUE(save_draws)) {
    result$draws <- list(
      fixed_weight_estimate = fixed_draws,
      weight_relearn_estimate = relearned_draws,
      fold_weights = weight_draws,
      fold_source_included = inclusion_draws,
      fold_wald_statistics = wald_draws
    )
  }
  result
}
