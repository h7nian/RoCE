# Diagnosis-only helpers for interpreting FACE-HD DGP gamma parameters.
#
# data$gamma_params stores multinomial-logit coefficients for
#   log P(R=s_j,A=a | X) / P(R=t | X).
# Source-assisted diagnostics average over the conditional source-site sample
# R=s_j and use I(A=a) exp(-phi^T gamma).  The corresponding calibration
# coefficient therefore needs an intercept shift log(P(R=t) / P(R=s_j)).

c2_source_index <- function(source_name) {
  source_idx <- suppressWarnings(as.integer(sub("^s", "", as.character(source_name))))
  if (!is.finite(source_idx) || is.na(source_idx) || source_idx <= 0L) {
    stop(sprintf("Invalid source name '%s'. Expected names like 's1'.", source_name),
         call. = FALSE)
  }
  source_idx
}

c2_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

c2_facehd_function <- function(name) {
  if (exists(name, mode = "function", inherits = TRUE)) {
    return(get(name, mode = "function", inherits = TRUE))
  }
  ns <- asNamespace("FACEHD")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(sprintf("FACEHD namespace does not contain required function '%s'.", name),
         call. = FALSE)
  }
  get(name, envir = ns, inherits = FALSE)
}

c2_true_source_joint_gamma <- function(data, source_name, A_val = 1L) {
  source_idx <- c2_source_index(source_name)
  gamma_key <- paste0("s", source_idx, "_", A_val)
  gamma <- as.numeric(data$gamma_params[[gamma_key]])
  if (length(gamma) == 0L) {
    stop(sprintf("Missing true joint gamma key '%s'.", gamma_key), call. = FALSE)
  }
  gamma
}

c2_true_source_log_target_over_source <- function(data, source_name,
                                                  ratio = c("population", "sample"),
                                                  target_n = NULL,
                                                  source_n = NULL) {
  ratio <- match.arg(ratio)
  source_idx <- c2_source_index(source_name)

  if (identical(ratio, "sample")) {
    if (is.null(target_n)) {
      target_n <- sum(data$R == "t")
    }
    if (is.null(source_n)) {
      source_n <- sum(data$R == paste0("s", source_idx))
    }
    if (!is.finite(target_n) || !is.finite(source_n) ||
        target_n <= 0L || source_n <= 0L) {
      stop("Sample true-gamma normalization requires positive target_n and source_n.",
           call. = FALSE)
    }
    return(log(as.numeric(target_n) / as.numeric(source_n)))
  }

  calculate_site_probabilities <- c2_facehd_function("calculate_site_probabilities")
  Z_true <- as.matrix(c2_null_coalesce(data$Z_site_true, data$Z_site))
  probs <- calculate_site_probabilities(Z_true, data$gamma_params, data$K)
  p_target <- mean(probs$p_target)
  p_source <- mean(probs[[paste0("p_s", source_idx, "_0")]] +
                     probs[[paste0("p_s", source_idx, "_1")]])
  if (!is.finite(p_target) || !is.finite(p_source) ||
      p_target <= 0 || p_source <= 0) {
    stop("Population true-gamma normalization produced non-positive site probability.",
         call. = FALSE)
  }
  log(p_target / p_source)
}

c2_true_source_calibration_gamma <- function(data, source_name, A_val = 1L,
                                             ratio = c("population", "sample"),
                                             gamma_raw = NULL,
                                             target_n = NULL,
                                             source_n = NULL) {
  gamma <- if (is.null(gamma_raw)) {
    c2_true_source_joint_gamma(data, source_name, A_val)
  } else {
    as.numeric(gamma_raw)
  }
  gamma[1L] <- gamma[1L] + c2_true_source_log_target_over_source(
    data = data,
    source_name = source_name,
    ratio = ratio,
    target_n = target_n,
    source_n = source_n
  )
  gamma
}

c2_true_source_treatment_propensity <- function(data, source_name, A_val = 1L,
                                                rows = NULL) {
  source_idx <- c2_source_index(source_name)
  calculate_site_probabilities <- c2_facehd_function("calculate_site_probabilities")
  Z_true <- as.matrix(c2_null_coalesce(data$Z_site_true, data$Z_site))
  probs <- calculate_site_probabilities(Z_true, data$gamma_params, data$K)

  source_prob <- probs[[paste0("p_s", source_idx, "_0")]] +
    probs[[paste0("p_s", source_idx, "_1")]]
  arm_prob <- probs[[paste0("p_s", source_idx, "_", A_val)]]
  propensity <- arm_prob / pmax(source_prob, .Machine$double.eps)
  if (!is.null(rows)) {
    propensity <- propensity[rows]
  }
  propensity
}
