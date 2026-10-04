# A bounded, sparse simulation with independently sampled sites. The source
# joint laws make merged-weight correctness explicit; target PS is specified
# separately. This is a new DGP, not a replacement for archived FACE draws.

BOUNDED_DGP_VERSION <- "bounded_joint_v3"

.bounded_dgp_control <- function(control = NULL) {
  defaults <- list(covariate_shift = 1, source_treatment_scale = 1)
  if (is.null(control) || (is.list(control) && length(control) == 0L)) return(defaults)
  if (!is.list(control) || is.null(names(control)) || anyNA(names(control)) || anyDuplicated(names(control)) ||
      any(!names(control) %in% names(defaults))) {
    stop("dgp_control supports covariate_shift and source_treatment_scale.", call. = FALSE)
  }
  for (name in names(control)) {
    value <- control[[name]]
    upper <- if (name == "source_treatment_scale") 1 else 2
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value < 0 || value > upper) {
      stop("dgp_control$", name, " must be one number in [0, ", upper, "].", call. = FALSE)
    }
    defaults[[name]] <- as.numeric(value)
  }
  defaults
}

.bounded_covariate_spec <- function() {
  limit <- 2.5
  probability <- 2 * stats::pnorm(limit) - 1
  scale <- sqrt(1 - 2 * limit * stats::dnorm(limit) / probability)
  list(limit = limit, probability = probability, scale = scale, bound = limit / scale)
}

.draw_bounded_covariates <- function(n, p) {
  spec <- .bounded_covariate_spec()
  probabilities <- stats::pnorm(-spec$limit) + stats::runif(n * p) * spec$probability
  matrix(stats::qnorm(probabilities) / spec$scale, nrow = n, ncol = p)
}

.bounded_active_features <- function(X, strength = 0) {
  active <- X[, 1:4, drop = FALSE]
  interaction <- strength
  active[, 1L] <- active[, 1L] * (1 + interaction * active[, 2L]) / sqrt(1 + interaction^2)
  active
}

.bounded_quadrature <- function(order = 24L) {
  index <- seq_len(order - 1L)
  jacobi <- matrix(0, order, order)
  off_diagonal <- index / sqrt(4 * index^2 - 1)
  jacobi[cbind(index, index + 1L)] <- off_diagonal
  jacobi[cbind(index + 1L, index)] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  spec <- .bounded_covariate_spec()
  z <- spec$limit * decomposition$values
  weights <- decomposition$vectors[1L, ]^2 * 2 * spec$limit * stats::dnorm(z) / spec$probability
  grid <- as.matrix(expand.grid(rep(list(seq_len(order)), 4L)))
  list(X = matrix(z[grid] / spec$scale, ncol = 4L),
       weights = weights[grid[, 1L]] * weights[grid[, 2L]] * weights[grid[, 3L]] * weights[grid[, 4L]])
}

.bounded_predictor_range <- function(coefficients, strength) {
  bound <- .bounded_covariate_spec()$bound
  corners <- as.matrix(expand.grid(rep(list(c(-bound, bound)), 4L)))
  # The predictor is affine in each coordinate separately, so its extrema
  # on this rectangular support occur at corners, including with interaction.
  range(drop(.bounded_active_features(corners, strength) %*% coefficients))
}

.draw_bounded_tilt <- function(n, coefficients, strength) {
  result <- matrix(0, n, 4L)
  if (n == 0L) return(result)
  maximum <- .bounded_predictor_range(coefficients, strength)[2L]
  filled <- 0L
  while (filled < n) {
    proposal <- .draw_bounded_covariates(max(32L, 2L * (n - filled)), 4L)
    log_acceptance <- drop(.bounded_active_features(proposal, strength) %*% coefficients) - maximum
    if (any(log_acceptance > 1e-10)) stop("Bounded tilt envelope is invalid.", call. = FALSE)
    accepted <- which(log(stats::runif(nrow(proposal))) <= log_acceptance)
    count <- min(length(accepted), n - filled)
    if (count > 0L) {
      result[filled + seq_len(count), ] <- proposal[accepted[seq_len(count)], , drop = FALSE]
      filled <- filled + count
    }
  }
  result
}

#' Generate a bounded sparse multisite experiment
#'
#' All p working features are directly observed standardized truncated-normal
#' coordinates. Four slopes are active; the intercept is additional. Sources
#' follow joint exponential tilts of the target law. A common tilt controls
#' covariate differences, while an arm contrast controls source treatment.
#' Zero common tilt alone does not imply identical site covariate laws when
#' the arm contrast is nonzero. Both controls zero give identical laws.
#'
#' C1 has correct outcome, target PS and source merged-weight models. C2 adds
#' an omitted outcome interaction; C3 adds an omitted assignment interaction.
#' C4 adds both. The bounded interaction depends only on the active coordinates,
#' preserving an independent centered inactive block in every site and arm.
#' The coefficient values are fixed before performance evaluation. This is a
#' literature-informed adaptation, not a reproduction of a published design.
#'
#' @inheritParams generate_face_data
#' @param dgp_control Optional list with covariate_shift and
#'   source_treatment_scale, in [0, 2] and [0, 1] respectively. Defaults are both 1.
#' @return A simulation dataset with oracle means, joint weights and explicit
#'   DGP version/control metadata. Only binary outcomes are supported.
#' @export
generate_bounded_data <- function(n_total = NULL, K = 2L, p = 10L, config = "C1",
                                  estimand_type = "superpopulation", outcome_type = "binary",
                                  ate_deviation = 0, n_deviated_sites = 0L,
                                  deviation_mechanism = c("treated_arm", "both_arms"),
                                  n_target = NULL, n_source_sizes = NULL,
                                  misspecification_strength = 0.75, dgp_control = NULL) {
  control <- .bounded_dgp_control(dgp_control)
  deviation_mechanism <- match.arg(deviation_mechanism)
  config <- match.arg(config, c("C1", "C2", "C3", "C4"))
  estimand_type <- match.arg(estimand_type, c("superpopulation", "sample"))
  if (!identical(outcome_type, "binary")) stop("The bounded DGP supports binary outcomes only.", call. = FALSE)
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p) || p < 4L || p != floor(p)) {
    stop("The bounded DGP requires integer p >= 4.", call. = FALSE)
  }
  if (length(ate_deviation) != 1L || !is.finite(ate_deviation) || ate_deviation < 0 ||
      length(n_deviated_sites) != 1L || !is.finite(n_deviated_sites) || n_deviated_sites < 0 ||
      n_deviated_sites != floor(n_deviated_sites)) stop("Invalid source outcome deviation.", call. = FALSE)
  strengths <- .face_misspecification_strengths(config, misspecification_strength)
  sizes <- resolve_face_site_sizes(n_total, n_target, n_source_sizes, K)
  K <- sizes$K
  n <- sizes$n_total
  site_names <- c("t", paste0("s", seq_len(K)))
  R <- rep(site_names, c(sizes$n_target, sizes$n_source_sizes))
  # Fixed uniforms pair treatment/outcome draws across source-shift settings.
  # In particular, changing source controls preserves the realized target data.
  treatment_uniform <- stats::runif(n)
  outcome0_uniform <- stats::runif(n)
  outcome1_uniform <- stats::runif(n)
  X <- .draw_bounded_covariates(n, p)
  colnames(X) <- paste0("X", seq_len(p))
  outcome_slopes <- c(.25, .20, .10, .05)
  propensity_slopes <- c(.35, -.25, .125, -.0625)
  shift_slopes <- c(.25, -.20, .15, -.10)
  quadrature <- .bounded_quadrature()
  assignment_grid <- .bounded_active_features(quadrature$X, strengths$propensity)
  source_parameters <- setNames(vector("list", K), site_names[-1L])
  for (j in seq_len(K)) {
    common <- control$covariate_shift * (j / K) * shift_slopes
    contrast <- control$source_treatment_scale * propensity_slopes / 2
    slopes <- cbind(mu0 = common - contrast, mu1 = common + contrast)
    normalizers <- colSums(exp(assignment_grid %*% slopes) * quadrature$weights)
    total <- sum(normalizers)
    rows <- which(R == site_names[j + 1L])
    # This component label samples the marginal covariate mixture. Actual
    # treatment is drawn afterward from its conditional propensity.
    component <- stats::rbinom(length(rows), 1L, normalizers[2L] / total)
    for (arm in 0:1) {
      selected <- rows[component == arm]
      X[selected, 1:4] <- .draw_bounded_tilt(length(selected), slopes[, arm + 1L], strengths$propensity)
    }
    source_parameters[[j]] <- list(slopes = slopes, normalizers = normalizers,
                                  log_normalizer = log(total))
  }
  assignment <- .bounded_active_features(X, strengths$propensity)
  propensity_predictor <- drop(assignment %*% propensity_slopes)
  propensity_predictor[R != "t"] <- propensity_predictor[R != "t"] * control$source_treatment_scale
  propensity <- stats::plogis(propensity_predictor)
  A <- as.integer(treatment_uniform < propensity)
  outcome_predictor <- drop(.bounded_active_features(X, strengths$outcome) %*% outcome_slopes)
  deviated <- R %in% utils::head(site_names[-1L], n_deviated_sites)
  shared_shift <- if (deviation_mechanism == "both_arms") ate_deviation * deviated else 0
  treatment_shift <- 1 + if (deviation_mechanism == "treated_arm") ate_deviation * deviated else 0
  mean0 <- stats::plogis(-.5 + outcome_predictor + shared_shift)
  mean1 <- stats::plogis(-.5 + outcome_predictor + shared_shift + treatment_shift)
  Y0 <- as.integer(outcome0_uniform < mean0)
  Y1 <- as.integer(outcome1_uniform < mean1)
  target <- R == "t"
  grid_predictor <- drop(.bounded_active_features(quadrature$X, strengths$outcome) %*% outcome_slopes)
  truth0 <- sum(quadrature$weights * stats::plogis(-.5 + grid_predictor))
  truth1 <- sum(quadrature$weights * stats::plogis(.5 + grid_predictor))
  oracle_weights <- matrix(NA_real_, n, 2L, dimnames = list(NULL, c("mu0", "mu1")))
  gamma_parameters <- list()
  for (site in names(source_parameters)) {
    rows <- which(R == site)
    model <- source_parameters[[site]]
    oracle_weights[rows, ] <- exp(model$log_normalizer - assignment[rows, , drop = FALSE] %*% model$slopes)
    for (arm in 0:1) gamma_parameters[[paste0(site, "_", arm)]] <-
      c(-model$log_normalizer, model$slopes[, arm + 1L], rep(0, p - 4L))
  }
  true_outcome_design <- true_assignment_design <- X
  true_outcome_design[, 1:4] <- .bounded_active_features(X, strengths$outcome)
  true_assignment_design[, 1:4] <- assignment
  list(n = n, K = K, p = p, X = X, X_dagger = true_outcome_design, R = R, A = A,
    Y = ifelse(A == 1L, Y1, Y0), Y_1 = Y1, Y_0 = Y0,
    p_treat_true = propensity, mu0_true = if (estimand_type == "sample") mean(mean0[target]) else truth0,
    mu1_true = if (estimand_type == "sample") mean(mean1[target]) else truth1,
    mu0_realized = mean(Y0[target]), mu1_realized = mean(Y1[target]),
    gamma_params = gamma_parameters, alpha0_true = c(-.5, outcome_slopes, rep(0, p - 4L)),
    alpha1_true = c(.5, outcome_slopes, rep(0, p - 4L)),
    Z_site = X, W_outcome = X, Z_site_true = true_assignment_design,
    W_outcome_true = true_outcome_design, config = config, dgp_type = "bounded",
    outcome_type = "binary", estimand_type = estimand_type,
    misspecification_strength = if (config == "C1") 0 else misspecification_strength,
    deviation_mechanism = deviation_mechanism, dgp_version = BOUNDED_DGP_VERSION,
    dgp_control = control, max_signal_slopes = 4L,
    truth_method = if (estimand_type == "sample") "sample_conditional_mean" else "gauss_legendre_24",
    oracle = list(mu0 = mean0, mu1 = mean1, source_weights = oracle_weights,
                  source_parameters = source_parameters))
}
