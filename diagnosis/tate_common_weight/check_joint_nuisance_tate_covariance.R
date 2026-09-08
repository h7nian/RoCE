#!/usr/bin/env Rscript

# Deterministic enumeration of actual joint (X,A,Y) outcomes, not independent
# artificial arm samples. Analytic score correction, no bootstrap/SE scaling.
.population_site_score <- function(details, events, theta = details$theta) {
  W <- details$W[events$x_index, , drop = FALSE]
  Z <- details$Z[events$x_index, , drop = FALSE]
  index <- details$indices
  initial_mean <- plogis(drop(W %*% theta[index$alpha_initial]))
  final_mean <- plogis(drop(W %*% theta[index$alpha_final]))
  initial_weight <- exp(-drop(Z %*% theta[index$gamma_initial]))
  final_weight <- exp(-drop(Z %*% theta[index$gamma_final]))
  derivative <- initial_mean*(1-initial_mean)
  in_arm <- as.numeric(events$A == details$arm_value)
  zero_outcome <- matrix(0, nrow(events), ncol(W))
  if (identical(unique(events$site), "target")) {
    moments <- cbind(W*(in_arm*(events$Y-initial_mean)), Z,
                     zero_outcome, Z*derivative)
    original_score <- final_mean
  } else if (identical(unique(events$site), "source")) {
    moments <- cbind(zero_outcome, -Z*(in_arm*initial_weight),
      W*(in_arm*initial_weight*(events$Y-final_mean)),
      -Z*(in_arm*final_weight*derivative))
    original_score <- in_arm*final_weight*(events$Y-final_mean)
  } else stop("expected exactly one known site")
  stopifnot(ncol(moments) == length(details$theta))
  list(score = original_score-drop(moments %*% details$projection), moments = moments)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% c(1L, 2L) ||
      (length(args) == 2L && args[[2L]] != "--sparse-equations")) {
    stop("usage: check_joint_nuisance_tate_covariance.R V19_LIBRARY [--sparse-equations]")
  }
  projection_solver <- NULL
  projection_calls <- 0L
  if (length(args) == 2L) {
    source("diagnosis/tate_common_weight/sparse_moment_projection.R")
    projection_solver <- function(jacobian, gradient, indices) {
      projection_calls <<- projection_calls+1L
      signs <- c(-1, 1, -1, 1)
      curvatures <- setNames(lapply(seq_along(indices), function(k) {
        index <- indices[[k]]
        signs[k]*jacobian[index, index, drop = FALSE]
      }), names(indices))
      couplings <- list(
        alpha_final_gamma_initial = jacobian[indices$alpha_final, indices$gamma_initial, drop = FALSE],
        gamma_final_alpha_initial = jacobian[indices$gamma_final, indices$alpha_initial, drop = FALSE])
      gradients <- list(alpha_final = gradient[indices$alpha_final],
                        gamma_final = gradient[indices$gamma_final])
      penalties <- setNames(as.list(rep(0, 4)), names(indices))
      fitted <- .solve_joint_nuisance_projection(curvatures, couplings, gradients,
                                                 penalties, tolerance = 1e-11)
      stopifnot(max(abs(drop(t(jacobian) %*% fitted$coefficients)-gradient)) < 1e-8)
      fitted$coefficients
    }
  }
  fixture <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/check_full_nuisance_orthogonalization.R", fixture)
  target_probability <- c(.3, .5, .7)
  treated <- fixture$main(args[1L], arm_value = 1L,
    target_treatment_probability = target_probability, return_details = TRUE, verbose = FALSE,
    projection_solver = projection_solver)
  control <- fixture$main(args[1L], arm_value = 0L,
    target_treatment_probability = target_probability, return_details = TRUE, verbose = FALSE,
    projection_solver = projection_solver)
  if (!is.null(projection_solver)) stopifnot(projection_calls == 6L)
  treated <- attr(treated, "case_details")
  control <- attr(control, "case_details")
  support <- c(-1, 0, 1)
  mean_treated <- plogis(-.2+.7*support+1.2*support^2)
  mean_control <- plogis(-1+.7*support+1.2*support^2)
  target_truth <- mean(mean_treated-mean_control)
  events <- lapply(c("target", "source"), function(site) {
    rows <- expand.grid(x_index = 1:3, A = 0:1, Y = 0:1)
    rows$site <- site
    x_mass <- if (site == "target") rep(1/3, 3) else c(.3, .4, .3)
    propensity <- if (site == "target") target_probability else
      c(.08, .28, .12)/x_mass
    arm_mean <- ifelse(rows$A == 1L, mean_treated[rows$x_index], mean_control[rows$x_index])
    treatment_mass <- ifelse(rows$A == 1L, propensity[rows$x_index], 1-propensity[rows$x_index])
    rows$probability <- x_mass[rows$x_index]*treatment_mass*
      ifelse(rows$Y == 1L, arm_mean, 1-arm_mean)
    stopifnot(all(rows$probability > 0), abs(sum(rows$probability)-1) < 1e-12)
    rows
  })
  names(events) <- c("target", "source")
  # Unequal sizes deliberately test the multisample/site-scaling convention.
  sample_sizes <- c(target = 1000, source = 750)
  weighted_cov <- function(x, y, probability) {
    sum(probability*(x-sum(probability*x))*(y-sum(probability*y)))
  }
  result <- lapply(names(treated), function(label) {
    arm1 <- treated[[label]]; arm0 <- control[[label]]
    score1 <- lapply(events, function(x) .population_site_score(arm1, x))
    score0 <- lapply(events, function(x) .population_site_score(arm0, x))
    for (pair in list(list(details = arm1, scores = score1), list(details = arm0, scores = score0))) {
      observed_moments <- Reduce(`+`, lapply(names(events), function(site)
        colSums(pair$scores[[site]]$moments*events[[site]]$probability)))
      stopifnot(max(abs(observed_moments-pair$details$population_moments)) < 1e-12)
    }
    site_results <- lapply(names(events), function(site) {
      probability <- events[[site]]$probability
      u <- score1[[site]]$score; v <- score0[[site]]$score
      contrast <- u-v
      data.frame(site = site, n = sample_sizes[[site]],
        mean_contrast = sum(probability*contrast),
        arm1_variance = weighted_cov(u, u, probability),
        arm0_variance = weighted_cov(v, v, probability),
        arm_covariance = weighted_cov(u, v, probability),
        contrast_variance = weighted_cov(contrast, contrast, probability))
    })
    sites <- do.call(rbind, site_results)
    mean_error <- abs(sum(sites$mean_contrast)-target_truth)
    direct_variance <- sum(sites$contrast_variance/sites$n)
    covariance_variance <- sum((sites$arm1_variance+sites$arm0_variance-2*sites$arm_covariance)/sites$n)
    naive_variance <- sum((sites$arm1_variance+sites$arm0_variance)/sites$n)
    total_n <- sum(sample_sizes)
    scaled_variance <- sum((sites$n/total_n)*(total_n/sites$n)^2*sites$contrast_variance)/total_n
    joint_theta <- c(arm1$theta, arm0$theta)
    joint_expectation <- function(theta) {
      theta1 <- theta[seq_along(arm1$theta)]
      theta0 <- theta[length(arm1$theta)+seq_along(arm0$theta)]
      sum(vapply(names(events), function(site) {
        u <- .population_site_score(arm1, events[[site]], theta1)$score
        v <- .population_site_score(arm0, events[[site]], theta0)$score
        sum(events[[site]]$probability*(u-v))
      }, numeric(1)))
    }
    gradient <- vapply(seq_along(joint_theta), function(j) {
      shift <- numeric(length(joint_theta)); shift[j] <- 1e-5
      (joint_expectation(joint_theta+shift)-joint_expectation(joint_theta-shift))/2e-5
    }, numeric(1))
    stopifnot(mean_error < 1e-8, max(abs(gradient)) < 1e-8,
      abs(direct_variance-covariance_variance) < 1e-12,
      abs(direct_variance-scaled_variance) < 1e-12,
      abs(naive_variance-direct_variance) > 1e-8, direct_variance > 0)
    data.frame(case = label, target_truth, target_identity_error = mean_error,
      full_tate_nuisance_gradient = max(abs(gradient)),
      analytic_variance = direct_variance,
      arm_covariance_identity_error = abs(direct_variance-covariance_variance),
      site_scaling_identity_error = abs(direct_variance-scaled_variance),
      naive_independent_arm_variance = naive_variance,
      target_arm_covariance = sites$arm_covariance[sites$site == "target"],
      source_arm_covariance = sites$arm_covariance[sites$site == "source"])
  })
  result <- do.call(rbind, result)
  cat("projection_solver=", if (is.null(projection_solver)) "dense_reference" else
        "blockwise_equations_zero_penalty_limit", "\n", sep = "")
  print(result, row.names = FALSE, digits = 12)
  cat("Exact joint-observation enumeration; no new Monte Carlo data.\n",
      "One target/source pair only. High-dimensional estimation, target anchor,\n",
      "adaptive common weights and full TATE coverage remain unvalidated.\n", sep = "")
  invisible(result)
}

if (sys.nframe() == 0L) main()
