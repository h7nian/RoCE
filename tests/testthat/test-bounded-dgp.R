test_that("bounded DGP has fixed site sizes and four active slopes", {
  set.seed(7101)
  data <- generate_simulation_data(K = 2L, p = 200L, config = "C1", dgp_type = "bounded",
    n_target = 1000L, n_source_sizes = c(1000L, 1000L), warn_ignored = FALSE)
  expect_identical(as.integer(table(data$R)), rep(1000L, 3L))
  expect_equal(dim(data$W_outcome), c(3000L, 200L))
  expect_identical(data$Z_site, data$W_outcome)
  expect_true(all(abs(data$X) <= .bounded_covariate_spec()$bound))
  expect_equal(data$Y, data$A * data$Y_1 + (1 - data$A) * data$Y_0)
  expect_equal(sum(data$alpha0_true[-1L] != 0), 4L)
  expect_equal(sum(data$alpha1_true[-1L] != 0), 4L)
  expect_true(all(vapply(data$gamma_params, function(x) sum(x[-1L] != 0) <= 4L, logical(1L))))
  expect_true(all(data$p_treat_true > 0 & data$p_treat_true < 1))
  expect_equal(data$oracle$mu0,
    plogis(drop(cbind(1, data$W_outcome_true) %*% data$alpha0_true)))
  for (site in c("s1", "s2")) for (arm in 0:1) {
    rows <- data$R == site
    reconstructed <- exp(-drop(cbind(1, data$Z_site_true[rows, ]) %*%
      data$gamma_params[[paste0(site, "_", arm)]]))
    expect_equal(data$oracle$source_weights[rows, arm + 1L], reconstructed, tolerance = 1e-12)
  }
})

test_that("bounded source laws normalize and have the advertised joint weights", {
  coarse <- .bounded_quadrature(16L)
  fine <- .bounded_quadrature(24L)
  expect_equal(sum(fine$weights), 1, tolerance = 1e-12)
  expect_equal(colSums(fine$X * fine$weights), rep(0, 4L), tolerance = 1e-12)
  expect_equal(crossprod(fine$X, fine$X * fine$weights), diag(4L), tolerance = 1e-12)
  for (config in c("C1", "C2", "C3")) {
    set.seed(7111)
    data <- generate_bounded_data(K = 1L, p = 4L, config = config,
      n_target = 1000L, n_source_sizes = 1000L,
      dgp_control = list(covariate_shift = 2))
    strengths <- .face_misspecification_strengths(config, .75)
    features <- .bounded_active_features(fine$X, strengths$propensity)
    model <- data$oracle$source_parameters$s1
    joint <- exp(features %*% model$slopes - model$log_normalizer)
    expect_equal(sum(fine$weights * rowSums(joint)), 1, tolerance = 1e-12)
    e <- plogis(drop(features %*% c(.35, -.25, .125, -.0625)))
    expect_equal(joint[, 2L] / rowSums(joint), e, tolerance = 1e-12)
    for (arm in 1:2) {
      q <- exp(model$log_normalizer - drop(features %*% model$slopes[, arm]))
      expect_equal(colSums(cbind(1, fine$X) * (fine$weights * joint[, arm] * q)),
        c(1, rep(0, 4L)), tolerance = 1e-12)
    }
    outcome <- .bounded_active_features(coarse$X, strengths$outcome)
    eta <- drop(outcome %*% c(.25, .20, .10, .05))
    expect_equal(data$mu1_true, sum(coarse$weights * plogis(.5 + eta)), tolerance = 1e-9)
    expect_equal(data$mu0_true, sum(coarse$weights * plogis(-.5 + eta)), tolerance = 1e-9)
  }
})

test_that("source shift controls preserve target draws and population truth", {
  make_data <- function(control, rho = 0, mechanism = "treated_arm") {
    set.seed(7121)
    generate_bounded_data(K = 2L, p = 8L, config = "C3", n_target = 1000L,
      n_source_sizes = c(1000L, 1000L), dgp_control = control,
      ate_deviation = rho, n_deviated_sites = 1L, deviation_mechanism = mechanism)
  }
  original <- make_data(NULL)
  changed <- make_data(list(covariate_shift = 2, source_treatment_scale = .5))
  target <- original$R == "t"
  for (field in c("A", "Y", "Y_0", "Y_1", "p_treat_true")) {
    expect_identical(original[[field]][target], changed[[field]][target])
  }
  expect_identical(original$X[target, ], changed$X[target, ])
  expect_identical(original$mu1_true, changed$mu1_true)
  expect_identical(original$mu0_true, changed$mu0_true)
  for (mechanism in c("treated_arm", "both_arms")) {
    deviated <- make_data(NULL, 2.5, mechanism)
    expect_identical(original$X, deviated$X)
    expect_identical(original$A, deviated$A)
    expect_identical(original$Y[target], deviated$Y[target])
    expect_identical(original$Y[original$R == "s2"], deviated$Y[original$R == "s2"])
    expect_identical(original$mu1_true, deviated$mu1_true)
    if (mechanism == "treated_arm") expect_identical(original$Y_0, deviated$Y_0)
  }
  zero <- make_data(list(covariate_shift = 0, source_treatment_scale = 0))
  expect_equal(zero$p_treat_true[zero$R != "t"], rep(.5, 2000L))
  expect_equal(as.numeric(zero$oracle$source_weights[zero$R != "t", ]), rep(2, 4000L), tolerance = 1e-12)
  expect_error(make_data(list(source_treatment_scale = 2)), "source_treatment_scale")
  expect_error(make_data(list(unknown = 1)), "dgp_control")
  expect_error(generate_bounded_data(p = 3), "p >= 4")
  expect_error(generate_bounded_data(outcome_type = "continuous"), "binary")
})

test_that("bounded DGP controls reach simulation output and separate summaries", {
  args <- list(sim_id = 2L, K = 1L, p = 4L, config = "C1", dgp_type = "bounded",
    n_target = 1000L, n_source_sizes = 1000L, methods = c("one_round_crossfit", "target_only"),
    estimate_ate = TRUE, n_folds = 4L, nlambda_init = 4L, verbose = FALSE,
    n_cores_internal = 1L, include_quadratic_bias_rule = FALSE,
    target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
    calibration_layout = "compact", nuisance_solver = "proximal_newton",
    M_tau = 12, M_tau_inference = 12,
    calibration_control = list(recipe = "score_derivative",
      target_propensity_initialization = "calibrated", target_radius = 12))
  a <- do.call(run_single_simulation, args)
  b <- do.call(run_single_simulation, c(args, list(dgp_control = list(covariate_shift = 2))))
  expect_identical(a$estimate[a$method == "target_only_ate"], b$estimate[b$method == "target_only_ate"])
  expect_identical(a$se[a$method == "target_only_ate"], b$se[b$method == "target_only_ate"])
  expect_true(all(a$dgp_covariate_shift == 1))
  expect_true(all(b$dgp_covariate_shift == 2))
  expect_true(all(a$working_dimension == 4L))
  summary <- summarize_results(rbind(a, b))
  expect_equal(nrow(summary), 2L * nrow(a))
  expect_true(all(summary$n_success == 1L))
  expect_error(do.call(run_single_simulation,
    modifyList(args, list(dgp_type = "face", dgp_control = list(covariate_shift = 1)))), "bounded")
})
