test_that("source-stage lambda rules validate and preserve inherited defaults", {
  expect_identical(.resolve_source_lambda_rules(NULL, "min"),
    list(initial_weight = "min", weight = "min", outcome = "min"))
  expect_identical(.resolve_source_lambda_rules(list(weight = "1se"), "min"),
    list(initial_weight = "min", weight = "1se", outcome = "min"))
  expect_identical(.resolve_source_lambda_rules(list(outcome = "min"), "1se"),
    list(initial_weight = "1se", weight = "1se", outcome = "min"))
  for (bad in list(list(), list(wieght = "1se"), list(weight = NA_character_),
                  list(weight = c("min", "1se")), list(weight = "other"),
                  setNames(list("min", "1se"), c("weight", "weight")))) {
    expect_error(.resolve_source_lambda_rules(bad), "source_lambda_rules")
  }
  expect_false("source_lambda_rules" %in% names(.validate_calibration_control()))
  control <- list(recipe = "score_derivative", source_lambda_rules = list(weight = "1se"))
  expect_identical(.validate_calibration_control(control, "lasso")$source_lambda_rules, list(weight = "1se"))
  expect_error(.validate_calibration_control(list(source_lambda_rules = list(weight = "1se"))), "requires")
  control$source_nuisance_method <- "standard"
  expect_error(.validate_calibration_control(control), "requires")
})

test_that("source final-weight and OR selection act on their own calibration losses", {
  set.seed(92701)
  make_site <- function(shift) {
    x <- matrix(rnorm(1600), 400, 4)
    x[, 1L] <- x[, 1L] + shift
    a <- rbinom(400, 1, plogis(.2 + .5*x[, 1L] - .2*x[, 2L]))
    list(n = 400L, W_outcome = x, Z_site = x, A = a,
         Y = rbinom(400, 1, plogis(-.6 + .5*a + .4*x[, 2L] - .3*x[, 3L])))
  }
  data <- list(t = make_site(0), s1 = make_site(.7))
  folds <- build_crossfit_folds(data, 4L)
  for (mode in c("one_round", "two_round")) for (arm in 0:1) {
    train <- function(rules = NULL) .fit_complete_source_program(
      folds$target_folds, folds$source_folds$s1, "s1", 2:4, mode, arm,
      "binomial", 5, 8L, 1000L, "min", layout = "compact", tol = 1e-10,
      calibration_control = list(recipe = "score_derivative", source_lambda_rules = rules))
    base <- train()
    weight <- train(list(weight = "1se"))
    outcome <- train(list(outcome = "1se"))
    initial_weight <- train(list(initial_weight = "1se"))
    all_source <- train(list(initial_weight = "1se", weight = "1se", outcome = "1se"))
    numeric_initial <- function(values) lapply(values, as.numeric)
    expect_equal(numeric_initial(weight$initial_outcome), numeric_initial(base$initial_outcome), tolerance = 1e-12)
    expect_equal(numeric_initial(weight$initial_weight), numeric_initial(base$initial_weight), tolerance = 1e-12)
    expect_equal(as.numeric(weight$outcome), as.numeric(base$outcome), tolerance = 1e-12)
    expect_equal(as.numeric(outcome$weight), as.numeric(base$weight), tolerance = 1e-12)
    expect_identical(attr(weight$weight, "lambda_rule"), "1se")
    expect_identical(attr(weight$outcome, "lambda_rule"), "min")
    expect_identical(attr(outcome$outcome, "lambda_rule"), "1se")
    expect_gte(attr(weight$weight, "lambda_used") + 1e-12, attr(base$weight, "lambda_used"))
    expect_gte(attr(outcome$outcome, "lambda_used") + 1e-12, attr(base$outcome, "lambda_used"))
    expect_true(all(vapply(initial_weight$initial_weight, function(value) identical(attr(value, "lambda_rule"), "1se"), logical(1L))))
    expect_identical(attr(initial_weight$weight, "lambda_rule"), "min")
    expect_identical(attr(initial_weight$outcome, "lambda_rule"), "min")
    expect_equal(numeric_initial(all_source$initial_outcome), numeric_initial(base$initial_outcome), tolerance = 1e-12)
    expect_true(all(vapply(all_source$initial_weight, function(value) identical(attr(value, "lambda_rule"), "1se"), logical(1L))))
    expect_identical(attr(all_source$weight, "lambda_rule"), "1se")
    expect_identical(attr(all_source$outcome, "lambda_rule"), "1se")
  }
})

test_that("source-stage rules respect nested validation boundaries", {
  set.seed(92702)
  make_site <- function() {
    x <- matrix(rnorm(2500), 500, 5)
    list(n = 500L, W_outcome = x, Z_site = x, A = rep(0:1, 250),
         Y = rbinom(500, 1, plogis(.3*x[, 1L] - .4*x[, 2L])))
  }
  data <- list(t = make_site(), s1 = make_site())
  folds <- build_crossfit_folds(data, 5L)
  changed <- data
  for (site in names(data)) {
    partition <- if (site == "t") folds$target_folds else folds$source_folds[[site]]
    rows <- unlist(lapply(partition[1:2], `[[`, "original_idx"), use.names = FALSE)
    changed[[site]]$Y[rows] <- 1 - changed[[site]]$Y[rows]
    changed[[site]]$W_outcome[rows, 1L] <- changed[[site]]$W_outcome[rows, 1L] + 20
    changed[[site]]$Z_site <- changed[[site]]$W_outcome
  }
  train <- function(partition) .fit_complete_source_program(
    partition$target_folds, partition$source_folds$s1, "s1", 3:5, "one_round", 1L,
    "binomial", 5, 8L, 1000L, "min", layout = "compact", tol = 1e-10,
    calibration_control = list(recipe = "score_derivative", source_lambda_rules = list(weight = "1se")))
  first <- train(folds)
  second <- train(build_crossfit_folds(changed, 5L))
  expect_equal(as.numeric(first$weight), as.numeric(second$weight), tolerance = 1e-12)
  expect_equal(as.numeric(first$outcome), as.numeric(second$outcome), tolerance = 1e-12)
  expect_identical(first$training_folds, 3:5)
  expect_identical(attr(first$weight, "lambda_rule"), "1se")
})

test_that("two-layer TATE source overrides leave target and the other final loss fixed", {
  set.seed(92703)
  generated <- generate_simulation_data(n_total = 600L, K = 1L, p = 4L,
    config = "C1", dgp_type = "bounded", outcome_type = "binary",
    estimand_type = "superpopulation", n_target = 300L, n_source_sizes = 300L)
  data <- split_data_by_site(generated)
  train <- function(rules = NULL, training_data = data) run_tate_crossfit(training_data, n_folds = 3L,
    communication_mode = "one_round", target_nuisance_method = "hou_calibrated",
    crossfit_layers = 2L, aggregation_mode = "joint_tate", nlambda_init = 6L,
    nuisance_lambda_rule = "min", nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    M_tau = 5, M_tau_inference = 5, lambda_selection = .5,
    calibration_layout = "compact", n_cores = 1L, verbose = FALSE,
    calibration_control = list(recipe = "score_derivative", target_radius = 5,
      source_lambda_rules = rules))
  base <- train()
  changed <- train(list(weight = "1se"))
  expect_equal(changed$target_only, base$target_only, tolerance = 1e-12)
  for (arm in c("mu1", "mu0")) for (fold in 1:3) {
    original <- base$arm_results[[arm]]$fold_results[[fold]]
    selected <- changed$arm_results[[arm]]$fold_results[[fold]]
    expect_equal(selected$target_only$varphi_ot, original$target_only$varphi_ot, tolerance = 1e-12)
    expect_equal(as.numeric(selected$source_results$s1$alpha_ts),
                 as.numeric(original$source_results$s1$alpha_ts), tolerance = 1e-12)
    expect_identical(attr(selected$source_results$s1$gamma_s, "lambda_rule"), "1se")
  }
  shifted <- data
  rows <- head(which(shifted$s1$A == 1L), 5L)
  shifted$s1$Y[rows] <- 1 - shifted$s1$Y[rows]
  reused <- .reuse_one_round_tate_across_rho(data, shifted, changed,
    changed_sources = "s1", lambda_selection = .5, verbose = FALSE, n_cores = 1L)
  fresh <- train(list(weight = "1se"), training_data = shifted)
  expect_equal(reused$estimate, fresh$estimate, tolerance = 1e-10)
  expect_equal(reused$se, fresh$se, tolerance = 1e-10)
  for (fold in 1:3) {
    expect_identical(attr(reused$arm_results$mu1$fold_results[[fold]]$source_results$s1$gamma_s, "lambda_rule"), "1se")
  }
})
