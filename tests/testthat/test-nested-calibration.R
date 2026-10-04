test_that("target calibration solves Hou's derivative-weighted empirical losses", {
  set.seed(907)
  x <- matrix(rnorm(3000), 1000, 3)
  a <- rbinom(1000, 1, plogis(.2 + .2 * x[, 1]))
  y <- rbinom(1000, 1, plogis(-.1 + .3 * a + .2 * x[, 2]))
  data <- list(W_outcome = x, Z_site = x, A = a, Y = y, n = 1000L)
  folds <- partition_into_folds(data, 4L, seed = 91)
  for (arm in 0:1) {
    fit <- .fit_target_calibration(folds, 3:4, arm, "binomial", 5, 6L,
                                   layout = "compact")
    blocks <- lapply(3:4, function(fold) materialize_fold(folds, fold))
    design <- do.call(rbind, lapply(blocks, function(block) cbind(1, block$Z_site)))
    outcome <- unlist(lapply(blocks, `[[`, "Y"))
    treatment <- unlist(lapply(blocks, `[[`, "A"))
    initial_or <- unlist(Map(function(block, coef) drop(cbind(1, block$W_outcome) %*% coef),
                            blocks, fit$initial_outcome), use.names = FALSE)
    initial_ps <- unlist(Map(function(block, coef) drop(cbind(1, block$Z_site) %*% coef),
                            blocks, fit$initial_weight), use.names = FALSE)
    h <- plogis(pmax(-5, pmin(5, initial_or)))
    h <- h * (1 - h)
    selected <- treatment == arm
    theta <- drop(design %*% fit$weight)
    beta <- drop(design %*% fit$outcome)
    ps_gradient <- colMeans(design * h *
      ((1 - selected) - selected * exp(-pmax(-5, pmin(5, theta)))))
    or_gradient <- colMeans(design * selected * exp(-pmax(-5, pmin(5, initial_ps))) *
                            (plogis(beta) - outcome))
    kkt <- function(gradient, coefficients) {
      penalty <- attr(coefficients, "lambda_used")
      slopes <- as.numeric(coefficients)[-1L]
      residual <- ifelse(abs(slopes) > 1e-8, gradient[-1L] + penalty * sign(slopes),
                         pmax(abs(gradient[-1L]) - penalty, 0))
      max(abs(c(gradient[1L], residual)))
    }
    expect_lt(kkt(ps_gradient, fit$weight), 2e-5)
    expect_lt(kkt(or_gradient, fit$outcome), 2e-5)
    evaluated <- .evaluate_target_calibration(fit, materialize_fold(folds, 2L), arm,
                                               "binomial", 5)
    evaluation <- materialize_fold(folds, 2L)
    expected <- evaluated$m_pred + (evaluation$A == arm) *
      (evaluation$Y - evaluated$m_pred) / evaluated$arm_propensity
    expect_equal(evaluated$estimate, mean(expected), tolerance = 1e-14)
    expect_equal(evaluated$varphi_ot, as.numeric(expected - mean(expected)), tolerance = 1e-14)
    expect_identical(attr(fit$initial_outcome$k2_3, "training_folds"), 4L)
    expect_identical(attr(fit$initial_weight$k2_4, "training_folds"), 3L)
  }
})

test_that("nested source and target training excludes both evaluation folds", {
  skip_on_cran()
  set.seed(918)
  site <- function() {
    x <- matrix(rnorm(3000), 1000, 3)
    a <- rep(0:1, 500)
    list(W_outcome = x, Z_site = x, A = a,
         Y = .2 * a + .3 * x[, 1] + rnorm(1000), n = 1000L)
  }
  data <- list(t = site(), s1 = site())
  folds <- build_crossfit_folds(data, 4L)
  changed <- data
  for (name in names(data)) {
    views <- if (name == "t") folds$target_folds else folds$source_folds[[name]]
    rows <- unlist(lapply(views[1:2], `[[`, "original_idx"))
    changed[[name]]$Y[rows] <- changed[[name]]$Y[rows] + 10
    changed[[name]]$W_outcome[rows, 1] <- changed[[name]]$W_outcome[rows, 1] + .5
    changed[[name]]$Z_site[rows, 1] <- changed[[name]]$Z_site[rows, 1] + .5
  }
  changed_folds <- build_crossfit_folds(changed, 4L)
  for (mode in c("one_round", "two_round")) for (arm in 0:1) {
    train <- function(partition, layout, cache) {
      messages <- .prepare_source_calibration_messages(
        partition$target_folds, partition$source_folds$s1, "s1", 3:4,
        mode, arm, "gaussian", 5, 6L, "min", cache
      )
      expect_identical(messages$k2_3$initial_training_folds, 4L)
      expect_identical(messages$k2_4$excluded_folds, 1:2)
      list(source = .fit_source_calibration(
        partition$source_folds$s1, "s1", 3:4,
        function(site, fold) messages[[paste0("k2_", fold)]],
        arm, 0L, 0L, 5, 6L, 500L, "min", cache, layout = layout),
        target = .fit_target_calibration(partition$target_folds, 3:4, arm,
          "gaussian", 5, 6L, fit_cache = cache, layout = layout))
    }
    cache <- new.env(parent = emptyenv())
    original <- train(folds, "block", cache)
    perturbed <- train(changed_folds, "block", cache)
    compact <- train(folds, "compact", NULL)
    for (site in c("source", "target")) for (model in c("weight", "outcome")) {
      expect_identical(as.numeric(original[[site]][[model]]), as.numeric(perturbed[[site]][[model]]))
      expect_equal(as.numeric(original[[site]][[model]]), as.numeric(compact[[site]][[model]]), tolerance = 1e-7)
      expect_equal(attr(original[[site]][[model]], "lambda_used"),
                   attr(compact[[site]][[model]], "lambda_used"), tolerance = 1e-12)
    }
  }
})
