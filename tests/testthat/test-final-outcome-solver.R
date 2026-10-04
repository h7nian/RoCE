test_that("explicit solver arguments apply to final fits and restore process state", {
  set.seed(91)
  x <- matrix(rnorm(600), 200, 3)
  a <- rep(0:1, 100)
  y <- x[, 1] + .5 * x[, 2] + rnorm(200)
  previous <- nuisance_solver_cpp()
  fits <- lapply(c("coordinate_descent", "proximal_newton"), function(solver) {
    fit <- fit_unified_outcome(
      x, y, a, A_val = 1L, gamma_s = numeric(4), lambda = .01,
      family = "gaussian", Z_site = x, tol = 1e-10,
      max_iter = 1000L, nuisance_solver = solver
    )
    expect_identical(attr(fit, "solver"), solver)
    expect_identical(nuisance_solver_cpp(), previous)
    fit
  })
  expect_equal(as.numeric(fits[[1]]), as.numeric(fits[[2]]), tolerance = 1e-7)
  expect_error(fit_unified_outcome(
    x, y, a, A_val = 2L, gamma_s = numeric(4), Z_site = x,
    nuisance_solver = "coordinate_descent"
  ), "A_val")
  expect_identical(nuisance_solver_cpp(), previous)
})

test_that("final outcome refits dispatch the requested solver in fresh R processes", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("ROCE_TEST_INSTALLED"), "1"),
              "Fresh processes require the current installed test library")
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(c(
    "suppressPackageStartupMessages(library(RoCE))",
    "set.seed(917)",
    "x <- matrix(rnorm(1200), 300, 4)",
    "x[,2] <- .9*x[,1] + .3*x[,2]",
    "weights <- runif(300, .5, 1.5)",
    "y <- 1 + x[,1] - x[,2] + .6*x[,3] + rnorm(300, sd=.2)",
    "fit <- RoCE:::fit_general_glm_cpp(x,y,weights,0L,0L,.005,1000L,1e-10,numeric())",
    "design <- cbind(1,x)",
    "gradient <- colMeans(design * (weights * (drop(design %*% fit$alpha)-y)))",
    "slopes <- as.numeric(fit$alpha)[-1L]",
    "kkt <- max(abs(gradient[1]), ifelse(abs(slopes)>1e-6, abs(gradient[-1]+.005*sign(slopes)), pmax(abs(gradient[-1])-.005,0)))",
    "stopifnot(fit$converged, kkt < 1e-7)",
    "stopifnot(identical(fit$solver, Sys.getenv('ROCE_NUISANCE_SOLVER')))",
    "saveRDS(fit,commandArgs(trailingOnly=TRUE)[1])"
  ), script)
  results <- lapply(c("coordinate_descent", "proximal_newton"), function(solver) {
    output <- tempfile(fileext = ".rds")
    log <- tempfile(fileext = ".log")
    on.exit(unlink(c(output, log)), add = TRUE)
    status <- system2(file.path(R.home("bin"), "Rscript"),
                      c("--vanilla", shQuote(script), shQuote(output)),
                      env = c(paste0("ROCE_NUISANCE_SOLVER=", solver),
                              paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))),
                      stdout = log, stderr = log)
    expect_identical(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
    if (status != 0L) return(NULL)
    readRDS(output)
  })
  if (any(vapply(results, is.null, logical(1L)))) return(invisible(NULL))
  expect_equal(as.numeric(results[[1]]$alpha), as.numeric(results[[2]]$alpha), tolerance = 1e-7)
  expect_false(identical(results[[1]]$iterations, results[[2]]$iterations))
})
