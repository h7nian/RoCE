library(testthat)
solver_path <- if (file.exists("sparse_moment_projection.R")) "sparse_moment_projection.R" else
  "diagnosis/tate_common_weight/sparse_moment_projection.R"
source(solver_path)

test_that("diagonal solution and zero-penalty limit are exact", {
  H <- diag(c(1, 2, 4))
  d <- c(-3, .1, 2)
  lambda <- c(.2, .2, 0)
  fitted <- .solve_sparse_moment_projection(H, d, lambda)
  expected <- sign(d)*pmax(abs(d)-lambda, 0)/diag(H)
  expect_equal(fitted$coefficients, expected, tolerance = 1e-12)
  expect_lte(fitted$max_kkt_residual, 1e-12)
  expect_equal(.solve_sparse_moment_projection(H, d, 0)$coefficients,
               drop(solve(H, d)), tolerance = 1e-12)
  expect_identical(fitted$numerical_psd_adjustment, 0)
  expect_identical(names(fitted), names(.solve_sparse_moment_projection(
    matrix(0, 3, 3), rep(0, 3), .1)))
})

test_that("nondiagonal zero-penalty solution matches independent dense solve", {
  H <- matrix(c(2, .3, .3, 1), 2)
  d <- c(.4, -.6)
  fitted <- .solve_sparse_moment_projection(H, d, 0, tolerance = 1e-11)
  expect_equal(fitted$coefficients, drop(solve(H, d)), tolerance = 1e-9)
  expect_lte(fitted$relative_kkt_residual, 1e-11)
})

test_that("common rescaling of the objective does not alter coefficients", {
  H <- matrix(c(2, .3, .3, 1), 2); d <- c(.4, -.6)
  reference <- .solve_sparse_moment_projection(H, d, .1)
  for (scale in c(1e-8, 1e8)) {
    fitted <- .solve_sparse_moment_projection(scale*H, scale*d, scale*.1)
    expect_equal(fitted$coefficients, reference$coefficients, tolerance = 1e-10)
    expect_lte(fitted$relative_kkt_residual, 1e-9)
  }
})

test_that("p greater than sample rank is supported without a positive ridge", {
  design <- cbind(diag(20), diag(20), diag(20))
  H <- crossprod(design)
  truth <- numeric(60); truth[c(1, 4, 9)] <- c(1, -.8, .6)
  d <- drop(H %*% truth)
  fitted <- .solve_sparse_moment_projection(H, d, .1)
  expect_equal(fitted$numerical_rank, 20L)
  expect_lte(fitted$relative_kkt_residual, 1e-9)
  expect_equal(drop(design %*% fitted$coefficients),
               sign(truth[1:20])*pmax(abs(truth[1:20])-.1, 0), tolerance = 1e-8)
  expect_lte(fitted$numerical_psd_adjustment, 1e-12)
  expect_lte(fitted$max_original_kkt_residual, 1e-8)
})

test_that("zero objective, unsupported directions and invalid inputs fail honestly", {
  expect_equal(.solve_sparse_moment_projection(matrix(0, 2, 2), c(.1, -.1), .2)$coefficients, c(0, 0))
  expect_error(.solve_sparse_moment_projection(matrix(0, 2, 2), c(1, 0), .2), "zero curvature")
  expect_error(.solve_sparse_moment_projection(matrix(1, 2, 2), c(1, -1), .1), "near-null")
  expect_error(.solve_sparse_moment_projection(diag(c(1, -.1)), c(0, 0), .1), "positive semidefinite")
  expect_error(.solve_sparse_moment_projection(matrix(c(1, .1, .2, 1), 2), c(0, 0), .1), "symmetric")
  expect_error(.solve_sparse_moment_projection(diag(2), c(NA, 0), .1), "gradient")
  expect_error(.solve_sparse_moment_projection(diag(2), c(1, 0), -.1), "penalty")
  expect_error(.solve_sparse_moment_projection(diag(2), c(1, 0), .1, max_iterations = 1.5), "iteration")
  error <- tryCatch(.solve_sparse_moment_projection(matrix(c(1, .99, .99, 1), 2),
    c(1, -1), .01, max_iterations = 1), error = identity)
  expect_s3_class(error, "roce_projection_error")
  expect_true(length(error$diagnostics$coefficients) == 2L)
  expect_gt(error$diagnostics$relative_kkt_residual, 1e-9)
})

test_that("four-block adjoint system retains initial-fit coupling and correct signs", {
  blocks <- c("alpha_initial", "gamma_initial", "alpha_final", "gamma_final")
  curvatures <- setNames(lapply(1:4, function(i) diag(c(i+1, i+2))), blocks)
  couplings <- list(alpha_final_gamma_initial = matrix(c(.1, .2, -.3, .1), 2),
                    gamma_final_alpha_initial = matrix(c(.2, -.1, .1, .4), 2))
  gradients <- list(alpha_final = c(.6, -.3), gamma_final = c(.2, .7))
  J <- matrix(0, 8, 8)
  for (i in 1:4) {
    index <- (2*i-1):(2*i)
    J[index, index] <- c(-1, 1, -1, 1)[i]*curvatures[[i]]
  }
  J[5:6, 3:4] <- couplings$alpha_final_gamma_initial
  J[7:8, 1:2] <- couplings$gamma_final_alpha_initial
  d <- c(0, 0, 0, 0, gradients$alpha_final, gradients$gamma_final)
  fitted <- .solve_joint_nuisance_projection(curvatures, couplings, gradients,
                                             setNames(as.list(rep(0, 4)), blocks))
  expect_equal(fitted$coefficients, drop(solve(t(J), d)), tolerance = 1e-10)
  expect_true(any(abs(fitted$coefficients[1:4]) > .001))
  regularized <- .solve_joint_nuisance_projection(curvatures, couplings, gradients,
                                                  setNames(as.list(rep(.02, 4)), blocks))
  expect_lte(max(abs(drop(t(J) %*% regularized$coefficients)-d)), .02+1e-9)
  couplings$alpha_final_gamma_initial <- matrix(0, 3, 2)
  expect_error(.solve_joint_nuisance_projection(curvatures, couplings, gradients,
    setNames(as.list(rep(.02, 4)), blocks)), "coupling")
})
