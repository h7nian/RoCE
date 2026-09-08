library(testthat)
source(if (file.exists("source_final_projection.R")) "source_final_projection.R" else
  "diagnosis/tate_common_weight/source_final_projection.R")
test_that("noise penalties retain separate site denominators", {
  target <- matrix(c(1, 2, 4, 0, 0, 0), 3, 2)
  source <- matrix(c(2, 4, 6, 8, 10, 1, 3, 5, 7, 9), 5, 2)
  expected <- sqrt(log(4)*(apply(target,2,var)/3+apply(source,2,var)/5))
  expect_equal(.source_projection_noise_penalty(target, source, 1), expected)
  expect_equal(.source_projection_noise_penalty(source, target, 1), expected)
  expect_equal(.source_projection_noise_penalty(target, source, 2), 2*expected)
  expect_equal(.source_projection_noise_penalty(target+7, source-4, 1), expected)
  expect_error(.source_projection_noise_penalty(target, source, 0), "penalty_scale")
  expect_error(.source_projection_noise_penalty(target[1,,drop=FALSE], source, 1), "derivative rows")
  expect_error(.source_projection_noise_penalty(target, source[,1,drop=FALSE], 1), "derivative rows")
  bad <- target; bad[1,1] <- NA_real_
  expect_error(.source_projection_noise_penalty(bad, source, 1), "derivative rows")
})
