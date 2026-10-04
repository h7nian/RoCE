test_that("real-data adapter preserves labels, rows and shared feature maps", {
  helper <- file.path(repo_root, "scripts/real_data/helpers.R")
  skip_if_not(file.exists(helper), "Repository-only collaborator adapter")
  scope <- new.env(parent = globalenv())
  sys.source(helper, envir = scope)
  frame <- data.frame(site = c("other", "target", "third", "target"),
                      A = c(1, 0, 0, 1), Y = c(0, 1, 0, 1), x1 = 1:4, x2 = 4:1)
  build <- function(data = frame, features = c("x2", "x1")) {
    scope$site_data_from_frame(data, "site", "A", "Y", features, "target")
  }
  result <- build()
  expect_identical(names(result), c("t", "s1", "s2"))
  expect_identical(unname(attr(result, "site_mapping")), c("target", "other", "third"))
  expect_identical(result$t$A, frame$A[c(2, 4)])
  expected <- as.matrix(frame[c(2, 4), c("x2", "x1")])
  rownames(expected) <- NULL
  expect_identical(result$t$W_outcome, expected)
  expect_identical(result$t$W_outcome, result$t$Z_site)
  expect_error(build(features = c("x1", "A")), "distinct")
  frame$x1[1] <- NA_real_
  expect_error(build(frame), "finite")
  frame$x1 <- factor(1:4)
  expect_error(build(frame), "numeric")
})
